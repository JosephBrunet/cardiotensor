import mmap
import math
import os
import shutil
import time
import warnings
from collections.abc import Sequence
from multiprocessing.pool import ThreadPool
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from alive_progress import alive_bar
from matplotlib import pyplot as plt

# from memory_profiler import profile
from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.orientation.orientation_computation_functions import (
    adjust_start_end_index,
    calculate_center_vector,
    calculate_structure_tensor,
    compute_azimuth_and_elevation,
    compute_fraction_anisotropy,
    compute_helical_and_intrusion_angles,
    interpolate_points,
    plot_images,
    rotate_vectors_to_new_axis,
    write_images,
)
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.image_io import (
    ZarrVectorFieldStore,
    initialize_zarr_vector_field,
    normalize_vector_format,
    vector_field_path,
    write_vector_field,
)
from cardiotensor.utils.utils import (
    get_available_cpu_count,
    get_available_memory_bytes,
    remove_corrupted_files,
)


# --- small helpers ---
def _resolve_colormap(colormap: str | None):
    """Resolve colormap names, including the project-specific helical angle map."""
    if colormap is None:
        return None
    cmap_name = colormap.strip()
    if not cmap_name:
        return None
    if cmap_name.lower() in {"helix_angle", "helix_angle_cmap"}:
        return helix_angle_cmap
    try:
        return plt.get_cmap(cmap_name)
    except ValueError as err:
        raise ValueError(
            f"Unknown colormap '{colormap}'. Use a valid matplotlib colormap name "
            "or 'helix_angle'."
        ) from err


def _centerline_neighborhood(
    center_line: np.ndarray, global_slice_idx: int, buffer: int = 5
) -> np.ndarray:
    """Return a chunk-independent centerline window around one global slice."""
    if not 0 <= global_slice_idx < len(center_line):
        raise IndexError(
            f"global_slice_idx={global_slice_idx} is outside centerline length "
            f"{len(center_line)}"
        )
    start = max(0, global_slice_idx - buffer)
    end = min(len(center_line), global_slice_idx + buffer + 1)
    return center_line[start:end]


def _normalize_vectors_in_place(vector_field: np.ndarray) -> np.ndarray:
    """Normalize one Z slice at a time, preserving zero-filled masked voxels."""
    for z in range(vector_field.shape[1]):
        vector_slice = vector_field[:, z]
        norms = np.linalg.norm(vector_slice, axis=0)
        np.divide(vector_slice, norms[None], out=vector_slice, where=norms[None] > 0)
        if isinstance(vector_field, np.memmap) and (
            (z + 1) % 4 == 0 or z + 1 == vector_field.shape[1]
        ):
            vector_field.flush()
            mapped = getattr(vector_field, "_mmap", None)
            if (
                mapped is not None
                and hasattr(mapped, "madvise")
                and hasattr(mmap, "MADV_DONTNEED")
            ):
                mapped.madvise(mmap.MADV_DONTNEED)
    return vector_field


def _estimate_low_memory_peak_bytes(
    volume_shape: tuple[int, int, int],
    input_dtype: np.dtype | type,
    *,
    include_eigenvalues: bool,
    has_mask: bool,
    write_angles: bool,
    sigma: float,
    rho: float,
    truncate: float,
    cpu_count: int | None = None,
    include_loaded_inputs: bool = True,
) -> int:
    """Conservative peak estimate for a disk-backed tensor calculation."""
    depth, height, width = volume_shape
    voxels = math.prod(volume_shape)
    output_components = 6 if include_eigenvalues else 3
    mapped_outputs = voxels * output_components * np.dtype(np.float32).itemsize
    input_bytes = voxels * np.dtype(input_dtype).itemsize
    mask_bytes = voxels if has_mask else 0

    # structure-tensor uses float64 blocks with overlap and several component arrays.
    block_size = 200
    radius = 2 * int(max(sigma, rho) * truncate + 0.5)
    block_shape = (
        min(depth, block_size + 2 * radius),
        min(height, block_size + 2 * radius),
        min(width, block_size + 2 * radius),
    )
    block_count = math.prod(math.ceil(size / block_size) for size in volume_shape)
    workers = min(
        max(1, get_available_cpu_count() if cpu_count is None else cpu_count),
        block_count,
    )
    worker_temporaries = workers * math.prod(block_shape) * 18 * 8

    # At least one angle-writing worker must fit. Vector-only writing is much lighter.
    slice_work = height * width * (128 if write_angles else 16)
    reserve = 512 * 1024**2

    if include_loaded_inputs:
        # Mask resizing can temporarily hold the resized mask and its fitted copy.
        load_phase = input_bytes + 2 * mask_bytes
        tensor_phase = mapped_outputs + input_bytes + mask_bytes + worker_temporaries
        post_phase = mapped_outputs + mask_bytes + slice_work
    else:
        # The runtime call is made after input/mask allocations, so available memory
        # already excludes them. Estimate only additional allocations from this point.
        load_phase = 0
        tensor_phase = mapped_outputs + worker_temporaries
        post_phase = mapped_outputs + slice_work

    return max(load_phase, tensor_phase, post_phase) + reserve


def _safe_low_memory_chunk_size(
    requested_chunk: int,
    total_slices: int,
    height: int,
    width: int,
    input_dtype: np.dtype | type,
    *,
    padding: int,
    include_eigenvalues: bool,
    has_mask: bool,
    write_angles: bool,
    sigma: float,
    rho: float,
    truncate: float,
    available_memory_bytes: int | None = None,
    cpu_count: int | None = None,
) -> tuple[int, int, int]:
    """Return safe payload slices plus requested/safe peak byte estimates."""
    if requested_chunk <= 0:
        raise ValueError("N_CHUNK must be > 0")
    available = (
        get_available_memory_bytes()
        if available_memory_bytes is None
        else available_memory_bytes
    )
    budget = int(available * 0.82)
    requested_chunk = min(requested_chunk, total_slices)

    def estimate(payload: int) -> int:
        padded_depth = min(total_slices, payload + 2 * padding)
        return _estimate_low_memory_peak_bytes(
            (padded_depth, height, width),
            input_dtype,
            include_eigenvalues=include_eigenvalues,
            has_mask=has_mask,
            write_angles=write_angles,
            sigma=sigma,
            rho=rho,
            truncate=truncate,
            cpu_count=cpu_count,
        )

    requested_peak = estimate(requested_chunk)
    safe_chunk = requested_chunk
    safe_peak = requested_peak
    while safe_chunk > 1 and safe_peak > budget:
        safe_chunk -= 1
        safe_peak = estimate(safe_chunk)

    if safe_peak > budget:
        raise MemoryError(
            "Even a one-slice low-memory chunk exceeds the conservative memory "
            f"budget ({safe_peak / 1024**3:.2f} GiB estimated, "
            f"{available / 1024**3:.2f} GiB available). Request more memory."
        )
    return safe_chunk, requested_peak, safe_peak


def _memory_aware_worker_count(
    num_slices: int,
    height: int,
    width: int,
    write_angles: bool,
    *,
    available_memory_bytes: int | None = None,
    cpu_count: int | None = None,
) -> int:
    """Choose a conservative worker count from slice size and free memory."""
    if num_slices <= 0:
        return 1
    available_memory_bytes = (
        get_available_memory_bytes()
        if available_memory_bytes is None
        else available_memory_bytes
    )
    cpu_count = get_available_cpu_count() if cpu_count is None else cpu_count

    # Angles create rotated vectors, three outputs, and temporary scalar arrays.
    # Vector-only workers mainly need codec and file-write buffers.
    bytes_per_pixel = 128 if write_angles else 12
    bytes_per_worker = max(1, height * width * bytes_per_pixel)
    memory_workers = max(1, (available_memory_bytes // 2) // bytes_per_worker)
    return max(1, min(num_slices, cpu_count, 32, memory_workers))


def _ensure_structure_tensor_memory(
    volume_shape: tuple[int, ...],
    *,
    include_eigenvalues: bool = True,
    low_memory: bool = False,
    available_memory_bytes: int | None = None,
    input_dtype: np.dtype | type = np.float32,
    has_mask: bool = False,
    write_angles: bool = False,
    sigma: float = 1.0,
    rho: float = 3.0,
    truncate: float = 4.0,
) -> None:
    """Fail before tensor allocation when a job cannot hold its output arrays."""
    voxel_count = math.prod(volume_shape)
    output_components = 6 if include_eigenvalues else 3
    output_bytes = voxel_count * output_components * np.dtype(np.float32).itemsize
    available_memory_bytes = (
        get_available_memory_bytes()
        if available_memory_bytes is None
        else available_memory_bytes
    )

    output_gib = output_bytes / 1024**3
    available_gib = available_memory_bytes / 1024**3
    print(
        f"Structure tensor outputs need {output_gib:.2f} GiB; "
        f"{available_gib:.2f} GiB is currently available"
    )

    if low_memory:
        additional_bytes = _estimate_low_memory_peak_bytes(
            tuple(volume_shape),
            input_dtype,
            include_eigenvalues=include_eigenvalues,
            has_mask=has_mask,
            write_angles=write_angles,
            sigma=sigma,
            rho=rho,
            truncate=truncate,
            include_loaded_inputs=False,
        )
        print(
            "Low-memory mode: tensor outputs will use temporary memory-mapped files; "
            f"conservative additional working set is {additional_bytes / 1024**3:.2f} GiB"
        )
        if additional_bytes > available_memory_bytes:
            raise MemoryError(
                f"This padded low-memory chunk still needs approximately "
                f"{additional_bytes / 1024**3:.2f} GiB beyond currently resident "
                f"inputs, but only {available_gib:.2f} GiB is available. Reduce "
                "N_CHUNK (or --chunk_size for SLURM) or request more memory."
            )
        return

    # Keep 30% for block temporaries, Python, codecs, and image-writing buffers.
    if output_bytes > available_memory_bytes * 0.7:
        raise MemoryError(
            f"Structure tensor outputs need {output_gib:.2f} GiB, leaving "
            f"insufficient working memory from {available_gib:.2f} GiB available. "
            "Reduce N_CHUNK (or --chunk_size for SLURM) or request more memory."
        )


def check_already_processed(
    output_dir: str,
    start_index: int,
    end_index: int,
    write_vectors: bool,
    write_angles: bool,
    output_format: str,
    vector_format: str = "zarr",
    zarr_store: ZarrVectorFieldStore | None = None,
    angle_names: tuple[str, str] = ("HA", "IA"),
    fa_name: str = "FA",
    extra_expected: Sequence[str] | None = None,
) -> bool:
    """
    Check whether all required output files already exist for every slice index.

    Parameters
    ----------
    output_dir : str
        Base output directory.
    start_index : int
        First global slice index to check (inclusive).
    end_index : int
        Last global slice index to check (exclusive).
    write_vectors : bool
        If True, expect completed vector output in the selected backend.
    write_angles : bool
        If True, expect angle images for angle_names[0], angle_names[1], and FA.
    output_format : str
        Image format/extension for angles, for example "jp2" or "tif".
    angle_names : tuple[str, str], optional
        Names of the two angle outputs, e.g. ("HA", "IA") or ("AZ", "EL").
    fa_name : str, optional
        Name of the FA subfolder, default "FA".
    extra_expected : sequence of str, optional
        Additional per-slice path templates to check. Each template must contain
        "{idx}" which will be formatted as a zero-padded integer (06d), and may
        also contain "{ext}" for the image extension.

    Returns
    -------
    bool
        True if all expected files for all indices exist (and pass the quick
        corruption filter), False otherwise.
    """
    if start_index >= end_index:
        # Nothing to check
        return True

    # Normalize extension
    ext = output_format.lstrip(".")
    if not ext:
        raise ValueError(
            "output_format must be a non-empty extension like 'jp2' or 'tif'."
        )

    # Prepare optional extras
    extra_expected = tuple(extra_expected or ())
    vector_format = normalize_vector_format(vector_format)
    zarr_completed = None
    if write_vectors and vector_format == "zarr":
        if zarr_store is None:
            print(f"Missing Zarr vector store in {output_dir}")
            return False
        zarr_completed = zarr_store.completed_range(start_index, end_index)

    for idx in range(start_index, end_index):
        expected_files = []

        if write_angles:
            a1, a2 = angle_names
            expected_files += [
                os.path.join(output_dir, a1, f"{a1}_{idx:06d}.{ext}"),
                os.path.join(output_dir, a2, f"{a2}_{idx:06d}.{ext}"),
                os.path.join(output_dir, fa_name, f"{fa_name}_{idx:06d}.{ext}"),
            ]

        if write_vectors and vector_format == "npy":
            expected_files.append(
                str(vector_field_path(output_dir, "npy") / f"eigen_vec_{idx:06d}.npy")
            )
        elif write_vectors and not zarr_completed[idx - start_index]:
            print(f"Missing or incomplete Zarr vector output for slice {idx:06d}")
            return False

        # User-specified extras, if any
        for tmpl in extra_expected:
            expected_files.append(tmpl.format(idx=f"{idx:06d}", ext=ext))

        # Remove small/corrupted files before checking (function defined elsewhere)
        corrupted_files = remove_corrupted_files(expected_files)
        if corrupted_files:
            print(
                f"Corrupted output file(s) found for slice {idx:06d}: "
                + ", ".join(corrupted_files)
            )

        # If any file is missing, we need to process
        missing_files = [path for path in expected_files if not os.path.exists(path)]
        if missing_files:
            if not corrupted_files:
                print(
                    f"Missing output file(s) for slice {idx:06d}: "
                    + ", ".join(missing_files)
                )
            return False

    print(f"Checking already processed files: all expected files exist in {output_dir}")
    return True


# --- main API ---
def compute_orientation(
    volume_path: str,
    mask_path: str | None = None,
    output_dir: str = "./output",
    output_format: str = "jp2",
    output_type: str = "8bit",
    sigma: float = 1.0,
    rho: float = 3.0,
    truncate: float = 4.0,
    axis_points: np.ndarray | None = None,
    vertical_padding: float | None = None,
    write_vectors: bool = False,
    angle_mode: str = "ha_ia",
    write_angles: bool = True,
    use_gpu: bool = True,
    gpu_workers_per_device: int = 1,
    is_test: bool = False,
    n_slice_test: int | None = None,
    show_quiver: bool = True,
    start_index: int = 0,
    end_index: int | None = None,
    colormap: str | None = None,
    colormap_angle1: str | None = None,
    colormap_angle2: str | None = None,
    projected: bool = False,
    vector_format: str = "zarr",
    low_memory: bool = False,
    low_memory_dir: str | os.PathLike[str] | None = None,
) -> None:
    """
    Compute the orientation for a volume dataset.

    Args:
        volume_path: Path to the 3D volume.
        mask_path: Optional binary mask path.
        output_dir: Output directory for results.
        output_format: Image format for results.
        output_type: Image type ("8bit" or "rgb").
        sigma: Noise scale for structure tensor.
        rho: Integration scale for structure tensor.
        truncate: Gaussian kernel truncation.
        axis_points: 3D points defining LV axis for cylindrical coordinates.
        vertical_padding: Padding slices for tensor computation.
        write_vectors: Whether to save eigenvectors. Ignored in test mode.
        vector_format: Vector storage backend, either "npy" or "zarr".
        low_memory: Store temporary tensor outputs as memory-mapped files.
        low_memory_dir: Optional scratch root. Defaults to
            OUTPUT_PATH/.cardiotensor_scratch, never the system /tmp directory.
        write_angles: Whether to save HA/IA/FA maps.
        projected: If True in ha_ia mode, write projected HA/IA legacy maps.
        use_gpu: Use GPU acceleration for tensor computation.
        gpu_workers_per_device: Worker processes assigned to each visible GPU.
        is_test: If True, runs in test mode and outputs plots.
        n_slice_test: Number of slices to process in test mode.
        show_quiver: If True, overlay the vector field on the test-slice figure.
        start_index: Start slice index.
        end_index: End slice index (None = last slice).
        colormap: Shared colormap name for RGB angle outputs.
        colormap_angle1: Colormap name for the first angle output.
        colormap_angle2: Colormap name for the second angle output.
    """

    # --- Sanity checks ---
    if sigma > rho:
        raise ValueError("sigma must be <= rho")
    if gpu_workers_per_device < 1:
        raise ValueError("gpu_workers_per_device must be at least 1")
    if not write_vectors and not write_angles:
        raise ValueError("At least one of write_vectors or write_angles must be True")
    vector_format = normalize_vector_format(vector_format)

    if angle_mode.lower() == "ha_ia":
        angle_names = ("HA_projected", "IA_projected") if projected else ("HA", "IA")
    elif angle_mode.lower() == "az_el":
        angle_names = ("AZ", "EL")
    else:
        raise ValueError("ANGLE_MODE must be 'ha_ia' or 'az_el'")

    if is_test and write_vectors:
        warnings.warn(
            "WRITE_VECTORS=True has no effect when is_test=True; "
            "vector fields are not written in test mode. "
            "Proceeding with write_vectors=False.",
            UserWarning,
            stacklevel=2,
        )
        write_vectors = False

    projected_status = projected if angle_mode.lower().strip() == "ha_ia" else "[n/a]"

    print(
        f"""
Parameters:
    - Volume path:    {volume_path}
    - Mask path:      {mask_path or "[None]"}
    - Output dir:     {output_dir}
    - Output format:  {output_format}
    - Output type:    {output_type}
    - sigma / rho:    {sigma} / {rho}
    - truncate:       {truncate}
    - Write angles:   {write_angles}
    - Angle mode:     {angle_mode}  -> {angle_names[0]}, {angle_names[1]}
    - Projected HA/IA:{projected_status}
    - Write vectors:  {write_vectors}
    - Vector format:  {vector_format}
    - Low memory:     {low_memory}
    - Use GPU:        {use_gpu}
    - GPU workers:    {gpu_workers_per_device} per device
    - Test mode:      {is_test}
    - Show quiver:    {show_quiver}
    - Colormap:       {colormap or "[default]"}
    - Colormap angle1:{colormap_angle1 or "[default]"}
    - Colormap angle2:{colormap_angle2 or "[default]"}
    """
    )

    print("\n" + "-" * 40)
    print("READING VOLUME INFORMATION")
    print("-" * 40 + "\n")

    print(f"Volume path: {volume_path}")

    data_reader = DataReader(volume_path)
    vector_store = None
    if write_vectors and vector_format == "zarr":
        vector_store = initialize_zarr_vector_field(
            output_dir, tuple(data_reader.shape)
        )

    if end_index is None:
        end_index = data_reader.shape[0]

    print(f"Number of slices: {data_reader.shape[0]}")

    # --- Check if already processed ---
    print("Check if file is already processed...")
    if (
        check_already_processed(
            output_dir,
            start_index,
            end_index,
            write_vectors,
            write_angles,
            output_format,
            vector_format=vector_format,
            zarr_store=vector_store,
            angle_names=angle_names,
        )
        and not is_test
    ):
        print("\nAll images are already processed. Skipping computation.\n")
        return

    print("\n---------------------------------")
    print("CALCULATE CENTER LINE\n")
    center_line = interpolate_points(axis_points, data_reader.shape[0])

    print("\n---------------------------------")
    print("CALCULATE PADDING START AND ENDING INDEXES\n")

    if vertical_padding is None:
        vertical_padding = int(sigma * truncate + 0.5) + int(rho * truncate + 0.5)

    padding_start = padding_end = math.ceil(vertical_padding)
    if not is_test:
        if padding_start > start_index:
            padding_start = start_index
        if padding_end > (data_reader.shape[0] - end_index):
            padding_end = data_reader.shape[0] - end_index
    if is_test:
        if n_slice_test is None:
            raise ValueError("n_slice_test is required in test mode")
        if not 0 <= n_slice_test < data_reader.shape[0]:
            raise ValueError(
                f"n_slice_test={n_slice_test} is outside volume bounds "
                f"[0, {data_reader.shape[0]})"
            )

    print(f"Padding start, Padding end : {padding_start}, {padding_end}")
    start_index_padded, end_index_padded = adjust_start_end_index(
        start_index,
        end_index,
        data_reader.shape[0],
        padding_start,
        padding_end,
        is_test,
        n_slice_test,
    )
    print(
        f"Start index padded, End index padded : {start_index_padded}, {end_index_padded}"
    )

    print("\n---------------------------------")
    print("LOAD DATASET\n")
    volume = data_reader.load_volume(start_index_padded, end_index_padded)
    print(f"Loaded volume shape {volume.shape}")
    invalid_mask = None
    if mask_path is not None:
        print("\n---------------------------------")
        print("LOAD MASK\n")
        mask_reader = DataReader(mask_path)

        mask = mask_reader.load_volume(
            start_index_padded, end_index_padded, unbinned_shape=data_reader.shape
        )

        if mask.shape != volume.shape:
            raise ValueError(
                f"Mask shape {mask.shape} does not match volume shape {volume.shape}"
            )

        invalid_mask = mask == 0
        volume[invalid_mask] = 0
        del mask

    print("\n" + "-" * 40)
    print("CALCULATING STRUCTURE TENSOR")
    print("-" * 40 + "\n")
    need_eigenvalues = write_angles or is_test
    _ensure_structure_tensor_memory(
        tuple(volume.shape),
        include_eigenvalues=need_eigenvalues,
        low_memory=low_memory,
        input_dtype=volume.dtype,
        has_mask=invalid_mask is not None,
        write_angles=write_angles,
        sigma=sigma,
        rho=rho,
        truncate=truncate,
    )
    tensor_scratch = None
    if low_memory:
        scratch_root = Path(
            low_memory_dir or Path(output_dir) / ".cardiotensor_scratch"
        ).expanduser().resolve()
        scratch_root.mkdir(parents=True, exist_ok=True)
        output_components = 6 if need_eigenvalues else 3
        required_scratch = (
            math.prod(volume.shape) * output_components * np.dtype(np.float32).itemsize
        )
        free_scratch = shutil.disk_usage(scratch_root).free
        print(
            f"Tensor scratch root: {scratch_root} | required: "
            f"{required_scratch / 1024**3:.2f} GiB | free: "
            f"{free_scratch / 1024**3:.2f} GiB"
        )
        if required_scratch > free_scratch * 0.9:
            raise OSError(
                f"Low-memory tensor outputs need {required_scratch / 1024**3:.2f} "
                f"GiB in {scratch_root}, but only {free_scratch / 1024**3:.2f} "
                "GiB is free. Choose LOW_MEMORY_DIR with more space."
            )
        tensor_scratch = TemporaryDirectory(
            prefix="cardiotensor_tensor_", dir=scratch_root
        )
        print(f"Tensor scratch directory: {tensor_scratch.name}")
    t1 = time.perf_counter()  # start time
    val, vec = calculate_structure_tensor(
        volume,
        sigma,
        rho,
        truncate=truncate,
        use_gpu=use_gpu,
        gpu_workers_per_device=gpu_workers_per_device,
        return_eigenvalues=need_eigenvalues,
        memmap_dir=None if tensor_scratch is None else tensor_scratch.name,
    )

    print(f"Vector field shape: {vec.shape}")
    if not is_test:
        # The intensity data is no longer needed once tensor calculation ends.
        del volume
        volume = None

    array_end = vec.shape[1] - padding_end
    interior = slice(padding_start, array_end)
    if low_memory:
        vec = vec[:, interior, :, :]
        if val is not None:
            val = val[:, interior, :, :]
    else:
        cropped_vec = vec[:, interior, :, :].copy(order="C")
        del vec
        vec = cropped_vec
        if val is not None:
            cropped_val = val[:, interior, :, :].copy(order="C")
            del val
            val = cropped_val

    if is_test:
        cropped_volume = volume[interior, :, :].copy(order="C")
        del volume
        volume = cropped_volume
    print(f"Vector shape after removing padding: {vec.shape}")

    if invalid_mask is not None:
        mask_end = invalid_mask.shape[0] - padding_end
        invalid_mask = invalid_mask[padding_start:mask_end].copy(order="C")
        print("Applying mask to cropped tensors and vectors...")
        if val is not None:
            val[:, invalid_mask] = 0
        vec[:, invalid_mask] = 0
        for mapped_output in (val, vec):
            if isinstance(mapped_output, np.memmap):
                mapped_output.flush()
                mapped = getattr(mapped_output, "_mmap", None)
                if (
                    mapped is not None
                    and hasattr(mapped, "madvise")
                    and hasattr(mmap, "MADV_DONTNEED")
                ):
                    mapped.madvise(mmap.MADV_DONTNEED)
        print("Masking complete")
    else:
        invalid_mask = None

    _normalize_vectors_in_place(vec)

    if write_vectors and vector_format == "zarr":
        vector_completed = vector_store.completed_range(
            start_index, start_index + vec.shape[1]
        )
    elif write_vectors:
        npy_dir = vector_field_path(output_dir, "npy")
        vector_completed = np.array(
            [
                (npy_dir / f"eigen_vec_{idx:06d}.npy").exists()
                for idx in range(start_index, start_index + vec.shape[1])
            ],
            dtype=bool,
        )
    else:
        vector_completed = np.zeros(vec.shape[1], dtype=bool)

    t2 = time.perf_counter()  # stop time
    print(f"finished calculating structure tensors in {t2 - t1} seconds")

    print("\n" + "-" * 40)
    print("ANGLE & ANISOTROPY CALCULATION")
    print("-" * 40 + "\n")

    if not is_test:
        num_slices = vec.shape[1]
        num_workers = _memory_aware_worker_count(
            num_slices,
            height=vec.shape[2],
            width=vec.shape[3],
            write_angles=write_angles,
        )

        print(f"Using {num_workers} memory-aware slice worker(s)")

        def update_bar(_):
            """Callback to tick the progress bar after each finished task."""
            bar()

        with ThreadPool(processes=num_workers) as pool:
            with alive_bar(
                num_slices, title="Processing slices (ThreadPool)", bar="smooth"
            ) as bar:
                results = []
                for z in range(num_slices):
                    global_slice_idx = start_index + z
                    result = pool.apply_async(
                        compute_slice_angles_and_anisotropy,
                        (
                            z,
                            vec[:, z, :, :],
                            None,
                            np.around(center_line[global_slice_idx]),
                            None if val is None else val[:, z, :, :],
                            center_line,
                            output_dir,
                            output_format,
                            output_type,
                            start_index,
                            write_vectors,
                            write_angles,
                            is_test,
                            show_quiver,
                            angle_mode,
                            colormap,
                            colormap_angle1,
                            colormap_angle2,
                            projected,
                            global_slice_idx,
                            None if invalid_mask is None else invalid_mask[z],
                            vector_format,
                            vector_store,
                            bool(vector_completed[z]),
                        ),
                        callback=update_bar,
                    )
                    results.append(result)

                # Ensure all tasks complete before leaving the pool context
                total_compute_time = 0.0
                total_write_time = 0.0
                skipped_slices = 0
                for r in results:
                    compute_time, write_time, skipped = r.get()
                    total_compute_time += compute_time
                    total_write_time += write_time
                    skipped_slices += int(skipped)

        processed_slices = num_slices - skipped_slices
        print(
            "Angle/FA timing: "
            f"compute={total_compute_time:.2f}s, "
            f"write={total_write_time:.2f}s, "
            f"processed={processed_slices}, skipped={skipped_slices}"
        )
    else:
        if vec.shape[1] != 1:
            raise ValueError(
                f"Test mode expected exactly 1 slice after padding removal, got {vec.shape[1]}"
            )

        z = 0
        global_slice_idx = n_slice_test
        compute_slice_angles_and_anisotropy(
            z,
            vec[:, z, :, :],
            # Test mode is the only path that retains the intensity image.
            volume[z, :, :],
            np.around(center_line[global_slice_idx]),
            None if val is None else val[:, z, :, :],
            center_line,
            output_dir,
            output_format,
            output_type,
            start_index,
            write_vectors,
            write_angles,
            is_test,
            show_quiver,
            angle_mode,
            colormap,
            colormap_angle1,
            colormap_angle2,
            projected,
            global_slice_idx,
            None if invalid_mask is None else invalid_mask[z],
            vector_format,
            vector_store,
            bool(vector_completed[z]),
        )

    if is_test:
        print(f"\nFinished processing test slice {n_slice_test}")
    else:
        end_index_local = start_index + vec.shape[1]
        print(f"\nFinished processing slices {start_index} to {end_index_local}")
    print("---------------------------------\n\n")

    if tensor_scratch is not None:
        if isinstance(vec, np.memmap):
            vec.flush()
        if isinstance(val, np.memmap):
            val.flush()
        del vec, val
        tensor_scratch.cleanup()
    return


def compute_slice_angles_and_anisotropy(
    z: int,
    vector_field_slice: np.ndarray,
    img_slice: np.ndarray | None,
    center_point: np.ndarray,
    eigen_val_slice: np.ndarray | None,
    center_line: np.ndarray,
    output_dir: str,
    output_format: str = "jp2",
    output_type: str = "8bit",
    start_index: int = 0,
    write_vectors: bool = False,
    write_angles: bool = True,
    is_test: bool = False,
    show_quiver: bool = True,
    angle_mode: str = "ha_ia",
    colormap: str | None = None,
    colormap_angle1: str | None = None,
    colormap_angle2: str | None = None,
    projected: bool = False,
    global_slice_idx: int | None = None,
    invalid_mask_slice: np.ndarray | None = None,
    vector_format: str = "zarr",
    zarr_store: ZarrVectorFieldStore | None = None,
    vector_already_written: bool = False,
) -> tuple[float, float, bool]:
    """
    Compute either HA/IA or Azimuth/Elevation plus FA for a single slice,
    then plot and/or write outputs depending on flags.
    """
    # Decide angle labels and ranges based on mode
    mode = angle_mode.lower().strip()
    if mode == "ha_ia":
        angle_names = ("HA_projected", "IA_projected") if projected else ("HA", "IA")
        angle_ranges = ((-90.0, 90.0), (-90.0, 90.0))
    elif mode == "az_el":
        angle_names = ("AZ", "EL")
        angle_ranges = ((-180.0, 180.0), (0.0, 90.0))
    else:
        raise ValueError("ANGLE_MODE must be 'ha_ia' or 'az_el'")

    ext = output_format.lstrip(".")
    idx = start_index + z if global_slice_idx is None else global_slice_idx
    shared_colormap = _resolve_colormap(colormap)
    colormap_angle = _resolve_colormap(colormap_angle1) or shared_colormap
    colormap_angle2 = _resolve_colormap(colormap_angle2) or shared_colormap
    if mode == "az_el":
        if colormap_angle is None:
            colormap_angle = helix_angle_cmap
        if colormap_angle2 is None:
            colormap_angle2 = plt.cm.viridis
    rows, cols = vector_field_slice.shape[1:]
    center_x, center_y = float(center_point[0]), float(center_point[1])
    if not (0 <= center_x < cols and 0 <= center_y < rows):
        print(
            f"⚠️ Warning: center point ({center_x:.1f}, {center_y:.1f}) is outside "
            f"image bounds for slice {idx}: x=[0, {cols - 1}], y=[0, {rows - 1}]"
        )

    # Expected outputs for skip logic
    expected_paths = []
    if write_angles:
        a1, a2 = angle_names
        expected_paths = [
            os.path.join(output_dir, a1, f"{a1}_{idx:06d}.{ext}"),
            os.path.join(output_dir, a2, f"{a2}_{idx:06d}.{ext}"),
            os.path.join(output_dir, "FA", f"FA_{idx:06d}.{ext}"),
        ]
    if write_vectors and vector_format == "npy":
        expected_paths.append(
            str(vector_field_path(output_dir, "npy") / f"eigen_vec_{idx:06d}.npy")
        )

    # Skip if all outputs are already present and we are not in test mode
    if (
        not is_test
        and (write_angles or write_vectors)
        and all(os.path.exists(p) for p in expected_paths)
        and (not write_vectors or vector_already_written)
    ):
        return 0.0, 0.0, True

    compute_t0 = time.perf_counter()

    # Compute FA and the chosen angle pair
    if write_angles or is_test:
        if eigen_val_slice is None:
            raise ValueError("Eigenvalues are required for angle and FA calculation")
        center_vec = calculate_center_vector(_centerline_neighborhood(center_line, idx))
        img_FA = compute_fraction_anisotropy(eigen_val_slice)
        vector_field_slice_rotated = rotate_vectors_to_new_axis(
            vector_field_slice, center_vec
        )

        if mode == "ha_ia":
            img_angle1, img_angle2 = compute_helical_and_intrusion_angles(
                vector_field_slice_rotated, center_point, projected=projected
            )
        else:  # "az_el"
            img_angle1, img_angle2 = compute_azimuth_and_elevation(
                vector_field_slice_rotated
            )

        if invalid_mask_slice is not None:
            img_angle1[invalid_mask_slice] = np.nan
            img_angle2[invalid_mask_slice] = np.nan
            img_FA[invalid_mask_slice] = np.nan

    compute_time = time.perf_counter() - compute_t0
    write_t0 = time.perf_counter()

    # Test mode: visualize a 2x2 figure and write to test subfolder
    if is_test:
        if img_slice is None:
            raise ValueError("The intensity image is required in test mode")
        if mode == "ha_ia":
            titles = (
                ("Projected Helical Angle", "Projected Intrusion Angle")
                if projected
                else ("Helical Angle", "Intrusion Angle")
            )
        else:
            titles = ("Azimuth", "Elevation")
        test_output_dir = os.path.join(output_dir, "test_slice")
        plot_images(
            img_slice,
            img_angle1,
            img_angle2,
            img_FA,
            center_point,
            vector_field_slice=vector_field_slice if show_quiver else None,
            colormap_angle=colormap_angle,
            colormap_angle2=colormap_angle2,
            angle1_title=titles[0],
            angle2_title=titles[1],
            angle_ranges=angle_ranges,
            save_path=os.path.join(test_output_dir, "result_test_slice.png"),
            overlay_scalar_map=img_angle1 if show_quiver else None,
        )
        write_images(
            img_angle1,
            img_angle2,
            img_FA,
            idx,
            test_output_dir,
            ext,
            output_type,
            0,
            colormap_angle=colormap_angle,
            colormap_angle2=colormap_angle2,
            angle_names=angle_names,
            angle_ranges=angle_ranges,
        )
        return compute_time, time.perf_counter() - write_t0, False

    # Persist outputs
    if write_angles:
        write_images(
            img_angle1,
            img_angle2,
            img_FA,
            idx,
            output_dir,
            ext,
            output_type,
            0,
            colormap_angle=colormap_angle,
            colormap_angle2=colormap_angle2,
            angle_names=angle_names,
            angle_ranges=angle_ranges,
        )
    if write_vectors and not vector_already_written:
        write_vector_field(
            vector_field_slice,
            idx,
            output_dir,
            0,
            vector_format=vector_format,
            zarr_store=zarr_store,
        )

    return compute_time, time.perf_counter() - write_t0, False
