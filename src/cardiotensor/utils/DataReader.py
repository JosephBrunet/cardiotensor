from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
import tifffile
import zarr
from alive_progress import alive_bar
from skimage.measure import block_reduce

from cardiotensor.utils.image_io import _validate_zarr_vector_group, read_image_file
from cardiotensor.utils.utils import (
    get_available_cpu_count,
    get_available_memory_bytes,
)


_MHD_DTYPES = {
    "MET_FLOAT": np.float32,
    "MET_DOUBLE": np.float64,
    "MET_CHAR": np.int8,
    "MET_UCHAR": np.uint8,
    "MET_SHORT": np.int16,
    "MET_USHORT": np.uint16,
    "MET_INT": np.int32,
    "MET_UINT": np.uint32,
}


# ---------------------------
# Integer-only upsampling util
# ---------------------------
def _upsample_mask_integer(mask: np.ndarray, zf: int, yf: int, xf: int) -> np.ndarray:
    """Pure-NumPy nearest-neighbor upsampling for integer factors."""
    out = mask
    if zf != 1:
        out = np.repeat(out, zf, axis=0)
    if yf != 1:
        out = np.repeat(out, yf, axis=1)
    if xf != 1:
        out = np.repeat(out, xf, axis=2)
    return out


def _fit(
    arr: np.ndarray, target: tuple[int, int, int], pad_value: int = 0
) -> np.ndarray:
    """
    Force `arr` (Z, Y, X) to exactly match `target` shape.
    Crops if too large, pads with `pad_value` if too small.
    """
    tz, ty, tx = target
    z, y, x = arr.shape

    if (z, y, x) == target:
        return arr

    if z >= tz and y >= ty and x >= tx:
        return arr[:tz, :ty, :tx]

    fitted = np.full(target, pad_value, dtype=arr.dtype)
    fitted[: min(z, tz), : min(y, ty), : min(x, tx)] = arr[:tz, :ty, :tx]
    return fitted


def _normalize_single_tiff_volume(volume: np.ndarray) -> np.ndarray:
    """
    Normalize a single TIFF file to a stack-like layout (Z, Y, X).

    Supported inputs:
    - 2D grayscale image -> (1, Y, X)
    - 3D grayscale stack -> (Z, Y, X)
    """
    if volume.ndim == 2:
        return volume[None, :, :]

    if volume.ndim == 3:
        if volume.shape[2] == 1:
            return volume[:, :, 0][None, :, :]
        if volume.shape[2] in {3, 4}:
            raise ValueError(
                "Single TIFF file appears to be a color image, not a grayscale stack."
            )
        return volume

    raise ValueError(
        f"Unsupported TIFF volume shape {volume.shape}. Expected 2D image or 3D stack."
    )


class DataReader:
    def __init__(self, path: str | Path):
        """
        Initializes the DataReader with a path to the volume.

        Args:
            path (str | Path): Path to the volume directory or file.
        """
        self.path = Path(path)
        self.supported_extensions = ["tif", "tiff", "jp2", "png", "npy"]
        self.volume_info = self._get_volume_info()

    # ---------------------------
    # Properties
    # ---------------------------
    @property
    def shape(self) -> tuple[int, ...]:
        """Returns the shape of the volume as (Z, Y, X) or (Z, Y, X, C)."""
        return self.volume_info["shape"]

    @property
    def dtype(self) -> np.dtype:
        """Returns the data type of the volume."""
        return self.volume_info["dtype"]

    @property
    def volume_size_gb(self) -> float:
        """Returns the total size of the volume in GB."""
        n_bytes = np.prod(self.shape) * np.dtype(self.dtype).itemsize
        return n_bytes / (1024**3)

    def _get_volume_info(self) -> dict:
        """
        Detects volume type, shape, and dtype.
        Returns a dict with keys: type, stack, file_list, shape, dtype
        """
        volume_info = {
            "type": "",
            "stack": False,
            "file_list": [],
            "shape": None,
            "dtype": None,
        }

        if not self.path.exists():
            raise ValueError(f"The path does not exist: {self.path}")

        # Case 1: Cardiotensor Zarr vector field
        if self.path.is_dir() and (self.path / "zarr.json").is_file():
            group = zarr.open_group(store=str(self.path), mode="r")
            _validate_zarr_vector_group(group)
            vectors = group["vectors"]
            volume_info["type"] = "zarr"
            volume_info["shape"] = tuple(vectors.shape)
            volume_info["dtype"] = np.dtype(vectors.dtype)

        # Case 2: Directory of images
        elif self.path.is_dir():
            volume_info["stack"] = True
            image_files = {
                ext: sorted(self.path.glob(f"*.{ext}"))
                for ext in self.supported_extensions
            }
            volume_info["type"], volume_info["file_list"] = max(
                image_files.items(), key=lambda item: len(item[1])
            )
            if not volume_info["file_list"]:
                raise ValueError(
                    "No supported image files found in the specified directory."
                )

            # Inspect first file
            first_image = self._custom_image_reader(volume_info["file_list"][0])
            volume_info["dtype"] = first_image.dtype

            # Shape: handle scalar vs vector
            if (
                volume_info["type"] == "npy"
                and first_image.ndim == 3
                and first_image.shape[0] == 3
            ):
                # Vector field stored as (3, Y, X)
                volume_info["shape"] = (
                    3,
                    len(volume_info["file_list"]),
                    first_image.shape[1],
                    first_image.shape[2],
                )
            elif first_image.ndim == 3:
                # 4D scalar stack (Z from files)
                volume_info["shape"] = (
                    len(volume_info["file_list"]),
                    *first_image.shape,
                )
            else:
                # Standard 3D stack
                volume_info["shape"] = (
                    len(volume_info["file_list"]),
                    first_image.shape[0],
                    first_image.shape[1],
                )

        # Case 3: Single MHD file
        elif self.path.is_file() and self.path.suffix.lower() == ".mhd":
            volume_info["type"] = "mhd"
            metadata = _read_mhd(self.path)
            if int(metadata.get("NDims", 0)) != 3:
                raise ValueError("Only 3D MHD volumes are supported")
            try:
                dtype = np.dtype(_MHD_DTYPES[metadata["ElementType"]])
            except KeyError as exc:
                raise NotImplementedError(
                    f"ElementType {metadata.get('ElementType')} is not supported"
                ) from exc
            byte_order = ">" if metadata.get("BinaryDataByteOrderMSB", False) else "<"
            dtype = dtype.newbyteorder(byte_order)
            shape = tuple(reversed(metadata["DimSize"]))
            channels = int(metadata.get("ElementNumberOfChannels", 1))
            volume_info["shape"] = shape if channels == 1 else (*shape, channels)
            volume_info["dtype"] = dtype

        # Case 4: Single TIFF file storing a stack
        elif self.path.is_file() and self.path.suffix.lower() in {".tif", ".tiff"}:
            volume_info["type"] = "tiff"
            with tifffile.TiffFile(self.path) as tif:
                series = tif.series[0]
                shape = tuple(series.shape)
                dtype = np.dtype(series.dtype)
            if len(shape) == 2:
                shape = (1, *shape)
            elif len(shape) != 3 or shape[-1] in {3, 4}:
                raise ValueError(
                    f"Unsupported TIFF volume shape {shape}. "
                    "Expected a 2D image or 3D grayscale stack."
                )
            volume_info["shape"] = shape
            volume_info["dtype"] = dtype

        else:
            raise ValueError(f"Unsupported volume type for path: {self.path}")

        return volume_info

    def load_volume(
        self,
        start_index: int = 0,
        end_index: int | None = None,
        unbinned_shape: tuple[int, int, int] | None = None,
        show_progress: bool = True,
    ) -> np.ndarray:
        """
        Loads the volume and resizes it to unbinned_shape if provided, using fast
        integer-only resampling:
        - np.repeat for upsampling
        - block_reduce (max) for downsampling

        Args:
            start_index (int): Start index for slicing (for stacks).
            end_index (int): End index for slicing (for stacks). If None, loads the entire stack.
            unbinned_shape (tuple): Desired shape (Z, Y, X). If None, no resizing is done.
            show_progress: Display the stack-loading progress bar.

        Returns:
            np.ndarray: Loaded volume.
        """
        if end_index is None:
            end_index = self.shape[1] if len(self.shape) == 4 else self.shape[0]

        # Check memory available is enough
        effective_shape = list(self.shape)
        if len(effective_shape) == 3:
            effective_shape[0] = end_index - start_index
        elif len(effective_shape) == 4:
            effective_shape[1] = end_index - start_index
        self.check_memory_requirement(
            tuple(effective_shape), self.dtype, verbose=show_progress
        )

        # Decide if resize is needed
        need_resize = False
        if unbinned_shape is not None and self.shape != unbinned_shape:
            need_resize = True
            zoom_factors = tuple(u / s for u, s in zip(unbinned_shape, self.shape))
            if show_progress:
                print(f"Resample factors: {zoom_factors}")
        else:
            zoom_factors = (1.0, 1.0, 1.0)

        if need_resize:
            start_index_ini, end_index_ini = start_index, end_index
            start_index = int(start_index_ini / zoom_factors[0]) - 1
            start_index = max(start_index, 0)
            end_index = int(end_index_ini / zoom_factors[0]) + 1
            end_index = min(end_index, self.shape[0])
            if show_progress:
                print(f"Volume start index padded: {start_index} - end: {end_index}")

        # Load volume from stack or mhd
        if not self.volume_info["stack"]:
            if self.volume_info["type"] == "zarr":
                group = zarr.open_group(store=str(self.path), mode="r")
                volume = np.asarray(group["vectors"][:, start_index:end_index])
            elif self.volume_info["type"] == "mhd":
                volume, _ = _load_raw_data_with_mhd(
                    self.path, start_index=start_index, end_index=end_index
                )
            elif self.volume_info["type"] == "tiff":
                with tifffile.TiffFile(self.path) as tif:
                    series = tif.series[0]
                    if len(series.pages) == self.shape[0]:
                        volume = np.array(
                            _normalize_single_tiff_volume(
                                tif.asarray(
                                    key=slice(start_index, end_index), series=series
                                )
                            ),
                            copy=True,
                            order="C",
                        )
                    else:
                        scratch_root = self.path.parent / ".cardiotensor_scratch"
                        scratch_root.mkdir(parents=True, exist_ok=True)
                        with TemporaryDirectory(
                            prefix="cardiotensor_tiff_", dir=scratch_root
                        ) as tmpdir:
                            mapped = series.asarray(
                                out=str(Path(tmpdir) / "volume.memmap")
                            )
                            try:
                                volume = np.array(
                                    _normalize_single_tiff_volume(
                                        mapped[start_index:end_index]
                                    ),
                                    copy=True,
                                    order="C",
                                )
                            finally:
                                mapped_file = getattr(mapped, "_mmap", None)
                                if mapped_file is not None:
                                    mapped_file.close()
            else:
                raise ValueError(f"Unsupported volume type for path: {self.path}")
        else:
            volume = self._load_image_stack(
                self.volume_info["file_list"],
                start_index,
                end_index,
                show_progress=show_progress,
            )

        if need_resize:
            if show_progress:
                print("Resizing with integer-only resampling...")

            _, y1, x1 = unbinned_shape
            z1 = end_index_ini - start_index_ini
            fz, fy, fx = zoom_factors

            def _as_int_factor_up(f: float) -> int | None:
                k = int(round(f))
                return k if (k >= 1 and abs(f - k) < 1e-1) else None

            def _as_int_factor_down(f: float) -> int | None:
                # f < 1 expected, so divisor ~ round(1/f)
                if f >= 1:
                    return None
                k = int(round(1.0 / f))
                return k if k >= 1 and abs(f - (1.0 / k)) < 1e-1 else None

            # init all factors to None
            kz = ky = kx = None
            dz = dy = dx = None

            # Z
            if fz > 1:
                kz = _as_int_factor_up(fz)
                if kz is None:
                    raise ValueError(f"Non-integer upsample factor on Z: {fz}")
                volume = np.repeat(volume, kz, axis=0)
            elif fz < 1:
                dz = _as_int_factor_down(fz)
                if dz is None:
                    raise ValueError(f"Non-integer downsample factor on Z: {fz}")
                volume = block_reduce(volume, block_size=(dz, 1, 1), func=np.max)

            # Y
            if fy > 1:
                ky = _as_int_factor_up(fy)
                if ky is None:
                    raise ValueError(f"Non-integer upsample factor on Y: {fy}")
                volume = np.repeat(volume, ky, axis=1)
            elif fy < 1:
                dy = _as_int_factor_down(fy)
                if dy is None:
                    raise ValueError(f"Non-integer downsample factor on Y: {fy}")
                volume = block_reduce(volume, block_size=(1, dy, 1), func=np.max)

            # X
            if fx > 1:
                kx = _as_int_factor_up(fx)
                if kx is None:
                    raise ValueError(f"Non-integer upsample factor on X: {fx}")
                volume = np.repeat(volume, kx, axis=2)
            elif fx < 1:
                dx = _as_int_factor_down(fx)
                if dx is None:
                    raise ValueError(f"Non-integer downsample factor on X: {fx}")
                volume = block_reduce(volume, block_size=(1, 1, dx), func=np.max)

            def _fmt(k, d):
                if k is not None:
                    return f"x{k}"
                if d is not None:
                    return f"/{d}"
                return "x1"

            if show_progress:
                print(
                    f"Applied integer resampling: Z {_fmt(kz, dz)}, "
                    f"Y {_fmt(ky, dy)}, X {_fmt(kx, dx)} -> "
                    f"resulting shape {volume.shape}"
                )

            # Enforce exact shape
            if show_progress:
                print(f"Fitting to exact unbinned shape: {unbinned_shape}")
            volume = _fit(volume, (z1, y1, x1), pad_value=0)

        return volume

    def load_region(
        self,
        start_index: int = 0,
        end_index: int | None = None,
        start_y: int = 0,
        end_y: int | None = None,
        start_x: int = 0,
        end_x: int | None = None,
    ) -> np.ndarray:
        """Load a Z/Y/X region, using direct chunk reads for Zarr vectors."""
        is_vector = len(self.shape) == 4 and self.shape[0] == 3
        depth = self.shape[1] if is_vector else self.shape[0]
        height = self.shape[-2]
        width = self.shape[-1]
        end_index = depth if end_index is None else end_index

        if not 0 <= start_index < end_index <= depth:
            raise ValueError(
                f"Invalid Z range [{start_index}, {end_index}) for depth {depth}"
            )

        y_slice = slice(start_y, end_y)
        x_slice = slice(start_x, end_x)
        y_start, y_stop, y_step = y_slice.indices(height)
        x_start, x_stop, x_step = x_slice.indices(width)
        if y_step != 1 or x_step != 1 or y_start >= y_stop or x_start >= x_stop:
            raise ValueError("Y and X regions must be non-empty contiguous ranges")

        if self.volume_info["type"] == "zarr":
            region_shape = (
                3,
                end_index - start_index,
                y_stop - y_start,
                x_stop - x_start,
            )
            self.check_memory_requirement(region_shape, self.dtype)
            group = zarr.open_group(store=str(self.path), mode="r")
            return np.asarray(
                group["vectors"][
                    :, start_index:end_index, y_start:y_stop, x_start:x_stop
                ]
            )

        volume = self.load_volume(start_index=start_index, end_index=end_index)
        if is_vector:
            return volume[:, :, y_start:y_stop, x_start:x_stop]
        return volume[:, y_start:y_stop, x_start:x_stop]

    def _custom_image_reader(self, file_path: Path) -> np.ndarray:
        """
        Reads an image from the given file path into a NumPy array.
        """
        return read_image_file(file_path)

    def _load_image_stack(
        self,
        file_list: list[Path],
        start_index: int,
        end_index: int,
        show_progress: bool = True,
    ) -> np.ndarray:
        """
        Efficiently loads a stack of images into a 3D (or 4D for vector .npy) NumPy array
        by preallocating the output and filling it in place. The number of
        in-flight reads is bounded so loaded slices are not duplicated in memory.
        """

        if end_index == 0:
            end_index = len(file_list)

        if start_index < 0 or end_index > len(file_list) or start_index >= end_index:
            raise ValueError(
                f"Invalid indices: start_index={start_index}, end_index={end_index}, total_files={len(file_list)}"
            )

        paths = sorted(file_list[start_index:end_index])
        total_files = len(paths)

        # Probe first slice to determine per slice shape and dtype
        first = self._custom_image_reader(paths[0])
        slice_shape = first.shape
        slice_dtype = first.dtype

        # Preallocate output with correct layout
        if (
            self.volume_info["type"] == "npy"
            and first.ndim == 3
            and first.shape[0] == 3
        ):
            # Vector field stored per file as (3, Y, X) -> target (3, Z, Y, X)
            out_shape = (3, total_files, slice_shape[1], slice_shape[2])
            volume = np.empty(out_shape, dtype=slice_dtype)

            # Place first slice
            volume[:, 0, :, :] = first
            start_fill_idx = 1

            def _assign(z_idx: int, arr: np.ndarray):
                if arr.shape != slice_shape:
                    raise ValueError(
                        f"Inconsistent file shape at z {z_idx}: expected {slice_shape}, got {arr.shape}"
                    )
                volume[:, z_idx, :, :] = arr

        else:
            # Scalar per file, typical 2D image -> target (Z, Y, X)
            if first.ndim == 3:
                # handle RGB or multi channel by squeezing single channel if needed
                # if you want to allow multi channel volumes, adapt here
                # for now, enforce 2D
                if first.shape[2] == 1:
                    first = first[:, :, 0]
                    slice_shape = first.shape
                else:
                    raise ValueError(
                        f"Expected 2D slice, got 3D with shape {first.shape}. Adjust reader if you need multi channel."
                    )
            out_shape = (total_files, slice_shape[0], slice_shape[1])
            volume = np.empty(out_shape, dtype=first.dtype)

            # Place first slice
            volume[0, :, :] = first
            start_fill_idx = 1

            def _assign(z_idx: int, arr: np.ndarray):
                if arr.ndim == 3 and arr.shape[2] == 1:
                    arr = arr[:, :, 0]
                if arr.shape != slice_shape:
                    raise ValueError(
                        f"Inconsistent file shape at z {z_idx}: expected {slice_shape}, got {arr.shape}"
                    )
                volume[z_idx, :, :] = arr

        # Keep concurrent decoded slices within a small, predictable buffer.
        # A single very large slice still gets one reader.
        read_buffer_bytes = 2 * 1024**3
        memory_workers = max(1, read_buffer_bytes // max(first.nbytes, 1))
        max_workers = min(
            32,
            get_available_cpu_count(default=8),
            memory_workers,
        )

        # Fill the rest with a progress bar
        with alive_bar(
            total_files,
            title="Loading Volume",
            length=40,
            disable=not show_progress,
        ) as bar:
            # We already placed index 0
            bar()  # account for the first already loaded slice

            if start_fill_idx < total_files:
                indexed_paths = iter(enumerate(paths[start_fill_idx:], start_fill_idx))
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    pending = {}
                    for _ in range(max_workers):
                        try:
                            z_idx, path = next(indexed_paths)
                        except StopIteration:
                            break
                        pending[executor.submit(read_image_file, path)] = z_idx

                    while pending:
                        done, _ = wait(pending, return_when=FIRST_COMPLETED)
                        for future in done:
                            z_idx = pending.pop(future)
                            _assign(z_idx, future.result())
                            bar()
                            try:
                                next_idx, next_path = next(indexed_paths)
                            except StopIteration:
                                continue
                            pending[executor.submit(read_image_file, next_path)] = (
                                next_idx
                            )

        return volume

    def check_memory_requirement(self, shape, dtype, safety_factor=0.8, verbose=True):
        """
        Check if the dataset can fit in available memory.

        Args:
            shape (tuple[int]): Shape of the array.
            dtype (np.dtype): NumPy dtype of the array.
            safety_factor (float): Fraction of available memory allowed to be used.
        """
        # Compute dataset size in bytes
        n_bytes = np.prod(shape) * np.dtype(dtype).itemsize
        size_gb = n_bytes / (1024**3)

        # Check available memory
        available_gb = get_available_memory_bytes() / (1024**3)

        if verbose:
            print(
                f"Dataset size: {size_gb:.2f} GB | "
                f"Available memory: {available_gb:.2f} GB"
            )

        if size_gb > available_gb * safety_factor:
            raise MemoryError(
                f"Dataset requires {size_gb:.2f} GB but only "
                f"{available_gb:.2f} GB is available"
            )


def _read_mhd(filename: PathLike[str]) -> dict[str, Any]:
    """
    Return a dictionary of meta data from an MHD meta header file.

    Args:
        filename (PathLike[str]): File path to the .mhd file.

    Returns:
        dict[str, Any]: A dictionary containing parsed metadata.
    """
    meta_dict: dict[str, Any] = {}
    tag_set = [
        "ObjectType",
        "NDims",
        "DimSize",
        "ElementType",
        "ElementDataFile",
        "ElementNumberOfChannels",
        "BinaryData",
        "BinaryDataByteOrderMSB",
        "CompressedData",
        "CompressedDataSize",
        "Offset",
        "CenterOfRotation",
        "AnatomicalOrientation",
        "ElementSpacing",
        "TransformMatrix",
        "Comment",
        "SeriesDescription",
        "AcquisitionDate",
        "AcquisitionTime",
        "StudyDate",
        "StudyTime",
    ]

    with open(filename) as fn:
        for line in fn:
            tags = line.split("=")
            if len(tags) < 2:
                continue
            key, content = tags[0].strip(), tags[1].strip()
            if key in tag_set:
                if key in [
                    "ElementSpacing",
                    "Offset",
                    "CenterOfRotation",
                    "TransformMatrix",
                ]:
                    # Parse as a list of floats
                    meta_dict[key] = [float(value) for value in content.split()]
                elif key in ["NDims", "ElementNumberOfChannels"]:
                    # Parse as an integer
                    meta_dict[key] = int(content)
                elif key == "DimSize":
                    # Parse as a list of integers
                    meta_dict[key] = [int(value) for value in content.split()]
                elif key in ["BinaryData", "BinaryDataByteOrderMSB", "CompressedData"]:
                    # Parse as a boolean
                    meta_dict[key] = content.lower() == "true"
                else:
                    # Parse as a string
                    meta_dict[key] = content
    return meta_dict


def _load_raw_data_with_mhd(
    filename: PathLike[str],
    start_index: int = 0,
    end_index: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """
    Load a MHD file

    :param filename: file of type .mhd that should be loaded
    :returns: tuple with raw data and dictionary of meta data
    """
    meta_dict = _read_mhd(filename)
    if int(meta_dict.get("NDims", 0)) != 3:
        raise ValueError("Only 3D MHD volumes are supported")
    if meta_dict.get("CompressedData", False):
        raise NotImplementedError("Compressed MHD data cannot be memory-mapped")
    if not meta_dict.get("BinaryData", True):
        raise NotImplementedError("ASCII MHD data is not supported")

    try:
        np_type = np.dtype(_MHD_DTYPES[meta_dict["ElementType"]])
    except KeyError as exc:
        raise NotImplementedError(
            f"ElementType {meta_dict.get('ElementType')} is not supported"
        ) from exc
    byte_order = ">" if meta_dict.get("BinaryDataByteOrderMSB", False) else "<"
    np_type = np_type.newbyteorder(byte_order)
    shape = tuple(reversed(meta_dict["DimSize"]))
    element_channels = int(meta_dict.get("ElementNumberOfChannels", 1))
    if element_channels > 1:
        shape = (*shape, element_channels)

    pwd = Path(filename).parents[0].resolve()
    data_file = Path(meta_dict["ElementDataFile"])
    if str(data_file).upper() == "LOCAL":
        raise NotImplementedError("MHD files with embedded LOCAL data are unsupported")
    if not data_file.is_absolute():
        data_file = pwd / data_file

    depth = shape[0]
    end_index = depth if end_index is None else end_index
    if not 0 <= start_index < end_index <= depth:
        raise ValueError(
            f"Invalid MHD Z range [{start_index}, {end_index}) for depth {depth}"
        )

    # Copy-on-write keeps reads lazy and lets masking modify pages in memory
    # without changing the source RAW file.
    data = np.memmap(data_file, dtype=np_type, mode="c", shape=shape)
    return data[start_index:end_index], meta_dict
