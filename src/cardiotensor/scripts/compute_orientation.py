import argparse
import math
import sys
import time

from cardiotensor.orientation.orientation_computation_pipeline import (
    _safe_low_memory_chunk_size,
    compute_orientation,
)
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import read_conf_file


def script() -> None:
    """
    Main entry point for the orientation computation pipeline.

    Supports:
    1. Interactive mode (no args): GUI file picker for .conf.
    2. CLI mode: Use arguments for configuration and processing range.
    """
    parser = argparse.ArgumentParser(
        description="Compute 3D cardiomyocyte orientation using a configuration file."
    )

    parser.add_argument(
        "conf_file_path",
        type=str,
        nargs="?",
        help="Path to the configuration file (.conf). Optional in interactive mode.",
    )
    parser.add_argument(
        "--start_index",
        type=int,
        default=0,
        help="Starting slice index for processing (default: 0).",
    )
    parser.add_argument(
        "--end_index",
        type=int,
        default=None,
        help="Ending slice index (default: None, processes all slices).",
    )
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="Process slices in reverse order starting from the end.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Force single-slice test mode even if TEST is False in the config.",
    )
    parser.add_argument(
        "--colormap",
        type=str,
        default=None,
        help="Shared colormap for both angle outputs, for example 'hsv'.",
    )
    parser.add_argument(
        "--colormap-angle1",
        type=str,
        default=None,
        help="Colormap for the first angle output, for example HA or AZ.",
    )
    parser.add_argument(
        "--colormap-angle2",
        type=str,
        default=None,
        help="Colormap for the second angle output, for example IA or EL.",
    )
    parser.add_argument(
        "--projected",
        action="store_true",
        help="In ha_ia mode, write legacy projected HA/IA maps instead of 3D HA/IA.",
    )

    args = parser.parse_args()
    conf_file_path = args.conf_file_path
    start_index = args.start_index
    end_index = args.end_index
    reverse = args.reverse
    force_test = args.test
    colormap = args.colormap
    colormap_angle1 = args.colormap_angle1
    colormap_angle2 = args.colormap_angle2
    projected = args.projected

    # --- Read configuration ---
    try:
        params = read_conf_file(conf_file_path)
    except Exception as err:
        print(f"⚠️  Error reading configuration file '{conf_file_path}': {err}")
        sys.exit(1)

    # Extract parameters
    volume_path = params.get("IMAGES_PATH", "")
    if not volume_path:
        print("❌ No volume path specified in the configuration file.")
        sys.exit(1)
    mask_path = params.get("MASK_PATH", None)
    output_dir = params.get("OUTPUT_PATH", "./output")
    output_format = params.get("OUTPUT_FORMAT", "jp2")
    output_type = params.get("OUTPUT_TYPE", "8bit")
    sigma = params.get("SIGMA", 1.0)
    rho = params.get("RHO", 3.0)
    truncate = params.get("TRUNCATE", 4.0)
    axis_points = params.get("AXIS_POINTS", None)
    vertical_padding = params.get("VERTICAL_PADDING", None)
    write_vectors = params.get("WRITE_VECTORS", False)
    vector_format = params.get("VECTOR_FORMAT", "zarr")
    write_angles = params.get("WRITE_ANGLES", True)
    angle_mode = params.get("ANGLE_MODE", "ha_ia")
    use_gpu = params.get("USE_GPU", True)
    low_memory = params.get("LOW_MEMORY", False)
    low_memory_dir = params.get("LOW_MEMORY_DIR", None)
    is_test = force_test or params.get("TEST", False)
    n_slice_test = params.get("N_SLICE_TEST", None)
    show_quiver = params.get("SHOW_QUIVER", True)
    projected = projected or params.get("PROJECTED_ANGLES", False)
    n_chunk = params.get("N_CHUNK", 100)
    colormap = colormap or params.get("COLORMAP", None)
    colormap_angle1 = colormap_angle1 or params.get("COLORMAP_ANGLE1", None)
    colormap_angle2 = colormap_angle2 or params.get("COLORMAP_ANGLE2", None)

    if not reverse:  # CLI --reverse overrides config
        reverse = params.get("REVERSE", False)

    # --- Validate Volume ---
    data_reader = DataReader(volume_path)
    total_slices = data_reader.shape[0]
    if end_index is None:
        end_index = total_slices

    # --- Test Mode ---
    if is_test:
        if n_slice_test is None:
            print("❌ Test mode requires N_SLICE_TEST in the configuration file.")
            sys.exit(1)
        if force_test and not params.get("TEST", False):
            print("⚙️  TEST mode forced from CLI (--test)")
        print(f"⚙️  TEST mode: processing slice {n_slice_test}")
        t0 = time.time()
        compute_orientation(
            volume_path=volume_path,
            mask_path=mask_path,
            output_dir=output_dir,
            output_format=output_format,
            output_type=output_type,
            sigma=sigma,
            rho=rho,
            truncate=truncate,
            axis_points=axis_points,
            vertical_padding=vertical_padding,
            write_vectors=write_vectors,
            vector_format=vector_format,
            angle_mode=angle_mode,
            write_angles=write_angles,
            use_gpu=use_gpu,
            low_memory=low_memory,
            low_memory_dir=low_memory_dir,
            is_test=is_test,
            n_slice_test=n_slice_test,
            show_quiver=show_quiver,
            start_index=start_index,
            end_index=end_index,
            colormap=colormap,
            colormap_angle1=colormap_angle1,
            colormap_angle2=colormap_angle2,
            projected=projected,
        )
        print(f"--- {time.time() - t0:.1f} seconds (TEST mode) ---")
        return

    if low_memory:
        padding = math.ceil(
            vertical_padding
            if vertical_padding is not None
            else int(sigma * truncate + 0.5) + int(rho * truncate + 0.5)
        )
        requested_chunk = min(n_chunk, end_index - start_index)
        safe_chunk, requested_peak, safe_peak = _safe_low_memory_chunk_size(
            requested_chunk,
            total_slices,
            data_reader.shape[-2],
            data_reader.shape[-1],
            data_reader.dtype,
            padding=padding,
            include_eigenvalues=write_angles,
            has_mask=mask_path is not None,
            write_angles=write_angles,
            sigma=sigma,
            rho=rho,
            truncate=truncate,
        )
        if safe_chunk < requested_chunk:
            print(
                "⚠️  LOW_MEMORY adjusted N_CHUNK from "
                f"{requested_chunk} to {safe_chunk}: conservative peak "
                f"{requested_peak / 1024**3:.2f} GiB -> "
                f"{safe_peak / 1024**3:.2f} GiB."
            )
        else:
            print(
                f"LOW_MEMORY chunk check: N_CHUNK={safe_chunk}, "
                f"conservative peak={safe_peak / 1024**3:.2f} GiB."
            )
        n_chunk = safe_chunk

    # --- Build Chunks ---
    chunks: list[tuple[int, int]] = []
    if reverse:
        for e in range(end_index, start_index, -n_chunk):
            s = max(e - n_chunk, start_index)
            chunks.append((s, e))
    else:
        for s in range(start_index, end_index, n_chunk):
            e = min(s + n_chunk, end_index)
            chunks.append((s, e))

    print(f"📦 Will process {len(chunks)} chunks{' in reverse' if reverse else ''}.\n")

    # --- Process Chunks ---
    for idx, (s, e) in enumerate(chunks, start=1):
        print("=" * 60)
        print(f"▶️  Chunk {idx}/{len(chunks)}: slices {s}–{e - 1}")
        print("=" * 60)
        t0 = time.time()

        compute_orientation(
            volume_path=volume_path,
            mask_path=mask_path,
            output_dir=output_dir,
            output_format=output_format,
            output_type=output_type,
            sigma=sigma,
            rho=rho,
            truncate=truncate,
            axis_points=axis_points,
            vertical_padding=vertical_padding,
            write_vectors=write_vectors,
            vector_format=vector_format,
            write_angles=write_angles,
            angle_mode=angle_mode,
            use_gpu=use_gpu,
            low_memory=low_memory,
            low_memory_dir=low_memory_dir,
            is_test=is_test,
            n_slice_test=n_slice_test,
            show_quiver=show_quiver,
            start_index=s,
            end_index=e,
            colormap=colormap,
            colormap_angle1=colormap_angle1,
            colormap_angle2=colormap_angle2,
            projected=projected,
        )

        elapsed = time.time() - t0
        print(f"✅ Finished chunk {idx}/{len(chunks)} in {elapsed:.1f}s\n")

    print("🎉 All chunks completed successfully!")


if __name__ == "__main__":
    script()
