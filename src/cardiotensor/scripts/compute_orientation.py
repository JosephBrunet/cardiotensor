import argparse
import sys
import time

from cardiotensor.orientation.orientation_computation_pipeline import (
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
        help="Colormap to apply, for example 'hsv'. If None, use helix_angle_cmap.",
    )

    args = parser.parse_args()
    conf_file_path = args.conf_file_path
    start_index = args.start_index
    end_index = args.end_index
    reverse = args.reverse
    force_test = args.test
    colormap = args.colormap

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
    write_angles = params.get("WRITE_ANGLES", True)
    angle_mode = params.get("ANGLE_MODE", "ha_ia")
    use_gpu = params.get("USE_GPU", True)
    is_test = force_test or params.get("TEST", False)
    n_slice_test = params.get("N_SLICE_TEST", None)
    show_quiver = params.get("SHOW_QUIVER", True)
    n_chunk = params.get("N_CHUNK", 100)

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
            print(
                "❌ Test mode requires N_SLICE_TEST in the configuration file."
            )
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
            angle_mode=angle_mode,
            write_angles=write_angles,
            use_gpu=use_gpu,
            is_test=is_test,
            n_slice_test=n_slice_test,
            show_quiver=show_quiver,
            start_index=start_index,
            end_index=end_index,
            colormap=colormap,
        )
        print(f"--- {time.time() - t0:.1f} seconds (TEST mode) ---")
        return

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
            write_angles=write_angles,
            angle_mode=angle_mode,
            use_gpu=use_gpu,
            is_test=is_test,
            n_slice_test=n_slice_test,
            show_quiver=show_quiver,
            start_index=s,
            end_index=e,
            colormap=colormap,
        )

        elapsed = time.time() - t0
        print(f"✅ Finished chunk {idx}/{len(chunks)} in {elapsed:.1f}s\n")

    print("🎉 All chunks completed successfully!")


if __name__ == "__main__":
    script()
