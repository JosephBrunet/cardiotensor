import argparse
import os
import sys
import time
from tkinter import Tk
from tkinter.filedialog import askopenfilename

from cardiotensor.orientation.orientation_computation_pipeline import compute_orientation
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

    # --- Interactive Mode ---
    if len(sys.argv) < 2:
        parser.print_help()
        print("\nNo configuration file provided. Launching interactive mode...")
        Tk().withdraw()
        conf_file_path: str = askopenfilename(
            initialdir=f"{os.getcwd()}/param_files",
            title="Select a configuration file",
            filetypes=[("Config files", "*.conf *.json *.yaml *.yml"), ("All files", "*.*")]
        )
        if not conf_file_path:
            sys.exit("❌ No file selected! Exiting.")
        start_index: int = 0
        end_index: int | None = None
        reverse: bool = False
    else:
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
        start_index = args.start_index
        end_index = args.end_index
        reverse = args.reverse

    # --- Read configuration ---
    try:
        params = read_conf_file(conf_file_path)
    except Exception as err:
        print(f"⚠️  Error reading configuration file '{conf_file_path}': {err}")
        sys.exit(1)

    # Extract parameters
    volume_path = params.get("IMAGES_PATH", "")
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
    use_gpu = params.get("USE_GPU", True)
    is_test = params.get("TEST", False)
    n_slice_test = params.get("N_SLICE_TEST", None)
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
        print(f"⚙️  TEST mode: processing slices {start_index}–{end_index - 1}")
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
            use_gpu=use_gpu,
            is_test=is_test,
            n_slice_test=n_slice_test,
            start_index=start_index,
            end_index=end_index,
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
            use_gpu=use_gpu,
            is_test=is_test,
            n_slice_test=n_slice_test,
            start_index=s,
            end_index=e,
        )

        elapsed = time.time() - t0
        print(f"✅ Finished chunk {idx}/{len(chunks)} in {elapsed:.1f}s\n")

    print("🎉 All chunks completed successfully!")


if __name__ == "__main__":
    script()
