import argparse
import os
import sys
import time
from tkinter import Tk
from tkinter.filedialog import askopenfilename

from cardiotensor.orientation.orientation_computation_pipeline import (
    compute_orientation,
)
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import read_conf_file


def script() -> None:
    """
    Executes the main pipeline for computing orientation using configuration parameters.

    Supports:
    1. Interactive mode: Opens a file dialog if no CLI args provided.
    2. CLI mode: Parses arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "This script computes orientation for a 3D volume based on the provided configuration file. "
            "You can run it in interactive mode (no arguments) or provide the configuration file path "
            "and other options as command-line arguments. It supports GPU acceleration and chunk-based processing."
        )
    )

    # Explicit Optional for end_index
    parser.add_argument(
        "conf_file_path",
        type=str,
        nargs="?",
        help="Path to the configuration file (optional in interactive mode).",
    )
    parser.add_argument(
        "--start_index",
        type=int,
        default=0,
        help="Starting index for processing (default: 0).",
    )
    parser.add_argument(
        "--end_index",
        type=int,
        default=None,  # Must be Optional[int]
        help="Ending index for processing (default: None, processes all data).",
    )
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="Start the orientation computation from the end of the volume",
    )

    # Print help and exit if no arguments are provided
    if len(sys.argv) < 2:
        parser.print_help()
        Tk().withdraw()
        conf_file_path: str = askopenfilename(
            initialdir=f"{os.getcwd()}/param_files", title="Select file"
        )
        if not conf_file_path:
            sys.exit("No file selected!")
        start_index: int = 0
        end_index: int | None = None
        reverse: bool = False
    else:
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
        start_index = args.start_index
        end_index = args.end_index  # Optional[int]
        reverse = args.reverse

    try:
        params = read_conf_file(conf_file_path)
    except Exception as err:
        print(f"⚠️  Error reading parameter file '{conf_file_path}': {err}")
        sys.exit(1)

    VOLUME_PATH = params.get("IMAGES_PATH", "")
    if reverse is False:
        reverse = params.get("REVERSE", False)
    N_CHUNK = params.get("N_CHUNK", 100)
    IS_TEST = params.get("TEST", False)

    data_reader = DataReader(VOLUME_PATH)
    total_slices = data_reader.shape[0]
    if end_index is None:
        end_index = total_slices

    # TEST mode
    if IS_TEST:
        print(f"⚙️  TEST mode: processing slices {start_index}–{end_index - 1}")
        t0 = time.time()
        compute_orientation(
            conf_file_path,
            start_index=start_index,
            end_index=end_index,
        )
        print(f"--- {time.time() - t0:.1f} seconds (TEST mode) ---")
        return

    # Build chunks list
    chunks = []
    if reverse:
        for e in range(end_index, start_index, -N_CHUNK):
            s = max(e - N_CHUNK, start_index)
            chunks.append((s, e))
    else:
        for s in range(start_index, end_index, N_CHUNK):
            e = min(s + N_CHUNK, end_index)
            chunks.append((s, e))

    print(f"Will process {len(chunks)} chunks{' in reverse' if reverse else ''}.")

    # Execute chunks
    for idx, (s, e) in enumerate(chunks, start=1):
        print("=" * 60)
        print(f"▶️  Chunk {idx}/{len(chunks)}: slices {s}–{e - 1}")
        print("=" * 60)
        t0 = time.time()

        compute_orientation(
            conf_file_path,
            start_index=s,
            end_index=e,
        )

        elapsed = time.time() - t0
        print(f"✅ Finished chunk {idx}/{len(chunks)} in {elapsed:.1f}s\n")

    print("🤖 All chunks complete!")

    return None
