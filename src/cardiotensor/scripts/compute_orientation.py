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

    This function supports two modes of operation:
    1. Interactive mode: Opens a file dialog to select a configuration file if no arguments are provided.
    2. Command-line mode: Parses command-line arguments to specify the configuration file, start and end indices,
       and whether to use GPU for computations.

    Returns:
        None
    """
    # Set up the argument parser with a description
    parser = argparse.ArgumentParser(
        description=(
            "This script computes orientation for a 3D volume based on the provided configuration file. "
            "You can run it in interactive mode (no arguments) or provide the configuration file path "
            "and other options as command-line arguments. It supports GPU acceleration and chunk-based processing."
        )
    )

    # Add arguments to the parser
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
        default=None,
        help="Ending index for processing (default: None, processes all data).",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Activate GPU acceleration for computations.",
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
        conf_file_path = askopenfilename(
            initialdir=f"{os.getcwd()}/param_files", title="Select file"
        )
        if not conf_file_path:
            sys.exit("No file selected!")
        start_index = 0
        end_index = None
        use_gpu = True
        reverse = False

    else:
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
        start_index = args.start_index
        end_index = args.end_index
        use_gpu = args.gpu
        reverse = args.reverse

    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(f"⚠️  Error reading parameter file '{conf_file_path}': {e}")
        sys.exit(1)

    # Extracting parameters safely using .get() with defaults where necessary
    VOLUME_PATH = params.get("IMAGES_PATH", "")
    if reverse is False:
        reverse = params.get("REVERSE", False)
    N_CHUNK = params.get("N_CHUNK", 100)
    IS_TEST = params.get("TEST", False)

    data_reader = DataReader(VOLUME_PATH)
    volume_shape = data_reader.shape

    # Set end_index to total_images if it's zero
    if end_index is None:
        end_index = volume_shape[0]

    if not IS_TEST:
        # If reverse is True, run the loop in reverse order; otherwise, run normally
        if reverse:
            for idx in range(end_index, start_index, -N_CHUNK):
                start_time = time.time()
                compute_orientation(
                    conf_file_path,
                    start_index=max(idx - N_CHUNK, 0),
                    end_index=idx,
                    use_gpu=use_gpu,
                )
                print("--- %s seconds ---" % (time.time() - start_time))
        else:
            for idx in range(start_index, end_index, N_CHUNK):
                start_time = time.time()
                compute_orientation(
                    conf_file_path,
                    start_index=idx,
                    end_index=min(idx + N_CHUNK, end_index),
                    use_gpu=use_gpu,
                )
                print("--- %s seconds ---" % (time.time() - start_time))

    else:
        start_time = time.time()
        compute_orientation(conf_file_path, start_index, end_index, use_gpu=use_gpu)
        print("--- %s seconds ---" % (time.time() - start_time))

    return None
