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

    else:
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
        start_index = args.start_index
        end_index = args.end_index
        use_gpu = args.gpu

    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(f"⚠️  Error reading parameter file: {conf_file_path}")
        print(f"\nError is {e}")
        sys.exit()

    (
        VOLUME_PATH,
        MASK_PATH,
        IS_FLIP,
        OUTPUT_DIR,
        OUTPUT_TYPE,
        SIGMA,
        RHO,
        N_CHUNK,
        PT_MV,
        PT_APEX,
        REVERSE,
        IS_TEST,
        N_SLICE_TEST,
    ) = (
        params[key]
        for key in [
            "IMAGES_PATH",
            "MASK_PATH",
            "FLIP",
            "OUTPUT_PATH",
            "OUTPUT_TYPE",
            "SIGMA",
            "RHO",
            "N_CHUNK",
            "POINT_MITRAL_VALVE",
            "POINT_APEX",
            "REVERSE",
            "TEST",
            "N_SLICE_TEST",
        ]
    )

    data_reader = DataReader(VOLUME_PATH)
    volume_shape = data_reader.shape

    # Set end_index to total_images if it's zero
    if end_index is None:
        end_index = volume_shape[0]

    if not IS_TEST:
        # If REVERSE is True, run the loop in reverse order; otherwise, run normally
        if REVERSE:
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
