"""
3D_Data_Processing
"""

import argparse
import os
import sys
from tkinter import Tk  # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename

from cardiotensor.launcher.slurm_launcher import slurm_launcher
from cardiotensor.utils.utils import read_conf_file


def script() -> None:
    """
    Main script to process 3D data. Reads configuration files, launches processing tasks,
    and logs processing time.
    """
    if len(sys.argv) < 2:
        # If no argument is passed, show file dialog to select a configuration file
        Tk().withdraw()
        conf_file_path: str = askopenfilename(
            initialdir=f"{os.getcwd()}/param_files", title="Select file"
        )  # Show an "Open" dialog box and return the path to the selected file
        if not conf_file_path:
            sys.exit("No file selected!")

    else:
        # Parse the configuration file path from command-line arguments
        parser = argparse.ArgumentParser(
            description="Process 3D data using the specified configuration file."
        )
        parser.add_argument(
            "conf_file_path", type=str, help="Path to the input configuration file."
        )
        args = parser.parse_args()
        conf_file_path = args.conf_file_path

    try:
        # Read configuration parameters from the selected file
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(e)
        sys.exit(f"⚠️  Error reading parameter file: {conf_file_path}")

    # Extract configuration parameters
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

    # Launch processing using slurm_launcher
    slurm_launcher(conf_file_path)


if __name__ == "__main__":
    script()
