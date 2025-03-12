"""
3D_Data_Processing
"""

import argparse
import os
import sys
from tkinter import Tk  # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename

from cardiotensor.launcher.slurm_launcher import slurm_launcher


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

    # Launch processing using slurm_launcher
    slurm_launcher(conf_file_path)


if __name__ == "__main__":
    script()
