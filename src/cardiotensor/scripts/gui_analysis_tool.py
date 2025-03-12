import argparse
import sys

from PyQt5.QtWidgets import (
    QApplication,
)

from cardiotensor.analysis.gui_analysis_tool import Window


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the image processing script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Open a GUI to interactively plot transmural profiles of the angles."
    )
    parser.add_argument(
        "conf_file_path", type=str, help="Path to the configuration file"
    )
    parser.add_argument("N_slice", type=int, help="Slice number")
    parser.add_argument("--N_line", type=int, default=5, help="Number of lines")
    parser.add_argument(
        "--angle_range", type=float, default=20, help="Angle range in degrees"
    )
    parser.add_argument(
        "--image_mode", type=str, default="HA", help="Output mode (HA, IA, or FA)"
    )

    return parser.parse_args()


def script() -> None:
    """
    Launch the GUI for analyzing image slices based on the provided configuration.

    This function initializes a PyQt5 application and opens the GUI window
    with the specified configuration parameters.

    Returns:
        None
    """
    args = parse_arguments()

    app = QApplication(sys.argv)
    w = Window(
        args.conf_file_path,
        args.N_slice,
        args.N_line,
        args.angle_range,
        args.image_mode,
    )
    w.show()
    app.exec()
