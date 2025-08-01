import argparse
import sys

from PyQt5.QtWidgets import QApplication

from cardiotensor.analysis.gui_analysis_tool import Window
from cardiotensor.utils.utils import read_conf_file


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the GUI analysis tool.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Open a GUI to interactively plot transmural profiles of the angles."
    )
    parser.add_argument(
        "conf_file_path", type=str, help="Path to the configuration file"
    )
    parser.add_argument(
        "--N_slice", type=int, default=None, help="Slice number (optional)"
    )
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
    """
    args = parse_arguments()

    # Load parameters from configuration file
    try:
        params = read_conf_file(args.conf_file_path)
    except Exception as e:
        print(f"⚠️ Error reading configuration file '{args.conf_file_path}': {e}")
        sys.exit(1)

    # Extract parameters safely
    mask_path = params.get("MASK_PATH", "")
    output_dir = params.get("OUTPUT_PATH", "./output")

    # Determine slice number: CLI overrides config
    N_slice = (
        args.N_slice if args.N_slice is not None else params.get("N_SLICE_TEST", 0)
    )

    # Initialize the PyQt5 application
    app = QApplication(sys.argv)
    w = Window(
        output_dir,
        mask_path,
        N_slice=N_slice,
        N_line=args.N_line,
        angle_range=args.angle_range,
        image_mode=args.image_mode,
    )
    w.show()
    app.exec()


if __name__ == "__main__":
    script()
