import argparse
import sys

from PyQt5.QtWidgets import (
    QApplication,
)

from cardiotensor.analysis.gui_analysis_tool import Window

# from cardiotensor.orientation import (
#     interpolate_points,
#     calculate_center_vector,
# )


def parse_arguments():
    parser = argparse.ArgumentParser(description="Image Processing Script")
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


def script():
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
