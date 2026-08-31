#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

from PyQt5.QtWidgets import QApplication

from cardiotensor.analysis.gui_analysis_tool import Window
from cardiotensor.utils.utils import read_conf_file


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for the GUI analysis tool.
    """
    parser = argparse.ArgumentParser(
        description="Open a GUI to interactively plot transmural profiles of FA or angle maps."
    )
    parser.add_argument(
        "conf_file_path",
        type=Path,
        help="Path to the configuration .conf file",
    )
    parser.add_argument(
        "--N-slice",
        type=int,
        default=None,
        help="Slice index to open, overrides config if provided",
    )
    parser.add_argument(
        "--N-line",
        type=int,
        default=5,
        help="Number of sampled lines around the main profile",
    )
    parser.add_argument(
        "--angle-range",
        type=float,
        default=20.0,
        help="Angular spread in degrees for the multi line sampling",
    )
    parser.add_argument(
        "--image-mode",
        type=str,
        default=None,
        choices=["HA", "IA", "EL", "AZ", "FA"],
        help=(
            "Which map to display first. Defaults to HA for ANGLE_MODE=ha_ia "
            "and AZ for ANGLE_MODE=az_el. Can be switched in the GUI"
        ),
    )
    return parser.parse_args()


def image_mode_from_config(params: dict) -> str:
    """Return the primary display map for the configured angle pair."""
    angle_mode = str(params.get("ANGLE_MODE", "ha_ia")).strip().lower()
    if angle_mode == "az_el":
        return "AZ"
    return "HA"


def script() -> None:
    """
    Launch the GUI for analyzing image slices based on the provided configuration.
    """
    args = parse_arguments()

    # Validate config path
    if not args.conf_file_path.exists():
        print(f"❌ Config file not found: {args.conf_file_path}")
        sys.exit(1)

    # Load parameters from configuration file
    try:
        params = read_conf_file(args.conf_file_path)
    except Exception as e:
        print(f"⚠️ Error reading configuration file '{args.conf_file_path}': {e}")
        sys.exit(1)

    # Extract parameters
    mask_path = params.get("MASK_PATH", "")
    output_dir = params.get("OUTPUT_PATH", "./output")
    image_mode = args.image_mode or image_mode_from_config(params)

    # Determine slice number, CLI overrides config
    N_slice = (
        args.N_slice if args.N_slice is not None else params.get("N_SLICE_TEST", 0)
    )

    # Launch Qt app
    app = QApplication(sys.argv)
    w = Window(
        output_dir=str(output_dir),
        mask_path=str(mask_path),
        N_slice=N_slice,
        N_line=args.N_line,
        angle_range=args.angle_range,
        image_mode=image_mode,
    )
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    script()
