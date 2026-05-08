#!/usr/bin/env python3
"""
visualize_streamlines.py
------------------------
CLI tool to visualize cardiac streamlines from a .trk file using FURY.

Input
  1) a .conf file with OUTPUT_PATH, will use OUTPUT_PATH/streamlines.trk
  2) a direct path to a .trk file

Color-by
  - elevation  computed from streamline geometry
  - any per-point scalar stored in the TRK data_per_point, e.g. HA, IA, AZ, EL
  - use --list-color-by to print all available options and exit
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib

from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.utils import read_conf_file
from cardiotensor.visualization.streamlines import visualize_streamlines


def _discover_trk_color_fields(trk_path: Path) -> list[str]:
    """Return list of available per-point scalar keys in TRK, case-preserving."""
    obj = nib.streamlines.load(str(trk_path))
    tg = obj.tractogram
    dpp = getattr(tg, "data_per_point", None)
    if not dpp:
        return []
    # keys can be e.g. "HA", "IA", "AZ", "EL", or custom
    return list(dpp.keys())


def _resolve_streamlines_path(input_path: Path) -> Path:
    if not input_path.exists():
        print(f"Input path not found: {input_path}")
        sys.exit(1)
    suf = input_path.suffix.lower()
    if suf == ".conf":
        params = read_conf_file(input_path)
        out_dir = Path(params.get("OUTPUT_PATH", "./output"))
        trk = out_dir / "streamlines.trk"
        if not trk.exists():
            print(f"streamlines.trk not found at: {trk}")
            sys.exit(1)
        return trk
    if suf == ".trk":
        return input_path
    print("Unsupported input. Provide either a .conf file or a .trk file.")
    sys.exit(1)


def script() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize cardiac streamlines from .trk, color by elevation or stored per-point fields",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help="Path to a .conf file with OUTPUT_PATH or directly a .trk file",
    )
    parser.add_argument(
        "--color-by",
        type=str,
        default="auto",
        help="Which scalar to color by. Options include 'elevation' or any per-point field present in the TRK, e.g. HA, IA, AZ, EL. Use --list-color-by to see what is available. 'auto' prefers HA, else any available field, else elevation.",
    )
    parser.add_argument(
        "--list-color-by",
        action="store_true",
        help="List available color-by options from the .trk and exit",
    )
    parser.add_argument(
        "--line-width",
        type=float,
        default=4.0,
        help="Tube width",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=1,
        help="Keep every Nth streamline for faster rendering",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="Filter out streamlines shorter than this number of vertices",
    )
    parser.add_argument(
        "--downsample-factor",
        type=int,
        default=1,
        help="Keep every Nth vertex along each streamline",
    )
    parser.add_argument(
        "--max-streamlines",
        type=int,
        default=None,
        help="Maximum number of streamlines to render after filtering",
    )
    parser.add_argument(
        "--crop-x",
        nargs=2,
        type=float,
        metavar=("XMIN", "XMAX"),
        help="Crop to X range",
    )
    parser.add_argument(
        "--crop-y",
        nargs=2,
        type=float,
        metavar=("YMIN", "YMAX"),
        help="Crop to Y range",
    )
    parser.add_argument(
        "--crop-z",
        nargs=2,
        type=float,
        metavar=("ZMIN", "ZMAX"),
        help="Crop to Z range",
    )
    parser.add_argument(
        "--no-interactive",
        action="store_true",
        help="Disable interactive window, useful when only saving a screenshot",
    )
    parser.add_argument(
        "--screenshot",
        type=str,
        default=None,
        help="Path to save a PNG screenshot, no file is saved if omitted",
    )
    parser.add_argument(
        "--spline-subdiv",
        type=int,
        default=16,
        help="Spline subdivisions for tube rendering (higher = smoother, heavier)",
    )
    parser.add_argument("--width", type=int, default=800, help="Window width in pixels")
    parser.add_argument(
        "--height", type=int, default=800, help="Window height in pixels"
    )
    parser.add_argument(
        "--colormap",
        type=str,
        default=None,
        help="Colormap name. By default, choose one from --color-by: helix_angle for HA/IA, viridis for EL, hsv for AZ.",
    )

    args = parser.parse_args()

    # Resolve .trk
    trk_path = _resolve_streamlines_path(args.input_path)

    # List available color-by and exit if requested
    available_dpp = _discover_trk_color_fields(trk_path)
    available = ["elevation"] + available_dpp
    if args.list_color_by:
        print("Available color-by options:")
        for k in available:
            print(f"  - {k}")
        return

    # Determine color-by
    color_by = args.color_by.strip().lower()
    dpp_lower = {k.lower(): k for k in available_dpp}  # map lower to original key

    if color_by == "auto":
        if "ha" in dpp_lower:
            color_by = "ha"
        elif available_dpp:
            # pick the first available stored field
            color_by = list(dpp_lower.keys())[0]
        else:
            color_by = "elevation"
        print(f"Auto color-by resolved to: {color_by}")

    # Validate choice
    if color_by != "elevation" and color_by not in dpp_lower:
        print(
            f"Requested color-by '{args.color_by}' not found. Available: {', '.join(available)}"
        )
        sys.exit(2)

    # Choose colormap only when explicitly requested; otherwise the visualizer
    # picks the default based on color_by.
    chosen_cmap = None
    if args.colormap is not None:
        if args.colormap.lower() == "helix_angle":
            chosen_cmap = helix_angle_cmap
        else:
            try:
                chosen_cmap = plt.get_cmap(args.colormap)
            except ValueError:
                print(
                    f"Unknown colormap '{args.colormap}', falling back to helix_angle"
                )
                chosen_cmap = helix_angle_cmap

    # Crop bounds
    crop_bounds = None
    if args.crop_x or args.crop_y or args.crop_z:
        crop_bounds = (
            tuple(args.crop_x or [-float("inf"), float("inf")]),
            tuple(args.crop_y or [-float("inf"), float("inf")]),
            tuple(args.crop_z or [-float("inf"), float("inf")]),
        )

    # Call the visualizer
    # Always render as tubes as requested earlier
    visualize_streamlines(
        streamlines_file=trk_path,
        color_by=(
            color_by if color_by == "elevation" else dpp_lower[color_by]
        ),  # pass original key if stored
        line_width=args.line_width,
        subsample_factor=args.subsample,
        filter_min_len=args.min_length,
        downsample_factor=args.downsample_factor,
        max_streamlines=args.max_streamlines,
        crop_bounds=crop_bounds,
        interactive=not args.no_interactive,
        screenshot_path=args.screenshot,
        window_size=(args.width, args.height),
        colormap=chosen_cmap,
        spline_subdiv=args.spline_subdiv,
    )


if __name__ == "__main__":
    script()
