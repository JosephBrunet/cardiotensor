#!/usr/bin/env python3
"""
visualize_streamlines.py

Command-line tool to visualize cardiac streamlines with color encoding.

This script loads streamlines from a .npz file (generated via a .conf pipeline),
applies optional filtering, downsampling, and subsampling, then renders the streamlines
as tubes or lines, colored either by helix angle (HA) or elevation angle.

Example usage:
    cardio-streamlines-visualize parameters.conf --color-by ha --mode tube --subsample 10 --screenshot heart.png

Options:
    --color-by            Choose streamline color mode: 'ha' (helix angle) or 'elevation' (angle w.r.t XY plane) [default: ha]
    --mode                Visual style: 'tube', 'fake_tube', or 'line' [default: tube]
    --line-width          Width of lines or tubes in pixels [default: 4]
    --subsample           Keep 1 in N streamlines [default: 1]
    --min-length          Minimum number of points per streamline [default: 10]
    --downsample-factor   Downsample streamline points by this factor [default: 1]
    --max-streamlines     Maximum number of streamlines to plot [default: 1000]
    --no-interactive      If set, disables GUI and requires --screenshot
    --screenshot          Save to PNG instead of displaying interactively
    --width / --height    Render window or screenshot size in pixels [default: 800 x 800]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from src.visualization.visualize_streamlines import show_streamlines
from cardiotensor.utils.utils import read_conf_file


def compute_streamline_bounds(streamlines):
    """
    Compute min and max bounds in X, Y, Z for a list of 3D streamlines.

    Returns:
        ((x_min, x_max), (y_min, y_max), (z_min, z_max))
    """
    all_points = np.concatenate(streamlines, axis=0)
    x_min, y_min, z_min = all_points.min(axis=0)
    x_max, y_max, z_max = all_points.max(axis=0)
    return (x_min, x_max), (y_min, y_max), (z_min, z_max)



def script():
    parser = argparse.ArgumentParser(
        description="Visualize streamlines from .npz file."
    )
    parser.add_argument(
        "conf_file", type=str, help="Path to the .conf file used for generation."
    )

    parser.add_argument(
        "--color-by",
        choices=["ha", "elevation"],
        default="ha",
        help="Color streamlines by helix angle (ha) or elevation angle (elevation). Default: ha.",
    )

    parser.add_argument(
        "--mode",
        choices=["tube", "fake_tube", "line"],
        default="tube",
        help="Rendering style: tube (default), fake_tube, or line.",
    )

    parser.add_argument(
        "--line-width",
        type=float,
        default=4,
        help="Line or tube thickness in pixels. Default: 4.",
    )

    parser.add_argument(
        "--subsample",
        type=int,
        default=1,
        help="Randomly keep 1 in every N streamlines. Default: 1 (keep all).",
    )

    parser.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="Minimum number of points per streamline. If not set, no filtering is applied.",
    )

    parser.add_argument(
        "--downsample-factor",
        type=int,
        default=1,
        help="Downsample points within each streamline by this factor. Default: 1 (no downsampling).",
    )

    parser.add_argument(
        "--max-streamlines",
        type=int,
        default=None,
        help="Maximum number of streamlines to render. If not set, no limit is applied.",
    )
    
    parser.add_argument(
        "--crop-z", nargs=2, type=float, metavar=("ZMIN", "ZMAX"),
        help="Crop streamlines to this Z range (after conversion to (x,y,z))."
    )
    parser.add_argument(
        "--crop-y", nargs=2, type=float, metavar=("YMIN", "YMAX"),
        help="Crop streamlines to this Y range."
    )
    parser.add_argument(
        "--crop-x", nargs=2, type=float, metavar=("XMIN", "XMAX"),
        help="Crop streamlines to this X range."
    )

    parser.add_argument(
        "--no-interactive",
        action="store_true",
        help="Run in headless mode (no GUI). Requires --screenshot to be set.",
    )

    parser.add_argument(
        "--screenshot",
        type=str,
        default=None,
        help="Path to PNG screenshot (used only if --no-interactive is set).",
    )

    parser.add_argument(
        "--width",
        type=int,
        default=800,
        help="Window or screenshot width in pixels. Default: 800.",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=800,
        help="Window or screenshot height in pixels. Default: 800.",
    )


    args = parser.parse_args()
    conf_path = Path(args.conf_file)

    if not conf_path.exists():
        print(f"❌ Config file not found: {conf_path}")
        sys.exit(1)

    try:
        params = read_conf_file(conf_path)
    except Exception as e:
        print(f"⚠️ Error reading config: {e}")
        sys.exit(1)

    output_dir = Path(params.get("OUTPUT_PATH", "./output"))
    npz_path = output_dir / "streamlines.npz"

    if not npz_path.exists():
        print(f"❌ Missing streamlines.npz in {output_dir}")
        sys.exit(1)

    data = np.load(npz_path, allow_pickle=True)
    raw_streamlines = data.get("streamlines")

    if raw_streamlines is None:
        print("❌ 'streamlines' array missing in .npz.")
        sys.exit(1)

    # Convert (z, y, x) to (x, y, z)
    streamlines_xyz = [
        np.array([(pt[2], pt[1], pt[0]) for pt in sl], dtype=np.float32)
        for sl in raw_streamlines.tolist()
    ]

    if args.color_by == "elevation":
        from src.visualization.visualize_streamlines import (
            compute_elevation_angles,
        )

        color_values = compute_elevation_angles(streamlines_xyz)
    else:
        ha_values = data.get("ha_values")

        if ha_values is None:
            print("❌ 'ha_values' array missing in .npz.")
            sys.exit(1)

        # Convert from [0, 255] to [-90, 90]
        ha_values = ha_values.astype(np.float32)
        color_values = (ha_values / 255.0) * 180.0 - 90.0
        
    # Build crop bounds with open-ended defaults
    crop_z = tuple(args.crop_z) if args.crop_z else (-float("inf"), float("inf"))
    crop_y = tuple(args.crop_y) if args.crop_y else (-float("inf"), float("inf"))
    crop_x = tuple(args.crop_x) if args.crop_x else (-float("inf"), float("inf"))

    if any([args.crop_z, args.crop_y, args.crop_x]):
        crop_bounds = (crop_x, crop_y, crop_z)
    else:
        crop_bounds = None

    print(f"Original bounds for streamlines:")
    x_bounds, y_bounds, z_bounds = compute_streamline_bounds(streamlines_xyz)
    print(f"X bounds: {int(x_bounds[0])} - {int(x_bounds[1])}")
    print(f"Y bounds: {int(y_bounds[0])} - {int(y_bounds[1])}")
    print(f"Z bounds: {int(z_bounds[0])} - {int(z_bounds[1])}")
        
    show_streamlines(
        streamlines_xyz=streamlines_xyz,
        color_values=color_values,
        mode=args.mode,
        line_width=args.line_width,
        interactive=not args.no_interactive,
        screenshot_path=args.screenshot,
        window_size=(args.width, args.height),
        downsample_factor=args.downsample_factor,
        max_streamlines=args.max_streamlines,
        filter_min_len=args.min_length,
        subsample_factor=args.subsample,
        crop_bounds=crop_bounds,
    )


if __name__ == "__main__":
    script()
