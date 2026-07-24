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
import json
import random
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib

from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.utils import read_conf_file
from cardiotensor.visualization.streamlines import visualize_streamlines


def _discover_trk_color_fields(trk_path: Path) -> list[str]:
    """Return list of available per-point scalar keys in TRK, case-preserving."""
    # Only the TRK header is needed here. A normal load would read the whole
    # tractogram once for discovery and then a second time for visualization.
    obj = nib.streamlines.load(str(trk_path), lazy_load=True)
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


def _apply_session_settings(args, settings: dict, provided_options: set[str]) -> None:
    """Use saved values as defaults while keeping explicit CLI options."""
    option_map = {
        "color_by": ("color_by", "--color-by"),
        "line_width": ("line_width", "--line-width"),
        "subsample_factor": ("subsample", "--subsample"),
        "filter_min_len": ("min_length", "--min-length"),
        "downsample_factor": ("downsample_factor", "--downsample-factor"),
        "max_streamlines": ("max_streamlines", "--max-streamlines"),
        "top_clusters": ("top_clusters", "--top-clusters"),
        "spline_subdiv": ("spline_subdiv", "--spline-subdiv"),
        "colormap": ("colormap", "--colormap"),
        "backend": ("backend", "--backend"),
        "mode": ("mode", "--mode"),
        "background_color": ("background_color", "--background-color"),
        "tube_sides": ("tube_sides", "--tube-sides"),
        "opacity": ("opacity", "--opacity"),
        "quality": ("quality", "--quality"),
        "random_seed": ("random_seed", "--random-seed"),
    }
    for key, (attribute, option) in option_map.items():
        if option not in provided_options and key in settings:
            setattr(args, attribute, settings[key])

    if "--line-width" not in provided_options and "tube_thickness" in settings:
        args.line_width = float(settings["tube_thickness"])

    if not {"--crop-x", "--crop-y", "--crop-z"} & provided_options:
        crop_bounds = settings.get("crop_bounds")
        if crop_bounds:
            args.crop_x, args.crop_y, args.crop_z = crop_bounds

    window_size = settings.get("window_size")
    if window_size:
        if "--width" not in provided_options:
            args.width = int(window_size[0])
        if "--height" not in provided_options:
            args.height = int(window_size[1])

    if "--hide-axes" not in provided_options and "show_axes" in settings:
        args.hide_axes = not bool(settings["show_axes"])
    if not {"--show-bounds", "--hide-bounds"} & provided_options:
        args.show_bounds = bool(settings.get("show_bounds", False))
        args.hide_bounds = False
    if not {"--shadows", "--no-shadows"} & provided_options:
        args.shadows = bool(settings.get("shadows", False))
        args.no_shadows = False


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
        "--session",
        type=Path,
        default=None,
        help=(
            "FURY session JSON to load and update. When omitted, use "
            "streamlines_session.json beside the TRK. Press Ctrl+S in the viewer "
            "to save the current camera, zoom, crop gizmos and controls"
        ),
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Seed used for reproducible streamline subsampling",
    )
    parser.add_argument(
        "--line-width",
        type=float,
        default=4.0,
        help="Tube radius in tube mode, or screen line width in line mode",
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
        "--top-clusters",
        type=int,
        default=None,
        help=(
            "Keep all streamlines from the N largest clusters; requires "
            "cluster_id metadata"
        ),
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
        help="Disable interactive window, useful when only saving a screenshot or video",
    )
    parser.add_argument(
        "--screenshot",
        type=str,
        default=None,
        help="Path to save a PNG screenshot, no file is saved if omitted",
    )
    parser.add_argument(
        "--video",
        type=str,
        default=None,
        metavar="PATH",
        help=(
            "Path for an orbit video (e.g. orbit.mp4 or orbit.gif). "
            "Implies --no-interactive. FURY MP4 uses OpenCV; PyVista MP4 may require imageio-ffmpeg."
        ),
    )
    parser.add_argument(
        "--video-fps",
        type=int,
        default=30,
        help="Frames per second for the orbit video (default: 30)",
    )
    parser.add_argument(
        "--video-frames",
        type=int,
        default=120,
        help="Total number of frames in the orbit video (default: 120 = 4 s at 30 fps)",
    )
    parser.add_argument(
        "--spline-subdiv",
        type=int,
        default=None,
        help=(
            "Number of cardinal-spline segments per streamline; 1 disables "
            "smoothing (quality preset default when omitted)"
        ),
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
    parser.add_argument(
        "--backend",
        choices=("fury", "pyvista"),
        default="fury",
        help="Rendering backend",
    )
    parser.add_argument(
        "--mode",
        choices=("tube", "line"),
        default="tube",
        help="Render explicit tube geometry or faster line geometry",
    )
    parser.add_argument(
        "--background-color",
        type=str,
        default=None,
        help="Renderer background color. Default: white for PyVista, black for FURY",
    )
    parser.add_argument(
        "--tube-sides",
        type=int,
        default=None,
        help="Number of sides for tube geometry (quality preset default when omitted)",
    )
    parser.add_argument(
        "--quality",
        choices=("interactive", "publication"),
        default="interactive",
        help="FURY rendering preset; individual options still override it",
    )
    parser.add_argument(
        "--opacity",
        type=float,
        default=1.0,
        help="Streamline opacity",
    )
    parser.add_argument(
        "--shadows",
        action="store_true",
        help="Enable PyVista shadow rendering at startup",
    )
    parser.add_argument(
        "--no-shadows",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--hide-axes",
        action="store_true",
        help="Hide the PyVista orientation axes",
    )
    parser.add_argument(
        "--show-bounds",
        action="store_true",
        help="Show the PyVista bounds grid at startup",
    )
    parser.add_argument(
        "--hide-bounds",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    provided_options = {
        token.split("=", 1)[0] for token in sys.argv[1:] if token.startswith("--")
    }
    args = parser.parse_args()

    # Resolve .trk
    trk_path = _resolve_streamlines_path(args.input_path)

    restore_session = args.session is not None
    if args.session is None:
        args.session = (
            trk_path.resolve().parent / "streamlines_session.json"
        ).resolve()
        restore_session = False
        if args.session.exists() and not args.list_color_by:
            if sys.stdin.isatty():
                answer = (
                    input(
                        f"Found default visualization session {args.session}. "
                        "Restart from it? [Y/n] "
                    )
                    .strip()
                    .lower()
                )
                restore_session = answer in {"", "y", "yes"}
            else:
                print(
                    f"Found default visualization session {args.session}, but input "
                    "is not interactive. Starting a fresh view; pass --session "
                    "explicitly to restore it."
                )
        elif not args.list_color_by:
            print(f"Visualization session will be saved by default to {args.session}")
    else:
        args.session = args.session.expanduser().resolve()

    if not args.list_color_by:
        if args.session.exists() and restore_session:
            try:
                session_data = json.loads(args.session.read_text())
                saved_file = session_data.get("streamlines_file")
                saved_size = session_data.get("streamlines_size")
                if saved_file and Path(saved_file).resolve() != trk_path.resolve():
                    print(
                        "Warning: session was created for a different TRK path: "
                        f"{saved_file}"
                    )
                if saved_size and int(saved_size) != trk_path.stat().st_size:
                    print("Warning: session TRK size differs from the current file.")
                _apply_session_settings(
                    args, session_data.get("settings", {}), provided_options
                )
                print(f"Loaded visualization settings from {args.session}")
            except (OSError, ValueError, TypeError, json.JSONDecodeError) as err:
                print(f"Warning: could not load session {args.session}: {err}")
                restore_session = False
        elif args.session.exists():
            print(f"Starting a fresh visualization session at {args.session}")
        elif "--session" in provided_options:
            print(f"New visualization session will be saved to {args.session}")

    if args.top_clusters is not None and args.top_clusters < 1:
        parser.error("--top-clusters must be at least 1")

    if args.random_seed is None:
        args.random_seed = random.SystemRandom().randrange(2**63)
    if args.spline_subdiv is None:
        args.spline_subdiv = 1 if args.quality == "interactive" else 2
    if args.tube_sides is None:
        args.tube_sides = 6 if args.quality == "interactive" else 12

    # List available color-by and exit if requested
    available_dpp = _discover_trk_color_fields(trk_path)
    computed_color_options = ["elevation", "azimuth", "az", "el"]
    available = computed_color_options + available_dpp
    if args.list_color_by:
        print("Available color-by options:")
        for k in available:
            print(f"  - {k}")
        return

    # Determine color-by
    color_by = args.color_by.strip().lower().replace("-", "_")
    dpp_lower = {
        key.lower().replace("-", "_"): key for key in available_dpp
    }

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
    if color_by not in computed_color_options and color_by not in dpp_lower:
        print(
            f"Requested color-by '{args.color_by}' not found. Available: {', '.join(available)}"
        )
        sys.exit(2)

    selected_color_by = dpp_lower[color_by] if color_by in dpp_lower else color_by

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

    # --video implies off-screen
    is_interactive = not args.no_interactive and not args.video

    session_settings = {
        "color_by": selected_color_by,
        "line_width": args.line_width,
        "subsample_factor": args.subsample,
        "filter_min_len": args.min_length,
        "downsample_factor": args.downsample_factor,
        "max_streamlines": args.max_streamlines,
        "top_clusters": args.top_clusters,
        "crop_bounds": crop_bounds,
        "window_size": [args.width, args.height],
        "colormap": args.colormap,
        "spline_subdiv": args.spline_subdiv,
        "backend": args.backend,
        "mode": args.mode,
        "background_color": args.background_color,
        "tube_sides": args.tube_sides,
        "opacity": args.opacity,
        "quality": args.quality,
        "show_axes": not args.hide_axes,
        "show_bounds": args.show_bounds and not args.hide_bounds,
        "shadows": args.shadows and not args.no_shadows,
        "random_seed": args.random_seed,
    }

    # Call the visualizer
    visualize_streamlines(
        streamlines_file=trk_path,
        color_by=selected_color_by,
        line_width=args.line_width,
        subsample_factor=args.subsample,
        filter_min_len=args.min_length,
        downsample_factor=args.downsample_factor,
        max_streamlines=args.max_streamlines,
        top_clusters=args.top_clusters,
        crop_bounds=crop_bounds,
        interactive=is_interactive,
        screenshot_path=args.screenshot,
        video_path=args.video,
        video_fps=args.video_fps,
        video_frames=args.video_frames,
        window_size=(args.width, args.height),
        colormap=chosen_cmap,
        spline_subdiv=args.spline_subdiv,
        backend=args.backend,
        mode=args.mode,
        background_color=args.background_color,
        tube_sides=args.tube_sides,
        pyvista_opacity=args.opacity,
        pyvista_show_axes=not args.hide_axes,
        pyvista_show_bounds=args.show_bounds and not args.hide_bounds,
        pyvista_shadows=args.shadows and not args.no_shadows,
        random_seed=args.random_seed,
        session_path=args.session,
        restore_session=restore_session,
        session_settings=session_settings,
        fury_quality=args.quality,
    )


if __name__ == "__main__":
    script()
