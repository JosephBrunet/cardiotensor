import argparse
from pathlib import Path

from cardiotensor.utils.utils import read_conf_file
from cardiotensor.visualization.streamlines import visualize_streamlines


def script():
    parser = argparse.ArgumentParser(
        description="Visualize cardiac streamlines from a .npz file."
    )
    parser.add_argument("conf_file", type=Path, help="Path to configuration file")
    parser.add_argument("--color-by", choices=["ha", "elevation"], default="ha")
    parser.add_argument("--mode", choices=["tube", "fake_tube", "line"], default="tube")
    parser.add_argument("--line-width", type=float, default=4)
    parser.add_argument("--subsample", type=int, default=1)
    parser.add_argument("--min-length", type=int, default=None)
    parser.add_argument("--downsample-factor", type=int, default=1)
    parser.add_argument("--max-streamlines", type=int, default=None)
    parser.add_argument("--crop-x", nargs=2, type=float)
    parser.add_argument("--crop-y", nargs=2, type=float)
    parser.add_argument("--crop-z", nargs=2, type=float)
    parser.add_argument("--no-interactive", action="store_true")
    parser.add_argument("--screenshot", type=str, default=None)
    parser.add_argument("--width", type=int, default=800)
    parser.add_argument("--height", type=int, default=800)

    args = parser.parse_args()

    # Read conf file only in CLI
    params = read_conf_file(args.conf_file)
    output_dir = Path(params.get("OUTPUT_PATH", "./output"))
    streamlines_file = output_dir / "streamlines.npz"

    # Compute crop bounds if provided
    crop_bounds = None
    if args.crop_x or args.crop_y or args.crop_z:
        crop_bounds = (
            tuple(args.crop_x or [-float("inf"), float("inf")]),
            tuple(args.crop_y or [-float("inf"), float("inf")]),
            tuple(args.crop_z or [-float("inf"), float("inf")]),
        )

    # Call Python function
    visualize_streamlines(
        streamlines_file=streamlines_file,
        color_by=args.color_by,
        mode=args.mode,
        line_width=args.line_width,
        subsample_factor=args.subsample,
        filter_min_len=args.min_length,
        downsample_factor=args.downsample_factor,
        max_streamlines=args.max_streamlines,
        crop_bounds=crop_bounds,
        interactive=not args.no_interactive,
        screenshot_path=args.screenshot,
        window_size=(args.width, args.height),
    )


if __name__ == "__main__":
    script()
