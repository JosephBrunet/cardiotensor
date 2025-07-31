import argparse

from src.cardiotensor.tractography.export_to_amira import amira_writer


def script() -> None:
    parser = argparse.ArgumentParser(
        description="Convert images between tif and jpeg2000 formats"
    )
    parser.add_argument("conf_file_path", type=str, help="Path to the input text file")
    parser.add_argument(
        "--start_index", type=int, default=None, help="Start index for volume subset"
    )
    parser.add_argument(
        "--end_index", type=int, default=None, help="End index for volume subset"
    )
    parser.add_argument("--bin", type=int, default=None, help="binning volume")
    parser.add_argument(
        "--num_ini_points",
        type=int,
        default=1000,
        help="Number of starting random points.",
    )
    parser.add_argument(
        "--num_steps",
        type=int,
        default=5000,
        help="Number of steps to follow vectors.",
    )
    parser.add_argument(
        "--max_length",
        type=int,
        default=500,
        help="Length of each segment in voxel.",
    )
    parser.add_argument(
        "--angle_threshold",
        type=float,
        default=90,
        help="Maximum allowed angle change.",
    )
    parser.add_argument(
        "--segment_min_length_threshold",
        type=int,
        default=5,
        help="Minimum length of valid fibers.",
    )

    args = parser.parse_args()

    conf_file_path = args.conf_file_path
    start_index = args.start_index
    end_index = args.end_index
    bin_factor = args.bin
    num_ini_points = args.num_ini_points
    num_steps = args.num_steps
    max_length = args.max_length
    angle_threshold = args.angle_threshold
    segment_min_length_threshold = args.segment_min_length_threshold

    amira_writer(
        conf_file_path,
        start_index,
        end_index,
        bin_factor,
        num_ini_points,
        num_steps,
        max_length,
        angle_threshold,
        segment_min_length_threshold,
    )
