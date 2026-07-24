"""
3D_Data_Processing
"""

import argparse

from cardiotensor.launcher.slurm_launcher import slurm_launcher


def script() -> None:
    """
    Submit cardiotensor processing as SLURM array jobs.
    """

    parser = argparse.ArgumentParser(
        description="Launch cardio-tensor jobs on a SLURM cluster.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "conf_file_path", type=str, help="Path to the input configuration file."
    )
    parser.add_argument(
        "--start_index",
        type=int,
        default=0,
        help="Starting slice index for processing (default: 0).",
    )
    parser.add_argument(
        "--end_index",
        type=int,
        default=None,
        help="Ending slice index (exclusive). Default: None (process until end).",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=None,
        help="Slices per SLURM task. If omitted, uses N_CHUNK from config.",
    )
    parser.add_argument(
        "--partition",
        type=str,
        default=None,
        help="SLURM partition name. If omitted, no partition directive is added.",
    )
    parser.add_argument(
        "--time_limit",
        type=str,
        default="2:00:00",
        help="SLURM time limit (e.g., 2:00:00).",
    )
    parser.add_argument(
        "--cpus_per_task",
        type=int,
        default=8,
        help="SLURM CPUs per task.",
    )
    parser.add_argument(
        "--mem_gb",
        type=int,
        default=64,
        help="SLURM memory per task in GB.",
    )
    parser.add_argument(
        "--array_parallel",
        type=int,
        default=100,
        help="Maximum concurrent tasks in SLURM array.",
    )
    parser.add_argument(
        "--log_dir",
        type=str,
        default=None,
        help="Directory for SLURM logs. If omitted, uses OUTPUT_PATH/slurm/log.",
    )
    parser.add_argument(
        "--submit_dir",
        type=str,
        default=None,
        help="Directory where .slurm scripts are stored. If omitted, uses OUTPUT_PATH/slurm/submit.",
    )
    parser.add_argument(
        "--no_monitor",
        action="store_true",
        help="Do not wait and monitor output files after submission.",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Generate scripts and print sbatch commands without submitting.",
    )
    args = parser.parse_args()
    conf_file_path = args.conf_file_path

    try:
        slurm_launcher(
            conf_file_path,
            start_index=args.start_index,
            end_index=args.end_index,
            chunk_size=args.chunk_size,
            partition=args.partition,
            time_limit=args.time_limit,
            cpus_per_task=args.cpus_per_task,
            mem_gb=args.mem_gb,
            array_parallel=args.array_parallel,
            log_dir=args.log_dir,
            submit_dir=args.submit_dir,
            monitor=not args.no_monitor,
            dry_run=args.dry_run,
        )
    except Exception as err:
        parser.exit(1, f"Error: {err}\n")


if __name__ == "__main__":
    script()
