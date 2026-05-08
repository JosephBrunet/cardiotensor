#!/usr/bin/env python3
"""
generate_streamlines.py
-----------------------
CLI script to generate streamlines from a 3D eigenvector field.
Reads a .conf file for paths and parameters, then calls the library function.

Behavior
- Writes .trk with all discovered per-point angle fields
- Writes .am with per-edge mean fields for all discovered angles
"""

import argparse
import sys
from pathlib import Path

from cardiotensor.tractography.generate_streamlines import (
    generate_streamlines_from_params,
)
from cardiotensor.utils.utils import read_conf_file


def script() -> None:
    parser = argparse.ArgumentParser(
        description="Generate streamlines from a 3D vector field and export .trk and .am",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("conf_file", type=Path, help="Path to configuration .conf file")
    parser.add_argument("--start-z", type=int, default=0, help="Start slice index in Z")
    parser.add_argument("--end-z", type=int, default=None, help="End slice index in Z")
    parser.add_argument("--start-y", type=int, default=0, help="Start slice index in Y")
    parser.add_argument("--end-y", type=int, default=None, help="End slice index in Y")
    parser.add_argument("--start-x", type=int, default=0, help="Start slice index in X")
    parser.add_argument("--end-x", type=int, default=None, help="End slice index in X")
    parser.add_argument("--bin", type=int, default=1, help="Downsampling factor")
    parser.add_argument("--seeds", type=int, default=20000, help="Number of seeds")
    parser.add_argument(
        "--fa-seed-min", type=float, default=0.2, help="Min FA for seeding"
    )
    parser.add_argument("--fa-threshold", type=float, default=0.1, help="FA threshold")
    parser.add_argument("--step", type=float, default=0.5, help="Step length in voxels")
    parser.add_argument(
        "--max-steps", type=int, default=None, help="Max steps per streamline"
    )
    parser.add_argument(
        "--angle", type=float, default=60.0, help="Max turning angle in degrees"
    )
    parser.add_argument(
        "--min-len", type=int, default=10, help="Minimum streamline length in points"
    )

    # Base angle folder selection for discovery
    parser.add_argument(
        "--angle-mode",
        choices=["ha_ia", "az_el"],
        help="Select the base angle folder used for discovery. "
        "ha_ia chooses HA/, az_el chooses AZ/. All siblings (HA IA AZ EL) found next to it are included.",
    )

    args = parser.parse_args()

    # Load config
    try:
        params = read_conf_file(args.conf_file)
    except Exception as e:
        print(f"Failed to read config file: {e}")
        sys.exit(1)

    angle_mode = (args.angle_mode or params.get("ANGLE_MODE", "ha_ia")).lower().strip()
    if angle_mode not in {"ha_ia", "az_el"}:
        print(f"Invalid ANGLE_MODE '{angle_mode}'. Use 'ha_ia' or 'az_el'.")
        sys.exit(2)

    # Pick base angle folder for discovery
    angle_folder = "HA" if angle_mode == "ha_ia" else "AZ"

    # Resolve standard output structure
    output_dir = Path(params.get("OUTPUT_PATH", "./output"))
    vector_field_dir = output_dir / "eigen_vec"
    fa_dir = output_dir / "FA"
    angle_dir = output_dir / angle_folder
    mask_path = params.get("MASK_PATH", None)

    # Call library
    generate_streamlines_from_params(
        vector_field_dir=vector_field_dir,
        output_dir=output_dir,
        fa_dir=fa_dir,
        angle_dir=angle_dir,  # base folder, siblings discovered automatically
        mask_path=mask_path,
        start_xyz=(args.start_z, args.start_y, args.start_x),
        end_xyz=(args.end_z, args.end_y, args.end_x),
        bin_factor=args.bin,
        num_seeds=args.seeds,
        fa_seed_min=args.fa_seed_min,
        fa_threshold=args.fa_threshold,
        step_length=args.step,
        max_steps=args.max_steps,
        angle_threshold=args.angle,
        min_length_pts=args.min_len,
        angle_mode=angle_mode,
    )

    print(
        f"Done! Streamlines generated and saved to {output_dir}."
    )


if __name__ == "__main__":
    script()
