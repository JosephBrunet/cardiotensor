#!/usr/bin/env python3
"""
trk_to_am.py
Convert a .trk file to Amira SpatialGraph .am with multiple edge scalars.

Examples
  # default, export all available per-edge scalars, compute elevation if missing
  trk_to_am.py heart.trk

  # export only HA and IA
  trk_to_am.py heart.trk --edge-scalar-sources ha,ia --edge-reduce median
"""

import argparse
from pathlib import Path

from cardiotensor.utils.am_utils import write_spatialgraph_am
from cardiotensor.utils.streamlines_io_utils import (
    compute_elevation_angles,  # derive EL from geometry if needed
    load_trk_streamlines,  # -> (streamlines_xyz, attrs_dict)
    normalize_attrs_to_degrees,  # normalize HA/IA/AZ/EL if byte-scaled
    reduce_per_edge,
)


def parse_sources(s: str) -> list[str]:
    return [t.strip().lower() for t in s.split(",") if t.strip()]


def script():
    parser = argparse.ArgumentParser(
        description="Convert .trk to Amira SpatialGraph .am with multiple edge scalars",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", type=Path, help="Path to .trk")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output .am path. If writer does not support multiple scalars, one file per scalar will be written with suffixes.",
    )
    parser.add_argument(
        "--edge-scalar-sources",
        type=str,
        default="all",
        help=(
            "Comma separated among {ha,ia,az,el,elevation} or 'all'. "
            "Default 'all' exports all angle fields present in the TRK and also computes elevation if EL is absent."
        ),
    )
    parser.add_argument(
        "--edge-reduce",
        choices=["mean", "median"],
        default="mean",
        help="Reduction across points per edge",
    )
    args = parser.parse_args()

    inp = args.input
    if inp.suffix.lower() != ".trk":
        raise ValueError("Only .trk is supported.")
    if not inp.exists():
        raise FileNotFoundError(f"Input file not found: {inp}")

    # Load streamlines and per-point attributes
    streamlines_xyz, attrs = load_trk_streamlines(
        inp
    )  # dict like {"HA":[...], "IA":[...], ...}
    keep = [
        index
        for index, streamline in enumerate(streamlines_xyz)
        if len(streamline) >= 2
    ]
    removed = len(streamlines_xyz) - len(keep)
    if removed:
        print(f"Skipping {removed} streamline(s) with fewer than 2 points")
        streamlines_xyz = [streamlines_xyz[index] for index in keep]
        attrs = {
            name: [values[index] for index in keep] for name, values in attrs.items()
        }
    if not streamlines_xyz:
        raise ValueError("No streamlines with at least 2 points to export")

    attrs_deg = normalize_attrs_to_degrees(
        attrs
    )  # cast to float32 degrees where needed

    # Decide which to export
    requested = parse_sources(args.edge_scalar_sources) or ["all"]

    targets = set()
    if "all" in requested:
        targets.update(k for k in attrs_deg.keys() if k in {"HA", "IA", "AZ", "EL"})
        if "EL" not in targets:
            # compute elevation from geometry
            attrs_deg["EL"] = compute_elevation_angles(streamlines_xyz)
            targets.add("EL")
    else:
        alias = {"elevation": "EL"}
        for item in requested:
            key = alias.get(item, item).upper()
            if key == "EL":
                if "EL" not in attrs_deg:
                    attrs_deg["EL"] = compute_elevation_angles(streamlines_xyz)
                targets.add("EL")
            elif key in {"HA", "IA", "AZ"}:
                if key not in attrs_deg:
                    raise ValueError(
                        f"Requested '{key}' is not present in TRK per-point data. "
                        f"Available: {sorted(attrs_deg.keys())}"
                    )
                targets.add(key)
            else:
                raise ValueError(
                    "Sources must be among {ha, ia, az, el, elevation} or 'all'"
                )

    if not targets:
        raise ValueError("No edge scalars to export after filtering.")

    # Reduce per-point to per-edge
    edge_scalar_dict = {}
    for key in sorted(targets):
        edge_scalar_dict[key] = reduce_per_edge(attrs_deg[key], how=args.edge_reduce)

    base_out = args.output if args.output is not None else inp.with_suffix(".am")

    # Try multi field write if your writer supports dict, else fallback per field
    try:
        write_spatialgraph_am(
            out_path=base_out,
            streamlines_xyz=streamlines_xyz,
            point_thickness=None,
            edge_scalar=edge_scalar_dict,  # dict[str, np.ndarray]
            edge_scalar_name=None,  # ignored when dict is provided
        )
        print(
            f"Wrote Amira SpatialGraph with fields {sorted(edge_scalar_dict.keys())}: {base_out}"
        )
    except TypeError:
        print(
            "Writer does not accept multiple edge scalars, writing one file per scalar..."
        )
        for key, vals in edge_scalar_dict.items():
            out_path = base_out.with_name(f"{base_out.stem}_{key}{base_out.suffix}")
            write_spatialgraph_am(
                out_path=out_path,
                streamlines_xyz=streamlines_xyz,
                point_thickness=None,
                edge_scalar=vals,
                edge_scalar_name=key,
            )
            print(f"  -> {out_path}")


if __name__ == "__main__":
    script()
