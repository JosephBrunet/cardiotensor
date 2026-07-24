#!/usr/bin/env python3
"""
streamline_compare.py — histogram by default, KDE optional

Distributions:
- Length (arc length), Curvature, Tortuosity

Features:
- Histogram (counts by default), KDE optional
- --normalize rescales each curve so y ∈ [0,1] (relative shape only)
- Percentile clipping via --clip PLOW PHIGH
- Axis scales selectable via --xscale/--yscale (linear/log)
- Saves PNG + PDF
"""

from __future__ import annotations

import argparse
import csv
import itertools
from pathlib import Path
from typing import Iterator

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update(
    {
        "savefig.dpi": 300,
        "savefig.transparent": True,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 120,
        "font.size": 14,
        "axes.labelsize": 14,
        "legend.fontsize": 14,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.frameon": False,
    }
)

# ---------- helpers ----------


def iter_streamlines(path: Path, key: str | None = None) -> Iterator[np.ndarray]:
    """
    Yield streamline geometry without loading unused TRK per-point fields.

    Nibabel applies the TRK affine lazily, so yielded TRK coordinates are in
    the physical space encoded in the file header. Legacy CardioTensor TRKs
    used an identity affine and can still be scaled with --voxel-size.
    """
    suffix = path.suffix.lower()

    if suffix == ".trk":
        import nibabel as nib

        tractogram = nib.streamlines.load(str(path), lazy_load=True).tractogram
        source = tractogram.streamlines
    elif suffix in {".npz", ".npy"}:
        data = np.load(path, allow_pickle=True)
        source = data[key or "streamlines"] if suffix == ".npz" else data
    else:
        raise ValueError(f"Unsupported file type: {path.suffix}")

    for streamline in source:
        points = np.asarray(streamline, dtype=np.float32)
        if points.ndim == 2 and points.shape[1] == 3 and len(points) >= 2:
            yield points


def load_streamlines(path: Path, key: str | None = None) -> list[np.ndarray]:
    """Compatibility wrapper for callers that explicitly need a list."""
    return list(iter_streamlines(path, key))


VoxelSize = float | tuple[float, float, float]


def _scale_phys(P: np.ndarray, voxel_size: VoxelSize) -> np.ndarray:
    Q = P.astype(np.float64).copy()
    if np.isscalar(voxel_size):
        Q *= float(voxel_size)
    else:
        arr = np.asarray(voxel_size, dtype=float)
        Q *= (arr.ravel()[:3])[None, :]
    return Q


def resample_streamline(P: np.ndarray, spacing: float) -> np.ndarray:
    """Linearly resample a polyline at approximately uniform arc-length spacing."""
    if len(P) < 2 or spacing <= 0:
        return P

    segment_lengths = np.linalg.norm(np.diff(P, axis=0), axis=1)
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    keep = np.concatenate(([True], np.diff(cumulative) > 0))
    cumulative = cumulative[keep]
    P = P[keep]
    if len(P) < 2 or cumulative[-1] <= 0:
        return P

    n_points = max(2, int(np.ceil(cumulative[-1] / spacing)) + 1)
    distances = np.linspace(0.0, cumulative[-1], n_points)
    return np.column_stack(
        [np.interp(distances, cumulative, P[:, axis]) for axis in range(3)]
    )


def streamline_length(sl_xyz: np.ndarray, voxel_size: VoxelSize) -> float:
    P = _scale_phys(sl_xyz, voxel_size)
    diffs = np.diff(P, axis=0)
    return float(np.linalg.norm(diffs, axis=1).sum())


def chord_length(sl_xyz: np.ndarray, voxel_size: VoxelSize) -> float:
    P = _scale_phys(sl_xyz, voxel_size)
    return float(np.linalg.norm(P[-1] - P[0]))


def curvature_discrete(sl_xyz: np.ndarray, voxel_size: VoxelSize) -> np.ndarray:
    """Menger curvature through each consecutive triplet of points."""
    if len(sl_xyz) < 3:
        return np.zeros((0,), dtype=np.float32)
    P = _scale_phys(sl_xyz, voxel_size)
    A, B, C = P[:-2], P[1:-1], P[2:]
    AB, AC = B - A, C - A
    cross = np.cross(AB, AC)
    area2 = np.linalg.norm(cross, axis=1)
    ab = np.linalg.norm(AB, axis=1)
    bc = np.linalg.norm(C - B, axis=1)
    ac = np.linalg.norm(AC, axis=1)
    denom = ab * bc * ac
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa = np.where(denom > 0, 2.0 * area2 / denom, 0.0)
    return np.nan_to_num(kappa).astype(np.float32)


def percentile_bounds(
    values_list: list[np.ndarray], p_lo: float, p_hi: float
) -> tuple[float, float]:
    finite_values = [v[np.isfinite(v)] for v in values_list if v.size > 0]
    if not finite_values:
        return 0.0, 1.0
    vals = np.concatenate(finite_values, axis=0)
    if vals.size == 0:
        return 0.0, 1.0
    lo = float(np.percentile(vals, p_lo))
    hi = float(np.percentile(vals, p_hi))
    if hi <= lo:
        lo, hi = lo - 0.5, hi + 0.5
    return lo, hi


def gaussian_kde_1d(
    samples: np.ndarray, grid: np.ndarray, bandwidth: float | None = None
) -> np.ndarray:
    """Approximate a Gaussian KDE using a fixed-memory smoothed histogram."""
    x = samples[np.isfinite(samples)].astype(np.float64)
    n = len(x)
    if n < 2 or len(grid) < 2:
        return np.zeros_like(grid)
    std = x.std(ddof=1)
    if bandwidth is None:
        h = (
            std * (n ** (-1.0 / 5.0))
            if std > 0
            else (np.ptp(x) / 100.0 if np.ptp(x) > 0 else 1.0)
        )
    else:
        h = float(bandwidth)
    step = float(np.median(np.diff(grid)))
    if not np.isfinite(h) or h <= 0 or not np.isfinite(step) or step <= 0:
        return np.zeros_like(grid)

    edges = np.concatenate((grid - step / 2.0, [grid[-1] + step / 2.0]))
    counts = np.histogram(x, bins=edges)[0].astype(np.float64)
    sigma_bins = h / step
    radius = max(1, int(np.ceil(4.0 * sigma_bins)))
    offsets = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (offsets / sigma_bins) ** 2)
    kernel /= kernel.sum()
    smoothed = np.convolve(
        np.pad(counts, radius, mode="constant"), kernel, mode="same"
    )[radius:-radius]
    return smoothed / (n * step)


def make_grid(lo: float, hi: float, n_points: int = 512) -> np.ndarray:
    if hi <= lo:
        hi = lo + 1.0
    return np.linspace(lo, hi, n_points)


# ---------- plotting ----------


def normalize_curve(y: np.ndarray, do_norm: bool) -> np.ndarray:
    if not do_norm or y.size == 0:
        return y
    m = np.max(y)
    return y / m if m > 0 else y


def plot_hist(ax, data, bins, xlabel, normalize, label):
    counts, _ = np.histogram(data, bins=bins)
    weights = None
    if normalize and counts.size and counts.max() > 0:
        weights = np.full(data.shape, 1.0 / counts.max(), dtype=float)
    ax.hist(
        data,
        bins=bins,
        weights=weights,
        alpha=0.4,
        label=label,
        edgecolor="black",
        linewidth=0.8,
    )

    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel(
        "Frequency" if not normalize else "Normalized (0–1)", fontweight="bold"
    )


def plot_kde(ax, data, grid, xlabel, bandwidth_mult, normalize, label):
    n = len(data)
    if n == 0:
        return
    std = data.std(ddof=1) if n > 1 else 0.0
    h_scott = (
        std * (n ** (-1.0 / 5.0))
        if std > 0
        else (np.ptp(data) / 100.0 if np.ptp(data) > 0 else 1.0)
    )
    bw = h_scott * bandwidth_mult if bandwidth_mult != 1.0 else None
    dens = gaussian_kde_1d(data, grid, bandwidth=bw)
    y = dens
    y = normalize_curve(y, normalize)
    ax.plot(grid, y, label=label)
    ax.fill_between(grid, 0, y, alpha=0.25)
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel("Density" if not normalize else "Normalized (0–1)", fontweight="bold")


# ---------- compute metrics ----------


def compute_metrics(
    path: Path,
    key: str | None,
    voxel_size: VoxelSize,
    min_points: int,
    resample_spacing: float | None = None,
):
    print(f"Load {path}")
    iterator = iter_streamlines(path, key)
    sample = list(itertools.islice(iterator, 512))
    if not sample:
        return {
            "length": np.array([]),
            "mean_curvature": np.array([]),
            "tortuosity": np.array([]),
            "total_streamlines": 0,
            "eligible_streamlines": 0,
            "space_size": np.zeros(3),
            "space_diagonal": 0.0,
        }

    if resample_spacing is None:
        sampled_lengths = []
        for sl in sample:
            P = _scale_phys(sl, voxel_size)
            lengths = np.linalg.norm(np.diff(P, axis=0), axis=1)
            sampled_lengths.append(lengths[np.isfinite(lengths) & (lengths > 0)])
        valid_lengths = [values for values in sampled_lengths if values.size]
        resample_spacing = (
            float(np.median(np.concatenate(valid_lengths))) if valid_lengths else 0.0
        )
    if resample_spacing <= 0:
        raise ValueError("resample_spacing must be positive")
    print(f"Uniform curvature sampling interval: {resample_spacing:g}")

    lengths, mean_curv, torts = [], [], []
    min_vals = np.full(3, np.inf)
    max_vals = np.full(3, -np.inf)
    total_streamlines = 0
    for sl in itertools.chain(sample, iterator):
        P = _scale_phys(sl, voxel_size)
        total_streamlines += 1
        min_vals = np.minimum(min_vals, np.min(P, axis=0))
        max_vals = np.maximum(max_vals, np.max(P, axis=0))
        if sl.shape[0] < min_points:
            continue
        L = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
        D = float(np.linalg.norm(P[-1] - P[0]))
        uniform_points = resample_streamline(P, resample_spacing)
        K = curvature_discrete(uniform_points, 1.0)
        tort = (L / D - 1.0) if D > 0 else np.nan
        lengths.append(L)
        mean_curv.append(K.mean() if K.size else 0.0)
        torts.append(tort)
    size = max_vals - min_vals
    print(
        "Coordinate-space size: "
        f"{size} (x={size[0]:.3g}, y={size[1]:.3g}, z={size[2]:.3g})"
    )
    if (
        path.suffix.lower() == ".trk"
        and np.allclose(np.asarray(voxel_size), 1.0)
        and np.max(size) > 1000
    ):
        print(
            "WARNING: This looks like a legacy TRK stored in voxel coordinates. "
            "QuickBundles thresholds would use voxels, not millimetres. Rerun with "
            "--voxel-size MICROMETRES/1000 (for example, 0.016495 for 16.495 µm).",
            flush=True,
        )
    return dict(
        length=np.array(lengths),
        mean_curvature=np.array(mean_curv),
        tortuosity=np.array(torts),
        total_streamlines=total_streamlines,
        eligible_streamlines=len(lengths),
        space_size=size,
        space_diagonal=float(np.linalg.norm(size)),
    )


def run_quickbundles(
    path: Path,
    key: str | None,
    label: str,
    outdir: Path,
    voxel_size: VoxelSize,
    min_points: int,
    thresholds: list[float],
    cluster_points: int,
    max_streamlines: int,
    seed: int,
    cluster_min_size: int = 1,
    max_clusters: int | None = None,
) -> dict:
    """Cluster a reproducible sample and save its centroids and members."""
    if max_clusters is not None and max_clusters < 1:
        raise ValueError("max_clusters must be at least 1")
    from dipy.segment.clustering import QuickBundles, QuickBundlesX
    from dipy.segment.featurespeed import ResampleFeature
    from dipy.segment.metricspeed import AveragePointwiseEuclideanMetric

    feature = ResampleFeature(nb_points=cluster_points)
    rng = np.random.default_rng(seed)
    reservoir: list[tuple[int, np.ndarray]] = []
    eligible = 0
    limit = None if max_streamlines == 0 else max_streamlines
    method = "QuickBundles" if len(thresholds) == 1 else "QuickBundlesX"

    print(f"\nPreparing {method} for {path}", flush=True)
    print(
        f"  thresholds: {thresholds} | points per streamline: {cluster_points}",
        flush=True,
    )
    if limit is None:
        print("  sampling: all eligible streamlines", flush=True)
    else:
        print(
            f"  sampling: reproducible reservoir of at most {limit:,} "
            f"streamlines (seed={seed})",
            flush=True,
        )

    for source_index, streamline in enumerate(iter_streamlines(path, key)):
        if source_index and source_index % 100_000 == 0:
            print(
                f"  scanned {source_index:,} | eligible {eligible:,} | "
                f"selected {len(reservoir):,}",
                flush=True,
            )
        if len(streamline) < min_points:
            continue
        if not np.all(np.isfinite(streamline)):
            continue

        eligible += 1
        replacement = None
        if limit is not None and len(reservoir) >= limit:
            replacement = int(rng.integers(eligible))
            if replacement >= limit:
                continue

        points = feature.extract(_scale_phys(streamline, voxel_size)).astype(
            np.float32, copy=False
        )
        item = (source_index, points)
        if replacement is None:
            reservoir.append(item)
        else:
            reservoir[replacement] = item

    if not reservoir:
        raise ValueError(f"No valid streamlines available for clustering in {path}")

    # QuickBundles is order-dependent. Restore source order after sampling so
    # that the fixed random seed gives reproducible results.
    reservoir.sort(key=lambda item: item[0])
    source_indices = np.asarray([item[0] for item in reservoir], dtype=np.int64)
    streamlines = [item[1] for item in reservoir]
    print(
        f"Sampling complete: selected {len(streamlines):,} of "
        f"{eligible:,} eligible streamlines",
        flush=True,
    )

    metric = AveragePointwiseEuclideanMetric(feature)
    print(f"Running {method}...", flush=True)
    if len(thresholds) == 1:
        clusters = QuickBundles(thresholds[0], metric=metric).cluster(streamlines)
    else:
        tree = QuickBundlesX(thresholds, metric=metric).cluster(streamlines)
        clusters = tree.get_clusters(len(thresholds))

    ranked_clusters = sorted(clusters, key=len, reverse=True)
    centroids = [
        np.asarray(cluster.centroid, dtype=np.float32) for cluster in ranked_clusters
    ]
    sizes = np.asarray([len(cluster) for cluster in ranked_clusters], dtype=np.int32)
    print(f"Clustering complete: {len(ranked_clusters):,} clusters", flush=True)
    if len(sizes):
        largest = ", ".join(f"{size:,}" for size in sizes[:5])
        print(f"  five largest cluster sizes: {largest}", flush=True)

    eligible_cluster_ids = np.flatnonzero(sizes >= cluster_min_size)
    if not len(eligible_cluster_ids):
        raise ValueError(
            f"No cluster has at least {cluster_min_size} streamlines; "
            f"the largest has {int(sizes.max())}"
        )
    kept_ids = eligible_cluster_ids[:max_clusters]
    saved_cluster_mask = np.zeros(len(sizes), dtype=bool)
    saved_cluster_mask[kept_ids] = True
    kept_centroids = [centroids[cluster_id] for cluster_id in kept_ids]
    kept_sizes = sizes[kept_ids]
    kept_streamlines = int(kept_sizes.sum())
    print(
        f"Keeping {len(kept_ids):,}/{len(eligible_cluster_ids):,} clusters with "
        f"at least {cluster_min_size:,} members; they contain "
        f"{kept_streamlines:,}/{len(streamlines):,} sampled streamlines "
        f"({100.0 * kept_streamlines / len(streamlines):.1f}%)",
        flush=True,
    )
    assignments = np.full(len(streamlines), -1, dtype=np.int32)
    for cluster_id, cluster in enumerate(ranked_clusters):
        assignments[np.asarray(cluster.indices, dtype=np.int64)] = cluster_id
    if np.any(assignments < 0):
        raise RuntimeError("QuickBundles did not assign every sampled streamline")

    kept_member_indices = np.flatnonzero(saved_cluster_mask[assignments])
    member_clusters = {
        int(source_indices[index]): int(assignments[index])
        for index in kept_member_indices
    }

    safe_label = "".join(
        char if char.isalnum() or char in "-_." else "_" for char in label
    ).strip("._") or path.stem
    prefix = outdir / f"{safe_label}_quickbundles"
    csv_path = prefix.with_name(f"{prefix.name}_clusters.csv")
    trk_path = prefix.with_name(f"{prefix.name}_centroids.trk")
    member_trk_path = prefix.with_name(f"{prefix.name}_members.trk")
    membership_path = prefix.with_name(f"{prefix.name}_membership.npz")

    print("Writing clustering results...", flush=True)
    with csv_path.open("w", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=(
                "cluster_id",
                "streamline_count",
                "percentage",
                "centroid_length",
                "method",
                "thresholds",
                "saved_centroid",
                "sampled_streamlines",
                "eligible_streamlines",
            ),
        )
        writer.writeheader()
        for cluster_id, (centroid, size) in enumerate(zip(centroids, sizes)):
            centroid_length = float(
                np.linalg.norm(np.diff(centroid, axis=0), axis=1).sum()
            )
            writer.writerow(
                {
                    "cluster_id": cluster_id,
                    "streamline_count": int(size),
                    "percentage": 100.0 * int(size) / len(streamlines),
                    "centroid_length": centroid_length,
                    "method": method,
                    "thresholds": ";".join(map(str, thresholds)),
                    "sampled_streamlines": len(streamlines),
                    "eligible_streamlines": eligible,
                    "saved_centroid": bool(saved_cluster_mask[cluster_id]),
                }
            )

    np.savez_compressed(
        membership_path,
        source_streamline_index=source_indices,
        cluster_id=assignments,
    )

    import nibabel as nib

    centroid_tractogram = nib.streamlines.Tractogram(
        kept_centroids,
        data_per_point={
            "cluster_id": [
                np.full((len(centroid), 1), cluster_id, dtype=np.float32)
                for cluster_id, centroid in zip(kept_ids, kept_centroids)
            ]
        },
        data_per_streamline={
            "cluster_size": kept_sizes.astype(np.float32)[:, None],
        },
        affine_to_rasmm=np.eye(4),
    )
    header = None
    if path.suffix.lower() == ".trk":
        header = nib.streamlines.load(str(path), lazy_load=True).header.copy()
    nib.streamlines.save(
        nib.streamlines.TrkFile(
            centroid_tractogram,
            header=header.copy() if header is not None else None,
        ),
        str(trk_path),
    )

    from nibabel.streamlines.tractogram import TractogramItem

    def iter_cluster_members():
        written = 0
        for source_index, streamline in enumerate(iter_streamlines(path, key)):
            cluster_id = member_clusters.get(source_index)
            if cluster_id is None:
                continue
            points = _scale_phys(streamline, voxel_size).astype(
                np.float32, copy=False
            )
            written += 1
            if written % 10_000 == 0:
                print(f"  wrote {written:,} cluster members", flush=True)
            yield TractogramItem(
                points,
                {"cluster_size": np.asarray([sizes[cluster_id]], dtype=np.float32)},
                {
                    "cluster_id": np.full(
                        (len(points), 1), cluster_id, dtype=np.float32
                    )
                },
            )

    print(
        f"Writing {len(member_clusters):,} full-resolution cluster members...",
        flush=True,
    )
    member_tractogram = nib.streamlines.LazyTractogram.from_data_func(
        iter_cluster_members
    )
    member_tractogram.affine_to_rasmm = np.eye(4)
    nib.streamlines.save(
        nib.streamlines.TrkFile(
            member_tractogram,
            header=header.copy() if header is not None else None,
        ),
        str(member_trk_path),
    )

    print(f"Cluster summary: {csv_path}", flush=True)
    print(f"Cluster membership: {membership_path}", flush=True)
    print(f"Cluster centroids: {trk_path}", flush=True)
    print(
        f"Cluster members: {member_trk_path} "
        f"({len(member_clusters):,} streamlines)",
        flush=True,
    )
    return {
        "method": method,
        "eligible_streamlines": eligible,
        "sampled_streamlines": len(streamlines),
        "clusters": len(centroids),
        "csv": csv_path,
        "membership": membership_path,
        "centroids": trk_path,
        "members": member_trk_path,
        "saved_clusters": len(kept_centroids),
        "cluster_sizes": sizes,
    }


def write_cohort_report(
    results: list[dict],
    outdir: Path,
    input_paths: list[Path] | None = None,
    cluster_results: list[dict] | None = None,
    max_samples: int = 20_000,
) -> None:
    """Write equal-sampled, heart-level comparison tables and figures."""
    if len(results) < 2:
        return

    print("\nCreating cohort comparison report...", flush=True)
    rng = np.random.default_rng(0)
    prepared = []
    for result in results:
        diagonal = float(result.get("space_diagonal", 0.0))
        if diagonal <= 0:
            raise ValueError(f"Cannot normalize {result['label']}: zero spatial extent")
        values = np.column_stack(
            (
                result["length"] / diagonal,
                result["mean_curvature"],
                result["tortuosity"],
            )
        )
        values = values[np.all(np.isfinite(values), axis=1)]
        if not len(values):
            raise ValueError(f"No finite streamline metrics for {result['label']}")
        prepared.append(values)

    sample_count = min(len(values) for values in prepared)
    if max_samples:
        sample_count = min(sample_count, max_samples)
    samples = []
    for values in prepared:
        indices = rng.choice(len(values), size=sample_count, replace=False)
        samples.append(values[indices])
    print(
        f"  using {sample_count:,} reproducibly sampled streamlines per heart",
        flush=True,
    )

    metric_specs = (
        (0, "normalized_length", "Normalized streamline length"),
        (1, "mean_curvature", "Mean curvature"),
        (2, "tortuosity", "Tortuosity"),
    )
    labels = [str(result["label"]) for result in results]
    paths = input_paths or [Path("")] * len(results)

    summary_path = outdir / "heart_summary.csv"
    summary_rows = []
    for index, (result, path) in enumerate(zip(results, paths)):
        diagonal = float(result["space_diagonal"])
        full_metrics = {
            "length": np.asarray(result["length"], dtype=float),
            "normalized_length": np.asarray(result["length"], dtype=float)
            / diagonal,
            "mean_curvature": np.asarray(result["mean_curvature"], dtype=float),
            "tortuosity": np.asarray(result["tortuosity"], dtype=float),
        }
        row = {
            "heart": result["label"],
            "input": str(path),
            "total_streamlines": result.get(
                "total_streamlines", len(result["length"])
            ),
            "eligible_streamlines": result.get(
                "eligible_streamlines", len(result["length"])
            ),
            "space_x": float(result["space_size"][0]),
            "space_y": float(result["space_size"][1]),
            "space_z": float(result["space_size"][2]),
            "space_diagonal": diagonal,
        }
        for name, values in full_metrics.items():
            values = values[np.isfinite(values)]
            q1, median, q3 = np.percentile(values, [25, 50, 75])
            row.update(
                {
                    f"{name}_mean": float(np.mean(values)),
                    f"{name}_q1": float(q1),
                    f"{name}_median": float(median),
                    f"{name}_q3": float(q3),
                    f"{name}_iqr": float(q3 - q1),
                }
            )

        if cluster_results:
            clustering = cluster_results[index]
            sizes = np.asarray(clustering["cluster_sizes"], dtype=float)
            fractions = sizes / sizes.sum()
            entropy = float(-np.sum(fractions * np.log(fractions)))
            row.update(
                {
                    "quickbundles_method": clustering["method"],
                    "quickbundles_sampled_streamlines": clustering[
                        "sampled_streamlines"
                    ],
                    "quickbundles_clusters": clustering["clusters"],
                    "quickbundles_top10_fraction": float(fractions[:10].sum()),
                    "quickbundles_entropy": entropy,
                    "quickbundles_normalized_entropy": (
                        entropy / np.log(len(sizes)) if len(sizes) > 1 else 0.0
                    ),
                    "quickbundles_effective_clusters": float(np.exp(entropy)),
                }
            )
        summary_rows.append(row)

    summary_path = outdir / "heart_summary.csv"
    with summary_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"  heart-level summary: {summary_path}", flush=True)

    colors = plt.get_cmap("tab10")(np.linspace(0, 1, len(results)))
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for metric_index, _, title in metric_specs:
        ax = axes[metric_index]
        for label, color, values in zip(labels, colors, samples):
            sorted_values = np.sort(values[:, metric_index])
            probability = np.arange(1, len(sorted_values) + 1) / len(sorted_values)
            ax.plot(sorted_values, probability, label=label, color=color, linewidth=2)
        ax.set_xlabel(title, fontweight="bold")
        ax.set_ylabel("Cumulative fraction" if metric_index == 0 else "")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.2)
    axes[-1].legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"geometry_distributions.{suffix}", bbox_inches="tight")
    plt.close(fig)
    print(
        f"  geometry distributions: {outdir / 'geometry_distributions.png'}",
        flush=True,
    )

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), sharex=True)
    x = np.arange(len(results))
    for metric_index, _, title in metric_specs:
        medians, lower, upper = [], [], []
        for values in prepared:
            q1, median, q3 = np.percentile(values[:, metric_index], [25, 50, 75])
            medians.append(median)
            lower.append(median - q1)
            upper.append(q3 - median)
        axes[metric_index].errorbar(
            x,
            medians,
            yerr=np.asarray([lower, upper]),
            fmt="o",
            capsize=4,
            color="black",
        )
        axes[metric_index].scatter(x, medians, c=colors, s=55, zorder=3)
        axes[metric_index].set_title(title, fontweight="bold")
        axes[metric_index].set_xticks(x, labels, rotation=45, ha="right")
        axes[metric_index].set_ylabel("Median and interquartile range")
        axes[metric_index].grid(axis="y", alpha=0.2)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"heart_metric_summary.{suffix}", bbox_inches="tight")
    plt.close(fig)
    print(f"  metric summary: {outdir / 'heart_metric_summary.png'}", flush=True)

    distance = np.zeros((len(results), len(results)), dtype=float)
    metric_distances = []
    for metric_index, _, _ in metric_specs:
        matrix = np.zeros_like(distance)
        pooled = np.concatenate([values[:, metric_index] for values in samples])
        q1, q3 = np.percentile(pooled, [25, 75])
        scale = q3 - q1
        if scale <= 0:
            scale = np.std(pooled) or 1.0
        sorted_values = [np.sort(values[:, metric_index]) for values in samples]
        for left in range(len(results)):
            for right in range(left + 1, len(results)):
                value = float(
                    np.mean(np.abs(sorted_values[left] - sorted_values[right]))
                )
                matrix[left, right] = matrix[right, left] = value / scale
        metric_distances.append(matrix)
    distance = np.mean(metric_distances, axis=0)

    matrix_path = outdir / "heart_distance_matrix.csv"
    with matrix_path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["heart", *labels])
        for label, row in zip(labels, distance):
            writer.writerow([label, *row])

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    image = ax.imshow(distance, cmap="magma", vmin=0)
    ax.set_xticks(np.arange(len(labels)), labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.set_title("Combined streamline-distribution distance", fontweight="bold")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Mean Wasserstein distance / pooled IQR")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"heart_distance_heatmap.{suffix}", bbox_inches="tight")
    plt.close(fig)
    print(f"  distance matrix: {matrix_path}", flush=True)
    print(
        f"  distance heatmap: {outdir / 'heart_distance_heatmap.png'}",
        flush=True,
    )

    if cluster_results:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        for label, color, clustering in zip(labels, colors, cluster_results):
            sizes = np.asarray(clustering["cluster_sizes"], dtype=float)
            fractions = sizes / sizes.sum()
            ranks = np.arange(1, len(sizes) + 1)
            axes[0].plot(ranks, 100 * fractions, label=label, color=color)
            axes[1].plot(
                ranks, 100 * np.cumsum(fractions), label=label, color=color
            )
        axes[0].set_xscale("log")
        axes[0].set_yscale("log")
        axes[0].set_xlabel("Cluster rank")
        axes[0].set_ylabel("Streamlines in cluster (%)")
        axes[0].set_title("Cluster rank-abundance", fontweight="bold")
        axes[1].set_xscale("log")
        axes[1].set_xlabel("Number of largest clusters")
        axes[1].set_ylabel("Cumulative streamlines (%)")
        axes[1].set_ylim(0, 100)
        axes[1].set_title("Cumulative cluster coverage", fontweight="bold")
        for ax in axes:
            ax.grid(alpha=0.2)
        axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            fig.savefig(
                outdir / f"cluster_rank_abundance.{suffix}", bbox_inches="tight"
            )
        plt.close(fig)
        print(
            f"  cluster rank-abundance: {outdir / 'cluster_rank_abundance.png'}",
            flush=True,
        )
    else:
        print("  cluster rank-abundance skipped (use --quickbundles to create it)")


# ---------- CLI ----------


def script():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", type=Path)
    ap.add_argument("--labels", nargs="+", type=str)
    ap.add_argument("--key", type=str, default=None)
    ap.add_argument(
        "--voxel-size",
        type=float,
        default=1.0,
        help=(
            "Coordinate multiplier for legacy TRK/NPY/NPZ inputs. Keep 1 for "
            "TRKs whose header already contains physical spacing"
        ),
    )
    ap.add_argument(
        "--resample-spacing",
        type=float,
        default=None,
        help=(
            "Physical spacing used for curvature. By default, estimate one common "
            "spacing per input from its first 512 streamlines"
        ),
    )
    ap.add_argument(
        "--quickbundles",
        type=float,
        nargs="+",
        metavar="THRESHOLD",
        help=(
            "Cluster streamlines in physical units. One threshold uses QuickBundles; "
            "multiple decreasing thresholds use QuickBundlesX"
        ),
    )
    ap.add_argument(
        "--cluster-points",
        type=int,
        default=12,
        help="Points used by the orientation-invariant clustering metric (default: 12)",
    )
    ap.add_argument(
        "--cluster-max-streamlines",
        type=int,
        default=50_000,
        help=(
            "Maximum reproducible sample used for clustering (default: 50000; "
            "use 0 for all streamlines)"
        ),
    )
    ap.add_argument(
        "--cluster-min-size",
        type=int,
        default=1,
        help=(
            "Write centroid and member streamlines only for clusters with at "
            "least this many sampled members"
        ),
    )
    ap.add_argument(
        "--max-clusters",
        type=int,
        default=None,
        help=(
            "Write only the N largest clusters that pass --cluster-min-size; "
            "the CSV and membership NPZ still contain all clusters"
        ),
    )
    ap.add_argument(
        "--cluster-seed",
        type=int,
        default=0,
        help="Random seed for reservoir sampling (default: 0)",
    )
    ap.add_argument("--min-points", type=int, default=2)
    ap.add_argument("--bins", type=int, default=60)
    ap.add_argument("--kde-points", type=int, default=512)
    ap.add_argument("--bandwidth-mult", type=float, default=1.0)
    ap.add_argument("--clip", type=float, nargs=2, default=(0.1, 99.9))
    ap.add_argument("--kde", action="store_true")
    ap.add_argument(
        "--normalize", action="store_true", help="Rescale y-axis between 0 and 1"
    )
    ap.add_argument("--xscale", type=str, choices=["linear", "log"], default="linear")
    ap.add_argument("--yscale", type=str, choices=["linear", "log"], default="linear")
    ap.add_argument("--outdir", type=Path, default=Path("./compare_out"))
    args = ap.parse_args()

    labels = args.labels or [p.stem for p in args.inputs]
    if len(labels) != len(args.inputs):
        ap.error("--labels must contain exactly one label per input")
    if len(set(labels)) != len(labels):
        ap.error("--labels must be unique for cohort comparison")
    if args.voxel_size <= 0:
        ap.error("--voxel-size must be positive")
    if args.min_points < 2:
        ap.error("--min-points must be at least 2")
    if args.resample_spacing is not None and args.resample_spacing <= 0:
        ap.error("--resample-spacing must be positive")
    if args.cluster_points < 2:
        ap.error("--cluster-points must be at least 2")
    if args.cluster_max_streamlines < 0:
        ap.error("--cluster-max-streamlines cannot be negative")
    if args.cluster_min_size < 1:
        ap.error("--cluster-min-size must be at least 1")
    if args.max_clusters is not None and args.max_clusters < 1:
        ap.error("--max-clusters must be at least 1")
    if args.cluster_seed < 0:
        ap.error("--cluster-seed cannot be negative")
    if args.quickbundles:
        if any(threshold <= 0 for threshold in args.quickbundles):
            ap.error("--quickbundles thresholds must be positive")
        if len(args.quickbundles) > 1 and any(
            coarse <= fine
            for coarse, fine in zip(args.quickbundles, args.quickbundles[1:])
        ):
            ap.error("QuickBundlesX thresholds must be strictly decreasing")

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("\nComparing:")
    for p, lab in zip(args.inputs, labels):
        print(f"  - {lab}: {p}")
    print("")

    results = []
    cluster_results = []
    for path, label in zip(args.inputs, labels):
        metrics = compute_metrics(
            path,
            args.key,
            args.voxel_size,
            args.min_points,
            args.resample_spacing,
        )
        metrics["label"] = label
        results.append(metrics)
        if args.quickbundles:
            cluster_results.append(
                run_quickbundles(
                    path=path,
                    key=args.key,
                    label=label,
                    outdir=args.outdir,
                    voxel_size=args.voxel_size,
                    min_points=args.min_points,
                    thresholds=args.quickbundles,
                    cluster_points=args.cluster_points,
                    max_streamlines=args.cluster_max_streamlines,
                    seed=args.cluster_seed,
                    cluster_min_size=args.cluster_min_size,
                    max_clusters=args.max_clusters,
                )
            )

    write_cohort_report(
        results,
        args.outdir,
        input_paths=args.inputs,
        cluster_results=cluster_results or None,
    )

    # Percentile bounds
    arr_len = [m["length"] for m in results if m["length"].size > 0]
    arr_curv = [m["mean_curvature"] for m in results if m["mean_curvature"].size > 0]
    arr_tort = [m["tortuosity"] for m in results if m["tortuosity"].size > 0]

    lo_len, hi_len = percentile_bounds(arr_len, *args.clip)
    lo_curv, hi_curv = percentile_bounds(arr_curv, *args.clip)
    lo_tort, hi_tort = percentile_bounds(arr_tort, *args.clip)

    grid_len = make_grid(lo_len, hi_len, args.kde_points)
    grid_curv = make_grid(lo_curv, hi_curv, args.kde_points)
    grid_tort = make_grid(lo_tort, hi_tort, args.kde_points)

    bins_len = np.linspace(lo_len, hi_len, args.bins + 1)
    bins_curv = np.linspace(lo_curv, hi_curv, args.bins + 1)
    bins_tort = np.linspace(lo_tort, hi_tort, args.bins + 1)

    # --- plotting ---
    for metric, grid, bins, xlabel in [
        ("length", grid_len, bins_len, "Streamline length"),
        ("mean_curvature", grid_curv, bins_curv, "Mean curvature"),
        ("tortuosity", grid_tort, bins_tort, "Tortuosity"),
    ]:
        fig, ax = plt.subplots(figsize=(6, 4))
        for m in results:
            data = m[metric]
            data = data[np.isfinite(data)]
            label = m["label"]
            if args.kde:
                plot_kde(
                    ax, data, grid, xlabel, args.bandwidth_mult, args.normalize, label
                )
            else:
                plot_hist(ax, data, bins, xlabel, args.normalize, label)

        ax.set_xscale(args.xscale)
        ax.set_yscale(args.yscale)

        # Force non-negative lower limits if linear
        if args.xscale == "linear":
            ax.set_xlim(left=0)
        if args.yscale == "linear":
            ax.set_ylim(bottom=0)

        ax.legend()
        fig.tight_layout()
        fig.savefig(args.outdir / f"{metric}.png", bbox_inches="tight")
        fig.savefig(args.outdir / f"{metric}.pdf", bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    script()
