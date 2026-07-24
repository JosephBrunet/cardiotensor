import csv
import importlib
import sys
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

analysis = importlib.import_module("cardiotensor.scripts.analysis_streamlines")
generator_cli = importlib.import_module("cardiotensor.scripts.generate_streamlines")


def test_menger_curvature_has_correct_scale():
    semicircle = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0]], dtype=float)

    curvature = analysis.curvature_discrete(semicircle, 1.0)

    np.testing.assert_allclose(curvature, [1.0])


def test_uniform_arc_length_resampling_handles_duplicate_points():
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 0, 0], [3, 0, 0]], dtype=float)

    resampled = analysis.resample_streamline(points, spacing=1.0)

    np.testing.assert_allclose(resampled[:, 0], [0, 1, 2, 3])
    np.testing.assert_allclose(np.diff(resampled, axis=0)[:, 0], 1.0)


def test_metric_computation_streams_and_uses_physical_scale(tmp_path):
    streamlines = np.empty(2, dtype=object)
    streamlines[0] = np.array([[0, 0, 0], [0.25, 0, 0], [3, 0, 0]])
    streamlines[1] = np.array([[0, 0, 0], [0, 2, 0]])
    path = tmp_path / "streamlines.npy"
    np.save(path, streamlines, allow_pickle=True)

    metrics = analysis.compute_metrics(
        path,
        key=None,
        voxel_size=2.0,
        min_points=2,
        resample_spacing=1.0,
    )

    np.testing.assert_allclose(metrics["length"], [6.0, 4.0])
    np.testing.assert_allclose(metrics["mean_curvature"], 0.0)
    np.testing.assert_allclose(metrics["tortuosity"], 0.0)
    np.testing.assert_allclose(metrics["space_size"], [6.0, 4.0, 0.0])
    assert np.isclose(metrics["space_diagonal"], np.sqrt(52))
    assert metrics["total_streamlines"] == 2
    assert metrics["eligible_streamlines"] == 2


def test_cohort_report_writes_tables_and_figures(tmp_path):
    results = [
        {
            "label": "heart_a",
            "length": np.array([2.0, 4.0, 6.0, 8.0]),
            "mean_curvature": np.array([0.1, 0.2, 0.3, 0.4]),
            "tortuosity": np.array([0.0, 0.1, 0.2, 0.3]),
            "total_streamlines": 5,
            "eligible_streamlines": 4,
            "space_size": np.array([6.0, 8.0, 0.0]),
            "space_diagonal": 10.0,
        },
        {
            "label": "heart_b",
            "length": np.array([3.0, 6.0, 9.0, 12.0]),
            "mean_curvature": np.array([0.2, 0.3, 0.4, 0.5]),
            "tortuosity": np.array([0.1, 0.2, 0.3, 0.4]),
            "total_streamlines": 4,
            "eligible_streamlines": 4,
            "space_size": np.array([0.0, 0.0, 12.0]),
            "space_diagonal": 12.0,
        },
    ]
    clusters = [
        {
            "method": "QuickBundles",
            "sampled_streamlines": 4,
            "clusters": 2,
            "cluster_sizes": np.array([3, 1]),
        },
        {
            "method": "QuickBundles",
            "sampled_streamlines": 4,
            "clusters": 2,
            "cluster_sizes": np.array([2, 2]),
        },
    ]

    analysis.write_cohort_report(
        results,
        tmp_path,
        input_paths=[tmp_path / "a.trk", tmp_path / "b.trk"],
        cluster_results=clusters,
        max_samples=3,
    )

    expected = [
        "heart_summary.csv",
        "heart_distance_matrix.csv",
        "geometry_distributions.png",
        "geometry_distributions.pdf",
        "heart_metric_summary.png",
        "heart_metric_summary.pdf",
        "heart_distance_heatmap.png",
        "heart_distance_heatmap.pdf",
        "cluster_rank_abundance.png",
        "cluster_rank_abundance.pdf",
    ]
    assert all((tmp_path / name).exists() for name in expected)

    with (tmp_path / "heart_summary.csv").open(newline="") as input_file:
        rows = list(csv.DictReader(input_file))
    assert [row["heart"] for row in rows] == ["heart_a", "heart_b"]
    assert np.isclose(float(rows[0]["normalized_length_median"]), 0.5)
    assert np.isclose(float(rows[0]["quickbundles_top10_fraction"]), 1.0)

    with (tmp_path / "heart_distance_matrix.csv").open(newline="") as input_file:
        matrix_rows = list(csv.reader(input_file))
    matrix = np.asarray([row[1:] for row in matrix_rows[1:]], dtype=float)
    np.testing.assert_allclose(matrix, matrix.T)
    np.testing.assert_allclose(np.diag(matrix), 0.0)
    assert matrix[0, 1] > 0


def test_trk_analysis_requests_lazy_geometry_only(monkeypatch, tmp_path):
    calls = []
    tractogram = SimpleNamespace(
        streamlines=iter([np.array([[0, 0, 0], [1, 0, 0]], dtype=float)])
    )

    def fake_load(path, **kwargs):
        calls.append((path, kwargs))
        return SimpleNamespace(tractogram=tractogram)

    monkeypatch.setattr("nibabel.streamlines.load", fake_load)
    path = tmp_path / "streamlines.trk"

    streamlines = list(analysis.iter_streamlines(path))

    assert len(streamlines) == 1
    assert calls == [(str(path), {"lazy_load": True})]


def test_fixed_memory_kde_and_normalized_histogram():
    rng = np.random.default_rng(123)
    samples = rng.normal(size=100_000)
    grid = np.linspace(-5, 5, 512)

    density = analysis.gaussian_kde_1d(samples, grid)

    assert density.shape == grid.shape
    assert np.all(np.isfinite(density))
    assert 0.98 < np.trapezoid(density, grid) < 1.02

    fig, ax = plt.subplots()
    analysis.plot_hist(
        ax,
        samples,
        bins=np.linspace(-5, 5, 51),
        xlabel="value",
        normalize=True,
        label="sample",
    )
    heights = [patch.get_height() for patch in ax.patches]
    assert np.isclose(max(heights), 1.0)
    plt.close(fig)


def test_quickbundlesx_saves_clusters_centroids_and_membership(tmp_path):
    streamlines = np.empty(5, dtype=object)
    streamlines[0] = np.array([[0, 0, 0], [5, 0, 0], [10, 0, 0]])
    streamlines[1] = np.array([[0, 0.2, 0], [5, 0.2, 0], [10, 0.2, 0]])
    streamlines[2] = np.array([[10, 0.1, 0], [5, 0.1, 0], [0, 0.1, 0]])
    streamlines[3] = np.array([[0, 10, 0], [5, 10, 0], [10, 10, 0]])
    streamlines[4] = np.array([[0, 10.2, 0], [5, 10.2, 0], [10, 10.2, 0]])
    path = tmp_path / "streamlines.npy"
    np.save(path, streamlines, allow_pickle=True)

    result = analysis.run_quickbundles(
        path=path,
        key=None,
        label="heart",
        outdir=tmp_path,
        voxel_size=1.0,
        min_points=2,
        thresholds=[5.0, 1.0],
        cluster_points=8,
        max_streamlines=0,
        seed=0,
        cluster_min_size=3,
    )

    assert result["method"] == "QuickBundlesX"
    assert result["eligible_streamlines"] == 5
    assert result["sampled_streamlines"] == 5
    assert result["clusters"] == 2
    assert result["saved_clusters"] == 1
    assert result["csv"].exists()
    with result["csv"].open(newline="") as input_file:
        rows = list(csv.DictReader(input_file))
    assert len(rows) == 2
    assert [row["saved_centroid"] for row in rows] == ["True", "False"]

    membership = np.load(result["membership"])
    np.testing.assert_array_equal(
        membership["source_streamline_index"], np.arange(5)
    )
    np.testing.assert_array_equal(
        np.sort(np.bincount(membership["cluster_id"])), [2, 3]
    )

    import nibabel as nib

    centroids = nib.streamlines.load(str(result["centroids"])).tractogram
    assert len(centroids.streamlines) == 1
    np.testing.assert_array_equal(
        centroids.data_per_streamline["cluster_size"].ravel(), [3]
    )
    for cluster_id, point_values in enumerate(
        centroids.data_per_point["cluster_id"]
    ):
        assert np.all(point_values == cluster_id)

    members = nib.streamlines.load(str(result["members"])).tractogram
    assert len(members.streamlines) == 3
    assert [len(streamline) for streamline in members.streamlines] == [3, 3, 3]
    np.testing.assert_array_equal(
        members.data_per_streamline["cluster_size"].ravel(), [3, 3, 3]
    )
    for point_values in members.data_per_point["cluster_id"]:
        assert np.all(point_values == 0)

    quickbundles = analysis.run_quickbundles(
        path=path,
        key=None,
        label="heart_qb",
        outdir=tmp_path,
        voxel_size=1.0,
        min_points=2,
        thresholds=[1.0],
        cluster_points=8,
        max_streamlines=0,
        seed=0,
        max_clusters=1,
    )

    assert quickbundles["method"] == "QuickBundles"
    assert quickbundles["clusters"] == 2
    assert quickbundles["saved_clusters"] == 1
    assert len(
        nib.streamlines.load(str(quickbundles["members"])).tractogram.streamlines
    ) == 3
    with quickbundles["csv"].open(newline="") as input_file:
        rows = list(csv.DictReader(input_file))
    assert [row["saved_centroid"] for row in rows] == ["True", "False"]


def test_generate_streamlines_writes_configured_spacing_in_mm(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(
        generator_cli,
        "read_conf_file",
        lambda path: {
            "OUTPUT_PATH": str(tmp_path),
            "MASK_PATH": None,
            "VECTOR_FORMAT": "zarr",
            "ANGLE_MODE": "ha_ia",
            "VOXEL_SIZE": 16.495,
        },
    )
    monkeypatch.setattr(
        generator_cli, "vector_field_path", lambda output, storage_format: output
    )
    monkeypatch.setattr(
        generator_cli,
        "generate_streamlines_from_params",
        lambda **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["cardio-generate-streamlines", str(tmp_path / "parameters.conf")],
    )

    generator_cli.script()

    np.testing.assert_allclose(captured["voxel_sizes_zyx"], [0.016495] * 3)
