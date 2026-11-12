import inspect
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
import zarr

import cardiotensor.orientation.orientation_computation_functions as orientation_functions
import cardiotensor.orientation.orientation_computation_pipeline as orientation_pipeline
from cardiotensor.orientation.orientation_computation_functions import (
    adjust_start_end_index,
    calculate_structure_tensor,
    compute_azimuth_and_elevation,
    compute_fraction_anisotropy,
    compute_helical_and_intrusion_angles,
    interpolate_points,
    orient_vectors_z_positive,
    plot_images,
    remove_padding,
    rotate_vectors_to_new_axis,
    write_images,
)
from cardiotensor.orientation.orientation_computation_pipeline import (
    _ensure_structure_tensor_memory,
    _memory_aware_worker_count,
    _safe_low_memory_chunk_size,
    check_already_processed,
    compute_orientation,
)
from cardiotensor.tractography.generate_streamlines import _select_seed_points
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.image_io import (
    ZARR_VECTOR_FORMAT_VERSION,
    ZarrVectorFieldStore,
    initialize_zarr_vector_field,
    open_zarr_vector_field,
    write_vector_field,
)


def test_interpolate_points_linear():
    p1 = (0, 0, 0)
    p2 = (10, 10, 10)

    # Interpolate 5 points between p1 and p2
    points = interpolate_points([p1, p2], N_img=20)

    # Validate shape (20 points, 3 coordinates)
    assert isinstance(points, np.ndarray)
    assert points.shape == (20, 3)

    # Ensure points are evenly spaced
    diffs = np.diff(points, axis=0)
    distances = np.linalg.norm(diffs, axis=1)
    assert np.allclose(distances, distances[0])


def test_adjust_start_end_index():
    # Case 1: Normal mode, no padding
    start, end = adjust_start_end_index(0, 5, 10, 0, 0, False, 0)
    assert start == 0
    assert end == 5

    # Case 2: Normal mode, with padding
    start, end = adjust_start_end_index(0, 5, 10, 2, 2, False, 0)
    assert start == 0  # cannot go below 0
    assert end == 7  # 5 + 2

    # Case 3: Normal mode, with end clamping
    start, end = adjust_start_end_index(8, 9, 10, 1, 5, False, 0)
    assert start == 7  # 8 - 1
    assert end == 10  # clamped to N_img

    # Case 4: Test mode, centered on n_slice with padding
    start, end = adjust_start_end_index(0, 0, 10, 1, 1, True, 5)
    assert start == 4  # 5 - 1
    assert end == 7  # 5 + 1 + 1


def test_structure_tensor_and_fa():
    """
    Test that structure tensor eigen decomposition and FA computation produce
    correct shapes and valid values.
    """
    # Generate a small random 3D volume
    volume = np.random.rand(5, 5, 5).astype(np.float32)

    # Compute eigenvalues and eigenvectors of the structure tensor
    val, vec = calculate_structure_tensor(volume, sigma=1.0, rho=1.0)

    # vec shape should be (3, Z, Y, X)
    assert isinstance(vec, np.ndarray)
    assert vec.shape[0] == 3, "Eigenvectors should have 3 components (v1,v2,v3)"
    assert vec.shape == (3,) + volume.shape

    # Compute Fractional Anisotropy (FA) from eigenvalues
    fa_map = compute_fraction_anisotropy(val[:, 0, :, :])
    assert isinstance(fa_map, np.ndarray)
    assert fa_map.shape == volume.shape[1:]
    assert np.all(np.isfinite(fa_map)), "FA map contains NaN or Inf"
    assert np.all(fa_map >= 0) and np.all(fa_map <= 1), "FA values out of range [0,1]"


def test_structure_tensor_can_skip_eigenvalue_allocation(monkeypatch):
    captured = {}

    def fake_parallel(volume, sigma, rho, **kwargs):
        captured.update(kwargs)
        vectors = np.zeros((3, *volume.shape), dtype=np.float32)
        return None, None, vectors

    monkeypatch.setattr(
        orientation_functions, "parallel_structure_tensor_analysis", fake_parallel
    )
    values, vectors = calculate_structure_tensor(
        np.ones((2, 3, 4), dtype=np.uint16),
        sigma=1.0,
        rho=1.0,
        devices=["cpu"],
        return_eigenvalues=False,
    )

    assert values is None
    assert captured["eigenvalues"] is None
    assert vectors.shape == (3, 2, 3, 4)


def test_structure_tensor_can_use_memory_mapped_outputs(tmp_path: Path, monkeypatch):
    captured = {}

    def fake_parallel(volume, sigma, rho, **kwargs):
        captured.update(kwargs)
        kwargs["eigenvectors"][:] = 1
        kwargs["eigenvalues"][:] = 2
        kwargs["progress_callback_fn"](16, 16)
        return None, kwargs["eigenvalues"], kwargs["eigenvectors"]

    monkeypatch.setattr(
        orientation_functions, "parallel_structure_tensor_analysis", fake_parallel
    )
    values, vectors = calculate_structure_tensor(
        np.ones((2, 3, 4), dtype=np.uint8),
        sigma=1.0,
        rho=1.0,
        devices=["cpu"],
        memmap_dir=tmp_path,
    )

    assert isinstance(values, np.memmap)
    assert isinstance(vectors, np.memmap)
    np.testing.assert_array_equal(values, 2)
    np.testing.assert_array_equal(vectors, 1)


def test_vector_only_pipeline_preserves_input_dtype_and_cleans_memmaps(
    tmp_path: Path, monkeypatch
):
    volume_dir = tmp_path / "volume"
    volume_dir.mkdir()
    for z in range(3):
        np.save(
            volume_dir / f"slice_{z:03d}.npy",
            np.full((4, 5), z + 1, dtype=np.uint16),
        )

    captured = {}

    def fake_structure_tensor(volume, sigma, rho, **kwargs):
        captured["input_dtype"] = volume.dtype
        captured.update(kwargs)
        scratch = Path(kwargs["memmap_dir"])
        captured["scratch"] = scratch
        vectors = np.lib.format.open_memmap(
            scratch / "eigenvectors.npy",
            mode="w+",
            dtype=np.float32,
            shape=(3, *volume.shape),
        )
        vectors[:] = 0
        vectors[0] = 1
        return None, vectors

    monkeypatch.setattr(
        orientation_pipeline, "calculate_structure_tensor", fake_structure_tensor
    )
    output_dir = tmp_path / "output"
    compute_orientation(
        volume_path=str(volume_dir),
        output_dir=str(output_dir),
        sigma=0.1,
        rho=0.1,
        axis_points=np.array([(2, 2, 0), (2, 2, 2)]),
        write_vectors=True,
        write_angles=False,
        use_gpu=False,
        low_memory=True,
    )

    assert captured["input_dtype"] == np.uint16
    assert captured["return_eigenvalues"] is False
    assert captured["scratch"].parent == (
        output_dir / ".cardiotensor_scratch"
    ).resolve()
    assert not captured["scratch"].exists()
    vector_store = open_zarr_vector_field(output_dir)
    assert vector_store.vectors.shape == (3, 3, 4, 5)
    assert np.all(vector_store.completed[:])


def test_vector_rotation_and_angles():
    """
    Test rotation of a 2D vector field slice (3, Y, X) to a new axis.
    - Rotating a field aligned with the target axis should leave it unchanged.
    - Rotating to a different axis should produce rotated vectors.
    """
    Y, X = 3, 3
    vector_slice = np.zeros((3, Y, X), dtype=np.float32)

    # Create vectors all pointing along +Z
    vector_slice[2, :, :] = 1.0  # all vectors = (0, 0, 1)

    # --- Test 1: Rotate to the same axis (Z-axis) ---
    rotated_same = rotate_vectors_to_new_axis(
        vector_slice, np.array([0, 0, 1], dtype=np.float32)
    )
    assert isinstance(rotated_same, np.ndarray)
    assert rotated_same.shape == vector_slice.shape
    assert np.allclose(
        rotated_same, vector_slice, atol=1e-6
    ), "Rotation to same axis should not change vectors"
    assert rotated_same.dtype == np.float32


def test_vector_rotation_preserves_input_and_vector_magnitude():
    vector_slice = np.zeros((3, 1, 1), dtype=np.float32)
    vector_slice[:, 0, 0] = [2.0, 0.0, 0.0]
    original = vector_slice.copy()

    rotated = rotate_vectors_to_new_axis(
        vector_slice, np.array([1.0, 0.0, 0.0], dtype=np.float32)
    )

    np.testing.assert_array_equal(vector_slice, original)
    np.testing.assert_allclose(np.linalg.norm(rotated[:, 0, 0]), 2.0, atol=1e-6)


def test_orient_vectors_z_positive_flips_negative_z_vectors():
    vector_slice = np.zeros((3, 2, 2), dtype=np.float32)
    vector_slice[0] = [[1, 1], [1, 1]]
    vector_slice[2] = [[-1, 1], [-0.5, 0.5]]

    oriented = orient_vectors_z_positive(vector_slice)

    assert np.all(oriented[2] >= 0)
    assert np.allclose(oriented[:, 0, 0], [-1, 0, 1])
    assert np.allclose(oriented[:, 0, 1], [1, 0, 1])


def test_compute_azimuth_and_elevation_uses_z_positive_convention():
    vector_slice = np.zeros((3, 1, 2), dtype=np.float32)
    vector_slice[:, 0, 0] = [1, 0, 1]
    vector_slice[:, 0, 1] = [-1, 0, -1]

    azimuth, elevation = compute_azimuth_and_elevation(vector_slice)

    assert np.allclose(elevation, [[45, 45]], atol=1e-6)
    assert np.allclose(azimuth, [[0, 0]], atol=1e-6)


def test_compute_helical_and_intrusion_angles_orients_circumferential_positive():
    vector_slice = np.zeros((3, 1, 2), dtype=np.float32)
    vector_slice[:, 0, 0] = [0, 1, 1]
    vector_slice[:, 0, 1] = [0, -1, -1]

    helical, intrusion = compute_helical_and_intrusion_angles(
        vector_slice, center_point=(0, 0, 0)
    )

    np.testing.assert_allclose(helical, [[45, 45]], atol=1e-6)
    np.testing.assert_allclose(intrusion, [[0, 0]], atol=1e-6)


def test_compute_angles_use_full_vector_without_projection():
    vector_slice = np.zeros((3, 1, 1), dtype=np.float32)
    vector_slice[:, 0, 0] = [1, 2, 3]

    helical, intrusion = compute_helical_and_intrusion_angles(
        vector_slice, center_point=(0, 0, 0)
    )

    expected_helical = np.rad2deg(np.arctan2(3, np.hypot(1, 2)))
    expected_intrusion = np.rad2deg(np.arctan2(1, np.hypot(2, 3)))
    np.testing.assert_allclose(helical, [[expected_helical]], atol=1e-6)
    np.testing.assert_allclose(intrusion, [[expected_intrusion]], atol=1e-6)
    assert helical.dtype == np.float32
    assert intrusion.dtype == np.float32


def test_compute_angles_can_return_projected_intrusion_values():
    vector_slice = np.zeros((3, 1, 1), dtype=np.float32)
    vector_slice[:, 0, 0] = [1, 2, 3]

    helical, intrusion = compute_helical_and_intrusion_angles(
        vector_slice, center_point=(0, 0, 0), projected=True
    )

    expected_helical = np.rad2deg(np.arctan2(3, 2))
    expected_intrusion = np.rad2deg(np.arctan2(1, 2))
    np.testing.assert_allclose(helical, [[expected_helical]], atol=1e-6)
    np.testing.assert_allclose(intrusion, [[expected_intrusion]], atol=1e-6)


def test_write_images_and_vectors(tmp_path: Path):
    # --- Prepare dummy 2D slices ---
    img_helical = np.ones((5, 5), dtype=np.float32)
    img_intrusion = np.ones((5, 5), dtype=np.float32)
    img_FA = np.ones((5, 5), dtype=np.float32)

    # --- Test write_images ---
    write_images(
        img_helical,
        img_intrusion,
        img_FA,
        start_index=0,
        output_dir=str(tmp_path),
        output_format="tif",
        output_type="8bit",
        z=0,
    )

    # Expect 3 TIFF files in HA, IA, FA folders
    ha_files = list((tmp_path / "HA").glob("*.tif"))
    ia_files = list((tmp_path / "IA").glob("*.tif"))
    fa_files = list((tmp_path / "FA").glob("*.tif"))

    assert ha_files, "HA images were not created"
    assert ia_files, "IA images were not created"
    assert fa_files, "FA images were not created"

    # --- Test write_vector_field ---
    vec_field_slice = np.ones((3, 5, 5), dtype=np.float32)  # shape (3, Y, X)
    write_vector_field(
        vec_field_slice,
        start_index=0,
        output_dir=str(tmp_path),
        slice_idx=0,
    )

    eigen_vec_files = list((tmp_path / "eigen_vec").glob("*.npy"))
    assert eigen_vec_files, "Vector field .npy file was not created"


def test_check_already_processed_prints_missing_file(tmp_path: Path, capsys):
    assert not check_already_processed(
        str(tmp_path),
        0,
        1,
        write_vectors=True,
        write_angles=False,
        output_format="tif",
        vector_format="npy",
    )

    output = capsys.readouterr().out
    assert "Missing output file(s) for slice 000000" in output
    assert "eigen_vec_000000.npy" in output


def test_check_already_processed_prints_corrupted_file(tmp_path: Path, capsys):
    vector_dir = tmp_path / "eigen_vec"
    vector_dir.mkdir()
    vector_path = vector_dir / "eigen_vec_000000.npy"
    vector_path.write_bytes(b"invalid")

    with pytest.warns(UserWarning, match="potentially corrupted"):
        assert not check_already_processed(
            str(tmp_path),
            0,
            1,
            write_vectors=True,
            write_angles=False,
            output_format="tif",
            vector_format="npy",
        )

    output = capsys.readouterr().out
    assert "Corrupted output file(s) found for slice 000000" in output
    assert str(vector_path) in output


def test_write_zarr_vector_field_roundtrip_and_restart_marker(tmp_path: Path):
    store = initialize_zarr_vector_field(tmp_path, (3, 5, 7))
    vector_slice = np.ones((3, 5, 7), dtype=np.float32)
    vector_slice[:, :2, :2] = 0

    assert not check_already_processed(
        str(tmp_path),
        1,
        2,
        write_vectors=True,
        write_angles=False,
        output_format="tif",
        vector_format="zarr",
        zarr_store=store,
    )

    write_vector_field(
        vector_slice,
        start_index=1,
        output_dir=str(tmp_path),
        slice_idx=0,
        vector_format="zarr",
        zarr_store=store,
    )

    assert store.vectors.chunks == (3, 1, 5, 7)
    assert store.vectors.shards == (3, 1, 5, 7)
    assert "BloscCodec" in type(store.vectors.compressors[0]).__name__
    assert bool(store.completed[1])
    assert check_already_processed(
        str(tmp_path),
        1,
        2,
        write_vectors=True,
        write_angles=False,
        output_format="tif",
        vector_format="zarr",
        zarr_store=store,
    )

    reader = DataReader(tmp_path / "eigen_vec.zarr")
    loaded = reader.load_volume(start_index=1, end_index=2)
    assert reader.shape == (3, 3, 5, 7)
    np.testing.assert_array_equal(loaded[:, 0], vector_slice)


def test_zarr_vector_slices_can_be_written_concurrently(tmp_path: Path):
    initialize_zarr_vector_field(tmp_path, (4, 8, 9))

    def write_slice(z: int) -> None:
        local_store = open_zarr_vector_field(tmp_path, mode="r+")
        data = np.full((3, 8, 9), z + 1, dtype=np.float32)
        local_store.write_slice(data, z)

    with ThreadPoolExecutor(max_workers=4) as executor:
        list(executor.map(write_slice, range(4)))

    store = open_zarr_vector_field(tmp_path)
    np.testing.assert_array_equal(store.completed[:], np.ones(4, dtype=bool))
    for z in range(4):
        np.testing.assert_array_equal(store.vectors[:, z], z + 1)


def test_interrupted_zarr_write_does_not_set_completion_marker(tmp_path: Path):
    store = initialize_zarr_vector_field(tmp_path, (2, 5, 7))
    vector_slice = np.ones((3, 5, 7), dtype=np.float32)

    class InterruptedVectors:
        shape = store.vectors.shape

        def __setitem__(self, key, value):
            z = key[1]
            store.vectors[:, z, :2, :] = value[:, :2, :]
            raise OSError("simulated interrupted write")

    interrupted_store = ZarrVectorFieldStore(
        store.path, InterruptedVectors(), store.completed
    )
    with pytest.raises(OSError, match="interrupted"):
        interrupted_store.write_slice(vector_slice, 0)

    assert not bool(store.completed[0])
    store.write_slice(vector_slice, 0)
    assert bool(store.completed[0])
    np.testing.assert_array_equal(store.vectors[:, 0], vector_slice)


def test_existing_zarr_store_rejects_incompatible_shards(tmp_path: Path):
    group = zarr.open_group(
        store=str(tmp_path / "eigen_vec.zarr"), mode="w", zarr_format=3
    )
    group.create_array(
        "vectors",
        shape=(3, 2, 1024, 1024),
        chunks=(3, 1, 512, 512),
        shards=(3, 1, 512, 512),
        dtype="float32",
        fill_value=0.0,
        compressors=zarr.codecs.BloscCodec(
            cname="zstd",
            clevel=3,
            shuffle=zarr.codecs.BloscShuffle.bitshuffle,
        ),
        dimension_names=("component", "z", "y", "x"),
    )
    group.create_array(
        "completed", shape=(2,), chunks=(1,), dtype="bool", fill_value=False
    )
    group.attrs.update(
        {
            "cardiotensor_format": "vector_field",
            "cardiotensor_format_version": ZARR_VECTOR_FORMAT_VERSION,
            "axis_order": ["component", "z", "y", "x"],
            "components": ["x", "y", "z"],
            "masked_fill_value": 0.0,
        }
    )

    with pytest.raises(ValueError, match="shards"):
        initialize_zarr_vector_field(tmp_path, (2, 1024, 1024))


def test_open_zarr_store_rejects_incompatible_compressor(tmp_path: Path):
    group = zarr.open_group(
        store=str(tmp_path / "eigen_vec.zarr"), mode="w", zarr_format=3
    )
    group.create_array(
        "vectors",
        shape=(3, 2, 4, 4),
        chunks=(3, 1, 4, 4),
        shards=(3, 1, 4, 4),
        dtype="float32",
        fill_value=0.0,
        compressors=zarr.codecs.ZstdCodec(level=3),
        dimension_names=("component", "z", "y", "x"),
    )
    group.create_array(
        "completed", shape=(2,), chunks=(1,), dtype="bool", fill_value=False
    )
    group.attrs.update(
        {
            "cardiotensor_format": "vector_field",
            "cardiotensor_format_version": ZARR_VECTOR_FORMAT_VERSION,
            "axis_order": ["component", "z", "y", "x"],
            "components": ["x", "y", "z"],
            "masked_fill_value": 0.0,
        }
    )

    with pytest.raises(ValueError, match="Blosc"):
        open_zarr_vector_field(tmp_path)


def test_memory_aware_workers_limit_massive_angle_slices():
    workers = _memory_aware_worker_count(
        num_slices=10,
        height=20_000,
        width=20_000,
        write_angles=True,
        available_memory_bytes=64 * 1024**3,
        cpu_count=32,
    )
    assert workers == 1


def test_memory_aware_workers_respect_slice_and_cpu_counts():
    workers = _memory_aware_worker_count(
        num_slices=3,
        height=64,
        width=64,
        write_angles=True,
        available_memory_bytes=64 * 1024**3,
        cpu_count=2,
    )
    assert workers == 2


def test_remove_padding_returns_views_unless_copy_is_requested():
    volume = np.arange(5 * 2 * 3, dtype=np.float32).reshape(5, 2, 3)
    values = np.broadcast_to(volume, (3, *volume.shape)).copy()
    vectors = values.copy()

    cropped_volume, cropped_values, cropped_vectors = remove_padding(
        volume, values, vectors, padding_start=1, padding_end=1
    )

    np.testing.assert_array_equal(cropped_volume, volume[1:4])
    np.testing.assert_array_equal(cropped_values, values[:, 1:4])
    np.testing.assert_array_equal(cropped_vectors, vectors[:, 1:4])
    assert np.shares_memory(cropped_volume, volume)
    assert np.shares_memory(cropped_values, values)
    assert np.shares_memory(cropped_vectors, vectors)

    copied_volume, copied_values, copied_vectors = remove_padding(
        volume, values, vectors, padding_start=1, padding_end=1, copy=True
    )
    assert not np.shares_memory(copied_volume, volume)
    assert not np.shares_memory(copied_values, values)
    assert not np.shares_memory(copied_vectors, vectors)


def test_structure_tensor_memory_check_fails_before_large_allocation():
    with pytest.raises(MemoryError, match="Reduce N_CHUNK"):
        _ensure_structure_tensor_memory(
            (10, 10, 10),
            available_memory_bytes=10_000,
        )


def test_low_memory_chunk_size_shrinks_failed_64_gib_dataset():
    safe_chunk, requested_peak, safe_peak = _safe_low_memory_chunk_size(
        50,
        10_920,
        7_660,
        7_385,
        np.uint16,
        padding=9,
        include_eigenvalues=True,
        has_mask=True,
        write_angles=True,
        sigma=0.6,
        rho=4.0,
        truncate=2.0,
        available_memory_bytes=64 * 1024**3,
        cpu_count=8,
    )

    assert safe_chunk == 16
    assert requested_peak > 64 * 1024**3
    assert safe_peak <= int(64 * 1024**3 * 0.82)


def test_low_memory_runtime_guard_rejects_oversized_padded_chunk():
    with pytest.raises(MemoryError, match="padded low-memory chunk"):
        _ensure_structure_tensor_memory(
            (59, 7_660, 7_385),
            include_eigenvalues=True,
            low_memory=True,
            available_memory_bytes=64 * 1024**3,
            input_dtype=np.uint16,
            has_mask=True,
            write_angles=True,
            sigma=0.6,
            rho=4.0,
            truncate=2.0,
        )


def test_compute_orientation_keeps_new_options_at_end_of_public_signature():
    assert list(inspect.signature(compute_orientation).parameters)[-3:] == [
        "vector_format",
        "low_memory",
        "low_memory_dir",
    ]


def test_compute_orientation_rejects_no_requested_outputs():
    with pytest.raises(ValueError, match="At least one"):
        compute_orientation("missing-volume", write_vectors=False, write_angles=False)


def test_seed_selection_is_reproducible_and_returns_valid_coordinates():
    seed_mask = np.ones((4, 5, 6), dtype=bool)
    first = _select_seed_points(seed_mask, num_seeds=12, random_seed=42)
    second = _select_seed_points(seed_mask, num_seeds=12, random_seed=42)
    different = _select_seed_points(seed_mask, num_seeds=12, random_seed=7)

    np.testing.assert_array_equal(first, second)
    assert not np.array_equal(first, different)
    assert first.shape == (12, 3)
    assert np.all(seed_mask[tuple(first.T)])


def test_plot_images_with_vector_overlay(tmp_path: Path):
    plt.switch_backend("Agg")

    img = np.arange(36, dtype=np.float32).reshape(6, 6)
    img_angle1 = np.linspace(-90, 90, 36, dtype=np.float32).reshape(6, 6)
    img_angle2 = np.linspace(90, -90, 36, dtype=np.float32).reshape(6, 6)
    img_fa = np.linspace(0, 1, 36, dtype=np.float32).reshape(6, 6)

    vector_field_slice = np.zeros((3, 6, 6), dtype=np.float32)
    vector_field_slice[0, :, :] = 1.0
    vector_field_slice[1, 2:, :] = 0.5

    out_path = tmp_path / "test_slice" / "result_test_slice.png"
    plot_images(
        img,
        img_angle1,
        img_angle2,
        img_fa,
        center_point=(3, 2, 0),
        vector_field_slice=vector_field_slice,
        save_path=str(out_path),
        show=False,
        quiver_step=2,
        angle_ranges=((-90.0, 90.0), (-90.0, 90.0)),
    )

    assert out_path.is_file(), "Diagnostic test-slice figure was not created"
    assert out_path.stat().st_size > 0, "Diagnostic test-slice figure is empty"
