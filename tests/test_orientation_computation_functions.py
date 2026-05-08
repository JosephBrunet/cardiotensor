from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cardiotensor.orientation.orientation_computation_functions import (
    adjust_start_end_index,
    calculate_structure_tensor,
    compute_azimuth_and_elevation,
    compute_fraction_anisotropy,
    compute_helix_and_transverse_angles,
    interpolate_points,
    orient_vectors_z_positive,
    plot_images,
    rotate_vectors_to_new_axis,
    write_images,
)
from cardiotensor.utils.image_io import write_vector_field


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
    assert np.allclose(rotated_same, vector_slice, atol=1e-6), (
        "Rotation to same axis should not change vectors"
    )


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


def test_compute_helix_and_transverse_angles_orients_circumferential_positive():
    vector_slice = np.zeros((3, 1, 2), dtype=np.float32)
    vector_slice[:, 0, 0] = [0, 1, 1]
    vector_slice[:, 0, 1] = [0, -1, -1]

    helix, transverse = compute_helix_and_transverse_angles(
        vector_slice, center_point=(0, 0, 0)
    )

    assert np.allclose(helix, [[45, 45]], atol=1e-6)
    assert np.allclose(transverse, [[0, 0]], atol=1e-6)


def test_write_images_and_vectors(tmp_path: Path):
    # --- Prepare dummy 2D slices ---
    img_helix = np.ones((5, 5), dtype=np.float32)
    img_intrusion = np.ones((5, 5), dtype=np.float32)
    img_FA = np.ones((5, 5), dtype=np.float32)

    # --- Test write_images ---
    write_images(
        img_helix,
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
