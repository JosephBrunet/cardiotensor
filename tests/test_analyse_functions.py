import numpy as np
import pytest

from cardiotensor.analysis.analysis_functions import (
    _calculate_angle_line,
    calculate_intensities,
    find_end_points,
    save_intensity,
)


def test_calculate_angle_line():
    # Test with a horizontal line
    assert _calculate_angle_line((0, 0), (10, 0)) == 0

    # Test with a vertical line
    assert _calculate_angle_line((0, 0), (0, 10)) == 90

    # Test with a 45-degree line
    assert _calculate_angle_line((0, 0), (10, 10)) == 45


def test_find_end_points():
    start_point = (0, 0)
    end_point = (10, 0)
    angle_range = 90
    N_line = 5

    # Generate endpoints for lines at different angles
    endpoints = find_end_points(start_point, end_point, angle_range, N_line)

    # Check the number of generated endpoints
    assert len(endpoints) == N_line

    # Ensure all endpoints are tuples with two elements
    assert all(len(pt) == 2 for pt in endpoints)


def test_calculate_intensities():
    img_helix = np.zeros((100, 100))
    img_helix[50, 50:] = np.linspace(0, 1, 50)

    start_point = (50, 50)
    end_point = (50, 99)

    # Calculate intensity profiles
    intensity_profiles = calculate_intensities(
        img_helix, start_point, end_point, angle_range=10, N_line=3
    )

    # Check the number of intensity profiles matches N_line
    assert len(intensity_profiles) == 3

    # Check the length of the intensity profiles matches the line length
    assert all(len(profile) == 50 for profile in intensity_profiles)


# def test_plot_intensity():
#     intensity_profiles = [np.linspace(0, 1, 50), np.linspace(1, 0, 50)]

#     # Plot intensity profiles and ensure no errors occur
#     try:
#         plot_intensity(intensity_profiles, label_y="Intensity")
#     except Exception as e:
#         pytest.fail(f"Plotting failed with error: {e}")


def test_save_intensity(tmp_path):
    intensity_profiles = [np.linspace(0, 1, 50), np.linspace(1, 0, 50)]
    save_path = tmp_path / "intensity_profiles.csv"

    # Save intensity profiles
    save_intensity(intensity_profiles, str(save_path))

    # Verify the file exists
    assert save_path.exists()

    # Verify the file content
    with open(save_path) as f:
        content = f.read()
        assert "Profile 1" in content
        assert "Profile 2" in content


if __name__ == "__main__":
    pytest.main(["-v", __file__])
