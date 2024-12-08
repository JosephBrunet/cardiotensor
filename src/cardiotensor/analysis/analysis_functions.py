import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skimage.measure import profile_line


def _calculate_angle_line(start, end):
    """
    Calculate the angle of a line defined by start and end points.

    Parameters:
    start (tuple): The starting point of the line (x1, y1).
    end (tuple): The ending point of the line (x2, y2).

    Returns:
    float: The angle of the line in degrees.
    """
    x1, y1 = start
    x2, y2 = end

    delta_x = x2 - x1
    delta_y = y2 - y1

    # Calculate the angle in radians
    angle_rad = math.atan2(delta_y, delta_x)

    # Convert the angle to degrees
    angle_deg = math.degrees(angle_rad)

    return angle_deg


def find_end_points(start_point, end_point, angle_range, N_line):
    theta = _calculate_angle_line(start_point, end_point)

    if N_line > 1:
        theta_list = np.linspace(
            theta - angle_range / 2, theta + angle_range / 2, N_line
        )
    else:
        theta_list = [theta]

    vector = np.array(end_point) - np.array(start_point)
    norm = np.linalg.norm(vector)

    # Calculate end points
    end_points = []
    for angle in theta_list:
        theta = np.deg2rad(angle)
        end_x = int(start_point[0] + norm * np.cos(theta))
        end_y = int(start_point[1] + norm * np.sin(theta))
        end_points.append((end_x, end_y))

    return np.array(end_points)


def calculate_intensities(
    img_helix,
    start_point,
    end_point,
    angle_range=5,
    N_line=10,
    max_value=None,
    min_value=None,
):
    end_points = find_end_points(start_point, end_point, angle_range, N_line)

    img_helix[np.isnan(img_helix)] = 0

    intensity_profiles = []
    for i, end in enumerate(end_points):
        print(f"Measure {i+1}/{len(end_points)}")
        intensity_profile = profile_line(img_helix, start_point, end, order=0)

        if min_value != None and max_value != None:
            intensity_profile = (
                intensity_profile * (max_value - min_value) / 255 + min_value
            )

        intensity_profiles.append(intensity_profile)

    return intensity_profiles


def plot_intensity(
    intensity_profiles,
    label_y="",
    x_max_lim=None,
    x_min_lim=None,
    y_max_lim=None,
    y_min_lim=None,
):
    plt.figure(figsize=(10, 6))

    # Get the minimum length of the profiles
    min_length = min(
        intensity_profile.shape[0] for intensity_profile in intensity_profiles
    )

    # Trim the arrays to the minimum length
    trimmed_arrays = [
        intensity_profile[:min_length] for intensity_profile in intensity_profiles
    ]

    # Convert list of trimmed arrays to a 2D NumPy array
    intensity_profiles = np.stack(trimmed_arrays)

    # Calculate mean and median arrays
    mean_array = np.mean(intensity_profiles, axis=0)

    # Calculate the 5th and 95th percentiles
    lower_percentile = np.percentile(intensity_profiles, 5, axis=0)
    upper_percentile = np.percentile(intensity_profiles, 95, axis=0)

    # Create a normalised x-axis ranging from 0 to 1
    x_axis_arr = np.linspace(0, 1, len(mean_array))

    # Plot the mean
    plt.plot(
        x_axis_arr, mean_array, label="Mean", linewidth=2, color="k"
    )  # Make the line thicker for better visibility

    # Add shaded area for the 5th to 95th percentiles
    plt.fill_between(
        x_axis_arr,
        lower_percentile,
        upper_percentile,
        color="gray",
        alpha=0.5,
        label="5%-95% Percentiles",
    )

    # Increase axis and label thickness
    ax = plt.gca()  # Get current axis
    ax.spines["top"].set_linewidth(2)
    ax.spines["right"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)

    # Increase tick width and font size
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    plt.xticks(fontsize=12, weight="bold")
    plt.yticks(fontsize=12, weight="bold")

    # Set thicker labels with larger font size
    plt.xlabel("Normalised Transmural Depth", fontsize=14, weight="bold")
    plt.ylabel(label_y, fontsize=14, weight="bold")

    # Set axis limits if provided
    if x_max_lim:
        plt.xlim(
            [x_min_lim, x_max_lim]
        )  # Set the x-axis limits, e.g., plt.xlim([0, 1])
    if y_max_lim:
        plt.ylim(
            [y_min_lim, y_max_lim]
        )  # Set the y-axis limits, e.g., plt.ylim([min_value, max_value])

    # Show legend and plot
    # plt.legend(fontsize=12)
    plt.show()


def save_intensity(intensity_profiles, save_path):
    # Convert intensity_profiles into a DataFrame
    df = pd.DataFrame(intensity_profiles)

    # Rename columns to start from "Value 1"
    df.columns = [f"Value {i+1}" for i in range(df.shape[1])]

    # Add a "Profile" row
    df.insert(0, "Profile", [f"Profile {i+1}" for i in range(df.shape[0])])

    # Transpose the DataFrame so profiles are in columns
    df = df.transpose()

    # Set the first row as the header
    df.columns = df.iloc[0]
    df = df[1:]

    # Save DataFrame to CSV with a semicolon delimiter
    df.to_csv(save_path, sep=";", index=False)

    print(f"Profile saved to {save_path}")
