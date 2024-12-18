import math
import multiprocessing as mp
import os
import sys
from pathlib import Path

import cv2
import numpy as np
from alive_progress import alive_bar
from cardiotensor.utils.utils import (
    convert_to_8bit,
    get_volume_shape,
    load_volume,
    read_conf_file,
)
from skimage.measure import block_reduce


def process_vector_block(block, bin_factor, h, w, output_dir, idx):
    """
    Function to process a single block of numpy files and save the downsampled output.
    """
    print(f"Processing block: {idx}")
    output_file = f"{output_dir}/eigen_vec/eigen_vec_{idx:06d}.npy"

    # Skip processing if the output already exists
    if os.path.exists(output_file):
        print(f"Skipping block {idx}, file already exists.")
        return

    array = np.empty((3, len(block), h, w))
    bin_array = np.empty((3, math.ceil(h / bin_factor), math.ceil(w / bin_factor)))

    for i, p in enumerate(block):
        print(f"Reading file: {p.name}")
        array[:, i, :, :] = np.load(p)  # Load the numpy data

    # Compute the mean across the stack
    array = array.mean(axis=1)

    # Define the block size for binning
    block_size = (bin_factor, bin_factor)

    # Apply block_reduce to downsample the volume
    bin_array[0, :, :] = block_reduce(
        array[0, :, :], block_size=block_size, func=np.mean
    )
    bin_array[1, :, :] = block_reduce(
        array[1, :, :], block_size=block_size, func=np.mean
    )
    bin_array[2, :, :] = block_reduce(
        array[2, :, :], block_size=block_size, func=np.mean
    )

    bin_array = bin_array.astype(np.float32)

    # Save the downsampled array
    np.save(output_file, bin_array)
    print(f"Saved block {idx} to {output_file}")


def downsample_vector_volume(input_npy, bin_factor, output_dir):
    """
    Downsamples a vector volume using multiprocessing.
    """
    output_dir = Path(output_dir) / f"bin{bin_factor}"
    os.makedirs(output_dir / "eigen_vec", exist_ok=True)

    npy_list = sorted(list(input_npy.glob("*.npy")))
    _, h, w = np.load(npy_list[0]).shape

    # Split files into blocks based on bin_factor
    blocks = [npy_list[i : i + bin_factor] for i in range(0, len(npy_list), bin_factor)]

    # Prepare the arguments for each block
    tasks = [
        (block, bin_factor, h, w, output_dir, idx) for idx, block in enumerate(blocks)
    ]

    # Use multiprocessing with apply_async
    with mp.Pool(processes=min(mp.cpu_count(), 16)) as pool:
        results = [
            pool.apply_async(
                process_vector_block, args=(block, bin_factor, h, w, output_dir, idx)
            )
            for block, bin_factor, h, w, output_dir, idx in tasks
        ]
        # Wait for all processes to finish
        [result.wait() for result in results]


def process_image_block(block, bin_factor, h, w, output_dir, idx):
    """
    Function to process a single block of image files and save the downsampled output.
    """
    print(f"Processing block: {idx}")
    output_file = f"{output_dir}/HA/HA_{idx:06d}.tif"

    # Skip processing if the output already exists
    if os.path.exists(output_file):
        print(f"Skipping block {idx}, file already exists.")
        return

    array = np.empty((len(block), h, w))
    bin_array = np.empty((math.ceil(h / bin_factor), math.ceil(w / bin_factor)))

    for i, p in enumerate(block):
        print(f"Reading file: {p.name}")
        array[i, :, :] = cv2.imread(
            str(p), cv2.IMREAD_UNCHANGED
        )  # Load the image as numpy array

    # Compute the mean across the stack
    array = array.mean(axis=0)

    # Define the block size for binning
    block_size = (bin_factor, bin_factor)

    # Apply block_reduce to downsample the volume
    bin_array[:, :] = block_reduce(array[:, :], block_size=block_size, func=np.mean)

    bin_array = convert_to_8bit(bin_array, output_min=-90, output_max=90)

    # Save the downsampled image
    cv2.imwrite(output_file, bin_array)
    print(f"Saved block {idx} to {output_file}")


def downsample_volume(input_path, bin_factor, output_dir):
    output_dir = output_dir / f"bin{bin_factor}"
    os.makedirs(output_dir / "HA", exist_ok=True)

    HA_list = sorted(list(input_path.glob("*.tif")))

    h, w = cv2.imread(str(HA_list[0]), cv2.IMREAD_UNCHANGED).shape

    blocks = [HA_list[i : i + bin_factor] for i in range(0, len(HA_list), bin_factor)]

    # Prepare the arguments for each block
    tasks = [
        (block, bin_factor, h, w, output_dir, idx) for idx, block in enumerate(blocks)
    ]

    print(f"Total blocks to process: {len(blocks)}")

    # Use multiprocessing with apply_async
    with mp.Pool(processes=min(mp.cpu_count(), 16)) as pool:
        results = [
            pool.apply_async(
                process_image_block, args=(block, bin_factor, h, w, output_dir, idx)
            )
            for block, bin_factor, h, w, output_dir, idx in tasks
        ]
        # Wait for all processes to finish
        [result.wait() for result in results]

    print("All blocks processed.")


def angle_between_vectors(vec1, vec2):
    """
    Calculate the element-wise angle between two vector fields.

    Parameters:
    vec1, vec2: np.ndarray
        Arrays of shape (3, z, y, x) representing the 3D vector fields.

    Returns:
    np.ndarray
        Array of shape (z, y, x) with the angle in radians between vec1 and vec2 at each point.
    """
    # Calculate the dot product along the vector component axis (axis=0)
    dot_product = np.sum(vec1 * vec2, axis=0)

    # Calculate the magnitudes of each vector along axis=0
    magnitude_vec1 = np.linalg.norm(vec1, axis=0)
    magnitude_vec2 = np.linalg.norm(vec2, axis=0)

    # Avoid division by zero by adding a small epsilon where magnitudes are zero
    epsilon = 1e-10
    magnitude_vec1 = np.maximum(magnitude_vec1, epsilon)
    magnitude_vec2 = np.maximum(magnitude_vec2, epsilon)

    # Calculate cosine of the angle between vectors
    cos_theta = dot_product / (magnitude_vec1 * magnitude_vec2)

    # Clamp values to the valid range for arccos to avoid numerical errors
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    # Calculate the angle in radians
    theta = np.arccos(cos_theta)

    return np.rad2deg(theta)


# Function to find consecutive points in the vector field direction
def find_consecutive_points(
    start_point, vector_field, num_steps=4, segment_length=10, angle_threshold=60
):
    """
    Given a starting point, find consecutive points in the direction specified by the vector field.

    Parameters:
        start_point (tuple): The starting point (z, y, x).
        vector_field (np.array): The vector field with shape (3, Z, Y, X).
        num_steps (int): Number of consecutive steps to take in the vector direction.

    Returns:
        list: A list of points (z, y, x) following the vector direction.
    """
    consecutive_points = [start_point]
    current_point = start_point
    direction_vector_tmp = np.array([])

    for _ in range(num_steps):
        # Get the vector at the current point
        z, y, x = current_point

        if (
            0 <= z < vector_field.shape[1]
            and 0 <= y < vector_field.shape[2]
            and 0 <= x < vector_field.shape[3]
        ):
            direction_vector = vector_field[:, z, y, x] * segment_length
            if np.isnan(direction_vector).any():
                break
            if direction_vector_tmp.any():
                # print(f"Angle: {angle_between_vectors(direction_vector_tmp, direction_vector)}")
                if (
                    angle_between_vectors(direction_vector_tmp, direction_vector)
                    > angle_threshold
                ):
                    break
            # Calculate the next point by moving in the direction of the vector
            next_point = (
                int(np.round(z + direction_vector[0])),
                int(np.round(y + direction_vector[1])),
                int(np.round(x + direction_vector[2])),
            )

            # Ensure the point is within bounds
            if (
                0 <= next_point[0] < vector_field.shape[1]
                and 0 <= next_point[1] < vector_field.shape[2]
                and 0 <= next_point[2] < vector_field.shape[3]
            ):
                consecutive_points.append(next_point)
                current_point = next_point
            else:
                break  # Stop if we go out of bounds
        else:
            break  # Stop if the current point is out of bounds

        direction_vector_tmp = direction_vector

    return consecutive_points


def write_am_file(consecutive_points_list, HA_angle, z_angle, file_path="output.am"):
    """
    Write an .am file with start and end vertices for each element in consecutive_points_list.

    Parameters:
        consecutive_points_list (list): List of lists of consecutive points. Each inner list contains tuples of points (x, y, z).
        file_path (str): Path to the output .am file.
    """

    N_point = 0
    for i in consecutive_points_list:
        N_point += len(i)

    with open(file_path, "w") as f:
        # Write header information
        f.write("# AmiraMesh 3D ASCII 3.0\n\n\n")
        f.write(
            f"define VERTEX {len(consecutive_points_list) * 2}\n"
        )  # Two vertices per line
        f.write(
            f"define EDGE {len(consecutive_points_list)}\n"
        )  # One edge per line segment
        f.write(f"define POINT {N_point}\n")  # One edge per line segment
        f.write("\nParameters {\n")
        f.write('    ContentType "HxSpatialGraph"\n')
        f.write("}\n\n")

        f.write("VERTEX { float[3] VertexCoordinates } @1\n")
        f.write("EDGE { int[2] EdgeConnectivity } @2\n")
        f.write("EDGE { int NumEdgePoints } @3\n")
        f.write("POINT { float[3] EdgePointCoordinates } @4\n")
        f.write("POINT { float thickness } @5\n")
        f.write("POINT { float HA_angle } @6 \n")
        f.write("POINT { float z_angle } @7 \n")

        f.write("\n# Data section follows\n")

        f.write("@1\n")
        vertex = []
        for segment in consecutive_points_list:
            if len(segment) >= 2:
                # Start and end points for each segment
                start = segment[0]
                end = segment[-1]
                f.write(f"{start[0]} {start[1]} {start[2]}\n")
                f.write(f"{end[0]} {end[1]} {end[2]}\n")
                vertex.append((start, end))

        f.write("\n")

        # Write edges connecting consecutive points
        f.write("@2\n")
        for i in range(len(vertex)):
            f.write(f"{i * 2} {i * 2 + 1}\n")  # Connect start to end for each segment

        f.write("\n")

        # Write number of points per edge
        f.write("@3\n")
        for segment in consecutive_points_list:
            f.write(f"{len(segment)}\n")

        f.write("\n")

        # Write number of points per edge
        f.write("@4\n")
        for segment in consecutive_points_list:
            for point in segment:
                f.write(" ".join(map(str, point)) + "\n")

        f.write("\n")

        # # Add an optional "Edge Data" section if needed (e.g., thickness of edges)
        # f.write("@5\n")
        # for _ in range(len(vertex)):
        #     f.write("1.0\n")  # Example thickness value, can be adjusted as needed

        # Add an optional "Edge Data" section if needed (e.g., thickness of edges)
        f.write("@5\n")
        for segment in consecutive_points_list:
            for point in segment:
                f.write("1.0\n")  # Example thickness value, can be adjusted as needed

        f.write("\n")

        # Add an optional "Edge Data" section if needed (e.g., thickness of edges)
        f.write("@6\n")
        for i, angle in enumerate(HA_angle):
            f.write(f"{angle}\n")  # Example thickness value, can be adjusted as needed

        f.write("\n")

        f.write("@7\n")
        for i, angle in enumerate(z_angle):
            f.write(f"{angle}\n")

    print(f"Amira file written to {file_path}")


def amira_writer(
    conf_file_path,
    pixel_size_um,
    start_index=None,
    end_index=None,
    bin_factor=None,
    num_ini_points=20000,
    num_steps=1000000,
    segment_length=20,
    angle_threshold=60,
    segment_min_length_threshold=30,
):
    if not start_index:
        start_index = 0

    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(e)
        sys.exit(f"⚠️  Error reading parameter file: {conf_file_path}")

    (
        VOLUME_PATH,
        MASK_PATH,
        IS_FLIP,
        OUTPUT_DIR,
        OUTPUT_TYPE,
        SIGMA,
        RHO,
        N_CHUNK,
        PT_MV,
        PT_APEX,
        IS_TEST,
        N_SLICE_TEST,
    ) = (
        params[key]
        for key in [
            "IMAGES_PATH",
            "MASK_PATH",
            "FLIP",
            "OUTPUT_PATH",
            "OUTPUT_TYPE",
            "SIGMA",
            "RHO",
            "N_CHUNK",
            "POINT_MITRAL_VALVE",
            "POINT_APEX",
            "TEST",
            "N_SLICE_TEST",
        ]
    )

    OUTPUT_DIR = Path(OUTPUT_DIR)

    w, h, N_img = get_volume_shape(VOLUME_PATH)
    if end_index is None:
        end_index = N_img

    output_npy = OUTPUT_DIR / "eigen_vec"
    output_HA = OUTPUT_DIR / "HA"

    if bin_factor:
        downsample_vector_volume(output_npy, bin_factor, OUTPUT_DIR)
        output_npy = OUTPUT_DIR / f"bin{bin_factor}/eigen_vec"

        downsample_volume(output_HA, bin_factor, OUTPUT_DIR)
        output_HA = OUTPUT_DIR / f"bin{bin_factor}/HA"

        start_index = int(start_index / bin_factor)
        end_index = int(end_index / bin_factor)

    npy_list = sorted(list(output_npy.glob("*.npy")))  # [start_index:end_index]

    vector_field = load_volume(npy_list, start_index=start_index, end_index=end_index)
    vector_field = np.moveaxis(vector_field, 0, 1)

    print("\nAlign vectors in same direction")
    # Flip the vectors where the z-component is negative
    vector_field[:, vector_field[0] > 0] *= -1

    print("Mask creation")
    # Update mask_volume to be 0 where any component of vector_field is NaN
    mask_volume = (~np.isnan(vector_field).any(axis=0)).astype(np.uint8)

    # ---------------------------------------------
    # HA

    HA_list = sorted(list(output_HA.glob("*.tif"))) + sorted(
        list(output_HA.glob("*.jp2"))
    )
    HA_volume = load_volume(HA_list, start_index=start_index, end_index=end_index)
    # mask_volume = np.where(HA_volume == 0, 0, 1)

    # HA_volume = HA_volume *90/255 - 90

    # fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    # # Display the original mask slice
    # axes[0].imshow(mask_volume[5, :, :], cmap='gray')
    # axes[0].set_title("Original Mask Volume Slice")

    # kernel = np.ones((3, 3), np.uint8)  # A 3x3 kernel; adjust as needed for the erosion effect
    # # Initialize an empty array to store the eroded volume
    # eroded_volume = np.zeros_like(mask_volume)
    # # Perform erosion slice by slice along the chosen axis (e.g., the z-axis)
    # for i in range(mask_volume.shape[0]):  # Loop through each 2D slice
    #     eroded_volume[i] = cv2.erode(mask_volume[i], kernel, iterations=1)
    # mask_volume = eroded_volume

    # # Plot the eroded mask_volume for the same slice
    # axes[1].imshow(eroded_volume[5, :, :], cmap='gray')
    # axes[1].set_title("Eroded Mask Volume Slice")
    # # Display the plots
    # plt.show()

    print("\nCreation of random points")
    # Get all indices where mask_volume is 1
    valid_indices = np.argwhere(mask_volume == 1)

    if len(valid_indices) < num_ini_points:
        print(
            "Not enough points with mask value 1. Adjust the number of points or check mask_volume."
        )
    else:
        # Randomly select num_ini_points indices from the valid indices
        random_points = valid_indices[
            np.random.choice(valid_indices.shape[0], num_ini_points, replace=False)
        ]

        print("Random 3D points within mask:\n", random_points)
    print("Random points created")

    # Iterate over each point in random_points to find consecutive points
    consecutive_points_list = []
    with alive_bar(len(random_points), title="Processing Points") as bar:
        for point in random_points:
            # print(f"\nCreation of fiber from point {point}")

            point = tuple(int(x) for x in point)
            consecutive_points = find_consecutive_points(
                point,
                vector_field,
                num_steps=num_steps,
                segment_length=segment_length,
                angle_threshold=angle_threshold,
            )
            # Add the start_index to offset correctly
            consecutive_points = [
                (point[0] + start_index, point[1], point[2])
                for point in consecutive_points
            ]

            if len(consecutive_points) >= segment_min_length_threshold:
                consecutive_points_list.append(consecutive_points)
                # print(f"\nCreated a fiber of {len(consecutive_points)} points")
            # else:
            #     print(f"\nFiber to small")

            # Update the progress bar after processing each point
            bar()

    # # Print results
    # for idx, points in enumerate(consecutive_points_list):
    #     print(f"Starting from point {random_points[idx]}:\nConsecutive points: {points}\n")

    print(f"{len(consecutive_points_list)}")

    HA_angle = []
    for point_list in consecutive_points_list:
        for point in point_list:
            HA_angle.append(
                float(HA_volume[point[0] - start_index, point[1], point[2]])
            )

    z_angle = []
    for point_list in consecutive_points_list:
        for point in point_list:
            vector = vector_field[:, point[0] - start_index, point[1], point[2]]

            # Calculate the angle in radians
            theta = np.arccos(abs(vector[0]) / np.linalg.norm(vector))

            # Convert the angle to degrees if needed
            theta_degrees = np.degrees(theta)
            z_angle.append(theta_degrees)

    if bin_factor:
        pixel_size_um = pixel_size_um * bin_factor

    print(f"Voxel size: {pixel_size_um}um")

    # Multiply each element by pixel size
    consecutive_points_list = scale_points(consecutive_points_list, pixel_size_um)

    # Reorder each point in each list from (z, y, x) to (x, y, z)
    consecutive_points_list = [
        [(point[2], point[1], point[0]) for point in point_list]
        for point_list in consecutive_points_list
    ]

    write_am_file(
        consecutive_points_list, HA_angle, z_angle, file_path=OUTPUT_DIR / "output.am"
    )


def scale_points(consecutive_points, pixel_size):
    """
    Scales each coordinate in a list of points by the specified pixel size.

    Parameters:
        consecutive_points (list of list of tuples): A nested list where each sublist contains tuples representing points (x, y, z).
        pixel_size (float): The scaling factor for each coordinate (e.g., pixel size in micrometers).

    Returns:
        list of list of tuples: A new list with each point's coordinates scaled by pixel_size.
    """
    scaled_points = []

    for point_list in consecutive_points:
        # Scale each point within the current list
        scaled_point_list = [
            tuple(coord * pixel_size for coord in point) for point in point_list
        ]
        scaled_points.append(scaled_point_list)

    return scaled_points
