import os
import sys
import glob
import numpy as np
import cv2
import dask
import matplotlib.pyplot as plt
from distutils.util import strtobool
import warnings
from pathlib import Path
import configparser
import tifffile

# Optional GPU support
try:
    import cupy as cp
    USE_GPU = True
    print("GPU support enabled.")
except ImportError:
    USE_GPU = False


from structure_tensor import parallel_structure_tensor_analysis




def convert_to_8bit(img, minimum=0, maximum=100):
    """
    Convert the input image to 8-bit format.

    Parameters:
    - img: numpy.ndarray
        The input image.

    Returns:
    - img_8bit: numpy.ndarray
        The 8-bit converted image.
    """
    # norm = np.zeros(img.shape[:2], dtype=np.float32)
    # img_8bit = cv2.normalize(img, norm, 0, 255, cv2.NORM_MINMAX).astype('uint8')
    
    minimum,maximum = np.nanpercentile(img, (minimum, maximum))

    print(f"Minimum, Maximum : {minimum}, {maximum}")
    
    minimum = int(minimum)
    maximum = int(maximum)
    
    img_normalized = (img + np.abs(minimum)) * (255 / (maximum - minimum))  # Normalize the image to the range [0, 255]
    img_8bit = img_normalized.astype('np.uint8')
    return img_8bit




def read_conf_file(file_path):
    config = configparser.ConfigParser()
    config.read(file_path)

    config_dict = {
        'IMAGES_PATH': config.get('DATASET', 'IMAGES_PATH', fallback=None).strip(),
        'MASK_PATH': config.get('DATASET', 'MASK_PATH', fallback=None).strip(),
        'FLIP': config.getboolean('DATASET', 'FLIP', fallback=True),
        'OUTPUT_PATH': config.get('OUTPUT', 'OUTPUT_PATH', fallback=None).strip(),
        'OUTPUT_TYPE': config.get('OUTPUT', 'OUTPUT_TYPE', fallback=None).strip(),
        'SIGMA': config.getfloat('STRUCTURE TENSOR CALCULATION', 'SIGMA', fallback=None),
        'RHO': config.getfloat('STRUCTURE TENSOR CALCULATION', 'RHO', fallback=None),
        'N_CHUNK': config.getint('STRUCTURE TENSOR CALCULATION', 'N_CHUNK', fallback=100),
        'POINT_MITRAL_VALVE': np.array(list(map(int, config.get('LV AXIS COORDINATES', 'POINT_MITRAL_VALVE', fallback=None).split(',')))),
        'POINT_APEX': np.array(list(map(int, config.get('LV AXIS COORDINATES', 'POINT_APEX', fallback=None).split(',')))),
        'TEST': config.getboolean('TEST', 'TEST', fallback=True),
        'N_SLICE_TEST': config.getint('TEST', 'N_SLICE_TEST', fallback=None)
    }

    return config_dict


def get_image_list(directory):
    """
    Get the list of image paths and their corresponding image types in the directory.
    
    Args:
    - directory (Path): Directory to scan.
    
    Returns:
    - list: List of tuples containing image paths and their corresponding image types.
    - str: Image type.
    """
    if isinstance(directory, str):
        directory = Path(directory)
        
    tif_files = sorted(list(directory.glob('*.tif*')))
    jp2_files = sorted(list(directory.glob('*.jp2')))
    png_files = sorted(list(directory.glob('*.png')))
    edf_files = sorted(list(directory.glob('*.edf')))
    
    # Find the largest list
    largest_list = max([tif_files, jp2_files, png_files, edf_files], key=len)
    img_list = largest_list
    if largest_list == tif_files:
        img_type = 'tif'
    elif largest_list == jp2_files:
        img_type = 'jp2'
    elif largest_list == png_files:
        img_type = 'png'
    elif largest_list == edf_files:
        img_type = 'edf'
    else:
        sys.exit("No image files found in the directory.")
    
    return img_list, img_type







def interpolate_points(point1, point2, N_img):
    """
    Generates interpolated points between two points.

    Parameters:
    - point1: tuple
        The first point (x, y, z).
    - point2: tuple
        The second point (x, y, z).
    - N_img: int
        The number of images or points to interpolate.

    Returns:
    - np.array
        Array of interpolated points.
    """
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    z_values = list(range(N_img))
    
    points = []
    for z in z_values:
        t = (z - z1) / (z2 - z1)  # Calculate the interpolation parameter
        x = x1 + (x2 - x1) * t
        y = y1 + (y2 - y1) * t
        points.append((x, y, z))
    return np.array(points)




def calculate_center_vector(pt_mv, pt_apex, is_flip=True):
    """
    Calculates the center vector between two points.

    Parameters:
    - pt_mv: numpy.ndarray
        The mitral valve point.
    - pt_apex: numpy.ndarray
        The apex point.

    Returns:
    - center_vec: numpy.ndarray
        The normalized center vector.
    """
    # Calculate center vector
    
    center_vec = pt_mv - pt_apex
    center_vec = center_vec / np.linalg.norm(center_vec)
    
    if is_flip:
        center_vec[2] = -center_vec[2]  # Invert z components (Don't know why but it works)

    # center_vec[0] = -center_vec[0]  # Invert z components (Don't know why but it works)
    # center_vec[1] = -center_vec[1]  # Invert z components (Don't know why but it works)


    # if center_vec[2] < 0:
    #     center_vec = -center_vec
        
    return center_vec




def adjust_start_end_index(start_index, end_index, N_img, padding_start, padding_end, is_test, n_slice):
    """
    Adjusts the start and end indices for processing.

    Parameters:
    - start_index: int
        The initial start index.
    - end_index: int
        The initial end index.
    - N_img: int
        Number of images in the volume data.
    - padding_start: int
        Padding to add at the start.
    - padding_end: int
        Padding to add at the end.
    - is_test: bool
        Flag to indicate if it is a test run.

    Returns:
    - start_index_padded: int
        The adjusted start index.
    - end_index_padded: int
        The adjusted end index.
    """
    
    # Validate input indices
    if start_index < 0:
        raise ValueError("start_index must be greater than 0.")
    if end_index < 0:
        raise ValueError("end_index must be greater than 0.")
    
    # Adjust indices for test condition
    if is_test:
        test_index = n_slice
        # else:
        #     test_index = int(N_img / 1.68)
        #     # test_index = 1723

        start_index_padded = max(test_index - padding_start, 0)
        end_index_padded = min(test_index + 1 + padding_end, N_img)
    else:
        # Adjust start and end indices considering padding
        start_index_padded = max(start_index - padding_start, 0)
        end_index_padded = min(end_index + padding_end, N_img)
    

    return start_index_padded, end_index_padded




def load_volume(file_list, start_index_padded, end_index_padded):
    """
    Loads the volume data from a list of file paths.

    Parameters:
    - file_list: list
        List of file paths to load.
    - start_index_padded: int
        The starting index with padding.
    - end_index_padded: int
        The ending index with padding.

    Returns:
    - volume: numpy.ndarray
        The loaded volume data.
    """
    
    count_img = 0
    
    def custom_image_reader(file_path: str) -> np.ndarray:
        """
        Read an image from the given file path.
        Args:
            file_path (str): The path to the image file.
        Returns:
            np.ndarray: The image data as a NumPy array.
        """
        nonlocal count_img  # Declare count_img as non-local to modify it
        print(f"{count_img}/{len(file_list)} - Reading image: {os.path.basename(file_path)}")
        sys.stdout.flush()
        image_data = cv2.imread(str(file_path), cv2.IMREAD_UNCHANGED)
        
        # Add unsharp mask filter
        # image_data = unsharp_mask(image_data, radius=3, amount=0.6)
        
        count_img += 1
        
        return image_data
    
    # Create a list of delayed tasks to read each image
    delayed_tasks = [dask.delayed(lambda x: np.array(custom_image_reader(x)))(file_path) for file_path in sorted(file_list)]
    
    # volume = volume_dask[start_index_padded:end_index_padded,:,:].compute()
    volume = np.array(dask.compute(*delayed_tasks[start_index_padded:end_index_padded]))
    return volume




def calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU):
    """
    Calculates the structure tensor for the given volume data.

    Parameters:
    - volume: numpy.ndarray
        The volume data.
    - SIGMA: float
        Sigma value for structure tensor calculation.
    - RHO: float
        Rho value for structure tensor calculation.
    - USE_GPU: bool
        Flag to use GPU if available.

    Returns:
    - tuple: (s, vec, val)
        The structure tensor (s), eigenvectors (vec), and eigenvalues (val).
    """
    # Filter or ignore specific warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    
    # flag_GPU=False
    # if USE_GPU:
    #     S, vec, val = parallel_structure_tensor_analysis(volume, 
    #                                                     SIGMA, 
    #                                                     RHO, 
    #                                                     devices=['cuda:0'] +['cuda:1']+64*['cpu'], 
    #                                                     block_size=200,
    #                                                     structure_tensor=True)                                                   
    # else:
    S, vec, val = parallel_structure_tensor_analysis(volume, SIGMA, RHO,
                                                        structure_tensor=True) # vec has shape =(3,x,y,z) in the order of (z,y,x)

    return S, vec, val







def remove_padding(volume, s, vec, val, padding_start, padding_end):
    """
    Removes padding from the processed data.

    Parameters:
    - volume: numpy.ndarray
        The volume data.
    - s: numpy.ndarray
        The structure tensor data.
    - vec: numpy.ndarray
        The eigenvectors.
    - val: numpy.ndarray
        The eigenvalues.
    - padding_start: int
        The start padding to remove.
    - padding_end: int
        The end padding to remove.

    Returns:
    - tuple: (volume, s, vec, val)
        The adjusted data without padding.
    """
    array_end = vec.shape[1] - padding_end
    volume = volume[padding_start:array_end,:,:]
    s = s[:,padding_start:array_end,:,:]
    vec = vec[:,padding_start:array_end,:,:]
    val = val[:,padding_start:array_end,:,:]
    
    return volume, s, vec, val


    
    
def calculate_fraction_anisotropy(eigenvalues_2d):
    """
    Calculates Fractional Anisotropy from eigenvalues of a structure tensor.

    Parameters:
    - eigenvalues_2d (numpy.ndarray): Eigenvalues of a slice of the 3D volume (l1, l2, l3).

    Returns:
    - numpy.ndarray: An array of Fractional Anisotropy values.
    """
    l1 = eigenvalues_2d[0, :, :]
    l2 = eigenvalues_2d[1, :, :]
    l3 = eigenvalues_2d[2, :, :]
    mean_eigenvalue  = (l1 + l2 + l3) / 3
    numerator = np.sqrt((l1 - mean_eigenvalue )**2 + (l2 - mean_eigenvalue )**2 + (l3 - mean_eigenvalue )**2)
    denominator = np.sqrt(l1**2 + l2**2 + l3**2)
    img_FA = np.sqrt(3 / 2) * (numerator / denominator)
    
    return img_FA







def rotate_vectors_to_new_axis(vector_field_slice, new_axis_vec):
    """
    Rotates vectors to align with a new axis.

    Parameters:
    - vector_field_slice (numpy.ndarray): Array of vectors to be rotated.
    - new_axis_vec (numpy.ndarray): The new axis to align the vectors with.

    Returns:
    - numpy.ndarray: Rotated vectors aligned with the new axis.
    """
    # Ensure new_axis_vec is normalized
    new_axis_vec = new_axis_vec / np.linalg.norm(new_axis_vec)
    
    # Calculate the rotation matrix
    vec1 = np.array([0, 0, 1])  # Initial vertical axis
    
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (new_axis_vec).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    if np.any(kmat):
        rotation_matrix = np.eye(3) + kmat + np.dot(kmat, kmat) * ((1 - c) / (np.linalg.norm(v) ** 2))
    else:
        rotation_matrix = np.eye(3)
    
    # Reshape vec_2D to (3, N) for matrix multiplication
    vec_2D_reshaped = np.reshape(vector_field_slice, (3, -1))  
    
    vec_2D_reshaped = vec_2D_reshaped / np.linalg.norm(vec_2D_reshaped, axis=0)

    # Rotate the vectors
    rotated_vecs = np.dot(rotation_matrix, vec_2D_reshaped)

    # Reshape back to the original shape
    rotated_vecs = rotated_vecs.reshape(vector_field_slice.shape)
    
    # print(f"Rotation matrix:\n{rotation_matrix}")
    
    return rotated_vecs







def compute_helix_and_transverse_angles(vector_field_2d, center_point):
    """
    Calculate the helix angles at each point of a 2D image using a vector field.

    Parameters:
    - image (numpy.ndarray): The 2D image array.
    - vector_field (numpy.ndarray): The 2D orientation vector field.
    - center (tuple): The center point (x, y) of the image.

    Returns:
    - tuple of numpy.ndarray: Arrays containing the helix and transverse angles.
    """
    center = center_point[0:2]  # Replace with actual values
    rows, cols = vector_field_2d.shape[1:3]

    reshaped_vector_field  = np.reshape(vector_field_2d, (3, -1))

    center_x, center_y  = center[0], center[1]
    
    X, Y = np.meshgrid(np.arange(cols) - center_x, np.arange(rows) - center_y)
    
    theta = -np.arctan2(Y.flatten(), X.flatten())
    cos_angle  = np.cos(theta)
    sin_angle = np.sin(theta)

    # Rotate the vector field
    rotated_vector_field  = np.copy(reshaped_vector_field )
    rotated_vector_field [0, :] = cos_angle * reshaped_vector_field[0, :] - sin_angle * reshaped_vector_field[1, :]
    rotated_vector_field [1, :] = sin_angle * reshaped_vector_field[0, :] + cos_angle * reshaped_vector_field[1, :]


    # Reshape rotated vector field to original image dimensions
    reshaped_rotated_vector_field  = np.zeros((3, rows, cols))
    for i in range(3):
        reshaped_rotated_vector_field[i] = rotated_vector_field[i].reshape(rows, cols)


    # Calculate helix and transverse angles    
    helix_angle = np.arctan(reshaped_rotated_vector_field [2, :, :] / reshaped_rotated_vector_field [1, :, :])
    transverse_angle = np.arctan(reshaped_rotated_vector_field [0, :, :] / reshaped_rotated_vector_field [1, :, :])

    helix_angle = np.rad2deg(helix_angle)
    transverse_angle = np.rad2deg(transverse_angle)

    return helix_angle, transverse_angle










def plot_images(img, img_helix, img_intrusion, img_FA, center_point, PT_MV, PT_APEX):
    """
    Plot images of the heart.

    Args:
        img (numpy.ndarray): The grayscale image of the heart.
        img_helix (numpy.ndarray): The helix image of the heart.
        img_intrusion (numpy.ndarray): The intrusion image of the heart.
        img_FA (numpy.ndarray): The FA (fractional anisotropy) image of the heart.
        center_point (tuple): The coordinates of the center point.
        PT_MV (tuple): The coordinates of the mitral valve point.
        PT_APEX (tuple): The coordinates of the apex point.

    Returns:
        None
    """
    
    
    print("\nPlotting images...")
    img_vmin,img_vmax = np.nanpercentile(img, (5, 95))
    orig_map=plt.get_cmap('hsv')
    reversed_map = orig_map.reversed()
    fig, axes = plt.subplots(2, 2, figsize=(8, 4))
    ax = axes
    ax[0,0].imshow(img, vmin=img_vmin, vmax=img_vmax ,cmap=plt.cm.gray)   
    # Draw a red point at the specified coordinate
    x, y = center_point[0:2]
    ax[0,0].scatter(x, y, c='red', s=50, marker='o')    
    # ax[0,0].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
    # ax[0,0].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
    tmp = ax[0,1].imshow(img_helix,cmap=orig_map)
    ax[1,0].imshow(img_intrusion,cmap=orig_map)  
    tmp2 = ax[1,1].imshow(img_FA,cmap=plt.get_cmap('inferno'))   
    fig.colorbar(tmp) 
                             
    plt.show()


    


def write_images(img_helix, img_intrusion, img_FA, start_index, OUTPUT_DIR, OUTPUT_TYPE, z):
    """
    Write images to the specified output directory.

    Args:
        img_helix (numpy.ndarray): The image data for helix.
        img_intrusion (numpy.ndarray): The image data for intrusion.
        img_FA (numpy.ndarray): The image data for FA (False Alarm).
        start_index (int): The starting index for the image filenames.
        OUTPUT_DIR (str): The output directory to save the images.

    Returns:
        None
    """
  
    # # Convert the float64 image to int8
    # img_helix = convert_to_8bit(img_helix, -90, 90)
    # img_intrusion = convert_to_8bit(img_intrusion, -90, 90)
    # img_FA = convert_to_8bit(img_FA, 0, 1)
    
    os.makedirs(OUTPUT_DIR + '/HA', exist_ok=True)
    os.makedirs(OUTPUT_DIR + '/IA', exist_ok=True)
    os.makedirs(OUTPUT_DIR + '/FA', exist_ok=True)        

    print(f"Saving image:")

    if '8bit' in OUTPUT_TYPE:

        cv2.imwrite(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.tif", img_helix)
        cv2.imwrite(f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.tif", img_intrusion)
        cv2.imwrite(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.tif", img_FA)
        
    elif 'rgb' in OUTPUT_TYPE:
        def write_img_rgb(img, output_path, cmap=plt.get_cmap('hsv')):
            minimum = np.nanmin(img)
            maximum = np.nanmax(img)
            img = (img + np.abs(minimum)) * (1 / (maximum - minimum)) 
            img = cmap(img)
            img = (img[:, :, :3] * 255).astype(np.uint8)
            
            # cv2.imwrite(output_path, img)
            print(f"Writing image to {output_path}")
            tifffile.imwrite(output_path, img)

        write_img_rgb(img_helix,f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.tif", cmap=plt.get_cmap('hsv'))
        write_img_rgb(img_intrusion,f"{OUTPUT_DIR}/IA/IA_{(start_index + z):06d}.tif", cmap=plt.get_cmap('hsv'))
        write_img_rgb(img_FA,f"{OUTPUT_DIR}/FA/FA_{(start_index + z):06d}.tif", cmap=plt.get_cmap('inferno'))

