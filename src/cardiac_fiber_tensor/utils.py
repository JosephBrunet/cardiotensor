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

# Optional GPU support
try:
    import cupy as cp
    USE_GPU = True
except ImportError:
    USE_GPU = False


from structure_tensor import parallel_structure_tensor_analysis




def convert_to_8bit(img, minimum=None, maximum=None):
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
    
    minimum,maximum = np.nanpercentile(img, (1, 99))
    if maximum == 0:
         minimum,maximum = np.nanpercentile(img, (0, 100))
    print(f"Minimum, Maximum : {minimum}, {maximum}")
    
    minimum = int(minimum)
    maximum = int(maximum)
    
    img_normalized = (img + np.abs(minimum)) * (255 / (maximum - minimum))  # Normalize the image to the range [0, 255]
    img_8bit = img_normalized.astype('uint8')
    return img_8bit


def read_parameter_file(para_file_path):
    # Define expected parameters with default values (if any)
    params = {}

    # Define parameter processing functions
    param_processors = {
        'SIGMA': float,
        'RHO': float,
        'POINT_MITRAL_VALVE': lambda x: np.fromstring(x, dtype=int, sep=','),
        'POINT_APEX': lambda x: np.fromstring(x, dtype=int, sep=','),
        'TEST': lambda x: bool(strtobool(x.strip()))
    }

    with open(para_file_path, 'r') as file:
        for line in file:
            line = line.split('#')[0].strip().replace(' ','')  # Remove comments and whitespace
            
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip().upper()  # Normalize key
                
                # Process specific parameters
                processor = param_processors.get(key, str.strip)
                params[key] = processor(value)

                print(f"{key} found ({params[key]})")

    return params




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




def adjust_start_end_index(start_index, end_index, N_img, padding_start, padding_end, is_test, N_SLICE=0):
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
    
    # Set end_index to total_images if it's zero
    if end_index == 0:
        end_index = N_img

    # Adjust indices for test condition
    if is_test:
        if N_SLICE != 0:
            test_index = N_SLICE
        else:
            test_index = int(N_img / 1.68)
            # test_index = 1723

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
    
    def custom_image_reader(file_path: str) -> np.ndarray:
        """
        Read an image from the given file path.
        Args:
            file_path (str): The path to the image file.
        Returns:
            np.ndarray: The image data as a NumPy array.
        """
        print(f"Reading image: {os.path.basename(file_path)}")
        sys.stdout.flush()
        image_data = cv2.imread(str(file_path), cv2.IMREAD_UNCHANGED)
        
        # Add unsharp mask filter
        # image_data = unsharp_mask(image_data, radius=3, amount=0.6)
        
        return image_data
    
    # Create a list of delayed tasks to read each image
    delayed_tasks = [dask.delayed(lambda x: np.array(custom_image_reader(x)))(file_path) for file_path in file_list]
    
    # volume = volume_dask[start_index_padded:end_index_padded,:,:].compute()
    volume = np.array(dask.compute(*delayed_tasks[start_index_padded:end_index_padded]))
    return volume





import numpy as np
from scipy import ndimage
from scipy.signal import convolve2d
from scipy.ndimage import convolve

def structure_tensor_3d_test2(volume):
    # Assuming I3d is your 3D image and FilterType and flip are defined
    I3d = np.transpose(volume, (1, 2, 0))
    
    h, w, d = I3d.shape
    wd = w * d
    hw = h * w
    k = (d + 1) // 2

    flip = 0
    FilterType = 1
    # Calculate gradients
    if FilterType == 1:  # Finite difference
        gx, gy, gz = np.gradient(I3d)
        if flip == 0:
            gz = -gz
    else:  # Sobel filter
        gx = ndimage.sobel(I3d, axis=0, mode='constant')
        gy = ndimage.sobel(I3d, axis=1, mode='constant')
        gz = ndimage.sobel(I3d, axis=2, mode='constant')
        if flip == 0:
            gz = -gz

    # Calculate sums of gradient products
    nT = d - 2
    sumMask3 = np.ones((1, 1, nT))
    sumMask22 = np.ones((nT, nT))

    def convolve_sum(gradient_product, h, w, d, sumMask22, sumMask3, k):
        # Reshape for 2D convolution
        tmp2d = np.reshape(gradient_product, (h, w * d))
        tmp2d = convolve(tmp2d, sumMask22, mode='constant', cval=0)

        # Reshape for 3D convolution and perform it
        tmp3d = np.reshape(tmp2d, (h, w, d))
        tmp3d = convolve(tmp3d, sumMask3, mode='constant', cval=0)

        # Extract the central slice
        return np.reshape(tmp3d[:, :, k], (h * w,))

    
    # Define convolution masks
    sumMask3 = np.ones((1, 1, d-2))
    sumMask22 = np.ones((d-2, d-2))

    # Convolve and compute sums for each gradient product
    Sgxgx = convolve_sum(gx * gx, h, w, d, sumMask22, sumMask3, k)
    Sgygy = convolve_sum(gy * gy, h, w, d, sumMask22, sumMask3, k)
    Sgzgz = convolve_sum(gz * gz, h, w, d, sumMask22, sumMask3, k)
    Sgxgy = convolve_sum(gx * gy, h, w, d, sumMask22, sumMask3, k)
    Sgxgz = convolve_sum(gx * gz, h, w, d, sumMask22, sumMask3, k)
    Sgygz = convolve_sum(gy * gz, h, w, d, sumMask22, sumMask3, k)
    
    
    
    
    
        
    
    
    d = Sgxgx * (-Sgygz**2 + Sgygy * Sgzgz) + Sgxgz * (Sgxgy * Sgygz - Sgxgz * Sgygy) + Sgxgy * (Sgxgz * Sgygz - Sgxgy * Sgzgz)
    c = Sgxgy**2 + Sgxgz**2 + Sgygz**2 - Sgxgx * (Sgygy + Sgzgz) - Sgygy * Sgzgz
    b = Sgxgx + Sgygy + Sgzgz
    a = -np.ones(d.shape)

    ind_a0 = a == 0
    ind_b0 = b == 0
    ind_c0 = c == 0

    hw = d.size  # Assuming hw is the number of elements in d
    Eigens = np.zeros((hw, 3))
    nn = np.full(hw, np.nan)

    ind_ord0 = ind_a0 & ind_b0 & ind_c0
    ind_ord1 = ind_a0 & ind_b0 & ~ind_c0
    ind_ord2 = np.logical_xor(np.logical_xor(ind_a0, ind_ord0), ind_ord1)

    delta = b * b - 4 * a * c
    sqdelta = np.sqrt(delta)
    ind_ord2_1 = ind_ord2 & (delta == 0)
    ind_ord2_2 = ind_ord2 & (delta != 0)
    inva = 0.5 / a[ind_ord2_2]

    ind_ord3 = ~ind_a0

    Eigens[ind_ord0, :] = np.column_stack((nn[ind_ord0], nn[ind_ord0], nn[ind_ord0]))
    Eigens[ind_ord1, :] = np.column_stack((-d[ind_ord1] / c[ind_ord1], nn[ind_ord1], nn[ind_ord1]))
    Eigens[ind_ord2_1, :] = -0.5 * np.column_stack((b[ind_ord2_1] / a[ind_ord2_1], b[ind_ord2_1] / a[ind_ord2_1], nn[ind_ord2_1]))
    Eigens[ind_ord2_2, :] = np.column_stack(((-b[ind_ord2_2] + sqdelta[ind_ord2_2]) * inva, (-b[ind_ord2_2] - sqdelta[ind_ord2_2]) * inva, nn[ind_ord2_2]))

    # Replace solveOrd3Eigens with the equivalent Python function or library call
    # Eigens[ind_ord3, :] = solveOrd3Eigens(a[ind_ord3], b[ind_ord3], c[ind_ord3], d[ind_ord3])

    Em = np.min(Eigens, axis=1)  # Equivalent to MATLAB's min(Eigens,[],2)

    import pdb; pdb.set_trace()


    # Reshape each set of eigenvalues into the 3D format
    for i in range(3):
        EV[:, :, i] = Eigens[:, i].reshape((h, w))




    # Assuming Eigens is a 2D array of shape (h*w, 3) where each row represents the three eigenvalues of a pixel
    # Reshape the eigenvalues to match the image dimensions
    Em = Eigens[:, 0].reshape((h, w))

    # Calculate the eigenvectors
    num = Em * Sgygz + Sgxgy * Sgxgz - Sgxgx * Sgygz
    den = Em**2 - Em * Sgxgx - Em * Sgygy - Sgxgy**2 + Sgxgx * Sgygy
    den[den == 0] = np.finfo(float).eps  # Avoid division by zero

    v3 = np.ones_like(Em)
    v2 = num / den
    v1 = (Sgxgy * v2 + Sgxgz) / (Em - Sgxgx)
    v1[Em == Sgxgx] = 0  # Handle special cases

    # Normalize the eigenvectors
    Ampl = np.sqrt(v1**2 + v2**2 + v3**2)
    FV = np.zeros((h, w, 3))
    FV[:, :, 0] = v1 / Ampl
    FV[:, :, 1] = v2 / Ampl
    FV[:, :, 2] = v3 / Ampl

    # Ensure the eigenvectors are only set for positive definite cases
    posdef = np.sum(Eigens < 0, axis=1).reshape((h, w)) == 0
    FV *= posdef[..., np.newaxis]

    import pdb; pdb.set_trace()
    return FV

    
    









import scipy.ndimage as ndi
def sobel_gradients(image):
    # Compute gradients using Sobel operator
    grad_x = ndi.sobel(image, axis=2, mode='constant')
    grad_y = ndi.sobel(image, axis=1, mode='constant')
    grad_z = ndi.sobel(image, axis=0, mode='constant')
    return grad_x, grad_y, grad_z

import numpy as np
from scipy.ndimage import gaussian_filter
def compute_gradients(array):
    # Compute gradients along each dimension
    grad_x = np.gradient(array, axis=2)
    grad_y = np.gradient(array, axis=1)
    grad_z = np.gradient(array, axis=0)
    
    #grad_z = -grad_z  # Invert z gradients to match the coordinate system
    
    return grad_x, grad_y, grad_z

def structure_tensor_3d_test(volume, sigma=1):
    # Compute gradients
    grad_x = np.gradient(volume, axis=2)
    grad_y = np.gradient(volume, axis=1)
    grad_z = np.gradient(volume, axis=0)
    
    grad_z = -grad_z  # Invert z gradients to match the coordinate system
    
    # grad_x, grad_y, grad_z = sobel_gradients(volume)


    # Initialize structure tensor components
    Jxx = grad_x * grad_x
    Jxy = grad_x * grad_y
    Jxz = grad_x * grad_z
    Jyy = grad_y * grad_y
    Jyz = grad_y * grad_z
    Jzz = grad_z * grad_z

    if sigma > 0:
        # Apply Gaussian filter
        Jxx = gaussian_filter(Jxx, sigma)
        Jxy = gaussian_filter(Jxy, sigma)
        Jxz = gaussian_filter(Jxz, sigma)
        Jyy = gaussian_filter(Jyy, sigma)
        Jyz = gaussian_filter(Jyz, sigma)
        Jzz = gaussian_filter(Jzz, sigma)


    J = np.array([Jxx, Jyy, Jzz, Jxy, Jxz, Jyz])

    # # Form the structure tensor
    # J = np.array([[Jxx, Jxy, Jxz],
    #               [Jxy, Jyy, Jyz],
    #               [Jxz, Jyz, Jzz]])

    return J

def eigen_decomposition(tensor):
    eigenvalues, eigenvectors = np.linalg.eigh(tensor)
    return eigenvalues, eigenvectors





# def eigen_decomposition(tensor):


#     tensor_reshaped = np.zeros((
#         tensor.shape[1], 
#         tensor.shape[2],
#         tensor.shape[3],
#         3, 
#         3))

#     # Assigning the elements to the 3x3 matrices
#     tensor_reshaped[..., 0, 0] = tensor[0]  # s_xx
#     tensor_reshaped[..., 1, 1] = tensor[1]  # s_yy
#     tensor_reshaped[..., 2, 2] = tensor[2]  # s_zz
#     tensor_reshaped[..., 0, 1] = tensor_reshaped[..., 1, 0] = tensor[3]  # s_xy
#     tensor_reshaped[..., 0, 2] = tensor_reshaped[..., 2, 0] = tensor[4]  # s_xz
#     tensor_reshaped[..., 1, 2] = tensor_reshaped[..., 2, 1] = tensor[5]  # s_yz

#     # Compute eigenvalues and eigenvectors for each 3x3 matrix
#     eigenvalues, eigenvectors = np.linalg.eigh(tensor_reshaped)

#     # Find the index of the smallest eigenvalue at each grid point
#     indices_of_smallest_eigenvalues = np.argmin(eigenvalues, axis=-1)

#     # Select the eigenvector corresponding to each smallest eigenvalue
#     smallest_eigenvectors = np.array([eigenvectors[i, j, k, :, indices_of_smallest_eigenvalues[i, j, k]] 
#                                     for i in range(tensor_reshaped.shape[0])
#                                     for j in range(tensor_reshaped.shape[1])
#                                     for k in range(tensor_reshaped.shape[2])])

#     # Reshape smallest_eigenvectors back to the grid shape
#     smallest_eigenvectors = smallest_eigenvectors.reshape(
#         tensor.shape[1], 
#         tensor.shape[2],
#         tensor.shape[3],
#         3)
    
#     return eigenvalues, smallest_eigenvectors


# def calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU):

#     s = structure_tensor_3d(volume, sigma=3)
    
    
#     val, vec = eig_special_3d(s)
#     # val, vec = eigen_decomposition(s)
    

#     plt.imshow(vec[2,4, :, :])
#     plt.show()
    
#     sys.exit()
    
#     import pdb; pdb.set_trace()
    

#     return s, vec, val



# def calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU):
#     """
#     Calculates the structure tensor for the given volume data.

#     Parameters:
#     - volume: numpy.ndarray
#         The volume data.
#     - SIGMA: float
#         Sigma value for structure tensor calculation.
#     - RHO: float
#         Rho value for structure tensor calculation.
#     - USE_GPU: bool
#         Flag to use GPU if available.

#     Returns:
#     - tuple: (s, vec, val)
#         The structure tensor (s), eigenvectors (vec), and eigenvalues (val).
#     """
    
#     from skimage import data, feature   
#     sigma = (1, 1, 1)
#     A_elems  = feature.structure_tensor(volume, sigma=sigma)
    
#     import pdb; pdb.set_trace()

#     #####################################################################
#     # We can then compute the eigenvalues of the structure tensor.

#     eigen = feature.structure_tensor_eigenvalues(A_elems)
#     print(eigen.shape)

#     #####################################################################
#     # Where is the largest eigenvalue?
#     coords = np.unravel_index(eigen.argmax(), eigen.shape)
#     assert coords[0] == 0  # by definition
#     coords
        
    
#     # import pdb; pdb.set_trace()


#     return s, vec, val











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
    
    # flag_GPU=True
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


    # S = structure_tensor_3d(volume,SIGMA,RHO)
    
    
    # # import scipy.io
    # # S = S[:,3,:,:]
    # # S= np.reshape(S, (6,-1))
    # # S = np.transpose(S, (1,0))
    # # print(S.shape)
    # # scipy.io.savemat('S.mat', {'S': S})
    # # sys.exit()
    
    
    # # from skimage import data, feature   
    # # sigma = (1, 1, 1)
    # # S  = feature.structure_tensor(volume, sigma=sigma)
    
    
    # # S = structure_tensor_3d_test(volume, sigma=3)
    # # S = structure_tensor_3d_test2(volume)
    
    # # reshape_eigenvalues
    
    # val, vec = eig_special_3d(S)
    

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
    
    
    print(rotation_matrix)
    
    
    # posdef = (vec_2D_reshaped.reshape(vector_field_slice.shape)[0,:] < 0) == 0
    # plt.imshow(posdef)
    # plt.figure(2)
    # posdef = (rotated_vecs[0, :] < 0) == 0
    # plt.imshow(posdef)
    # # plt.figure(3)
    # # plt.imshow(vector_field_slice[0,:])
    # # plt.figure(4)
    # # plt.imshow(rotated_vecs[0, :])
    # plt.show()
    
    
    
    # print(rotation_matrix)
    # posdef = (vector_field_slice[1,:] < 0) == 0
    # plt.imshow(posdef)
    # plt.figure(2)
    # posdef = (rotated_vecs[1, :] < 0) == 0
    # plt.imshow(posdef)
    # plt.show()
    
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
    
    orig_map=plt.get_cmap('hsv')
    print("\nPlotting images...")
    plt.imshow(img_helix, cmap=orig_map)                      
    plt.show()

    
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
    ax[0,0].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
    ax[0,0].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
    tmp = ax[0,1].imshow(img_helix,cmap=orig_map)
    ax[1,0].imshow(img_intrusion,cmap=orig_map)  
    tmp2 = ax[1,1].imshow(img_FA,cmap=plt.get_cmap('inferno'))   
    fig.colorbar(tmp) 
                             
    plt.show()


    


def write_images(img_helix, img_intrusion, img_FA, start_index, OUTPUT_DIR, z):
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

    print(f"Saving images: {start_index + z}")

    # cv2.imwrite(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.tif", img_helix)
    # cv2.imwrite(f"{OUTPUT_DIR}/IA/IA_{(start_index + z):05d}.tif", img_intrusion)
    # cv2.imwrite(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):05d}.tif", img_FA)
    

    def write_img_rgb(img, output_path, cmap=plt.get_cmap('hsv')):
        minimum = np.nanmin(img)
        maximum = np.nanmax(img)
        img = (img + np.abs(minimum)) * (1 / (maximum - minimum)) 
        img = cmap(img)
        img = (img[:, :, :3] * 255).astype(np.uint8)
        
        cv2.imwrite(output_path, img)
    
    write_img_rgb(img_helix,f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.tif", cmap=plt.get_cmap('hsv'))
    write_img_rgb(img_intrusion,f"{OUTPUT_DIR}/IA/IA_{(start_index + z):05d}.tif", cmap=plt.get_cmap('hsv'))
    write_img_rgb(img_FA,f"{OUTPUT_DIR}/FA/FA_{(start_index + z):05d}.tif", cmap=plt.get_cmap('inferno'))



def main():

    print("FINISH ! ")
    
if __name__ == '__main__':  
    main()