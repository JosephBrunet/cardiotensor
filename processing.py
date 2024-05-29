import os
import sys
import time
import math
import glob
import argparse
import numpy as np
import cv2
import dask
import dask_image.imread
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from memory_profiler import profile
from skimage.filters import unsharp_mask
from distutils.util import strtobool
import warnings

# Optional GPU support
try:
    import cupy as cp
    USE_GPU = True
except ImportError:
    USE_GPU = False


from structure_tensor import parallel_structure_tensor_analysis




def convert_to_8bit(img, minimum, maximum):
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
    
    minimum = int(minimum)
    maximum = int(maximum)
    
    img_normalized = (img + np.abs(minimum)) * (255 / (maximum - minimum))  # Normalize the image to the range [0, 255]
    img_8bit = img_normalized.astype('uint8')
    return img_8bit


def read_parameter_file(para_file_path):
    # Define expected parameters with default values (if any)
    params = {
        'IMAGES_PATH': None,
        'OUTPUT_PATH': None,
        'SIGMA': None,
        'RHO': None,
        'POINT_MITRAL_VALVE': None,
        'POINT_APEX': None,
        'TEST': False
    }

    with open(para_file_path, 'r') as file:
        for line in file:
            line = line.split('#')[0].strip().replace(' ','')  # Remove comments and whitespace
            
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip().upper()  # Normalize key
                
                # Process specific parameters
                if key in ['SIGMA', 'RHO']:
                    params[key] = float(value)
                elif key in ['POINT_MITRAL_VALVE', 'POINT_APEX']:
                    params[key] = np.fromstring(value, dtype=int, sep=',')
                elif key == 'TEST':
                    params[key] = bool(strtobool(value))
                else:
                    params[key] = value.strip()

                print(f"{key} found")

    # Validate required parameters
    for param, value in params.items():
        if value is None:
            raise ValueError(f"Missing parameter: {param}")

    return [params[key] for key in ['IMAGES_PATH', 'OUTPUT_PATH', 'SIGMA', 'RHO', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'TEST']]




def get_file_list(volume_path):
    """
    Retrieves a sorted list of file paths from a specified directory.

    Parameters:
    - volume_path: str
        Path to the directory containing the files.

    Returns:
    - file_type: str
        The type of files found (tif or jp2).
    - file_list: list
        A sorted list of file paths.
    """
    # Check number of files to convert 
    file_list_tif = sorted(glob.glob(os.path.join(volume_path, "*.tif*")))
    file_list_jp2 = sorted(glob.glob(os.path.join(volume_path, "*.jp2*")))

    if file_list_tif and file_list_jp2:
        raise ValueError('Both tif and jp2 files were found. Check your folder path.')
    elif file_list_tif:
        file_list = file_list_tif      
        file_type = 'tif'
    elif file_list_jp2:
        file_list = file_list_jp2
        file_type = 'jp2'
    else:  
        sys.exit('No files were found (check your folder path)')
    
    print(f"File type: {file_type}")        
    return file_type, file_list


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


def calculate_center_vector(pt_mv, pt_apex):
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
    center_vec = pt_apex - pt_mv
    
    center_vec = pt_mv - pt_apex
    center_vec = center_vec / np.linalg.norm(center_vec)
    
    center_vec[2] = -center_vec[2]  # Invert z components (Don't know why but it works)

    if center_vec[2] < 0:
        center_vec = -center_vec
        
    return center_vec



def adjust_start_end_index(start_index, end_index, N_img, padding_start, padding_end, is_test):
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
        test_index = int(N_img / 1.68)
        start_index_padded = max(test_index - padding_start, 0)
        end_index_padded = min(test_index + 1 + padding_end, N_img)
    else:
        # Adjust start and end indices considering padding
        start_index_padded = max(start_index - padding_start, 0)
        end_index_padded = min(end_index + padding_end, N_img)
    
    print(f"Start index padded, End index padded : {start_index_padded}, {end_index_padded}")
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
        image_data = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)
        
        # Add unsharp mask filter
        # image_data = unsharp_mask(image_data, radius=3, amount=0.6)
        
        return image_data
    
    # Create a list of delayed tasks to read each image
    delayed_tasks = [dask.delayed(lambda x: np.array(custom_image_reader(x)))(file_path) for file_path in file_list]
    
    # volume = volume_dask[start_index_padded:end_index_padded,:,:].compute()
    volume = np.array(dask.compute(*delayed_tasks[start_index_padded:end_index_padded]))
    return volume.astype('float32')




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
    if USE_GPU:
        s, vec, val = parallel_structure_tensor_analysis(volume, 
                                                        SIGMA, 
                                                        RHO, 
                                                        devices=['cuda:0'] +['cuda:1']+64*['cpu'], 
                                                        block_size=200,
                                                        structure_tensor=True)                                                   
    else:
        s, vec, val = parallel_structure_tensor_analysis(volume, SIGMA, RHO,
                                                        structure_tensor=True) # vec has shape =(3,x,y,z) in the order of (z,y,x)

    return s, vec, val



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


def compute_helix_and_transverse_angles(img, vector_field_2d, center_point):
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
    rows, cols = img.shape

    reshaped_vector_field  = np.reshape(vector_field_2d, (3, -1))

    center_x, center_y  = center[0], center[1]
    # X, Y = np.meshgrid(np.arange(1, cols + 1) - center_x, np.arange(1, rows + 1) - center_y )
    X, Y = np.meshgrid(np.arange(cols) - center_x, np.arange(rows) - center_y)
    gamma = -np.arctan2(Y.flatten(), X.flatten())
    cos_angle  = np.cos(gamma)
    sin_angle = np.sin(gamma)

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
    
    # FI_unwrapped = np.arctan2(reshaped_rotated_vector_field [2, :, :], reshaped_rotated_vector_field [1, :, :])
    # ind_wrapped = (FI_unwrapped - FI) != 0
    # VF_ind = np.tile(ind_wrapped[np.newaxis, :, :], (3, 1, 1))

    # vector_field_2d[VF_ind] *= -1

    helix_angle = np.rad2deg(-helix_angle)
    transverse_angle = np.rad2deg(-transverse_angle)

    return helix_angle, transverse_angle


def rotate_vectors_to_new_axis(vector_field_3d, new_axis_vec):
    """
    Rotates vectors to align with a new axis.

    Parameters:
    - vector_field_2d (numpy.ndarray): Array of vectors to be rotated.
    - new_axis_vec (numpy.ndarray): The new axis to align the vectors with.

    Returns:
    - numpy.ndarray: Rotated vectors aligned with the new axis.
    """
    # Ensure new_axis_vec is normalized
    new_axis_vec = new_axis_vec / np.linalg.norm(new_axis_vec)
    
    # Calculate the rotation matrix
    vec1 = np.array([0, 0, 1])  # Initial vertical axis
    
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (new_axis_vec / np.linalg.norm(new_axis_vec)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    if np.any(kmat):
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    else:
        rotation_matrix = np.eye(3)
        
    # Reshape vec_2D to (3, 100) for matrix multiplication
    vec_2D_reshaped = np.reshape(vector_field_3d, (3, -1))
    # Rotate the vectors
    rotated_vecs = np.dot(rotation_matrix, vec_2D_reshaped)
    
    # Reshape back to the original shape (3, 10, 10)
    rotated_vecs = rotated_vecs.reshape(vector_field_3d.shape)
    
    return rotated_vecs       



# @profile  
def process_3d_data(para_file_path, start_index=0, end_index=0, IS_TEST=False):
    # function to process data
    
    print(para_file_path)
    import numpy as np

    
    print(f'Start index, End index : {start_index}, {end_index}')
    
    VOLUME_PATH, OUTPUT_DIR, SIGMA, RHO, PT_MV, PT_APEX, IS_TEST = read_parameter_file(para_file_path)
    print(f'PARAMETERS : ',VOLUME_PATH, OUTPUT_DIR, SIGMA, RHO, PT_MV, PT_APEX)
    
    
    
    file_type, file_list = get_file_list(VOLUME_PATH)
    N_img = len(file_list)
    print(f"{N_img} {file_type} files found\n")  
    
    print('Reading images with Dask...')      
    volume_dask = dask_image.imread.imread(f'{VOLUME_PATH}/*.{file_type}')
    print('Finished reading images with Dask\n')
    print(f"Dask volume: {volume_dask}")
    
    center_line = interpolate_points(PT_MV, PT_APEX, N_img)
    center_vec = calculate_center_vector(PT_MV, PT_APEX)
    
    padding_start = math.ceil(RHO)
    if padding_start > start_index:
        padding_start = start_index
    padding_end = math.ceil(RHO)
    if padding_end > (N_img - end_index):
        padding_end = N_img - end_index
    
    start_index_padded, end_index_padded = adjust_start_end_index(start_index, end_index, N_img, padding_start, padding_end, IS_TEST)
    print(f"Start index padded, End index padded : {start_index_padded}, {end_index_padded}")
    
    volume = load_volume(file_list, start_index_padded, end_index_padded)
    print(f"Loaded volume shape {volume.shape}")    
    print('Dask size used by chunk: ',volume.size*4*12.5/1000000000, ' GB') # x4 because float32 / x 12.5 for structure tensore calculation 
     
     
    print(f'\ncalculating structure tensor')
    t1 = time.perf_counter()  # start time

    s, vec, val = calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU)
    print(vec.shape)

    volume, s, vec, val = remove_padding(volume, s, vec, val, padding_start, padding_end)
    print(vec.shape)

    center_line = center_line[start_index_padded:end_index_padded]
    
    t2 = time.perf_counter()  # stop time
    print(f'finished calculating structure tensors in {t2 - t1} seconds')

    print(vec.shape)

    print(f'\nCalculating helix/intrusion angle and fractional anisotropy:')
    for z in range(vec.shape[1]):
        print(f"Processing image: {start_index + z}")
        
        if os.path.exists(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.tif") and IS_TEST == False:
            print(f"File {(start_index + z):05d} already exist")
            continue
        
        center_point = np.around(center_line[z])        
        img = volume[z, :, :]
        vec_2D = vec[:,z, :, :]
        val_2D = val[:,z, :, :]
        

        img_FA = calculate_fraction_anisotropy(val_2D)

        vec_2D = rotate_vectors_to_new_axis(vec_2D, center_vec)     
        
        img_helix,img_intrusion = compute_helix_and_transverse_angles(img,vec_2D,center_point)
        
  

        if IS_TEST:
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

            tmp = ax[0,1].imshow(img_helix,cmap=reversed_map)
            tmp = ax[1,0].imshow(img_intrusion,cmap=reversed_map)                                
            tmp2 = ax[1,1].imshow(img_FA,cmap=plt.get_cmap('inferno'))                                
            fig.colorbar(tmp) 
            plt.show()
            

            # plt.imshow(img,vmin=img_vmin, vmax=img_vmax ,cmap=plt.cm.gray) 
            # plt.figure(2)
            # plt.imshow(img_helix,cmap=reversed_map) 
            # plt.figure(3)
            # plt.imshow(img_intrusion,cmap=reversed_map)   
            # plt.figure(4)
            # plt.imshow(img_FA_8bits,cmap=plt.get_cmap('inferno'))                                 
            # plt.show()

        
        #img_vmin,img_vmax = np.percentile(img_helix, (0.5, 99.5))
        #std_img=np.interp(std_img, (std_img.min(), std_img.max()), (img_vmin, img_vmax))
        # img_helix = np.interp(img_helix, (img_vmin, img_vmax), (0, 1))
        
        # img_helix = np.interp(img_helix, (np.nanmin(img_helix), np.nanmax(img_helix)), (0, 1))


        if not IS_TEST:
                
            # Convert the float64 image to int8
            # img_helix = img_helix.astype(np.float32)
            # img_intrusion = img_intrusion.astype(np.float32)
            # img_FA = img_FA.astype(np.float32)
         
            img_helix = convert_to_8bit(img_helix, -90, 90)
            img_intrusion = convert_to_8bit(img_intrusion, -90, 90)
            img_FA = convert_to_8bit(img_FA, 0, 1)
            
            
            os.makedirs(OUTPUT_DIR + '/HA',exist_ok=True)
            os.makedirs(OUTPUT_DIR + '/IA',exist_ok=True)
            os.makedirs(OUTPUT_DIR + '/FA',exist_ok=True)        

            # plt.imsave(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.png", img_helix, cmap=reversed_map)
            # plt.imsave(f"{OUTPUT_DIR}/IA/IA_{(start_index + z):05d}.png", img_intrusion, cmap=reversed_map)
            # plt.imsave(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):05d}.png", img_FA, cmap=plt.get_cmap('inferno'))

            print(f"Saving images: {start_index + z}")
            print(img_helix.dtype)

            cv2.imwrite(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.tif", img_helix)
            cv2.imwrite(f"{OUTPUT_DIR}/IA/IA_{(start_index + z):05d}.tif", img_intrusion)
            cv2.imwrite(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):05d}.tif", img_FA)

    
    print('finished - program end')





def main():

    if len(sys.argv) < 2:
        Tk().withdraw()
        para_file_path = askopenfilename(initialdir=f"{os.getcwd()}/param_files", title="Select file") # show an "Open" dialog box and return the path to the selected file
        # para_file_path = "/data/bm18/inhouse/JOSEPH/python/orientation_heart/parameters_317.6um_LADAF-2021-17_heart.txt"
        if not para_file_path:
            sys.exit("No file selected!")
        start_index = 0
        end_index = 0
        
    elif len(sys.argv) >= 2:
        parser = argparse.ArgumentParser(description='Convert images between tif and jpeg2000 formats.')
        parser.add_argument('para_file_path', type=str, help='Path to the input text file.')
        parser.add_argument('--start_index', type=int, default=0, help='Starting index for processing.')
        parser.add_argument('--end_index', type=int, default=0, help='Ending index for processing.')
        args = parser.parse_args()
        para_file_path = args.para_file_path
        start_index = args.start_index
        end_index = args.end_index
    
    start_time = time.time()
    process_3d_data(para_file_path, start_index, end_index) 
    print("--- %s seconds ---" % (time.time() - start_time))

    print("FINISH ! ")
    
if __name__ == '__main__':  
    main()