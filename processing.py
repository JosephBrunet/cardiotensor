"""
3D_Data_Processing
"""
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import dask_image.imread
import dask
import glob
import time
import argparse
from distutils.util import strtobool
import warnings
import cv2
from skimage.filters import unsharp_mask

# Optional GPU support
try:
    import cupy as cp
    USE_GPU = True
except ImportError:
    USE_GPU = False

from structure_tensor import parallel_structure_tensor_analysis, eig_special_3d, structure_tensor_3d
# Import structure tensor and utility functions based on available hardware

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename, askdirectory


from memory_profiler import profile






def convert_to_8bit(img):
    """
    Convert the input image to 8-bit format.

    Parameters:
    - img: numpy.ndarray
        The input image.

    Returns:
    - img_interp_8bits: numpy.ndarray
        The 8-bit converted image.
    """
    norm = np.zeros([img.shape[0], img.shape[1]])
    img_normalized = cv2.normalize(img, norm, 0, 255, cv2.NORM_MINMAX)
    img_8bit = img_normalized.astype('uint8')
    return img_8bit





def read_parameter_file(para_file_path):
    with open(para_file_path, 'r') as file:
        for line in file:
            line = line.replace(' ','').split('#')[0].strip()
            
            if line.startswith('IMAGES_PATH'):
                images_path = line.split('=')[1]
                print("IMAGES_PATH found")
            if line.startswith('OUTPUT_PATH'):
                output_dir = line.split('=')[1]
                print("OUTPUT_PATH found")
            if line.startswith('SIGMA'):
                sigma = float(line.split('=')[1])
                print("SIGMA found")
            if line.startswith('RHO'):
                rho = float(line.split('=')[1])
                print("RHO found")
            if line.startswith('POINT_MITRAL_VALVE'):
                pt_MV = line.split('=')[1]
                pt_MV = np.fromstring(pt_MV, dtype=int, sep=',')
                print("POINT_MITRAL_VALVE found")
            if line.startswith('POINT_APEX'):
                pt_apex = line.split('=')[1]
                pt_apex = np.fromstring(pt_apex, dtype=int, sep=',')
                print("POINT_APEX found")
            if line.startswith('TEST'):
                is_test = strtobool(line.split('=')[1])
                if is_test == 0:
                    is_test = False
                elif is_test == 1:
                    is_test = True
                print("IS_TEST found")
                
    return images_path, output_dir, sigma, rho, pt_MV, pt_apex, is_test






def get_file_list(volume_path):
    # Check number of files to convert 
    file_list_tif = sorted(glob.glob(volume_path + "/*.tif*"))
    file_list_jp2 = sorted(glob.glob(volume_path + "/*.jp2*"))

    if file_list_tif and file_list_jp2:
        sys.exit('Both tif and jp2 files were found (check your folder path)')
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
    # Calculate center vector
    center_vec = pt_apex - pt_mv
    center_vec = center_vec / np.linalg.norm(center_vec)
    center_vec[0:2] = -center_vec[0:2]  # Invert x and y components (check how to correct this in the future)
    return center_vec



def adjust_start_end_index(start_index, end_index, volume_dask, padding_start, padding_end, is_test):
    if start_index < 0:
        raise ValueError("start_index must be greater than 0.")
    if end_index < 0:
        raise ValueError("end_index must be greater than 0.")
    if end_index == 0:
        end_index = volume_dask.shape[0]
    
    if is_test == True:
        n_index = int(volume_dask.shape[0] / 1.68)
        center_line = center_line[n_index:n_index + 1, :]
        start_index_padded = n_index - padding_start
        end_index_padded = n_index + 1 + padding_end
    else:
        center_line = center_line[start_index:end_index, :]
        
        if start_index > padding_start:
            start_index_padded = start_index - padding_start
        else:
            start_index_padded = 0
            padding_start = 0
        if end_index + padding_end > volume_dask.shape[0]:
            padding_end = 0
            end_index_padded = end_index
        else:
            end_index_padded = end_index + padding_end
    
    print(f"Start index padded, End index padded : {start_index_padded}, {end_index_padded}")
    return start_index_padded, end_index_padded


def load_volume(file_list, start_index_padded, end_index_padded):
    
    
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
    volume = volume.astype('float32')
    print(f"Loaded volume shape {volume.shape}")    
    print('Dask size used by chunk: ',volume.size*4*12.5/1000000000, ' GB') # x4 because float32 / x 12.5 for structure tensore calculation 
    return volume




def calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU):
    # Calculates the structure tensor, eigenvalues, and eigenvectors based on the given volume, SIGMA, RHO, and USE_GPU
    # Returns the structure tensor (s), eigenvectors (vec), and eigenvalues (val)
    
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
    # Removes the padding from the volume, structure tensor, eigenvectors, and eigenvalues based on the given padding start
    # Returns the updated volume, structure tensor, eigenvectors, and eigenvalues
    
    array_end = vec.shape[1] - padding_end
    volume = volume[padding_start:array_end,:,:]
    s = s[:,padding_start:array_end,:,:]
    vec = vec[:,padding_start:array_end,:,:]
    val = val[:,padding_start:array_end,:,:]
    
    print(padding_start,array_end)


    
    return volume, s, vec, val

    
# @profile  
def processing(para_file_path, start_index, end_index, IS_TEST=False):
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
    print(f"Dask volume: {volume_dask}")
        
    center_line = interpolate_points(PT_MV, PT_APEX, N_img)
    center_vec = calculate_center_vector(PT_MV, PT_APEX)
    
    padding_start = math.ceil(RHO)
    padding_end = math.ceil(RHO)
    
    start_index_padded, end_index_padded = adjust_start_end_index(start_index, end_index, volume_dask, padding_start, padding_end, IS_TEST)
    
    volume = load_volume(volume_dask, start_index_padded, end_index_padded)
       
    
     
    print(f'\ncalculating structure tensor')
    t1 = time.perf_counter()  # start time

    s, vec, val = calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU)
    volume, s, vec, val = remove_padding(volume, s, vec, val, padding_start)
    
    t2 = time.perf_counter()  # stop time
    print(f'finished calculating structure tensors in {t2 - t1} seconds')





    print(f'\nCalculating helix/intrusion angle and fractional anisotropy\n')
    for z in range(vec.shape[1]):
        print(f"Processing image: {start_index + z}")
        
        if os.path.exists(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.tif"):
            print(f"File {(start_index + z):05d} already exist")
            continue
        
        center_point = np.around(center_line[z])        
        img = volume[z, :, :]
        vec_2D = vec[:,z, :, :]
        val_2D = val[:,z, :, :]
        

        img_FA = calculate_fraction_anisotropy(val_2D)

        vec_2D = rotate_vectors_to_new_axis(vec_2D, center_vec)     
        
        img_helix,img_intrusion,_ = calculate_helix_angle(img,vec_2D,center_point)
        


        # Convert the float64 image to int8
        # img_helix = convert_to_8bit(img_helix)
        # img_intrusion = convert_to_8bit(img_intrusion)
         
        img_helix = img_helix.astype(np.int8)
        img_intrusion = img_intrusion.astype(np.int8)

        img_FA_8bits = convert_to_8bit(img_FA)
        
  
  

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
            tmp2 = ax[1,1].imshow(img_FA_8bits,cmap=plt.get_cmap('inferno'))                                
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
            os.makedirs(OUTPUT_DIR + '/HA',exist_ok=True)
            os.makedirs(OUTPUT_DIR + '/IA',exist_ok=True)
            os.makedirs(OUTPUT_DIR + '/FA',exist_ok=True)        

            # plt.imsave(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.png", img_helix, cmap=reversed_map)
            # plt.imsave(f"{OUTPUT_DIR}/IA/IA_{(start_index + z):05d}.png", img_intrusion, cmap=reversed_map)
            # plt.imsave(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):05d}.png", img_FA_8bits, cmap=plt.get_cmap('inferno'))

            cv2.imwrite(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.tif", img_helix)
            cv2.imwrite(f"{OUTPUT_DIR}/IA/IA_{(start_index + z):05d}.tif", img_intrusion)
            cv2.imwrite(f"{OUTPUT_DIR}/FA/FA_{(start_index + z):05d}.tif", img_FA_8bits)

        
                
    

 
 
    # import plotly.graph_objects as go

    # import pandas as pd

    # df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/streamtube-wind.csv').drop(['Unnamed: 0'],axis=1)
    
    
    
    # x, y, z = np.mgrid[0:volume.shape[2], 0:volume.shape[1], 0:volume.shape[0]]
    # x = x.flatten()
    # y = y.flatten()
    # z = z.flatten()
    
    # u = vec[2].flatten()
    # v = vec[1].flatten()
    # w = vec[0].flatten()
    
    
    
    # fig = go.Figure(data=go.Streamtube(
    #     x = x,
    #     y = y,
    #     z = z,
    #     u = u,
    #     v = v,
    #     w = w,
    #     starts = dict(
    #         x = [290,300,325,330] * 4,
    #         y = [340] * 16,
    #         # z = [340,350,360,370] * 4,
    #         z = [1] * 16,

    #     ),
    #     sizeref = 1,
    #     colorscale = 'Portland',
    #     showscale = False,
    #     maxdisplayed = 20000
    # ))

    # # fig.update_layout(
    # #     scene = dict(
    # #         aspectratio = dict(
    # #             x = 2,
    # #             y = 1,
    # #             z = 0.3
    # #         )
    # #     ),
    # #     margin = dict(
    # #         t = 20,
    # #         b = 20,
    # #         l = 20,
    # #         r = 20
    # #     )
    # # )

    # fig.show()
 
 

    
    

 
 
 
    
    print('finished - program end')


    
    
    
    
def calculate_fraction_anisotropy(val_2D):
    
    # Fraction Anisotropy Calculation
    l1 = val_2D[0, :, :]
    l2 = val_2D[1, :, :]
    l3 = val_2D[2, :, :]
    lm = (l1 + l2 + l3) / 3
    numerator = np.sqrt((l1 - lm)**2 + (l2 - lm)**2 + (l3 - lm)**2)
    denominator = np.sqrt(l1**2 + l2**2 + l3**2)
    img_FA = np.sqrt(3 / 2) * (numerator / denominator)
    
    return img_FA


def calculate_helix_angle(img, vec_2D, center_point):
    """
    Calculate the signed helix angles at each point of a 2D image using a 2D orientation vector matrix 
    and a center point defining the center of the radial coordinate system.

    Parameters:
        img (numpy.ndarray): 2D image array.
        vec_2D (numpy.ndarray): 2D vector.
        center_point (tuple): Tuple (x, y) representing the center point of the image grid.

    Returns:
        numpy.ndarray: An array containing the signed angles between the 2D vector and the normal vectors of the image grid.
    """
    # Given inputs
    center = center_point[0:2]  # Replace with actual values
    rows, cols = img.shape

    FV = vec_2D  # Replace with actual value
    FV = np.reshape(FV, (3, -1))

    VF = vec_2D

    Cx, Cy = center[0], center[1]
    X, Y = np.meshgrid(np.arange(1, cols + 1) - Cx, np.arange(1, rows + 1) - Cy)
    R = np.sqrt(X ** 2 + Y ** 2)
    gama = -np.arctan2(Y.flatten(), X.flatten())
    cg = np.cos(gama)
    sg = np.sin(gama)

    FVrot = np.copy(FV)
    FVrot[2, :] = FV[2, :]
    FVrot[0, :] = cg * FV[0, :] - sg * FV[1, :]
    FVrot[1, :] = sg * FV[0, :] + cg * FV[1, :]

    VFrot = np.zeros((3, rows, cols))

    VFrot[2, :, :] = FVrot[2, :].reshape(rows, cols)
    VFrot[0, :, :] = FVrot[0, :].reshape(rows, cols)
    VFrot[1, :, :] = FVrot[1, :].reshape(rows, cols)

    FT = np.arctan(VFrot[0, :, :] / VFrot[1, :, :])
    FI = np.arctan(VFrot[2, :, :] / VFrot[1, :, :])
    FI_unwrapped = np.arctan2(VFrot[2, :, :], VFrot[1, :, :])
    ind_wrapped = (FI_unwrapped - FI) != 0
    VF_ind = np.tile(ind_wrapped[np.newaxis, :, :], (3, 1, 1))

    VF[VF_ind] *= -1

    FI = np.rad2deg(-FI)
    FT = np.rad2deg(-FT)

    return FI, FT, VF


def rotate_vectors_to_new_axis(vec_3D, center_vec):
    # Ensure center_vec is normalized
    center_vec = center_vec / np.linalg.norm(center_vec)
    
    # Calculate the rotation matrix
    vec1 = np.array([0, 0, 1])  # Initial vertical axis
    vec2=center_vec
    
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    if np.any(kmat):
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    else:
        rotation_matrix = np.eye(3)
        
    # Reshape vec_2D to (3, 100) for matrix multiplication
    vec_2D_reshaped = np.reshape(vec_3D, (3, -1))
    # Rotate the vectors
    rotated_vecs = np.dot(rotation_matrix, vec_2D_reshaped)
    
    # Reshape back to the original shape (3, 10, 10)
    rotated_vecs = rotated_vecs.reshape(vec_3D.shape)
    
    return rotated_vecs       








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
    processing(para_file_path, start_index, end_index) 
    print("--- %s seconds ---" % (time.time() - start_time))

    print("FINISH ! ")
    
if __name__ == '__main__':  
    main()