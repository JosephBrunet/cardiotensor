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
from skimage.filters import unsharp_mask
from distutils.util import strtobool
import warnings
import numpy as np

# Optional GPU support
try:
    import cupy as cp
    USE_GPU = True
except ImportError:
    USE_GPU = False
USE_GPU = False

from MyoTensor.utils import *




# @profile  
def process_3d_data(conf_file_path, start_index=0, end_index=0, IS_TEST=False):
    # function to process data
    
    print(f"\n---------------------------------")
    print(f'READING PARAMETER FILE : {conf_file_path}\n')
    
    print(f'Start index, End index : {start_index}, {end_index}\n')
    
    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(e)
        sys.exit(f'⚠️  Error reading parameter file: {conf_file_path}')
    
    
    VOLUME_PATH, MASK_PATH, IS_FLIP, OUTPUT_DIR, OUTPUT_TYPE, SIGMA, RHO, N_CHUNK, PT_MV, PT_APEX, IS_TEST, N_SLICE_TEST = [params[key] for key in ['IMAGES_PATH', 'MASK_PATH', 'FLIP', 'OUTPUT_PATH', 'OUTPUT_TYPE', 'SIGMA', 'RHO', 'N_CHUNK', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'TEST', 'N_SLICE_TEST']]

    if not IS_TEST:
        is_already_done = True
        if end_index == 0:
            is_already_done = False
        for idx in range(start_index, end_index):
            if not os.path.exists(f"{OUTPUT_DIR}/HA/HA_{(idx):06d}.tif"):
                is_already_done = False
                break

        if is_already_done:            
            print('All images are already done')
            return
         
    if MASK_PATH:
        is_mask = True
    else:
        print("Mask path not provided")
        is_mask = False
        
    import pdb; pdb.set_trace()
       
        
    def read_raw_volume(file_path, shape, dtype):
        """
        Reads a raw volume file and returns it as a numpy array.

        Parameters:
        - file_path (str): Path to the raw volume file.
        - shape (tuple of ints): Dimensions of the volume, e.g., (z, y, x).
        - dtype (numpy dtype): Data type of the volume, e.g., np.uint8, np.float32.

        Returns:
        - volume (numpy.ndarray): The volume as a numpy array.
        """
        # Calculate the total number of elements
        num_elements = np.prod(shape)
        
        # Read the binary data from the file
        with open(file_path, 'rb') as file:
            data = np.fromfile(file, dtype=dtype, count=num_elements)
        
        # Reshape the flat array to the desired 3D shape
        volume = data.reshape(shape)
        
        return volume
import numpy as np

# Define the dimensions and data type of the volume
width, height, depth = 512, 512, 256  # Example dimensions
data_type = np.uint8  # Example data type

# Path to the .raw file
file_path = 'path/to/your/file.raw'

# Read the binary data from the .raw file
volume_data = np.fromfile(file_path, dtype=data_type)

# Reshape the data to the volume dimensions
volume_data = volume_data.reshape((depth, height, width))

# Verify the shape
print(volume_data.shape)




        
    print(f"\n---------------------------------")
    print(f"READING VOLUME INFORMATION\n")
    print(f"Volume path: {VOLUME_PATH}")
    img_list, img_type = get_image_list(VOLUME_PATH)
    N_img = len(img_list)
    print(f"{N_img} {img_type} files found\n")  
    
    print('Reading images with Dask...')      
    volume_dask = dask_image.imread.imread(f'{VOLUME_PATH}/*.{img_type}')
    print(f"Dask volume: {volume_dask}")
    print('\nAll information gathered')

    
    if is_mask:
        print(f"\n---------------------------------")
        print(f"READING MASK INFORMATION FROM {MASK_PATH}...\n")
        mask_list, mask_type = get_image_list(MASK_PATH)
        print(f"{len(mask_list)} {mask_type} files found\n")  
        # mask_dask = dask_image.imread.imread(f'{MASK_PATH}/*.{mask_type}')
        # print(f"Dask mask: {mask_dask}")
        N_mask = len(mask_list)
        binning_factor = N_img / N_mask
        print(f"Mask bining factor: {binning_factor}\n")
        print('\nAll information gathered')

    
    print(f"\n---------------------------------")
    print("CALCULATE CENTER LINE AND CENTER VECTOR\n")
    center_line = interpolate_points(PT_MV, PT_APEX, N_img)
    center_vec = calculate_center_vector(PT_MV, PT_APEX, IS_FLIP)
    print(f"Center vector: {center_vec}")
    
    print(f"\n---------------------------------")
    print("CALCULATE PADDING START AND ENDING INDEXES\n")
    padding_start = math.ceil(RHO)
    padding_end = math.ceil(RHO)
    if not IS_TEST:
        if padding_start > start_index:
            padding_start = start_index
        if padding_end > (N_img - end_index):
            padding_end = N_img - end_index
            
    print(f"Padding start, Padding end : {padding_start}, {padding_end}")
    start_index_padded, end_index_padded = adjust_start_end_index(start_index, end_index, N_img, padding_start, padding_end, IS_TEST, N_SLICE_TEST)
    print(f"Start index padded, End index padded : {start_index_padded}, {end_index_padded}")
    

    
    print(f"\n---------------------------------")
    print("LOAD VOLUMES\n")
    volume = load_volume(img_list, start_index_padded, end_index_padded).astype('float32')
    print(f"Loaded volume shape {volume.shape}")    

    if is_mask:
        start_index_padded_mask = int((start_index_padded / binning_factor)-1)#-binning_factor/2)
        if start_index_padded_mask < 0:
            start_index_padded_mask = 0
        end_index_padded_mask = int((end_index_padded / binning_factor)+1)#+binning_factor/2)
        if end_index_padded_mask > N_mask:
            end_index_padded_mask = N_mask
        
        
        print(f"Mask start index padded: {start_index_padded_mask} - Mask end index padded : {end_index_padded_mask}")
        
        mask = load_volume(mask_list, start_index_padded_mask, end_index_padded_mask)
        print(f"Mask volume loaded of shape {mask.shape}")
        
        
        from scipy.ndimage import zoom
        mask = zoom(mask, binning_factor, order=0)
        
        
        # start_index_mask_upscaled = int(mask.shape[0]/2)-int(volume.shape[0]/2)
        # end_index_mask_upscaled = start_index_mask_upscaled+volume.shape[0]
            
        start_index_mask_upscaled = int(np.abs(start_index_padded_mask * binning_factor - start_index_padded))
        end_index_mask_upscaled = start_index_mask_upscaled + volume.shape[0]
        
        print(volume.shape)
        print(mask.shape)
        print(f"Start index: {start_index_mask_upscaled}, End index: {end_index_mask_upscaled}")

        
        if start_index_mask_upscaled < 0:
            start_index_mask_upscaled = 0
        if end_index_mask_upscaled > mask.shape[0]:
            end_index_mask_upscaled = mask.shape[0]
        
        mask = mask[start_index_mask_upscaled:end_index_mask_upscaled,:]
        
        print(f"Start index: {start_index_mask_upscaled}, End index: {end_index_mask_upscaled}")
        print(mask.shape)

        
        mask_resized = np.empty_like(volume)
        for i in range(mask.shape[0]):
            # Resize the slice to match the corresponding slice of the volume
            mask_resized[i] = cv2.resize(mask[i], (volume.shape[2], volume.shape[1]), interpolation = cv2.INTER_LINEAR)
            
            kernel_size = kernel_size = RHO * 2
            # Ensure kernel_size is odd
            if kernel_size % 2 == 0:
                kernel_size += 1
            kernel_size = int(kernel_size)
            
            assert kernel_size % 2 == 1, "Kernel size has to be an odd number"      

        #     kernel = np.ones((kernel_size,kernel_size),np.uint8)
        #     mask_resized[i] = cv2.dilate(mask_resized[i], kernel, iterations = 1)
        # print(f"Applying a dilation to mask with kernel = [{kernel_size},{kernel_size}]")

        mask = mask_resized
        
        assert mask.shape == volume.shape, f"Mask shape {mask.shape} does not match volume shape {volume.shape}"           

        volume[mask == 0] = np.nan
        
        print(f"Mask applied to image volume")  
        del mask_resized    
        del mask

    # num_cpus= 64
    # aa = ["num_cpus*['cpu']", 
    #       "32*['cuda:0'] + 32*['cuda:1'] + num_cpus*['cpu']", 
    #       "16*['cuda:0'] + 16*['cuda:1'] + num_cpus*['cpu']",
    #       "8*['cuda:0'] + 8*['cuda:1'] + num_cpus*['cpu']",
    #       "8*['cuda:0'] + num_cpus*['cpu']",
    #       "16*['cuda:0'] + num_cpus*['cpu']",
    #       "32*['cuda:0'] + num_cpus*['cpu']"]
    # bb = [200,300,400]
    # for a in aa:
    #     for b in bb:
            
    #         print(f"\n---------------------------------")
    #         print(f'CALCULATING STRUCTURE TENSOR')
    #         t1 = time.perf_counter()  # start time
    #         s, val, vec  = calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU,device=eval(a), block_size=b)    
    #         t2 = time.perf_counter()  # stop time
    #         print(f'{t2 - t1} seconds - device: {a}, block_size: {b}')
    #         time.sleep(2)
            
            
    # sys.exit()
    
    
    



    print(f"\n---------------------------------")
    print(f'CALCULATING STRUCTURE TENSOR')
    t1 = time.perf_counter()  # start time
    s, val, vec  = calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU)    
    print(f"Vector shape: {vec.shape}")
    volume, _, val, vec = remove_padding(volume, s, val, vec, padding_start, padding_end)
    del s
    print(f"Vector shape after removing padding: {vec.shape}")
    
    # print( np.isnan(vec[1,0,:,:]).all())
    # sys.exit()

    center_line = center_line[start_index_padded:end_index_padded]  # adjust the center line to the new start and end index
    
    
    # # plt.imshow(vec[0,0, :, :]);plt.show()
    # import pdb; pdb.set_trace()
    
    
    # Putting all the vectors in positive direction
    posdef = np.all(val >= 0, axis=0)  # Check if all elements are non-negative along the first axis
    vec_norm = np.linalg.norm(vec, axis=0)  # Compute the norm of the vectors along the first axis

    # Use broadcasting to normalize the vectors
    # vec = vec / (vec_norm * posdef)
    
    vec = vec / (vec_norm)


    # Check for negative z component and flip if necessary
    # negative_z = vec[2, :] < 0
    # vec[:, negative_z] *= -1
    


    t2 = time.perf_counter()  # stop time
    print(f'finished calculating structure tensors in {t2 - t1} seconds')

    IS_NP_SAVE = False
    if IS_NP_SAVE:
        os.makedirs(OUTPUT_DIR,exist_ok=True)
        np.save(f"{OUTPUT_DIR}/eigen_vec.npy", vec)


    print(f'\nCalculating helix/intrusion angle and fractional anisotropy:')
    for z in range(vec.shape[1]):
        print(f"Processing image: {start_index + z}")
        
        if os.path.exists(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):06d}.tif") and IS_TEST == False:
            print(f"File {(start_index + z):06d} already exist")
            continue
        
        center_point = np.around(center_line[z])        
        img = volume[z, :, :]
        vector_field_slice = vec[:,z, :, :]
        
        eigen_val_slice = val[:,z, :, :]
        
        img_FA = calculate_fraction_anisotropy(eigen_val_slice)

        vector_field_slice = rotate_vectors_to_new_axis(vector_field_slice, center_vec) 
        

        img_helix,img_intrusion = compute_helix_and_transverse_angles(vector_field_slice,center_point)
          


                
        if IS_TEST:
            plot_images(img, img_helix, img_intrusion, img_FA, center_point, PT_MV, PT_APEX)
        if not IS_TEST:           
            write_images(img_helix, img_intrusion, img_FA, start_index, OUTPUT_DIR, OUTPUT_TYPE, z)

    return

def main():

    if len(sys.argv) < 2:
        Tk().withdraw()
        conf_file_path = askopenfilename(initialdir=f"{os.getcwd()}/param_files", title="Select file") # show an "Open" dialog box and return the path to the selected file
        if not conf_file_path:
            sys.exit("No file selected!")
        start_index = 0
        end_index = 0
        
    elif len(sys.argv) >= 2:
        parser = argparse.ArgumentParser(description='Convert images between tif and jpeg2000 formats.')
        parser.add_argument('conf_file_path', type=str, help='Path to the input text file.')
        parser.add_argument('--start_index', type=int, default=0, help='Starting index for processing.')
        parser.add_argument('--end_index', type=int, default=0, help='Ending index for processing.')
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
        start_index = args.start_index
        end_index = args.end_index
    
    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        sys.exit(f'⚠️  Error reading parameter file: {conf_file_path}')
    
    VOLUME_PATH, MASK_PATH, IS_FLIP, OUTPUT_DIR, OUTPUT_TYPE, SIGMA, RHO, N_CHUNK, PT_MV, PT_APEX, IS_TEST, N_SLICE_TEST = [params[key] for key in ['IMAGES_PATH', 'MASK_PATH', 'FLIP', 'OUTPUT_PATH', 'OUTPUT_TYPE', 'SIGMA', 'RHO', 'N_CHUNK', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'TEST', 'N_SLICE_TEST']]
    
    img_list, img_type = get_image_list(VOLUME_PATH)
    
    # Set end_index to total_images if it's zero
    if end_index == 0:
        end_index = len(img_list)
    
    if not IS_TEST:
        for idx in range(start_index, end_index, N_CHUNK):
            print(f"Processing slices {idx} to {idx+N_CHUNK}")
            start_time = time.time()
            process_3d_data(conf_file_path, idx, idx + N_CHUNK) 
            print("--- %s seconds ---" % (time.time() - start_time))
            
            if time.time() - start_time > 5:
                sys.exit()

    else:
        start_time = time.time()
        process_3d_data(conf_file_path, start_index, end_index) 
        print("--- %s seconds ---" % (time.time() - start_time))

    print("FINISH ! ")
    
if __name__ == '__main__':  
    main()