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

from utils import *



N_SLICE = 1269



# @profile  
def process_3d_data(para_file_path, start_index=0, end_index=0, IS_TEST=False):
    # function to process data
    
    print(para_file_path)
    import numpy as np

    
    print(f'Start index, End index : {start_index}, {end_index}')
    
    params = read_parameter_file(para_file_path)
    VOLUME_PATH, OUTPUT_DIR, SIGMA, RHO, IS_TEST = [params[key] for key in ['IMAGES_PATH', 'OUTPUT_PATH', 'SIGMA', 'RHO', 'TEST']] 
    print(f'PARAMETERS : ',VOLUME_PATH, OUTPUT_DIR, SIGMA, RHO)
    
    
    
    file_type, file_list = get_file_list(VOLUME_PATH)
    N_img = len(file_list)
    print(f"{N_img} {file_type} files found\n")  
    
    print('Reading images with Dask...')      
    volume_dask = dask_image.imread.imread(f'{VOLUME_PATH}/*.{file_type}')
    print('Finished reading images with Dask\n')
    print(f"Dask volume: {volume_dask}")
    

    
    padding_start = math.ceil(RHO)
    if padding_start > start_index:
        padding_start = start_index
    padding_end = math.ceil(RHO)
    if padding_end > (N_img - end_index):
        padding_end = N_img - end_index
    
    start_index_padded, end_index_padded = adjust_start_end_index(start_index, end_index, N_img, padding_start, padding_end, IS_TEST, N_SLICE)
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
    
    
    
    t2 = time.perf_counter()  # stop time
    print(f'finished calculating structure tensors in {t2 - t1} seconds')



    IS_NP_SAVE = False
    if IS_NP_SAVE:
        os.makedirs(OUTPUT_DIR,exist_ok=True)
        np.save(f"{OUTPUT_DIR}/eigen_vec.npy", vec)



    for z in range(vec.shape[1]):
        print(f"Processing image: {start_index + z}")
        
        if os.path.exists(f"{OUTPUT_DIR}/HA/HA_{(start_index + z):05d}.tif") and IS_TEST == False:
            print(f"File {(start_index + z):05d} already exist")
            continue
        
        img = volume[z, :, :]
        vec_2D = vec[:,z, :, :]
        val_2D = val[:,z, :, :]
        
        
        n = np.array([0,0,1])
        
        def calculate_angles(vec_2D, n):
            # Calculate the magnitude of each vector in the grid
            vec_magnitudes = np.linalg.norm(vec_2D, axis=0)
            
            angles = 90 - np.abs(90 * (vec_2D[2,:,:] / vec_magnitudes))
             
            
            # # Calculate the dot product of each vector in the grid with the normal vector
            # dot_products = np.einsum('ijk, i->jk', vec_2D, n)
            
            # # Calculate the cosine of the angle between each vector and the normal vector
            # cosines = dot_products / vec_magnitudes
            
            # # Use arccos to find the angle. Clip values to the range [-1, 1] to avoid NaNs due to floating point errors.
            # angles = np.arccos(np.clip(cosines, -1, 1))
            # #angles = np.arccos(dot_product / (vec_2D_magnitude * n_magnitude))
            
            # angles = np.degrees(angles)
            
            return angles
                
        angles = calculate_angles(vec_2D, n)
        print(angles.shape)
        

        img_FA = calculate_fraction_anisotropy(val_2D)
        
        
        
        img = convert_to_8bit(img)
        
        
        if IS_TEST:
            # print("\nPlotting images...")
            # img_vmin,img_vmax = np.nanpercentile(img, (20, 80))
            # orig_map=plt.get_cmap('hsv')
            # reversed_map = orig_map.reversed()
            # fig, axes = plt.subplots(2, 2, figsize=(8, 4))
            # ax = axes
            # ax[0,0].imshow(img, vmin=img_vmin, vmax=img_vmax ,cmap=plt.cm.gray)   
            # # Draw a red point at the specified coordinate
            
            # tmp = ax[0,1].imshow(angles,cmap='inferno') 
            # fig.colorbar(tmp) 
            # tmp = ax[1,0].imshow(img_FA,cmap='inferno')                               
            # plt.show()


            print("\nPlotting images...")

            tmp = plt.imshow(angles,cmap='inferno')                                
            plt.colorbar(tmp) 
            plt.show()

            
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