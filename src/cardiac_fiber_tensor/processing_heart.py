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

from structure_tensor import parallel_structure_tensor_analysis

from utils import *




# @profile  
def process_3d_data(para_file_path, start_index=0, end_index=0, IS_TEST=False):
    # function to process data
    
    print(f"\n---------------------------------")
    print(f'READING PARAMETER FILE : {para_file_path}\n')

    
    print(f'Start index, End index : {start_index}, {end_index}\n')
    
    params = read_parameter_file(para_file_path)
    VOLUME_PATH, MASK_PATH, OUTPUT_DIR, SIGMA, RHO, PT_MV, PT_APEX, IS_FLIP, IS_TEST = [params[key] for key in ['IMAGES_PATH', 'MASK_PATH', 'OUTPUT_PATH', 'SIGMA', 'RHO', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'FLIP', 'TEST']]

    if not IS_TEST:
        is_already_done = True
        if end_index == 0:
            is_already_done = False
        for idx in range(start_index, end_index):
            if not os.path.exists(f"{OUTPUT_DIR}/HA/HA_{(start_index + idx):05d}.tif"):
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
        
    print(f"\n---------------------------------")
    print(f"READING VOLUME INFORMATION FROM {VOLUME_PATH}...\n")
    img_list, img_type = get_image_list(VOLUME_PATH)
    N_img = len(img_list)
    print(f"{N_img} {img_type} files found\n")  
    
    print('Reading images with Dask...')      
    volume_dask = dask_image.imread.imread(f'{VOLUME_PATH}/*.{img_type}')
    print(f"Dask volume: {volume_dask}")
    print('\nAll information gathered\n')

    
    if is_mask:
        print(f"\n---------------------------------")
        print(f"READING MASK INFORMATION FROM {MASK_PATH}...\n")
        mask_type, mask_list = get_image_list(MASK_PATH)
        print(f"{len(mask_list)} {mask_type} files found\n")  
        mask_dask = dask_image.imread.imread(f'{MASK_PATH}/*.{mask_type}')
        print(f"Dask mask: {mask_dask}")
        N_mask = len(mask_list)
        binning_factor = N_img / N_mask
        print(f"Mask bining factor: {binning_factor}\n")
        print('\nAll information gathered\n')

    
    print(f"\n---------------------------------")
    print("CALCULATE CENTER LINE AND CENTER VECTOR\n")
    center_line = interpolate_points(PT_MV, PT_APEX, N_img)
    center_vec = calculate_center_vector(PT_MV, PT_APEX, IS_FLIP)
    print(f"\nCenter vector: {center_vec}")
    
    print(f"\n---------------------------------")
    print("CALCULATE PADDING START AND ENDING INDEXES\n")
    padding_start = math.ceil(RHO)
    padding_end = math.ceil(RHO)
    if padding_start > start_index:
        padding_start = start_index
    if padding_end > (N_img - end_index):
        padding_end = N_img - end_index
    start_index_padded, end_index_padded = adjust_start_end_index(start_index, end_index, N_img, padding_start, padding_end, IS_TEST)
    print(f"Start index padded, End index padded : {start_index_padded}, {end_index_padded}")
    

    
    print(f"\n---------------------------------")
    print("LOAD VOLUMES\n")
    volume = load_volume(img_list, start_index_padded, end_index_padded).astype('float64')
    print(f"Loaded volume shape {volume.shape}")    
    print('Dask size used by chunk: ',volume.size*4*12.5/1000000000, ' GB') # x4 because float32 / x 12.5 for structure tensore calculation 

    if is_mask:
        start_index_padded_mask = int((start_index_padded / binning_factor)-binning_factor/2)
        if start_index_padded_mask < 0:
            start_index_padded_mask = 0
        end_index_padded_mask = int((end_index_padded / binning_factor)+binning_factor/2)
        if end_index_padded_mask > N_mask:
            end_index_padded_mask = N_mask
        
        
        print(f"Mask start index padded, Mask end index padded : {start_index_padded_mask}, {end_index_padded_mask}")
        
        mask = load_volume(mask_list, start_index_padded_mask, end_index_padded_mask)
        print(f"Mask volume loaded of shape {mask.shape}")
        
        
        from scipy.ndimage import zoom
        mask = zoom(mask, binning_factor, order=0)
        
        # start_index_mask_upscaled = int(mask.shape[0]/2)-int(volume.shape[0]/2)
        # end_index_mask_upscaled = start_index_mask_upscaled+volume.shape[0]
            
        start_index_mask_upscaled = int(np.abs(start_index_padded_mask*binning_factor - start_index_padded))
        end_index_mask_upscaled = start_index_mask_upscaled+volume.shape[0]
        
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
            
            kernel = np.ones((2,2),np.uint8)
            mask_resized[i] = cv2.dilate(mask_resized[i], kernel, iterations = 1)

        mask = mask_resized
        
        assert mask.shape == volume.shape, f"Mask shape {mask.shape} does not match volume shape {volume.shape}"

        volume[mask == 0] = np.nan
        
        print(f"Mask applied to image volume")  
        del mask_resized    
        del mask





    print(f"\n---------------------------------")
    print(f'CALCULATING STRUCTURE TENSOR')
    t1 = time.perf_counter()  # start time
    s, vec, val = calculate_structure_tensor(volume, SIGMA, RHO, USE_GPU)    
    print(vec.shape)
    volume, _, vec, val = remove_padding(volume, s, vec, val, padding_start, padding_end)
    del s
    print(vec.shape)
    

    center_line = center_line[start_index_padded:end_index_padded]  # adjust the center line to the new start and end index
    
    posdef = np.sum(val < 0, axis=0) == 0
    vec_norm = np.linalg.norm(vec,axis=0)
    for i in range(3):
        vec[i] = vec[i] / (vec_norm*posdef)
        
    #vec[2,:] = -vec[2,:]
    
    negative_z = vec[2, :] < 0
    vec[:, negative_z] *= -1

    t2 = time.perf_counter()  # stop time
    print(f'finished calculating structure tensors in {t2 - t1} seconds')

    IS_NP_SAVE = True
    if IS_NP_SAVE:
        os.makedirs(OUTPUT_DIR,exist_ok=True)
        np.save(f"{OUTPUT_DIR}/eigen_vec.npy", vec)




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
        
        print(f"Center point: {center_point}")
        print(f"Center vector: {center_vec}")
        

        # plt.imshow(vec_2D[2,:])
        # plt.figure(2)
        # plt.imshow(vec_2D[0,:])
        # plt.show()


        img_FA = calculate_fraction_anisotropy(val_2D)



        # print("\nPlotting hivkkihgivhv...")
        # img_vmin,img_vmax = np.nanpercentile(img, (5, 95))
        # orig_map=plt.get_cmap('hsv')
        # reversed_map = orig_map.reversed()
        # fig, axes = plt.subplots(2, 2, figsize=(8, 4))
        # ax = axes
        # ax[0,0].imshow(img, vmin=img_vmin, vmax=img_vmax ,cmap=plt.cm.gray)   
        # # Draw a red point at the specified coordinate
        # x, y = center_point[0:2]
        # ax[0,0].scatter(x, y, c='red', s=50, marker='o')
        # ax[0,0].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
        # ax[0,0].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
        # tmp = ax[0,1].imshow(vec_2D[0,:,:],cmap=orig_map)
        # ax[1,0].imshow(vec_2D[1,:,:],cmap=orig_map)
        # ax[1,1].imshow(vec_2D[2,:,:],cmap=orig_map)
        # fig.colorbar(tmp) 
        # plt.show()

        # posdef = np.sum(vec_2D < 0, axis=0) == 0
        
        # plt.imshow(posdef)
        # plt.show()



        # posdef = (img_helix < 0) == 0
        # plt.imshow(posdef)
        # plt.show()

        
        # posdef = (vec_2D[0,:] < 0) == 0
        # plt.imshow(posdef)
        # plt.figure(2)
        # posdef = (vec_2D[1,:] < 0) == 0
        # plt.imshow(posdef)
        # plt.show()
        
        # center_vec = np.array([0.8,0,0.2])
        
        
        
        
        # import scipy.io
        # mat = scipy.io.loadmat('FV.mat')
        
        # vec_2D = mat["VF"]#.reshape(val_2D.shape)
        # vec_2D = vec_2D.transpose(2, 0, 1)
        


        # import scipy.io
        # scipy.io.savemat('test.mat', {'vec_2D': vec_2D})
        # sys.exit()
        
        
        
        
        vec_2D = rotate_vectors_to_new_axis(vec_2D, center_vec) 
        
        
        
        # posdef = (vec_2D[0,:] < 0) == 0
        # plt.imshow(posdef)
        # plt.figure(2)
        # posdef = (vec_2D[1,:] < 0) == 0
        # plt.imshow(posdef)
        # plt.show()
           
           
           
        img_helix,img_intrusion = compute_helix_and_transverse_angles(vec_2D,center_point)
        





        # print("\nPlotting images...")
        # img_vmin,img_vmax = np.nanpercentile(img, (5, 95))
        # orig_map=plt.get_cmap('hsv')
        # reversed_map = orig_map.reversed()
        # fig, axes = plt.subplots(2, 2, figsize=(8, 4))
        # ax = axes
        # ax[0,0].imshow(img, vmin=img_vmin, vmax=img_vmax ,cmap=plt.cm.gray)   
        # # Draw a red point at the specified coordinate
        # x, y = center_point[0:2]
        # ax[0,0].scatter(x, y, c='red', s=50, marker='o')   
                
        # tmp = ax[0,1].imshow(img_helix,cmap=orig_map)
        # ax[1,0].imshow(img_intrusion,cmap=orig_map)  
        # tmp2 = ax[1,1].imshow(img_FA,cmap=plt.get_cmap('inferno'))  
        
        # ax[0,0].scatter(x, y, c='red', s=50, marker='o')   
        # ax[0,0].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
        # ax[0,0].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
        # ax[0,1].scatter(x, y, c='red', s=50, marker='o')   
        # ax[0,1].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
        # ax[0,1].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
        # ax[1,0].scatter(x, y, c='red', s=50, marker='o')   
        # ax[1,0].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
        # ax[1,0].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
        # fig.colorbar(tmp)                     
        # plt.show()        
        
        
        
        # plt.imshow(img_helix,cmap=orig_map)
        # plt.show()
        # plt.imsave(f"testingH.jpg", img_helix, cmap=orig_map)
        
        # plt.imshow(img_intrusion,cmap=orig_map)
        # plt.show()
        # plt.imsave(f"testingI.jpg", img_intrusion, cmap=orig_map)
        
        if IS_TEST:
            plot_images(img, img_helix, img_intrusion, img_FA, center_point, PT_MV, PT_APEX)
        if not IS_TEST:           
            write_images(img_helix, img_intrusion, img_FA, start_index, OUTPUT_DIR, z)
            
            
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
    
    # start_index = 7500
    # end_index = 7520
    
    start_time = time.time()
    process_3d_data(para_file_path, start_index, end_index) 
    print("--- %s seconds ---" % (time.time() - start_time))

    print("FINISH ! ")
    
if __name__ == '__main__':  
    main()