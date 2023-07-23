"""
3D_Data_Processing
Version 6.0
Date: 31/08/2022
Created by Simon Burt in accordance with the MSc Individual Project
Key feature of this version - Ability to read nifti files for processing
"""
import time
import sys
import os
import math
import numpy as np
from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename, askdirectory
import glob 
import multiprocessing

import skimage.io
import skimage.transform
import multiprocessing as mp
from tqdm import tqdm
from structure_tensor import parallel_structure_tensor_analysis, eig_special_3d, structure_tensor_3d
import matplotlib.pyplot as plt
import tifffile

import dask_image
import dask_image.imread

import cv2

#test for GPU
flag_GPU = 0
try:
    import cupy as cp
except:
    flag_GPU = 0
#import 3D structure tensor and utils (helper functions)
if not(flag_GPU):
    from structure_tensor import eig_special_3d, structure_tensor_3d    #CPU version
else:
    from structure_tensor.cp import eig_special_3d, structure_tensor_3d #GPU version      
import utilsST_3D













def main():



    #Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    #folder_path = askdirectory(initialdir=os.getcwd(), title="Select folder") # show an "Open" dialog box and return the path to the selected file
    
    folder_path = 'W:/LADAF-2021-17/heart/19.85um_complete-organ_bm18/158.8um_LADAF-2021-17_heart_pag-0.04_0.12_'
    folder_path = 'W:/LADAF-2021-17/heart/19.85um_complete-organ_bm18/test'
    
    folder_path = '/data/projects/hop/data_repository/LADAF-2021-17/heart/19.85um_complete-organ_bm18/158.8um_LADAF-2021-17_heart_pag-0.04_0.12_'
    folder_path = '/data/projects/hop/data_repository/LADAF-2021-17/heart/19.85um_complete-organ_bm18/79.40um_LADAF-2021-17_heart_pag-0.04_0.12_'
    folder_path = '/data/projects/hop/data_repository/LADAF-2021-17/heart/19.85um_complete-organ_bm18/39.70um_LADAF-2021-17_heart_pag-0.04_0.12_'
    
    folder_path = '/data/projects/hop/data_repository/LADAF-2021-17/heart/19.85um_complete-organ_bm18/19.85um_LADAF-2021-17_heart_pag-0.04_0.12_jp2__tif'

    # folder_path = '/data/projects/hop/data_repository/biology/pig/19.59um_complete-organ/39.18um_pig-fresh_heart_pag-0.07_0.07_'
    #folder_path = '/data/projects/hop/data_repository/biology/pig/19.59um_complete-organ/19.59um_pig-fresh_heart_pag-0.07_0.07_'


    folder_path = "/data/projects/md1290/EDF/8um_LADAF-2021-17_heart/.ZZ_UNDISTORTED_D0.0025_nobackup/PROCESSING/voltif/32.04um_LADAF-2021-17_heart-c__pag-0.03_0.07_"
    #folder_path = "/data/projects/md1290/EDF/8um_LADAF-2021-17_heart/.ZZ_UNDISTORTED_D0.0025_nobackup/PROCESSING/voltif/16.02um_LADAF-2021-17_heart-c__pag-0.03_0.07_"
    #folder_path = "/data/projects/md1290/EDF/8um_LADAF-2021-17_heart/.ZZ_UNDISTORTED_D0.0025_nobackup/PROCESSING/voltif/8.01um_LADAF-2021-17_heart-c__pag-0.03_0.07_"
    

    
    flag_bin_mask = False
    path_bin_mask = '/data/bm18/inhouse/JOSEPH/python/orientation/bin_mask_39um.png'
    path_bin_mask_heart = '/data/bm18/inhouse/JOSEPH/python/orientation/bin_mask_39um_heart.png'

    
    crop_x = 0
    crop_y = 0
    
    
    sigma = 1 # noise scale
    rho = 3 # integration scale


    
    
    

    save_as_tvk = False    # True - if want to save results as tvk and nifti file (program takes longer to run)

    
    #Check number of files to convert 
    file_list = glob.glob(folder_path+"/*.tif*")
    file_list=sorted(file_list)
    if not file_list:
        print("No files were found (check your folder path)")
        sys.exit('No files were found (check your folder path)')
    else:
        N = len(file_list)
        print(str(N)+" tif files found\n")
        
    
    crop_z = True
    # Reads tiff files    
    img_array_dask = dask_image.imread.imread(f'{folder_path}/*.tif')
    
    print(img_array_dask)
    
    if crop_z:
        N = rho*10+1
        #volume = img_array_dask[int(len(file_list)/1.57):int(len(file_list)/1.57+N),:,:].compute()
        N_img = img_array_dask.shape[0]/1.68
        print(N_img)
        
        volume = img_array_dask[int(N_img-N/2):int(N_img+N/2),
                                int(img_array_dask.shape[1]*crop_x):int(img_array_dask.shape[1]*(1-crop_x)),
                                int(img_array_dask.shape[2]*crop_y):int(img_array_dask.shape[2]*(1-crop_y))].compute()
        #volume = img_array_dask[int(11180-15):int(11180+15),:,:].compute()

    else:
        volume = img_array_dask[:,:,:].compute()
    
    volume = volume.astype('float32')
    print(volume.shape)

    
    
    # #volume = volume.transpose(1, 2, 0)   # reshape to (x,y,z)
    # print(volume.shape)


    # skimage.io.imshow(volume[0, :, :])    # z-axis
    # plt.show()
        
    #skimage.io.imsave('test.tif',volume[:, :, 0])   
    
    # check sigma and rho is no more than 10% of data size
    smallest_size = min(volume.shape[0], volume.shape[1], volume.shape[2])
    # check sigma value is smaller than 10% of image size
    if not smallest_size * 0.1 >= sigma:
        sys.exit('sigma value must be smaller than 10% of volume size in all planes')
    # check rho value is smaller than 10% of image size
    if not smallest_size * 0.1 >= rho:
        sys.exit('rho value must be smaller than 10% of volume size in all planes')

    # plot 2d slice of data in 3-planes
    
    plot= False
    if plot:
        show_slice= True
        if show_slice:
            fig, axes = plt.subplots(1, 3, figsize=(8, 4))
            fig.suptitle('2d slice of data in 3-planes')
            ax = axes.ravel()
            ax[0].imshow(volume[int(volume.shape[0]/2), :, :], cmap=plt.cm.gray)
            ax[1].imshow(volume[:, int(volume.shape[1]/2), :], cmap=plt.cm.gray)
            ax[2].imshow(volume[:, :, int(volume.shape[2]/2)], cmap=plt.cm.gray)
            fig.tight_layout()
            plt.show()



    
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # ask user location to save results
    save_filename = "'U:/inhouse/JOSEPH/python/joseph/orientation/test'"
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if save_filename == '':  # asksaveasfile return `none` if dialog closed with "cancel".
        sys.exit('no file name or location was selected to save results')

    # check to see if data needs splitting into chunks for processing
    if volume.size >= 160000000000:
        # split data into smaller chunks for processing + add padding to avoid boundary errors
        chunked_data = True
        split_amount = math.ceil(volume.size / 16000000)  # number of chunks to split the data into
        print('\033[93m'+ f'Warning: Large dataset loaded. Data will be split into {split_amount} chunks for processing.' + '\033[0m')
        volume_chunks, chunk_amount_z, split_axis = chunk_split(volume, split_amount, rho)    # calls function to split data into manageable chunks
    else:
        print('No chuncked')
        chunked_data = False
        split_amount = 1
        order = 0
        chunk_amount_z = volume.shape[2]
        split_axis = 0


    # start pool for multicore processing
    if mp.cpu_count()>=60:
        pool = mp.Pool(mp.cpu_count()-10)   # error occurs if pc has more than 60-cores
    else:
        pool = mp.Pool(mp.cpu_count()-5)

    # process data
    final_data = [None] * split_amount  # create empty list
    for inds in range(split_amount):
        print(inds)
        if chunked_data:
            chunk = volume_chunks[inds]['data']
            order = volume_chunks[inds]['order']
            padding_start = volume_chunks[inds]['padding_start']
            padding_end = volume_chunks[inds]['padding_end']
        else:
            chunk = volume.astype('float32')

        # calculate structure tensor,eigenvalues and eigenvectors
        print(f'calculating structure tensor for chunk {order+1}...')
        t1 = time.perf_counter()  # start time
        
        
        flag_GPU=0
        if flag_GPU == 1:
            s, vec, val = parallel_structure_tensor_analysis(chunk, 
                                                            sigma, 
                                                            rho, 
                                                            devices=4*['cuda:0'] +4*['cuda:1']+64*['cpu'], 
                                                            block_size=100,
                                                            structure_tensor=True)                                                   
        else:
            s, vec, val = parallel_structure_tensor_analysis(chunk, sigma, rho,
                                                         structure_tensor=True) # vec has shape =(3,x,y,z) in the order of (z,y,x)
        # s = structure_tensor_3d(chunk, sigma, rho)
        # val,vec = eig_special_3d(s)
        

        # remove padding from st, e'values & e'vectors if data is chunked
        if chunked_data:
            array_end = vec.shape[split_axis+1] - padding_end
            s = array_slice(s, split_axis+1, padding_start, array_end)
            vec = array_slice(vec, split_axis+1, padding_start, array_end)
            val = array_slice(val, split_axis+1, padding_start, array_end)

        t2 = time.perf_counter()  # stop time
        print(f'finished calculating structure tensors in {t2 - t1} seconds')
        
        
        
                
        
        
        
        
        show_slice= False
        if show_slice:
            import utilsST_3D_2
    
            #volume = volume.transpose(1, 2, 0)   # reshape to (x,y,z)
    
            utilsST_3D_2.show_vol_flow(volume, vec[0:2], s=20, double_arrow = True) 
            
            
            # utilsST_3D_2.show_vol_orientation(volume, 
            #                                 vec, 
            #                                 coloring = utilsST_3D_2.fan_coloring)

        
        

            
        
        def get_click_point(img):
            # Open the image
            img_vmin,img_vmax = np.percentile(img, (5, 95))
            plt.imshow(img, vmin=img_vmin, vmax=img_vmax, cmap=plt.cm.gray)
            plt.title('CHOOSE THE CENTER OF LV')

            # Wait for the user to click on the image
            point = plt.ginput(1, timeout=-1)
            
            # Close the image
            plt.close()
            
            # Return the point that was clicked
            return point[0]
        
        

        
        def calculate_angles(v,n):
            # Define the normal vector of the plane p
            p = np.array([0, 0, 1])
            
            # Calculate the dot product between each vector in v and the normal vector p
            dot_products = np.sum(v * p.reshape(3, 1, 1), axis=0)
            
            # Calculate the magnitudes of each vector in v and the normal vector p
            v_magnitudes = np.sqrt(np.sum(v**2, axis=0))
            p_magnitude = np.linalg.norm(p)
            
            # Calculate the angles between each vector in v and the normal vector p
            angles = np.arccos(dot_products / (v_magnitudes * p_magnitude))* 180 / np.pi
            
            angles = 90 - angles
            # Calculate the signed angles by checking the z component of each vector in v
            signed_angles = np.where(v[2] > 90, angles - 180, angles) 
            
        
            # c = np.cross(n,np.array([0,0,1]),axis=0)
            # dot_product = np.sum(v * c, axis=0)
            
            # # Calculate the sign of each scalar in b
            # signs = np.sign(dot_product)
            # signs = np.where(signs==0, 1, signs)
            
            
            # # Multiply a by signs elementwise
            # signed_angles = signed_angles * signs
            
            return signed_angles



        def calculate_helix_angle(img,vec_2D):
            
            rows, cols = img.shape
            
            
            
            #center_point=np.array([100,500])

            # #Vector axis inversed ! x,y and not row,column

            # vec_2D[:,100,510]=np.array([1, 2, 3])
            # vec_2D[:,100,490]=np.array([1,2, 3])

            # normal_matrix[:,100,490]
            # normal_matrix[:,100,510]
            # circ_matrix[:,100,490]
            # circ_matrix[:,100,510]
            
            

            y, x = np.meshgrid(np.arange(rows), np.arange(cols))
            normal_matrix = np.stack((y - center_point[1], x - center_point[0], np.zeros_like(x)), axis=0)
            normal_matrix = normal_matrix / np.linalg.norm(normal_matrix, axis=0, keepdims=True)
            
            v_proj = np.einsum('ijk,ijk->jk', vec_2D, normal_matrix) / np.einsum('ijk,ijk->jk', normal_matrix, normal_matrix) * normal_matrix
            v_orth_proj = vec_2D - v_proj
            magnitudes = np.linalg.norm(v_orth_proj, axis=0, keepdims=True)

            v_orth_proj = v_orth_proj / magnitudes
            
            
            
            circ_matrix = np.cross(normal_matrix,np.array([0,0,1]),axis=0)
            dot_product = np.sum(v_orth_proj * circ_matrix, axis=0)
           
            # Calculate the sign of each scalar in b
            signs = np.sign(dot_product)
            signs = np.where(signs==0, 1, signs)
           
           
            # Multiply a by signs elementwise
            v_orth_proj = v_orth_proj * signs
           

            # v_orth_proj[:,100,510]
            # v_orth_proj[:,100,490]

            

            
            
    


    
            # np.rad2deg(np.arccos(np.clip(np.dot(v, np.array([1,0])), -1.0, 1.0)))  
            #helix_angle = np.rad2deg(angle_between(proj_vec, circ_vec))
            
            #helix_angle = np.rad2deg(signed_angle(v, np.array([1,0])))
            
            

            
            # v = v_orth_proj
            # n = normal_matrix
            signed_angles = calculate_angles(v_orth_proj,normal_matrix)
            
            return signed_angles
            
            
            
        z = int(volume.shape[0]/2)
        img = volume[z, :, :]
        center_point = get_click_point(img)
        print(center_point)
        center_point = np.asarray(center_point)
        plt.close()
        
        img = volume[z, :, :]
        vec_2D = vec[:,z, :, :]
        img_helix = calculate_helix_angle(img,vec_2D)
        
        
    
        if flag_bin_mask:
            bin_mask = cv2.imread(path_bin_mask, cv2.IMREAD_GRAYSCALE) 
            bin_mask = bin_mask[int(bin_mask.shape[0]*crop_x):int(bin_mask.shape[0]*(1-crop_x)),
                                int(bin_mask.shape[1]*crop_y):int(bin_mask.shape[1]*(1-crop_y))]
            print('bin_mask: ',bin_mask.shape)
            print('helix img: ',img_helix.shape)
            img_helix = img_helix * bin_mask
            #img_helix = cv2.bitwise_and(img_helix, bin_mask)

        img_vmin,img_vmax = np.percentile(img, (5, 95))
        # plot 2d slice of data in 3-planes
        show_slice= True
        if show_slice:
            orig_map=plt.get_cmap('hsv')
            reversed_map = orig_map.reversed()
            fig, axes = plt.subplots(1, 2, figsize=(8, 4))
            ax = axes.ravel()
            ax[0].imshow(img, vmin=img_vmin, vmax=img_vmax ,cmap=plt.cm.gray)      
            tmp = ax[1].imshow(img_helix,cmap=reversed_map)                
            fig.colorbar(tmp) 
            #plt.savefig('test1.png', dpi=1200)           
            plt.show()
            
            
            
            #img_vmin,img_vmax = np.percentile(img_helix, (0.5, 99.5))
            #std_img=np.interp(std_img, (std_img.min(), std_img.max()), (img_vmin, img_vmax))
            # img_helix = np.interp(img_helix, (img_vmin, img_vmax), (0, 1))
            img_helix = np.interp(img_helix, (np.nanmin(img_helix), np.nanmax(img_helix)), (0, 1))



            
            # #plt.imsave('test2.png', img, cmap=plt.cm.gray)
            # tifffile.imsave('test2.png', img)

            # plt.imsave('test3.png', img_helix, cmap=reversed_map)
            plt.figure()
            plt.imshow(img_helix, cmap=reversed_map)
            plt.show()
            
            
        #img_helix = img_helix.astype('uint16')  

        #tifffile.imsave('orien.tif', img_helix)
 





        sys.exit()







        # import SimpleITK as sitk
        # vol = sitk.GetImageFromArray(volume)
        # sitk.Show(vol, '3')


                
        
        
        # call function to calculate coherency
        print(f'calculating coherence values for chunk {order+1} ...')
        t1 = time.perf_counter()  # start time
        val_chunks = np.array_split(val, val.shape[3], axis=3)      # split array into chunks to be parallel processed
        coh = list(tqdm(pool.imap(coherence, [chunk_coh for chunk_coh in val_chunks]),
                        total=len(val_chunks)))      # calculates coherency from eigenvalues using parallel processing
        
        # coh = coherence(val)
        # coh.shape

        coh_result = np.concatenate(coh, axis=2)    # joins array chunks back into one
        coh_result[coh_result > 1] = 1
        coh_result[coh_result < 0] = 0

        del coh     # saves some memory
        t2 = time.perf_counter()        # stop time
        print(f'finished calculating coherence in {t2 - t1} seconds')

        tifffile.imsave('coh.tif', coh_result[int(coh_result.shape[0]/2), :, :])

        # plot 2d slice of data in 3-planes
        show_slice= True
        if show_slice:
            fig, axes = plt.subplots(1, 3, figsize=(8, 4))
            fig.suptitle('Coherence (2d slice of data in 3-planes)')
            ax = axes.ravel()
            ax[0].imshow(coh_result[int(coh_result.shape[0]/2), :, :], cmap=plt.cm.gray)
            ax[1].imshow(coh_result[:, int(coh_result.shape[1]/2), :], cmap=plt.cm.gray)
            ax[2].imshow(coh_result[:, :, int(coh_result.shape[2]/2)], cmap=plt.cm.gray)
            fig.tight_layout()
            plt.show(block=False)
            plt.pause(0.001)


        #mask = np.zeros_like(img_helix)
        img_helix[coh_result[z,:,:]<0.1] = None
        # plot 2d slice of data in 3-planes
        show_slice= True
        if show_slice:
            fig, axes = plt.subplots(1, 2, figsize=(8, 4))
            ax = axes.ravel()
            ax[0].imshow(img, cmap=plt.cm.gray) 
            cmap = plt.get_cmap('hsv')
            cmap.set_bad(color='black')
            tmp = ax[1].imshow(img_helix,cmap=plt.get_cmap('hsv'))                
            fig.colorbar(tmp)            
            plt.show()
        






        # call function to calculate fractional anisotropy (fa)
        print(f'calculating fa values for chunk {order + 1} ...')
        t1 = time.perf_counter()  # start time
        fa = list(tqdm(pool.imap(fractional_anisotropy, [chunk_coh for chunk_coh in val_chunks]),
                        total=len(val_chunks)))  # calculates coherency from eigenvalues using parallel processing
        fa_result = np.concatenate(fa, axis=2)  # joins array chunks back into one
        del fa      # saves some memory
        t2 = time.perf_counter()  # stop time
        print(f'finished calculating fa in {t2 - t1} seconds')
        
        
        tifffile.imsave('fa.tif', fa_result[int(coh_result.shape[0]/2), :, :])

        # plot 2d slice of data in 3-planes
        show_slice= True
        if show_slice:
            fig, axes = plt.subplots(1, 3, figsize=(8, 4))
            fig.suptitle('fractional anisotropy (2d slice of data in 3-planes)')
            ax = axes.ravel()
            ax[0].imshow(fa_result[int(fa_result.shape[0]/2), :, :], cmap=plt.cm.gray)
            ax[1].imshow(fa_result[:, int(fa_result.shape[1]/2), :], cmap=plt.cm.gray)
            ax[2].imshow(fa_result[:, :, int(fa_result.shape[2]/2)], cmap=plt.cm.gray)
            fig.tight_layout()
            plt.show()


        sys.exit()

    # write results to file
    print('saving results to file...')
    
    print(final_data)    
    

    pool.close()
    pool.join()
    print('finished - program end')
    
    
    
    
    


def coherence(val):
    # function calculates the coherence value of every voxel
    # Based on the following formula:
    # c_s = (3 * eigval3) / (eigval1 ** 2 + eigval2 ** 2 + eigval3 ** 2) ** 0.5

    # Flatten S.
    input_shape = val.shape
    val = val.reshape(3, -1)
    c_s = np.empty((1,) + val.shape[1:], dtype="float32")
    c_a = np.empty((1,) + val.shape[1:], dtype="float32")
    tmp = np.empty((4,) + val.shape[1:], dtype="float32")

    # compute c_s
    np.multiply(val[2], 3, out=tmp[0])         # 3 * eigval3
    np.multiply(val[0], val[0], out=tmp[1])   # eigval1^2
    np.multiply(val[1], val[1], out=tmp[2])   # eigval2^2
    np.multiply(val[2], val[2], out=tmp[3])   # eigval3^2
    a = np.add(tmp[1], tmp[2])
    a += tmp[3]
    a **= 0.5
    np.divide(tmp[0], a, out=c_s, where=a != 0)
    # compute c_a
    np.subtract(1, c_s, out=c_a)
    del tmp
    return c_a.reshape(input_shape[1:])

def fractional_anisotropy(val):
    # function calculates fractional anisotropy
    # Based on the following formula:
    # fa = (1/2)**0.5 * ((eigval1 - eigval2)**2 + (eigval2-eigval3)**2 + (eigval3 - eigval1)**2)** 0.5/(eigval1 ** 2 + eigval2 ** 2 + eigval3 ** 2)) ** 0.5

    # Flatten val
    input_shape = val.shape
    val = val.reshape(3, -1)
    fa = np.empty((1,) + val.shape[1:], dtype="float32")
    tmp = np.empty((5,) + val.shape[1:], dtype="float32")

    # compute fa
    np.multiply(val[0], val[0], out=tmp[0])   # eigval1^2
    np.multiply(val[1], val[1], out=tmp[1])   # eigval2^2
    np.multiply(val[2], val[2], out=tmp[2])   # eigval3^2
    np.add(tmp[0], tmp[1], out=tmp[3])
    tmp[3] += tmp[2]

    a = np.subtract(val[0], val[1])
    a **= 2
    b = np.subtract(val[1], val[2])
    b **= 2
    c = np.subtract(val[2], val[0])
    c **= 2
    np.add(a, b, out=tmp[4])
    tmp[4] += c

    d = np.divide(tmp[4], tmp[3], where=tmp[3] != 0)
    d **= 0.5
    e = 1/(2**0.5)

    np.multiply(e, d, out=fa)
    del tmp
    return fa.reshape(input_shape[1:])


def colouring(vec, fa):
    # function to assign rgba value to data
    vec = abs(vec)  # negative numbers don't matter for assigning colour
    vec[~np.isfinite(vec)] = np.nan  # Replaces any inf values with nan

    input_shape = vec.shape
    vec = vec.reshape(3, -1)
    tmp = np.empty((6,) + vec.shape[1:], dtype="float32")
    # find red
    np.subtract(vec[2], np.nanmin(vec[2]), out=tmp[0])
    np.subtract(np.nanmax(vec[2]), np.nanmin(vec[2]), out=tmp[1])
    r = np.divide(tmp[0], tmp[1])  # x-axis component is the red channel
    # find blue
    np.subtract(vec[1], np.nanmin(vec[1]), out=tmp[2])
    np.subtract(np.nanmax(vec[1]), np.nanmin(vec[1]), out=tmp[3])
    b = np.divide(tmp[2], tmp[3])  # x-axis component is the red channel
    # find green
    np.subtract(vec[0], np.nanmin(vec[0]), out=tmp[4])
    np.subtract(np.nanmax(vec[0]), np.nanmin(vec[0]), out=tmp[5])
    g = np.divide(tmp[4], tmp[5])  # x-axis component is the red channel
    a = fa  # use fa as alpha channel

    rgba = np.empty((input_shape[1], input_shape[2], input_shape[3], 4))
    rgba[..., 0] = r.reshape(input_shape[1:])
    rgba[..., 1] = g.reshape(input_shape[1:])
    rgba[..., 2] = b.reshape(input_shape[1:])
    rgba[..., 3] = a
    # rgba *= 255  # normalise data
    rgba = rgba.astype(np.uint8)  # convert to uint8 to reduce file size
    del tmp
    return rgba

def array_slice(a, axis, start, end, step=1):
    return a[(slice(None),) * (axis % a.ndim) + (slice(start, end, step),)]

def chunk_split(volume, split_amount, padding):
    # function to split large datasets into smaller manageable chunks for processing
    # input:
    #       volume = numpy array of the dataset to be split into chunks.
    #       split_amount = integer of the number of chunks to split the data into.
    #       padding = integer of the number of layers to add to each chunk. takes the data from the neighbouring chunk.
    #                 if first or last chunk, padding will only be added to one side of the chunk.
    # output:
    #       volume_chunk = list containing multiple numpy arrays of the chunks.
    #       chunk_amount_z = tuple containing the z-dimension of the every chunk the data was split into.

    # Determine largest axis and split data along it
    # if volume.shape[0] >= volume.shape[1] and volume.shape[0] >= volume.shape[2]:
    #     split_axis = 0
    # elif volume.shape[1] >= volume.shape[0] and volume.shape[1] >= volume.shape[2]:
    #     split_axis = 1
    # else:
    #     split_axis = 2
    
    split_axis = 0

    chunk_amount_z = tuple([volume.shape[split_axis] // split_amount + int(x < volume.shape[split_axis] % split_amount) for x in range(split_amount)])
    volume_chunks = [None] * split_amount  # initialise variable
    for iz in range(split_amount):
        if iz == 0:
            # deals with initial chunk
            volume_chunks[iz] = {
                                'data': array_slice(volume, split_axis, 0, (chunk_amount_z[iz] + padding)),
                                'order': int(iz),
                                'padding_start': 0,
                                'padding_end': padding
                                }
            previous_finish = chunk_amount_z[iz]
        elif iz == (split_amount - 1):
            # deals with last chunk
            volume_chunks[iz] = {
                                'data': array_slice(volume, split_axis, (previous_finish - padding), (previous_finish + chunk_amount_z[iz])),
                                'order': int(iz),
                                'padding_start': padding,
                                'padding_end': 0
                                }
        else:
            # deals with all other chunks
            volume_chunks[iz] = {
                                'data': array_slice(volume, split_axis, (previous_finish - padding),(previous_finish + chunk_amount_z[iz] + padding)),
                                'order': int(iz),
                                'padding_start': padding,
                                'padding_end': padding
                                }
            previous_finish = previous_finish + chunk_amount_z[iz]
    return(volume_chunks, chunk_amount_z, split_axis)








# def main(folder_path,sigma,rho):
    
    
if __name__ == '__main__':  
    main()
    
    
    
