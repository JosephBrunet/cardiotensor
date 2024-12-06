import os
import sys
import math
import numpy as np
import cv2
import matplotlib.pyplot as plt
from pathlib import Path
import math
import argparse
from skimage.measure import block_reduce


from cardiotensor.utils import read_conf_file, get_volume_shape, load_volume





def writeStructuredVTK(
    aspectRatio=[1.0, 1.0, 1.0],
    origin=[0.0, 0.0, 0.0],
    cellData={},
    pointData={},
    fileName="spam.vtk",
):
    """Write a plain text regular grid vtk from

    - 3D arrays for 3D scalar fields
    - 4D arrays for 3D vector fields

    Parameters
    ----------
        aspectRatio : size 3 array, float
            Length between two nodes in every direction `e.i.` size of a cell
            Default = [1, 1, 1]

        origin : size 3 array float
            Origin of the grid
            Default = [0, 0, 0]

        cellData : dict ``{"field1name": field1, "field2name": field2, ...}``
            Cell fields, not interpolated by paraview.
            The field values are reshaped into a flat array in the lexicographic order.
            ``field1`` and ``field2`` are ndimensional array
                (3D arrays are scalar fields and 4D array are vector valued fields).

        pointData : dict ``{"field1name": field1, "field2name": field2, ...}``
            Nodal fields, interpolated by paraview. ``pointData`` has the same shape as ``cellData``.

        fileName : string
            Name of the output file.
            Default = 'spam.vtk'

    WARNING
    -------
        This function deals with structured mesh thus ``x`` and ``z`` axis are swapped **in python**.
    """

    dimensions = []

    # Check dimensions
    if len(cellData) + len(pointData) == 0:
        print("spam.helpers.writeStructuredVTK() Empty files. Not writing {}".format(fileName))
        return 0

    if len(cellData):
        dimensionsCell = list(cellData.values())[0].shape[:3]
        for k, v in cellData.items():
            if set(dimensionsCell) != set(v.shape[:3]):
                print("spam.helpers.writeStructuredVTK() Inconsistent cell field sizes {} != {}".format(dimensionsCell, v.shape[:3]))
                return 0
        dimensions = [n + 1 for n in dimensionsCell]

    if len(pointData):
        dimensionsPoint = list(pointData.values())[0].shape[:3]
        for k, v in pointData.items():
            if set(dimensionsPoint) != set(v.shape[:3]):
                print("spam.helpers.writeStructuredVTK() Inconsistent point field sizes {} != {}".format(dimensionsPoint, v.shape[:3]))
                return 0
        dimensions = dimensionsPoint

    if len(cellData) and len(pointData):
        if {n + 1 for n in dimensionsCell} != set(dimensionsPoint):
            print(
                "spam.helpers.writeStructuredVTK() Inconsistent point VS cell field sizes.\
                 Point size should be +1 for each axis."
            )

    with open(fileName, "w") as f:
        # header
        f.write("# vtk DataFile Version 2.0\n")
        f.write("VTK file from spam: {}\n".format(fileName))
        f.write("ASCII\n\n")
        f.write("DATASET STRUCTURED_POINTS\n")

        f.write("DIMENSIONS {} {} {}\n".format(*reversed(dimensions)))
        f.write("ASPECT_RATIO {} {} {}\n".format(*reversed(aspectRatio)))
        f.write("ORIGIN {} {} {}\n\n".format(*reversed(origin)))

        # pointData
        if len(pointData) == 1:
            f.write("POINT_DATA {}\n\n".format(dimensions[0] * dimensions[1] * dimensions[2]))
            _writeFieldInVtk(pointData, f)
        elif len(pointData) > 1:
            f.write("POINT_DATA {}\n\n".format(dimensions[0] * dimensions[1] * dimensions[2]))
            for k in pointData:
                _writeFieldInVtk({k: pointData[k]}, f)

        # cellData
        if len(cellData) == 1:
            f.write("CELL_DATA {}\n\n".format((dimensions[0] - 1) * (dimensions[1] - 1) * (dimensions[2] - 1)))
            _writeFieldInVtk(cellData, f)
        elif len(cellData) > 1:
            f.write("CELL_DATA {}\n\n".format((dimensions[0] - 1) * (dimensions[1] - 1) * (dimensions[2] - 1)))
            for k in cellData:
                _writeFieldInVtk({k: cellData[k]}, f)

        f.write("\n")


def _writeFieldInVtk(data, f, flat=False):
    """
    Private helper function for writing vtk fields
    """
    
    for key in data:
        field = data[key]

        if flat:
            # SCALAR flatten (n by 1)
            if len(field.shape) == 1:
                f.write("SCALARS {} float\n".format(key.replace(" ", "_")))
                f.write("LOOKUP_TABLE default\n")
                for item in field:
                    f.write("    {}\n".format(item))
                f.write("\n")

            # VECTORS flatten (n by 3)
            elif len(field.shape) == 2 and field.shape[1] == 3:
                f.write("VECTORS {} float\n".format(key.replace(" ", "_")))
                for item in field:
                    f.write("    {} {} {}\n".format(*reversed(item)))
                f.write("\n")

        else:
            # SCALAR not flatten (n1 by n2 by n3)
            if len(field.shape) == 3:
                f.write("SCALARS {} float\n".format(key.replace(" ", "_")))
                f.write("LOOKUP_TABLE default\n")
                for item in field.reshape(-1):
                    f.write("    {}\n".format(item))
                f.write("\n")

            # VECTORS (n1 by n2 by n3 by 3)
            elif len(field.shape) == 4 and field.shape[3] == 3:                
                f.write("VECTORS {} float\n".format(key.replace(" ", "_")))
                for item in field.reshape((field.shape[0] * field.shape[1] * field.shape[2], field.shape[3])):
                    f.write("    {} {} {}\n".format(*reversed(item)))
                f.write("\n")

            # TENSORS (n1 by n2 by n3 by 3 by 3)
            elif len(field.shape) == 5 and field.shape[3] * field.shape[4] == 9:
                f.write("TENSORS {} float\n".format(key.replace(" ", "_")))
                for item in field.reshape(
                    (
                        field.shape[0] * field.shape[1] * field.shape[2],
                        field.shape[3] * field.shape[4],
                    )
                ):
                    f.write("    {} {} {}\n    {} {} {}\n    {} {} {}\n\n".format(*reversed(item)))
                f.write("\n")
            else:
                print("spam.helpers.vtkio._writeFieldInVtk(): I'm in an unknown condition!")









def vtk_writer(conf_file_path, bin_factor=1, start_index=None, end_index=None):

    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(e)
        sys.exit(f'⚠️  Error reading parameter file: {conf_file_path}')
    
    VOLUME_PATH, MASK_PATH, IS_FLIP, OUTPUT_DIR, OUTPUT_TYPE, SIGMA, RHO, N_CHUNK, PT_MV, PT_APEX, IS_TEST, N_SLICE_TEST = [params[key] for key in ['IMAGES_PATH', 'MASK_PATH', 'FLIP', 'OUTPUT_PATH', 'OUTPUT_TYPE', 'SIGMA', 'RHO', 'N_CHUNK', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'TEST', 'N_SLICE_TEST']]
    
    OUTPUT_DIR = Path(OUTPUT_DIR)
    
    w, h, N_img = get_volume_shape(VOLUME_PATH)
    
    if start_index==None:
        start_index = 0
    if end_index==None:
        end_index = N_img

    output_npy = OUTPUT_DIR / 'eigen_vec'
    npy_list = sorted(list(output_npy.glob('*.npy')))[start_index:end_index]

    shape = (end_index - start_index, h, w)
    vector_field = np.empty((3,) + shape)
      
      
    blocks = [npy_list[i:i + bin_factor] for i in range(0, len(npy_list), bin_factor)]

    bin_array = np.empty((3, len(blocks), math.ceil(h/bin_factor), math.ceil(w/bin_factor)))
    # bin_array = np.empty((3, len(blocks),h, w))
    # Load and assign data in chunks based on bin factor
    for i, b in enumerate(blocks):
        print(f"Processing block: {i}/{len(blocks)}")
        array = np.empty((3, len(b), h, w))
        for idx, p in enumerate(b):
            print(f"Reading file: {p.name}")
            
            # Load the numpy data
            array[:,idx,:,:] = np.load(p)  # Shape should match the expected volume slice
            
            
        array = array.mean(axis=1)
        
        # Define the block size for each axis
        block_size = (bin_factor, bin_factor)
        
        # Use block_reduce to bin the volume        
        bin_array[0,i,:,:] = block_reduce(array[0,:,:], block_size=block_size, func=np.mean) 
        bin_array[1,i,:,:] = block_reduce(array[1,:,:], block_size=block_size, func=np.mean) 
        bin_array[2,i,:,:] = block_reduce(array[2,:,:], block_size=block_size, func=np.mean) 

        # bin_array[:,i,:,:] = array[:,0,:,:]

    vector_field = bin_array
    shape = bin_array.shape[1:]


    # Check where the z-component (index 2) is negative
    negative_z_mask = vector_field[0, :, :, :] < 0

    # Flip the vectors where the z-component is negative
    vector_field[:, negative_z_mask] *= -1
    
    
    
    mask_volume = np.where(np.isnan(vector_field[0,:,:,:]), 0, 1)

    # for i in range(0,3):
    #     mask_volume[vector_field[i, :, :, :] == 0] = 0

    
    mask_volume = mask_volume.astype(np.uint8)


    # for i in range(0,3):
    #     vector_field[i,:,:,:][vector_field[0,:,:,:] <= 0] = 0
    
    
    # #---------------------------------------------
    # # Binning
    # print("Binning...")
    
    # # Define the binning factor
    # bin_factor = 32  # Adjust this as needed

    # # Calculate new dimensions by cropping to the nearest multiple of bin_factor
    # z_new = vector_field.shape[1] - (vector_field.shape[1] % bin_factor)
    # y_new = vector_field.shape[2] - (vector_field.shape[2] % bin_factor)
    # x_new = vector_field.shape[3] - (vector_field.shape[3] % bin_factor)

    # # Crop the array if necessary to make dimensions multiples of bin_factor
    # vector_field_cropped = vector_field[:, :z_new, :y_new, :x_new]

    # # Reshape and bin by averaging
    # vector_field_binned = vector_field_cropped.reshape(3, z_new // bin_factor, bin_factor, y_new // bin_factor, bin_factor, x_new // bin_factor, bin_factor).mean(axis=(2, 4, 6))
    
    # vector_field = vector_field_binned
    # shape = vector_field_binned.shape[1:]





    #---------------------------------------------
    # HA

    output_HA = OUTPUT_DIR / 'HA'
    HA_list = sorted(list(output_HA.glob('*.tif'))) + sorted(list(output_HA.glob('*.jp2')))
    HA_volume = load_volume(HA_list[start_index:end_index])
    # mask_volume = np.where(HA_volume == 0, 0, 1)

    # Define the block size for each axis
    block_size = (bin_factor, bin_factor, bin_factor)
    # block_size = (bin_factor, 1, 1)

    # Use block_reduce to bin the volume
    HA_volume = block_reduce(HA_volume, block_size=block_size, func=np.mean)  
    
    # mask_volume = block_reduce(mask_volume, block_size=block_size, func=np.mean)  
    # mask_volume = np.where(mask_volume < 0.5, 0, 1)

    # HA_volume = HA_volume *90/255 - 90


    
    kernel = np.ones((3, 3), np.uint8)  # A 3x3 kernel; adjust as needed for the erosion effect

    # Initialize an empty array to store the eroded volume
    eroded_volume = np.zeros_like(mask_volume)

    # Perform erosion slice by slice along the chosen axis (e.g., the z-axis)
    for i in range(mask_volume.shape[0]):  # Loop through each 2D slice
        eroded_volume[i] = cv2.erode(mask_volume[i], kernel, iterations=1)
      
    mask_volume = eroded_volume

        
        
    
    
    cellData = {}
    # cellData["eigenVectors"] = vector_field.reshape((shape[0], h, w, 3))
    cellData["eigenVectors"] = np.moveaxis(vector_field, 0, -1)
        
    # cellData["eigenVectors"] = vector_field.reshape((shape[0], shape[1], shape[2], 3))
    cellData["HA_angles"] = HA_volume.reshape(shape)
    cellData["mask"] = mask_volume.reshape(shape)

    # Overwrite nans and infs with 0, rubbish I know
    cellData["eigenVectors"][np.logical_not(np.isfinite(cellData["eigenVectors"]))] = 0
    # cellData["eigenVectors"] = np.nan_to_num(cellData["eigenVectors"])

    cellData["HA_angles"][np.logical_not(np.isfinite(cellData["HA_angles"]))] = 0
    cellData["mask"][np.logical_not(np.isfinite(cellData["mask"]))] = 0


    
    # print("eigenVectors shape:", cellData["eigenVectors"].shape)
    # print("HA_angles shape:", cellData["HA_angles"].shape)
    # print("mask shape:", cellData["mask"].shape)
    
    
    try:
        
        for idx in range(0,cellData["eigenVectors"].shape[0]):
            fig, axes = plt.subplots(2, 3, figsize=(15, 5))
            # Plot each slice in a separate subplot
            axes[0,0].imshow(cellData["eigenVectors"][idx, :, :, 0])
            axes[0,1].imshow(cellData["eigenVectors"][idx, :, :, 1])
            axes[0,2].imshow(cellData["eigenVectors"][idx, :, :, 2])
            axes[1,0].imshow(cellData["mask"][idx, :, :])
            axes[1,1].imshow(cellData["HA_angles"][idx])       
            plt.show()
            
    except:
        print("\n/!\ C'ant plot graph\n")

               

    
    
    vtf_name = OUTPUT_DIR / "paraview.vtk"
    print(f"Writing the .vtk file: {vtf_name}")
    writeStructuredVTK(cellData=cellData, fileName=vtf_name)

