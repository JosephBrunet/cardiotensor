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
from typing import Dict, Tuple, Any
from os import PathLike
from alive_progress import alive_bar

import SimpleITK as sitk






def read_conf_file(file_path):
    
    if not os.path.exists(file_path):
        sys.exit(f"The configuration file {file_path} does not exist.")
    
    config = configparser.ConfigParser()
    config.read(file_path)

    config_dict = {
        'IMAGES_PATH': config.get('DATASET', 'IMAGES_PATH', fallback=None).strip(),
        'MASK_PATH': config.get('DATASET', 'MASK_PATH', fallback=None).strip(),
        'FLIP': config.getboolean('DATASET', 'FLIP', fallback=True),
        'OUTPUT_PATH': config.get('OUTPUT', 'OUTPUT_PATH', fallback=None).strip(),
        'OUTPUT_FORMAT': config.get('OUTPUT', 'OUTPUT_FORMAT', fallback='jp2').strip(),
        'OUTPUT_TYPE': config.get('OUTPUT', 'OUTPUT_TYPE', fallback=None).strip(),
        'VECTORS': config.getboolean('OUTPUT', 'VECTORS', fallback=False),
        'SIGMA': config.getfloat('STRUCTURE TENSOR CALCULATION', 'SIGMA', fallback=None),
        'RHO': config.getfloat('STRUCTURE TENSOR CALCULATION', 'RHO', fallback=None),
        'N_CHUNK': config.getint('STRUCTURE TENSOR CALCULATION', 'N_CHUNK', fallback=100),
        'POINT_MITRAL_VALVE': np.array(list(map(float, config.get('LV AXIS COORDINATES', 'POINT_MITRAL_VALVE', fallback=None).split(',')))),
        'POINT_APEX': np.array(list(map(float, config.get('LV AXIS COORDINATES', 'POINT_APEX', fallback=None).split(',')))),
        'REVERSE': config.getboolean('RUN', 'REVERSE', fallback=False),
        'TEST': config.getboolean('TEST', 'TEST', fallback=True),
        'N_SLICE_TEST': config.getint('TEST', 'N_SLICE_TEST', fallback=None)
    }

    return config_dict







def convert_to_8bit(img, perc_min = 0, perc_max = 100, output_min = None, output_max = None):
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
    
    minimum,maximum = np.nanpercentile(img, (perc_min, perc_max))

    # print(f"Minimum, Maximum : {minimum}, {maximum}")
    
    # minimum = int(minimum)
    # maximum = int(maximum)

    if output_min and output_max:   
        minimum = output_min
        maximum = output_max
        
    
    img_normalized = (img + np.abs(minimum)) * (255 / (maximum - minimum))  # Normalize the image to the range [0, 255]
    img_8bit = img_normalized.astype(np.uint8)
    
    return img_8bit





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




def get_volume_shape(VOLUME_PATH) -> Tuple[int, int, int]:
    """
    Determines the volume shape (x, y, z) based on the VOLUME_PATH.

    Parameters:
    - VOLUME_PATH: str
        The path to either a .mhd file or a folder containing slices.
    
    Returns:
    - tuple: The volume shape as (x, y, z).
    """
    # Case 1: If it's a .mhd file
    if os.path.isfile(VOLUME_PATH) and VOLUME_PATH.endswith('.mhd'):
        # Use SimpleITK to read the .mhd file and get the image size
        image = sitk.ReadImage(VOLUME_PATH)
        size = image.GetSize()  # returns (x, y, z)
        return size  # Return the full shape as (x, y, z)

    # Case 2: If it's a folder with image slices
    elif os.path.isdir(VOLUME_PATH):
        # Count the number of image slices in the folder (assuming all are in the same format)
        supported_extensions = ('.tif', '.tiff', 'jp2')
        slice_files = [f for f in os.listdir(VOLUME_PATH) if f.lower().endswith(supported_extensions)]
        z_dim = len(slice_files)  # Z-dimension is the number of slices
        
        if z_dim > 0:
            # Load the first image to get x and y dimensions
            first_image = cv2.imread(os.path.join(VOLUME_PATH, slice_files[0]), cv2.IMREAD_UNCHANGED)
            y_dim, x_dim = first_image.shape[:2]
            return (x_dim, y_dim, z_dim)
        else:
            raise ValueError("No image slices found in the specified directory.")
    
    else:
        raise ValueError("The VOLUME_PATH must either be a .mhd file or a folder containing slices.")



def load_volume(file_list, start_index=0, end_index=0):
    """
    Loads the volume data from a list of file paths with a progress bar.

    Parameters:
    - file_list: list
        List of file paths to load.
    - start_index: int
        The starting index with padding.
    - end_index: int
        The ending index with padding. If 0, loads the entire volume.

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
        if os.path.basename(file_path)[-4:] == '.npy':
            image_data = np.load(str(file_path))
        elif os.path.basename(file_path)[-4:] == '.jp2':
            image_data = cv2.imread(str(file_path), cv2.IMREAD_UNCHANGED)
        else:
            image_data = cv2.imread(str(file_path), cv2.IMREAD_UNCHANGED)
        return image_data

    # If end_index is 0, set it to the length of the file list to load all files
    if end_index == 0:
        end_index = len(file_list)

    # Calculate the total number of files to be loaded
    total_files = end_index - start_index
    print(f"Loading {total_files} files...")

    # Create a progress bar with alive-progress
    with alive_bar(total_files, title="Loading Volume", length=40) as bar:
        # Create a list of delayed tasks to read each image and update the progress bar
        def wrapped_reader(file_path):
            result = custom_image_reader(file_path)
            bar()  # Update the progress bar
            return result

        delayed_tasks = [dask.delayed(lambda x: np.array(wrapped_reader(x)))(file_path) 
                         for file_path in sorted(file_list)]

        # Compute the volume
        volume = np.array(dask.compute(*delayed_tasks[start_index:end_index]))
    
    return volume





def read_mhd(filename: PathLike) -> Dict[str, Any]:
    """
    Return a dictionary of meta data from MHD meta header file
    
    :param filename: file of type .mhd that should be loaded
    :returns: dictionary of meta data
    """
    
    # Define tags
    meta_dict = {}
    tag_set = []
    tag_set.extend(['ObjectType', 'NDims', 'DimSize', 'ElementType', 'ElementDataFile', 'ElementNumberOfChannels'])
    tag_set.extend(['BinaryData', 'BinaryDataByteOrderMSB', 'CompressedData', 'CompressedDataSize'])
    tag_set.extend(['Offset', 'CenterOfRotation', 'AnatomicalOrientation', 'ElementSpacing', 'TransformMatrix'])
    tag_set.extend(['Comment', 'SeriesDescription', 'AcquisitionDate', 'AcquisitionTime', 'StudyDate', 'StudyTime'])

    tag_flag = [False] * len(tag_set)

    with open(filename, "r") as fn:
        line = fn.readline()
        while line:
            tags = str.split(line, '=')
            # print(tags[0])
            for i in range(len(tag_set)):
                tag = tag_set[i]
                if (str.strip(tags[0]) == tag) and (not tag_flag[i]):
                    # print(tags[1])
                    content = str.strip(tags[1])
                    if tag in ['ElementSpacing', 'Offset', 'CenterOfRotation', 'TransformMatrix']:
                        meta_dict[tag] = [float(s) for s in content.split()]
                    elif tag in ['NDims', 'ElementNumberOfChannels']:
                        meta_dict[tag] = int(content)
                    elif tag in ['DimSize']:
                        meta_dict[tag] = [int(s) for s in content.split()]
                    elif tag in ['BinaryData', 'BinaryDataByteOrderMSB', 'CompressedData']:
                        if content == "True":
                            meta_dict[tag] = True
                        else:
                            meta_dict[tag] = False
                    else:
                        meta_dict[tag] = content
                    tag_flag[i] = True
            line = fn.readline()
    return meta_dict




def load_raw_data_with_mhd(filename: PathLike) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Load a MHD file

    :param filename: file of type .mhd that should be loaded
    :returns: tuple with raw data and dictionary of meta data
    """
    meta_dict = read_mhd(filename)
    dim = int(meta_dict['NDims'])
    if "ElementNumberOfChannels" in meta_dict:
        element_channels = int(meta_dict["ElementNumberOfChannels"])
    else:
        element_channels = 1

    if meta_dict['ElementType'] == 'MET_FLOAT':
        np_type = np.float32
    elif meta_dict['ElementType'] == 'MET_DOUBLE':
        np_type = np.float64
    elif meta_dict['ElementType'] == 'MET_CHAR':
        np_type = np.byte
    elif meta_dict['ElementType'] == 'MET_UCHAR':
        np_type = np.ubyte
    elif meta_dict['ElementType'] == 'MET_SHORT':
        np_type = np.int16
    elif meta_dict['ElementType'] == 'MET_USHORT':
        np_type = np.ushort
    elif meta_dict['ElementType'] == 'MET_INT':
        np_type = np.int32
    elif meta_dict['ElementType'] == 'MET_UINT':
        np_type = np.uint32
    else:
        raise NotImplementedError("ElementType " + meta_dict['ElementType'] + " not understood.")
    arr = list(meta_dict['DimSize'])



    volume = np.prod(arr[0:dim - 1])

    pwd = Path(filename).parents[0].resolve()
    data_file = Path(meta_dict['ElementDataFile'])
    if not data_file.is_absolute():
        data_file = pwd / data_file

    shape = (arr[dim - 1], volume, element_channels)
    with open(data_file,'rb') as f:
        data = np.fromfile(f, count=np.prod(shape), dtype=np_type)
    data.shape = shape

    # Adjust byte order in numpy array to match default system byte order
    if 'BinaryDataByteOrderMSB' in meta_dict:
        sys_byteorder_msb = sys.byteorder == 'big'
        file_byteorder_ms =  meta_dict['BinaryDataByteOrderMSB']
        if sys_byteorder_msb != file_byteorder_ms:
            data = data.byteswap()

    # Begin 3D fix
    arr.reverse()
    if element_channels > 1:
        data = data.reshape(arr + [element_channels])
    else:
        data = data.reshape(arr)
    # End 3D fix
    
    return (data, meta_dict)



