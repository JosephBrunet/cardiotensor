import os
import sys
import math
import numpy as np
import cv2
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from skimage.measure import profile_line
from pathlib import Path
import math
import argparse
import csv
import pandas as pd


from MyoTensor.utils import *


def calculate_angle_line(start, end):
    """
    Calculate the angle of a line defined by start and end points.
    
    Parameters:
    start (tuple): The starting point of the line (x1, y1).
    end (tuple): The ending point of the line (x2, y2).
    
    Returns:
    float: The angle of the line in degrees.
    """
    x1, y1 = start
    x2, y2 = end
    
    delta_x = x2 - x1
    delta_y = y2 - y1
    
    # Calculate the angle in radians
    angle_rad = math.atan2(delta_y, delta_x)
    
    # Convert the angle to degrees
    angle_deg = math.degrees(angle_rad)
    
    return angle_deg



def find_end_points(start_point, end_point, angle_range, N_line):
    
    theta = calculate_angle_line(start_point, end_point)
    
    if N_line > 1:
        theta_list = np.linspace(theta-angle_range/2,theta+angle_range/2,N_line)
    else:
        theta_list = [theta]
        
    vector = np.array(end_point) - np.array(start_point)
    norm = np.linalg.norm(vector)
    
    # Calculate end points 
    end_points = []
    for angle in theta_list:
        theta = np.deg2rad(angle)
        end_x = int(start_point[0] + norm * np.cos(theta))
        end_y = int(start_point[1] + norm * np.sin(theta))
        end_points.append((end_x, end_y))  
        
    return np.array(end_points)


def calculate_intensities(img_helix, start_point, end_point, angle_range=5, N_line=10):
      
    end_points = find_end_points(start_point, end_point, angle_range, N_line)
    
    img_helix[np.isnan(img_helix)] = 0
    
    intensity_profiles = []
    for i, end in enumerate(end_points):
        print(f"Measure {i+1}/{len(end_points)}")
        intensity_profile = profile_line(img_helix, start_point, end, order=0) * 180/255 - 90
        intensity_profiles.append(intensity_profile)
   
    return intensity_profiles


def plot_intensity(intensity_profiles):
    plt.figure(figsize=(10, 6))

    min_length = min(intensity_profile.shape[0] for intensity_profile in intensity_profiles)
    
    # Trim the arrays to the minimum length
    trimmed_arrays = [intensity_profile[:min_length] for intensity_profile in intensity_profiles]

    # Convert list of trimmed arrays to a 2D NumPy array
    intensity_profiles = np.stack(trimmed_arrays)
    
    # for idx,i in enumerate(intensity_profiles):
    #     plt.plot(i, label=f'{idx}', linewidth=0.5, c='k')

    mean_array = np.mean(intensity_profiles, axis=0)
    median_array = np.median(intensity_profiles, axis=0)
    
    # Calculate the 2.5th and 97.5th percentiles
    lower_percentile = np.percentile(intensity_profiles, 5, axis=0)
    upper_percentile = np.percentile(intensity_profiles, 95, axis=0)
    
    # Plot the mean and median
    plt.plot(mean_array, label='Mean')
    plt.plot(median_array, label='Median')
    
    # Add shaded area for the 95% centiles
    plt.fill_between(range(min_length), lower_percentile, upper_percentile, color='gray', alpha=0.5, label='95% Centiles')
    
    # plt.ylim(-90, 90)
    plt.title('Intensity Profiles along the Lines')
    plt.xlabel('Pixel position along the line')
    plt.ylabel('Intensity')
    plt.legend()
    plt.show()

def save_intensity(intensity_profiles, save_path):
    
    # Convert intensity_profiles into a DataFrame
    df = pd.DataFrame(intensity_profiles)
    
    # Rename columns to start from "Value 1"
    df.columns = [f"Value {i+1}" for i in range(df.shape[1])]
    
    # Add a "Profile" row
    df.insert(0, "Profile", [f"Profile {i+1}" for i in range(df.shape[0])])
    
    # Transpose the DataFrame so profiles are in columns
    df = df.transpose()
    
    # Set the first row as the header
    df.columns = df.iloc[0]
    df = df[1:]
    
    # Save DataFrame to CSV with a semicolon delimiter
    df.to_csv(save_path, sep=';', index=False)
    
    print(f"Profile saved to {save_path}")
    
    
    

def main():
    if len(sys.argv) < 2:
        Tk().withdraw()
        conf_file_path = askopenfilename(initialdir=f"{os.getcwd()}/param_files", title="Select file") # show an "Open" dialog box and return the path to the selected file
        if not conf_file_path:
            sys.exit("No file selected!")
        N_slice = ''
        
    elif len(sys.argv) >= 2:
        parser = argparse.ArgumentParser(description='Convert images between tif and jpeg2000 formats.')
        parser.add_argument('conf_file_path', type=str, help='Path to the input text file.')
        parser.add_argument('--slice', type=int, default=0, help='Starting index for processing.')
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
        N_slice = args.slice
    
    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        sys.exit(f'⚠️  Error reading parameter file: {conf_file_path}')
    
    VOLUME_PATH, MASK_PATH, IS_FLIP, OUTPUT_DIR, OUTPUT_TYPE, SIGMA, RHO, N_CHUNK, PT_MV, PT_APEX, IS_TEST, N_SLICE_TEST = [params[key] for key in ['IMAGES_PATH', 'MASK_PATH', 'FLIP', 'OUTPUT_PATH', 'OUTPUT_TYPE', 'SIGMA', 'RHO', 'N_CHUNK', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'TEST', 'N_SLICE_TEST']]
    
    img_list, img_type = get_image_list(VOLUME_PATH)
    
    if not N_slice:
        N_slice = int(input(f"Slice number (0 - {len(img_list)}) ? "))
        
    slice_path = Path(OUTPUT_DIR) / "HA" / f"HA_{(N_slice):06d}.tif"
        
    if not slice_path.exists():
        print(f"File {slice_path} doesn't exist")
        return
    else:
        print(f"Slice found ({slice_path})")
        
    img = cv2.imread(img_list[N_slice], -1)
    img_helix = cv2.imread(str(slice_path), -1)
    
    center_line = interpolate_points(PT_MV, PT_APEX, len(img_list))
    center_point = np.around(center_line[N_slice])
    
    
    
    start = (7500, 7500)  # Starting point of the line (row, column)
    end = (5380, 11340)    # Ending point of the line (row, column)
    angle_range = 5 #degree
    N_line = 5
    
    
    find_end_points(start_point, end_point, angle_range, N_line)
        
    
    
        
    print("Measure the intensity along the lines")
    intensity_profiles = calculate_intensities(img_helix, start_point, end_points)
    plot_intensity(intensity_profiles)
    


    # # Plot the intensity profile
    # plt.figure(figsize=(10, 4))
    # plt.plot(intensity_profile)
    # plt.title('Intensity Profile along the Line')
    # plt.xlabel('Pixel position along the line')
    # plt.ylabel('Intensity')
    # plt.show()
    
    
    
    
    
    print("\nPlotting images...")
    img_vmin,img_vmax = np.nanpercentile(img, (5, 95))
    orig_map=plt.get_cmap('hsv')
    reversed_map = orig_map.reversed()
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    ax = axes
    ax[0].imshow(img, vmin=img_vmin, vmax=img_vmax ,cmap=plt.cm.gray)   
    # Draw a red point at the specified coordinate
    x, y = center_point[0:2]
    ax[0].scatter(x, y, c='red', s=50, marker='o')    
    # ax[0,0].scatter(PT_MV[0],PT_MV[1], c='green', s=50, marker='o')
    # ax[0,0].scatter(PT_APEX[0],PT_APEX[1], c='blue', s=50, marker='o')
    tmp = ax[1].imshow(img_helix,cmap=orig_map)

    ax[0].plot([start[1], end[1]], [start[0], end[0]], marker = 'o', linewidth=3)
    ax[1].plot([start[1], end[1]], [start[0], end[0]], marker = 'o', linewidth=3)
    for e in end_points:
        ax[1].plot([start[1], e[1]], [start[0], e[0]], marker = 'o', linewidth=3)
    fig.colorbar(tmp) 
                             
    plt.show()
    # show slice
    # analysis ?
    
    
    print("FINISH ! ")
    
if __name__ == '__main__':  
    main()