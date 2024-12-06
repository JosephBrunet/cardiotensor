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
import multiprocessing as mp

from cardiotensor.export import amira_writer




def script():
    parser = argparse.ArgumentParser(description='Convert images between tif and jpeg2000 formats')
    parser.add_argument('conf_file_path', type=str, help='Path to the input text file')
    parser.add_argument('--start_index', type=int, default=None, help='Start index for volume subset')
    parser.add_argument('--end_index', type=int, default=None, help='End index for volume subset')
    parser.add_argument('--bin', type=int, default=None, help='binning volume')
    parser.add_argument("--num_points", type=int, default=20000, help="Number of starting random points.")
    parser.add_argument("--num_steps", type=int, default=1000000, help="Number of steps to follow vectors.")
    parser.add_argument("--segment_length", type=int, default=20, help="Length of each segment.")
    parser.add_argument("--angle_threshold", type=float, default=60, help="Maximum allowed angle change.")
    parser.add_argument("--segment_min_length_threshold", type=int, default=30, help="Minimum length of valid fibers.")

    args = parser.parse_args()
    
    conf_file_path = args.conf_file_path
    start_index = args.start_index
    end_index = args.end_index
    bin_factor = args.bin
    num_points = args.num_points
    num_steps = args.num_steps
    segment_length = args.segment_length
    angle_threshold = args.angle_threshold
    segment_min_length_threshold = args.segment_min_length_threshold
    
    amira_writer(conf_file_path,start_index,end_index, bin_factor)
