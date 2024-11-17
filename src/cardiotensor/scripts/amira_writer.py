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
    args = parser.parse_args()
    
    conf_file_path = args.conf_file_path
    start_index = args.start_index
    end_index = args.end_index
    bin_factor = args.bin
    
    amira_writer(conf_file_path,start_index,end_index, bin_factor)
