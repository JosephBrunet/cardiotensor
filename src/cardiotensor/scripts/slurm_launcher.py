"""
3D_Data_Processing
"""
import os
import sys
import math
import argparse
import numpy as np
import dask_image.imread
import glob
import time
from distutils.util import strtobool
import subprocess
import random
import math
import glob
import inspect

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename, askdirectory



from cardiotensor.launcher import slurm_launcher
from cardiotensor.utils import read_conf_file

            





def script():

    if len(sys.argv) < 2:
        Tk().withdraw()
        conf_file_path = askopenfilename(initialdir=f"{os.getcwd()}/param_files", title="Select file") # show an "Open" dialog box and return the path to the selected file
        if not conf_file_path:
            sys.exit("No file selected!")
        
    elif len(sys.argv) >= 2:
        parser = argparse.ArgumentParser(description='Convert images between tif and jpeg2000 formats.')
        parser.add_argument('conf_file_path', type=str, help='Path to the input text file.')
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
    
    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(e)
        sys.exit(f'⚠️  Error reading parameter file: {conf_file_path}')
    
    VOLUME_PATH, MASK_PATH, IS_FLIP, OUTPUT_DIR, OUTPUT_TYPE, SIGMA, RHO, N_CHUNK, PT_MV, PT_APEX, REVERSE, IS_TEST, N_SLICE_TEST = [params[key] for key in ['IMAGES_PATH', 'MASK_PATH', 'FLIP', 'OUTPUT_PATH', 'OUTPUT_TYPE', 'SIGMA', 'RHO', 'N_CHUNK', 'POINT_MITRAL_VALVE', 'POINT_APEX', 'REVERSE', 'TEST', 'N_SLICE_TEST']]
    

    start_time = time.time()
    slurm_launcher(conf_file_path) 
    print("--- %s seconds ---" % (time.time() - start_time))

    print("FINISH ! ")

    
    
if __name__ == '__main__':  
    main()
    
    
    