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

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename, askdirectory





def read_parameter_file(para_file_path):
    with open(para_file_path, 'r') as file:
        for line in file:
            line = line.replace(' ','').split('#')[0].strip()
            
            if line.startswith('IMAGES_PATH'):
                images_path = line.split('=')[1]
                print("IMAGES_PATH found")
            if line.startswith('OUTPUT_PATH'):
                output_dir = line.split('=')[1]
                print("OUTPUT_PATH found")
            if line.startswith('SIGMA'):
                sigma = float(line.split('=')[1])
                print("SIGMA found")
            if line.startswith('RHO'):
                rho = float(line.split('=')[1])
                print("RHO found")
            if line.startswith('POINT_MITRAL_VALVE'):
                pt_MV = line.split('=')[1]
                pt_MV = np.fromstring(pt_MV, dtype=int, sep=',')
                print("POINT_MITRAL_VALVE found")
            if line.startswith('POINT_APEX'):
                pt_apex = line.split('=')[1]
                pt_apex = np.fromstring(pt_apex, dtype=int, sep=',')
                print("POINT_APEX found")
            if line.startswith('TEST'):
                is_test = strtobool(line.split('=')[1])
                if is_test == 0:
                    is_test = False
                elif is_test == 1:
                    is_test = True
                print("IS_TEST found")
                
    return images_path, output_dir, sigma, rho, pt_MV, pt_apex, is_test






def submit_job_to_slurm(executable_path: str, txt_file_path: str, start_image: int, end_image: int, mem_needed=20, additional_args: str = '') -> int:
    """
    Submit a Slurm job and return its job ID.

    Parameters:
        executable (str): Path to the executable script.
        folder_path (str): Directory path containing the images.
        start_image (int): Index of the first image.
        end_image (int): Index of the last image.
        additional_args (str, optional): Additional arguments for the executable script.

    Returns:
        int: The Slurm job ID.
    """
    # log_dir = '/tmp_14_days/bm18/slurm/log'
    # submit_dir = '/tmp_14_days/bm18/slurm/submit'
    log_dir = '/data/projects/md1290/tmp/slurm/log'
    submit_dir = '/data/projects/md1290/tmp/slurm/submit'

    
    executable_path = executable_path.split('.py')[0]
    
    executable = os.path.basename(executable_path)
    print(f'Script to start: {executable}')

    random_suffix = random.randint(10000, 99999)
    job_name = f"{executable}"
    job_filename = f"{submit_dir}/{job_name}.slurm"

    
    slurm_script_content = f"""#!/bin/bash -l
#SBATCH --output={log_dir}/slurm-%x-%j.out
#SBATCH --partition=nice
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
###SBATCH --gres=gpu:2
#SBATCH --mem={math.ceil(mem_needed)}G
#SBATCH --job-name={job_name}
#SBATCH --time=4:00:00
###SBATCH --exclude=gpbm18-01

echo ------------------------------------------------------                                                     
echo SLURM_NNODES: $SLURM_NNODES                                                     
echo SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST                                                     
echo SLURM_SUBMIT_DIR: $SLURM_SUBMIT_DIR                                                     
echo SLURM_SUBMIT_HOST: $SLURM_SUBMIT_HOST                                                     
echo SLURM_JOB_ID: $SLURM_JOB_ID                                                     
echo SLURM_JOB_NAME: $SLURM_JOB_NAME                                                     
echo SLURM_JOB_PARTITION: $SLURM_JOB_PARTITION                                                     
echo SLURM_NTASKS: $SLURM_NTASKS 
echo SLURM_CPUS-PER-TASK: $SLURM_CPUS_PER_TASK                                                    
echo SLURM_TASKS_PER_NODE: $SLURM_TASKS_PER_NODE                                                     
echo SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE
echo SLURM_MEM_PER_CPU: $SLURM_MEM_PER_CPU
echo SLURM_MEM_PER_NODE: $SLURM_MEM_PER_NODE                                              
echo ------------------------------------------------------ 

module load cuda/12

###echo Activate python environment
###source /home/esrf/joseph08091994/python/pyEnv/bin/activate 

echo mem used {math.ceil(mem_needed)}G

# Starting python script 
echo python3 {executable_path}.py {txt_file_path} --start_index {start_image} --end_index {end_image}

python3 {executable_path}.py {txt_file_path} --start_index {start_image} --end_index {end_image}
"""
    print(f"python3 {executable_path}.py {txt_file_path} --start_index {start_image} --end_index {end_image}")
    print(job_filename)
    with open(job_filename, 'w') as file:
        file.write(slurm_script_content)

    try:
        result = subprocess.run(['sbatch', job_filename], capture_output=True, text=True, check=True)
        job_id = result.stdout.split()[-1]
        print(f"sbatch {job_id} - Index {start_image} to {end_image}")
        return job_id
    except subprocess.CalledProcessError:
        print(f"Failed to submit Slurm job with script {job_filename}")
        return None







def monitor_job_output(output_directory: str, total_images: int, file_extension: str) -> None:
    """
    Monitor the output directory until all images are processed.

    Parameters:
        output_directory (str): Directory path to monitor for output images.
        total_images (int): Total number of images expected.
        file_extension (str): Extension of output files to look for.
    """
    time.sleep(1)
    tmp_count = len(glob.glob(f'{output_directory}/HA/*'))
    while True:
        current_files_count = len(glob.glob(f'{output_directory}/HA/*'))

        print(f"{current_files_count}/{total_images} processed")
        
        if current_files_count > tmp_count:
            rate = (current_files_count - tmp_count) * 60 / 10 # images per minute
            remaining_time = (total_images - current_files_count) / rate # minutes
            print(f"{current_files_count - tmp_count} images processed in 10sec. Approximately {remaining_time:.2f} minutes remaining")
        tmp_count = current_files_count

        if current_files_count >= total_images:
            break
                
        print("\nWaiting 10 seconds...\n")
        time.sleep(10)



    


def run_heart_orientation_processing(file_path):

    images_path, output_dir, sigma, rho, pt_MV, pt_apex, is_test = read_parameter_file(file_path)
    print('PARAMETERS : ',images_path, output_dir, sigma, rho, pt_MV, pt_apex, is_test)
    
    if is_test == True:
        sys.exit('Test mode activated, run directly 3D_processing.py or deactivate test mode in the parameter file')
    
    
    #Check number of files to convert 
    file_list_tif = sorted(glob.glob(images_path + "/*.tif*"))
    file_list_jp2 = sorted(glob.glob(images_path + "/*.jp2*"))

    if file_list_tif and file_list_jp2:
        sys.exit('Both tif and jp2 files were found (check your folder path)')
    elif file_list_tif:
        file_list = file_list_tif      
        file_type = 'tif'
    elif file_list_jp2:
        file_list = file_list_jp2
        file_type = 'jp2'
    else:  
        sys.exit('No files were found (check your folder path)')
    
    print(f"File type: {file_type}")        
    N_img = len(file_list)
    print(f"{N_img} {file_type} files found\n")  
       
    
    if is_test == True:
        N_img = 1

    # Reads files
    print('Reading files...')    
    volume_dask = dask_image.imread.imread(f'{images_path}/*.{file_type}')
    print('Dask volume: ', volume_dask)
    print('Dask size: ',volume_dask.size*2)


    def chunk_split(num_images,n):

        # Calculate the number of images per interval
        images_per_interval = num_images // n
        
        print(f"Number of images per interval: {images_per_interval}")

        # Create a list of intervals
        intervals = []
        start_index = 0
        while start_index < num_images:
            end_index = start_index + images_per_interval
            intervals.append([start_index, end_index])
            start_index = end_index
        intervals[-1][1] = num_images

        print(intervals)

        return intervals

    padding_start = math.ceil(rho)
    padding_end = math.ceil(rho)
    print(f"Padding: {padding_start} and {padding_end}")

    MAX_SIZE = 256 #GB
    MAX_SIZE = MAX_SIZE * 1000000000 
    
    print(f"The size of the volume is : {volume_dask.size/(1024*1024*1024)} GB")
    
    safety_coef = 4
          
    # check to see if data needs splitting into chunks for processing
    if volume_dask.size*4*12.5 * safety_coef >= MAX_SIZE:   # x4 because float32 / x 12.5 for structure tensore calculation 
        # split data into smaller chunks for processing + add padding to avoid boundary errors
        
        mem_image = volume_dask[0,:,:].size*4 *12.5 * safety_coef # x4 because float32 / x 12.5 for structure tensore calculation
        print(f"Memory per image: {mem_image/(1024*1024*1024)} GB")
        
        pad_size = mem_image*(padding_end+padding_start)
        print(f"Padding size: {pad_size/(1024*1024*1024)} GB")
        
        max_N_images = math.floor((MAX_SIZE-pad_size)/mem_image)
        print(f"Max number of images per chunk: {max_N_images}")
        
        split_amount = math.ceil(volume_dask.shape[0]/max_N_images)
        print(f"Splitting data into {split_amount} chunks for processing")
        
        # split_amount = math.ceil(volume_dask.size*4*12.5*1.2 / MAX_SIZE)  # number of chunks to split the data into
        print('\033[93m'+ f'Warning: Large dataset loaded. Data will be split into {split_amount} chunks for processing.' + '\033[0m')
        
        index_intervals = chunk_split(N_img,split_amount)
        
        # if max_N_images < padding_start+padding_end:
        #     sys.exit('Warning: Padding is too large for the chunk size. Not enough memory per chunk to process the data.')
        
        if max_N_images < 5:
            sys.exit('Warning: Chunk size is too small ({max_N_images} images)})')
        
        mem_needed = volume_dask[index_intervals[1][0]-padding_start:index_intervals[1][1]+padding_end,:,:].size*4 *12.5 *safety_coef/1000000000 # x4 because float32 / x 12.5 for structure tensore calculation 
        print(mem_needed)
        mem_needed = MAX_SIZE / 1000000000

        print(mem_needed)
    else:
        print('No chuncked')
        split_amount = 1
        index_intervals = [[0, N_img]]
        
        
    ans = input('Do you want to continue? [y]')
    if 'n' in ans.lower():
        sys.exit('Aborted by user')
    
    for i in index_intervals:
        
        print('Sending to SLURM...')
        job_id = submit_job_to_slurm('processing_heart', file_path, i[0], i[1], mem_needed=mem_needed)
        # sys.exit()
                
    monitor_job_output(output_dir, N_img, file_type)

    print("\n🤖 Beep beep boop! Binning complete, human! 🤖\n")
        

    
        
    return
            










def main():

    if len(sys.argv) < 2:
        Tk().withdraw()
        para_file_path = askopenfilename(initialdir=os.getcwd(), title="Select file") # show an "Open" dialog box and return the path to the selected file
        if not para_file_path:
            sys.exit("No file selected!")

    elif len(sys.argv) >= 2:
        parser = argparse.ArgumentParser(description='Convert images between tif and jpeg2000 formats.')
        parser.add_argument('para_file_path', type=str, help='Path to the input text file.')
        args = parser.parse_args()
        para_file_path = args.para_file_path
        print(para_file_path)

    start_time = time.time()
    run_heart_orientation_processing(para_file_path) 
    print("--- %s seconds ---" % (time.time() - start_time))

    print("FINISH ! ")

    
    
if __name__ == '__main__':  
    main()