import glob
import inspect
import math
import os
import subprocess
import sys
import time
from datetime import datetime

from cardiotensor.orientation.orientation_computation_pipeline import (
    compute_orientation,
)
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import (
    read_conf_file,
)


def submit_job_to_slurm(
    executable_path: str,
    conf_file_path: str,
    start_image: int,
    end_image: int,
    N_chunk: int = 10,
    mem_needed: int = 64,
) -> int:
    """
    Submit a Slurm job and return its job ID.

    Args:
        executable_path (str): Path to the executable script.
        conf_file_path (str): Path to the configuration file.
        start_image (int): Index of the first image to process.
        end_image (int): Index of the last image to process.
        N_chunk (int, optional): Number of chunks for the job. Default is 10.
        mem_needed (int, optional): Memory required in GB. Default is 64.

    Returns:
        int: The Slurm job ID.
    """
    log_dir = "/tmp_14_days/bm18/slurm/log/"
    submit_dir = "/tmp_14_days/bm18/slurm/submit"

    executable_path = executable_path.split(".py")[0]
    executable = os.path.basename(executable_path)
    print(f"Script to start: {executable}")

    # Get the current date in YYYY-MM-DD format
    current_date = datetime.now().strftime("%Y-%m-%d")

    job_name = f"{executable}_{current_date}"
    # Generate the job filename
    job_filename = f"{submit_dir}/{job_name}.slurm"

    # Calculate the total number of images to process
    total_images = end_image - start_image + 1
    IMAGES_PER_JOB = N_chunk
    N_jobs = math.ceil(total_images / IMAGES_PER_JOB)

    print(f"\nN_jobs = {N_jobs}, IMAGES_PER_JOB = {IMAGES_PER_JOB}")

    slurm_script_content = f"""#!/bin/bash -l
#SBATCH --output={log_dir}/slurm-%x-%j.out
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
###SBATCH --gres=gpu:2
#SBATCH --mem={math.ceil(mem_needed)}G
#SBATCH --job-name={job_name}
#SBATCH --time=2:00:00
#SBATCH --array=1-{N_jobs}%400

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

START_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) * {IMAGES_PER_JOB} + {start_image}))
END_INDEX=$(( SLURM_ARRAY_TASK_ID * {IMAGES_PER_JOB} + {start_image}))
if [ $END_INDEX -ge {total_images} ]; then END_INDEX={total_images} - 1; fi
echo Start index, End index : $START_INDEX: $END_INDEX

echo mem used {math.ceil(mem_needed)}G

# Starting python script
echo cardio-tensor {conf_file_path} --start_index $START_INDEX --end_index $END_INDEX

# cardio-tensor {executable_path}.py {conf_file_path} --start_index $START_INDEX --end_index $END_INDEX
cardio-tensor {conf_file_path} --start_index $START_INDEX --end_index $END_INDEX

"""

    with open(job_filename, "w") as file:
        file.write(slurm_script_content)

    try:
        result = subprocess.run(
            ["sbatch", job_filename], capture_output=True, text=True, check=True
        )
        job_id = result.stdout.split()[-1]
        print(f"sbatch {job_id} - Index {start_image} to {end_image}")
        return int(job_id)
    except subprocess.CalledProcessError:
        print(f"⚠️ - Failed to submit Slurm job with script {job_filename}")
        sys.exit()


def monitor_job_output(
    output_directory: str, total_images: int, file_extension: str
) -> None:
    """
    Monitor the output directory until all images are processed.

    Args:
        output_directory (str): Path to the directory to monitor.
        total_images (int): Total number of expected images.
        file_extension (str): File extension to monitor for.

    Returns:
        None
    """
    time.sleep(1)
    tmp_count = len(glob.glob(f"{output_directory}/HA/*"))
    while True:
        current_files_count = len(glob.glob(f"{output_directory}/HA/*"))

        print(f"{current_files_count}/{total_images} processed")

        if current_files_count > tmp_count:
            rate = current_files_count - tmp_count  # images per minute
            remaining_time = (total_images - current_files_count) / rate  # minutes
            print(
                f"{current_files_count - tmp_count} images processed in 60sec. Approximately {remaining_time:.2f} minutes remaining"
            )
        tmp_count = current_files_count

        if current_files_count >= total_images:
            break

        print("\nWaiting 60 seconds...\n")
        time.sleep(60)


def slurm_launcher(conf_file_path: str) -> None:
    """
    Launch Slurm jobs for 3D data processing.

    Args:
        conf_file_path (str): Path to the configuration file.

    Returns:
        None
    """
    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(f"⚠️  Error reading parameter file '{conf_file_path}': {e}")
        sys.exit(1)

    # Extracting parameters safely using .get() with defaults where necessary
    VOLUME_PATH = params.get("IMAGES_PATH", "")
    OUTPUT_DIR = params.get("OUTPUT_PATH", "./output")
    N_CHUNK = params.get("N_CHUNK", 100)
    IS_TEST = params.get("TEST", False)

    if IS_TEST == True:
        sys.exit(
            "Test mode activated, run directly 3D_processing.py or deactivate test mode in the parameter file"
        )

    data_reader = DataReader(VOLUME_PATH)

    mem_needed = 128

    def chunk_split(num_images: int, n: int) -> list[tuple[int, int]]:
        """
        Split a range of images into smaller intervals (chunks) for processing.

        Args:
            num_images (int): Total number of images.
            n (int): Number of chunks to divide the images into.

        Returns:
            List[Tuple[int, int]]: List of tuples where each tuple represents
                                    the start and end indices of a chunk.
        """
        # Calculate the number of images per interval
        images_per_interval = num_images // n

        print(f"Number of images per interval: {images_per_interval}")

        # Create a list of intervals
        intervals = []
        start_index = 0
        while start_index < num_images:
            end_index = start_index + images_per_interval
            intervals.append((start_index, end_index))  # Use tuple here
            start_index = end_index
        intervals[-1] = (
            intervals[-1][0],
            num_images,
        )  # Ensure the last chunk ends at num_images

        return intervals

    n_jobs = math.ceil(data_reader.shape[0] / N_CHUNK)

    index_intervals = chunk_split(data_reader.shape[0], n_jobs)

    print(
        f"Splitting data into {n_jobs} chunks of {N_CHUNK} slices each for processing"
    )

    # Split the index_intervals into batches of max length 1000
    N_job_max_per_array = 999
    batched_intervals = [
        index_intervals[i : i + N_job_max_per_array]
        for i in range(0, len(index_intervals), N_job_max_per_array)
    ]

    print(
        f"Launching {len(batched_intervals)} batches of jobs (Number of jobs: {[len(b) for b in batched_intervals]})"
    )

    answer = input("Do you want to continue? [y]\n")
    if "n" in answer.lower():
        sys.exit("Aborted by user")

    python_file_path = os.path.abspath(inspect.getfile(compute_orientation))

    # Launch each batch
    for batch in batched_intervals:
        start, end = batch[0][0], batch[-1][1]
        job_id = submit_job_to_slurm(
            python_file_path,
            conf_file_path,
            start,
            end,
            N_chunk=N_CHUNK,
            mem_needed=mem_needed,
        )
        print(
            f"Submitted job for batch starting at {start} and ending at {end} (job ID: {job_id})"
        )

    monitor_job_output(OUTPUT_DIR, data_reader.shape[0], conf_file_path)

    return
