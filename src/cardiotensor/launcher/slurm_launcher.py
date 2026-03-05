import glob
import math
import os
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import read_conf_file


def submit_job_to_slurm(
    conf_file_path: str,
    start_index: int,
    end_index_exclusive: int,
    chunk_size: int,
    *,
    partition: str | None,
    time_limit: str,
    cpus_per_task: int,
    mem_gb: int,
    array_parallel: int,
    log_dir: str,
    submit_dir: str,
    dry_run: bool = False,
) -> int | None:
    """
    Submit a SLURM array job for a contiguous index window [start_index, end_index_exclusive).
    """
    total_images = end_index_exclusive - start_index
    if total_images <= 0:
        raise ValueError("end_index_exclusive must be greater than start_index.")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0")

    n_jobs = math.ceil(total_images / chunk_size)
    if n_jobs <= 0:
        raise ValueError("No jobs to submit.")

    log_dir_path = Path(log_dir).resolve()
    submit_dir_path = Path(submit_dir).resolve()
    log_dir_path.mkdir(parents=True, exist_ok=True)
    submit_dir_path.mkdir(parents=True, exist_ok=True)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_name = f"cardiotensor_{ts}_{start_index:06d}_{end_index_exclusive - 1:06d}"
    job_filename = submit_dir_path / f"{job_name}.slurm"
    conf_file_abs = str(Path(conf_file_path).resolve())
    conf_file_quoted = shlex.quote(conf_file_abs)

    print(
        f"Preparing SLURM array: n_jobs={n_jobs}, chunk_size={chunk_size}, "
        f"window=[{start_index}, {end_index_exclusive})",
        flush=True,
    )

    partition_line = (
        f"#SBATCH --partition={partition}\n"
        if partition is not None and partition.strip() != ""
        else ""
    )

    slurm_script_content = f"""#!/bin/bash -l
#SBATCH --output={log_dir_path}/slurm-%x-%A_%a.out
{partition_line}#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem_gb}G
#SBATCH --job-name={job_name}
#SBATCH --time={time_limit}
#SBATCH --array=0-{n_jobs - 1}%{array_parallel}

set -euo pipefail

echo ------------------------------------------------------
echo SLURM_NNODES: ${{SLURM_NNODES:-<unset>}}
echo SLURM_JOB_NODELIST: ${{SLURM_JOB_NODELIST:-<unset>}}
echo SLURM_SUBMIT_DIR: ${{SLURM_SUBMIT_DIR:-<unset>}}
echo SLURM_SUBMIT_HOST: ${{SLURM_SUBMIT_HOST:-<unset>}}
echo SLURM_JOB_ID: ${{SLURM_JOB_ID:-<unset>}}
echo SLURM_JOB_NAME: ${{SLURM_JOB_NAME:-<unset>}}
echo SLURM_JOB_PARTITION: ${{SLURM_JOB_PARTITION:-<unset>}}
echo SLURM_NTASKS: ${{SLURM_NTASKS:-<unset>}}
echo SLURM_CPUS_PER_TASK: ${{SLURM_CPUS_PER_TASK:-<unset>}}
echo SLURM_TASKS_PER_NODE: ${{SLURM_TASKS_PER_NODE:-<unset>}}
echo SLURM_NTASKS_PER_NODE: ${{SLURM_NTASKS_PER_NODE:-<unset>}}
echo SLURM_MEM_PER_CPU: ${{SLURM_MEM_PER_CPU:-<unset>}}
echo SLURM_MEM_PER_NODE: ${{SLURM_MEM_PER_NODE:-<unset>}}
echo ------------------------------------------------------

IMAGES_PER_JOB={chunk_size}
START_INDEX_BASE={start_index}
END_INDEX_EXCLUSIVE={end_index_exclusive}

START_INDEX=$(( SLURM_ARRAY_TASK_ID * IMAGES_PER_JOB + START_INDEX_BASE ))
END_INDEX=$(( START_INDEX + IMAGES_PER_JOB ))

if [ "$START_INDEX" -ge "$END_INDEX_EXCLUSIVE" ]; then
  echo "Task index out of range for requested window. Skipping."
  exit 0
fi

if [ "$END_INDEX" -gt "$END_INDEX_EXCLUSIVE" ]; then
  END_INDEX=$END_INDEX_EXCLUSIVE
fi

echo "Start index (inclusive): $START_INDEX"
echo "End index (exclusive):   $END_INDEX"
echo "Requested memory (GB):   {mem_gb}"

# Fix Qt headless error
export QT_QPA_PLATFORM=offscreen

echo cardio-tensor {conf_file_quoted} --start_index "$START_INDEX" --end_index "$END_INDEX"
cardio-tensor {conf_file_quoted} --start_index "$START_INDEX" --end_index "$END_INDEX"
"""

    job_filename.write_text(slurm_script_content)

    if dry_run:
        print(f"[dry-run] Generated script: {job_filename}", flush=True)
        print(f"[dry-run] sbatch {job_filename}", flush=True)
        return None

    try:
        result = subprocess.run(
            ["sbatch", str(job_filename)],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        print(f"⚠️ Failed to submit SLURM script {job_filename}", flush=True)
        if exc.stdout:
            print(exc.stdout, flush=True)
        if exc.stderr:
            print(exc.stderr, flush=True)
        sys.exit(1)

    stdout = result.stdout.strip()
    match = re.search(r"Submitted batch job\s+(\d+)", stdout)
    if match is None:
        numbers = re.findall(r"\d+", stdout)
        if not numbers:
            raise RuntimeError(f"Could not parse SLURM job id from: {stdout}")
        job_id = int(numbers[-1])
    else:
        job_id = int(match.group(1))

    print(f"Submitted sbatch job {job_id} for [{start_index}, {end_index_exclusive})")
    return job_id


def slurm_launcher(
    conf_file_path: str,
    start_index: int = 0,
    end_index: int | None = None,
    chunk_size: int | None = None,
    partition: str | None = None,
    time_limit: str | None = None,
    cpus_per_task: int | None = None,
    mem_gb: int | None = None,
    array_parallel: int | None = None,
    log_dir: str | None = None,
    submit_dir: str | None = None,
    monitor: bool = True,
    dry_run: bool = False,
) -> None:
    """
    Launch SLURM array jobs for a subset [start_index, end_index) of the volume.

    Parameters
    ----------
    conf_file_path : str
        Path to cardiotensor configuration file.
    start_index : int
        Global start index (inclusive).
    end_index : int | None
        Global end index (exclusive). If None, process until the last slice.
    monitor : bool
        If True, wait and monitor outputs after submission.
    dry_run : bool
        If True, generate scripts but do not submit jobs.
    """
    try:
        params = read_conf_file(conf_file_path)
    except Exception as err:
        print(f"⚠️ Error reading parameter file '{conf_file_path}': {err}", flush=True)
        sys.exit(1)

    volume_path = params.get("IMAGES_PATH", "")
    output_dir = params.get("OUTPUT_PATH", "./output")
    output_format = params.get("OUTPUT_FORMAT", "jp2")
    angle_mode = str(params.get("ANGLE_MODE", "ha_ia")).strip().lower()
    write_angles = bool(params.get("WRITE_ANGLES", True))
    write_vectors = bool(params.get("WRITE_VECTORS", False))
    is_test = bool(params.get("TEST", False))

    if is_test:
        sys.exit(
            "Test mode is enabled. Disable TEST in the configuration for SLURM runs."
        )

    default_chunk_size = int(params.get("N_CHUNK", 100))
    chunk_size = default_chunk_size if chunk_size is None else int(chunk_size)
    partition = None if partition is None else str(partition).strip()
    if partition == "":
        partition = None
    time_limit = "2:00:00" if time_limit is None else str(time_limit).strip()
    cpus_per_task = 8 if cpus_per_task is None else int(cpus_per_task)
    mem_gb = 64 if mem_gb is None else int(mem_gb)
    array_parallel = 100 if array_parallel is None else int(array_parallel)
    default_log_dir = os.path.join(output_dir, "slurm", "log")
    default_submit_dir = os.path.join(output_dir, "slurm", "submit")
    log_dir = default_log_dir if log_dir is None else str(log_dir)
    submit_dir = default_submit_dir if submit_dir is None else str(submit_dir)

    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0")
    if cpus_per_task <= 0:
        raise ValueError("cpus_per_task must be > 0")
    if mem_gb <= 0:
        raise ValueError("mem_gb must be > 0")
    if array_parallel <= 0:
        raise ValueError("array_parallel must be > 0")

    data_reader = DataReader(volume_path)
    total_slices = int(data_reader.shape[0])
    if total_slices <= 0:
        raise ValueError(f"Dataset at {volume_path} contains no slices.")

    first = max(0, int(start_index))
    last_exclusive = total_slices if end_index is None else int(end_index)
    last_exclusive = max(0, min(last_exclusive, total_slices))

    if last_exclusive <= first:
        raise ValueError(
            f"Invalid range: start_index={first}, end_index={last_exclusive} "
            f"(must satisfy start_index < end_index <= {total_slices})."
        )

    window_len = last_exclusive - first
    print(
        f"Processing slice window [{first}, {last_exclusive}) "
        f"(len={window_len}) out of 0..{total_slices}",
        flush=True,
    )
    print(
        "SLURM settings: "
        f"partition={partition}, time_limit={time_limit}, cpus_per_task={cpus_per_task}, "
        f"mem_gb={mem_gb}, array_parallel={array_parallel}, chunk_size={chunk_size}",
        flush=True,
    )
    print(f"SLURM paths: log_dir={log_dir}, submit_dir={submit_dir}", flush=True)

    def build_intervals(
        first_idx: int, last_idx_exclusive: int, step: int
    ) -> list[tuple[int, int]]:
        out: list[tuple[int, int]] = []
        s = first_idx
        while s < last_idx_exclusive:
            e = min(s + step, last_idx_exclusive)
            out.append((s, e))
            s = e
        return out

    intervals = build_intervals(first, last_exclusive, chunk_size)
    n_jobs_total = len(intervals)
    print(
        f"Splitting data into {n_jobs_total} jobs of up to {chunk_size} slices each",
        flush=True,
    )

    max_tasks_per_array = 999
    batched = [
        intervals[i : i + max_tasks_per_array]
        for i in range(0, n_jobs_total, max_tasks_per_array)
    ]
    print(
        f"Launching {len(batched)} array batch(es) "
        f"(tasks per batch: {[len(batch) for batch in batched]})",
        flush=True,
    )

    job_ids: list[int] = []
    start_t = time.time()
    for batch in batched:
        batch_start = batch[0][0]
        batch_end_exclusive = batch[-1][1]
        job_id = submit_job_to_slurm(
            conf_file_path,
            batch_start,
            batch_end_exclusive,
            chunk_size,
            partition=partition,
            time_limit=time_limit,
            cpus_per_task=cpus_per_task,
            mem_gb=mem_gb,
            array_parallel=array_parallel,
            log_dir=log_dir,
            submit_dir=submit_dir,
            dry_run=dry_run,
        )
        if job_id is not None:
            job_ids.append(job_id)

    if dry_run:
        print("Dry run complete. No jobs were submitted.", flush=True)
        return

    if not monitor:
        print("Submission complete. Monitoring disabled by user.", flush=True)
        return

    monitor_job_output(
        output_directory=output_dir,
        start_index=first,
        end_index_exclusive=last_exclusive,
        output_format=output_format,
        angle_mode=angle_mode,
        write_angles=write_angles,
        write_vectors=write_vectors,
    )

    elapsed = time.time() - start_t
    print(
        f"SLURM jobs submitted: {job_ids if job_ids else '[unknown ids]'}\n"
        f"Elapsed submit+monitor time (s): {elapsed:.1f}",
        flush=True,
    )


def _count_files_in_range(
    folder: str,
    prefix: str,
    extension: str,
    start_index: int,
    end_index_exclusive: int,
) -> int:
    """
    Count files named like {prefix}_{index:06d}.{extension} within [start, end).
    """
    if not os.path.isdir(folder):
        return 0

    ext = extension.lstrip(".")
    pattern = os.path.join(folder, f"{prefix}_*.{ext}")
    rgx = re.compile(rf"^{re.escape(prefix)}_(\d+)\.{re.escape(ext)}$")

    count = 0
    for file_path in glob.glob(pattern):
        name = os.path.basename(file_path)
        match = rgx.match(name)
        if match is None:
            continue
        idx = int(match.group(1))
        if start_index <= idx < end_index_exclusive:
            count += 1
    return count


def monitor_job_output(
    output_directory: str,
    start_index: int,
    end_index_exclusive: int,
    output_format: str = "jp2",
    angle_mode: str = "ha_ia",
    write_angles: bool = True,
    write_vectors: bool = False,
    poll_interval_sec: int = 60,
) -> None:
    """
    Monitor output progress for the requested index range.

    The monitor tracks one representative output per slice:
    - HA or AZ image when `write_angles=True`
    - eigen_vec when only vectors are written
    """
    total_images = end_index_exclusive - start_index
    if total_images <= 0:
        return

    if write_angles:
        if angle_mode == "az_el":
            mode_prefix = "AZ"
        else:
            mode_prefix = "HA"
        folder = os.path.join(output_directory, mode_prefix)
        extension = output_format
        prefix = mode_prefix
    elif write_vectors:
        folder = os.path.join(output_directory, "eigen_vec")
        extension = "npy"
        prefix = "eigen_vec"
    else:
        print(
            "No angle/vector outputs requested; skipping monitor (nothing to track).",
            flush=True,
        )
        return

    print(
        f"Monitoring outputs in {folder} for range [{start_index}, {end_index_exclusive})",
        flush=True,
    )

    start_t = time.time()
    prev_count = _count_files_in_range(
        folder, prefix, extension, start_index, end_index_exclusive
    )
    while True:
        current_count = _count_files_in_range(
            folder, prefix, extension, start_index, end_index_exclusive
        )
        processed = current_count
        remaining = max(total_images - processed, 0)
        print(f"{processed}/{total_images} processed", flush=True)

        if processed >= total_images:
            break

        delta = current_count - prev_count
        if delta > 0:
            rate_per_min = delta * (60.0 / max(poll_interval_sec, 1))
            eta_min = remaining / rate_per_min if rate_per_min > 0 else float("inf")
            print(
                f"{delta} new files in last {poll_interval_sec}s. "
                f"Approx. {eta_min:.2f} minutes remaining.",
                flush=True,
            )

        prev_count = current_count
        print(f"Elapsed time (s): {time.time() - start_t:.1f}", flush=True)
        print(f"Waiting {poll_interval_sec} seconds...\n", flush=True)
        time.sleep(poll_interval_sec)


def is_chunk_done(
    output_dir: str,
    start: int,
    end: int,
    output_format: str = "jp2",
    angle_mode: str = "ha_ia",
) -> bool:
    """
    Check if all output files for a given chunk [start, end) are already present.
    """
    ext = output_format.lstrip(".")
    mode = angle_mode.lower().strip()
    if mode == "az_el":
        angle1, angle2 = "AZ", "EL"
    else:
        angle1, angle2 = "HA", "IA"

    for idx in range(start, end):
        expected_files = [
            f"{output_dir}/{angle1}/{angle1}_{idx:06d}.{ext}",
            f"{output_dir}/{angle2}/{angle2}_{idx:06d}.{ext}",
            f"{output_dir}/FA/FA_{idx:06d}.{ext}",
        ]
        if not all(os.path.exists(path) for path in expected_files):
            return False
    return True
