import glob
import math
import os
import re
import shlex
import subprocess
import time
from datetime import datetime
from pathlib import Path

from cardiotensor.orientation.orientation_computation_pipeline import (
    _safe_low_memory_chunk_size,
)
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.image_io import (
    initialize_zarr_vector_field,
    normalize_vector_format,
    open_zarr_vector_field,
    vector_field_path,
)
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
echo TMPDIR: ${{TMPDIR:-/tmp}}
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
        details = (exc.stderr or exc.stdout or str(exc)).strip()
        raise RuntimeError(
            f"Failed to submit SLURM script {job_filename}: {details}"
        ) from exc

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
    params = read_conf_file(conf_file_path)

    volume_path = params.get("IMAGES_PATH", "")
    output_dir = params.get("OUTPUT_PATH", "./output")
    output_format = params.get("OUTPUT_FORMAT", "jp2")
    angle_mode = str(params.get("ANGLE_MODE", "ha_ia")).strip().lower()
    write_angles = bool(params.get("WRITE_ANGLES", True))
    write_vectors = bool(params.get("WRITE_VECTORS", False))
    vector_format = normalize_vector_format(params.get("VECTOR_FORMAT", "zarr"))
    projected = bool(params.get("PROJECTED_ANGLES", False))
    low_memory = bool(params.get("LOW_MEMORY", False))
    is_test = bool(params.get("TEST", False))

    if is_test:
        raise ValueError(
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
    if write_vectors and vector_format == "zarr" and not dry_run:
        initialize_zarr_vector_field(output_dir, tuple(data_reader.shape))

    first = max(0, int(start_index))
    last_exclusive = total_slices if end_index is None else int(end_index)
    last_exclusive = max(0, min(last_exclusive, total_slices))

    if last_exclusive <= first:
        raise ValueError(
            f"Invalid range: start_index={first}, end_index={last_exclusive} "
            f"(must satisfy start_index < end_index <= {total_slices})."
        )

    window_len = last_exclusive - first
    if low_memory:
        sigma = float(params.get("SIGMA", 1.0))
        rho = float(params.get("RHO", 3.0))
        truncate = float(params.get("TRUNCATE", 4.0))
        vertical_padding = params.get("VERTICAL_PADDING", None)
        padding = math.ceil(
            vertical_padding
            if vertical_padding is not None
            else int(sigma * truncate + 0.5) + int(rho * truncate + 0.5)
        )
        requested_chunk = min(chunk_size, window_len)
        safe_chunk, requested_peak, safe_peak = _safe_low_memory_chunk_size(
            requested_chunk,
            total_slices,
            data_reader.shape[-2],
            data_reader.shape[-1],
            data_reader.dtype,
            padding=padding,
            include_eigenvalues=write_angles,
            has_mask=params.get("MASK_PATH") is not None,
            write_angles=write_angles,
            sigma=sigma,
            rho=rho,
            truncate=truncate,
            available_memory_bytes=mem_gb * 1024**3,
            cpu_count=cpus_per_task,
        )
        if safe_chunk < requested_chunk:
            print(
                "LOW_MEMORY adjusted SLURM chunk_size from "
                f"{requested_chunk} to {safe_chunk} for {mem_gb} GiB: "
                f"conservative peak {requested_peak / 1024**3:.2f} GiB -> "
                f"{safe_peak / 1024**3:.2f} GiB.",
                flush=True,
            )
            chunk_size = safe_chunk
        else:
            print(
                f"LOW_MEMORY validated chunk_size={chunk_size}: conservative peak "
                f"{safe_peak / 1024**3:.2f} GiB for {mem_gb} GiB requested.",
                flush=True,
            )
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
    all_intervals = intervals
    intervals = [
        interval
        for interval in all_intervals
        if not is_chunk_done(
            output_dir,
            *interval,
            output_format=output_format,
            angle_mode=angle_mode,
            write_angles=write_angles,
            write_vectors=write_vectors,
            vector_format=vector_format,
            projected=projected,
        )
    ]
    skipped_jobs = len(all_intervals) - len(intervals)
    print(
        f"Output check: {skipped_jobs} complete job(s) skipped, "
        f"{len(intervals)} job(s) needed",
        flush=True,
    )
    if not intervals:
        print("All requested outputs already exist. Nothing to submit.", flush=True)
        return

    n_jobs_total = len(intervals)
    print(
        f"Splitting data into {n_jobs_total} jobs of up to {chunk_size} slices each",
        flush=True,
    )

    max_tasks_per_array = 999
    batched: list[list[tuple[int, int]]] = []
    for interval in intervals:
        if (
            not batched
            or len(batched[-1]) >= max_tasks_per_array
            or batched[-1][-1][1] != interval[0]
        ):
            batched.append([interval])
        else:
            batched[-1].append(interval)
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
        vector_format=vector_format,
        projected=projected,
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
    vector_format: str = "zarr",
    poll_interval_sec: int = 60,
    projected: bool = False,
) -> None:
    """
    Monitor every requested output for the requested index range.

    Completion requires both angle outputs, FA, and vectors when enabled.
    """
    total_images = end_index_exclusive - start_index
    if total_images <= 0:
        return

    vector_format = normalize_vector_format(vector_format)
    mode = angle_mode.strip().lower()
    if mode == "az_el":
        angle_names = ("AZ", "EL")
    elif mode == "ha_ia":
        angle_names = ("HA_projected", "IA_projected") if projected else ("HA", "IA")
    else:
        raise ValueError("ANGLE_MODE must be 'ha_ia' or 'az_el'")

    if not write_angles and not write_vectors:
        print(
            "No angle/vector outputs requested; skipping monitor (nothing to track).",
            flush=True,
        )
        return

    print(
        f"Monitoring all requested outputs for range "
        f"[{start_index}, {end_index_exclusive})",
        flush=True,
    )

    def count_completed() -> dict[str, int]:
        counts = {}
        if write_angles:
            for name in (*angle_names, "FA"):
                counts[name] = _count_files_in_range(
                    os.path.join(output_directory, name),
                    name,
                    output_format,
                    start_index,
                    end_index_exclusive,
                )

        if write_vectors:
            if vector_format == "zarr":
                store = open_zarr_vector_field(output_directory)
                counts["Zarr vectors"] = int(
                    store.completed_range(start_index, end_index_exclusive).sum()
                )
            else:
                counts["NPY vectors"] = _count_files_in_range(
                    str(vector_field_path(output_directory, "npy")),
                    "eigen_vec",
                    "npy",
                    start_index,
                    end_index_exclusive,
                )
        return counts

    start_t = time.time()
    previous_total = None
    while True:
        counts = count_completed()
        print(
            " | ".join(
                f"{name}: {min(count, total_images)}/{total_images}"
                for name, count in counts.items()
            ),
            flush=True,
        )

        if all(count >= total_images for count in counts.values()):
            break

        current_total = sum(counts.values())
        target_total = total_images * len(counts)
        remaining = max(target_total - current_total, 0)
        delta = 0 if previous_total is None else max(current_total - previous_total, 0)
        if delta:
            rate_per_min = delta * (60.0 / max(poll_interval_sec, 1))
            eta_min = remaining / rate_per_min if rate_per_min > 0 else float("inf")
            print(
                f"{delta} new outputs in last {poll_interval_sec}s. "
                f"Approx. {eta_min:.2f} minutes remaining.",
                flush=True,
            )

        previous_total = current_total
        print(f"Elapsed time (s): {time.time() - start_t:.1f}", flush=True)
        print(f"Waiting {poll_interval_sec} seconds...\n", flush=True)
        time.sleep(poll_interval_sec)


def is_chunk_done(
    output_dir: str,
    start: int,
    end: int,
    output_format: str = "jp2",
    angle_mode: str = "ha_ia",
    write_angles: bool = True,
    write_vectors: bool = False,
    vector_format: str = "zarr",
    projected: bool = False,
) -> bool:
    """Return True when every requested output exists for [start, end)."""
    ext = output_format.lstrip(".")
    mode = angle_mode.lower().strip()
    if mode == "az_el":
        angle1, angle2 = "AZ", "EL"
    else:
        angle1, angle2 = (
            ("HA_projected", "IA_projected") if projected else ("HA", "IA")
        )

    vector_format = normalize_vector_format(vector_format)
    zarr_completed = None
    if write_vectors and vector_format == "zarr":
        try:
            zarr_completed = open_zarr_vector_field(output_dir).completed_range(
                start, end
            )
        except (FileNotFoundError, KeyError, ValueError):
            return False

    for idx in range(start, end):
        expected_files = []
        if write_angles:
            expected_files.extend(
                [
                    f"{output_dir}/{angle1}/{angle1}_{idx:06d}.{ext}",
                    f"{output_dir}/{angle2}/{angle2}_{idx:06d}.{ext}",
                    f"{output_dir}/FA/FA_{idx:06d}.{ext}",
                ]
            )
        if write_vectors and vector_format == "npy":
            expected_files.append(
                str(vector_field_path(output_dir, "npy") / f"eigen_vec_{idx:06d}.npy")
            )
        if not all(os.path.exists(path) for path in expected_files):
            return False
        if zarr_completed is not None and not zarr_completed[idx - start]:
            return False
    return True
