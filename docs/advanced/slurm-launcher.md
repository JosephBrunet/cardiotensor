# SLURM Launcher

Use `cardio-tensor-slurm` to submit chunked orientation jobs as SLURM array tasks.

## What It Does

- Splits a dataset into chunks.
- Generates a `.slurm` script.
- Submits one or more SLURM arrays with `sbatch`.
- Monitors output progress.

The launcher calls `cardio-tensor` on each task with:

```bash
cardio-tensor <conf_file> --start_index <start> --end_index <end>
```

## Basic Usage

```bash
cardio-tensor-slurm path/to/parameters.conf
```

## Useful Options

```bash
cardio-tensor-slurm path/to/parameters.conf \
  --start_index 0 \                     # First slice index (inclusive)
  --end_index 5000 \                    # Last slice index (exclusive)
  --chunk_size 100 \                    # Slices processed per SLURM task
  --time_limit 4:00:00 \                # Wall-time limit for each task
  --cpus_per_task 8 \                   # CPU cores requested per task
  --mem_gb 64 \                         # Memory requested per task (GB)
  --array_parallel 50 \                 # Max concurrent tasks in the SLURM array
  --partition nice \                    # Optional partition/queue name
  --log_dir /path/to/logs \             # Where SLURM stdout/stderr logs are written
  --submit_dir /path/to/slurm_scripts   # Where generated .slurm scripts are saved
```

Other useful flags:

- `--no_monitor`: submit jobs and return immediately (do not watch output progress).
- `--dry_run`: generate the `.slurm` script and print `sbatch` command without submitting.

Run help for the full list:

```bash
cardio-tensor-slurm -h
```

## Troubleshooting

- **`sbatch` fails**
  Verify partition/account/time/memory policy on your cluster.

- **Job fails before running `cardio-tensor`**
  Check logs in `OUTPUT_PATH/slurm/log/`.

- **No progress in monitor**
  Confirm output folder permissions and that tasks can write outputs.

- **Want to debug without submitting jobs**
  Use `--dry_run`.
