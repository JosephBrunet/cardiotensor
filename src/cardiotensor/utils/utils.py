import ast
import configparser
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import psutil


def get_available_cpu_count(default: int = 1) -> int:
    """Return the CPU count available to this process, respecting SLURM limits."""
    for env_name in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        value = os.environ.get(env_name)
        if value:
            try:
                count = int(value)
            except ValueError:
                continue
            if count > 0:
                return count

    return os.cpu_count() or default


def _memory_available_from_cgroup(
    proc_cgroup_path: Path = Path("/proc/self/cgroup"),
    cgroup_root: Path = Path("/sys/fs/cgroup"),
) -> int | None:
    """Return unused bytes in the process memory cgroup, when available."""

    def read_pair(directory: Path, limit_name: str, usage_name: str) -> int | None:
        try:
            limit_text = (directory / limit_name).read_text().strip()
            usage = int((directory / usage_name).read_text().strip())
        except (OSError, ValueError):
            return None
        if limit_text == "max":
            return None
        try:
            limit = int(limit_text)
        except ValueError:
            return None
        return max(limit - usage, 0)

    try:
        entries = proc_cgroup_path.read_text().splitlines()
    except OSError:
        entries = []

    for entry in entries:
        parts = entry.split(":", maxsplit=2)
        if len(parts) != 3:
            continue
        _, controllers, relative_path = parts
        relative_path = relative_path.lstrip("/")

        if not controllers:  # cgroup v2
            available = read_pair(
                cgroup_root / relative_path, "memory.max", "memory.current"
            )
            if available is not None:
                return available
        elif "memory" in controllers.split(","):  # cgroup v1
            for base in (cgroup_root / "memory", cgroup_root):
                available = read_pair(
                    base / relative_path,
                    "memory.limit_in_bytes",
                    "memory.usage_in_bytes",
                )
                if available is not None:
                    return available

    # Common when the process sees a cgroup namespace rooted at its own group.
    available = read_pair(cgroup_root, "memory.max", "memory.current")
    if available is not None:
        return available
    return read_pair(
        cgroup_root / "memory",
        "memory.limit_in_bytes",
        "memory.usage_in_bytes",
    )


def _memory_available_from_slurm() -> int | None:
    """Estimate unused allocation from SLURM variables when cgroups are hidden."""
    memory_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if memory_mb is None:
        memory_per_cpu_mb = os.environ.get("SLURM_MEM_PER_CPU")
        if memory_per_cpu_mb is not None:
            try:
                memory_mb = str(
                    int(memory_per_cpu_mb) * get_available_cpu_count(default=1)
                )
            except ValueError:
                return None
    if memory_mb is None:
        return None
    try:
        limit = int(memory_mb) * 1024**2
    except ValueError:
        return None
    return max(limit - psutil.Process().memory_info().rss, 0)


def get_available_memory_bytes() -> int:
    """Return free memory visible to this process, respecting job limits."""
    host_available = int(psutil.virtual_memory().available)
    job_available = _memory_available_from_cgroup()
    if job_available is None:
        job_available = _memory_available_from_slurm()
    return (
        host_available if job_available is None else min(host_available, job_available)
    )


def get_gpu_count() -> int:
    """Return the number of usable visible NVIDIA GPUs.

    CuPy is deliberately exercised in a subprocess. Initializing CUDA in this
    process before ``structure_tensor`` forks its GPU workers leaves the workers
    with an invalid inherited CUDA context on Linux.
    """

    def _cupy_gpu_count(max_count: int | None = None) -> int:
        probe = """
import sys
import cupy as cp

limit = int(sys.argv[1])
count = int(cp.cuda.runtime.getDeviceCount())
if limit >= 0:
    count = min(count, limit)

usable = 0
for device_id in range(count):
    try:
        with cp.cuda.Device(device_id):
            cp.arange(2, dtype=cp.float32).sum().item()
        usable += 1
    except Exception:
        pass
print(usable)
"""
        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    probe,
                    str(max_count if max_count is not None else -1),
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=60,
            )
            return int(result.stdout.strip().splitlines()[-1])
        except (OSError, subprocess.SubprocessError, ValueError, IndexError):
            return 0

    if os.environ.get("SLURM_JOB_ID"):
        allocated = os.environ.get("SLURM_JOB_GPUS") or os.environ.get(
            "SLURM_STEP_GPUS"
        )
        if not allocated:
            return 0
        count = len([x for x in allocated.split(",") if x.strip()])
        return _cupy_gpu_count(max_count=count)

    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if visible:
        ids = [x for x in visible.split(",") if x.strip()]
        if ids:
            return _cupy_gpu_count(max_count=len(ids))

    try:
        if platform.system() == "Windows":
            result = subprocess.run(
                [r"C:\Program Files\NVIDIA Corporation\NVSMI\nvidia-smi.exe", "-L"],
                capture_output=True,
                text=True,
                check=True,
            )
        else:
            result = subprocess.run(
                ["nvidia-smi", "-L"],
                capture_output=True,
                text=True,
                check=True,
            )
        count = len(
            [line for line in result.stdout.strip().splitlines() if "GPU" in line]
        )
        return _cupy_gpu_count(max_count=count)
    except Exception:
        return 0


def read_conf_file(file_path: str) -> dict[str, Any]:
    """
    Reads and parses a configuration file into a dictionary.

    Args:
        file_path (str): Path to the configuration file.

    Returns:
        Dict[str, Any]: Parsed configuration parameters.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
        ValueError: If expected numerical or array values are incorrectly formatted.
    """

    config_path = Path(file_path).resolve()
    file_path = str(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"The configuration file {file_path} does not exist.")

    if not file_path.endswith(".conf"):
        raise ValueError("The file is not a .conf file")

    config = configparser.ConfigParser()
    config.read(file_path)

    def resolve_path(value: str) -> str:
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = config_path.parent / path
        return str(path.resolve())

    def parse_coordinates(section: str, option: str, fallback: str = ""):
        """
        Parses a coordinate string from a configuration file into a list of tuples.

        Parameters:
            config (ConfigParser): The configuration parser object.
            section (str): The section in the config.
            option (str): The option name.
            fallback (str): A fallback value if the option is missing.

        Returns:
            list[tuple[int, int, int]]: List of 3D coordinate tuples or an empty list.
        """
        value = config.get(section, option, fallback=fallback).strip().replace(" ", "")

        if not value:  # Return an empty list if the value is empty
            return []

        try:
            value = f"[{value}]"  # Ensure it's formatted as a list
            points_list = ast.literal_eval(value)  # Safely evaluate the string
            return [tuple(point) for point in points_list]  # Convert to list of tuples
        except (SyntaxError, ValueError) as e:
            raise ValueError(
                f"Invalid coordinate format for {option} in [{section}]: {value}"
            ) from e

    # Read the two paths
    images_path = resolve_path(config.get("DATASET", "IMAGES_PATH").strip())
    mask_path = config.get("DATASET", "MASK_PATH", fallback=None)
    if mask_path == "":
        mask_path = None
    if mask_path is not None:
        mask_path = resolve_path(mask_path.strip())

    # Existence check (file or directory)
    if not os.path.exists(images_path):
        raise FileNotFoundError(f"The IMAGES_PATH '{images_path}' does not exist.")
    if mask_path and not os.path.exists(mask_path):
        raise FileNotFoundError(f"The MASK_PATH '{mask_path}' does not exist.")

    # OUTPUT
    output_dir = config.get("OUTPUT", "OUTPUT_PATH", fallback="").strip()
    if not output_dir:
        output_dir = "./output"  # Default output directory if not specified
    output_dir = resolve_path(output_dir)
    low_memory_dir = config.get(
        "STRUCTURE TENSOR CALCULATION", "LOW_MEMORY_DIR", fallback=""
    ).strip()
    low_memory_dir = (
        resolve_path(os.path.expandvars(low_memory_dir)) if low_memory_dir else None
    )

    params = {
        # DATASET
        "IMAGES_PATH": images_path,
        "MASK_PATH": mask_path,
        "VOXEL_SIZE": config.getfloat("DATASET", "VOXEL_SIZE", fallback=1.0),
        # STRUCTURE TENSOR CALCULATION
        "SIGMA": config.getfloat("STRUCTURE TENSOR CALCULATION", "SIGMA", fallback=1.0),
        "RHO": config.getfloat("STRUCTURE TENSOR CALCULATION", "RHO", fallback=3.0),
        "TRUNCATE": config.getfloat(
            "STRUCTURE TENSOR CALCULATION", "TRUNCATE", fallback=4.0
        ),
        "VERTICAL_PADDING": config.getfloat(
            "STRUCTURE TENSOR CALCULATION", "VERTICAL_PADDING", fallback=None
        ),
        "N_CHUNK": config.getint(
            "STRUCTURE TENSOR CALCULATION", "N_CHUNK", fallback=100
        ),
        "USE_GPU": config.getboolean(
            "STRUCTURE TENSOR CALCULATION", "USE_GPU", fallback=False
        ),
        "GPU_WORKERS_PER_DEVICE": config.getint(
            "STRUCTURE TENSOR CALCULATION",
            "GPU_WORKERS_PER_DEVICE",
            fallback=1,
        ),
        "LOW_MEMORY": config.getboolean(
            "STRUCTURE TENSOR CALCULATION", "LOW_MEMORY", fallback=False
        ),
        "LOW_MEMORY_DIR": low_memory_dir,
        "WRITE_VECTORS": config.getboolean(
            "STRUCTURE TENSOR CALCULATION", "WRITE_VECTORS", fallback=False
        ),
        "VECTOR_FORMAT": config.get(
            "STRUCTURE TENSOR CALCULATION", "VECTOR_FORMAT", fallback="zarr"
        )
        .strip()
        .lower(),
        "REVERSE": config.getboolean(
            "STRUCTURE TENSOR CALCULATION", "REVERSE", fallback=False
        ),
        # ANGLE CALCULATION
        "WRITE_ANGLES": config.getboolean(
            "ANGLE CALCULATION", "WRITE_ANGLES", fallback=True
        ),
        "ANGLE_MODE": config.get(
            "ANGLE CALCULATION", "ANGLE_MODE", fallback="ha_ia"
        ).strip(),
        "PROJECTED_ANGLES": config.getboolean(
            "ANGLE CALCULATION",
            "PROJECTED_ANGLES",
            fallback=config.getboolean(
                "ANGLE CALCULATION", "PROJECTED", fallback=False
            ),
        ),
        "AXIS_POINTS": parse_coordinates("ANGLE CALCULATION", "AXIS_POINTS"),
        # TEST
        "TEST": config.getboolean("TEST", "TEST", fallback=False),
        "N_SLICE_TEST": config.getint("TEST", "N_SLICE_TEST", fallback=None),
        "SHOW_QUIVER": config.getboolean("TEST", "SHOW_QUIVER", fallback=True),
        # OUTPUT
        "OUTPUT_PATH": output_dir,
        "OUTPUT_FORMAT": config.get("OUTPUT", "OUTPUT_FORMAT", fallback="jp2").strip(),
        "OUTPUT_TYPE": config.get("OUTPUT", "OUTPUT_TYPE", fallback="8bit").strip(),
        "COLORMAP": config.get("OUTPUT", "COLORMAP", fallback=None),
        "COLORMAP_ANGLE1": config.get("OUTPUT", "COLORMAP_ANGLE1", fallback=None),
        "COLORMAP_ANGLE2": config.get("OUTPUT", "COLORMAP_ANGLE2", fallback=None),
    }
    if params["GPU_WORKERS_PER_DEVICE"] < 1:
        raise ValueError("GPU_WORKERS_PER_DEVICE must be at least 1")
    if not params["WRITE_VECTORS"] and not params["WRITE_ANGLES"]:
        raise ValueError("At least one of WRITE_VECTORS or WRITE_ANGLES must be True")
    return params


def remove_corrupted_files(file_paths, size_threshold=200):
    import warnings

    removed_files = []
    for file_path in file_paths:
        if os.path.exists(file_path) and os.path.getsize(file_path) < size_threshold:
            warnings.warn(
                f"Removing potentially corrupted file (size < {size_threshold} bytes): {file_path}",
                UserWarning,
                stacklevel=2,
            )
            os.remove(file_path)
            removed_files.append(file_path)

    return removed_files


def convert_to_8bit(
    img: np.ndarray,
    perc_min: int = 0,
    perc_max: int = 100,
    min_value: float | None = None,
    max_value: float | None = None,
) -> np.ndarray:
    """
    Converts a NumPy array to an 8-bit image.

    Args:
        img (np.ndarray): Input image array.
        perc_min (int): Minimum percentile for normalization. Default is 0.
        perc_max (int): Maximum percentile for normalization. Default is 100.
        min_value (Optional[float]): Optional explicit minimum value.
        max_value (Optional[float]): Optional explicit maximum value.

    Returns:
        np.ndarray: 8-bit converted image.
    """
    if min_value is not None and max_value is not None:
        minimum, maximum = min_value, max_value
    else:
        minimum, maximum = np.nanpercentile(img, (perc_min, perc_max))
        if min_value is not None:
            minimum = min_value
        if max_value is not None:
            maximum = max_value

    # Avoid division by zero
    if maximum == minimum:
        return np.zeros_like(img, dtype=np.uint8)

    # Normalize and scale to 8-bit range
    img_normalized = (img - minimum) / (maximum - minimum) * 255

    # Convert NaN and Inf values to numbers
    img_cleaned = np.nan_to_num(img_normalized, nan=0.0, posinf=255.0, neginf=0.0)

    # Clip values to ensure they are in 8-bit range
    img_clipped = np.clip(img_cleaned, 0, 255)

    return img_clipped.astype(np.uint8)
