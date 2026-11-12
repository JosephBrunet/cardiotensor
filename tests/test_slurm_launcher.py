from pathlib import Path

import numpy as np

from cardiotensor.launcher.slurm_launcher import monitor_job_output, slurm_launcher
from cardiotensor.utils.image_io import initialize_zarr_vector_field


def _create_angle_outputs(output_dir: Path, names: tuple[str, ...], count: int):
    for name in names:
        folder = output_dir / name
        folder.mkdir(parents=True)
        for index in range(count):
            (folder / f"{name}_{index:06d}.tif").touch()


def test_monitor_waits_for_zarr_when_angles_are_complete(
    tmp_path: Path, monkeypatch, capsys
):
    _create_angle_outputs(tmp_path, ("HA", "IA", "FA"), count=2)
    store = initialize_zarr_vector_field(tmp_path, (2, 4, 5))
    store.completed[0] = True

    sleep_calls = []

    def complete_zarr_after_first_poll(seconds):
        sleep_calls.append(seconds)
        store.completed[1] = True

    monkeypatch.setattr(
        "cardiotensor.launcher.slurm_launcher.time.sleep",
        complete_zarr_after_first_poll,
    )

    monitor_job_output(
        output_directory=str(tmp_path),
        start_index=0,
        end_index_exclusive=2,
        output_format="tif",
        write_angles=True,
        write_vectors=True,
        vector_format="zarr",
        poll_interval_sec=1,
    )

    output = capsys.readouterr().out
    assert "HA: 2/2 | IA: 2/2 | FA: 2/2 | Zarr vectors: 1/2" in output
    assert "HA: 2/2 | IA: 2/2 | FA: 2/2 | Zarr vectors: 2/2" in output
    assert sleep_calls == [1]


def test_monitor_supports_projected_angle_names(tmp_path: Path, capsys):
    _create_angle_outputs(tmp_path, ("HA_projected", "IA_projected", "FA"), count=1)

    monitor_job_output(
        output_directory=str(tmp_path),
        start_index=0,
        end_index_exclusive=1,
        output_format="tif",
        write_angles=True,
        write_vectors=False,
        projected=True,
    )

    output = capsys.readouterr().out
    assert "HA_projected: 1/1 | IA_projected: 1/1 | FA: 1/1" in output


def test_slurm_low_memory_dry_run_reduces_unsafe_chunk(
    tmp_path: Path, monkeypatch, capsys
):
    params = {
        "IMAGES_PATH": "large-volume",
        "MASK_PATH": "large-mask",
        "OUTPUT_PATH": str(tmp_path / "output"),
        "OUTPUT_FORMAT": "jp2",
        "ANGLE_MODE": "ha_ia",
        "WRITE_ANGLES": True,
        "WRITE_VECTORS": True,
        "VECTOR_FORMAT": "npy",
        "LOW_MEMORY": True,
        "SIGMA": 0.6,
        "RHO": 4.0,
        "TRUNCATE": 2.0,
        "VERTICAL_PADDING": 9,
        "N_CHUNK": 50,
        "TEST": False,
    }

    class FakeReader:
        shape = (10_920, 7_660, 7_385)
        dtype = np.dtype(np.uint16)

        def __init__(self, path):
            pass

    monkeypatch.setattr(
        "cardiotensor.launcher.slurm_launcher.read_conf_file", lambda _: params
    )
    monkeypatch.setattr("cardiotensor.launcher.slurm_launcher.DataReader", FakeReader)

    slurm_launcher(
        "parameters.conf",
        start_index=0,
        end_index=100,
        chunk_size=50,
        cpus_per_task=8,
        mem_gb=64,
        array_parallel=4,
        log_dir=str(tmp_path / "logs"),
        submit_dir=str(tmp_path / "submit"),
        monitor=False,
        dry_run=True,
    )

    output = capsys.readouterr().out
    assert "adjusted SLURM chunk_size from 50 to 16" in output
    script = next((tmp_path / "submit").glob("*.slurm")).read_text()
    assert "IMAGES_PER_JOB=16" in script
    assert "echo TMPDIR: ${TMPDIR:-/tmp}" in script


def test_slurm_skips_complete_chunks(tmp_path: Path, monkeypatch, capsys):
    output_dir = tmp_path / "output"
    _create_angle_outputs(output_dir, ("HA", "IA", "FA"), count=2)
    params = {
        "IMAGES_PATH": "volume",
        "OUTPUT_PATH": str(output_dir),
        "OUTPUT_FORMAT": "tif",
        "ANGLE_MODE": "ha_ia",
        "WRITE_ANGLES": True,
        "WRITE_VECTORS": False,
        "N_CHUNK": 2,
        "TEST": False,
    }

    class FakeReader:
        shape = (4, 8, 8)
        dtype = np.dtype(np.uint16)

        def __init__(self, path):
            pass

    monkeypatch.setattr(
        "cardiotensor.launcher.slurm_launcher.read_conf_file", lambda _: params
    )
    monkeypatch.setattr("cardiotensor.launcher.slurm_launcher.DataReader", FakeReader)

    slurm_launcher(
        "parameters.conf",
        chunk_size=2,
        log_dir=str(tmp_path / "logs"),
        submit_dir=str(tmp_path / "submit"),
        monitor=False,
        dry_run=True,
    )

    output = capsys.readouterr().out
    assert "1 complete job(s) skipped, 1 job(s) needed" in output
    scripts = list((tmp_path / "submit").glob("*.slurm"))
    assert len(scripts) == 1
    assert "START_INDEX_BASE=2" in scripts[0].read_text()


def test_slurm_gpu_request_is_written_to_script(tmp_path, monkeypatch, capsys):
    params = {
        "IMAGES_PATH": "volume",
        "OUTPUT_PATH": str(tmp_path / "output"),
        "OUTPUT_FORMAT": "tif",
        "ANGLE_MODE": "ha_ia",
        "WRITE_ANGLES": True,
        "WRITE_VECTORS": False,
        "USE_GPU": True,
        "N_CHUNK": 2,
        "TEST": False,
    }

    class FakeReader:
        shape = (4, 8, 8)
        dtype = np.dtype(np.uint16)

        def __init__(self, path):
            pass

    monkeypatch.setattr(
        "cardiotensor.launcher.slurm_launcher.read_conf_file", lambda _: params
    )
    monkeypatch.setattr(
        "cardiotensor.launcher.slurm_launcher.DataReader", FakeReader
    )

    slurm_launcher(
        "parameters.conf",
        chunk_size=2,
        gpus=2,
        log_dir=str(tmp_path / "logs"),
        submit_dir=str(tmp_path / "submit"),
        monitor=False,
        dry_run=True,
    )

    script = next((tmp_path / "submit").glob("*.slurm")).read_text()
    assert "#SBATCH --gres=gpu:2" in script
    assert "echo SLURM_JOB_GPUS: ${SLURM_JOB_GPUS:-<unset>}" in script
    assert 'echo "Requested GPUs:          2"' in script
    assert "gpus=2" in capsys.readouterr().out


def test_slurm_zero_gpus_omits_gpu_resource(tmp_path, monkeypatch):
    params = {
        "IMAGES_PATH": "volume",
        "OUTPUT_PATH": str(tmp_path / "output"),
        "OUTPUT_FORMAT": "tif",
        "ANGLE_MODE": "ha_ia",
        "WRITE_ANGLES": True,
        "WRITE_VECTORS": False,
        "N_CHUNK": 2,
        "TEST": False,
    }

    class FakeReader:
        shape = (2, 8, 8)
        dtype = np.dtype(np.uint16)

        def __init__(self, path):
            pass

    monkeypatch.setattr(
        "cardiotensor.launcher.slurm_launcher.read_conf_file", lambda _: params
    )
    monkeypatch.setattr(
        "cardiotensor.launcher.slurm_launcher.DataReader", FakeReader
    )

    slurm_launcher(
        "parameters.conf",
        gpus=0,
        log_dir=str(tmp_path / "logs"),
        submit_dir=str(tmp_path / "submit"),
        monitor=False,
        dry_run=True,
    )

    script = next((tmp_path / "submit").glob("*.slurm")).read_text()
    assert "#SBATCH --gres=gpu:" not in script
