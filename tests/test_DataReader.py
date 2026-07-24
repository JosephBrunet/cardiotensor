from pathlib import Path

import cv2
import numpy as np
import pytest
import tifffile

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.image_io import initialize_zarr_vector_field


@pytest.fixture
def test_stack(tmp_path: Path):
    # Create a fake image stack
    for i in range(5):
        img = np.full((10, 10), i * 50, dtype=np.uint8)
        cv2.imwrite(str(tmp_path / f"img_{i:03d}.tif"), img)
    return tmp_path


def test_read_volume_shape(test_stack):
    reader = DataReader(test_stack)
    vol = reader.load_volume()
    assert vol.shape[0] == 5
    assert vol.ndim == 3


def test_mask_optional(test_stack):
    reader = DataReader(test_stack)
    # Should not fail without mask
    vol = reader.load_volume()
    assert isinstance(vol, np.ndarray)


def test_read_single_tiff_stack(tmp_path: Path):
    stack = np.stack(
        [np.full((10, 10), i * 50, dtype=np.uint8) for i in range(5)], axis=0
    )
    stack_path = tmp_path / "stack.tif"
    tifffile.imwrite(stack_path, stack)

    reader = DataReader(stack_path)
    vol = reader.load_volume()

    assert vol.shape == (5, 10, 10)
    assert np.array_equal(vol[2], stack[2])


def test_read_single_tiff_stack_loads_only_requested_pages(tmp_path: Path, monkeypatch):
    stack = np.arange(6 * 8 * 9, dtype=np.uint16).reshape(6, 8, 9)
    stack_path = tmp_path / "stack.tif"
    tifffile.imwrite(stack_path, stack, photometric="minisblack")

    requested_keys = []
    original_asarray = tifffile.TiffFile.asarray

    def record_key(self, key=None, **kwargs):
        requested_keys.append(key)
        return original_asarray(self, key=key, **kwargs)

    monkeypatch.setattr(tifffile.TiffFile, "asarray", record_key)
    volume = DataReader(stack_path).load_volume(2, 5)

    np.testing.assert_array_equal(volume, stack[2:5])
    assert requested_keys == [slice(2, 5)]


def test_single_tiff_fallback_memmap_uses_data_scratch(tmp_path: Path, monkeypatch):
    stack = np.arange(3 * 4 * 5, dtype=np.uint16).reshape(3, 4, 5)
    stack_path = tmp_path / "stack.tif"
    tifffile.imwrite(stack_path, stack, photometric="minisblack")
    reader = DataReader(stack_path)
    captured = {}

    class FakeSeries:
        pages = [object()]

        def asarray(self, *, out):
            captured["out"] = Path(out)
            mapped = np.memmap(
                out, mode="w+", dtype=stack.dtype, shape=stack.shape
            )
            mapped[:] = stack
            return mapped

    class FakeTiff:
        series = [FakeSeries()]

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

    monkeypatch.setattr(tifffile, "TiffFile", lambda path: FakeTiff())

    volume = reader.load_volume(1, 3)

    np.testing.assert_array_equal(volume, stack[1:3])
    assert captured["out"].parent.parent == tmp_path / ".cardiotensor_scratch"
    assert not captured["out"].parent.exists()


def test_mhd_reads_requested_z_range_as_copy_on_write_memmap(tmp_path: Path):
    volume = np.arange(5 * 6 * 7, dtype=np.int16).reshape(5, 6, 7)
    raw_path = tmp_path / "volume.raw"
    raw_path.write_bytes(volume.astype("<i2").tobytes())
    mhd_path = tmp_path / "volume.mhd"
    mhd_path.write_text(
        "ObjectType = Image\n"
        "NDims = 3\n"
        "BinaryData = True\n"
        "BinaryDataByteOrderMSB = False\n"
        "DimSize = 7 6 5\n"
        "ElementType = MET_SHORT\n"
        "ElementDataFile = volume.raw\n"
    )

    reader = DataReader(mhd_path)
    loaded = reader.load_volume(1, 4)

    assert reader.shape == volume.shape
    assert isinstance(loaded, np.memmap)
    np.testing.assert_array_equal(loaded, volume[1:4])

    loaded[0] = 0
    unchanged = np.memmap(raw_path, dtype="<i2", mode="r", shape=volume.shape)
    np.testing.assert_array_equal(unchanged, volume)


def test_zarr_region_read_does_not_load_full_slices(tmp_path: Path, monkeypatch):
    store = initialize_zarr_vector_field(tmp_path, (3, 6, 7))
    full = np.empty((3, 3, 6, 7), dtype=np.float32)
    for z in range(3):
        full[:, z] = z + np.arange(3, dtype=np.float32)[:, None, None]
        store.write_slice(full[:, z], z)

    reader = DataReader(store.path)

    def fail_full_load(*args, **kwargs):
        raise AssertionError("load_region should read Zarr directly")

    monkeypatch.setattr(reader, "load_volume", fail_full_load)
    region = reader.load_region(1, 3, start_y=2, end_y=5, start_x=1, end_x=6)

    np.testing.assert_array_equal(region, full[:, 1:3, 2:5, 1:6])
