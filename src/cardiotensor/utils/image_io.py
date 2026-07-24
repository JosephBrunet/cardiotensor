import math
import os
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import cv2
import glymur
import numpy as np
import tifffile
import zarr

VECTOR_FORMATS = {"npy", "zarr"}
ZARR_VECTOR_FORMAT_VERSION = 1
ZARR_VECTOR_CHUNK_SIZE = 512
ZARR_VECTOR_SHARD_CHUNKS = 8


def normalize_vector_format(vector_format: str) -> str:
    """Validate and normalize a vector output format."""
    normalized = vector_format.lower().strip()
    if normalized not in VECTOR_FORMATS:
        choices = ", ".join(sorted(VECTOR_FORMATS))
        raise ValueError(f"VECTOR_FORMAT must be one of: {choices}")
    return normalized


def vector_field_path(output_dir: str | Path, vector_format: str) -> Path:
    """Return the vector-field path for one configured storage format."""
    vector_format = normalize_vector_format(vector_format)
    name = "eigen_vec" if vector_format == "npy" else "eigen_vec.zarr"
    return Path(output_dir) / name


def _zarr_vector_layout(height: int, width: int) -> tuple[tuple[int, ...], ...]:
    """Return the fixed chunk and shard layout for a vector field."""
    chunk_y = min(ZARR_VECTOR_CHUNK_SIZE, height)
    chunk_x = min(ZARR_VECTOR_CHUNK_SIZE, width)
    shard_y = chunk_y * min(ZARR_VECTOR_SHARD_CHUNKS, math.ceil(height / chunk_y))
    shard_x = chunk_x * min(ZARR_VECTOR_SHARD_CHUNKS, math.ceil(width / chunk_x))
    return (3, 1, chunk_y, chunk_x), (3, 1, shard_y, shard_x)


def _codec_value(value: Any) -> Any:
    """Return the plain value stored by a codec enum."""
    return getattr(value, "value", value)


def _validate_zarr_vector_group(
    group: Any,
    volume_shape: tuple[int, int, int] | None = None,
    *,
    allow_missing_version: bool = False,
) -> None:
    """Validate arrays and metadata used by a Cardiotensor vector store."""
    if "vectors" not in group or "completed" not in group:
        raise ValueError("Zarr vector store must contain 'vectors' and 'completed'")

    vectors = group["vectors"]
    completed = group["completed"]
    if len(vectors.shape) != 4 or vectors.shape[0] != 3:
        raise ValueError(
            f"Zarr vectors must have shape (3, Z, Y, X), got {vectors.shape}"
        )

    z_count, height, width = (int(size) for size in vectors.shape[1:])
    if volume_shape is not None and (z_count, height, width) != tuple(volume_shape):
        raise ValueError(
            f"Existing Zarr volume shape {(z_count, height, width)} does not match "
            f"requested shape {tuple(volume_shape)}"
        )
    if np.dtype(vectors.dtype) != np.dtype(np.float32):
        raise ValueError(f"Zarr vector dtype must be float32, got {vectors.dtype}")

    expected_chunks, expected_shards = _zarr_vector_layout(height, width)
    if tuple(vectors.chunks) != expected_chunks:
        raise ValueError(
            f"Zarr vector chunks {vectors.chunks} do not match {expected_chunks}"
        )
    if vectors.shards is None or tuple(vectors.shards) != expected_shards:
        raise ValueError(
            f"Zarr vector shards {vectors.shards} do not match {expected_shards}"
        )
    if vectors.fill_value != 0:
        raise ValueError(f"Zarr vector fill value must be 0, got {vectors.fill_value}")
    dimension_names = getattr(vectors.metadata, "dimension_names", None)
    if dimension_names is None or tuple(dimension_names) != (
        "component",
        "z",
        "y",
        "x",
    ):
        raise ValueError("Zarr vector dimensions must be ('component', 'z', 'y', 'x')")

    compressors = vectors.compressors
    if len(compressors) != 1 or not isinstance(compressors[0], zarr.codecs.BloscCodec):
        raise ValueError("Zarr vectors must use one Blosc compressor")
    compressor = compressors[0]
    if (
        _codec_value(compressor.cname) != "zstd"
        or compressor.clevel != 3
        or _codec_value(compressor.shuffle) != "bitshuffle"
    ):
        raise ValueError("Zarr vectors must use Blosc Zstd level 3 with bit-shuffle")

    if tuple(completed.shape) != (z_count,):
        raise ValueError(
            f"Zarr completion shape {completed.shape} does not match {(z_count,)}"
        )
    if tuple(completed.chunks) != (1,) or np.dtype(completed.dtype) != np.dtype(bool):
        raise ValueError("Zarr completion markers must be bool with chunks=(1,)")
    if completed.fill_value not in (False, 0):
        raise ValueError("Zarr completion fill value must be False")

    expected_attrs = {
        "cardiotensor_format": "vector_field",
        "axis_order": ["component", "z", "y", "x"],
        "components": ["x", "y", "z"],
        "masked_fill_value": 0.0,
    }
    for name, expected in expected_attrs.items():
        if group.attrs.get(name) != expected:
            raise ValueError(
                f"Invalid Zarr metadata {name!r}: expected {expected!r}, "
                f"got {group.attrs.get(name)!r}"
            )
    version = group.attrs.get("cardiotensor_format_version")
    if version is None and allow_missing_version:
        return
    if version != ZARR_VECTOR_FORMAT_VERSION:
        raise ValueError(
            f"Unsupported Zarr vector format version {version!r}; "
            f"expected {ZARR_VECTOR_FORMAT_VERSION}"
        )


@contextmanager
def _zarr_metadata_lock(lock_path: Path):
    """Serialize Zarr metadata creation across local or SLURM processes."""
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a+b") as lock_file:
        try:
            import fcntl
        except ImportError:  # pragma: no cover - Windows fallback
            fcntl = None

        if fcntl is not None:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


@dataclass(frozen=True)
class ZarrVectorFieldStore:
    """Open arrays belonging to one Cardiotensor Zarr vector field."""

    path: Path
    vectors: Any
    completed: Any

    def write_slice(self, vector_slice: np.ndarray, global_slice_idx: int) -> None:
        expected_shape = (3, *self.vectors.shape[2:])
        if vector_slice.shape != expected_shape:
            raise ValueError(
                f"Vector slice shape {vector_slice.shape} does not match "
                f"Zarr slice shape {expected_shape}"
            )
        if not 0 <= global_slice_idx < self.vectors.shape[1]:
            raise IndexError(
                f"Vector slice {global_slice_idx} is outside Zarr bounds "
                f"[0, {self.vectors.shape[1]})"
            )

        self.vectors[:, global_slice_idx, :, :] = np.asarray(
            vector_slice, dtype=np.float32
        )
        # Completion is written last, so interrupted slices are recomputed.
        self.completed[global_slice_idx] = True

    def completed_range(self, start_index: int, end_index: int) -> np.ndarray:
        return np.asarray(self.completed[start_index:end_index], dtype=bool)


def initialize_zarr_vector_field(
    output_dir: str | Path,
    volume_shape: tuple[int, int, int],
) -> ZarrVectorFieldStore:
    """Create or validate the sharded Zarr vector field for a full volume."""
    z_count, height, width = (int(size) for size in volume_shape)
    if min(z_count, height, width) <= 0:
        raise ValueError(f"Invalid vector volume shape: {volume_shape}")

    store_path = vector_field_path(output_dir, "zarr")
    lock_path = Path(output_dir) / ".eigen_vec.zarr.lock"
    vector_shape = (3, z_count, height, width)
    chunks, shards = _zarr_vector_layout(height, width)

    with _zarr_metadata_lock(lock_path):
        group = zarr.open_group(
            store=str(store_path),
            mode="a",
            zarr_format=3,
        )
        existing_store = "vectors" in group and "completed" in group
        if existing_store:
            _validate_zarr_vector_group(group, volume_shape, allow_missing_version=True)

        if "vectors" not in group:
            vectors = group.create_array(
                "vectors",
                shape=vector_shape,
                chunks=chunks,
                shards=shards,
                dtype="float32",
                fill_value=0.0,
                compressors=zarr.codecs.BloscCodec(
                    cname="zstd",
                    clevel=3,
                    shuffle=zarr.codecs.BloscShuffle.bitshuffle,
                ),
                dimension_names=("component", "z", "y", "x"),
                config={"write_empty_chunks": False},
            )
        else:
            vectors = group["vectors"]

        if "completed" not in group:
            completed = group.create_array(
                "completed",
                shape=(z_count,),
                chunks=(1,),
                dtype="bool",
                fill_value=False,
                config={"write_empty_chunks": False},
            )
        else:
            completed = group["completed"]

        group.attrs.update(
            {
                "cardiotensor_format": "vector_field",
                "axis_order": ["component", "z", "y", "x"],
                "components": ["x", "y", "z"],
                "masked_fill_value": 0.0,
                "cardiotensor_format_version": ZARR_VECTOR_FORMAT_VERSION,
            }
        )
        _validate_zarr_vector_group(group, volume_shape)

    return ZarrVectorFieldStore(store_path, vectors, completed)


def open_zarr_vector_field(
    output_dir: str | Path,
    mode: str = "r",
) -> ZarrVectorFieldStore:
    """Open an existing Cardiotensor Zarr vector field."""
    store_path = vector_field_path(output_dir, "zarr")
    group = zarr.open_group(store=str(store_path), mode=mode)
    _validate_zarr_vector_group(group)
    return ZarrVectorFieldStore(store_path, group["vectors"], group["completed"])


def read_image_file(file_path: str | Path) -> np.ndarray:
    """Read one image file to a NumPy array using the preferred backend."""
    file_path = Path(file_path)
    suffix = file_path.suffix.lower()

    if suffix == ".npy":
        return np.load(file_path)

    if suffix in {".tif", ".tiff"}:
        try:
            return tifffile.imread(file_path)
        except Exception as e:
            img = cv2.imread(str(file_path), -1)
            if img is not None:
                print(
                    f"⚠ Falling back to OpenCV for TIFF '{file_path}' "
                    f"(tifffile error: {e})"
                )
                return img
            raise RuntimeError(
                "Failed to read TIFF with tifffile and OpenCV. "
                "If the file uses LZW/other codecs, install 'imagecodecs' "
                "to enable decoding."
            ) from e

    if suffix == ".jp2":
        return glymur.Jp2k(str(file_path))[:]

    img = cv2.imread(str(file_path), -1)
    if img is None:
        raise RuntimeError(
            f"cv2.imread failed for '{file_path}'. "
            "File may be missing or in an unsupported/invalid format."
        )
    return img


def write_atomic(out_path: str | Path, writer: Callable[[str], None]) -> None:
    """Write to a temporary file first, then replace the final output."""
    out_path = Path(out_path)
    if out_path.suffix:
        tmp_path = out_path.with_name(f"{out_path.stem}.tmp{out_path.suffix}")
    else:
        tmp_path = out_path.with_name(f"{out_path.name}.tmp")

    if tmp_path.exists():
        tmp_path.unlink()

    try:
        writer(str(tmp_path))
        os.replace(tmp_path, out_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def write_jp2(out_path: str | Path, data: np.ndarray) -> None:
    write_atomic(
        out_path,
        lambda tmp_path: glymur.Jp2k(
            tmp_path, data=data, cratios=[10], numres=8, irreversible=True
        ),
    )


def write_tif(out_path: str | Path, data: np.ndarray) -> None:
    write_atomic(out_path, lambda tmp_path: tifffile.imwrite(tmp_path, data))


def write_image(out_path: str | Path, data: np.ndarray) -> None:
    """Write an image using the backend selected from the output extension."""
    out_path = Path(out_path)
    suffix = out_path.suffix.lower()

    if suffix == ".jp2":
        write_jp2(out_path, data)
    elif suffix in {".tif", ".tiff"}:
        write_tif(out_path, data)
    else:
        if not cv2.imwrite(str(out_path), data):
            raise RuntimeError(f"cv2.imwrite failed for '{out_path}'")


def write_vector_field(
    vector_field_slice: np.ndarray,
    start_index: int,
    output_dir: str,
    slice_idx: int,
    vector_format: str = "npy",
    zarr_store: ZarrVectorFieldStore | None = None,
) -> None:
    """Save one vector field slice using the configured vector backend."""
    vector_format = normalize_vector_format(vector_format)
    global_slice_idx = start_index + slice_idx
    if vector_format == "zarr":
        if zarr_store is None:
            raise ValueError("zarr_store is required when vector_format='zarr'")
        zarr_store.write_slice(vector_field_slice, global_slice_idx)
        return

    eig_dir = vector_field_path(output_dir, "npy")
    eig_dir.mkdir(parents=True, exist_ok=True)
    out_path = eig_dir / f"eigen_vec_{global_slice_idx:06d}.npy"

    def _writer(tmp_path: str) -> None:
        with open(tmp_path, "wb") as f:
            np.save(f, vector_field_slice)

    write_atomic(out_path, _writer)
