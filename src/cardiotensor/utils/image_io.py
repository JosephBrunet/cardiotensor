import os
from pathlib import Path
from typing import Callable

import cv2
import glymur
import numpy as np
import tifffile


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
    out_path = str(out_path)
    tmp_path = f"{out_path}.tmp"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)

    try:
        writer(tmp_path)
        os.replace(tmp_path, out_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


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
    vector_field_slice: np.ndarray, start_index: int, output_dir: str, slice_idx: int
) -> None:
    """Save one vector field slice as eigen_vec/eigen_vec_XXXXXX.npy."""
    eig_dir = Path(output_dir) / "eigen_vec"
    eig_dir.mkdir(parents=True, exist_ok=True)
    out_path = eig_dir / f"eigen_vec_{(start_index + slice_idx):06d}.npy"

    def _writer(tmp_path: str) -> None:
        with open(tmp_path, "wb") as f:
            np.save(f, vector_field_slice)

    write_atomic(out_path, _writer)
