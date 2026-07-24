#!/usr/bin/env python3
"""Create an RGB movie from a scalar image stack."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterator
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
from alive_progress import alive_bar

from cardiotensor.colormaps.helix_angle import helix_angle_cmap
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.image_io import read_image_file
from cardiotensor.utils.utils import read_conf_file

ANGLE_LIMITS = {
    "HA": (-90.0, 90.0),
    "IA": (-90.0, 90.0),
    "EL": (-90.0, 90.0),
    "AZ": (0.0, 360.0),
    "FA": (0.0, 1.0),
}


def _resolve_stack_path(input_path: Path, volume: str) -> Path:
    input_path = input_path.resolve()
    volume = volume.upper()

    if input_path.suffix.lower() == ".conf":
        params = read_conf_file(str(input_path))
        return (Path(params.get("OUTPUT_PATH", "./output")) / volume).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    if input_path.is_dir() and (input_path / volume).is_dir():
        return (input_path / volume).resolve()

    return input_path


def _iter_slices(
    reader: DataReader, start: int, end: int, step: int
) -> Iterator[np.ndarray]:
    if reader.volume_info["stack"]:
        files = reader.volume_info["file_list"]
        for path in files[start:end:step]:
            frame = read_image_file(path)
            if frame.ndim == 3 and frame.shape[2] == 1:
                frame = frame[:, :, 0]
            yield frame
        return

    for z in range(start, end, step):
        yield reader.load_volume(start_index=z, end_index=z + 1)[0]


def _normalise_scalar(
    frame: np.ndarray,
    *,
    volume: str,
    vmin: float | None,
    vmax: float | None,
) -> np.ndarray:
    original_dtype = frame.dtype
    frame = frame.astype(np.float32, copy=False)
    finite = np.isfinite(frame)

    if vmin is None or vmax is None:
        if original_dtype == np.uint8 or (
            finite.any() and float(np.nanmax(frame)) > 1.5
        ):
            norm = frame / 255.0
            return np.clip(np.nan_to_num(norm, nan=0.0), 0.0, 1.0)

        default_limits = ANGLE_LIMITS.get(volume.upper())
        if default_limits is not None:
            auto_vmin, auto_vmax = default_limits
        elif finite.any():
            auto_vmin = float(np.nanpercentile(frame[finite], 1.0))
            auto_vmax = float(np.nanpercentile(frame[finite], 99.0))
        else:
            auto_vmin, auto_vmax = 0.0, 1.0

        vmin = auto_vmin if vmin is None else vmin
        vmax = auto_vmax if vmax is None else vmax

    if vmax <= vmin:
        return np.zeros(frame.shape, dtype=np.float32)

    norm = (frame - vmin) / (vmax - vmin)
    return np.clip(np.nan_to_num(norm, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)


def _to_rgb_frame(
    frame: np.ndarray,
    *,
    volume: str,
    colormap: str,
    vmin: float | None,
    vmax: float | None,
) -> np.ndarray:
    if frame.ndim == 3 and frame.shape[2] in {3, 4}:
        rgb = frame[:, :, :3]
        if rgb.dtype != np.uint8:
            rgb = _normalise_scalar(rgb, volume=volume, vmin=vmin, vmax=vmax) * 255
        return rgb.astype(np.uint8)

    norm = _normalise_scalar(frame, volume=volume, vmin=vmin, vmax=vmax)
    cmap = helix_angle_cmap if colormap == "helix_angle" else plt.get_cmap(colormap)
    return (cmap(norm)[:, :, :3] * 255).astype(np.uint8)


def _resize_frame(
    frame: np.ndarray,
    *,
    width: int | None,
    height: int | None,
    scale: float,
) -> np.ndarray:
    if scale <= 0:
        raise ValueError("--scale must be > 0")

    src_h, src_w = frame.shape[:2]
    dst_w = int(round(src_w * scale)) if width is None else width
    dst_h = int(round(src_h * scale)) if height is None else height

    if width is None and height is not None:
        dst_w = int(round(src_w * (height / src_h)))
    elif height is None and width is not None:
        dst_h = int(round(src_h * (width / src_w)))

    if (dst_w, dst_h) == (src_w, src_h):
        return frame

    return cv2.resize(frame, (dst_w, dst_h), interpolation=cv2.INTER_AREA)


def create_movie(
    input_path: Path,
    *,
    volume: str = "HA",
    output: Path | None = None,
    fps: float = 30.0,
    start: int = 0,
    end: int | None = None,
    step: int = 1,
    width: int | None = None,
    height: int | None = None,
    scale: float = 1.0,
    colormap: str = "helix_angle",
    vmin: float | None = None,
    vmax: float | None = None,
) -> Path:
    stack_path = _resolve_stack_path(input_path, volume)
    reader = DataReader(stack_path)

    if len(reader.shape) not in {3, 4}:
        raise ValueError(
            f"Expected a 3D stack, got shape {reader.shape} at {stack_path}"
        )

    n_slices = reader.shape[0]
    start = max(0, start)
    end = n_slices if end is None else min(end, n_slices)
    if step < 1:
        raise ValueError("--step must be >= 1")
    if start >= end:
        raise ValueError(f"Empty slice range: start={start}, end={end}")

    if output is None:
        output = stack_path.parent / f"{stack_path.name}_movie.mp4"
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    first_slice = next(_iter_slices(reader, start, start + 1, 1))
    first_rgb = _to_rgb_frame(
        first_slice, volume=volume, colormap=colormap, vmin=vmin, vmax=vmax
    )
    first_rgb = _resize_frame(first_rgb, width=width, height=height, scale=scale)
    out_h, out_w = first_rgb.shape[:2]

    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(str(output), fourcc, fps, (out_w, out_h), isColor=True)
    if not writer.isOpened():
        raise RuntimeError(f"Could not open video writer for {output}")

    print(f"Creating RGB movie from {stack_path}")
    print(f"Slices: {start}:{end}:{step} | Size: {out_w}x{out_h} | FPS: {fps:g}")
    try:
        total_frames = len(range(start, end, step))
        with alive_bar(total_frames, title="Writing movie", length=40) as bar:
            for frame in _iter_slices(reader, start, end, step):
                rgb = _to_rgb_frame(
                    frame, volume=volume, colormap=colormap, vmin=vmin, vmax=vmax
                )
                rgb = _resize_frame(rgb, width=width, height=height, scale=scale)
                writer.write(cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR))
                bar()
    finally:
        writer.release()

    print(f"Movie created: {output}")
    return output


def script() -> None:
    parser = argparse.ArgumentParser(
        description="Create an RGB MP4 movie from a scalar stack such as HA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input",
        type=Path,
        help="A .conf file, an OUTPUT_PATH directory, or a direct stack directory/file.",
    )
    parser.add_argument(
        "--volume",
        default="HA",
        help="Subfolder to use when input is a .conf or output directory.",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output movie path.")
    parser.add_argument("--fps", type=float, default=30.0, help="Frames per second.")
    parser.add_argument("--start", type=int, default=0, help="First slice index.")
    parser.add_argument("--end", type=int, default=None, help="Stop slice index.")
    parser.add_argument("--step", type=int, default=1, help="Slice step.")
    parser.add_argument(
        "--width", type=int, default=None, help="Output width in pixels."
    )
    parser.add_argument(
        "--height", type=int, default=None, help="Output height in pixels."
    )
    parser.add_argument("--scale", type=float, default=1.0, help="Output scale factor.")
    parser.add_argument(
        "--colormap",
        default="helix_angle",
        help="Matplotlib colormap name, or 'helix_angle' for the project HA map.",
    )
    parser.add_argument("--vmin", type=float, default=None, help="Scalar value for 0.")
    parser.add_argument("--vmax", type=float, default=None, help="Scalar value for 1.")
    args = parser.parse_args()

    try:
        create_movie(
            args.input,
            volume=args.volume,
            output=args.output,
            fps=args.fps,
            start=args.start,
            end=args.end,
            step=args.step,
            width=args.width,
            height=args.height,
            scale=args.scale,
            colormap=args.colormap,
            vmin=args.vmin,
            vmax=args.vmax,
        )
    except Exception as exc:
        sys.exit(f"Error: {exc}")


if __name__ == "__main__":
    script()
