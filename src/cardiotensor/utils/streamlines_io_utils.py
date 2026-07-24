from __future__ import annotations

from pathlib import Path

import numpy as np


def load_npz_streamlines(
    p: Path,
) -> tuple[list[np.ndarray], dict[str, list[np.ndarray]]]:
    """
    Load streamlines from a .npz file saved as object arrays.
    Expects 'streamlines' in (z, y, x). Converts to (x, y, z).
    Collects any per-point arrays whose keys end with '_values' and
    exposes them as uppercase names without the suffix, e.g. 'ha_values' -> 'HA'.

    Returns:
        streamlines_xyz: list[np.ndarray], each (N_i, 3) in (x, y, z)
        per_point: dict[str, list[np.ndarray]] keyed by field, each list aligned to streamlines
    """
    data = np.load(p, allow_pickle=True)
    raw_streamlines = data.get("streamlines")
    if raw_streamlines is None:
        raise ValueError("'streamlines' array missing in .npz")

    # stored as (z, y, x) convert to (x, y, z)
    streamlines_xyz: list[np.ndarray] = [
        np.asarray([(pt[2], pt[1], pt[0]) for pt in sl], dtype=np.float32)
        for sl in raw_streamlines.tolist()
    ]

    per_point: dict[str, list[np.ndarray]] = {}
    for key in data.files:
        if key == "streamlines":
            continue
        if key.endswith("_values"):
            base = key[:-7].upper()  # remove '_values'
            obj = data[key]
            vals = [np.asarray(a).reshape(-1) for a in obj.tolist()]
            if len(vals) != len(streamlines_xyz):
                raise ValueError(
                    f"{key} length {len(vals)} does not match streamlines {len(streamlines_xyz)}"
                )
            per_point[base] = vals

    return streamlines_xyz, per_point


def load_trk_streamlines(
    p: Path,
    include_per_streamline: bool = False,
):
    """
    Load streamlines and all per-point fields from a TrackVis .trk file.
    Returns streamlines in (x, y, z) voxel/world space (as stored in the TRK),
    and a dict of per-point fields, one list per field aligned with streamlines.
    """
    import nibabel as nib

    obj = nib.streamlines.load(str(p))
    tg = obj.tractogram

    # Prefer RAS mm if available
    try:
        tg = tg.to_rasmm()
    except AttributeError:
        try:
            tg = tg.to_world()
        except Exception:
            # Neither to_rasmm() nor to_world() is available; use coordinates as-is
            pass

    streamlines_xyz = [np.asarray(sl, dtype=np.float32) for sl in tg.streamlines]

    per_point: dict[str, list[np.ndarray]] = {}
    dpp = getattr(tg, "data_per_point", None)
    if dpp:
        for name, arrseq in dpp.items():
            # arrseq is an ArraySequence with shape (Ni, C). We flatten to 1D per point if C==1.
            vals: list[np.ndarray] = []
            for a in arrseq:
                a = np.asarray(a)
                if a.ndim == 2 and a.shape[1] == 1:
                    vals.append(a.reshape(-1))
                else:
                    # keep last axis if multi-channel, but make it contiguous
                    vals.append(a.astype(np.float32))
            if len(vals) != len(streamlines_xyz):
                raise ValueError(
                    f"data_per_point['{name}'] length {len(vals)} does not match "
                    f"number of streamlines {len(streamlines_xyz)}"
                )
            per_point[name] = vals

    if not include_per_streamline:
        return streamlines_xyz, per_point

    per_streamline: dict[str, np.ndarray] = {}
    dps = getattr(tg, "data_per_streamline", None)
    if dps:
        for name, values in dps.items():
            values = np.asarray(values)
            if len(values) != len(streamlines_xyz):
                raise ValueError(
                    f"data_per_streamline['{name}'] length {len(values)} does not "
                    f"match number of streamlines {len(streamlines_xyz)}"
                )
            per_streamline[name] = values

    return streamlines_xyz, per_point, per_streamline


def ha_to_degrees_per_streamline(ha_list: list[np.ndarray]) -> list[np.ndarray]:
    """
    Convert HA values that might be byte-scaled (0..255) to degrees (-90..90).
    Leaves values unchanged if they already look like degrees.
    """
    out: list[np.ndarray] = []
    for ha in ha_list:
        ha = np.asarray(ha)
        if ha.size > 0 and np.nanmax(ha) > 1.5:  # likely 0..255
            ha_deg = (ha.astype(np.float32) / 255.0) * 180.0 - 90.0
        else:
            ha_deg = ha.astype(np.float32)
        out.append(ha_deg)
    return out


def normalize_attrs_to_degrees(attrs: dict | None) -> dict[str, list[np.ndarray]]:
    """
    Normalize HA, IA, AZ, EL fields to degrees if stored as 0–255.
    If already in degrees or unit vectors, returns unchanged except cast to float32.

    Input:
      attrs: dict[str, list[np.ndarray]] from TRK (e.g., {"HA":[...], "IA":[...], ...})

    Returns:
      normalized: same keys, each entry = list of np.ndarray (float32) in degrees.
    """
    if not attrs:
        return {}

    angle_ranges = {
        "HA": (-90.0, 90.0),
        "IA": (-90.0, 90.0),
        "AZ": (-180.0, 180.0),
        "EL": (0.0, 90.0),
    }
    normalized = {}
    for key, arr_list in attrs.items():
        name = key.upper()
        source_arrays = [np.asarray(arr) for arr in arr_list]
        nonempty = [arr.reshape(-1) for arr in source_arrays if arr.size]
        byte_scaled = any(arr.dtype == np.uint8 for arr in source_arrays)
        if nonempty and name in angle_ranges and not byte_scaled:
            field_min = min(float(np.nanmin(arr)) for arr in nonempty)
            field_max = max(float(np.nanmax(arr)) for arr in nonempty)
            _, physical_max = angle_ranges[name]
            byte_scaled = (
                field_min >= 0.0 and field_max <= 255.0 and field_max > physical_max
            )

        out_list = []
        for arr in source_arrays:
            a = np.asarray(arr, dtype=np.float32).reshape(-1)
            if byte_scaled and name in angle_ranges:
                angle_min, angle_max = angle_ranges[name]
                a = angle_min + (a / 255.0) * (angle_max - angle_min)
            out_list.append(a.astype(np.float32))

        normalized[name] = out_list

    return normalized


def compute_elevation_angles(streamlines_xyz: list[np.ndarray]) -> list[np.ndarray]:
    """
    Compute per-vertex elevation angle from streamline geometry:
    elevation = arcsin(z-component of unit tangent) in degrees.
    The last vertex copies the previous value to keep lengths aligned.
    """
    all_angles: list[np.ndarray] = []
    for pts in streamlines_xyz:
        pts = np.asarray(pts, dtype=np.float32)
        n = len(pts)
        if n < 2:
            all_angles.append(np.zeros((n,), dtype=np.float32))
            continue
        vecs = np.diff(pts, axis=0)
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        normalized = np.divide(vecs, norms, out=np.zeros_like(vecs), where=norms != 0)
        elev = np.arcsin(np.clip(normalized[:, 2], -1.0, 1.0)) * (180.0 / np.pi)
        elev = np.concatenate([elev, [elev[-1]]]).astype(np.float32)
        all_angles.append(elev)
    return all_angles


def reduce_per_edge(
    values_per_point: list[np.ndarray], how: str = "mean"
) -> np.ndarray:
    """
    Reduce per-point values along each streamline to a single scalar per edge.
    """
    if how == "mean":
        return np.array(
            [float(np.nanmean(v)) if v.size else np.nan for v in values_per_point],
            dtype=float,
        )
    if how == "median":
        return np.array(
            [float(np.nanmedian(v)) if v.size else np.nan for v in values_per_point],
            dtype=float,
        )
    raise ValueError("how must be 'mean' or 'median'")
