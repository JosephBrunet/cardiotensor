import numpy as np
from alive_progress import alive_bar
from pathlib import Path
import math

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.downsampling import downsample_vector_volume, downsample_volume



def trilinear_interpolate_vector(
    vector_field: np.ndarray, pt: tuple[float, float, float]
) -> np.ndarray:
    """
    Given a fractional (z,y,x), returns the trilinearly‐interpolated 3‐vector
    from `vector_field` (shape = (3, Z, Y, X)). Clamps to nearest voxel if out‐of‐bounds.
    """
    zf, yf, xf = pt
    _, Z, Y, X = vector_field.shape

    # Clamp floor and ceil to valid ranges
    z0 = max(min(int(np.floor(zf)), Z - 1), 0)
    z1 = min(z0 + 1, Z - 1)
    y0 = max(min(int(np.floor(yf)), Y - 1), 0)
    y1 = min(y0 + 1, Y - 1)
    x0 = max(min(int(np.floor(xf)), X - 1), 0)
    x1 = min(x0 + 1, X - 1)

    dz = zf - z0
    dy = yf - y0
    dx = xf - x0

    # 8 corner vectors
    c000 = vector_field[:, z0, y0, x0]
    c001 = vector_field[:, z0, y0, x1]
    c010 = vector_field[:, z0, y1, x0]
    c011 = vector_field[:, z0, y1, x1]
    c100 = vector_field[:, z1, y0, x0]
    c101 = vector_field[:, z1, y0, x1]
    c110 = vector_field[:, z1, y1, x0]
    c111 = vector_field[:, z1, y1, x1]

    # Interpolate along X
    c00 = c000 * (1 - dx) + c001 * dx
    c01 = c010 * (1 - dx) + c011 * dx
    c10 = c100 * (1 - dx) + c101 * dx
    c11 = c110 * (1 - dx) + c111 * dx

    # Interpolate along Y
    c0 = c00 * (1 - dy) + c01 * dy
    c1 = c10 * (1 - dy) + c11 * dy

    # Interpolate along Z
    c = c0 * (1 - dz) + c1 * dz
    return c  # shape (3,)


def trilinear_interpolate_scalar(
    volume: np.ndarray, pt: tuple[float, float, float]
) -> float:
    """
    Trilinearly interpolate a scalar volume at fractional point (z, y, x).
    Clamps to valid range.
    """
    zf, yf, xf = pt
    Z, Y, X = volume.shape

    z0 = int(np.floor(zf))
    z1 = min(z0 + 1, Z - 1)
    y0 = int(np.floor(yf))
    y1 = min(y0 + 1, Y - 1)
    x0 = int(np.floor(xf))
    x1 = min(x0 + 1, X - 1)

    dz = zf - z0
    dy = yf - y0
    dx = xf - x0

    c000 = volume[z0, y0, x0]
    c001 = volume[z0, y0, x1]
    c010 = volume[z0, y1, x0]
    c011 = volume[z0, y1, x1]
    c100 = volume[z1, y0, x0]
    c101 = volume[z1, y0, x1]
    c110 = volume[z1, y1, x0]
    c111 = volume[z1, y1, x1]

    c00 = c000 * (1 - dx) + c001 * dx
    c01 = c010 * (1 - dx) + c011 * dx
    c10 = c100 * (1 - dx) + c101 * dx
    c11 = c110 * (1 - dx) + c111 * dx

    c0 = c00 * (1 - dy) + c01 * dy
    c1 = c10 * (1 - dy) + c11 * dy

    c = c0 * (1 - dz) + c1 * dz
    return float(c)


def trace_streamline(
    start_pt: tuple[float, float, float],
    vector_field: np.ndarray,
    fa_volume: np.ndarray | None = None,
    fa_threshold: float = 0.1,
    step_length: float = 0.5,
    max_steps: int | None = 1000,
    angle_threshold: float = 60.0,
    eps: float = 1e-10,
    direction: int = 1,
) -> list[tuple[float, float, float]]:
    """
    Trace one streamline from `start_pt` (z,y,x) in the continuous vector_field.
    - Interpolate & normalize each sub‐step
    - Move forward by `step_length` voxels each step using RK4
    - `direction` = +1 (default) or -1 to reverse integration direction
    """
    Z, Y, X = vector_field.shape[1:]
    coords: list[tuple[float, float, float]] = [
        (float(start_pt[0]), float(start_pt[1]), float(start_pt[2]))
    ]
    current_pt = np.array(start_pt, dtype=np.float64)
    prev_dir: np.ndarray | None = None  # previous unit vector

    def interp_unit(pt: np.ndarray) -> np.ndarray | None:
        """Return a normalized direction vector at fractional pt, or None if invalid."""
        vec = trilinear_interpolate_vector(vector_field, (pt[0], pt[1], pt[2]))
        if np.isnan(vec).any():
            return None
        norm = np.linalg.norm(vec)
        if norm < eps:
            return None
        return np.array([vec[2], vec[1], vec[0]]) / norm * direction  # flip to (z,y,x)

    step_count = 0
    while max_steps is None or step_count < max_steps:
        step_count += 1

        if fa_volume is not None:
            fa_value = trilinear_interpolate_scalar(fa_volume, tuple(current_pt))
            if fa_value < fa_threshold:
                break

        k1 = interp_unit(current_pt)
        if k1 is None:
            break
        if prev_dir is not None:
            angle = np.degrees(np.arccos(np.clip(np.dot(prev_dir, k1), -1.0, 1.0)))
            if angle > angle_threshold:
                break

        mid1 = current_pt + 0.5 * step_length * k1
        k2 = interp_unit(mid1)
        if k2 is None:
            break

        mid2 = current_pt + 0.5 * step_length * k2
        k3 = interp_unit(mid2)
        if k3 is None:
            break

        end_pt = current_pt + step_length * k3
        k4 = interp_unit(end_pt)
        if k4 is None:
            break

        angle4 = np.degrees(np.arccos(np.clip(np.dot(k1, k4), -1.0, 1.0)))
        if angle4 > angle_threshold:
            break

        increment = (step_length / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        next_pt = current_pt + increment
        next_dir = k1

        zn, yn, xn = next_pt
        if not (0 <= zn < Z and 0 <= yn < Y and 0 <= xn < X):
            break

        coords.append((float(zn), float(yn), float(xn)))
        current_pt = next_pt
        prev_dir = next_dir

    return coords



def generate_streamlines_from_vector_field(
    vector_field: np.ndarray,
    seed_points: np.ndarray,
    fa_volume: np.ndarray | None = None,
    fa_threshold: float = 0.1,
    step_length: float = 0.5,
    max_steps: int | None = None,
    angle_threshold: float = 60.0,
    min_length_pts: int = 10,
    bidirectional: bool = True,
) -> list[list[tuple[float, float, float]]]:
    """
    Given a 3D vector_field (shape = (3, Z, Y, X)) and a set of integer‐seed voxels,
    returns a list of streamlines (each streamline = a list of float (z,y,x) points),
    filtered so that only those longer than `min_length_pts` are kept.
    Optionally trace in both directions from each seed point.
    """
    all_streamlines: list[list[tuple[float, float, float]]] = []
    total_seeds = len(seed_points)

    with alive_bar(total_seeds, title="Tracing Streamlines") as bar:
        for zi, yi, xi in seed_points:
            start = (float(zi), float(yi), float(xi))

            # Forward tracing
            forward_pts = trace_streamline(
                start_pt=start,
                vector_field=vector_field,
                fa_volume=fa_volume,
                fa_threshold=fa_threshold,
                step_length=step_length,
                max_steps=max_steps,
                angle_threshold=angle_threshold,
                direction=1,
            )

            # Backward tracing if enabled
            if bidirectional:
                backward_pts = trace_streamline(
                    start_pt=start,
                    vector_field=vector_field,
                    fa_volume=fa_volume,
                    fa_threshold=fa_threshold,
                    step_length=step_length,
                    max_steps=max_steps,
                    angle_threshold=angle_threshold,
                    direction=-1,
                )
                # Remove duplicate seed point and reverse
                backward_pts = backward_pts[::-1][:-1] if len(backward_pts) > 1 else []
                full_streamline = backward_pts + forward_pts
            else:
                full_streamline = forward_pts

            if len(full_streamline) >= min_length_pts:
                all_streamlines.append(full_streamline)

            bar()

    return all_streamlines




def generate_streamlines_from_params(
    vector_field_dir: str | Path,
    output_dir: str | Path,
    fa_dir: str | Path,
    ha_dir: str | Path,
    mask_path: str | Path | None = None,
    start_xyz: tuple[int, int, int] = (0, 0, 0),
    end_xyz: tuple[int | None, int | None, int | None] = (None, None, None),
    bin_factor: int = 1,
    num_seeds: int = 20000,
    fa_threshold: float = 0.1,
    step_length: float = 0.5,
    max_steps: int | None = None,
    angle_threshold: float = 60.0,
    min_length_pts: int = 10,
    bidirectional: bool = True,
) -> None:
    """
    Generate streamlines from a vector field and save to NPZ.

    Parameters
    ----------
    vector_field_dir : str or Path
        Directory containing eigenvector volumes.
    output_dir : str or Path
        Directory where output NPZ will be saved.
    fa_dir : str or Path
        Directory containing FA volumes.
    ha_dir : str or Path
        Directory containing HA volumes.
    mask_path : str or Path, optional
        Path to mask volume to filter vectors.
    start_xyz, end_xyz : tuple[int]
        Cropping bounds in Z, Y, X (None = full dimension).
    bin_factor : int
        Spatial downsampling factor.
    num_seeds : int
        Number of random seed points for streamline generation.
    fa_threshold : float
        Minimum FA value to continue tracing.
    step_length : float
        Streamline step length in voxels.
    max_steps : int, optional
        Maximum number of steps per streamline.
    angle_threshold : float
        Maximum angle (degrees) allowed between steps.
    min_length_pts : int
        Minimum number of points for a valid streamline.
    bidirectional : bool
        Trace in both directions from each seed.
    """
    vector_field_dir = Path(vector_field_dir)
    fa_dir = Path(fa_dir)
    ha_dir = Path(ha_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    start_z, start_y, start_x = start_xyz
    end_z, end_y, end_x = end_xyz

    # --- Load vector field for shape ---
    vec_reader = DataReader(vector_field_dir)
    full_shape = vec_reader.shape  # (3, Z, Y, X)
    if end_z is None:
        end_z = full_shape[1]
    if end_y is None:
        end_y = full_shape[2]
    if end_x is None:
        end_x = full_shape[3]

    # --- Handle binning ---
    if bin_factor > 1:
        downsample_vector_volume(vector_field_dir, bin_factor, output_dir)
        vec_load_dir = output_dir / f"bin{bin_factor}" / vector_field_dir.name

        downsample_volume(fa_dir, bin_factor, output_dir, subfolder="FA", out_ext="tif")
        fa_load_dir = output_dir / f"bin{bin_factor}" / "FA"

        downsample_volume(ha_dir, bin_factor, output_dir, subfolder="HA", out_ext="tif")
        ha_load_dir = output_dir / f"bin{bin_factor}" / "HA"

        start_z_binned = start_z // bin_factor
        end_z_binned = math.ceil(end_z / bin_factor)
        start_y_binned = start_y // bin_factor
        end_y_binned = math.ceil(end_y / bin_factor)
        start_x_binned = start_x // bin_factor
        end_x_binned = math.ceil(end_x / bin_factor)
    else:
        vec_load_dir = vector_field_dir
        fa_load_dir = fa_dir
        ha_load_dir = ha_dir
        start_z_binned, end_z_binned = start_z, end_z
        start_y_binned, end_y_binned = start_y, end_y
        start_x_binned, end_x_binned = start_x, end_x

    # --- Load vector field ---
    print("📥 Loading vector field...")
    vector_field = DataReader(vec_load_dir).load_volume(
        start_index=start_z_binned, end_index=end_z_binned
    )[:, :, start_y_binned:end_y_binned, start_x_binned:end_x_binned]

    # Ensure channel order (3, Z, Y, X)
    if vector_field.ndim == 4 and vector_field.shape[-1] == 3:
        vector_field = np.moveaxis(vector_field, -1, 0)

    # Flip vectors for consistency
    neg_mask = vector_field[0] < 0
    vector_field[:, neg_mask] *= -1

    # --- Mask if provided ---
    if mask_path:
        print(f"🩹 Applying mask: {mask_path}")
        mask_reader = DataReader(mask_path)
        
        mask = mask_reader.load_volume(
            start_index=start_z_binned,
            end_index=end_z_binned,
            unbinned_shape=vec_reader.shape[1:],  # (Z, Y, X)
        )

        mask = mask[:, start_y_binned:end_y_binned, start_x_binned:end_x_binned]
        mask = (mask > 0).astype(np.uint8)
        
        vector_field[:, mask == 0] = np.nan

    # --- Load FA and create seed points ---
    print("📥 Loading FA volume...")
    fa_volume = DataReader(fa_load_dir).load_volume(
        start_index=start_z_binned, end_index=end_z_binned
    )
    fa_volume = fa_volume[:, start_y_binned:end_y_binned, start_x_binned:end_x_binned]

    seed_mask = fa_volume > (fa_threshold * 255)
    valid_indices = np.argwhere(seed_mask)
    if len(valid_indices) < num_seeds:
        print("⚠️ Not enough valid seed points; using all available voxels above threshold.")
        chosen_indices = valid_indices
    else:
        chosen_indices = valid_indices[
            np.random.choice(valid_indices.shape[0], num_seeds, replace=False)
        ]

    # --- Generate streamlines ---
    streamlines = generate_streamlines_from_vector_field(
        vector_field=vector_field,
        seed_points=chosen_indices,
        fa_volume=fa_volume,
        fa_threshold=fa_threshold,
        step_length=step_length,
        max_steps=max_steps,
        angle_threshold=angle_threshold,
        min_length_pts=min_length_pts,
        bidirectional=bidirectional,
    )

    # --- Load HA for sampling ---
    print("📥 Loading HA volume for sampling...")
    ha_volume = DataReader(ha_load_dir).load_volume(
        start_index=start_z_binned, end_index=end_z_binned
    )
    ha_volume = ha_volume[:, start_y_binned:end_y_binned, start_x_binned:end_x_binned]

    def sample_ha_along(streamline):
        values = []
        Z, Y, X = ha_volume.shape
        for z, y, x in streamline:
            zi, yi, xi = map(int, [round(z), round(y), round(x)])
            zi, yi, xi = max(0, min(zi, Z - 1)), max(0, min(yi, Y - 1)), max(0, min(xi, X - 1))
            values.append(float(ha_volume[zi, yi, xi]))
        return values

    all_ha = [val for sl in streamlines for val in sample_ha_along(sl)]

    # --- Save output ---
    out_path = output_dir / "streamlines.npz"
    np.savez_compressed(
        out_path,
        streamlines=np.array(streamlines, dtype=object),
        ha_values=np.array(all_ha, dtype=np.float32),
    )
    print(f"✅ Saved {len(streamlines)} streamlines to {out_path}")

