import numpy as np
from alive_progress import alive_bar


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
    z0 = int(np.floor(zf))
    z1 = min(z0 + 1, Z - 1)
    y0 = int(np.floor(yf))
    y1 = min(y0 + 1, Y - 1)
    x0 = int(np.floor(xf))
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
) -> list[tuple[float, float, float]]:
    """
    Trace one streamline from `start_pt` (z,y,x) in the continuous vector_field.
    - Interpolate & normalize each sub‐step
    - Move forward by `step_length` voxels each step (Euler or RK4)
    - Stop if turning angle > angle_threshold or out of bounds or `vec` too small.
    - If max_steps is None, trace until a stopping condition is hit (no hard limit).
    """
    Z, Y, X = vector_field.shape[1:]
    coords: list[tuple[float, float, float]] = [
        (float(start_pt[0]), float(start_pt[1]), float(start_pt[2]))
    ]
    current_pt = np.array(start_pt, dtype=np.float64)
    prev_dir: np.ndarray | None = None  # previous “unit vector”

    def interp_unit(pt: np.ndarray) -> np.ndarray | None:
        """Return a normalized direction vector at fractional pt, or None if invalid."""
        vec = trilinear_interpolate_vector(vector_field, (pt[0], pt[1], pt[2]))
        if np.isnan(vec).any():
            return None
        norm = np.linalg.norm(vec)
        if norm < eps:
            return None
        return np.array([vec[2], vec[1], vec[0]]) / norm  # flip to (z,y,x) order

    step_count = 0
    while max_steps is None or step_count < max_steps:
        step_count += 1

        if fa_volume is not None:
            fa_value = trilinear_interpolate_scalar(fa_volume, tuple(current_pt))
            if fa_value < fa_threshold:
                break

        # ------
        # Runge-Kutta 4th order integration

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

        # ------

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
) -> list[list[tuple[float, float, float]]]:
    """
    Given a 3D vector_field (shape = (3, Z, Y, X)) and a set of integer‐seed voxels,
    returns a list of streamlines (each streamline = a list of float (z,y,x) points),
    filtered so that only those longer than `min_length_pts` are kept.
    Displays a progress bar during processing.
    """
    all_streamlines: list[list[tuple[float, float, float]]] = []
    total_seeds = len(seed_points)

    with alive_bar(total_seeds, title="Tracing Streamlines") as bar:
        for zi, yi, xi in seed_points:
            start = (float(zi), float(yi), float(xi))
            pts = trace_streamline(
                start_pt=start,
                vector_field=vector_field,
                fa_volume=fa_volume,
                fa_threshold=fa_threshold,
                step_length=step_length,
                max_steps=max_steps,
                angle_threshold=angle_threshold,
            )
            if len(pts) >= min_length_pts:
                all_streamlines.append(pts)
            bar()

    return all_streamlines
