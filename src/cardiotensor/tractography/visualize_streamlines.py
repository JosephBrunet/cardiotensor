import numpy as np
import fury
import random


def compute_elevation_angles(streamlines):
    """Compute per-vertex elevation angle."""
    all_angles = []
    for pts in streamlines:
        if len(pts) < 2:
            all_angles.append(np.zeros((len(pts),), dtype=np.float32))
            continue
        vecs = np.diff(pts, axis=0)
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        normalized = np.divide(vecs, norms, where=norms != 0)
        z_components = normalized[:, 2]
        elev = np.arcsin(z_components) * 180.0 / np.pi
        elev = np.concatenate([elev, [elev[-1]]])  # repeat last value to match point count
        all_angles.append(elev.astype(np.float32))
    return all_angles


def downsample_streamline(streamline, factor=2):
    return streamline if len(streamline) < 3 else streamline[::factor]


def show_streamlines(
    streamlines_xyz: list[np.ndarray],
    color_values: list[np.ndarray],
    mode: str = "tube",
    line_width: int = 4,
    interactive: bool = True,
    screenshot_path: str | None = None,
    window_size: tuple[int, int] = (800, 800),
    downsample_factor: int = 2,
    max_streamlines: int | None = None,
    filter_min_len: int | None = None,
    subsample_factor: int = 1,
):
    print(f"Initial number of streamlines: {len(streamlines_xyz)}")
    if filter_min_len:
        print(f"Filtering out streamlines shorter than {filter_min_len} points")
    if downsample_factor > 1:
        print(f"Downsampling each streamline by factor {downsample_factor}")
    if subsample_factor > 1:
        print(f"Subsampling: keeping 1 in every {subsample_factor} streamlines")
    if max_streamlines:
        print(f"Limiting to max {max_streamlines} streamlines")

    # Downsample and filter
    downsampled_streamlines = []
    downsampled_colors = []
    idx = 0
    for sl in streamlines_xyz:
        color_slice = color_values[idx:idx + len(sl)]
        ds_sl = downsample_streamline(sl, downsample_factor)
        ds_cl = downsample_streamline(color_slice, downsample_factor)

        if filter_min_len is None or len(ds_sl) >= filter_min_len:
            downsampled_streamlines.append(ds_sl)
            downsampled_colors.append(ds_cl)

        idx += len(sl)

    streamlines_xyz = downsampled_streamlines
    color_values = downsampled_colors

    if not streamlines_xyz:
        raise ValueError("❌ No streamlines left after downsampling and filtering.")

    # Subsample
    if subsample_factor > 1:
        total = len(streamlines_xyz)
        selected_idx = sorted(random.sample(range(total), total // subsample_factor))
        streamlines_xyz = [streamlines_xyz[i] for i in selected_idx]
        color_values = [color_values[i] for i in selected_idx]

    # Cap max
    if max_streamlines is not None and len(streamlines_xyz) > max_streamlines:
        selected_idx = sorted(random.sample(range(len(streamlines_xyz)), max_streamlines))
        streamlines_xyz = [streamlines_xyz[i] for i in selected_idx]
        color_values = [color_values[i] for i in selected_idx]

    print(f"Final number of streamlines to render: {len(streamlines_xyz)}")

    flat_colors = np.concatenate(color_values)
    print(f"Coloring mode: min={flat_colors.min():.2f}, max={flat_colors.max():.2f}")
    print(f"Rendering mode: {mode}")

    lut = fury.actor.colormap_lookup_table(
        scale_range=(flat_colors.min(), flat_colors.max()),
        hue_range=(0.7, 0.0),
        saturation_range=(0.5, 1.0),
    )

    scene = fury.window.Scene()

    if mode == "tube":
        actor = fury.actor.streamtube(streamlines_xyz, flat_colors, linewidth=line_width, opacity=1, lookup_colormap=lut)
    elif mode == "fake_tube":
        actor = fury.actor.line(streamlines_xyz, flat_colors, linewidth=line_width, fake_tube=True, depth_cue=True)
    elif mode == "line":
        actor = fury.actor.line(streamlines_xyz, flat_colors, linewidth=line_width, fake_tube=False, depth_cue=False, lookup_colormap=lut)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    scene.add(actor)
    scene.add(fury.actor.scalar_bar(lut))
    scene.reset_camera()

    if interactive:
        print("🕹️ Opening interactive window...")
        fury.window.show(scene, size=window_size, reset_camera=False)
    else:
        if not screenshot_path:
            raise ValueError("Must specify screenshot_path when interactive=False.")
        print(f"📸 Saving screenshot to: {screenshot_path}")
        fury.window.record(scene=scene, out_path=screenshot_path, size=window_size)
