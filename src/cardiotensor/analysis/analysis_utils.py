# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 11:37:23 2025

@author: A Cook
"""


import os

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt


import cv2

from matplotlib.patches import Wedge, Circle
from matplotlib.patches import Patch

from scipy.stats import circmean, circstd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy import ndimage as ndi
from scipy.ndimage import map_coordinates

from skimage import measure
from skimage.morphology import remove_small_objects, remove_small_holes, closing, opening, disk, erosion
from skimage.io import imread
from skimage import exposure

from astropy import stats

from sklearn.metrics import r2_score

import pyvista as pv


def get_septum_point(lv_mask):
    """
    Manually select a septum point from an LV mask.

    Parameters
    ----------
    lv_mask : str
        Directory containing LV mask slices.

    Returns
    -------
    np.ndarray
        Septum point as (z, y, x).
    """

    # Load volume
    files = sorted(
        f for f in os.listdir(lv_mask)
        if f.lower().endswith((".tif", ".tiff", ".jp2"))
    )

    if not files:
        raise ValueError("Directory contains no TIFF or JP2 images")

    slices = []
    for f in files:
        img = Image.open(os.path.join(lv_mask, f))
        slices.append(np.array(img))

    vol = np.stack(slices, axis=0)

    # Find slice with largest LV area
    zsum = vol.sum(axis=(1, 2))
    mid = int(np.argmax(zsum))

    # Display slice and collect click
    clicked = []

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(vol[mid], cmap="gray")
    ax.set_title(
        f"Click once on the septum (slice {mid})"
    )

    def onclick(event):
        if event.xdata is None or event.ydata is None:
            return

        clicked.append(
            np.array([
                mid,
                int(round(event.ydata)),
                int(round(event.xdata))
            ])
        )

        plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", onclick)
    plt.show()

    if not clicked:
        raise RuntimeError("No click recorded")

    return clicked[0]




def get_RV_mask(MV, AX, heart_mask, lv_mask):

    """Return the right-ventricle mask (heart_mask minus lv_mask),
    clipped by the plane perpendicular to the long-axis vector that
    intersects the mitral valve point `MV`.

    Parameters
    ----------
    MV, AX : array-like
        Mitral valve and apex points (z,y,x) in same coordinate system as masks.
    heart_mask : ndarray or path-like
        Whole-heart mask volume or path to folder (will be resampled to lv_mask shape).
    lv_mask : ndarray
        Left-ventricle mask volume.

    Returns
    -------
    rv_mask : ndarray
        Binary mask of the RV clipped at the MV plane (voxels on the "apex"
        side of the MV are kept).
    """

    axis_vec = np.asarray(MV, dtype=float) - np.asarray(AX, dtype=float)
    axis_vec = axis_vec / np.linalg.norm(axis_vec)

    # Ensure masks have same shape
    if heart_mask.shape != lv_mask.shape:
        heart_mask = resample(heart_mask, target_shape=lv_mask.shape, order=0, TwoD=False)

    # RV = heart - LV (assumes binary masks)
    rv_mask = (heart_mask > 0).astype(np.uint8) - (lv_mask > 0).astype(np.uint8)
    rv_mask = np.clip(rv_mask, 0, 1).astype(np.uint8)

    # Clip RV mask at plane through MV with normal axis_vec
    # Build coordinates of all voxels (z,y,x)
    nz, ny, nx = rv_mask.shape
    zz, yy, xx = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing='ij')
    coords = np.stack([zz, yy, xx], axis=-1).reshape(-1, 3)

    # Compute signed distance along axis_vec from MV
    MV = np.asarray(MV, dtype=float)
    rel = coords - MV
    signed = rel @ axis_vec

    # Keep voxels on apex side of MV (signed <= 0)
    keep = signed <= 0

    rv_flat = rv_mask.ravel()
    rv_flat[~keep] = 0
    rv_mask = rv_flat.reshape(rv_mask.shape)

    # ------------------------------------------------------------
    """
    # 3D mesh of mask
    vol_pad = np.pad(rv_mask.astype(np.uint8), 1, mode='constant')
    verts, faces, norms, vals = measure.marching_cubes(vol_pad, 0.5)

    faces = np.hstack([np.full((faces.shape[0],1),3), faces]).astype(np.int64)
    mesh = pv.PolyData(verts, faces)

    # ------------------------------------------------------------
    # Plot
    pl = pv.Plotter()
    pl.add_mesh(mesh, color="lightgray", opacity=0.6)

    pl.add_axes()
    pl.show_grid()
    pl.show()
    """

    return rv_mask



def load_volume(img_dir):
    #Load TIFF/JP2 stack into 3D numpy array
    files = sorted([f for f in os.listdir(img_dir) if f.lower().endswith(('.tif','.tiff','.jp2'))])
    if not files:
        raise ValueError(f"No images found in {img_dir}")
    slices = []
    for f in files:
        img = Image.open(os.path.join(img_dir, f))
        slices.append(np.array(img))
    vol = np.stack(slices, axis=0)
    return vol


def create_aha_mask(
    lv_mask_path,
    axis_points,
    septum_point,
    images_path,
    return_radial_map=False,
    additional_downsample=1,
    rv=False,
    rv_mask_path=None, 
    heart_mask_path=None,
    headless=False,
    savepath=None
):
    """
    Create a 17-segment AHA mask from an LV mask.

    Parameters
    ----------
    lv_mask_path : str
        Path to LV mask volume.
    axis_points : list
        [MV_point, Apex_point] from config file.
    septum_point : array-like
        Septum reference point.
    images_path : str
        Path to image stack (used to determine downsampling).
    return_radial_map : bool, optional
        If True, also returns the radial subdivision map.

    Returns
    -------
    segment_map : ndarray
        17-segment AHA mask.
    segment_radial_map : ndarray, optional
        Segment map with radial quartiles encoded as seg*10 + radial_group.
    """

    # --------------------------------------------------
    # Determine downsampling
    # --------------------------------------------------
    files_lv = sorted(
        f for f in os.listdir(lv_mask_path)
        if f.lower().endswith((".tif", ".tiff", ".jp2"))
    )
    files_img = sorted(
        f for f in os.listdir(images_path)
        if f.lower().endswith((".tif", ".tiff", ".jp2"))
    )

    # Avoid zero division or zero scaling when image and mask slice counts differ greatly.
    # Compute ratio and ensure integer scaling factor >= 1
    mask_ratio = float(len(files_img)) / max(1, float(len(files_lv)))
    mask_factor = max(1, int(round(mask_ratio)))

    # Compute downsample as aha_factor / mask_factor, at least 1
    try:
        downsample = int(max(1, round(float(additional_downsample) * float(mask_factor))))
    except Exception:
        downsample = mask_factor

    # --------------------------------------------------
    # Load and resample LV mask
    # --------------------------------------------------
    lv_mask_vol = load_volume(lv_mask_path)
    lv_mask_vol_d = resample(
        lv_mask_vol,
        zoom=1 / additional_downsample,
        order=0,
    )

    # --------------------------------------------------
    # Long-axis definition
    # --------------------------------------------------
    MV_d = np.array(axis_points[0])[::-1] / downsample
    Apex_d = np.array(axis_points[1])[::-1] / downsample

    axis_vec = MV_d - Apex_d
    axis_vec /= np.linalg.norm(axis_vec)

    coords = np.array(np.nonzero(lv_mask_vol_d)).T
    centroid = (MV_d + Apex_d) / 2

    septum_point = np.asarray(septum_point) / additional_downsample

    # --------------------------------------------------
    # Create RV mask if requested
    # --------------------------------------------------
    if rv== True:
        if rv_mask_path is None:
            if heart_mask_path is None:
                raise ValueError("If rv=True, either rv_mask or heart_mask must be provided.")
            else:
                heart_mask = load_volume(heart_mask_path)
                rv_mask = get_RV_mask(MV_d, Apex_d, heart_mask, lv_mask_vol_d)
        else:
            rv_mask = load_volume(rv_mask_path)
            # Ensure rv_mask has same shape as lv_mask_vol_
            if rv_mask.shape != lv_mask_vol_d.shape:
                rv_mask = resample(rv_mask, target_shape=lv_mask_vol_d.shape, order=0, TwoD=False, is_mask=True)


    # --------------------------------------------------
    # Long-axis projections
    # --------------------------------------------------
    projections = coords @ axis_vec

    # Discard extreme percentiles to guard againts noise and outliers in the mask
    pmin = np.percentile(projections, 2)
    pmax = np.percentile(projections, 98)
    length = pmax - pmin

    basal_thresh = pmin + 2 * length / 3.0
    mid_thresh = pmin + length / 3.0

    # --------------------------------------------------
    # Septum reference angle
    # --------------------------------------------------
    septum_vec = septum_point - centroid

    pts2d, u, v = project_to_plane(
        coords,
        centroid,
        axis_vec,
    )

    septum_point2d = np.array([
        np.dot(septum_vec, u),
        np.dot(septum_vec, v)
    ])

    ref_angle = np.arctan2(
        septum_point2d[1],
        septum_point2d[0]
    )

    angles = np.arctan2(pts2d[:, 1], pts2d[:, 0])

    angles_rel = (angles - ref_angle) % (2 * np.pi)
    angles_rel = (angles_rel + 2 * np.pi / 3) % (2 * np.pi)

    # --------------------------------------------------
    # Create segment map
    # --------------------------------------------------
    segment_map = np.zeros(
        lv_mask_vol_d.shape,
        dtype=np.uint16
    )

    def assign_segment(mask_idx, seg):
        segment_map[tuple(coords[mask_idx].T)] = seg

    basal_idx = np.where(projections >= basal_thresh)[0]

    mid_idx = np.where(
        (projections >= mid_thresh) &
        (projections < basal_thresh)
    )[0]

    apical_idx = np.where(
        (projections < mid_thresh) &
        (projections >= pmin + 0.02 * length)
    )[0]

    apex_idx = np.where(
        projections < (pmin + 0.02 * length)
    )[0]

    def sectors_from_indices(indices, n_sectors, start_seg):
        sector_width = 2 * np.pi / n_sectors

        for s in range(n_sectors):
            a0 = s * sector_width
            a1 = (s + 1) * sector_width

            mask_idx = indices[
                (angles_rel[indices] >= a0) &
                (angles_rel[indices] < a1)
            ]

            assign_segment(mask_idx, start_seg + s)

    # Basal: 1–6
    sectors_from_indices(basal_idx, 6, 1)

    # Mid: 7–12
    sectors_from_indices(mid_idx, 6, 7)

    # Apical: 13–16
    sector_width = 2 * np.pi / 4

    for s in range(4):
        a0 = s * sector_width
        a1 = (s + 1) * sector_width

        mask_idx = apical_idx[
            (angles_rel[apical_idx] >= a0) &
            (angles_rel[apical_idx] < a1)
        ]

        assign_segment(mask_idx, 13 + s)

    # Apex: 17
    if apex_idx.size > 0:
        assign_segment(apex_idx, 17)

    
    # ==================================================
    # RV segmentation (segments 18-26)
    # ==================================================

    if rv:

        rv_coords = np.array(np.nonzero(rv_mask)).T

        if len(rv_coords) > 0:

            # ------------------------------------------
            # Long-axis coordinate
            # ------------------------------------------

            rv_proj = rv_coords @ axis_vec

            rv_pmin = np.percentile(rv_proj, 2)
            rv_pmax = np.percentile(rv_proj, 98)

            mid_proj = 0.5 * (rv_pmin + rv_pmax)

            rv_length = rv_pmax - rv_pmin

            rv_basal_thresh = rv_pmin + 2 * rv_length / 3
            rv_mid_thresh = rv_pmin + rv_length / 3

            # ------------------------------------------
            # Short-axis anatomical basis
            # ------------------------------------------

            long_axis = axis_vec 

            lv_vec = septum_vec / np.linalg.norm(septum_vec)

            # project LV direction into short-axis plane
            e0 = (
                lv_vec
                - np.dot(lv_vec, long_axis) * long_axis
            )

            e0 /= np.linalg.norm(e0)

            # perpendicular direction in short-axis plane
            e1 = np.cross(long_axis, e0)
            e1 /= np.linalg.norm(e1)

            # ------------------------------------------
            # RV angular coordinate
            # ------------------------------------------

            #rv_centroid = rv_coords.mean(axis=0)
            
            rv_centroid = (
                septum_point
                + (
                    mid_proj
                    - np.dot(septum_point, long_axis)
                ) * long_axis
            )

            rv_rel = rv_coords - rv_centroid

            x = rv_rel @ e1
            y = rv_rel @ e0

            rv_angles = np.arctan2(y, x)
            rv_angles = rv_angles % (2 * np.pi)

            # ------------------------------------------
            # Longitudinal regions
            # ------------------------------------------

            rv_basal_idx = np.where(
                rv_proj >= rv_basal_thresh
            )[0]

            rv_mid_idx = np.where(
                (rv_proj >= rv_mid_thresh) &
                (rv_proj < rv_basal_thresh)
            )[0]

            rv_apical_idx = np.where(
                rv_proj < rv_mid_thresh
            )[0]

            # ------------------------------------------
            # Sector assignment
            # ------------------------------------------

            def assign_rv_region(
                indices,
                start_seg
            ):
                """
                Six 60° sectors are first defined.

                sector 0 -> seg 1
                sector 1 -> seg 2
                sector 2 -> seg 3

                sector 3 + half sector 4 -> seg 3

                half sector 4 + sector 5 -> seg 1
                """

                if len(indices) == 0:
                    return

                theta = rv_angles[indices]

                seg_local = np.full(
                    theta.shape,
                    -1,
                    dtype=np.int32
                )

                pi3 = np.pi / 3


                # sector 0
                seg_local[
                    (theta >= 0) &
                    (theta < pi3)
                ] = 2

                # sector 1
                seg_local[
                    (theta >= pi3) &
                    (theta < 2*pi3)
                ] = 1

                # sector 2
                seg_local[
                    (theta >= 2*pi3) &
                    (theta < 3*pi3)
                ] = 0

                # sector 3 + first half sector 4
                seg_local[
                    (theta >= 3*pi3) &
                    (theta < 4.5*pi3)
                ] = 0

                # second half sector 4 + sector 5
                seg_local[
                    (theta >= 4.5*pi3) &
                    (theta < 6*pi3)
                ] = 2
                
                for local_seg in range(3):

                    mask_idx = indices[
                        seg_local == local_seg
                    ]

                    if len(mask_idx) == 0:
                        continue

                    segment_map[
                        tuple(rv_coords[mask_idx].T)
                    ] = start_seg + local_seg

            # basal: 18-20
            assign_rv_region(
                rv_basal_idx,
                18
            )

            # mid: 21-23
            assign_rv_region(
                rv_mid_idx,
                21
            )

            # apical: 24-26
            assign_rv_region(
                rv_apical_idx,
                24
            )
            
    # --------------------------------------------------
    # Radial quartile subdivision
    # --------------------------------------------------
    
    if return_radial_map:
        dists = np.linalg.norm(coords - centroid, axis=1)

        lv_segment_radial_map = np.zeros_like(
            segment_map,
            dtype=np.uint16
        )

        for seg_id in range(1, 18):

            seg_idx = np.where(
                segment_map[tuple(coords.T)] == seg_id
            )[0]

            if len(seg_idx) == 0:
                continue

            seg_dists = dists[seg_idx]

            q25, q50, q75 = np.percentile(
                seg_dists,
                [25, 50, 75]
            )

            for voxel_idx, d in zip(seg_idx, seg_dists):

                if d <= q25:
                    radial_group = 1
                elif d <= q50:
                    radial_group = 2
                elif d <= q75:
                    radial_group = 3
                else:
                    radial_group = 4

                lv_segment_radial_map[
                    tuple(coords[voxel_idx])
                ] = seg_id * 10 + radial_group
        
        if rv:
            
            rv_dists = np.linalg.norm(rv_coords - rv_centroid, axis=1)

            rv_segment_radial_map = np.zeros_like(
                segment_map,
                dtype=np.uint16
            )

            for seg_id in range(18, 27):

                seg_idx = np.where(
                    segment_map[tuple(rv_coords.T)] == seg_id
                )[0]

                if len(seg_idx) == 0:
                    continue

                seg_dists = rv_dists[seg_idx]

                q25, q50, q75 = np.percentile(
                    seg_dists,
                    [25, 50, 75]
                )

                for voxel_idx, d in zip(seg_idx, seg_dists):

                    if d <= q25:
                        radial_group = 1
                    elif d <= q50:
                        radial_group = 2
                    elif d <= q75:
                        radial_group = 3
                    else:
                        radial_group = 4

                    rv_segment_radial_map[
                        tuple(rv_coords[voxel_idx])
                    ] = seg_id * 10 + radial_group

            # Combine LV and RV radial maps
            segment_radial_map = lv_segment_radial_map + rv_segment_radial_map

        else: 
            segment_radial_map = lv_segment_radial_map

    #----------------------------------------------
    # Visualisation 
    #----------------------------------------------

    if rv:
        max_seg = 26
    else:
        max_seg = 17

    def _get_segment_colors(n):
        # Prefer a categorical colormap for up to 20 colors, fall back to HSV for more
        try:
            if n <= 20:
                cmap = plt.get_cmap('tab20', n)
                colors = cmap(np.linspace(0, 1, n, endpoint=False))[:, :3]
            else:
                colors = plt.cm.hsv(np.linspace(0, 1, n, endpoint=False))[:, :3]
        except Exception:
            colors = plt.cm.hsv(np.linspace(0, 1, n, endpoint=False))[:, :3]
        return colors

    segment_colors = _get_segment_colors(max_seg)
    
    n_slices = 20

    # Step 1: project voxels along long axis
    apex_dir = -axis_vec

    # combine LV and RV coordinates if RV exists so projections and
    # slice extraction operate on a single coordinate array
    if rv:
        combined_coords = np.vstack([coords, rv_coords])
        combined_centroid = combined_coords.mean(axis=0)
    else:
        combined_coords = coords
        combined_centroid = centroid

    projections_norm = (combined_coords - combined_centroid) @ apex_dir
    slice1 = np.percentile(projections_norm, 2)
    slicen = np.percentile(projections_norm, 98)
    slice_positions = np.linspace(slice1, slicen, n_slices)
    
    # Step 2: find two orthonormal vectors perpendicular to apex_dir
    v = apex_dir / np.linalg.norm(apex_dir)
    if abs(v[0]) < abs(v[1]):
        u1 = np.array([0, -v[2], v[1]])
    else:
        u1 = np.array([-v[2], 0, v[0]])
    u1 /= np.linalg.norm(u1)
    u2 = np.cross(v, u1)
    u2 /= np.linalg.norm(u2)
    
    # Step 3: extract slices
    slices_rgb = []
    for s in slice_positions:
        
        # select voxels near the slice plane (from combined coords)
        dist = np.abs(projections_norm - s)
        plane_idx = np.where(dist < 2)[0]
        plane_coords = combined_coords[plane_idx]
    
        if len(plane_coords) == 0:
            continue
    
        # project to 2D plane
        pts2d = np.column_stack([
            (plane_coords - combined_centroid) @ u1,
            (plane_coords - combined_centroid) @ u2
        ])
        x = np.round(pts2d[:,0] - pts2d[:,0].min()).astype(int)
        y = np.round(pts2d[:,1] - pts2d[:,1].min()).astype(int)
    
        # create RGB slice (use float colors in 0-1 range)
        shape = (y.max()+1, x.max()+1, 3)
        rgb = np.zeros(shape, dtype=float)
        # ensure integer indexing into the segment map
        seg_idx0 = plane_coords[:, 0].astype(int)
        seg_idx1 = plane_coords[:, 1].astype(int)
        seg_idx2 = plane_coords[:, 2].astype(int)
        seg_ids = segment_map[seg_idx0, seg_idx1, seg_idx2]
        for seg_id in range(1, max_seg + 1):
            mask = seg_ids == seg_id
            if not np.any(mask):
                continue
            # use only RGB channels (ignore alpha) and values are 0-1 from colormap
            rgb[y[mask], x[mask]] = segment_colors[seg_id-1][:3]
        slices_rgb.append(rgb)
    
    # Step 4: plot slices in 2 rows with legend
    rows = 2
    cols = (len(slices_rgb)+1)//2
    fig, axes = plt.subplots(rows, cols, figsize=(3*cols, 3*rows))
    axes = axes.flatten()
    
    for i, rgb in enumerate(slices_rgb):
        axes[i].imshow(rgb)
        axes[i].set_title(f"Slice {i+1}")
        axes[i].axis('off')
    
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    
    # Legend
    # segment_colors values are floats in 0-1 from the colormap; use RGB channels
    patches = [Patch(facecolor=segment_colors[i][:3], label=f"Seg {i+1}") for i in range(max_seg)]
    fig.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.1, 0.95), title="AHA Segments")

    if not headless:
        plt.tight_layout()
        plt.show()
    
    if savepath:
        out_pdf = savepath / f"aha_map_slices.pdf"
        fig.savefig(out_pdf, bbox_inches="tight")

    if not return_radial_map:
        return segment_map, None

    return segment_map, segment_radial_map


def create_aha_mask_17(
    lv_mask_path,
    axis_points,
    septum_point,
    images_path,
    return_radial_map=False,
    additional_downsample=1
):
    """
    Create a 17-segment AHA mask from an LV mask.

    Parameters
    ----------
    lv_mask_path : str
        Path to LV mask volume.
    axis_points : list
        [MV_point, Apex_point] from config file.
    septum_point : array-like
        Septum reference point.
    images_path : str
        Path to image stack (used to determine downsampling).
    return_radial_map : bool, optional
        If True, also returns the radial subdivision map.

    Returns
    -------
    segment_map : ndarray
        17-segment AHA mask.
    segment_radial_map : ndarray, optional
        Segment map with radial quartiles encoded as seg*10 + radial_group.
    """

    # --------------------------------------------------
    # Determine downsampling
    # --------------------------------------------------
    files_lv = sorted(
        f for f in os.listdir(lv_mask_path)
        if f.lower().endswith((".tif", ".tiff", ".jp2"))
    )
    files_img = sorted(
        f for f in os.listdir(images_path)
        if f.lower().endswith((".tif", ".tiff", ".jp2"))
    )

    # Avoid zero division or zero scaling when image and mask slice counts differ greatly.
    # Compute ratio and ensure integer scaling factor >= 1
    mask_ratio = float(len(files_img)) / max(1, float(len(files_lv)))
    mask_factor = max(1, int(round(mask_ratio)))

    # Compute downsample as aha_factor / mask_factor, at least 1
    try:
        downsample = int(max(1, round(float(additional_downsample) * float(mask_factor))))
    except Exception:
        downsample = mask_factor

    # --------------------------------------------------
    # Load and resample LV mask
    # --------------------------------------------------
    lv_mask_vol = load_volume(lv_mask_path)
    lv_mask_vol_d = resample(
        lv_mask_vol,
        zoom=1 / additional_downsample,
        order=0,
    )

    # --------------------------------------------------
    # Long-axis definition
    # --------------------------------------------------
    MV_d = np.array(axis_points[0])[::-1] / downsample
    Apex_d = np.array(axis_points[1])[::-1] / downsample

    axis_vec = MV_d - Apex_d
    axis_vec /= np.linalg.norm(axis_vec)

    coords = np.array(np.nonzero(lv_mask_vol_d)).T
    centroid = (MV_d + Apex_d) / 2

    septum_point = np.asarray(septum_point) / additional_downsample

    # --------------------------------------------------
    # Long-axis projections
    # --------------------------------------------------
    projections = coords @ axis_vec

    # Discard extreme percentiles to guard againts noise and outliers in the mask
    pmin = np.percentile(projections, 2)
    pmax = np.percentile(projections, 98)
    length = pmax - pmin

    basal_thresh = pmin + 2 * length / 3.0
    mid_thresh = pmin + length / 3.0

    # --------------------------------------------------
    # Septum reference angle
    # --------------------------------------------------
    septum_vec = septum_point - centroid

    pts2d, u, v = project_to_plane(
        coords,
        centroid,
        axis_vec,
    )

    septum_point2d = np.array([
        np.dot(septum_vec, u),
        np.dot(septum_vec, v)
    ])

    ref_angle = np.arctan2(
        septum_point2d[1],
        septum_point2d[0]
    )

    angles = np.arctan2(pts2d[:, 1], pts2d[:, 0])

    angles_rel = (angles - ref_angle) % (2 * np.pi)
    angles_rel = (angles_rel + 2 * np.pi / 3) % (2 * np.pi)

    # --------------------------------------------------
    # Create segment map
    # --------------------------------------------------
    segment_map = np.zeros(
        lv_mask_vol_d.shape,
        dtype=np.uint16
    )

    def assign_segment(mask_idx, seg):
        segment_map[tuple(coords[mask_idx].T)] = seg

    basal_idx = np.where(projections >= basal_thresh)[0]

    mid_idx = np.where(
        (projections >= mid_thresh) &
        (projections < basal_thresh)
    )[0]

    apical_idx = np.where(
        (projections < mid_thresh) &
        (projections >= pmin + 0.02 * length)
    )[0]

    apex_idx = np.where(
        projections < (pmin + 0.02 * length)
    )[0]

    def sectors_from_indices(indices, n_sectors, start_seg):
        sector_width = 2 * np.pi / n_sectors

        for s in range(n_sectors):
            a0 = s * sector_width
            a1 = (s + 1) * sector_width

            mask_idx = indices[
                (angles_rel[indices] >= a0) &
                (angles_rel[indices] < a1)
            ]

            assign_segment(mask_idx, start_seg + s)

    # Basal: 1–6
    sectors_from_indices(basal_idx, 6, 1)

    # Mid: 7–12
    sectors_from_indices(mid_idx, 6, 7)

    # Apical: 13–16
    sector_width = 2 * np.pi / 4

    for s in range(4):
        a0 = s * sector_width
        a1 = (s + 1) * sector_width

        mask_idx = apical_idx[
            (angles_rel[apical_idx] >= a0) &
            (angles_rel[apical_idx] < a1)
        ]

        assign_segment(mask_idx, 13 + s)

    # Apex: 17
    if apex_idx.size > 0:
        assign_segment(apex_idx, 17)

    if not return_radial_map:
        return segment_map

    # --------------------------------------------------
    # Radial quartile subdivision
    # --------------------------------------------------
    dists = np.linalg.norm(coords - centroid, axis=1)

    segment_radial_map = np.zeros_like(
        segment_map,
        dtype=np.uint16
    )

    for seg_id in range(1, 18):

        seg_idx = np.where(
            segment_map[tuple(coords.T)] == seg_id
        )[0]

        if len(seg_idx) == 0:
            continue

        seg_dists = dists[seg_idx]

        q25, q50, q75 = np.percentile(
            seg_dists,
            [25, 50, 75]
        )

        for voxel_idx, d in zip(seg_idx, seg_dists):

            if d <= q25:
                radial_group = 1
            elif d <= q50:
                radial_group = 2
            elif d <= q75:
                radial_group = 3
            else:
                radial_group = 4

            segment_radial_map[
                tuple(coords[voxel_idx])
            ] = seg_id * 10 + radial_group

    n_slices = 20
    
    # 17-segment discrete colours (RGB)
    segment_colors = np.array([
        [166,206,227],[31,120,180],[178,223,138],[51,160,44],[251,154,153],
        [227,26,28],[253,191,111],[255,127,0],[202,178,214],[106,61,154],
        [255,255,153],[177,89,40],[141,211,199],[255,255,179],[190,186,218],
        [251,128,114],[128,177,211]
    ], dtype=np.uint8)
    
    # Step 1: project voxels along long axis
    
    apex_dir = -axis_vec
    
    projections_norm = (coords - centroid) @ apex_dir
    slice1 = np.percentile(projections_norm, 2)
    slicen = np.percentile(projections_norm, 98)
    slice_positions = np.linspace(slice1, slicen, n_slices)
    
    # Step 2: find two orthonormal vectors perpendicular to apex_dir
    v = apex_dir / np.linalg.norm(apex_dir)
    if abs(v[0]) < abs(v[1]):
        u1 = np.array([0, -v[2], v[1]])
    else:
        u1 = np.array([-v[2], 0, v[0]])
    u1 /= np.linalg.norm(u1)
    u2 = np.cross(v, u1)
    u2 /= np.linalg.norm(u2)
    
    # Step 3: extract slices
    slices_rgb = []
    for s in slice_positions:
        
        # select voxels near the slice plane
        dist = np.abs(projections_norm - s)
        plane_idx = np.where(dist < 2)[0]
        plane_coords = coords[plane_idx]
    
        if len(plane_coords) == 0:
            continue
    
        # project to 2D plane
        pts2d = np.column_stack([
            (plane_coords - centroid) @ u1,
            (plane_coords - centroid) @ u2
        ])
        x = np.round(pts2d[:,0] - pts2d[:,0].min()).astype(int)
        y = np.round(pts2d[:,1] - pts2d[:,1].min()).astype(int)
    
        # create RGB slice
        shape = (y.max()+1, x.max()+1, 3)
        rgb = np.zeros(shape, dtype=np.uint16)
        seg_ids = segment_map[plane_coords[:,0], plane_coords[:,1], plane_coords[:,2]]
        for seg_id in range(1, 18):
            mask = seg_ids == seg_id
            rgb[y[mask], x[mask]] = segment_colors[seg_id-1]
        slices_rgb.append(rgb)
    
    # Step 4: plot slices in 2 rows with legend
    rows = 2
    cols = (len(slices_rgb)+1)//2
    fig, axes = plt.subplots(rows, cols, figsize=(3*cols, 3*rows))
    axes = axes.flatten()
    
    for i, rgb in enumerate(slices_rgb):
        axes[i].imshow(rgb)
        axes[i].set_title(f"Slice {i+1}")
        axes[i].axis('off')
    
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    
    # Legend
    patches = [Patch(facecolor=segment_colors[i]/255, label=f"Seg {i+1}") for i in range(17)]
    fig.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.1, 0.95), title="AHA Segments")
    
    plt.tight_layout()
    plt.show()

    return segment_map, segment_radial_map



def resample(
    vol,
    target_shape=None,
    zoom=None,
    order=0,
    TwoD=False,
    is_mask=False
):
    """
    Resample 2D or 3D volume.

    Parameters
    ----------
    vol : ndarray
        Input array.
    target_shape : tuple, optional
        Desired output shape.
    zoom : float or tuple, optional
        Zoom factor(s).
    order : int
        Interpolation order (ignored for masks).
    TwoD : bool
        If True, treat input as 2D.
    is_mask : bool
        If True, use nearest-neighbour mask-safe resampling.

    Returns
    -------
    ndarray
        Resampled volume.
    """
    dtype = vol.dtype
    vol = np.asarray(vol).astype(dtype)

    # ------------------ Determine zoom factors ------------------
    if TwoD:
        if target_shape is not None:
            zoom_factors = (
                target_shape[0] / vol.shape[0],
                target_shape[1] / vol.shape[1]
            )
        elif zoom is not None:
            zoom_factors = (zoom, zoom)
        else:
            raise ValueError("resample requires target_shape or zoom")
    else:
        if target_shape is not None:
            zoom_factors = (
                target_shape[0] / vol.shape[0],
                target_shape[1] / vol.shape[1],
                target_shape[2] / vol.shape[2]
            )
        elif zoom is not None:
            zoom_factors = (zoom, zoom, zoom)
        else:
            raise ValueError("resample requires target_shape or zoom")

    # ------------------ MASK SAFE PATH ------------------
    if is_mask:
        if not TwoD:
            raise ValueError("Mask resampling currently supports 2D only")

        new_h = int(round(vol.shape[0] * zoom_factors[0]))
        new_w = int(round(vol.shape[1] * zoom_factors[1]))

        res = cv2.resize(
            vol,
            (new_w, new_h),
            interpolation=cv2.INTER_NEAREST
        )

        return res.astype(vol.dtype)

    # ------------------ IMAGE PATH ------------------
    #vol = vol.astype(np.float32)
    vol = vol.astype(dtype)
    return ndi.zoom(vol, zoom_factors, order=order)




def project_to_plane(points, origin, normal):
    #Project 3D points to plane defined by two vectors. Returns 2D coordinates in plane basis.
    
    # build orthonormal basis u,v for plane
    normal = normal / np.linalg.norm(normal)
    
    # choose arbitrary vector not parallel to normal
    arbitrary = np.array([1.0,0.0,0.0])
    if np.allclose(np.abs(np.dot(arbitrary, normal)), 1.0, atol=1e-3):
        arbitrary = np.array([0.0,1.0,0.0])
    u = np.cross(normal, arbitrary)
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    v = v / np.linalg.norm(v)
    rel = points - origin
    x = rel @ u
    y = rel @ v
    return np.column_stack([x, y]), u, v


def newAxes(vector1, vector2=None):
    vector1 = np.asarray(vector1, dtype=float)
    vector1 /= np.linalg.norm(vector1)

    if vector2 is None:
        # Choose a safe arbitrary axis
        if abs(vector1[1]) < 0.9:
            tmp = np.array([0, 1, 0])
        else:
            tmp = np.array([1, 0, 0])

        perp1 = np.cross(vector1, tmp)
    else:
        perp1 = np.asarray(vector2, dtype=float)

    perp1 /= np.linalg.norm(perp1)

    perp2 = np.cross(vector1, perp1)
    perp2 /= np.linalg.norm(perp2)

    return vector1, perp1, perp2



def newDimensions(vector, dimensions):
    depth, height, width = dimensions

    # Compute perpendicular axes
    yax = np.array([0, 1, 0])
    perp1 = np.cross(vector, yax) 
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(vector, perp1)
    perp2 /= np.linalg.norm(perp2)

    # Original volume corners (z, y, x)
    corners = np.array([
        [0, 0, 0],
        [0, 0, width],
        [0, height, 0],
        [0, height, width],
        [depth, 0, 0],
        [depth, 0, width],
        [depth, height, 0],
        [depth, height, width]
    ])

    # Center the volume
    centre = np.array([depth // 2, height // 2, width // 2])
    shifted = corners - centre  # shape (8, 3)

    # Project into new coordinate system
    new_coords = np.dot(shifted, np.vstack([vector, perp2, perp1]).T)  # shape (8, 3)
    
    # Find bounding box in new coordinates
    mins = np.min(new_coords, axis=0)
    maxs = np.max(new_coords, axis=0)
    sizes = maxs - mins

    # Output dimensions: depth, height, width
    out_d = int(np.ceil(sizes[0]))
    out_h = int(np.ceil(sizes[1]))
    out_w = int(np.ceil(sizes[2]))

    print(out_d, out_h, out_w)
    return np.array([out_d, out_h, out_w])



def stackSubset(folder, index, vector):
    # Open the first file in the folder to work out the dimensions
    files = sorted([f for f in os.listdir(folder) if f.lower().endswith(('.png', '.jpg', '.tif', '.tiff'))])
    im = imread(os.path.join(folder, files[0]))#
    depth = len(files)
    height, width = im.shape
    dimensions = [depth, height, width]
    newDims = newDimensions(vector, dimensions)
    
    # Work out size of subset required
    rang = newDims[0] - dimensions[0]
    
    half = int(round(rang / 1.2))

    p1 = index - half
    p2 = index + half
    
    start = min(p1, p2)
    end = max(p1, p2)
    
    # clamp to valid index range
    start = max(start, 0)
    end = min(end, len(files) - 1)
    
    print(start, end)
    print(rang)
    
    # Open necessary images
    stack = [imread(os.path.join(folder, f)) for f in files[start:end]]
    stack = np.stack(stack, axis=0)
    print('Stack loaded')
    
    return stack, dimensions, newDims



def getRange(folder, index, vector):
    # Open the first file in the folder to work out the dimensions
    files = sorted([f for f in os.listdir(folder) if f.lower().endswith(('.png', '.jpg', '.tif', '.tiff'))])
    im = imread(os.path.join(folder, files[0]))#
    depth = len(files)
    height, width = im.shape
    dimensions = [depth, height, width]
    newDims = newDimensions(vector, dimensions)
    
    # Work out size of subset required
    rang = newDims[0] - dimensions[0]
    
    half = int(round(rang / 1.2))

    p1 = index - half
    p2 = index + half
    
    start = min(p1, p2)
    end = max(p1, p2)

    
    print(start, end)
    print(rang)
    
    # Open necessary images
    stack = [imread(os.path.join(folder, f)) for f in files[start:end]]
    stack = np.stack(stack, axis=0)
    print('Stack loaded')
    
    return np.array([start, end]), dimensions, newDims




def ResliceOneSlice(vector, c_new, in_folder=None, stack=None):
    """
    Reslice an oblique slice through a given 3D point c_new (z, y, x),
    loading only a Z-subset where possible.
    """

    # --- build slice frame ---
    n, perp1, perp2 = newAxes(vector)

    if in_folder is not None: 

        # --- load one image to get volume shape ---
        files = sorted(
            f for f in os.listdir(in_folder)
            if f.lower().endswith(('.png', '.jpg', '.tif', '.tiff'))
        )
        im0 = imread(os.path.join(in_folder, files[0]))
        depth = len(files)
        height, width = im0.shape
        dimensions = [depth, height, width]
    
        # --- output slice size ---
        out_d, out_h, out_w = newDimensions(n, dimensions)
    
        # --- compute required Z slab ---
        z_radius = (
            abs(perp1[0]) * (out_w / 2) +
            abs(perp2[0]) * (out_h / 2)
        )
    
        z_min = int(np.floor(c_new[0] - z_radius))
        z_max = int(np.ceil(c_new[0] + z_radius))
    
        z_min_clamped = max(0, z_min)
        z_max_clamped = min(depth - 1, z_max)
    
        stack = np.stack(
            [
                imread(os.path.join(in_folder, files[z]))
                for z in range(z_min_clamped, z_max_clamped + 1)
            ],
            axis=0
        )
        z_offset = z_min_clamped
    
    else: 
        depth, height, width = stack.shape
        dimensions = [depth, height, width]
    
        # --- output slice size ---
        out_d, out_h, out_w = newDimensions(n, dimensions)
    
        # --- compute required Z slab ---
        z_radius = (
            abs(perp1[0]) * (out_w / 2) +
            abs(perp2[0]) * (out_h / 2)
        )
    
        z_min = int(np.floor(c_new[0] - z_radius))
        z_max = int(np.ceil(c_new[0] + z_radius))
    
        z_min_clamped = max(0, z_min)
        z_max_clamped = min(depth - 1, z_max)
    
        z_offset = z_min_clamped
        z_offset = 0
        

    # --- slice center in stack coordinates ---
    s_centre = np.array([
        c_new[0] - z_offset,
        c_new[1],
        c_new[2]
    ], dtype=float)

    # --- sample coordinates ---
    coords = []
    for y in range(out_h):
        for x in range(out_w):
            offset = (
                (x - out_w // 2) * perp1 +
                (y - out_h // 2) * -perp2
            )
            coords.append(s_centre + offset)

    coords = np.array(coords).T  # (3, N)

    slice_vals = map_coordinates(
        stack,
        coords,
        order=0,
        mode="constant",
        cval=0.0
    ).reshape(out_h, out_w)

    # --- geometry (volume coordinates!) ---
    slice_geometry = {
        "slice_center_3d": np.asarray(c_new, dtype=float),
        "axis_x_3d": perp1.astype(float),
        "axis_y_3d": (-perp2).astype(float),
        "normal_3d": n.astype(float),
        "width": out_w,
        "height": out_h,
        "z_slab": (z_min_clamped, z_max_clamped),
    }

    return slice_vals, slice_geometry


def ResliceMultipleSlices(vector, centers, in_folder=None, stack=None):
    """
    Reslice multiple oblique slices through given 3D points (z, y, x),
    using the same axis vector. Loads one volume slab containing all slices.
    
    Parameters
    ----------
    vector : array-like
        Normal vector of the slices.
    centers : list of array-like
        List of slice centers (z, y, x) in volume coordinates.
    in_folder : str, optional
        Folder containing image stack (if stack not provided).
    stack : np.ndarray, optional
        Preloaded volume (depth, height, width).

    Returns
    -------
    slices : list of np.ndarray
        List of resliced 2D slices.
    geometries : list of dict
        List of slice geometry dictionaries for each slice.
    """

    n, perp1, perp2 = newAxes(vector)
    centers = np.array(centers, dtype=float)

    # --- Load volume slab if folder is provided ---
    if in_folder is not None:
        files = sorted(
            f for f in os.listdir(in_folder)
            if f.lower().endswith(('.png', '.jpg', '.tif', '.tiff'))
        )
        im0 = imread(os.path.join(in_folder, files[0]))
        depth, height, width = len(files), *im0.shape
        dimensions = [depth, height, width]

        # Determine slice output size
        out_d, out_h, out_w = newDimensions(n, dimensions)

        # Determine Z slab needed to cover all slices
        z_offsets = (
            abs(perp1[0]) * (out_w / 2) +
            abs(perp2[0]) * (out_h / 2)
        )

        z_min = int(np.floor(np.min(centers[:, 0]) - z_offsets))
        z_max = int(np.ceil(np.max(centers[:, 0]) + z_offsets))

        z_min_clamped = max(0, z_min)
        z_max_clamped = min(depth - 1, z_max)

        stack = np.stack(
            [imread(os.path.join(in_folder, files[z]))
             for z in range(z_min_clamped, z_max_clamped + 1)],
            axis=0
        )
        z_offset = z_min_clamped

    else:
        depth, height, width = stack.shape
        dimensions = [depth, height, width]
        out_d, out_h, out_w = newDimensions(n, dimensions)
        z_offset = 0

    slices = []
    geometries = []

    # --- Reslice each center ---
    for c_new in centers:
        s_centre = np.array([
            c_new[0] - z_offset,
            c_new[1],
            c_new[2]
        ], dtype=float)

        # Compute sample coordinates
        coords = []
        for y in range(out_h):
            for x in range(out_w):
                offset = (x - out_w // 2) * perp1 + (y - out_h // 2) * -perp2
                coords.append(s_centre + offset)
        coords = np.array(coords).T  # (3, N)

        slice_vals = map_coordinates(
            stack,
            coords,
            order=0,
            mode="constant",
            cval=0.0
        ).reshape(out_h, out_w)

        slices.append(slice_vals)

        slice_geometry = {
            "slice_center_3d": c_new.astype(float),
            "axis_x_3d": perp1.astype(float),
            "axis_y_3d": (-perp2).astype(float),
            "normal_3d": n.astype(float),
            "width": out_w,
            "height": out_h,
            "z_slab": (z_min_clamped, z_max_clamped),
        }
        geometries.append(slice_geometry)

    return slices, geometries



def GetResliceGeometry(vector, centers, in_folder=None, stack=None):
    n, perp1, perp2 = newAxes(vector)
    centers = np.array(centers, dtype=float)

    # --- Load volume slab if folder is provided ---
    if in_folder is not None:
        files = sorted(
            f for f in os.listdir(in_folder)
            if f.lower().endswith(('.png', '.jpg', '.tif', '.tiff'))
        )
        im0 = imread(os.path.join(in_folder, files[0]))
        depth, height, width = len(files), *im0.shape
        dimensions = [depth, height, width]

        # Determine slice output size
        out_d, out_h, out_w = newDimensions(n, dimensions)

        # Determine Z slab needed to cover all slices
        z_offsets = (
            abs(perp1[0]) * (out_w / 2) +
            abs(perp2[0]) * (out_h / 2)
        )

        z_min = int(np.floor(np.min(centers[:, 0]) - z_offsets))
        z_max = int(np.ceil(np.max(centers[:, 0]) + z_offsets))

        z_min_clamped = max(0, z_min)
        z_max_clamped = min(depth - 1, z_max)

    else:
        depth, height, width = stack.shape
        dimensions = [depth, height, width]
        out_d, out_h, out_w = newDimensions(n, dimensions)

    geometries = []

    # --- Reslice each center ---
    for c_new in centers:
        slice_geometry = {
            "slice_center_3d": c_new.astype(float),
            "axis_x_3d": perp1.astype(float),
            "axis_y_3d": (-perp2).astype(float),
            "normal_3d": n.astype(float),
            "width": out_w,
            "height": out_h,
            "z_slab": (z_min_clamped, z_max_clamped),
        }
        geometries.append(slice_geometry)
    
    return geometries



# Two vector reslice for long axis - septum point and MV axis 
def Reslice3Points(p1, p2, p3, in_folder=None, stack=None):
    """
    Reslice using either 3 points (p1, p2 define main axis, p3 defines plane)

    Points are (z, y, x)
    """

    # --- build slice frame ---
    
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    p3 = np.asarray(p3, dtype=float)

    # main axis (replaces vector)
    perp1 = p2 - p1
    perp1 = perp1 / np.linalg.norm(perp1)

    # second in-plane vector from p3
    v2 = p3 - p1

    # normal
    n = np.cross(perp1, v2)
    norm_n = np.linalg.norm(n)
    if norm_n < 1e-6:
        raise ValueError("Points are collinear or too close")

    n = n / norm_n

    # orthogonal second axis
    perp2 = np.cross(n, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    # center = midpoint of p1->p2 projected into plane
    c_new = (p1 + p2) / 2.0


    # --- load data / dimensions ---
    if in_folder is not None:

        files = sorted(
            f for f in os.listdir(in_folder)
            if f.lower().endswith(('.png', '.jpg', '.tif', '.tiff'))
        )

        im0 = imread(os.path.join(in_folder, files[0]))
        depth = len(files)
        height, width = im0.shape
        dimensions = [depth, height, width]

        out_d, out_h, out_w = newDimensions(n, dimensions)

        # --- compute required Z slab ---
        z_radius = (
            abs(perp1[0]) * (out_w / 2) +
            abs(perp2[0]) * (out_h / 2)
        )

        z_min = int(np.floor(c_new[0] - z_radius))
        z_max = int(np.ceil(c_new[0] + z_radius))

        z_min_clamped = max(0, z_min)
        z_max_clamped = min(depth - 1, z_max)

        stack = np.stack(
            [
                imread(os.path.join(in_folder, files[z]))
                for z in range(z_min_clamped, z_max_clamped + 1)
            ],
            axis=0
        )

        z_offset = z_min_clamped

    else:
        depth, height, width = stack.shape
        dimensions = [depth, height, width]

        out_d, out_h, out_w = newDimensions(n, dimensions)

        z_radius = (
            abs(perp1[0]) * (out_w / 2) +
            abs(perp2[0]) * (out_h / 2)
        )

        z_min = int(np.floor(c_new[0] - z_radius))
        z_max = int(np.ceil(c_new[0] + z_radius))

        z_min_clamped = max(0, z_min)
        z_max_clamped = min(depth - 1, z_max)

        z_offset = z_min_clamped

    # --- slice center in stack coords ---
    s_centre = np.array([
        c_new[0] - z_offset,
        c_new[1],
        c_new[2]
    ], dtype=float)

    # --- sample coordinates ---
    coords = []
    for y in range(out_h):
        for x in range(out_w):
            offset = (
                (x - out_w // 2) * perp1 +
                (y - out_h // 2) * -perp2
            )
            coords.append(s_centre + offset)

    coords = np.array(coords).T

    slice_vals = map_coordinates(
        stack,
        coords,
        order=0,
        mode="constant",
        cval=0.0
    ).reshape(out_h, out_w)

    # --- geometry ---
    slice_geometry = {
        "slice_center_3d": np.asarray(c_new, dtype=float),
        "axis_x_3d": perp1.astype(float),
        "axis_y_3d": (-perp2).astype(float),
        "normal_3d": n.astype(float),
        "width": out_w,
        "height": out_h,
        "z_slab": (z_min_clamped, z_max_clamped),
    }

    return slice_vals, slice_geometry



def line_xy_at_z(p0, p1, z):
    """
    p0, p1: (x, y, z) endpoints of the line
    z: target z value

    returns: (x, y) at that z
    """
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)

    z0, z1 = p0[2], p1[2]

    if z0 == z1:
        raise ValueError("Line is parallel to z-plane (z0 == z1)")

    t = (z - z0) / (z1 - z0)

    x = p0[0] + t * (p1[0] - p0[0])
    y = p0[1] + t * (p1[1] - p0[1])

    return x, y


def project_3d_to_slice(p3d, geom):
    """
    Project a 3D point into 2D slice coordinates using slice geometry.
    """
    p3d = np.asarray(p3d, dtype=float)

    centre = geom["slice_center_3d"]
    axis_x = geom["axis_x_3d"]
    axis_y = geom["axis_y_3d"]
    out_w = geom["width"]
    out_h = geom["height"]

    # Vector from slice center to point
    v = p3d - centre

    # Project onto slice axes
    x = np.dot(v, axis_x) + out_w / 2
    y = np.dot(v, axis_y) + out_h / 2

    return x, y



def crop_from_mask(mask, border=5):
    """
    Compute crop slices that remove all-zero rows/cols,
    keeping a small border.
    """
    rows = np.any(mask > 0, axis=1)
    cols = np.any(mask > 0, axis=0)

    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    rmin = max(rmin - border, 0)
    cmin = max(cmin - border, 0)
    rmax = min(rmax + border + 1, mask.shape[0])
    cmax = min(cmax + border + 1, mask.shape[1])

    return slice(rmin, rmax), slice(cmin, cmax)

def normalise_percentile(img, p_low=1, p_high=99, out_range=(0, 1)):
    img = img.astype(np.float32)

    lo = np.nanpercentile(img, p_low)
    hi = np.nanpercentile(img, p_high)

    img_n = (img - lo) / (hi - lo)
    img_n = np.clip(img_n, 0, 1)

    out_lo, out_hi = out_range
    return img_n * (out_hi - out_lo) + out_lo

################################# Transumral PLots



def get_epi_contour(mask2d, close_radius=2, shrink_radius=2):
    se_close = disk(close_radius)
    closed = closing(mask2d, se_close)
    closed = remove_small_objects(closed.astype(bool), 500).astype(np.uint8)
    
    # Optional shrink
    if shrink_radius > 0:
        se_shrink = disk(shrink_radius)
        closed = erosion(closed, se_shrink).astype(np.uint8)

    # Find contour
    contours, _ = cv2.findContours(closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None, closed

    cnt = max(contours, key=cv2.contourArea)
    return cnt.squeeze(), closed


def get_epi_contour(mask2d):
    H, W = mask2d.shape
    scale = max(int(0.002 * min(H, W)), 2)
    
    closed = closing(mask2d, disk(6 * scale))
    closed = remove_small_objects(closed.astype(bool), 0.001 * H * W)
    closed = closed.astype(np.uint8)

    contours, _ = cv2.findContours(closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None, closed

    cnt = max(contours, key=cv2.contourArea)
    return cnt.squeeze(), closed


def get_endo_contour(mask2d, epi_cnt):
    H, W = mask2d.shape
    scale = max(int(0.002 * min(H, W)), 2)

    # Epicardial mask
    epi_mask = np.zeros_like(mask2d, dtype=np.uint8)
    cv2.fillPoly(epi_mask, [epi_cnt.astype(int)], 1)

    # Inside myocardium
    inside = mask2d * epi_mask

    # Candidate cavity
    mask2d[mask2d>0] = 1
    inv = (1 - mask2d)
    inv[epi_mask == 0] = 0
    inv = remove_small_objects(inv.astype(bool), 0.0005 * H * W)
    inv = inv.astype(np.uint8)

    # Smooth cavity
    inv = closing(inv, disk(3 * scale))
    inv = opening(inv, disk(scale))
    inv = remove_small_holes(inv.astype(bool), area_threshold=0.0005 * H * W)
    inv = inv.astype(np.uint8)
    

    # Extract contours
    contours, _ = cv2.findContours(inv, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None, inv

    
    # Epicardial centroid
    epi_center = epi_cnt.mean(axis=0)

    # Choose cavity: closest centroid to epi center
    best_cnt = None
    best_dist = np.inf

    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area < 0.001 * cv2.contourArea(epi_cnt):
            continue  # too small to be LV

        center = cnt.squeeze().mean(axis=0)
        dist = np.linalg.norm(center - epi_center)

        # Reject anything touching epicardium
        #if any(cv2.pointPolygonTest(epi_cnt, tuple(pt[0]), False) >= 0 for pt in cnt):
        #    continue

        if dist < best_dist:
            best_dist = dist
            best_cnt = cnt

    if best_cnt is None:
        return None, inv

    return best_cnt.squeeze(), inv


def trim_lines(lines, epi_pts, endo_pts=None, centroid2d=None):
    trimmed_lines = []

    # Prepare epicardial contour
    epi_cnt = epi_pts.astype(np.float32).reshape(-1, 1, 2)
    
    for start, end in lines:
        xs = np.linspace(start[0], end[0], 1000)
        ys = np.linspace(start[1], end[1], 1000)
        points = np.column_stack([xs, ys])

        # Determine start and end points for trimming
        if endo_pts is not None:
            endo_cnt = endo_pts.astype(np.float32).reshape(-1, 1, 2)
            # First point inside endo
            mask_endo = np.array([cv2.pointPolygonTest(endo_cnt, (y,x), False)>=0 for y,x in points])
            if np.any(mask_endo):
                start_pt = points[np.argmin(mask_endo)]
            else:
                start_pt = points[0] if centroid2d is None else centroid2d
        else:
            start_pt = points[0] if centroid2d is None else centroid2d

        # Last point inside epi
        mask_epi = np.array([cv2.pointPolygonTest(epi_cnt, (y,x), False)>=0 for y,x in points])
        if np.any(mask_epi):
            end_pt = points[np.argmax(mask_epi[::-1] != 0)]  # last True point
            end_pt = points[len(points)-1 - np.argmax(mask_epi[::-1])]
        else:
            continue  # skip line if it doesn't reach epicardium

        trimmed_lines.append((tuple(start_pt), tuple(end_pt)))

    return trimmed_lines



def trim_lines(lines, epi_pts, endo_pts=None, centroid2d=None):
    trimmed_lines = []

    epi_cnt = epi_pts.astype(np.float32).reshape(-1, 1, 2)

    for start, end in lines:
        xs = np.linspace(start[0], end[0], 1000)
        ys = np.linspace(start[1], end[1], 1000)
        points = np.column_stack([xs, ys])

        # Check inside epicardium
        mask_epi = np.array([
            cv2.pointPolygonTest(epi_cnt, (float(y), float(x)), False) >= 0
            for x, y in points
        ])

        if not np.any(mask_epi):
            continue  # no intersection at all

        # --- Find continuous True segments ---
        segments = []
        in_segment = False
        seg_start = 0

        for i, val in enumerate(mask_epi):
            if val and not in_segment:
                in_segment = True
                seg_start = i
            elif not val and in_segment:
                segments.append((seg_start, i))
                in_segment = False

        if in_segment:
            segments.append((seg_start, len(mask_epi)))

        # --- Choose best segment ---
        # Typically the longest segment inside epi
        seg_lengths = [end - start for start, end in segments]
        best_seg_idx = np.argmax(seg_lengths)
        seg_start, seg_end = segments[best_seg_idx]

        trimmed_segment = points[seg_start:seg_end]

        if len(trimmed_segment) < 2:
            continue

        start_pt = tuple(trimmed_segment[0])
        end_pt = tuple(trimmed_segment[-1])

        # Optional: enforce endocardium trimming (inner boundary)
        if endo_pts is not None:
            endo_cnt = endo_pts.astype(np.float32).reshape(-1, 1, 2)
            mask_endo = np.array([
                cv2.pointPolygonTest(endo_cnt, (float(y), float(x)), False) >= 0
                for x, y in trimmed_segment
            ])

            if np.any(mask_endo):
                first_inside = np.argmax(mask_endo)
                start_pt = tuple(trimmed_segment[first_inside])

        trimmed_lines.append((start_pt, end_pt))

    return trimmed_lines

def trim_lines(lines, epi_pts, endo_pts=None, centroid2d=None):
    trimmed_lines = []

    # Prepare epicardial contour
    epi_cnt = epi_pts.astype(np.float32).reshape(-1, 1, 2)
    
    for start, end in lines:
        xs = np.linspace(start[0], end[0], 1000)
        ys = np.linspace(start[1], end[1], 1000)
        points = np.column_stack([xs, ys])
        
        mask_epi = np.array([cv2.pointPolygonTest(epi_cnt, (y,x), False)>=0 for y,x in points])
        
        # Determine start and end points for trimming
        if endo_pts is not None:
            endo_cnt = endo_pts.astype(np.float32).reshape(-1, 1, 2)
            mask_endo = np.array([cv2.pointPolygonTest(endo_cnt, (y,x), False)>=0 for y,x in points[mask_epi]])
            # First point inside endo
            if np.any(mask_endo):
                start_pt = points[np.argmin(mask_endo)]
            else:
                start_pt = points[np.argmax(mask_epi)] 
                #start_pt = points[0] if centroid2d is None else centroid2d
        else:
            start_pt = points[np.argmax(mask_epi)] 
            #start_pt = points[0] if centroid2d is None else centroid2d

        # Last point inside epi
        if np.any(mask_epi):
            end_pt = points[np.argmax(mask_epi[::-1] != 0)]  # last True point
            end_pt = points[len(points)-1 - np.argmax(mask_epi[::-1])]
        else:
            continue  # skip line if it doesn't reach epicardium

        trimmed_lines.append((tuple(start_pt), tuple(end_pt)))

    return trimmed_lines


def filter_lines_by_length(lines, tol_frac=0.5):
    """
    Remove lines whose length is outside median ± tol_frac * median.

    Parameters
    ----------
    lines : list of ((x1,y1), (x2,y2))
        Trimmed line segments.
    tol_frac : float
        Fraction of median length to tolerate (default 0.5).

    Returns
    -------
    filtered_lines : list
        Lines within allowed length range.
    """

    if len(lines) == 0:
        return []

    # Compute lengths
    lengths = np.array([
        np.linalg.norm(np.asarray(p2) - np.asarray(p1))
        for p1, p2 in lines
    ])

    med = np.median(lengths)
    lo = med * (1.0 - tol_frac)
    hi = med * (1.0 + tol_frac)

    # Filter
    filtered_lines = [
        line for line, L in zip(lines, lengths)
        if lo <= L <= hi
    ]

    return filtered_lines




def profile_intensities(img_helix, radial_lines, factor=1.0):
    """
    Profiles intensities from helical image along radial lines.

    Parameters:
        img_helix (ndarray): The helical angle image.
        radial_lines(list): List of (start, end) tuples corresponding to radial lines
        factor (float): Scaling factor for line coordinates.

    Returns:
        dict: Grouped intensity profiles by region name.
    """
    #img_helix = np.nan_to_num(img_helix, nan=0.0)

    intensity_profiles = []
    good_lines = []
    for i, (start_pt, end_pt) in enumerate(radial_lines):
        print(f"Measuring line {i+1}/{len(radial_lines)}")

        # Scale and flip (x, y) → (row, col) = (y, x)
        start = [start_pt[1] * factor, start_pt[0] * factor]  # (row, col)
        end   = [end_pt[1]   * factor, end_pt[0]   * factor]

        intensity_profile = measure.profile_line(img_helix, start, end, order=0, mode='constant', cval=np.nan)
        
        n_points = len(intensity_profile)
        if n_points == 0 or (np.count_nonzero(np.isnan( intensity_profile)) / n_points) > 0.1:
            continue
        else: 
            good_lines.append([start_pt, end_pt])

        # Convert to helical angle
        #angle_profile = (intensity_profile * 180 / 255) - 90
        
        intensity_profile = np.asarray(intensity_profile, dtype=float)
        
        angle_profile = np.where(
            (intensity_profile != 0) & (~np.isnan(intensity_profile)),
            (intensity_profile * 180 / 255) - 90,
            np.nan
        )
        
        intensity_profiles.append(angle_profile)

    return intensity_profiles, good_lines





def transmuralSampling(centres, vector, slices=None, seg_map=None, centroids = None, mask_factor=None, n_radial=50, route_slices=None):
    
    all_profiles = []
    all_good_lines = []
    all_segments_assigned = []
    all_route_slices = []
    
    if centroids is None: 
        centroids = np.zeros(len(centres))
        
    if route_slices is None: 
        route_slices = np.zeros(len(centres))
    
    for slic, C, centroid, route_idx in zip(slices, centres, centroids, route_slices): 
    
        # ------------------- Extract 2D Mask Slice -------------------
        if seg_map is None:
            continue
        else:
            mf = mask_factor if mask_factor is not None else 1
            C_new = C / mf
            
            mask2d_sml, tmp = ResliceOneSlice(vector, C_new, stack=seg_map)
            mask2d = resample(mask2d_sml, target_shape=slic.shape, TwoD=True, order=0, is_mask=True)
            
            plt.imshow(mask2d_sml)
            plt.show()
            
            plt.imshow(mask2d)
            plt.show()
            
        mask_slice = slic > 0 
        mask2d[~mask_slice] = 0
        
        
        # ------------------- Epicardium -------------------
        epi_cnt, epi_closed = get_epi_contour(mask2d)
        if epi_cnt is None:
            continue
        
        # Convert epi contour to same 2D coordinate system as pts2d
        epi_pts = np.column_stack([epi_cnt[:,0],
                                   epi_cnt[:,1]])
        
        # ------------------- Endocardium -------------------
        endo_cnt, endo_clean = get_endo_contour(mask2d, epi_cnt)
        if endo_cnt is not None:
            endo_pts = np.column_stack([endo_cnt[:,0],
                                        endo_cnt[:,1]])
        else:
            endo_pts = None
        
        # ------------------- Verification -------------------
        
        fig, ax = plt.subplots(1,3, figsize=(15,5))
        ax[0].imshow(mask2d, cmap='gray')
        ax[0].set_title("Raw Mask Slice")
        
        ax[1].imshow(epi_closed, cmap='gray')
        ax[1].plot(epi_cnt[:,0], epi_cnt[:,1], 'r-')
        ax[1].set_title("Epicardium")
        
        ax[2].imshow(endo_clean, cmap='gray')
        if endo_pts is not None:
            ax[2].plot(endo_cnt[:,0], endo_cnt[:,1], 'b-')
        ax[2].set_title("Endocardium")
        plt.show()
        
        # ------------------- Compute Centroid if not given-------------------
        if centroid == 0: 
            if endo_pts is not None:
                centroid2d = (epi_pts.mean(axis=0) + endo_pts.mean(axis=0)) / 2
            else: 
                centroid2d = epi_pts.mean(axis=0)
            
            #centroid2d = centroid2d[::-1]
        else:
            centroid2d = np.array(centroid)
    
        # ------------------- Radial lines from centroid -------------------
            
        n_radial = 48  # number of radial lines
        angles = np.linspace(0, 2*np.pi, n_radial, endpoint=False)
        
        mean_radius = np.mean(np.linalg.norm(epi_pts - centroid2d[None,:], axis=1))
        L = mean_radius * 2  # extend a bit beyond epicardium
        
        p1_list, p2_list = [], []
        for theta in angles:
            direction = np.array([np.cos(theta), np.sin(theta)])
            p1_list.append(centroid2d)             # start at center
            p2_list.append(centroid2d + direction*L)  # extend outwards
        
        
        # ------------------- Trim using Epicardium -------------------
        trimmed_lines = trim_lines(list(zip(p1_list,p2_list)), epi_pts, endo_pts, centroid2d)
        trimmed_lines = filter_lines_by_length(trimmed_lines, tol_frac=0.5)
        
        # ------------------- Sample HA along radial lines -------------------
        profiles, good_lines = profile_intensities(slic, trimmed_lines)            
        
        all_profiles.extend(profiles)
        all_good_lines.extend(good_lines)
        
        # Segment assignment
        segments_assigned = []
        
        mx = []
        my = []
        
        p1x = []
        p1y = []
        
        p2x = [] 
        p2y = []
        
        route_slices = []
        for (pt1, pt2) in good_lines:
            pt1 = np.asarray(pt1, dtype=float)
            pt2 = np.asarray(pt2, dtype=float)
            mid_pt = 0.5 * (pt1 + pt2)   # [x, y]
            
            p1x.append(pt1[0])
            p1y.append(pt1[1])
            
            p2x.append(pt2[0])
            p2y.append(pt2[1])
            
            scale_x = mask2d_sml.shape[1] / slic.shape[1]
            scale_y = mask2d_sml.shape[0] / slic.shape[0]
            
            mid_pt_sml = np.array([
                mid_pt[0] * scale_x,
                mid_pt[1] * scale_y
            ])
            
            mx.append(mid_pt_sml[0])
            my.append(mid_pt_sml[1])
            
            seg = mask2d_sml[int(mid_pt_sml[1]), int(mid_pt_sml[0])]            
            segments_assigned.append(seg)
            
            #seg = mask2d[int(mid_pt[1]), int(mid_pt[0])] 
            #segments_assigned.append(seg)
            
            route_slices.append(route_idx)
        
        all_segments_assigned.extend(segments_assigned)
        all_route_slices.extend(route_slices)
        
        plt.imshow(mask2d_sml)
        plt.scatter(mx, my, color='g')
        plt.show()
        
        plt.imshow(mask2d)
        plt.scatter(p1x, p1y, color='r')
        plt.scatter(p2x, p2y, color='b')
        plt.show()

        # ------------------- Plot 2D slice -------------------
        
        plt.figure(figsize=(7,7))
        img = plt.imshow(slic, origin='upper', aspect='equal')
        for p1,p2 in good_lines:
            plt.plot([p1[0],p2[0]],[p1[1],p2[1]],'r-')
        plt.plot(epi_pts[:,0], epi_pts[:,1],'b-', lw=2)
        cbar = plt.colorbar(img)
        
        # Set ticks in LUT space (0–255)
        lut_ticks = np.linspace(0, 255, 7)
        
        # Convert ticks to physical angle space (-90..90)
        angle_ticks = np.linspace(-90, 90, 7)
        
        cbar.set_ticks(lut_ticks)
        cbar.set_ticklabels([f"{a:.0f}" for a in angle_ticks])
        cbar.set_label("HA (deg)")
        
        plt.title(f"Slic {np.round(C_new[0])} HA with Radial Lines")
        plt.show()
    
        
    return all_profiles, all_good_lines, all_segments_assigned, all_route_slices



def plot_transmural_slice(profiles, savepath=None):
    # --------------------------
    # Pad profiles to equal size
    # --------------------------
    max_len = max(len(p) for p in profiles)
    # Pad front to line up with epicardium
    #padded = [np.pad(p, (max_len - len(p), 0), mode="constant", constant_values=np.nan) for p in profiles]
    padded = [np.pad(p, (max_len - len(p), 0), mode="constant", constant_values=np.nan) for p in profiles]
    data = np.vstack(padded)    # shape (Nprofiles × Ndepth)

    # --------------------------
    # Circular stats on [-90°, 90°]
    # --------------------------
    data_rad = np.deg2rad(data)

    mean_rad = circmean(
        data_rad, axis=0,
        low=-np.pi/2, high=np.pi/2,
        nan_policy="omit"
    )

    std_rad = circstd(
        data_rad, axis=0,
        low=-np.pi/2, high=np.pi/2,
        nan_policy="omit"
    )

    # Range 0–1 for transmural depth
    r = np.linspace(0, 1, data.shape[1])

    # -------------------------------------------------
    # Convert mean/stdev to polar angles in [0, 2π]:
    #   -90° → 0°
    #   +90° → 180°
    #   then map [0°, 180°] to full polar circle by doubling
    # -------------------------------------------------
    mean_deg = np.rad2deg(mean_rad)
    theta = np.deg2rad((mean_deg + 90) * 2)

    data2pi = np.deg2rad((data + 90) * 2)

    upper = np.deg2rad((mean_deg + np.rad2deg(std_rad) + 90) * 2)
    lower = np.deg2rad((mean_deg - np.rad2deg(std_rad) + 90) * 2)

    # ===========================
    #       POLAR PLOT
    # ===========================
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))

    # Tick labels for -90 → +90
    ticks_360 = np.arange(0, 360, 30)
    ticks_labels = [f"{int((t/2)-90)}°" for t in ticks_360]

    ax.set_thetagrids(ticks_360, labels=ticks_labels)

    # Mean
    ax.plot(theta, r, "k", label="Mean Angle")
    for p in data2pi:
        ax.scatter(p, r, color='g', marker="x", alpha= 0.2)

    # Std band
    #ax.fill_betweenx(r, lower, upper, alpha=0.3, color="gray")

    ax.set_title("Polar Projection of Helical Angles", fontsize=14)
    ax.set_ylim(0, 1)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    fig.tight_layout()
    if savepath != None:
        plt.savefig(savepath + '_Ptransmural_lateral.pdf', dpi=1000)
    plt.show()


    # ===========================
    #       LINEAR PLOT
    # ===========================
    mean_deg = np.rad2deg(mean_rad)
    std_deg = np.rad2deg(std_rad)

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.plot(r, mean_deg, "k-x", label="Mean")
    ax.fill_between(r, mean_deg - std_deg, mean_deg + std_deg,
                    color="gray", alpha=0.4, label="±1 SD")

    ax.set_xlabel("Transmural depth (0 = endo, 1 = epi)")
    ax.set_ylabel("Helical angle (°)")
    ax.set_title("Transmural Helical Angle (Linear Plot)")
    ax.grid(True)
    ax.set_ylim(-90,90)
    ax.legend()

    fig.tight_layout()
    if savepath != None:
        plt.savefig(savepath + '_Ltransmural_lateral.pdf', dpi=1000)
    plt.show()
    
    mask = ~np.isnan(mean_deg)

    r = r[mask]
    mean_deg = mean_deg[mask]
    
    return r, mean_deg




# Models
def flat_model(x, c):
    return c * np.ones_like(x)

def linear_model(x, m, b):
    return m*x + b

def find_best_breakpoint(x, y, min_size=20):
    """
    Finds breakpoint where data transitions from flat to downward slope.
    min_size = minimum number of points allowed in each segment
    """
    n = len(x)
    best_sse = np.inf
    best_idx = None
    best_params = None

    # try each possible breakpoint
    for i in range(min_size, n - min_size):
        # segment 1: flat
        x1, y1 = x[:i], y[:i]
        popt1, _ = curve_fit(flat_model, x1, y1)
        y1_pred = flat_model(x1, *popt1)
        sse1 = np.sum((y1 - y1_pred)**2)

        # segment 2: linear
        x2, y2 = x[i:], y[i:]
        popt2, _ = curve_fit(linear_model, x2, y2)
        y2_pred = linear_model(x2, *popt2)
        sse2 = np.sum((y2 - y2_pred)**2)

        total_sse = sse1 + sse2

        # keep best candidate
        if total_sse < best_sse:
            best_sse = total_sse
            best_idx = i
            best_params = (popt1, popt2)

    return best_idx, best_params






def wrap_angles(y):
    """
    Wrap angles (in degrees) to [-180, 180)
    so that 90° == -90° in distance calculations.
    """
    return np.degrees(np.angle(np.exp(1j * np.radians(y))))





def fast_linear_fit(x, y):
    """Analytical least squares fit for y = mx + b."""
    n = len(x)
    sx = np.sum(x)
    sy = np.sum(y)
    sxx = np.sum(x*x)
    sxy = np.sum(x*y)

    denom = n*sxx - sx*sx
    if denom == 0:
        return 0.0, np.mean(y)

    m = (n*sxy - sx*sy) / denom
    b = (sy - m*sx) / n
    return m, b






############### Introduce constrained wrapping 



def constrained_unwrap(
    p,
    period=np.pi,
    max_range=np.deg2rad(350),
    step=0.05,
    max_frac=0.95
):
    p = np.asarray(p)

    mask = ~np.isnan(p)
    p_valid = p[mask]

    if len(p_valid) == 0:
        return p.copy()

    frac = 0.5

    while frac <= max_frac:
        thresh = frac * period
        out = np.copy(p_valid)

        # incremental unwrap
        for i in range(1, len(out)):
            dp = out[i] - out[i - 1]

            if dp > thresh:
                out[i] -= period
            elif dp < -thresh:
                out[i] += period

        # range check
        if (np.nanmax(out) - np.nanmin(out)) <= max_range:
            result = np.full_like(p, np.nan)
            result[mask] = out
            return result

        frac += step

    # fallback
    result = np.full_like(p, np.nan)
    result[mask] = out
    return result

def wrap_pi_deg(y):
    """
    Wrap angles (degrees) to [-90, 90)
    consistent with π-periodicity.
    """
    return (y + 90) % 180 - 90


def find_two_breakpoints(
    x, y,
    min_slope_size=100,
    w_slope=1.5,
    w_flat=0.1,
    penalty_strength=0.2
):

    n = len(x)
    max_tail = int(0.5 * n)  # 50% limit for 3rd flat section

    best_loss = np.inf
    best_breaks = None
    best_params = None

    i_candidates = np.linspace(1, n-2, 100, dtype=int)

    for i in i_candidates:

        valid_j_range = np.arange(i + min_slope_size, n)
        if len(valid_j_range) == 0:
            continue

        j_candidates = np.linspace(valid_j_range[0], valid_j_range[-1], 100, dtype=int)

        for j in j_candidates:

            # enforce max size of segment 3
            if (n - j) > max_tail:
                continue

            #### SEGMENT 1 (flat)
            y1 = y[:i]
            if len(y1):
                c1 = np.degrees(np.angle(np.mean(np.exp(1j*np.radians(y1)))))
                SSE1 = np.sum((wrap_pi_deg(y1 - c1))**2)
                SST1 = np.sum((wrap_pi_deg(y1 - np.mean(y1)))**2)
            else:
                SSE1 = 0
                SST1 = 0
                c1 = 0

            #### SEGMENT 2 (slope)
            x2 = x[i:j]
            y2 = y[i:j]

            if len(x2) < 2:
                continue

            # unwrap angles for linear fit
            #y2_unwrapped = np.rad2deg(constrained_unwrap(np.deg2rad(y2)))

            y2_unwrapped = y2
            
            m, b = fast_linear_fit(x2, y2_unwrapped)

            # enforce negative slope
            if m > 0:
                m = 0

            y2_pred = m*x2 + b

            SSE2 = np.sum((y2_unwrapped - y2_pred)**2)
            SST2 = np.sum((y2_unwrapped - np.mean(y2_unwrapped))**2)

            r2 = 1 - SSE2/SST2 if SST2 > 0 else -np.inf

            # refit unconstrained if fit is terrible
            if r2 < 0:
                m, b = fast_linear_fit(x2, y2_unwrapped)
                y2_pred = m*x2 + b
                SSE2 = np.sum((y2_unwrapped - y2_pred)**2)

            #### SEGMENT 3 (flat)
            y3 = y[j:]
            if len(y3):
                c3 = np.degrees(np.angle(np.mean(np.exp(1j*np.radians(y3)))))
                SSE3 = np.sum((wrap_pi_deg(y3 - c3))**2)
                SST3 = np.sum((wrap_pi_deg(y3 - np.mean(y3)))**2)
            else:
                SSE3 = 0
                SST3 = 0
                c3 = 0

            #### slope length penalty

            
            R2_2 = 1 - SSE2 / SST2
            
            R2_1 = 1 - SSE1 / SST1 if SST1 > 0 else 1
            R2_3 = 1 - SSE3 / SST3 if SST3 > 0 else 1
            
            L = j-i
            penalty = penalty_strength * L * (n / L)**2
            
            L1 = len(y1)
            L2 = len(y2)
            L3 = len(y3)
            
            loss = (
                w_slope * L2 * (1 - R2_2) +
                w_flat * (L1 * (1 - R2_1) + L3 * (1 - R2_3)) +
                penalty
            )
            
            print(f'E Slope: {w_slope * L2 * (1 - R2_2)}, E Flat: {w_flat * (L1 * (1 - R2_1) + L3 * (1 - R2_3))}, Penalty: {penalty},')

            if loss < best_loss:
                best_loss = loss
                best_breaks = (i, j)
                best_params = ((c1,), (m, b), (c3,))

    return best_breaks, best_params


def process_transmural_profiles(
    profiles,
    target_len=None,
    plot=False
):
    
    """
    Resample, unwrap, fit flat→slope→flat model, and circularly average profiles.
    Returns transmural coordinate, circular mean/std, and R² values for slope fits.
    """
    
    max_len = max(len(p) for p in profiles)
    padded_profiles = [np.pad(p, (max_len - len(p), 0), mode="constant", constant_values=np.nan) for p in profiles]
    
    pp_rad = np.deg2rad(padded_profiles)
    
    # --- circular statistics ---
    mean_rad = circmean(
        pp_rad,
        axis=0,
        low=-np.pi/2,
        high=np.pi/2,
        nan_policy="omit"
    )

    std_rad = circstd(
        pp_rad,
        axis=0,
        low=-np.pi/2,
        high=np.pi/2,
        nan_policy="omit"
    )

    if target_len != None:
        
        x = np.linspace(0, 1, target_len)
        x_orig = np.linspace(0, 1, len(mean_rad))
        
        mask_o = ~np.isnan(mean_rad) #interp can't handle nan
        mask_s = ~np.isnan(std_rad)

        p_unwrapped = constrained_unwrap(mean_rad[mask_o], period=np.pi)
        
        mask_u = ~np.isnan(p_unwrapped)
        
        p_resampled = np.interp(x, x_orig[mask_o], p_unwrapped[mask_u])
        st_resampled = np.interp(x, x_orig[mask_o], std_rad[mask_s])

        rewrapped = (p_resampled + np.pi/2) % (np.pi) - np.pi/2 # Rewrap to [-π/2, +π/2]

        mean_deg = np.rad2deg(rewrapped)
        std_deg = np.rad2deg(st_resampled)
    
        unwrapped_mean = np.rad2deg(p_resampled)
    
    else: 
        x = np.linspace(0, 1, len(mean_rad))
        
        mean_deg = np.rad2deg(mean_rad)
        std_deg = np.rad2deg(std_rad)
        
        unwrapped_mean = np.rad2deg(constrained_unwrap(np.deg2rad(mean_deg), period=np.pi))
    
    mask = ~np.isnan(unwrapped_mean)
    
    y = unwrapped_mean
    
    y_mid = y[int(0.4*len(y)):int(0.8*len(y))]
    
    while y_mid.mean() < -90 or y_mid.mean() > 90:
        y_mid = y[int(0.4*len(y)):int(0.8*len(y))]
        if y_mid.mean() < -90:
            y += 180
        elif y_mid.mean() > 90:
            y -= 180

    # --- breakpoint fitting ---

    (b1, b2), (flat1_params, slope_params, flat2_params) = find_two_breakpoints(
        x, 
        y[mask],
    )
    
    mask = mask = ~np.isnan(mean_deg)
    
    if plot != False: 
                
        c1 = flat1_params[0]
        m, b = slope_params
        c3 = flat2_params[0]
        
        # Generate fitted segments
        y_fit1 = flat_model(x[:b1], c1)
        y_fit2 = linear_model(x[b1:b2], m, b)
        y_fit3 = flat_model(x[b2:], c3)
        
        r2 = r2_score(y[b1:b2], y_fit2)
        
        plt.figure(figsize=(12,6))
        
        plt.plot(x, y, 'k.', markersize=3)
        
        plt.plot(x[:b1], y_fit1, 'r-', linewidth=2)
        plt.plot(x[b1:b2], y_fit2, 'b-', linewidth=2, label=f"Slope gradient {m:.3f}")
        plt.plot(x[b2:], y_fit3, 'g-', linewidth=2)
        
        plt.axvline(x[b1], color='r', linestyle='--', label=f"Break 1 = {x[b1]:.3f}")
        plt.axvline(x[b2], color='g', linestyle='--', label=f"Break 2 = {x[b2]:.3f}")
        
        plt.xlabel("Transmural depth")
        plt.ylabel("Helical angle")
        plt.title(f"{plot}, R2 = {r2}")
        plt.legend()
        plt.grid(True)
        plt.show()

    return (
        x[mask],
        mean_deg[mask],
        std_deg[mask],
        (b1, b2), 
        (flat1_params, slope_params, flat2_params)
    )


def inspect_segment_fraction(
    heart_id,
    segment,
    df_Mega,
    period=np.pi
):

    # --- load data ---
    mask_global = (
        (df_Mega["heart"] == heart_id) &
        (df_Mega["segment"] == segment)
    )

    df = df_Mega[mask_global].copy()
    if df.empty:
        print("No data found")
        return

    df = df.sort_values("depth").reset_index()

    t_depth = df["depth"].to_numpy()
    y = df["HA"].to_numpy()


    # --- fractional init ---
    b1_frac = df["s_compact_start"].iloc[0]
    b2_frac = df["s_compact_end"].iloc[0]   
    
    n = len(t_depth)
    
    shift_mask = np.zeros(n, dtype=bool)
    shift_sign = 0
    

    while True:
        
        y_plot = y.copy()

        n = len(t_depth)
        b1 = int(b1_frac * n)
        b2 = int(b2_frac * n)


        # apply shift safely (mapped to local df index)
        if shift_sign != 0:
            local_mask = np.zeros(n, dtype=bool)
            local_mask[:] = shift_mask[df.index]
            y_plot[local_mask] += shift_sign * 180

        # --- fit ---
        c1 = np.mean(y_plot[:b1])
        c3 = np.mean(y_plot[b2:])

        m, c = fast_linear_fit(t_depth[b1:b2], y_plot[b1:b2])

        y_fit = linear_model(t_depth[b1:b2], m, c)

        r2 = r2_score(y_plot[b1:b2], y_fit)  

        # --- plot ---
        plt.figure(figsize=(12, 6))

        plt.plot(t_depth, y_plot, 'k.', markersize=3)
        plt.plot(t_depth[:b1], np.full(b1, c1), 'r-')
        plt.plot(t_depth[b1:b2], y_fit, 'b-')
        plt.plot(t_depth[b2:], np.full(len(t_depth)-b2, c3), 'g-')

        plt.axvline(t_depth[b1], color='r', linestyle='--')
        plt.axvline(t_depth[b2], color='g', linestyle='--')

        plt.title(f"Heart {heart_id}, Segment {segment}, R2={r2:.3f}")
        plt.xlabel("Depth")
        plt.ylabel("HA")
        plt.minorticks_on()
        plt.grid(True, which='both')
        plt.show()

        cmd = input("Command (b1/b2/shift+/shift-/clear/commit/exit): ")

        if cmd == "exit":
            break

        elif cmd == "b1":
            b1_frac = float(input("b1 fraction (0-1): "))

        elif cmd == "b2":
            b2_frac = float(input("b2 fraction (0-1): "))

        elif cmd == "shift+":
            a, b = map(float, input("a:b (0-1): ").split(":"))
            shift_mask[:] = False
            idx = slice(int(a*n), int(b*n))
            shift_mask[df.index[idx]] = True
            shift_sign = +1

        elif cmd == "shift-":
            a, b = map(float, input("a:b (0-1): ").split(":"))
            shift_mask[:] = False
            idx = slice(int(a*n), int(b*n))
            shift_mask[df.index[idx]] = True
            shift_sign = -1

        elif cmd == "clear":
            shift_mask[:] = False
            shift_sign = 0

        elif cmd == "commit":

            # --- recompute final y ---
            y_commit = y.copy()
        
            if shift_sign != 0:
                y_commit[shift_mask] += shift_sign * 180
        
            # --- write back using original indices ---
            idx = df["index"].values  # original df_Mega indices
        
            df_Mega.loc[idx, "HA"] = y_commit
        
            df_Mega.loc[idx, "s_compact_start"] = b1_frac
            df_Mega.loc[idx, "s_compact_end"] = b2_frac
            df_Mega.loc[idx, "s_m"] = m
            df_Mega.loc[idx, "s_c"] = c
            df_Mega.loc[idx, "s_r2"] = r2
        
            print(f"Committed heart {heart_id}, segment {segment}")
                    

            
################################ Histograms


def histogram(folder=None, files=None, start=0, end=None, mask=None, seg_map=None, factor=1, rdr=None):
    """
    Compute histograms either from a folder of image files or from a DataReader `rdr`.

    Parameters
    ----------
    folder : str or None
        Path to folder with image slices (used if `rdr` is None).
    files : list or None
        Precomputed list of filenames (optional).
    start, end : int
        Slice range to process.
    mask : ndarray or None
        Binary mask volume (indexed by slice) to apply.
    seg_map : ndarray or None
        Segment map volume (indexed by slice) to compute per-segment histograms.
    factor : int
        Scale factor when mapping mask/seg_map indices to image indices.
    rdr : DataReader or None
        If provided, stream slices from this reader instead of reading files from disk.

    Returns
    -------
    bins : ndarray
        Histogram counts. If `seg_map` provided returns shape (n_segments, 256),
        else returns shape (256,).
    unique_segments : ndarray or None
        Array of segment labels when `seg_map` provided, else None.
    """

    unique_segments = None

    # Determine slice loader
    def load_slice(idx):
        if rdr is not None:
            return rdr.load_volume(start_index=idx, end_index=idx + 1)[0]
        else:
            return imread(os.path.join(folder, files[idx]))

    # If using folder, build file list
    if rdr is None:
        if files is None:
            files = sorted([f for f in os.listdir(folder) if f.lower().endswith(('.tif', '.tiff', '.jp2'))])
        if end is None:
            end = len(files)
    else:
        # rdr provided
        Z = rdr.shape[0]
        if end is None:
            end = Z

    # Branch by mask/segment/none
    if mask is not None:
        bins = np.zeros(256, dtype=np.int64)
        for i in range(start, end):
            im = load_slice(i).astype(np.uint8)
            im_mask = resample(mask[int(i // factor)], target_shape=im.shape, TwoD=True)
            vals = im[im_mask > 0]
            vals = vals[~np.isnan(vals)]
            if vals.size:
                bins += np.bincount(vals, minlength=256)
    elif seg_map is not None:
        unique_segments = np.unique(seg_map)
        unique_segments = unique_segments[unique_segments != 0]
        bins = np.zeros((len(unique_segments), 256), dtype=np.int64)
        for i in range(start, end):
            im = load_slice(i).astype(np.uint8)
            im_mask = resample(seg_map[int(i // factor)], target_shape=im.shape, TwoD=True)
            for si, s in enumerate(unique_segments):
                mask_loc = (im_mask == s)
                if not mask_loc.any():
                    continue
                vals = im[mask_loc]
                vals = vals[~np.isnan(vals)]
                if vals.size:
                    bins[si] += np.bincount(vals, minlength=256)
    else:
        bins = np.zeros(256, dtype=np.int64)
        for i in range(start, end):
            im = load_slice(i).astype(np.uint8)
            vals = im[~np.isnan(im)]
            if vals.size:
                bins += np.bincount(vals, minlength=256)

    return bins, unique_segments

  

def plot_segment_histogram(
    bins,
    segments,
    value_range=(-90, 90),
    smooth_sigma=0,
    normalize=False,
    show_mean=False,
    circular_mean=False,
    xlab = 'Value',
    ylab = 'Count',
    title=None, 
    divide_radial=False,
    headless= False
):

    # --------------------
    # Combine selected segments
    # --------------------
    if divide_radial:
        ##divide radial must be true
        seg_idx = []
        for s in segments:
            base = s // 10
            r = s % 10
            idx = (base - 1) * 4 + r -1
            seg_idx.append(idx)

    else:
        seg_idx = [s - 1 for s in segments]
    
    hist = bins[:, seg_idx].sum(axis=1)
    
    # Crop histogram to negate masking artifacts
    hist[0] = hist[1]
    hist[-1] = hist[-2]

    # --------------------
    # X axis mapping
    # --------------------
    x = np.linspace(value_range[0], value_range[1], hist.size)

    # --------------------
    # Normalize
    # --------------------
    if normalize and hist.sum() > 0:
        hist = hist / hist.sum()

    # --------------------
    # Smoothing
    # --------------------
    if smooth_sigma > 0:
        hist = gaussian_filter1d(hist, sigma=smooth_sigma)

    # --------------------
    # Plot
    # --------------------
    fig, ax = plt.subplots()
    ax.plot(x, hist, linewidth=2)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    if title:
        ax.set_title(title)

    # --------------------
    # Mean line
    # --------------------
    if show_mean:
        if circular_mean:
            angles = np.linspace(value_range[0], value_range[1], hist.size, endpoint=False)
            theta = np.deg2rad(2*angles)
            
            z = np.sum(hist * np.exp(1j*theta))
            cmean = np.rad2deg(np.angle(z))/2
            
            # wrap into [-90,90)
            if cmean < -90:
                cmean += 180
            elif cmean >= 90:
                cmean -= 180

            ax.axvline(cmean, linestyle='--', linewidth=2, label=f"Circlar Mean = {cmean:.2f}")
            mean_val = np.average(x, weights=hist)
            ax.axvline(mean_val, linestyle='-', linewidth=2, label=f"Mean = {mean_val:.2f}")
        else:
            mean_val = np.average(x, weights=hist)
            ax.axvline(mean_val, linestyle='-', linewidth=2, label=f"Mean = {mean_val:.2f}")
        ax.legend()
    
    if not headless: 
        plt.show()
    
    return fig 

            
def plot_segment_histogram_labels(
    bins,
    segments,
    value_range=(-90, 90),
    smooth_sigma=0,
    normalize=False,
    show_mean=False,
    circular_mean=False,
    xlab='Value',
    ylab='Count',
    title=None,
    ax=None   # <-- add this
):

    # --------------------
    # Combine selected segments
    # --------------------
    rows = bins[np.isin(bins[:, 0], segments)]
    hist = rows[:, 1:].sum(axis=0)

    # Crop histogram to negate masking artifacts
    hist[0] = hist[1]
    hist[-1] = hist[-2]

    # --------------------
    # X axis mapping
    # --------------------
    x = np.linspace(value_range[0], value_range[1], hist.size)

    # --------------------
    # Normalize
    # --------------------
    if normalize and hist.sum() > 0:
        hist = hist / hist.sum()

    # --------------------
    # Smoothing
    # --------------------
    if smooth_sigma > 0:
        hist = gaussian_filter1d(hist, sigma=smooth_sigma)

    # --------------------
    # Plot (use provided axis)
    # --------------------
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    ax.plot(x, hist, linewidth=2)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    if title:
        ax.set_title(title)

    # --------------------
    # Mean line
    # --------------------
    if show_mean and hist.sum() > 0:
        if circular_mean:
            scale_rad = np.linspace(0, 2*np.pi, 256)
            scale_mean_val = np.rad2deg(stats.circmean(scale_rad, weights=hist))
            f = (value_range[1] - value_range[0]) / 360
            cmean_val = scale_mean_val * f + value_range[0]
            if cmean_val < value_range[0]:
                cmean_val += f * 360

            ax.axvline(cmean_val, linestyle='--', linewidth=2, label=f"Circ Mean = {cmean_val:.2f}")

        mean_val = np.average(x, weights=hist)
        ax.axvline(mean_val, linestyle='-', linewidth=2, label=f"Mean = {mean_val:.2f}")
        ax.legend()

    return fig
            
            
def segment_means_from_histogram(
    bins,
    value_range=(-90, 90),
    circular=False
):
    """
    bins : ndarray (256, no_segs)
        Histogram counts per segment

    value_range : tuple
        Mapping of bin index -> physical values

    circular : bool
        Use circular mean (for angles)

    Returns
    -------
    means : dict
    """
    means = {}

    # Physical x-axis
    x = np.linspace(value_range[0], value_range[1], bins.shape[0])

    for s in range(bins.shape[1]):
        hist = bins[:,s]
        
        # Crop histogram to negate masking artifacts
        hist[0] = hist[1]
        hist[-1] = hist[-2]

        if hist.sum() == 0:
            continue

        if circular:
            scale_rad = np.linspace(0, 2*np.pi, len(hist), endpoint=False)
            scale_mean_val = np.rad2deg(stats.circmean(scale_rad, weights=hist))
            f = (value_range[1] - value_range[0])  / 360
            cmean_val = scale_mean_val*f + value_range[0]
            if cmean_val < value_range[0]:
                cmean_val += f*360
            means[s+1] = cmean_val
        else:
            means[s+1] = np.average(x, weights=hist)

    return means

