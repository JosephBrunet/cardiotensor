
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import read_conf_file
import cardiotensor.analysis.analysis_utils as utils
import csv
from cardiotensor.colormaps.helix_angle import helix_angle_cmap


def script() -> None:
    ap = argparse.ArgumentParser(description="AHA-segmented histograms")
    ap.add_argument("input", type=Path, help=".conf file or directory of images to sample")
    ap.add_argument("--sample", type=str, default="HA", help="Data for transmural analysis (e.g., 'HA', 'FA', 'IA)")
    ap.add_argument("--lv-mask", type=Path, required=True, help="LV mask directory (slices)")
    ap.add_argument("--axis-points", type=str, default=None, help="Override axis points as 'z1,y1,x1;z2,y2,x2'")
    ap.add_argument("--septum", type=str, default=None, help="Septum point 'z,y,x' (optional interactive if omitted)")
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--write_csv", type=bool, default=True, help="If true record intensity profiles in a CSV file in the output directory")
    ap.add_argument("--plot", type=bool, default=True, help="Plot mean_intensity profile")
    ap.add_argument("--mask_downsample", type=int, default=1, help="Downsampling factor for the LV mask to speed up processing (default: 1, no downsampling)")
    ap.add_argument("--nslices", type=int, default=5, help="Number of slices to sample")
    ap.add_argument("--slice_idxs", type=int, default=None, help="Specify slice indices for sampling. If not provided, slices will be evenly spaced across the LV mask.")
    args = ap.parse_args()

    # Resolve config/base
    if args.input.suffix.lower() == ".conf":
        params = read_conf_file(args.input)
        # When a .conf is provided, image data (histograms) live under OUTPUT_PATH/<sample>
        output_base = Path(params.get("OUTPUT_PATH"))
        images_path = str(output_base / args.sample)
        cfg_axis = params.get("AXIS_POINTS", [])
        # Place output in existing file structure
        output_base = Path(params.get("OUTPUT_PATH"))
        default_out = output_base / "Analysis"
    elif args.input.is_dir():
        images_path = str(args.input)
        cfg_axis = []
        default_out = Path("./transmral_analysis")
    else:
        raise ValueError("Input must be a .conf or a directory")

    ## Use requested output directory or default
    outdir = args.outdir if args.outdir is not None else default_out
    outdir.mkdir(parents=True, exist_ok=True)

    # HA directory discovery
    ha_dir = Path(images_path)  # fallback
    ha_rdr = DataReader(ha_dir)

    # Axis points
    if args.axis_points:
        try:
            p1s, p2s = args.axis_points.split(";")
            axis_points = [tuple(map(float, p1s.split(","))), tuple(map(float, p2s.split(",")))]
        except Exception:
            raise ValueError("--axis-points must be 'z1,y1,x1;z2,y2,x2'")
    else:
        axis_points = cfg_axis if cfg_axis else None
    
    # Septum
    if args.septum:
        sept = tuple(map(float, args.septum.split(",")))
    else:
        # interactive picker (will open a matplotlib window)
        sept = utils.get_septum_point(str(args.lv_mask))

    if axis_points is None or len(axis_points) < 2:
        raise ValueError("Axis points must be provided either in .conf or via --axis-points")

    print("Creating AHA mask")
    seg_map, _ = utils.create_aha_mask(
        str(args.lv_mask), 
        axis_points, 
        sept, 
        str(images_path), 
        return_radial_map=False, 
        additional_downsample=args.mask_downsample,
    )
    
    # Resample segment map to HA reader shape
    target_shape = tuple(ha_rdr.shape)
    seg_resampled = utils.resample(seg_map, target_shape=target_shape, order=0, TwoD=False, is_mask=False)

    axis_vec = np.asarray(params['AXIS_POINTS'][0], dtype=float) - np.asarray(params['AXIS_POINTS'][1], dtype=float)
    axis_vec = axis_vec / np.linalg.norm(axis_vec)

    mask_factor = np.round(seg_resampled.shape[0] / seg_map.shape[0])

    if not slice_idxs: 
        slice_idxs = np.linspace(params['AXIS_POINTS'][0][2], params['AXIS_POINTS'][1][2], 12)[1:-1]
    
    ## Search for existing reslices 
    search_dir = outdir / "Analysis" / "Reslices"
    search_dir.mkdir(parents=True, exist_ok=True)

    # Ensure slice_idxs is a list
    if np.isscalar(slice_idxs):
        slice_idxs = [slice_idxs]
    slice_idxs = list(slice_idxs)

    # Look for files named like reslice_<idx>* in search_dir. Load any found
    # and build `slices` as the list of indices missing from search_dir.
    slices = []  # indices to gather (not present on disk)
    reslices = {}  # idx -> loaded array
    for idx in slice_idxs:
        # Format index for filename matching: prefer integer if close
        try:
            fidx = int(idx) if float(idx).is_integer() else str(idx)
        except Exception:
            fidx = str(idx)

        matches = list(search_dir.glob(f"reslice_{fidx}*"))
        if matches:
            p = matches[0]
            try:
                if p.suffix.lower() == ".npy":
                    reslices[idx] = np.load(p)
                else:
                    reslices[idx] = plt.imread(p)
            except Exception:
                try:
                    reslices[idx] = np.load(p)
                except Exception:
                    reslices[idx] = None
        else:
            slices.append(idx)

    
    C_news = []
    C_masks = []
    C_reslice = []
    C_done = []
    xm = seg_resampled.shape[2] //2
    ym = seg_resampled.shape[1] //2
    
    for i in slice_idxs:
        x, y = utils.line_xy_at_z(params['AXIS_POINTS'][0], params['AXIS_POINTS'][1], i)
    
        C = np.array([i, ym, xm], dtype=float)
        P = np.array([i, y, x], dtype=float)
    
        t = np.dot(P - C, axis_vec)
        C_new = C + t * axis_vec
        
        C_news.append(C_new)
        C_masks.append(C_new/mask_factor)


        if i in slices:
            C_reslice.append(C_new)
        else:
            C_done.append(C_new)

    new_slices, geom = utils.ResliceMultipleSlices(axis_vec, C_reslice, images_path)


    ########################### To do: Match up geometry 
    for C_done in zip(C_done):
        geo = utils.GetResliceGeometry(axis_vec, C_reslice, images_path)

    for idx, slice in zip(slices, new_slices):
        reslices[idx] = slice
    
    

    centroids = []
    for s, geom in zip(slice_idxs, geom):
        centroid_slice = utils.line_xy_at_z(params['AXIS_POINTS'][0], params['AXIS_POINTS'][1], s)
        centroids.append(utils.project_3d_to_slice([s, centroid_slice[1], centroid_slice[0]], geom))

    t_profiles = {}
    t_profiles['PROFILE'], t_profiles['LINES'], t_profiles['SEGMENTS'], t_profiles['INDICES'] = utils.transmuralSampling(C_news, axis_vec, slices=list(reslices.values()), centroids=centroids, mask_factor=mask_factor, seg_map=seg_map, route_slices=slice_idxs)
    
        