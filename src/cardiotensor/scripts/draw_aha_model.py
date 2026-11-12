#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone drawing/test script for an AHA-style 17- or 26-segment model.

Run:
    python draw_aha26_model.py

Import:
    from draw_aha26_model import draw_aha_model

The geometry constants near the top are intentionally easy to tweak while you
match the RV shape to a reference image.
"""


import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib.path import Path as MplPath
from matplotlib.patches import Circle, PathPatch, Wedge
import cardiotensor.analysis.analysis_utils as utils

import pandas as pd


# =============================================================================
# TWEAK THESE TO MATCH YOUR REFERENCE
# =============================================================================
LV_CENTER = (0.55, 0.0)
LV_OUTER_RADIUS = 4.0

# Segment order starts at the top and goes counter-clockwise, matching the
# common AHA bullseye layout: 1 top, 2 upper-left, 3 lower-left, 4 bottom, etc.
LV_RING_DEFS = [
    (range(1, 7), 3.0, 4.0, np.linspace(60, 420, 7)),
    (range(7, 13), 2.0, 3.0, np.linspace(60, 420, 7)),
    (range(13, 17), 1.0, 2.0, np.linspace(45, 405, 5)),
]

# Angle (degrees from the +x axis) at which the outermost RV free wall attaches to the LV.
RV_INSERTION_ANGLE = 90.0

# How far the outermost RV free wall extends to the LEFT (horizontal semi-axis).
RV_OUTER_RX = 7

# Four concentric half-ellipse boundaries, outermost → innermost.
# Each entry: (insertion_angle_deg, rx)
#   insertion_angle_deg : where the ellipse tips touch the LV outer circle
#   rx                  : horizontal semi-axis (width of the half-ellipse)
# Decrease insertion_angle and rx for each successive inner boundary so the
# nested half-ellipses have smaller, lower insertion points — matching the picture.
RV_HALF_ELLIPSES = [
    (90.0, 6.5),   # free-wall outer boundary  (outer edge of band 18-20)
    (60.0, 7.4),   # divider between bands 18-20 and 21-23
    (50.0, 7.2),   # divider between bands 21-23 and 24-26
    (50.0, 4),   # septal inner boundary      (inner edge of band 24-26)
]

# Per-band row separators: each entry is (segments, outer_lvl, inner_lvl, (angle1, angle2))
# angle1/angle2 are radial directions from LV center (degrees from +x) that cut the band
# into top / middle / bottom rows.  Tune each band independently.
RV_SEPARATOR_SAMPLES = 2400
RV_BAND_DEFS = [
    ((18, 19, 20), 0, 1, (162.0, 198.0)),   # outer band
    ((21, 22, 23), 1, 2, (167.0, 193.0)),   # middle band
    ((24, 25, 26), 2, 3, (173.0, 187.0)),   # inner band
]
RV_LABEL_LV_CLEARANCE = 0.10

FIGSIZE = (9.0, 6.5)
DEFAULT_MODEL = 26
OUTPUT_PATH = None
SHOW_FIGURE = True
# =============================================================================


DEFAULT_SEGMENT_COLORS = {
    1: "#4154b3",
    2: "#ed4b52",
    3: "#ef6a45",
    4: "#233b9f",
    5: "#112f91",
    6: "#1b3c99",
    7: "#1e8d88",
    8: "#223fb3",
    9: "#4bc4d1",
    10: "#28d84e",
    11: "#b5f0b7",
    12: "#1e8d88",
    13: "#d6e53a",
    14: "#e8df2f",
    15: "#b9df32",
    16: "#c7e83e",
    17: "#f2de36",
    18: "#f16f2a",
    19: "#f7a53a",
    20: "#f7ce82",
    21: "#aa541b",
    22: "#bf7623",
    23: "#d79b43",
    24: "#6d3d18",
    25: "#83511e",
    26: "#a56828",
}



def _rv_insertion_xy() -> tuple[float, float]:
    """Top RV-LV insertion point (x, y) relative to LV_CENTER. Bottom is (x, -y)."""
    theta = np.deg2rad(RV_INSERTION_ANGLE)
    return LV_OUTER_RADIUS * np.cos(theta), LV_OUTER_RADIUS * np.sin(theta)


def _rv_ellipse_curve(t: np.ndarray | float, angle_deg: float, rx: float) -> np.ndarray:
    """Left half-ellipse whose tips lie on the LV outer circle at ±angle_deg from +x.
    t=0 → top tip, t=0.5 → leftmost point, t=1 → bottom tip."""
    t = np.asarray(t, dtype=float)
    angle = np.deg2rad(angle_deg)
    cx = LV_OUTER_RADIUS * np.cos(angle)   # x of the ellipse centre
    ry = LV_OUTER_RADIUS * np.sin(angle)   # vertical semi-axis = y of top tip
    phi = np.pi / 2.0 + t * np.pi          # sweeps π/2 → π → 3π/2
    x = cx + rx * np.cos(phi)
    y = ry * np.sin(phi)
    return np.column_stack((np.atleast_1d(x), np.atleast_1d(y)))


def _lv_arc(from_deg: float, to_deg: float, n: int = 18) -> np.ndarray:
    """Arc of the LV outer circle from from_deg to to_deg (degrees from +x axis)."""
    thetas = np.linspace(np.deg2rad(from_deg), np.deg2rad(to_deg), n)
    x = LV_OUTER_RADIUS * np.cos(thetas)
    y = LV_OUTER_RADIUS * np.sin(thetas)
    return np.column_stack((x, y))


def _rv_separator_t(
    ellipse_level: int,
    separator_angle_deg: float,
) -> float:
    """Find where a straight RV row separator intersects one half-ellipse."""
    angle_deg, rx = RV_HALF_ELLIPSES[ellipse_level]
    t_values = np.linspace(0.0, 1.0, RV_SEPARATOR_SAMPLES)
    points = _rv_ellipse_curve(t_values, angle_deg, rx)
    theta = np.deg2rad(separator_angle_deg)
    direction = np.array([np.cos(theta), np.sin(theta)])
    cross = np.abs(direction[0] * points[:, 1] - direction[1] * points[:, 0])
    along = points @ direction
    cross = np.where(along > 0.0, cross, np.inf)
    return float(t_values[np.argmin(cross)])


def _rv_separator_t_pair(
    outer_lvl: int,
    inner_lvl: int,
    separator_angle_deg: float,
) -> tuple[float, float]:
    """Return (t_out, t_in): the separator hits the outer ellipse at the point
    closest to separator_angle_deg from the LV center, then follows the inward
    normal of the outer ellipse to find where it crosses the inner ellipse."""
    t_out = _rv_separator_t(outer_lvl, separator_angle_deg)

    angle_out, rx_out = RV_HALF_ELLIPSES[outer_lvl]
    angle_in,  rx_in  = RV_HALF_ELLIPSES[inner_lvl]
    angle_out_rad = np.deg2rad(angle_out)
    angle_in_rad  = np.deg2rad(angle_in)
    cx_out = LV_OUTER_RADIUS * np.cos(angle_out_rad)
    ry_out = LV_OUTER_RADIUS * np.sin(angle_out_rad)
    cx_in  = LV_OUTER_RADIUS * np.cos(angle_in_rad)
    ry_in  = LV_OUTER_RADIUS * np.sin(angle_in_rad)

    phi = np.pi / 2.0 + t_out * np.pi
    x0 = cx_out + rx_out * np.cos(phi)
    y0 = ry_out * np.sin(phi)

    # Inward normal to the outer ellipse at (x0, y0): -∇f / |∇f|
    nx_raw = -np.cos(phi) / rx_out
    ny_raw = -np.sin(phi) / ry_out
    mag = np.hypot(nx_raw, ny_raw)
    nx, ny = nx_raw / mag, ny_raw / mag

    # Intersect the normal line with the inner ellipse
    A = (nx / rx_in) ** 2 + (ny / ry_in) ** 2
    B = 2.0 * ((x0 - cx_in) * nx / rx_in ** 2 + y0 * ny / ry_in ** 2)
    C = ((x0 - cx_in) / rx_in) ** 2 + (y0 / ry_in) ** 2 - 1.0
    disc = B ** 2 - 4.0 * A * C
    if disc < 0:
        return t_out, _rv_separator_t(inner_lvl, separator_angle_deg)

    s1 = (-B + np.sqrt(disc)) / (2.0 * A)
    s2 = (-B - np.sqrt(disc)) / (2.0 * A)
    positives = sorted(s for s in (s1, s2) if s > 1e-9)
    if not positives:
        return t_out, _rv_separator_t(inner_lvl, separator_angle_deg)

    s_hit = positives[0]
    x_hit = x0 + s_hit * nx
    y_hit = y0 + s_hit * ny

    # Convert hit point to t on inner ellipse
    cos_phi_in = np.clip((x_hit - cx_in) / rx_in, -1.0, 1.0)
    sin_phi_in = y_hit / ry_in
    phi_in = np.arctan2(sin_phi_in, cos_phi_in)
    if phi_in < np.pi / 2.0 - 1e-9:
        phi_in += 2.0 * np.pi
    t_in = float(np.clip((phi_in - np.pi / 2.0) / np.pi, 0.0, 1.0))
    return t_out, t_in


def _rv_segment_path(
    center: tuple[float, float],
    row1: int,
    row2: int,
    outer_lvl: int,
    inner_lvl: int,
    sep_angles: tuple,
) -> MplPath:
    """Closed path for one RV segment bounded by two nested half-ellipses.
    Row separators are perpendicular to the outer ellipse at the cut point."""
    N = 50
    angle_out, rx_out = RV_HALF_ELLIPSES[outer_lvl]
    angle_in,  rx_in  = RV_HALF_ELLIPSES[inner_lvl]

    if row1 == 0:
        t1_out, t1_in = 0.0, 0.0
    else:
        t1_out, t1_in = _rv_separator_t_pair(outer_lvl, inner_lvl, sep_angles[row1 - 1])

    if row2 == len(sep_angles) + 1:
        t2_out, t2_in = 1.0, 1.0
    else:
        t2_out, t2_in = _rv_separator_t_pair(outer_lvl, inner_lvl, sep_angles[row2 - 1])

    outer_arc = _rv_ellipse_curve(np.linspace(t1_out, t2_out, N), angle_out, rx_out)
    inner_arc = _rv_ellipse_curve(np.linspace(t2_in, t1_in, N), angle_in, rx_in)

    # Bottom edge: at t=1 both tips are on the LV circle -> connect via LV arc
    bottom_conn = (
        _lv_arc(-angle_out, -angle_in)[1:]
        if row2 == len(sep_angles) + 1
        else None
    )
    # Top edge: same at t=0
    top_conn = _lv_arc(angle_in, angle_out)[1:] if row1 == 0 else None

    parts = [outer_arc]
    if bottom_conn is not None:
        parts.append(bottom_conn)
    parts.append(inner_arc)
    if top_conn is not None:
        parts.append(top_conn)

    cx, cy = center
    vertices = np.vstack(parts + [outer_arc[:1]])
    vertices = vertices + np.array([cx, cy])
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(vertices) - 2) + [MplPath.CLOSEPOLY]
    return MplPath(vertices, codes)


def _outside_lv(
    points: np.ndarray,
    center: tuple[float, float],
    clearance: float = 0.0,
) -> np.ndarray:
    """Return boolean mask: True where each point lies outside (or on) the LV outer circle.

    Parameters
    ----------
    points : (N, 2) array of (x, y) coordinates in data space.
    center : (x, y) of the LV center, matching LV_CENTER.
    clearance : extra margin added to LV_OUTER_RADIUS before the distance test.
    """
    center_xy = np.asarray(center, dtype=float)
    distances = np.linalg.norm(points - center_xy, axis=1)
    return distances >= LV_OUTER_RADIUS + clearance


def _rv_label_position(
    path: MplPath,
    center: tuple[float, float],
) -> tuple[float, float]:
    """Find a good (x, y) position for a text label inside an RV segment patch.

    Samples a grid inside the segment bounding box, keeps only points that are
    (a) inside the segment path, (b) on the RV side (left of LV center), and
    (c) outside the LV outer circle.  Returns the grid point closest to the
    centroid of those visible candidates.  Falls back to progressively looser
    criteria, and ultimately to the raw vertex centroid, if no candidate passes.

    Parameters
    ----------
    path   : closed MplPath describing one RV segment (in data coordinates).
    center : (x, y) of the LV center, used to exclude the LV interior.

    Returns
    -------
    (x, y) in data coordinates.
    """
    xmin, ymin = path.vertices.min(axis=0)
    xmax, ymax = path.vertices.max(axis=0)
    x_values = np.linspace(xmin, xmax, 80)
    y_values = np.linspace(ymin, ymax, 80)
    xx, yy = np.meshgrid(x_values, y_values)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    in_segment = path.contains_points(points)
    rv_side = points[:, 0] <= center[0] - 0.05
    visible = in_segment & rv_side & _outside_lv(points, center, RV_LABEL_LV_CLEARANCE)
    if not np.any(visible):
        visible = in_segment & rv_side & _outside_lv(points, center)
    if not np.any(visible):
        visible = in_segment & _outside_lv(points, center)

    if np.any(visible):
        visible_points = points[visible]
        centroid = visible_points.mean(axis=0)
        closest = np.argmin(np.sum((visible_points - centroid) ** 2, axis=1))
        return tuple(visible_points[closest])

    centroid = path.vertices[:-1].mean(axis=0)
    return tuple(centroid)


def _value_color(
    segment: int,
    values: dict[int, float] | None,
    cmap,
    norm,
    segment_colors: dict[int, str],
):
    """Resolve the fill color for one segment.

    If *values* contains an entry for *segment*, maps the scalar through *cmap*/*norm*.
    Otherwise falls back to the hex string in *segment_colors*, or grey if absent.
    """
    if values is not None and segment in values:
        return cmap(norm(values[segment]))
    return segment_colors.get(segment, "#d9d9d9")


def _text_color(facecolor) -> str:
    """Return 'black' or 'white' for legible text on top of *facecolor*."""
    rgba = mcolors.to_rgba(facecolor)
    luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
    return "black" if luminance > 0.58 else "white"


def _polar_to_xy(
    center: tuple[float, float],
    radius: float,
    theta_deg: float,
) -> tuple[float, float]:
    """Convert polar coordinates to Cartesian (x, y) offset from *center*.

    Parameters
    ----------
    center    : (x, y) origin in data coordinates.
    radius    : radial distance from center.
    theta_deg : angle in degrees measured counter-clockwise from the +x axis.
    """
    theta = np.deg2rad(theta_deg)
    return center[0] + radius * np.cos(theta), center[1] + radius * np.sin(theta)


def draw_lv(
    ax,
    values: dict[int, float] | None,
    cmap,
    norm,
    segment_colors: dict[int, str],
    center: tuple[float, float],
    show_labels: bool = True,
):
    """Draw all LV segments (1–17) onto *ax*.

    Segments 1–16 are Wedge patches arranged in three concentric rings defined
    by LV_RING_DEFS.  Segment 17 (apex) is a filled Circle at *center*.

    Parameters
    ----------
    ax            : Matplotlib Axes to draw on.
    values        : mapping {segment: scalar} used for colormap coloring, or
                    None to use *segment_colors*.
    cmap          : Matplotlib colormap (returned by plt.get_cmap).
    norm          : Matplotlib Normalize instance matching *cmap*.
    segment_colors: fallback {segment: hex_color} dict when *values* is None
                    or a segment has no entry in *values*.
    center        : (x, y) of the LV center in data coordinates (= LV_CENTER).
    show_labels   : if True, draw the segment number text inside each patch.
    """
    for segments, r_in, r_out, angles in LV_RING_DEFS:
        for i, segment in enumerate(segments):
            theta1 = angles[i]
            theta2 = angles[i + 1]
            facecolor = _value_color(segment, values, cmap, norm, segment_colors)
            ax.add_patch(
                Wedge(
                    center,
                    r_out,
                    theta1=theta1,
                    theta2=theta2,
                    width=r_out - r_in,
                    facecolor=facecolor,
                    edgecolor="#202020",
                    linewidth=1.3,
                )
            )
            if show_labels:
                x, y = _polar_to_xy(center, (r_in + r_out) / 2, (theta1 + theta2) / 2)
                ax.text(
                    x,
                    y,
                    str(segment),
                    ha="center",
                    va="center",
                    fontsize=9,
                    weight="bold",
                    color=_text_color(facecolor),
                )

    facecolor = _value_color(17, values, cmap, norm, segment_colors)
    ax.add_patch(
        Circle(center, 1.0, facecolor=facecolor, edgecolor="#202020", linewidth=1.3)
    )
    if show_labels:
        ax.text(
            center[0],
            center[1],
            "17",
            ha="center",
            va="center",
            fontsize=10,
            weight="bold",
            color=_text_color(facecolor),
        )


def draw_rv(
    ax,
    values: dict[int, float] | None,
    cmap,
    norm,
    segment_colors: dict[int, str],
    center: tuple[float, float],
    show_labels: bool = True,
):
    """Draw all RV segments (18–26) onto *ax*.

    The RV is represented as three depth bands (outer 18-20, middle 21-23,
    inner 24-26), each bounded by a pair of concentric half-ellipses defined
    in RV_HALF_ELLIPSES.  Each band is divided into top/middle/bottom rows
    by perpendicular-to-ellipse separators whose angles are specified per band
    in RV_BAND_DEFS.

    Parameters
    ----------
    ax            : Matplotlib Axes to draw on.
    values        : mapping {segment: scalar} used for colormap coloring, or
                    None to use *segment_colors*.
    cmap          : Matplotlib colormap (returned by plt.get_cmap).
    norm          : Matplotlib Normalize instance matching *cmap*.
    segment_colors: fallback {segment: hex_color} dict when *values* is None
                    or a segment has no entry in *values*.
    center        : (x, y) of the LV center in data coordinates (= LV_CENTER).
    show_labels   : if True, draw the segment number text inside each patch,
                    placed outside the LV circle using _rv_label_position.
    """
    for segments, outer_lvl, inner_lvl, sep_angles in RV_BAND_DEFS:
        for i, segment in enumerate(segments):
            facecolor = _value_color(segment, values, cmap, norm, segment_colors)
            path = _rv_segment_path(center, i, i + 1, outer_lvl, inner_lvl, sep_angles)
            ax.add_patch(
                PathPatch(
                    path,
                    facecolor=facecolor,
                    edgecolor="#202020",
                    linewidth=1.3,
                    joinstyle="round",
                )
            )
            if show_labels:
                x, y = _rv_label_position(path, center)
                ax.text(
                    x,
                    y,
                    str(segment),
                    ha="center",
                    va="center",
                    fontsize=8,
                    weight="bold",
                    color=_text_color(facecolor),
                    zorder=4,
                )


def _validate_model(model: int | str) -> int:
    model = int(model)
    if model not in (17, 26):
        raise ValueError(f"model must be 17 or 26, got {model}")
    return model


def draw_aha_model(
    values: dict[int, float] | None = None,
    csv: str = None,
    ax=None,
    *,
    model: int | str = DEFAULT_MODEL,
    cmap_name: str = "viridis",
    cmap=None,
    value_range: tuple[float, float] = (0.0, 1.0),
    segment_colors: dict[int, str] | None = None,
    show_labels: bool = True,
    show_lv_rv_labels: bool = True,
    show_colorbar: bool = True,
    colorbar_label: str | None = None,
    title: str | None = "AHA 26 Model",
):
    """Draw an AHA 17- or 26-segment bullseye model.

    This is the main entry point.  It creates (or reuses) a Matplotlib figure,
    draws the LV (and optionally the RV), and returns the figure and axes.

    Parameters
    ----------
    values : dict[int, float] or None
        Optional mapping of segment number → scalar value.  When provided the
        segments are colored by mapping the scalar through *cmap_name* and
        *value_range*.  Segments absent from *values* fall back to
        *segment_colors*.  Pass None (default) to use only *segment_colors*.
    csv: string or None
        Path to csv file containing histogram data. Assumes data encoded as 8 bit.
    ax : matplotlib.axes.Axes or None
        Axes to draw on.  If None a new figure of size FIGSIZE is created.
    model : 17 or 26
        ``17`` draws only the LV (segments 1–17).
        ``26`` draws the LV plus the RV (segments 18–26).
    cmap_name : str
        Name of any Matplotlib colormap, used when *values* is provided.
        Default ``"viridis"``.  For helix-angle data ``"twilight_shifted"``
        works well.
    cmap : colourmap
        Provide a colourmap directly instead of using one form matplotlib
    value_range : (float, float)
        (vmin, vmax) passed to ``plt.Normalize``.  Scalars outside this range
        are clipped by the colormap.  Default ``(0.0, 1.0)``.
    segment_colors : dict[int, str] or None
        Per-segment fallback colors as ``{segment_number: "#rrggbb"}``.
        Merged with DEFAULT_SEGMENT_COLORS (caller values take precedence).
        Pass None to use only DEFAULT_SEGMENT_COLORS.
    show_labels : bool
        Draw segment-number text inside each patch.  Default True.
    show_lv_rv_labels : bool
        Draw the "LV" and "RV" annotations (only meaningful for model 26).
    show_colorbar : bool
        Attach a colorbar to the figure when *values* is provided.
    colorbar_label : str or None
        Label shown next to the colorbar axis.
    title : str or None
        Figure title.  Pass None to suppress.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes

    Examples
    --------
    # Plain color map with default colors
    fig, ax = draw_aha_model(model=26)

    # Color each segment by a fractional anisotropy value
    fa_values = {1: 0.42, 2: 0.55, ..., 26: 0.38}
    fig, ax = draw_aha_model(fa_values, model=26,
                             cmap_name="viridis", value_range=(0, 1),
                             colorbar_label="FA")
    plt.show()
    """
    model = _validate_model(model)
    segment_colors = dict(DEFAULT_SEGMENT_COLORS if segment_colors is None else segment_colors)
    if not cmap:
        cmap = plt.get_cmap(cmap_name)
    norm = plt.Normalize(*value_range)

    if ax is None:
        fig, ax = plt.subplots(figsize=FIGSIZE)
    else:
        fig = ax.figure

    if not values:
        ## Load a histogram from a csv file 
        if csv: 
            hist= pd.read_csv(csv)
            hist.drop(columns=['VALUE'], inplace=True)

            hist.columns = [
                int(c.removeprefix('SEGMENT-')) if c.startswith('SEGMENT-') else c
                for c in hist.columns
            ]

            # Concatenate radial division
            if max(hist.columns) > model: 
                raw_values = {}

                for s in range(1, model + 1):
                    raw_values[s] = np.zeros(len(hist))

                    base = s * 10
                    for r in range(1, 5):
                        label = base + r
                        if label in hist.columns:
                            raw_values[s] += hist[label].to_numpy()
            else: 
                raw_values = {k: hist[k].to_numpy() for k in hist.columns}
            
            hist_array = np.column_stack([raw_values[k] for k in sorted(raw_values)])
            values = utils.segment_means_from_histogram(hist_array, circular=True)



    if model == 26:
        draw_rv(ax, values, cmap, norm, segment_colors, LV_CENTER, show_labels=show_labels)
    draw_lv(ax, values, cmap, norm, segment_colors, LV_CENTER, show_labels=show_labels)

    if show_lv_rv_labels and model == 26:
        x_ins, y_ins = _rv_insertion_xy()
        ax.text(
            LV_CENTER[0] + x_ins - RV_OUTER_RX * 0.45,
            LV_CENTER[1] + y_ins * 0.95,
            "RV", fontsize=12, weight="bold",
        )
        ax.text(LV_CENTER[0] + 3.35, LV_CENTER[1] + 3.35, "LV", fontsize=12, weight="bold")

    if title:
        ax.set_title(title, fontsize=15, weight="bold")
    ax.set_aspect("equal")
    if model == 26:
        x_ins, _ = _rv_insertion_xy()
        rv_left = x_ins - RV_OUTER_RX
        ax.set_xlim(LV_CENTER[0] + rv_left - 0.6, LV_CENTER[0] + 4.35)
        ax.set_ylim(LV_CENTER[1] - 4.25, LV_CENTER[1] + 4.25)
    else:
        ax.set_xlim(LV_CENTER[0] - 4.25, LV_CENTER[0] + 4.25)
        ax.set_ylim(LV_CENTER[1] - 4.25, LV_CENTER[1] + 4.25)
    ax.axis("off")

    if values is not None and show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.04)
        if colorbar_label:
            cbar.set_label(colorbar_label)

    fig.tight_layout()
    return fig, ax



def fake_values_for_model(model: int | str, kind: str = "fa") -> dict[int, float]:
    """Create deterministic fake scalar values for testing segment coloring.

    Parameters
    ----------
    model : 17 or 26
        Which model to generate values for (determines segment count).
    kind  : ``"fa"`` or ``"ha"``
        ``"fa"`` → linearly spaced values in [0, 1] (fractional anisotropy range).
        ``"ha"`` → linearly spaced values in [-90, 90] (helix-angle range).

    Returns
    -------
    dict[int, float] mapping each segment number to a scalar.
    """
    model = _validate_model(model)
    segments = list(range(1, model + 1))
    if kind == "ha":
        values = np.linspace(-90.0, 90.0, len(segments))
    elif kind == "fa":
        values = np.linspace(0.0, 1.0, len(segments))
    else:
        raise ValueError(f"kind must be 'fa' or 'ha', got {kind}")
    return dict(zip(segments, values))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draw/test the AHA 17/26 model geometry.")
    parser.add_argument("--model", choices=("17", "26"), default=str(DEFAULT_MODEL))
    parser.add_argument("--fake-values", choices=("none", "fa", "ha"), default="none")
    parser.add_argument("--output", default=OUTPUT_PATH)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--no-save", action="store_true")
    parser.add_argument("--no-title", action="store_true")
    parser.add_argument("--no-colorbar", action="store_true")
    parser.add_argument("--csv", default=False)
    parser.add_argument("--value-range", nargs = 2, type=float, default=[0.0, 1.0])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model = _validate_model(args.model)
    values = None
    cmap_name = "viridis"
    value_range = (0.0, 1.0)
    colorbar_label = None
    title = None if args.no_title else f"AHA {model} Model"

    if args.fake_values != "none":
        values = fake_values_for_model(model, args.fake_values)
        title = None if args.no_title else f"AHA {model} Model - fake {args.fake_values.upper()}"
        colorbar_label = f"Fake {args.fake_values.upper()} value"
        if args.fake_values == "ha":
            cmap_name = "twilight_shifted"
            value_range = (-90.0, 90.0)

    fig, _ = draw_aha_model(
        values=values,
        csv = args.csv,
        model=model,
        cmap_name=cmap_name,
        value_range=args.value_range,
        show_colorbar=not args.no_colorbar,
        colorbar_label=colorbar_label,
        title=title,
    )

    if not args.no_save:
        output = Path(args.output or f"aha{model}_model_test.png")
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved: {output.resolve()}")

    if SHOW_FIGURE and not args.no_show:
        print("Opening figure. Close the plot window to finish.")
        plt.show()

    plt.close(fig)


if __name__ == "__main__":
    main()
