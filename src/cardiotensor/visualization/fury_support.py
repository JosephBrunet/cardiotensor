"""Shared FURY 2 / PyGfx rendering helpers."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pygfx as gfx
from fury import ui


class ScalarColors:
    """Keep numeric limits together with their display colormap."""

    def __init__(self, cmap, value_range):
        self.cmap = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
        self.value_range = tuple(map(float, value_range))

    def map(self, values):
        low, high = self.value_range
        normalized = (np.asarray(values) - low) / (high - low or 1.0)
        return np.asarray(self.cmap(normalized), dtype=np.float32)


def materials(obj):
    """Yield materials even when FURY batches a large actor into a group."""
    if getattr(obj, "material", None) is not None:
        yield obj.material
    for child in obj.children:
        yield from materials(child)


def set_opacity(obj, opacity):
    for material in materials(obj):
        material.opacity = float(opacity)


def scalar_bar(colors, title):
    panel = ui.Panel2D(size=(265, 70), position=(12, 12), opacity=0.7)
    label = ui.TextBlock2D(text=title, font_size=14, size=(250, 22))
    panel.add_element(label, (5, 45))
    for i, color in enumerate(colors.cmap(np.linspace(0, 1, 64))):
        strip = ui.Rectangle2D(size=(4, 12), color=tuple(color[:3]))
        panel.add_element(strip, (5 + i * 4, 28))
    low, high = colors.value_range
    limits = ui.TextBlock2D(
        text=f"{low:.3g}                         {high:.3g}",
        font_size=13,
        size=(250, 22),
    )
    panel.add_element(limits, (5, 3))
    return panel


def capture(showm, path, scale=1):
    """Capture the current camera and overlays at the requested resolution."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    previous = showm.renderer.pixel_ratio
    try:
        showm.renderer.pixel_ratio = previous * scale
        showm.render()
        # request_draw schedules a frame; force it before reading the GPU texture.
        if hasattr(showm.window, "force_draw"):
            showm.window.force_draw()
        else:
            showm.window.draw()
        pixels = showm.snapshot(fname=str(path))
        return pixels.shape[1], pixels.shape[0]
    finally:
        showm.renderer.pixel_ratio = previous


def camera_for_scene(scene, size):
    camera = gfx.PerspectiveCamera(45, aspect=size[0] / size[1])
    light = gfx.DirectionalLight(intensity=1.0)
    camera.add(light)
    scene.add(camera)
    camera.show_object(scene.main_scene)
    return camera
