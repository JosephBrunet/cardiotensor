import numpy as np
import pytest

from cardiotensor.visualization import fury_plotting_vectors as plotting


@pytest.mark.parametrize("value", [0.0, np.nan, np.inf])
@pytest.mark.parametrize("with_colors", [False, True])
def test_empty_vectors_skip_rendering(value, with_colors, monkeypatch, capsys):
    def unexpected_scene():
        pytest.fail("An empty selection must not open a rendering window")

    monkeypatch.setattr(plotting.window, "Scene", unexpected_scene)
    vectors = np.full((2, 2, 2, 3), value)
    colors = np.ones((2, 2, 2)) if with_colors else None
    plotting.plot_vector_field_fury(vectors, color_volume=colors, downsample=1)
    assert "No finite, nonzero vectors" in capsys.readouterr().out


@pytest.mark.parametrize("mode", ["arrow", "cylinder"])
def test_fury2_glyphs_preserve_positions_and_mask_opacity(mode, monkeypatch, tmp_path):
    from types import SimpleNamespace

    from cardiotensor.visualization.fury_support import materials

    captured = {}

    def manager(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(close=lambda: None)

    monkeypatch.setattr(plotting.window, "ShowManager", manager)
    monkeypatch.setattr(plotting, "capture", lambda *args: None)
    vectors = np.zeros((1, 1, 2, 3), dtype=np.float32)
    vectors[..., 0] = 1
    plotting.plot_vector_field_fury(
        vectors,
        mode=mode,
        downsample=1,
        size=0.01,
        color_volume=np.array([[[0.0, 90.0]]]),
        opacity_mask=np.array([[[False, True]]]),
        outside_opacity=0.1,
        save_path=tmp_path / "vectors.png",
    )
    glyphs = [
        obj
        for obj in captured["scene"].main_scene.children
        if getattr(obj, "geometry", None) is not None
    ]
    assert len(glyphs) == 2
    assert next(materials(glyphs[0])).opacity == pytest.approx(0.1)
    assert next(materials(glyphs[1])).opacity == 1.0
    centers = [obj.get_world_bounding_box().mean(axis=0) for obj in glyphs]
    np.testing.assert_allclose(centers[1] - centers[0], [1, 0, 0], atol=1e-5)
    assert not np.allclose(
        glyphs[0].geometry.colors.data[0], glyphs[1].geometry.colors.data[0]
    )
