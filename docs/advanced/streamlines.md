# Streamlines

This page covers advanced streamline generation and interpretation using `cardio-generate-streamlines` and `cardio-visualize-streamlines`.

## Goal

Reconstruct coherent 3D myocyte trajectories from the orientation vector field and analyze regional organization patterns.

## Prerequisites

Before generating streamlines:

- Run `cardio-tensor` with `WRITE_VECTORS = True`.
- Ensure orientation outputs are computed on the full volume.
- Verify FA and vector consistency on representative slices.

## Generate Streamlines

Basic example:

```console
$ cardio-generate-streamlines ./parameters_example.conf --seeds 10000 --start-z 150
```

Commonly tuned arguments:

- `--seeds`: number of initial seed points (higher values improve coverage but increase runtime).
- `--start-z`: start slice for generation when focusing on a sub-volume.
- `--fa-threshold`: minimum FA to continue tracking.
- `--angle`: maximum allowed turning angle between steps.
- `--min_len`: minimum streamline length to keep.

## Visualize Streamlines

```console
$ cardio-visualize-streamlines ./parameters_example.conf --line-width 1
```

Use thinner lines for dense datasets and thicker lines for sparse exploratory views.

## Robust Tuning Workflow

1. Start with moderate seeds and default thresholds.
2. Inspect streamline coverage and continuity.
3. Increase seeds to fill undersampled regions.
4. Tighten turning-angle and FA thresholds if many anatomically implausible curves appear.
5. Remove very short tracks using `--min_len`.

## Quality Control

Check the following before interpretation:

- Streamlines stay within myocardium boundaries.
- Trajectories remain smooth without abrupt direction flips.
- Regional transmural transitions are coherent with HA maps.
- Results are stable when repeating with nearby parameter values.

## Interpretation Notes

- High curvature clusters may indicate true architectural transitions or tracking instability.
- Large empty regions can come from low FA, low seeding density, or masking issues.
- Use transmural profiles (`cardio-analysis`) together with streamlines for stronger conclusions.

## Export and Reuse

Generated streamlines are saved in `.trk` format and can be:

- Reloaded for additional filtering/quantification.
- Visualized in 3D tools such as ParaView (after conversion if needed).
- Compared across samples using matched tracking settings.
