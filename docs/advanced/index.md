# Advanced

This section focuses on advanced workflows in `cardiotensor` for users who already run the standard orientation pipeline and want deeper analysis or larger 3D exploration.

## What You Can Do Here

- Perform detailed transmural quantification with the `cardio-analysis` GUI.
- Generate and inspect 3D streamlines from the computed orientation field.
- Apply practical quality-control checks before interpreting biological patterns.

## Recommended Prerequisites

Before starting these workflows, make sure you have:

- Run `cardio-tensor` on your full volume (not only test mode).
- Saved orientation vectors with `WRITE_VECTORS = True`.
- Generated angle maps (`HA`, `IA`) and `FA` outputs.

See:

- [Configuration file](../reference/configuration.md)
- [CLI commands](../reference/cli.md)
- [Example workflow](../getting-started/examples.md)

## Advanced Topics

- [Cardio Analysis](./cardio-analysis.md): in-depth usage of the transmural profile GUI and reproducible measurement strategy.
- [Streamlines](./streamlines.md): robust streamline generation, filtering, and interpretation in 3D.

## Suggested Order

1. Validate orientation and angle outputs on representative slices.
2. Extract transmural profiles in key anatomical regions.
3. Build 3D streamline visualizations and compare with 2D transmural trends.
