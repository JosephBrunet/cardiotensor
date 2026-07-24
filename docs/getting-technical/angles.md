# Angle Definitions

Cardiotensor computes helical angle (HA) and intrusion angle (IA) from the local myocyte-axis eigenvector.

By default, these are **unprojected 3D angles**. The vector is not flattened onto a 2D plane before measuring HA or IA, because projection removes one component of the vector and can bias the result.

## Local Components

Cardiotensor expresses the vector in a centerline-based cylindrical frame for the left ventricle:

- **R**: radial component, outward from the LV centerline
- **C**: circumferential component, tangent around the ventricle
- **L**: longitudinal component, along the LV axis

This is a practical LV cylindrical approximation. It should not be read as a full reproduction of anatomically corrected methods that explicitly model epicardial curvature.

## Helical Angle

The helical angle measures how far the vector points out of the local horizontal, or short-axis, plane. That plane is spanned by the radial and circumferential directions.

```text
HA = atan2(L, sqrt(R^2 + C^2))
```

Typical values are around -60 degrees at the epicardium, 0 degrees in the mid-wall, and +60 degrees at the endocardium.

## Intrusion Angle

The intrusion angle measures radial deviation from the local tangential plane. That plane is spanned by the circumferential and longitudinal directions.

```text
IA = atan2(R, sqrt(C^2 + L^2))
```

## Projected Legacy Angles

Set `PROJECTED_ANGLES = True` only when comparing with legacy projection-based maps:

```text
HA_projected = atan2(L, C)
IA_projected = atan2(R, C)
```

`IA_projected` corresponds to the projected transverse angle terminology often used in the literature.

Projection-based angles can differ from the true 3D orientation when the discarded component is large. This projection bias is discussed in Agger et al., "Anatomically correct assessment of the orientation of the cardiomyocytes using diffusion tensor imaging", *NMR in Biomedicine* (2020), https://doi.org/10.1002/nbm.4205.

## Angle Ranges

Both angles are reported in degrees:

- **HA**: -90 to +90 degrees
- **IA**: -90 to +90 degrees
