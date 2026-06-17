# Angle Definitions

This page explains how helical and intrusion angles are calculated from the 3D eigenvector field derived by `cardiotensor`.

By default, Cardiotensor reports **unprojected 3D angles**. The primary eigenvector is not first flattened onto a 2D plane before the angle is measured, because projection discards one component of the vector and can bias the resulting angle. Projection-based angles are available only for legacy comparison.

## Coordinate System

A transformation to a cylindrical coordinate system is defined for each voxel based on an approximation of the left ventricle (LV) centerline.

- **Radial (r)**: outward from the LV center
- **Circumferential (θ)**: tangential around the ventricle
- **Longitudinal (z)**: base to apex direction

To compute local fiber angles consistently, all eigenvectors are first rotated into this cylindrical coordinate frame. This alignment is performed using the Rodrigues rotation formula, which computes the minimal-angle rotation that maps the global reference axis (here the z-axis) onto the local longitudinal axis at each point. This allows a robust comparison of orientations across the myocardium.

## Helical Angle (HA)

The helical angle is defined as the angle between the primary myocyte-axis eigenvector \\( \vec{v}_1 \\) and the local circumferential plane.

In local cylindrical coordinates, with radial component \\(R\\), circumferential component \\(C\\), and longitudinal component \\(L\\), the unprojected helical angle is:

\\[
\mathrm{HA} = \arctan2\left(L, \sqrt{R^2 + C^2}\right)
\\]

It captures the transmural variation of fiber orientation from epicardium to endocardium.

Typical pattern:
- ~−60° at epicardium
- ~0° in mid-wall
- ~+60° at endocardium

## Intrusion Angle (IA)

The intrusion angle is the angle between the primary myocyte-axis eigenvector \\( \vec{v}_1 \\) and the **tangential plane** (longitudinal + circumferential).

Using the same local components, the unprojected intrusion angle is:

\\[
\mathrm{IA} = \arctan2\left(R, \sqrt{C^2 + L^2}\right)
\\]

It captures radial deviation of fiber aggregates and can help identify wall thickening or microstructural disruptions.

## Projection Bias

Conventional projected angles are computed after removing one vector component, for example \\(\arctan2(L, C)\\) for projected helical angle or \\(\arctan2(R, C)\\) for projected intrusion angle. These projected quantities can differ from the true 3D orientation when the discarded component is large.

Set `PROJECTED_ANGLES = True` only when you need legacy projected `HA_projected` and `IA_projected` maps for comparison with literature.

This convention follows the projection-error discussion in Agger et al., "Anatomically correct assessment of the orientation of the cardiomyocytes using diffusion tensor imaging", *NMR in Biomedicine* (2020), https://doi.org/10.1002/nbm.4205.

## Angle Ranges

Both angles are reported in degrees:
- **HA**: −90° to +90°
- **IA**: −90° to +90°

Angles are defined in a left-handed cylindrical coordinate system aligned to the LV.
