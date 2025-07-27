# Angle Definitions

This page explains how helix and intrusion angles are calculated from the 3D eigenvector field derived by `cardiotensor`.

## Coordinate System

A transformation to cylindrical coordinate system is defined for each voxel based on an approximation of the left ventricle (LV) centerline.

- **Radial (r)**: outward from the LV center
- **Circumferential (θ)**: tangential around the ventricle
- **Longitudinal (z)**: base to apex direction

## Helix Angle (HA)

The helix angle is defined as the angle between the third eigenvector \\( \vec{v}_3 \\) (smallest eigenvalue direction) and the **local circumferential plane**.

It captures the transmural variation of fiber orientation from epicardium to endocardium.

Typical pattern:
- ~−60° at epicardium
- ~0° in mid-wall
- ~+60° at endocardium

## Intrusion Angle (IA)

The intrusion angle is the angle between \\( \vec{v}_3 \\) and the **tangential plane** (longitudinal + circumferential).

It captures radial deviation of fiber aggregates and can help identify wall thickening or microstructural disruptions.

## Angle Ranges

Both angles are reported in degrees:
- **HA**: −90° to +90°
- **IA**: −90° to +90°

Angles are defined in a left-handed cylindrical coordinate system aligned to the LV.
