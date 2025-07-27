# Fractional Anisotropy (FA)

Fractional Anisotropy (FA) is a scalar value that quantifies the degree of anisotropy (directional dependence) in a tensor field. In `cardiotensor`, FA is derived from the structure tensor eigenvalues and helps identify regions with organized myocardial fiber orientation.

---

## Definition

Given the ordered eigenvalues of the structure tensor $\lambda_1 \leq \lambda_2 \leq \lambda_3$, the FA is computed as:

$$
FA = \sqrt{\frac{3}{2}} \cdot \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}
$$

Where:

* $\bar{\lambda} = (\lambda_1 + \lambda_2 + \lambda_3) / 3$: mean of the eigenvalues

---

## Interpretation

* **FA = 0**: Isotropic region, where diffusion or structural orientation is equal in all directions (e.g. cavity or noise).
* **FA = 1**: Perfectly anisotropic region, structure is highly aligned in one direction (e.g. dense aligned fibers).
* **0 < FA < 1**: Varying degrees of anisotropy.

This scalar map can be used to:

* Mask low-confidence areas for streamline seeding
* Identify anatomical regions of high or low alignment
* Visualize structural organization within the myocardium

---

## Output

The FA volume is saved as a single 3D scalar field with the same shape as the input image:

```python
fa.shape == (Z, Y, X)
```

It can be saved as TIFF or JP2 and visualized with standard volume rendering tools.

---

## Thresholding in Tractography

FA values can be used to restrict streamline propagation:

* **Minimum FA threshold**: ensures streamlines are seeded and propagated only in structured regions
* Default values range from 0.1 to 0.2 depending on noise level and resolution

---
