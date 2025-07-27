# Structure Tensor Computation

The structure tensor is a fundamental tool in 3D image analysis used to estimate local orientation and anisotropy. In `cardiotensor`, it is computed per voxel to quantify the principal direction of cardiomyocyte aggregates.

---

## What is a Structure Tensor?

The structure tensor $S$ is a 3-by-3 symmetric matrix that summarizes how intensity varies in a local neighborhood of a 3D image. It is defined as:

$$
S = K_\rho * (\nabla V_\sigma \cdot \nabla V_\sigma^T)
$$

Where:

* $V$: 3D image volume
* $\nabla V_\sigma$: smoothed gradient (via Gaussian derivative with noise scale $\sigma$)
* $K_\rho$: Gaussian kernel with standard deviation $\rho$ (integration scale)
* $*$: convolution over a local neighborhood

This matrix encodes orientation and intensity variation.

---

## Implementation

Cardiotensor computes the structure tensor using the [structure-tensor](https://github.com/Skielex/structure-tensor) Python library.

* It supports **parallel computation** and automatically uses the **GPU** (via CuPy) if available.
* CPU-based execution falls back to NumPy for full compatibility.
* This makes orientation estimation fast and scalable on large datasets.

---

## Step-by-Step Computation

1. **Noise Filtering**

   * Compute gradients $V_x, V_y, V_z$ using Gaussian derivative filters (standard deviation $\sigma$).

2. **Outer Product**

   * Form tensor components:

     * $V_x^2, V_y^2, V_z^2, V_xV_y, V_xV_z, V_yV_z$

3. **Smoothing**

   * Convolve each component with a Gaussian kernel (standard deviation $\rho$) to compute final tensor:

$$
S = \begin{bmatrix}
\langle V_x^2 \rangle & \langle V_x V_y \rangle & \langle V_x V_z \rangle \\
\langle V_x V_y \rangle & \langle V_y^2 \rangle & \langle V_y V_z \rangle \\
\langle V_x V_z \rangle & \langle V_y V_z \rangle & \langle V_z^2 \rangle
\end{bmatrix}
$$

---

## Eigen Decomposition

Each structure tensor $S$ is decomposed into eigenvalues $\lambda_1 \leq \lambda_2 \leq \lambda_3$ and corresponding eigenvectors $\vec{v}_1, \vec{v}_2, \vec{v}_3$.

* $\vec{v}_1$: direction with least intensity change (principal fiber direction)
* $\vec{v}_3$: direction of greatest intensity change

---

## Output and Interpretation

* The third eigenvector $\vec{v}_1$ is stored as the local myocyte orientation.
* Helix and intrusion angles are computed from $\vec{v}_1$ after transforming it into cylindrical coordinates.
* Fractional Anisotropy (FA) is computed using the three eigenvalues (see [FA section](./fractional-anisotropy.md)).

---

## Parameters

* **$\sigma$** *(noise scale)*: Controls how much local variation is smoothed in gradient calculation.
* **$\rho$** *(integration scale)*: Defines the neighborhood size for averaging the tensor field.

Higher $\sigma$ and $\rho$ values increase robustness to noise but reduce spatial precision.

---

