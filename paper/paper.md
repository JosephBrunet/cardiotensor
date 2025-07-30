---
title: 'Cardiotensor: A Python Library for Orientation Analysis and Tractography in 3D Cardiac Imaging'
tags:
  - Python
  - cardiac imaging
  - structure tensor
  - orientation analysis
  - histoanatomy
  - fiber architecture
  - heart
authors:
  - name: Joseph Brunet
    orcid: 0000-0002-8424-9510
    affiliation: "1, 2"
  - name: Andrew Cook
    orcid: 0000-0001-5079-7546
    affiliation: 3
  - name: Lisa Chestnutt
    orcid: 0009-0008-0783-7299
    affiliation: 3
  - name: Hector Dejea
    orcid: 0000-0003-2584-9812
    affiliation: 2
  - name: Vaishnavi Sabarigirivasan
    orcid: 0000-0003-2550-6262
    affiliation: 2
  - name: Matthieu Chourrout
    orcid: 0000-0002-2282-6976
    affiliation: 1
  - name: David Stansby
    orcid: 0000-0002-1365-1908
    affiliation: 1
  - name: Peter D. Lee
    orcid: 0000-0002-3898-8881
    affiliation: "1, 4"

affiliations:
  - name: Department of Mechanical Engineering, University College London, London, UK
    index: 1
  - name: European Synchrotron Radiation Facility, Grenoble, France
    index: 2
  - name: UCL Institute of Cardiovascular Science, London, UK
    index: 3
  - name: Research Complex at Harwell, Didcot, UK
    index: 4
date: 29 July 2025
bibliography: paper.bib
---

# Summary

Understanding the architecture of the human heart requires analyzing its microstructural organization across scales. With the advent of high-resolution imaging techniques such as synchrotron-based tomography, it has become possible to visualize entire hearts at micron-scale resolution. However, translating these large, complex volumetric datasets into interpretable, quantitative descriptors of cardiac organization remains a major challenge. Cardiotensor is an open-source Python package designed to quantify 3D cardiomyocyte orientation in whole-heart imaging datasets. It provides efficient, scalable implementations of structure tensor analysis, enabling extraction of directional metrics such as helix angle (HA), transverse angle (TA), and fractional anisotropy (FA). The package supports datasets reaching teravoxel scale and is optimized for high-performance computing environments, including parallel and chunk-based processing pipelines. In addition, cardiotensor includes tractography functionality to reconstruct continuous cardiomyocyte trajectories. This enables fiber-level visualization and structural mapping of cardiac tissue, allowing detailed assessments of anatomical continuity and regional organization.

# Statement of Need

Despite major advances in high-resolution 3D imaging, there is a lack of open-source tools to analyze cardiomyocyte orientation in large volumetric datasets. Most established frameworks were developed for diffusion tensor MRI (DT-MRI), where orientation is inferred from water diffusion. Examples include MRtrix3 [@tournier_mrtrix3_2019], DIPY [@garyfallidis_dipy_2014], and DSI Studio [@yeh_dsi_2010]. While powerful for diffusion-based neuroimaging and cardiac applications [@mekkaoui_diffusion_2017], these packages are not designed to handle direct image-gradient–based orientation estimation or the teravoxel-scale datasets produced by synchrotron tomography, micro-CT, or optical imaging.

For non-diffusion imaging modalities, researchers have historically relied on custom structure tensor implementations to estimate fiber orientation directly from image intensity gradients. However, most of these are in-house codes, often unpublished or not generalizable. For example, structure tensor analysis has been applied in the heart using micro-CT [@reichardt_fiber_2020], microscopy [@dileep_cardiomyocyte_2023; @garcia-canadilla_detailed_2022], and synchrotron tomography [@dejea_comprehensive_2019], but these methods were tailored to specific datasets and lacked scalability or public availability.

Cardiotensor addresses this gap by providing an open-source Python package specifically tailored to structure tensor analysis of large cardiac volumes. Rather than relying on diffusion modeling, cardiotensor infers tissue orientation directly from image intensity gradients, making it applicable across a wide range of modalities. Previous studies have demonstrated strong agreement between structure tensor–based orientation and DT-MRI–derived metrics when applied to the same human hearts [@teh_validation_2016]. The package supports full pipelines from raw image stacks to fiber orientation maps and tractography. Its architecture is optimized for large datasets, using chunked and parallel processing suitable for high-performance computing environments.

Cardiotensor has already been successfully applied in published work to characterize 3D cardiomyocyte architecture in healthy and diseased human hearts using synchrotron tomography [@brunet_multidimensional_2024] to datasets over a terabyte in size. While cardiotensor was conceived for cardiac imaging, the package is modality‑ and tissue‑agnostic. Any volumetric dataset exhibiting coherent fibrous or laminar microstructure can be analyzed, including brain white matter, skeletal muscle, and tendon. This generality makes the library useful for both cardiovascular and broader anatomical or histological studies.

![Helix angle map computed from a human heart dataset using `cardiotensor`.](figs/pipeline.png)

**Figure 1**: **Overview of the `cardiotensor` pipeline for 3D cardiac orientation analysis and tractography.**  
(a) Input data consist of a whole‑ or partial-heart image volume and, optionally, a binary mask to restrict analysis to myocardial tissue.
(b) Local cardiomyocyte orientation is estimated by computing the 3D structure tensor and performing eigenvector decomposition. The third eigenvector field (smallest eigenvalue) is visualized as arrows color‑coded by helix angle (HA); the inset shows a zoomed view of the ventricular septum highlighting the transmural fiber rotation.
(c) After transforming to a cylindrical coordinate system aligned with the left ventricle, voxel‑wise helical angle (HA), transverse angle (TA), and fractional anisotropy (FA) maps are computed for quantitative tissue analysis.
(d) Streamline tractography is generated from the vector field, revealing continuous cardiomyocyte bundles trajectories throughout the heart, color‑coded by HA.

## Implementation

Cardiotensor is implemented in pure Python and designed to efficiently process very large 3D cardiac imaging datasets. It relies primarily on NumPy [@van_der_walt_numpy_2011] for numerical computation, with I/O accelerated by tifffile [@gohlke_cgohlketifffile_2025], Glymur for JPEG 2000 volumes [@evans_quintusdiasglymur_2025], and OpenCV [@bradski_opencv_2000]. Dask [@rocklin_dask_2015] is used exclusively to parallelize file reading, while the core computations rely on Python’s multiprocessing module for local parallelism. The package builds on the structure-tensor library [@jeppesen_quantifying_2021] to calculate the 3D structure tensor and eigenvector decomposition.

The package supports multiple use cases:

- Command‑line workflows, which automate batch processing of terabyte‑scale heart volumes and produce results as live plots or files saved to disk.

- Embedded use in larger Python pipelines for specific cardiac imaging analysis.

Efficient computation is achieved through a chunk‑based processing strategy with padding, which avoids edge artifacts. This architecture allows cardiotensor to process whole‑heart volumes in hours rather than days while maintaining practical memory requirements, and can be parallelized across a computing cluster by splitting volumes into independent jobs.


# Architecture

The core functionality of cardiotensor is organized into four main modules:

- **`orientation`**: Implements the chunked-based structure tensor framework for estimating local cardiomyocyte orientation. This includes structure-tensor computation, eigenvalue decomposition, rotation to cylindrical coordinate system, and quantification of helix and transverse angles, as well as fractional anisotropy (FA).

- **`tractography`**: Provides tools for generating and filtering streamlines that trace cardiomyocyte trajectories based on the orientation field. This module enables fiber-level reconstruction and visualization.

- **`analysis`**: Contains higher-level methods that integrate orientation and tractography data for regional or statistical analysis. It also supports interpolation, centerline alignment, and cardiac anatomical mapping.

- **`utils`**: Includes general-purpose utilities such as I/O functions, image preprocessing, vector math, and configuration parsing. These functions support the broader package infrastructure. A general high-level volume reader is implemented using Dask for parallelisation.

This modular structure promotes clarity, reproducibility, and ease of extension for cardiac imaging researchers working with large 3D datasets.


## Documentation and Usage

The documentation for cardiotensor is available online at:

**[https://josephbrunet.github.io/cardiotensor](https://josephbrunet.github.io/cardiotensor)**

The main components of the documentation are:

* Step-by-step walkthroughs for installation, first steps, and a guided example covering all available commands. A small example dataset and its corresponding mask are provided with the package.
* In-depth explanations of the core algorithms used in cardiotensor, including structure tensor theory, helix angle calculation, fractional anisotropy (FA), and tractography integration.
* Reference guides for the command-line interface, configuration file format, and public API.

# Acknowledgements

This project has been made possible in part by grant number 2022-316777 from the Chan Zuckerberg Initiative DAF, an advised fund of Silicon Valley Community Foundation.

The authors also acknowledge ESRF beamtimes md1290 and md1389 as sources of the data.

AC’s research is enabled through the Noé Heart Centre Laboratories, which are gratefully supported by the Rachel Charitable Trust via Great Ormond Street Hospital Children’s Charity (GOSH Charity). The Noé Heart Centre Laboratories are based in The Zayed Centre for Research into Rare Disease in Children, which was made possible thanks to Her Highness Sheikha Fatima bint Mubarak, wife of the late Sheikh Zayed bin Sultan Al Nahyan, founding father of the United Arab Emirates, as well as other generous funders.

# References

