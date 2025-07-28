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
    orcid: 
    affiliation: "1, 2"
  - name: Lisa Chestnutt
    orcid: 
    affiliation: 3
  - name: Andrew Cook
    orcid: 
    affiliation: 3
  - name: Hector Dejea
    orcid: 
    affiliation: 2
  - name: Matthieu Chourrout
    orcid: 
    affiliation: 1
  - name: David Stansby
    orcid: 
    affiliation: 1
  - name: Peter D. Lee
    orcid: 
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
date: 26 June 2025
bibliography: paper.bib
---

# Summary

Understanding the architecture of the human heart requires analyzing its microstructural organization across scales. With the advent of high-resolution imaging techniques such as synchrotron-based tomography, it has become possible to visualize entire hearts at micron-scale resolution. However, translating these large, complex volumetric datasets into interpretable, quantitative descriptors of cardiac organization remains a major bottleneck.

Cardiotensor is an open-source Python package designed to quantify 3D cardiomyocyte orientation in whole-heart imaging datasets. It provides efficient, scalable implementations of structure tensor analysis, enabling extraction of directional metrics such as helix angle (HA), transverse angle (TA), and fractional anisotropy (FA). The package supports datasets reaching teravoxel scale and is optimized for high-performance computing environments, including parallel and chunk-based processing pipelines. In addition, cardiotensor includes tractography functionality to reconstruct continuous cardiomyocyte trajectories. This enables fiber-level visualization and structural mapping of cardiac tissue, supporting detailed assessments of anatomical continuity and regional organization.

By enabling scalable and reproducible analysis of cardiac microstructure, cardiotensor helps researchers study heart development, disease, and anatomy in 3D.



# Statement of Need

Despite major advances in high-resolution 3D imaging, there is a lack of scalable, open-source tools to analyze cardiomyocyte orientation in large volumetric datasets. Most established frameworks were developed for diffusion tensor MRI (DT-MRI), where orientation is inferred from water diffusion. Examples include MRtrix3 [@tournier_mrtrix3_2019], DIPY [@garyfallidis_dipy_2014], and DSI Studio [@yeh_dsi_2010]. While powerful for diffusion-based neuroimaging and cardiac applications [@mekkaoui_diffusion_2017], these packages are not designed to handle direct image-gradient–based orientation estimation or the teravoxel-scale datasets produced by synchrotron tomography, micro-CT, or optical imaging.

For non-diffusion imaging modalities, researchers have historically relied on custom structure tensor implementations to estimate fiber orientation directly from image intensity gradients. However, most of these are in-house codes, often unpublished or not generalizable. For example, structure tensor analysis has been applied in the heart using micro-CT [@], optical projection tomography [@], confocal microscopy [@dileep2023], episcopic microscopy [@], and synchrotron tomography [@], but these methods were tailored to specific datasets and lacked scalability or public availability.

These datasets present unique computational and analytical challenges, including memory constraints, limited processing throughput, and the need for spatially coherent orientation quantification across large fields of view. Moreover, the diversity of contrast mechanisms in non-diffusion imaging modalities requires algorithms that do not rely on water diffusion but instead exploit local image gradient patterns.

`cardiotensor` addresses this methodological gap by offering an open-source Python package specifically tailored to structure tensor analysis of high-resolution cardiac volumes. Rather than relying on diffusion modeling, `cardiotensor` infers tissue orientation directly from image intensity gradients, making it applicable across a wide range of modalities. Previous studies have demonstrated strong agreement between structure tensor–based orientation and DT-MRI–derived metrics when applied to the same human hearts [@teh_validation_2016].

The package supports full pipelines from raw image stacks to fiber orientation maps, HA and TA computation, FA, and tractography. Its architecture is optimized for large datasets, using chunked and parallel processing suitable for high-performance computing environments.

`cardiotensor` has already been successfully applied in published work to characterize 3D cardiomyocyte architecture in healthy and diseased human hearts using synchrotron tomography [@brunet_multidimensional_2024] to datasets over a terabyte in size, demonstrating its robustness and scalability.

 <img src="https://github.com/JosephBrunet/cardiotensor/raw/main/paper/figs/pipeline.jpg" alt="Helix angle map computed from a human heart dataset using `cardiotensor`." style="max-width: 80%">

**Figure 1**: Helix angle map computed from a human heart dataset using `cardiotensor`.

Image Stack + rendeirng -> vector field -> HA/IA/FA + Fiber tracing

The package also supports centerline interpolation and alignment to anatomical axes, which is useful for regional analysis of the heart.

# Architecture

The core functionality of `cardiotensor` is organized into four main modules:

- **`orientation/`**: Implements the structure tensor framework for estimating local cardiomyocyte orientation. This includes structure-tensor computation, eigenvalue decomposition, rotation to cylindrical coordinate system, and quantification of helix and transverse angles, as well as fractional anisotropy (FA).

- **`tractography/`**: Provides tools for generating and filtering streamlines that trace cardiomyocyte trajectories based on the orientation field. This module enables fiber-level reconstruction and visualization.

- **`analysis/`**: Contains higher-level methods that integrate orientation and tractography data for regional or statistical analysis. It also supports interpolation, centerline alignment, and cardiac anatomical mapping.

- **`utils/`**: Includes general-purpose utilities such as I/O functions, image preprocessing, vector math, and configuration parsing. These functions support the broader package infrastructure.

This modular structure promotes clarity, reproducibility, and ease of extension for cardiac imaging researchers working with large 3D datasets.


## Online Documentation

The documentation for cardiotensor is available online at:

**[https://josephbrunet.github.io/cardiotensor](https://josephbrunet.github.io/cardiotensor)**

The main components of the documentation are:

* Step-by-step walkthroughs for installation, first steps, and a guided example covering all available commands. A small example dataset and its corresponding mask are provided with the package.
* In-depth explanations of the core algorithms used in cardiotensor, including structure tensor theory, helix angle calculation, fractional anisotropy (FA), and tractography integration.
* Reference guides for the command-line interface, configuration file format, and public API.

# Acknowledgements

This project has been made possible in part by grant number 2022-316777 from the Chan Zuckerberg Initiative DAF, an advised fund of Silicon Valley Community Foundation.

The authors also acknowledge ESRF beamtimes md1252, md1290, and md1389 as sources of the data.

AC’s research is enabled through the Noé Heart Centre Laboratories, which are gratefully supported by the Rachel Charitable Trust via Great Ormond Street Hospital Children’s Charity (GOSH Charity). The Noé Heart Centre Laboratories are based in The Zayed Centre for Research into Rare Disease in Children, which was made possible thanks to Her Highness Sheikha Fatima bint Mubarak, wife of the late Sheikh Zayed bin Sultan Al Nahyan, founding father of the United Arab Emirates, as well as other generous funders.

# References

