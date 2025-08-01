<h1 align="center">Cardiotensor</h1>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/JosephBrunet/cardiotensor/raw/main/assets/logos/heart_logo_dark.png">
    <source media="(prefers-color-scheme: light)" srcset="https://github.com/JosephBrunet/cardiotensor/raw/main/assets/logos/heart_logo_light.png">
    <img alt="Cardiotensor logo" src="https://github.com/JosephBrunet/cardiotensor/raw/main/assets/logos/heart_logo_light.png" width="200px">
  </picture>
</p>
<br />

<p align="center">A Python package to quantify and visualize 3D cardiomyocyte orientation in heart imaging datasets</p>

[![CI](https://github.com/JosephBrunet/cardiotensor/actions/workflows/ci.yml/badge.svg)](https://github.com/JosephBrunet/cardiotensor/actions/workflows/ci.yml)
[![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://JosephBrunet.github.io/cardiotensor/)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/b0e80972e3104ffa890532738882f42e)](https://app.codacy.com?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![License](https://img.shields.io/github/license/JosephBrunet/cardiotensor)](https://github.com/JosephBrunet/cardiotensor/blob/main/LICENSE)
[![Python Version](https://img.shields.io/pypi/pyversions/cardiotensor.svg)](https://pypi.python.org/pypi/cardiotensor)
[![PyPI version](https://img.shields.io/pypi/v/cardiotensor.svg)](https://pypi.org/project/cardiotensor/)


## Introduction

**Cardiotensor** is a user-friendly and memory-efficient toolkit designed for analyzing the orientation of cardiomyocyte fibers in large heart imaging datasets. It uses advanced image processing techniques to help researchers to accurately quantify 3D cardiomyocyte orientations with high efficiency.



## Installation

Clone the repository and install the package using pip:

```bash
git clone https://github.com/JosephBrunet/cardiotensor.git
cd cardiotensor
pip install .
```

⚠️ Require python 3.10 or newer


## Documentation

cardiotensor's documentation is available at [josephbrunet.fr/cardiotensor/](https://www.josephbrunet.fr/cardiotensor/)

## Getting Started

Have a look at our [simple example](https://www.josephbrunet.fr/cardiotensor/getting-started/examples/) that runs you through all the commands of the package

<p align="center">
    <img src="https://github.com/JosephBrunet/cardiotensor/raw/main/assets/images/pipeline.png"
         alt="Cardiotensor pipeline for 3D cardiac orientation analysis"
         style="max-width: 100%; border-radius: 6px;">
    <br>
    <em>
        <strong>Overview of the <code>cardiotensor</code> pipeline for 3D cardiac orientation analysis and tractography.</strong>
        <strong>(a)</strong> Input data consist of a whole‑ or partial‑heart volume and, optionally, a binary mask to restrict analysis to myocardial tissue.
        <strong>(b)</strong> Local cardiomyocyte orientation is derived by 3D structure tensor computation and eigenvector decomposition.
        The third eigenvector (smallest eigenvalue) is visualized as arrows, color‑coded by helix angle (HA); inset shows a zoom of the ventricular septum highlighting transmural fiber rotation.
        <strong>(c)</strong> After transforming to a cylindrical coordinate system aligned with the left ventricle, voxel‑wise HA, transverse angle (TA), and fractional anisotropy (FA) maps are computed for quantitative analysis.
        <strong>(d)</strong> Streamline tractography generated from the eigenvector field reveals continuous cardiomyocyte bundles throughout the heart, color‑coded by HA.
    </em>
</p>


## More Information

This package uses the [structure-tensor](https://github.com/Skielex/structure-tensor) package to calculate the structure tensor, extending its capabilities for cardiac imaging.

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/JosephBrunet/cardiotensor/blob/main/LICENSE) file for details.

## Contributing

Contributions are welcome! If you encounter a bug or have suggestions for new features:

- **Report an Issue**: Open an issue in the repository.
- **Submit a Pull Request**: Fork the repository, make changes, and submit a pull request.

For major changes, please discuss them in an issue first.

## Contact

For questions, feedback, or support, please contact the maintainers at [j.brunet@ucl.ac.uk].

## Reference

Brunet, J., Cook, A. C., Walsh, C. L., Cranley, J., Tafforeau, P., Engel, K., Arthurs, O., Berruyer, C., Burke O’Leary, E., Bellier, A., et al. (2024). Multidimensional analysis of the adult human heart in health and disease using hierarchical phase-contrast tomography. *Radiology, 312*(1), e232731. https://doi.org/10.1148/radiol.232731. [[PDF](https://pubs.rsna.org/doi/epdf/10.1148/radiol.232731)]

```bibtex
@article{brunet2024multidimensional,
  title={Multidimensional analysis of the adult human heart in health and disease using hierarchical phase-contrast tomography},
  author={Brunet, Joseph and Cook, Andrew C and Walsh, Claire L and Cranley, James and Tafforeau, Paul and Engel, Klaus and Arthurs, Owen and Berruyer, Camille and Burke O’Leary, Emer and Bellier, Alexandre and others},
  journal={Radiology},
  volume={312},
  number={1},
  pages={e232731},
  year={2024},
  publisher={Radiological Society of North America}
}
```
