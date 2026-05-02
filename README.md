<h1 align="center">CardioTensor</h1>

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
[![Doc](https://img.shields.io/badge/docs-online-blue.svg)](https://JosephBrunet.github.io/cardiotensor/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/JosephBrunet/cardiotensor/blob/main/LICENSE)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.09720/status.svg)](https://doi.org/10.21105/joss.09720)
[![PyPI version](https://img.shields.io/pypi/v/cardiotensor.svg)](https://pypi.org/project/cardiotensor/)


## Introduction

**Cardiotensor** is a user-friendly and memory-efficient toolkit designed for analyzing the orientation of cardiomyocyte fibers in large heart imaging datasets. It uses advanced image processing techniques to help researchers to accurately quantify 3D cardiomyocyte orientations with high efficiency.



## Installation

**Cardiotensor** is published as a [Python package](https://pypi.org/project/cardiotensor/) and can be installed with
`pip`, ideally by using a [virtual environment](https://realpython.com/what-is-pip/). Open up a terminal and install
cardiotensor with:

```bash
pip install cardiotensor
```

## Documentation

**Cardiotensor**'s documentation is available at [josephbrunet.fr/cardiotensor/](https://www.josephbrunet.fr/cardiotensor/)

## Getting Started

Have a look at our [simple example](https://www.josephbrunet.fr/cardiotensor/getting-started/examples/) that runs you through all the commands of the package

<p align="center">
    <img src="https://github.com/JosephBrunet/cardiotensor/raw/main/assets/images/pipeline.png"
         alt="Cardiotensor pipeline for 3D cardiac orientation analysis"
         style="max-width: 100%; border-radius: 6px;">
    <br>
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

If you use **Cardiotensor** in your research, please cite:

### Primary Reference

> Brunet, J., Chestnutt, L., Chourrout, M., Dejea, H., Sabarigirivasan, V., Lee, P. D., & Cook, A. C. "Cardiotensor: A Python Library for Orientation Analysis and Tractography in 3D Cardiac Imaging." *Journal of Open Source Software* 11(121) (2026): 9720. [🗎](https://doi.org/10.21105/joss.09720)

``` bibtex
@article{Brunet2026,
  doi={10.21105/joss.09720},
  url={https://doi.org/10.21105/joss.09720},
  year={2026},
  publisher={The Open Journal},
  volume={11},
  number={121},
  pages={9720},
  author={Brunet, Joseph and Chestnutt, Lisa and Chourrout, Matthieu and Dejea, Hector and Sabarigirivasan, Vaishnavi and Lee, Peter D. and Cook, Andrew C.},
  title={Cardiotensor: A Python Library for Orientation Analysis and Tractography in 3D Cardiac Imaging},
  journal={Journal of Open Source Software}
}
```

### Other Papers

> Brunet, J., Cook, A. C., Walsh, C. L., Cranley, J., Tafforeau, P., Engel, K., Arthurs, O., Berruyer, C., Burke O'Leary, E., Bellier, A., et al. "Multidimensional analysis of the adult human heart in health and disease using hierarchical phase-contrast tomography." *Radiology* 312(1) (2024): e232731. [🗎](https://doi.org/10.1148/radiol.232731)
