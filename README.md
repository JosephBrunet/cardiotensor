[![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://JosephBrunet.github.io/cardiotensor/)
[![License](https://img.shields.io/github/license/JosephBrunet/cardiotensor)](https://github.com/JosephBrunet/cardiotensor/blob/main/LICENSE)
![Python Version](https://img.shields.io/badge/python->3.11-blue.svg)
[![PyPI version](https://img.shields.io/pypi/v/caridotensor.svg)](https://pypi.org/project/cardiotensor/)

# CardioTensor

A Python package for calculating 3D cardiomyocyte orientations in heart images.

![Example Result](./examples/result_HA_slice.png)

## Overview

**cardiotensor** is a powerful and user-friendly toolkit for analyzing the orientation of cardiomyocites fibers in the heart. This package leverages advanced image processing techniques, enabling researchers to quantify 3D cardiomyocyte orientations efficiently. It is designed to support a wide range of datasets and provide accurate results for cardiac research.

## Features

- **Automated Orientation Analysis**: Extract 3D orientations of cardiac fibers with ease.
- **Configurable Workflow**: Customize the process using configuration files.
- **Integration Ready**: Designed for seamless integration with downstream analysis tools.

## Installation

Clone the repository and install the package using pip:

```bash
git clone https://github.com/JosephBrunet/cardiotensor.git
cd cardiotensor
pip install .
```

## Requirements

- Python 3.10 or newer

## Getting Started

1. **Prepare a Configuration File**:

   - Use the provided template (`/param_files/parameters_template.conf`) to create a configuration file tailored to your dataset.

2. **Run the orientation computation**:

   - Execute the following command in your terminal:
     ```bash
     cardio-tensor
     ```
   - When prompted, select your `.conf` file to start the analysis.

3. **View Results**:
   - The results will be saved in the specified output directory as defined in your configuration file.

## More Information

This package builds upon the [Structure Tensor package](https://github.com/Skielex/structure-tensor), extending its capabilities for cardiac imaging.

## License

This project is licensed under the MIT License. See the [LICENSE.md](./LICENSE.md) file for details.

## Contributing

Contributions are welcome! If you encounter a bug or have suggestions for new features:

- **Report an Issue**: Open an issue in the repository.
- **Submit a Pull Request**: Fork the repository, make changes, and submit a pull request.

For major changes, please discuss them in an issue first.

## Contact

For questions, feedback, or support, please contact the maintainers at [j.brunet@ucl.ac.uk].

---

Thank you for using **CardioTensor**! Your contributions and feedback are invaluable in improving the package.
