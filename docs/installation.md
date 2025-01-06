# Installation

cardiotensor is a powerful and user-friendly toolkit for analyzing the orientation of cardiomyocites fibers in the heart

## Prerequisites

- Python 3.11 or newer
- Required libraries (see `pyproject.toml`)

## With pip <small>recommended</small>

cardiotensor is published as a [Python package] and can be installed with
`pip`, ideally by using a [virtual environment]. Open up a terminal and install
cardiotensor with:

``` sh
pip install cardiotensor 
```
        
  [Python package]: https://pypi.org/project/cardiotensor/
  [virtual environment]: https://realpython.com/what-is-pip/#using-pip-in-a-python-virtual-environment

## With git 

Clone the repository:

   ```bash
   git clone https://github.com/JosephBrunet/cardiotensor.git 
   cd cardiotensor
   ```

Install the package:

   ```bash
   pip install -e .  # (1)!
   ```

   1.  The `-e` flag in `pip install -e .` installs the package in editable mode, allowing changes to the source code to be immediately reflected without reinstallation.

Verify installation:

   ```bash
   cardio-tensor --help
   ```

---

  [GitHub]: https://github.com/squidfunk/mkdocs-material