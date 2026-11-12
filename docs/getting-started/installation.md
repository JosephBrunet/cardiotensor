# Installation

Cardiotensor is a powerful and user-friendly toolkit for analyzing the orientation of cardiomyocites fibers in the heart

## Prerequisites

- Python 3.11 or higher (tested on 3.11–3.13)
- A WebGPU-compatible graphics driver for FURY visualization

Cardiotensor uses FURY 2 with the PyGfx/WGPU renderer. VTK remains a
dependency for the PyVista backend and VTK export.

Upgrade the package and its dependencies with `python -m pip install --upgrade cardiotensor`.

## Installing with pip <small>recommended</small>

cardiotensor is published as a [Python package] and can be installed with
`pip`, ideally by using a [virtual environment]. Open up a terminal and install
cardiotensor with:

``` sh
pip install cardiotensor
```

  [Python package]: https://pypi.org/project/cardiotensor/
  [virtual environment]: https://realpython.com/what-is-pip/#using-pip-in-a-python-virtual-environment

## Installing from Source

To install cardiotensor from source, follow these steps:

1. Clone the repository from GitHub:

    ```console
    $ git clone https://github.com/JosephBrunet/cardiotensor.git
    ```

2. Navigate to the cloned repository directory:

    ```console
    $ cd cardiotensor
    ```

3. Install the package using pip:

    ```console
    $ pip install -e .  # (1)!
    ```

      1.  The `-e` flag in `pip install -e .` installs the package in editable mode, allowing changes to the source code to be immediately reflected without reinstallation.

## Uninstallation

If you need to remove cardiotensor from your system, follow these steps:

   ```console
   $ pip uninstall cardio-tensor
   ```

This should remove cardiotensor from your environment.
