# Installation

cardiotensor is a powerful and user-friendly toolkit for analyzing the orientation of cardiomyocites fibers in the heart

## Prerequisites

- Python 3.10 or higher

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

    ```bash
    git clone https://github.com/JosephBrunet/cardiotensor.git
    ```

2. Navigate to the cloned repository directory:

    ```bash
    cd cardiotensor
    ```

3. Install the package using pip:

    ```bash
    pip install -e .  # (1)!
    ```

      1.  The `-e` flag in `pip install -e .` installs the package in editable mode, allowing changes to the source code to be immediately reflected without reinstallation.

## Verify the installation

To verify that cardiotensor is installed correctly, you can import it in Python:

   ```bash
   cardio-tensor --help
   ```

If no errors occur, the installation is successful.


---