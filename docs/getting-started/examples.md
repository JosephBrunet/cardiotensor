# Example

This example demonstrates how to use the **cardiotensor** package with the example provided in the repository. This example will help you understand how to process input data, compute fiber orientations, and visualize results.

## Example Directory Overview

The `./examples/` directory contains:

- **Datasets**: Cropped 3D volumes and binary masks for testing and experimentation.
- **Configuration File**: A pre-filled `parameters_example.conf` file for running the workflows.

### Data

1. **Heart Volume**:
    - Path: `./data/635.2um_LADAF-2021-17_heart_overview_/`
    - Description: A downsampled 3D heart image volume designed for testing. The full-resolution dataset is available at the [Human Organ Atlas](https://human-organ-atlas.esrf.fr/datasets/1659197537).

2. **Binary Mask**:
    - Path: `./data/mask/`
    - Description: A binary mask used for segmenting the heart from the background.

## Running the Examples

### Installation

!!! note

    See [Installation](installation.md)

1. Clone the repository and install the package:
    ```console
    $ pip install cardiotensor
    ```

2. Navigate to the `examples` directory:
    ```console
    $ cd examples
    ```

### Processing a Test Slice

1. Open `parameters_example.conf` and set `TEST = True` in the `[TEST]` section.

!!! note

    For information about conf file see the section [Configuration file](../reference/configuration.md)

2. Run the following command:
    ```console
    $ cardio-tensor ./parameters_example.conf
    ```
3. The output will be displayed as a plot for a single slice:

<figure markdown="span">
![Test Slice Result](../assets/images/result_test_slice.png){ width="70%" }
<figcaption>Result from processing a single test slice.</figcaption>
</figure>

### Processing the Entire Volume

1. Set `TEST = False` in `parameters_example.conf`.

!!! note

    For information about conf file see the section [Configuration file](../reference/configuration.md)

2. Run the command:
    ```console
    $ cardio-tensor ./parameters_example.conf
    ```
3. Outputs will be saved in the `./output` directory with the following structure:
    ```
    ./output
    ├── HA          # Helix angle results
    ├── IA          # Intrusion angle results
    ├── FA          # Fractional anisotropy results
    └── eigen_vec   # 3rd Eigenvectors
    ```

### Visualizing Transmural Profiles

1. Use the `cardio-analysis` command:
    ```console
    $ cardio-analysis ./parameters_example.conf 150
    ```
    Replace `150` with the slice number you wish to analyze.

2. The GUI will appear, allowing you to:
    - Define a transmural profile line.
    - Adjust parameters like `Angle range` and `Number of lines`.
    - Plot and export the profile.
    <figure markdown="span">
    ![Analyse GUI](../assets/images/analyse_GUI.png){ width="100%" }
    <figcaption>Graphical interface for defining transmural profiles.</figcaption>
    </figure>

    The generated profile will resemble:

    <figure markdown="span">
    ![Transmural profile](../assets/images/transmural_profile.png){ width="70%" }
    <figcaption>Example of a generated transmural profile.</figcaption>
    </figure>

### Visualizing Vector field

!!! note

    WRITE_VECTORS must be equal to True in the [Configuration file](../reference/configuration.md)

Once the structure tensor and eigenvectors have been calculated on the whole volume.

1. Use the `cardio-visualize-vector` command:
    ```console
    $ cardio-visualize-vector parameters_example.conf --start 120 --stride 6 
    ```

    - `--start`
    : Show vectors only from slice 120 to the last slice.

    - `--stride 6`
    : Show only 1 vector out of every 6.


2. The plot will appear and you will be able to rotate it by dragging left click

    <figure markdown="span">
    ![Vector visualization](../assets/images/vectors_example.png){ width="70%" }
    <figcaption>Vector field visualization from computed structure tensors using Fury.</figcaption>
    </figure>


### Generating and visualizing streamlines

1. Use the `cardio-visualize-vector` command to generate the streamlines:
    ```console
    $ cardio-generate-streamlines parameters_example.conf --seeds 10000 --start 150
    ```

    - `--start`
    : Show vectors only from slice 120 to the last slice.

    - `--seeds`
    : Number of seeds to start streamlines

    !!! note

        The streamlines are generated in ```output/streamlines.npz```


2. Use the `cardio-visualize-streamlines` command to plot the streamlines:
    ```console
    $ cardio-visualize-streamlines parameters_example.conf --line-width 1
    ```

    - `--line-width`
    : The width of the streamline in the plot.


<figure markdown="span">
![streamline visualization](../assets/images/streamline_example.png){ width="70%" }
<figcaption>Streamline visualization using Fury.</figcaption>
</figure>



## Notes

- The provided dataset is for demonstration purposes only.
- Modify parameters in `parameters_example.conf` (e.g., `SIGMA`, `RHO`) to suit your data.

