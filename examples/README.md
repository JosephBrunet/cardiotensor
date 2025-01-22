# Examples

This directory contains example scripts and datasets to help you understand and utilize the functionality of the **cardiotensor** package. The examples demonstrate how to load input data, compute fiber orientations, and visualize the results.

## 1. Data

### Heart Volume

- **`./data/635.2um_LADAF-2021-17_heart_overview_/`**

This folder contains a small 3D image volume designed for testing and experimentation. It is a cropped and downsampled dataset that maintains key features while reducing file size.

The original heart data comes from a synchrotron X-ray scan of a complete human heart. The dataset's original voxel size is 19.85 µm/voxel. To accommodate GitHub's file size limitations, the 3D volume was downsampled by a factor of 32. The full-resolution dataset is available at [Human Organ Atlas](https://human-organ-atlas.esrf.fr/datasets/1659197537).

### Binary Mask

- **`./data/mask/`**

This folder contains the corresponding binary mask used to segment the heart from the background.

## 2. Configuration File

The configuration file `parameters_example.conf` is pre-filled and ready to execute. It contains all necessary parameters for running both test and full-volume processing workflows. You can modify it to suit your specific requirements.

## 3. Running the Examples

### Installation

1. Clone the repository and install the package:

   ```bash
   git clone https://github.com/JosephBrunet/cardiotensor.git
   cd cardiotensor
   pip install .
   ```

2. Navigate to the `examples` directory:

   ```bash
   cd examples
   ```

### Test Slice

To process and visualize a single test slice:

1. Open the configuration file `parameters_example.conf` and in the `[TEST]` section set:
   **`TEST = True`**

2. Run the following command:

   ```bash
   cardio-tensor ./parameters_example.conf
   ```

3. The processed result will be displayed as a plot:

   ![Example Result](../assets/images/result_test_slice.png)

### Processing the Whole Volume

To process the entire volume:

1. Open the configuration file `parameters_example.conf` and in the `[TEST]` section set:
   **`TEST = False`**

2. Run the following command:

   ```bash
   cardio-tensor ./parameters_example.conf
   ```

3. The results will be saved in the `./output` directory with the following structure:
   ```
   ./output
   ├── HA   # Helix angle results
   ├── IA   # Intrusion angle results
   └── FA   # Fractional anisotropy results
   ```

### Analysis

After processing the whole volume, you can plot a transmural profile by running:

```bash
cardio-analysis ./parameters_example.conf 150
```

Here, `150` refers to the slice number. The GUI will appear as follows:

![Analysis GUI](../assets/images/analyse_GUI.png)

1. Drag and drop on the image to define a transmural profile line.
2. Adjust parameters like `Angle range` and `Number of lines` as needed.
3. Click on `Plot profile` to generate the profile plot. The result will resemble this:

![Transmural Plot Profile](../assets/images/transmural_profile.png)

You can save the plot or click `Save Profile` to export the profile as a `.csv` file.

## Notes

- The dataset in `./data/635.2um_LADAF-2021-17_heart_overview_/` is intended for demonstration only and does not represent large-scale datasets.
- Modify parameters in the configuration file (e.g., `SIGMA`, `RHO`, `N_CHUNK`) to optimize processing for your data.

## Contributing

We welcome contributions! If you have suggestions for new examples or encounter any issues, please open a pull request or create an issue on the repository.
