import math
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from skimage.measure import block_reduce

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import read_conf_file


def writeStructuredVTK(
    aspectRatio: list[float] = [1.0, 1.0, 1.0],
    origin: list[float] = [0.0, 0.0, 0.0],
    cellData: dict[str, np.ndarray] = {},
    pointData: dict[str, np.ndarray] = {},
    fileName: str = "spam.vtk",
) -> None:
    """
    Write a plain text regular grid VTK file from 3D or 4D arrays.

    Args:
        aspectRatio (list[float]): Length between two nodes in every direction.
            Default is [1.0, 1.0, 1.0].
        origin (list[float]): Origin of the grid. Default is [0.0, 0.0, 0.0].
        cellData (dict[str, np.ndarray]): Cell fields; 3D arrays for scalar fields, 4D arrays for vector fields.
        pointData (dict[str, np.ndarray]): Nodal fields interpolated by Paraview.
        fileName (str): Name of the output file. Default is 'spam.vtk'.

    Returns:
        None
    """

    dimensions = []

    # Check dimensions
    if len(cellData) + len(pointData) == 0:
        print(f"spam.helpers.writeStructuredVTK() Empty files. Not writing {fileName}")
        return

    if len(cellData):
        dimensionsCell = list(cellData.values())[0].shape[:3]
        for k, v in cellData.items():
            if set(dimensionsCell) != set(v.shape[:3]):
                print(
                    f"spam.helpers.writeStructuredVTK() Inconsistent cell field sizes {dimensionsCell} != {v.shape[:3]}"
                )
                return
        dimensions = [n + 1 for n in dimensionsCell]

    if len(pointData):
        dimensionsPoint = list(pointData.values())[0].shape[:3]
        for k, v in pointData.items():
            if set(dimensionsPoint) != set(v.shape[:3]):
                print(
                    f"spam.helpers.writeStructuredVTK() Inconsistent point field sizes {dimensionsPoint} != {v.shape[:3]}"
                )
                return
        dimensions = dimensionsPoint

    if len(cellData) and len(pointData):
        if {n + 1 for n in dimensionsCell} != set(dimensionsPoint):
            print(
                "spam.helpers.writeStructuredVTK() Inconsistent point VS cell field sizes.\
                 Point size should be +1 for each axis."
            )

    with open(fileName, "w") as f:
        # header
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"VTK file from spam: {fileName}\n")
        f.write("ASCII\n\n")
        f.write("DATASET STRUCTURED_POINTS\n")

        f.write("DIMENSIONS {} {} {}\n".format(*reversed(dimensions)))
        f.write("ASPECT_RATIO {} {} {}\n".format(*reversed(aspectRatio)))
        f.write("ORIGIN {} {} {}\n\n".format(*reversed(origin)))

        # pointData
        if len(pointData) == 1:
            f.write(f"POINT_DATA {dimensions[0] * dimensions[1] * dimensions[2]}\n\n")
            _writeFieldInVtk(pointData, f)
        elif len(pointData) > 1:
            f.write(f"POINT_DATA {dimensions[0] * dimensions[1] * dimensions[2]}\n\n")
            for k in pointData:
                _writeFieldInVtk({k: pointData[k]}, f)

        # cellData
        if len(cellData) == 1:
            f.write(
                f"CELL_DATA {(dimensions[0] - 1) * (dimensions[1] - 1) * (dimensions[2] - 1)}\n\n"
            )
            _writeFieldInVtk(cellData, f)
        elif len(cellData) > 1:
            f.write(
                f"CELL_DATA {(dimensions[0] - 1) * (dimensions[1] - 1) * (dimensions[2] - 1)}\n\n"
            )
            for k in cellData:
                _writeFieldInVtk({k: cellData[k]}, f)

        f.write("\n")


def _writeFieldInVtk(data: dict[str, np.ndarray], f: Any, flat: bool = False) -> None:
    """
    Helper function to write fields into a VTK file.

    Args:
        data (dict[str, np.ndarray]): Data fields to write.
        f (Any): File object to write to.
        flat (bool): Whether to flatten the data before writing. Default is False.

    Returns:
        None
    """

    for key in data:
        field = data[key]

        if flat:
            # SCALAR flatten (n by 1)
            if len(field.shape) == 1:
                f.write("SCALARS {} float\n".format(key.replace(" ", "_")))
                f.write("LOOKUP_TABLE default\n")
                for item in field:
                    f.write(f"    {item}\n")
                f.write("\n")

            # VECTORS flatten (n by 3)
            elif len(field.shape) == 2 and field.shape[1] == 3:
                f.write("VECTORS {} float\n".format(key.replace(" ", "_")))
                for item in field:
                    f.write("    {} {} {}\n".format(*reversed(item)))
                f.write("\n")

        else:
            # SCALAR not flatten (n1 by n2 by n3)
            if len(field.shape) == 3:
                f.write("SCALARS {} float\n".format(key.replace(" ", "_")))
                f.write("LOOKUP_TABLE default\n")
                for item in field.reshape(-1):
                    f.write(f"    {item}\n")
                f.write("\n")

            # VECTORS (n1 by n2 by n3 by 3)
            elif len(field.shape) == 4 and field.shape[3] == 3:
                f.write("VECTORS {} float\n".format(key.replace(" ", "_")))
                for item in field.reshape(
                    (field.shape[0] * field.shape[1] * field.shape[2], field.shape[3])
                ):
                    f.write("    {} {} {}\n".format(*reversed(item)))
                f.write("\n")

            # TENSORS (n1 by n2 by n3 by 3 by 3)
            elif len(field.shape) == 5 and field.shape[3] * field.shape[4] == 9:
                f.write("TENSORS {} float\n".format(key.replace(" ", "_")))
                for item in field.reshape(
                    (
                        field.shape[0] * field.shape[1] * field.shape[2],
                        field.shape[3] * field.shape[4],
                    )
                ):
                    f.write(
                        "    {} {} {}\n    {} {} {}\n    {} {} {}\n\n".format(
                            *reversed(item)
                        )
                    )
                f.write("\n")
            else:
                print(
                    "spam.helpers.vtkio._writeFieldInVtk(): I'm in an unknown condition!"
                )


def vtk_writer(
    conf_file_path: str,
    bin_factor: int = 1,
    start_index: int | None = None,
    end_index: int | None = None,
) -> None:
    """
    Process volume data and write results to a VTK file.

    Args:
        conf_file_path (str): Path to the configuration file.
        bin_factor (int, optional): Binning factor for data reduction. Default is 1.
        start_index (Optional[int], optional): Starting index for processing. Default is None.
        end_index (Optional[int], optional): Ending index for processing. Default is None.

    Returns:
        None
    """
    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(f"⚠️  Error reading parameter file '{conf_file_path}': {e}")
        sys.exit(1)

    # Extracting parameters safely using .get() with defaults where necessary
    VOLUME_PATH = params.get("IMAGES_PATH", "")
    OUTPUT_DIR = params.get("OUTPUT_PATH", "./output")
    OUTPUT_DIR = Path(OUTPUT_DIR)

    data_reader_volume = DataReader(VOLUME_PATH)

    _, h, w = data_reader_volume.shape

    if start_index is None:
        start_index = 0
    if end_index is None:
        end_index = data_reader_volume.shape[0]

    output_npy = OUTPUT_DIR / "eigen_vec"
    npy_list = sorted(list(output_npy.glob("*.npy")))[start_index:end_index]

    shape = (end_index - start_index, h, w)
    vector_field = np.empty((3,) + shape)

    blocks = [npy_list[i : i + bin_factor] for i in range(0, len(npy_list), bin_factor)]

    bin_array = np.empty(
        (3, len(blocks), math.ceil(h / bin_factor), math.ceil(w / bin_factor))
    )
    # bin_array = np.empty((3, len(blocks),h, w))
    # Load and assign data in chunks based on bin factor
    for i, b in enumerate(blocks):
        print(f"Processing block: {i}/{len(blocks)}")
        array = np.empty((3, len(b), h, w))
        for idx, p in enumerate(b):
            print(f"Reading file: {p.name}")

            # Load the numpy data
            array[:, idx, :, :] = np.load(
                p
            )  # Shape should match the expected volume slice

        array = array.mean(axis=1)

        # Define the block size for each axis
        block_size_vec = (bin_factor, bin_factor)

        # Use block_reduce to bin the volume
        bin_array[0, i, :, :] = block_reduce(
            array[0, :, :], block_size=block_size_vec, func=np.mean
        )
        bin_array[1, i, :, :] = block_reduce(
            array[1, :, :], block_size=block_size_vec, func=np.mean
        )
        bin_array[2, i, :, :] = block_reduce(
            array[2, :, :], block_size=block_size_vec, func=np.mean
        )

        # bin_array[:,i,:,:] = array[:,0,:,:]

    vector_field = bin_array
    shape = bin_array.shape[1:]

    # Check where the z-component (index 2) is negative
    negative_z_mask = vector_field[0, :, :, :] < 0

    # Flip the vectors where the z-component is negative
    vector_field[:, negative_z_mask] *= -1

    mask_volume = np.where(np.isnan(vector_field[0, :, :, :]), 0, 1)

    # for i in range(0,3):
    #     mask_volume[vector_field[i, :, :, :] == 0] = 0

    mask_volume = mask_volume.astype(np.uint8)

    # for i in range(0,3):
    #     vector_field[i,:,:,:][vector_field[0,:,:,:] <= 0] = 0

    # #---------------------------------------------
    # # Binning
    # print("Binning...")

    # # Define the binning factor
    # bin_factor = 32  # Adjust this as needed

    # # Calculate new dimensions by cropping to the nearest multiple of bin_factor
    # z_new = vector_field.shape[1] - (vector_field.shape[1] % bin_factor)
    # y_new = vector_field.shape[2] - (vector_field.shape[2] % bin_factor)
    # x_new = vector_field.shape[3] - (vector_field.shape[3] % bin_factor)

    # # Crop the array if necessary to make dimensions multiples of bin_factor
    # vector_field_cropped = vector_field[:, :z_new, :y_new, :x_new]

    # # Reshape and bin by averaging
    # vector_field_binned = vector_field_cropped.reshape(3, z_new // bin_factor, bin_factor, y_new // bin_factor, bin_factor, x_new // bin_factor, bin_factor).mean(axis=(2, 4, 6))

    # vector_field = vector_field_binned
    # shape = vector_field_binned.shape[1:]

    # ---------------------------------------------
    # HA

    output_HA = OUTPUT_DIR / "HA"
    data_reader_HA = DataReader(output_HA)
    HA_volume = data_reader_HA.load_volume(start_index=start_index, end_index=end_index)
    # mask_volume = np.where(HA_volume == 0, 0, 1)

    block_size = (bin_factor, bin_factor, bin_factor)

    # Use block_reduce to bin the volume
    HA_volume = block_reduce(HA_volume, block_size=block_size, func=np.mean)

    # mask_volume = block_reduce(mask_volume, block_size=block_size, func=np.mean)
    # mask_volume = np.where(mask_volume < 0.5, 0, 1)

    # HA_volume = HA_volume *90/255 - 90

    cellData = {}
    # cellData["eigenVectors"] = vector_field.reshape((shape[0], h, w, 3))
    cellData["eigenVectors"] = np.moveaxis(vector_field, 0, -1)

    # cellData["eigenVectors"] = vector_field.reshape((shape[0], shape[1], shape[2], 3))
    cellData["HA_angles"] = HA_volume.reshape(shape)
    cellData["mask"] = mask_volume.reshape(shape)

    # Overwrite nans and infs with 0, rubbish I know
    cellData["eigenVectors"][np.logical_not(np.isfinite(cellData["eigenVectors"]))] = 0
    # cellData["eigenVectors"] = np.nan_to_num(cellData["eigenVectors"])

    cellData["HA_angles"][np.logical_not(np.isfinite(cellData["HA_angles"]))] = 0
    cellData["mask"][np.logical_not(np.isfinite(cellData["mask"]))] = 0

    # print("eigenVectors shape:", cellData["eigenVectors"].shape)
    # print("HA_angles shape:", cellData["HA_angles"].shape)
    # print("mask shape:", cellData["mask"].shape)

    try:
        for idx in range(0, cellData["eigenVectors"].shape[0]):
            fig, axes = plt.subplots(2, 3, figsize=(15, 5))
            # Plot each slice in a separate subplot
            axes[0, 0].imshow(cellData["eigenVectors"][idx, :, :, 0])
            axes[0, 1].imshow(cellData["eigenVectors"][idx, :, :, 1])
            axes[0, 2].imshow(cellData["eigenVectors"][idx, :, :, 2])
            axes[1, 0].imshow(cellData["mask"][idx, :, :])
            axes[1, 1].imshow(cellData["HA_angles"][idx])
            plt.show()

    except:
        print("\n⚠ C'ant plot graph\n")

    vtf_name = OUTPUT_DIR / "paraview.vtk"
    print(f"Writing the .vtk file: {vtf_name}")
    writeStructuredVTK(cellData=cellData, fileName=vtf_name)
