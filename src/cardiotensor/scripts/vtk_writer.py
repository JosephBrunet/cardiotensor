import argparse

from cardiotensor.export.vtk_writer import vtk_writer


def script() -> None:
    parser = argparse.ArgumentParser(
        description="Convert images between tif and jpeg2000 formats."
    )
    parser.add_argument("conf_file_path", type=str, help="Path to the input text file.")
    parser.add_argument(
        "--start_index", type=int, default=0, help="Start index for volume subset."
    )
    parser.add_argument(
        "--end_index", type=int, default=None, help="End index for volume subset."
    )
    parser.add_argument(
        "--bin", type=int, default=1, help="Binning factor for volume reduction."
    )
    args = parser.parse_args()

    conf_file_path = args.conf_file_path
    bin_factor = args.bin
    start_index = args.start_index
    end_index = args.end_index

    vtk_writer(conf_file_path, bin_factor, start_index, end_index)
