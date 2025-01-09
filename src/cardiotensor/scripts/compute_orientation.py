import argparse

# from line_profiler import LineProfiler
import os
import sys
import time
from tkinter import Tk
from tkinter.filedialog import askopenfilename

from cardiotensor.orientation.orientation_computation_pipeline import (
    compute_orientation,
)
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import read_conf_file


def script() -> None:
    if len(sys.argv) < 2:
        Tk().withdraw()
        conf_file_path = askopenfilename(
            initialdir=f"{os.getcwd()}/param_files", title="Select file"
        )  # show an "Open" dialog box and return the path to the selected file
        if not conf_file_path:
            sys.exit("No file selected!")
        start_index = 0
        end_index = None
        use_gpu = True

    elif len(sys.argv) >= 2:
        parser = argparse.ArgumentParser(
            description="Convert images between tif and jpeg2000 formats."
        )
        parser.add_argument(
            "conf_file_path", type=str, help="Path to the input text file."
        )
        parser.add_argument(
            "--start_index", type=int, default=0, help="Starting index for processing."
        )
        parser.add_argument(
            "--end_index", type=int, default=None, help="Ending index for processing."
        )
        parser.add_argument("--gpu", action="store_true", help="Activate the gpu")
        args = parser.parse_args()
        conf_file_path = args.conf_file_path
        start_index = args.start_index
        end_index = args.end_index
        use_gpu = args.gpu

    try:
        params = read_conf_file(conf_file_path)
    except Exception as e:
        print(f"⚠️  Error reading parameter file: {conf_file_path}")
        print(f"\nError is {e}")
        sys.exit()

    (
        VOLUME_PATH,
        MASK_PATH,
        IS_FLIP,
        OUTPUT_DIR,
        OUTPUT_TYPE,
        SIGMA,
        RHO,
        N_CHUNK,
        PT_MV,
        PT_APEX,
        REVERSE,
        IS_TEST,
        N_SLICE_TEST,
    ) = (
        params[key]
        for key in [
            "IMAGES_PATH",
            "MASK_PATH",
            "FLIP",
            "OUTPUT_PATH",
            "OUTPUT_TYPE",
            "SIGMA",
            "RHO",
            "N_CHUNK",
            "POINT_MITRAL_VALVE",
            "POINT_APEX",
            "REVERSE",
            "TEST",
            "N_SLICE_TEST",
        ]
    )

    data_reader = DataReader(VOLUME_PATH)
    volume_shape = data_reader.shape

    # Set end_index to total_images if it's zero
    if end_index is None:
        end_index = volume_shape[0]

    if not IS_TEST:
        # If REVERSE is True, run the loop in reverse order; otherwise, run normally
        if REVERSE:
            for idx in range(end_index, start_index, -N_CHUNK):
                start_time = time.time()
                compute_orientation(
                    conf_file_path,
                    start_index=max(idx - N_CHUNK, 0),
                    end_index=idx,
                    use_gpu=use_gpu,
                )
                print("--- %s seconds ---" % (time.time() - start_time))
        else:
            for idx in range(start_index, end_index, N_CHUNK):
                start_time = time.time()
                compute_orientation(
                    conf_file_path,
                    start_index=idx,
                    end_index=min(idx + N_CHUNK, end_index),
                    use_gpu=use_gpu,
                )
                print("--- %s seconds ---" % (time.time() - start_time))

    else:
        start_time = time.time()
        compute_orientation(conf_file_path, start_index, end_index, use_gpu=use_gpu)
        print("--- %s seconds ---" % (time.time() - start_time))

    return None
