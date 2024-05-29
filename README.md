#  CardiacFiberTensor

Utils for calculating the orientation of fibers in a heart. Godspeed

## How to use:

Install all the dependencies

Create a new parameter file in `./param_file/` similar to `param_template.txt` and fill all the parameters

Run `processing` with and indicate the param file (run first with the option TEST activated to have a preview of the result)

Use `python3 processing.py --help` to look at the optional parameters


dependencies = [  
    "structure_tensor",  
    "dask_image>=2023.3.0",  
    "dask",  
    "numpy,  
    "opencv-python",  
    "matplotlib",  
    "cupy",  
]  