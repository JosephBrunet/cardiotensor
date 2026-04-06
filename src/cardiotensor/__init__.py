"""
cardiotensor — 3D cardiomyocyte orientation analysis from volumetric imaging.

Core usage
----------
>>> from cardiotensor import compute_orientation, DataReader
>>> reader = DataReader("path/to/images")
>>> compute_orientation("path/to/images", output_dir="output", sigma=3.0, rho=1.0)

Sub-module imports
------------------
For tractography::

    from cardiotensor.tractography import (
        generate_streamlines_from_vector_field,
        generate_streamlines_from_params,
    )

For I/O utilities::

    from cardiotensor.utils import (
        read_conf_file,
        load_npz_streamlines,
        load_trk_streamlines,
        write_spatialgraph_am,
        export_vector_field_to_vtk,
    )

For analysis::

    from cardiotensor.analysis import calculate_intensities, find_end_points

For visualization (requires a display)::

    from cardiotensor.visualization import visualize_streamlines, visualize_vector_field
"""

# --- Core ---
from cardiotensor.orientation.orientation_computation_pipeline import compute_orientation
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.utils import convert_to_8bit, read_conf_file

# --- Orientation ---
from cardiotensor.orientation.orientation_computation_functions import (
    calculate_structure_tensor,
    compute_fraction_anisotropy,
    compute_helix_and_transverse_angles,
)

# --- Tractography ---
from cardiotensor.tractography.generate_streamlines import (
    generate_streamlines_from_params,
    generate_streamlines_from_vector_field,
)

# --- I/O ---
from cardiotensor.utils.am_utils import write_spatialgraph_am
from cardiotensor.utils.streamlines_io_utils import (
    load_npz_streamlines,
    load_trk_streamlines,
)
from cardiotensor.utils.vector_vtk_export import export_vector_field_to_vtk

# --- Analysis ---
from cardiotensor.analysis.analysis_functions import (
    calculate_intensities,
    find_end_points,
    save_intensity,
)

__all__ = [
    # Core
    "compute_orientation",
    "DataReader",
    "read_conf_file",
    "convert_to_8bit",
    # Orientation
    "calculate_structure_tensor",
    "compute_fraction_anisotropy",
    "compute_helix_and_transverse_angles",
    # Tractography
    "generate_streamlines_from_vector_field",
    "generate_streamlines_from_params",
    # I/O
    "write_spatialgraph_am",
    "load_npz_streamlines",
    "load_trk_streamlines",
    "export_vector_field_to_vtk",
    # Analysis
    "calculate_intensities",
    "find_end_points",
    "save_intensity",
]
