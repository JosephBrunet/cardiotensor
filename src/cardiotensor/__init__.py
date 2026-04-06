"""
cardiotensor — 3D cardiomyocyte orientation analysis from volumetric imaging.

Recommended imports
-------------------
>>> from cardiotensor import DataReader, compute_orientation, read_conf_file
>>> reader = DataReader("path/to/images")
>>> compute_orientation("path/to/images", output_dir="output", sigma=1.0, rho=3.0)

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

from cardiotensor.analysis import (
    calculate_intensities,
    find_end_points,
    plot_intensity,
    save_intensity,
)
from cardiotensor.orientation import (
    calculate_structure_tensor,
    compute_azimuth_and_elevation,
    compute_fraction_anisotropy,
    compute_helix_and_transverse_angles,
    compute_orientation,
)
from cardiotensor.tractography import (
    generate_streamlines_from_params,
    generate_streamlines_from_vector_field,
)
from cardiotensor.utils import (
    DataReader,
    convert_to_8bit,
    export_vector_field_to_vtk,
    load_npz_streamlines,
    load_trk_streamlines,
    read_conf_file,
    write_spatialgraph_am,
)

__all__ = [
    # Core
    "compute_orientation",
    "DataReader",
    "read_conf_file",
    "convert_to_8bit",
    # Orientation
    "calculate_structure_tensor",
    "compute_azimuth_and_elevation",
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
    "plot_intensity",
    "save_intensity",
]
