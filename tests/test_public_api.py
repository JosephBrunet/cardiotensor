from cardiotensor import (
    DataReader,
    calculate_intensities,
    calculate_structure_tensor,
    compute_azimuth_and_elevation,
    compute_fraction_anisotropy,
    compute_helical_and_intrusion_angles,
    compute_orientation,
    convert_to_8bit,
    export_vector_field_to_vtk,
    find_end_points,
    generate_streamlines_from_params,
    generate_streamlines_from_vector_field,
    load_npz_streamlines,
    load_trk_streamlines,
    plot_intensity,
    read_conf_file,
    save_intensity,
    write_spatialgraph_am,
)
from cardiotensor.analysis import (
    calculate_intensities as analysis_calculate_intensities,
)
from cardiotensor.analysis import find_end_points as analysis_find_end_points
from cardiotensor.analysis import plot_intensity as analysis_plot_intensity
from cardiotensor.analysis import save_intensity as analysis_save_intensity
from cardiotensor.orientation import (
    calculate_structure_tensor as orientation_calculate_structure_tensor,
)
from cardiotensor.orientation import (
    compute_azimuth_and_elevation as orientation_compute_azimuth_and_elevation,
)
from cardiotensor.orientation import (
    compute_fraction_anisotropy as orientation_compute_fraction_anisotropy,
)
from cardiotensor.orientation import (
    compute_helical_and_intrusion_angles as orientation_compute_helical_and_intrusion_angles,
)
from cardiotensor.orientation import (
    compute_orientation as orientation_compute_orientation,
)
from cardiotensor.orientation import (
    rotate_vectors_to_new_axis as orientation_rotate_vectors_to_new_axis,
)
from cardiotensor.tractography import (
    generate_streamlines_from_params as tractography_generate_streamlines_from_params,
)
from cardiotensor.tractography import (
    generate_streamlines_from_vector_field as tractography_generate_streamlines_from_vector_field,
)
from cardiotensor.utils import DataReader as utils_data_reader
from cardiotensor.utils import convert_to_8bit as utils_convert_to_8bit
from cardiotensor.utils import (
    export_vector_field_to_vtk as utils_export_vector_field_to_vtk,
)
from cardiotensor.utils import load_npz_streamlines as utils_load_npz_streamlines
from cardiotensor.utils import load_trk_streamlines as utils_load_trk_streamlines
from cardiotensor.utils import read_conf_file as utils_read_conf_file
from cardiotensor.utils import write_spatialgraph_am as utils_write_spatialgraph_am
from cardiotensor.visualization import (
    visualize_streamlines,
    visualize_vector_field,
)


def test_root_public_api_exports():
    assert DataReader is utils_data_reader
    assert calculate_intensities is analysis_calculate_intensities
    assert calculate_structure_tensor is orientation_calculate_structure_tensor
    assert compute_azimuth_and_elevation is orientation_compute_azimuth_and_elevation
    assert compute_fraction_anisotropy is orientation_compute_fraction_anisotropy
    assert compute_helical_and_intrusion_angles is orientation_compute_helical_and_intrusion_angles
    assert compute_orientation is orientation_compute_orientation
    assert convert_to_8bit is utils_convert_to_8bit
    assert export_vector_field_to_vtk is utils_export_vector_field_to_vtk
    assert find_end_points is analysis_find_end_points
    assert generate_streamlines_from_params is tractography_generate_streamlines_from_params
    assert (
        generate_streamlines_from_vector_field
        is tractography_generate_streamlines_from_vector_field
    )
    assert load_npz_streamlines is utils_load_npz_streamlines
    assert load_trk_streamlines is utils_load_trk_streamlines
    assert plot_intensity is analysis_plot_intensity
    assert read_conf_file is utils_read_conf_file
    assert save_intensity is analysis_save_intensity
    assert write_spatialgraph_am is utils_write_spatialgraph_am


def test_grouped_public_imports_expose_visualization_helpers():
    assert callable(visualize_streamlines)
    assert callable(visualize_vector_field)


def test_grouped_public_imports_expose_orientation_helper():
    assert callable(orientation_rotate_vectors_to_new_axis)
