"""
    data.jl

Data fetching and manipulation operations for common geospatial datasets.
"""
module Data

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using GeoData.GeoDataCoordinates: bounding_box, slice_indices
using Downloads
using NCDatasets
using JSON3
using HTTP
using Base64
using Statistics
using LinearAlgebra
using Interpolations
using Dates
using Random

# Re-export core functions
export
    geoload,
    geosave,
    geoslice,
    geoselect,
    geovalues,
    geonearest,
    geoindex,
    georegrid,
    georegrid_regular,
    geosummary,
    variable_stats,
    bounding_box,
    data_stats

# Data-specific exports
export
    # Coastline
    fetch_natural_earth_coastline,
    load_coastline_geodata,
    load_coastline_polygons_geodata,
    is_point_on_land_geodata,
    is_marine_water_geodata,
    # Bathymetry (ETOPO, ERDDAP)
    fetch_erddap_bathymetry,
    load_bathymetry_geodata,
    save_bathymetry_geodata,
    get_bathymetry_interpolator,
    regrid_bathymetry_from_etopo,
    etopo_bathymetry_field,
    # Surface winds (Open-Meteo ERA5)
    fetch_open_meteo_winds,
    load_wind_stress_geodata,
    build_bulk_surface_flux_geodata,
    wind_speed_to_kinematic_stress,
    wind_speed_from_stress,
    # WOA23 climatology
    fetch_woa23,
    load_woa23_interpolators,
    # Boundary hydrography
    fetch_boundary_hydrography_geodata,
    build_boundary_tracer_interpolators_geodata,
    # Regridding utilities
    regrid_2d_field,
    regrid_3d_field,
    slice_bathymetry_geodata,
    slice_wind_geodata,
    extract_grid_coordinates_geodata,
    # Bathymetry processing
    smooth_bathymetry,
    extract_marine_cells,
    sample_marine_coordinates,
    # Geospatial utilities
    buffer_distance_to_degrees,
    expand_domain_with_buffer,
    # Copernicus Marine and CDS
    copernicusmarine_executable,
    project_python,
    copernicus_credentials,
    copernicus_login_reminder,
    fetch_copernicus_physics_subset,
    fetch_copernicus_hydrography_with_fallback,
    fetch_copernicus_surface_winds

include("coastline.jl")
include("bathymetry.jl")
include("winds.jl")
include("woa23.jl")
include("boundary.jl")
include("regrid_utils.jl")
include("manifest.jl")
include("copernicus.jl")
include("regional_cube.jl")
using .GeoDataManifest
using .Copernicus

# Re-export GeoDataManifest exports
export DATA_SOURCES, DataSource, fetch_input, input_dir, file_digest, data_provenance, data_source, describe_data_sources

# Re-export regional cube ingestion symbols
export RegionalCubeConfig,
       load_cube_config,
       standard_ocean_depths,
       assimilate_regional_cube,
       build_cube_hsi_evaluator,
       evaluate_bioenergetic_scope

end # module Data