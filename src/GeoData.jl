"""
    GeoData.jl

Backend-agnostic geospatial data access layer.

Provides a uniform interface for loading, saving, querying, slicing, indexing, and joining
geospatial data across multiple backends (Zarr, YAXArrays, GeoParquet, NetCDF, etc.).

# Design Principles

- **Backend-agnostic**: Same code works with any supported backend
- **Composable**: Operations chain naturally (load → slice → query → save)
- **Lazy by default**: Operations return views/iterators until materialized
- **Coordinate-aware**: Understands lon/lat/depth/time dimensions and CRS
- **Extensible**: New backends via `GeoDataBackend` interface

# Quick Start

```julia
using GeoData

# Load from any backend (auto-detected from path/URI)
ds = geoload("data/temperature.zarr")

# Slice by coordinate ranges
ds_slice = geoslice(ds, lon=(-70, -55), lat=(40, 50), depth=(0, 100))

# Query at specific points
vals = geovalues(ds_slice, lon=-65.0, lat=45.0, depth=50.0)

# Save to any backend
geosave("output/temperature_subset.zarr", ds_slice)
```

# Supported Backends

| Backend | Package | URI Scheme | Strengths |
|---------|---------|------------|-----------|
| Zarr | `Zarr.jl` | `zarr://`, `.zarr` | Chunked, compressed, cloud-native |
| NetCDF/Zarr | `NCDatasets.jl` | `nczarr://`, `.nc` | CF conventions, interoperability |
| YAXArrays | `YAXArrays.jl` | `yaxarray://` | Lazy, labeled arrays, Dask-like |
| GeoParquet | `GeoParquet.jl` | `geoparquet://`, `.parquet` | Vector/tabular, Arrow ecosystem |
| NetCDF | `NCDatasets.jl` | `netcdf://`, `.nc` | Classic, CF-compliant |

# Core Operations

- `geoload(uri; backend, kwargs...)` - Load dataset
- `geosave(uri, data; backend, kwargs...)` - Save dataset
- `geoslice(data; dims...)` - Slice by coordinate ranges
- `geoselect(data; dims...)` - Select by coordinate values
- `geovalues(data; dims...)` - Interpolate/query at points
- `geojoin(datasets...; on, how)` - Join datasets
- `geoindex(data; dims...)` - Build spatial/temporal index
- `georegrid(data, target_grid; method)` - Regrid/interpolate
"""
module GeoData

using Reexport

# Core types and interfaces - include first to define submodules
include("core/types.jl")
include("core/interfaces.jl")
include("core/coordinates.jl")
include("core/operations.jl")

# Bring core modules into scope
using .GeoDataCoreTypes
using .GeoDataInterfaces
using .GeoDataCoordinates
using .GeoDataOperations

# Registry - include and bring into scope
include("backends/registry.jl")
using .GeoDataRegistry

# Backend implementations
include("backends/abstract_backend.jl")
include("backends/zarr_backend.jl")
include("backends/ncdatasets_backend.jl")
include("backends/yaxarray_backend.jl")
include("backends/geoparquet_backend.jl")

# High-level API
include("api/utils.jl")
include("api/load.jl")
include("api/save.jl")
include("api/slice.jl")
include("api/select.jl")
include("api/values.jl")
include("api/join.jl")
include("api/index.jl")
include("api/regrid.jl")
include("api/storage.jl")

# Data modules
include("data/data.jl")
using .Data
# Import manifest functions for re-export
import .Data: DATA_SOURCES, DataSource, fetch_input, input_dir, file_digest, data_provenance, data_source, describe_data_sources

# Data catalog and lakehouse
include("data/catalog.jl")
using .CatalogModule

include("data/lakehouse.jl")
using .LakehouseModule


# Register default backends
function _register_default_backends()
    register_backend(:zarr, ZarrBackend())
    register_backend(:nczarr, NCDatasetsBackend(; format="nczarr"))
    register_backend(:ncdatasets, NCDatasetsBackend())
    register_backend(:yaxarray, YAXArraysBackend())
    register_backend(:geoparquet, GeoParquetBackend())
end
_register_default_backends()

# Re-export public API
export
    # Types
    GeoDataset,
    GeoArray,
    GeoBackend,
    CoordinateSystem,
    Dimension,
    is_cartesian,
    # Core functions
    geoload,
    geosave,
    geoslice,
    geoselect,
    geovalues,
    geojoin,
    geoindex,
    georegrid,
    # Stats functions
    data_stats,
    variable_stats,
    geosummary,
    bounding_box,
    # Backend registration
    register_backend,
    get_backend,
    list_backends,
    # Utilities
    infer_backend,
    open_dataset,
    close_dataset,
    # Storage abstraction
    GeoStorage,
    open_geostorage,
    close_geostorage,
    create_storage_group,
    has_storage_group,
    write_storage_variable!,
    read_storage_variable,
    # Data functions (from Data module)
    # Coastline
    fetch_natural_earth_coastline,
    load_coastline_geodata,
    load_coastline_polygons_geodata,
    is_point_on_land_geodata,
    is_marine_water_geodata,
    # Bathymetry
    fetch_erddap_bathymetry,
    load_bathymetry_geodata,
    save_bathymetry_geodata,
    get_bathymetry_interpolator,
    regrid_bathymetry_from_etopo,
    etopo_bathymetry_field,
    # Winds
    fetch_open_meteo_winds,
    load_wind_stress_geodata,
    build_bulk_surface_flux_geodata,
    wind_speed_to_kinematic_stress,
    wind_speed_from_stress,
    # WOA23
    fetch_woa23,
    load_woa23_interpolators,
    # Boundary
    fetch_boundary_hydrography_geodata,
    build_boundary_tracer_interpolators_geodata,
    # Regridding
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
    # Copernicus Marine & Climate (GLORYS / ERA5)
    copernicusmarine_executable,
    project_python,
    copernicus_credentials,
    copernicus_login_reminder,
    fetch_copernicus_physics_subset,
    fetch_copernicus_hydrography_with_fallback,
    fetch_copernicus_surface_winds,
    # Manifest (provenance registry)
    DATA_SOURCES,
    DataSource,
    data_source,
    describe_data_sources,
    fetch_input,
    input_dir,
    file_digest,
    data_provenance,
     # Data catalog and lakehouse
     GeoDataCatalog,
     DatasetEntry,
     CatalogStats,
     default_catalog_path,
     lakehouse_root_dir,
     lakehouse_path,
     lakehouse_tiers,
     compute_lakehouse_checksum,
     geopublish_dataset!,
     geofetch_dataset,
     load_catalog!,
     save_catalog,
     register_dataset!,
     unregister_dataset!,
     update_dataset!,
     geopublish!,
     geofetch,
     get_dataset,
     get_dataset_by_key,
     get_datasets_by_name,
     find_datasets,
     get_best_dataset,
     access_dataset!,
     record_download!,
     get_stats,
     export_catalog_to_csv,
     get_global_catalog,
     init_global_catalog!,
     # Regional data cube assimilation
     RegionalCubeConfig,
     load_cube_config,
     standard_ocean_depths,
     assimilate_regional_cube,
     build_cube_hsi_evaluator,
     evaluate_bioenergetic_scope

# Convenience re-exports from dependencies
@reexport using Zarr, GeoParquet, NCDatasets, GeoInterface

end # module GeoData