"""
    GeoData.jl

Backend-agnostic geospatial data access layer.

Provides a unified interface for loading, saving, querying, slicing, indexing, and joining
geospatial datasets across multiple storage backends (Zarr, GeoParquet, NetCDF, etc.).

# Design Principles

- **Backend-Agnostic**: Identical functional API across supported storage formats.
- **Composable**: Functional operations compose directly (`geoload` → `geoslice` → `geovalues` → `geosave`).
- **Lazy Evaluation**: Where supported by the backend (e.g. Zarr), chunked arrays are
  maintained lazily without premature memory materialization.
- **Coordinate-Aware**: Preserves dimension semantics (longitude, latitude, depth, time)
  and Coordinate Reference Systems (CRS).
- **Explicit Failure**: Missing coordinate metadata, unsupported formats, or dimension
  mismatches raise informative errors rather than silently coercing defaults.

# Quick Start

```julia
using GeoData

# Load from backend (inferred from URI scheme or extension)
ds = geoload("data/temperature.zarr")

# Slice by coordinate intervals
ds_slice = geoslice(ds, lon=(-70.0, -55.0), lat=(40.0, 50.0), depth=(0.0, 100.0))

# Query at specific coordinates
vals = geovalues(ds_slice, lon=-65.0, lat=45.0, depth=50.0)

# Save to destination backend
geosave("output/temperature_subset.zarr", ds_slice)
```

# Supported Backends

| Backend | Package | URI Scheme | Features |
|---------|---------|------------|----------|
| Zarr | `Zarr.jl` | `zarr://`, `.zarr` | Chunked, compressed, cloud-native storage |
| NetCDF/Zarr | `NCDatasets.jl` | `nczarr://`, `.nc` | CF-compliant Zarr mapping |
| GeoParquet | `GeoParquet.jl` | `geoparquet://`, `.parquet` | Columnar vector/tabular geometry |
| NetCDF | `NCDatasets.jl` | `netcdf://`, `.nc` | Standard CF-compliant NetCDF4 |
"""
module GeoData

using Reexport

# Core types and interfaces
include("core/GeoDataTypes.jl")
include("core/GeoDataBackends.jl")
include("core/GeoDataCoordinates.jl")
include("core/GeoDataOperations.jl")

using .GeoDataTypes
using .GeoDataBackends
using .GeoDataCoordinates
using .GeoDataOperations

# Pure computational analysis: grid interpolation, land/water classification, bathymetry
include("core/GeoDataOpenGrid.jl")
include("core/GeoDataLandmask.jl")
include("core/GeoDataBathy.jl")

# Registry
include("backends/GeoDataRegistry.jl")
using .GeoDataRegistry

# Backend implementations
include("backends/common.jl")
include("backends/zarr_backend.jl")
include("backends/ncdatasets_backend.jl")
include("backends/geoparquet_backend.jl")

# Default backend capability fallback
function backend_capabilities(backend::GeoBackend)
    return BackendCapabilities(read = true, write = true)
end

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
include("api/aggregate.jl")
include("api/storage.jl")
include("api/run_state.jl")



# Register default backends
function _register_default_backends()
    register_backend(:zarr, ZarrBackend())
    register_backend(:nczarr, NCDatasetsBackend(; format="nczarr"))
    register_backend(:ncdatasets, NCDatasetsBackend())
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
    geoget,
    register_source_provider!,
    unregister_source_provider!,
    list_source_providers,
    geosave,
    geoslice,
    geovalues,
    # A Zarr store created with its variables' final shape and no data in them, for a
    # caller that fills it a slice at a time.
    create_zarr_template,
    geojoin,
    geoaggregate,
    geoindex,
    georegrid,
      # Stats functions
      data_stats,
      variable_stats,
      geosummary,
      bounding_box,
      # Coordinate resolution: an absent axis is a fact about the dataset, not to invent
      axis,
      dim_permutation,
      normalize_depth,
      standardize_dimension_name,
      is_regular_grid,
      grid_spacing,
      find_coord_indices,
      slice_indices,
    # Backend registration
    register_backend,
    get_backend,
    list_backends,
    # Time parsing: every source reports "seconds since ..." in its own units
    parse_time_units,
    datetime_to_time,
    # Grid-to-grid interpolation and geographic geometry
    regrid_2d_field,
    regrid_3d_field,
    buffer_distance_to_degrees,
    expand_domain_with_buffer,
    point_in_polygon,
    is_point_on_land_geodata,
    is_marine_water_geodata,
    load_coastline_polygons,
    smooth_bathymetry,
    get_bathymetry_interpolator,
    extract_marine_cells,
    sample_marine_cells,
    is_cached,
    variables_like,
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
    list_storage_variables,
    list_storage_groups,
    lakehouse_dir,
    set_lakehouse_dir!,
    project_lakehouse_path,
    # Run state: what a run leaves behind, and how it is asked to clear it
    RunStateTarget,
    run_state_targets,
    confirm_run_state_purge,
    apply_run_state_purge,
    purge_run_state

# Explicit package imports
using Zarr
using GeoParquet
using NCDatasets
using GeoInterface

end # module GeoData
