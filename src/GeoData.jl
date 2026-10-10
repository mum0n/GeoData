"""
    GeoData.jl

Backend-agnostic geospatial data access layer.

Provides a uniform interface for loading, saving, querying, slicing, indexing, and joining
geospatial data across multiple backends (Zarr, GeoParquet, NetCDF, etc.).

# Design Principles

- **Backend-agnostic**: Same code works with any supported backend
- **Composable**: Operations chain naturally (load → slice → query → save)
- **Lazy where the format allows it**: the Zarr backend keeps its chunked arrays as the
  dataset's data, so opening a store reads nothing and indexing materialises only the
  elements it touches; NetCDF and GeoParquet materialise on open, because their Julia
  readers do
- **Coordinate-aware**: Understands lon/lat/depth/time dimensions and CRS
- **Extensible**: New backends via the `GeoBackend` interface
- **Explicit failure**: An unknown backend, URI, dimension, or missing coordinate is an error,
  never a silent default

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
| GeoParquet | `GeoParquet.jl` | `geoparquet://`, `.parquet` | Vector/tabular, Arrow ecosystem |
| NetCDF | `NCDatasets.jl` | `netcdf://`, `.nc` | Classic, CF-compliant |

# Core Operations

- `geoload(uri; backend, kwargs...)` - Load dataset
- `geosave(uri, data; backend, kwargs...)` - Save dataset
- `geoslice(data; dims...)` - Slice by coordinate ranges
- `geoselect(data; dims...)` - Select by coordinate values
- `geovalues(data; dims...)` - Interpolate/query at points
- `geojoin(datasets...; on, how)` - Join datasets
- `geoaggregate(ds; dim, func)` - Reduce along a dimension (mean, maximum, ...)
- `geoindex(data; dims...)` - Build spatial/temporal index
- `georegrid(data, target_grid; method)` - Regrid/interpolate
"""
module GeoData

using Reexport

# Core types and interfaces - include first to define submodules
    # Core types and the backend contract: included first to define the submodules.
    include("core/GeoDataTypes.jl")
    include("core/GeoDataBackends.jl")
    include("core/GeoDataCoordinates.jl")
    include("core/GeoDataOperations.jl")

# Bring core modules into scope
using .GeoDataTypes
using .GeoDataBackends
using .GeoDataCoordinates
using .GeoDataOperations

# Analysis over data already in hand: grid-to-grid interpolation, land/ocean tests, and
# bathymetry interrogation. Included after the core types so their signatures resolve.
    include("core/GeoDataOpenGrid.jl")
    include("core/GeoDataLandmask.jl")
    include("core/GeoDataBathy.jl")

# Registry - include and bring into scope
    include("backends/GeoDataRegistry.jl")
using .GeoDataRegistry

# Backend implementations
include("backends/common.jl")
include("backends/zarr_backend.jl")
include("backends/ncdatasets_backend.jl")
include("backends/geoparquet_backend.jl")

# Default backend capabilities, defined here beside the four implementations it covers.
# A default declared inside `GeoDataBackends` never registered on the function the
# backends extend, so an unlisted backend raised MethodError. Define it in this module,
# after the concrete methods, so the fallback is always present.
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
    geosave,
    geoslice,
    geoselect,
    geovalues,
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
    list_storage_groups

# Convenience re-exports from dependencies
# The backend readers stay private to this package. `@reexport`-ing them used to give
# every consumer Zarr, GeoParquet, NCDatasets and GeoInterface transitively, so a caller
# that never mentioned any of them still pulled them in and could not tell which names
# were ours. Name what you use.
using Zarr
using GeoParquet
using NCDatasets
using GeoInterface

end # module GeoData
