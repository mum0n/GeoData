"""
Backend interface definitions for GeoData.

Each backend must implement these methods to be compatible with the GeoData API.
"""
module GeoDataInterfaces

using ..GeoDataCoreTypes
using Base: @kwdef

"""
    GeoBackendInterface

The complete interface a backend must implement.

Methods are grouped by capability - a backend only needs to implement
what it supports (checked via `BackendCapabilities`).
"""
abstract type GeoBackendInterface end

# ============================================================
# Core dataset operations
# ============================================================

"""
    backend_open(backend::GeoBackend, uri::String; mode::String = "r", kwargs...) -> GeoDataset

Open a dataset from the given URI.

# Arguments
- `backend`: Backend instance
- `uri`: Source URI (file path, URL, or backend-specific identifier)
- `mode`: "r" (read), "r+" (read/write), "w" (write/new), "a" (append)
- `kwargs`: Backend-specific options

# Returns
A `GeoDataset` with variables, coordinates, and metadata loaded.
"""
function backend_open end

"""
    backend_close(backend::GeoBackend, dataset::GeoDataset) -> Nothing

Close a dataset and release resources.
"""
function backend_close end

"""
    backend_create(backend::GeoBackend, uri::String; dims::Dict{Symbol, Dimension},
                   variables::Dict{String, <:AbstractArray}, coords::Dict{Symbol, <:AbstractArray},
                   crs::CoordinateSystem, attrs::Dict{String, Any}, kwargs...) -> GeoDataset

Create a new dataset at the given URI.

# Arguments
- `backend`: Backend instance
- `uri`: Destination URI
- `dims`: Dimension definitions
- `variables`: Data variables (name => array)
- `coords`: Coordinate variables
- `crs`: Coordinate reference system
- `attrs`: Global attributes
- `kwargs`: Backend-specific options (chunks, compression, etc.)
"""
function backend_create end

"""
    backend_write(backend::GeoBackend, dataset::GeoDataset; variables::Dict{String, <:AbstractArray} = Dict(),
                  coords::Dict{Symbol, <:AbstractArray} = Dict(), attrs::Dict{String, Any} = Dict(),
                  mode::String = "update", kwargs...) -> Nothing

Write data to an existing dataset.

# Arguments
- `backend`: Backend instance
- `dataset`: Open dataset (from `backend_open` with mode "r+" or "a")
- `variables`: Variables to write/update
- `coords`: Coordinates to write/update
- `attrs`: Global attributes to update
- `mode`: "update" (modify in place), "append" (append along unlimited dim)
- `kwargs`: Backend-specific options
"""
function backend_write end

# ============================================================
# Query and slicing operations
# ============================================================

"""
    backend_slice(backend::GeoBackend, dataset::GeoDataset; kwargs...) -> GeoDataset

Return a sliced view of the dataset.

# Arguments
- `backend`: Backend instance
- `dataset`: Source dataset
- `kwargs`: Dimension slices as `dim = (min, max)` or `dim = value` or `dim = [vals...]`

# Returns
New `GeoDataset` with sliced variables and coordinates.
"""
function backend_slice end

"""
    backend_select(backend::GeoBackend, dataset::GeoDataset; kwargs...) -> GeoDataset

Select specific coordinate values (exact match or nearest).

# Arguments
- `backend`: Backend instance
- `dataset`: Source dataset
- `kwargs`: Dimension selections as `dim = value` or `dim = [val1, val2...]`

# Returns
New `GeoDataset` with selected data.
"""
function backend_select end

"""
    backend_values(backend::GeoBackend, dataset::GeoDataset, varnames::Vector{String}; kwargs...) -> Dict{String, Any}

Query values at specific coordinate points.

# Arguments
- `backend`: Backend instance
- `dataset`: Source dataset
- `varnames`: Variables to query
- `kwargs`: Coordinate points as `lon = -65.0, lat = 45.0, depth = 50.0, time = ...`

# Returns
Dict mapping variable names to queried values (scalars or arrays).
"""
function backend_values end

# ============================================================
# Indexing and joining
# ============================================================

"""
    backend_index(backend::GeoBackend, dataset::GeoDataset; spatial::Bool = true, temporal::Bool = true) -> Any

Build or retrieve a spatial/temporal index for fast queries.

# Returns
Backend-specific index object (e.g., RTree, KDTree, interval tree).
"""
function backend_index end

"""
    backend_join(backend::GeoBackend, datasets::Vector{GeoDataset}; on::Vector{Symbol}, how::Symbol = :inner) -> GeoDataset

Join multiple datasets along shared dimensions.

# Arguments
- `backend`: Backend instance (from first dataset)
- `datasets`: Datasets to join
- `on`: Dimension names to join on (e.g., `[:time]`, `[:lon, :lat]`)
- `how`: Join type - `:inner`, `:outer`, `:left`, `:right`

# Returns
Joined `GeoDataset`.
"""
function backend_join end

# ============================================================
# Regridding and transformation
# ============================================================

"""
    backend_regrid(backend::GeoBackend, dataset::GeoDataset, target_grid::GeoDataset; method::Symbol = :bilinear,
                   vars::Vector{String} = String[]) -> GeoDataset

Regrid dataset to a target grid.

# Arguments
- `backend`: Backend instance
- `dataset`: Source dataset
- `target_grid`: Dataset defining target coordinates
- `method`: Interpolation method (`:bilinear`, `:nearest`, `:conservative`, `:bicubic`)
- `vars`: Variables to regrid (empty = all data variables)

# Returns
Regridded `GeoDataset`.
"""
function backend_regrid end

# ============================================================
# Introspection
# ============================================================

"""
    backend_info(backend::GeoBackend, dataset::GeoDataset) -> Dict{String, Any}

Return backend-specific information about the dataset.
"""
function backend_info end

"""
    backend_capabilities(backend::GeoBackend) -> BackendCapabilities

Return what this backend supports.
"""
function backend_capabilities end

# ============================================================
# Optional: Advanced operations
# ============================================================

"""
    backend_aggregate(backend::GeoBackend, dataset::GeoDataset; dim::Symbol, func::Function, kwargs...) -> GeoDataset

Aggregate along a dimension (mean, sum, min, max, etc.).
"""
function backend_aggregate end

"""
    backend_resample(backend::GeoBackend, dataset::GeoDataset; time_rule::String, func::Function) -> GeoDataset

Temporal resampling (e.g., daily to monthly).
"""
function backend_resample end

"""
    backend_subset(backend::GeoBackend, dataset::GeoDataset; bbox::Tuple, time_range::Tuple) -> GeoDataset

Fast spatial/temporal subsetting (uses index if available).
"""
function backend_subset end

export GeoBackendInterface,
    backend_open, backend_close, backend_create, backend_write,
    backend_slice, backend_select, backend_values,
    backend_index, backend_join, backend_regrid,
    backend_info, backend_capabilities,
    backend_aggregate, backend_resample, backend_subset

end # module GeoDataInterfaces