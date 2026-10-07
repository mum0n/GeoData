"""
Core type definitions for GeoData.

Generic n-dimensional array containers with coordinate metadata.
No geographic assumptions in the core - CRS defaults to cartesian.
"""
module GeoDataCoreTypes

using Base: @kwdef

"""
    DimensionType

Type of a dimension for generic handling.
"""
@enum DimensionType begin
    DIM_SPATIAL  # x, y, z (cartesian)
    DIM_TEMPORAL # t
    DIM_PARAMETRIC # parameter index (e.g., MCMC chain, draw)
    DIM_GENERIC  # any other dimension
end

"""
    Dimension

Describes a single dimension of a dataset.

# Fields
- `name::Symbol`: Dimension name (e.g., `:x`, `:y`, `:t`, `:chain`, `:draw`, `:param`)
- `size::Union{Int, Nothing}`: Size (number of points) or `nothing` if unknown/unbounded
- `coords::Union{AbstractVector, Nothing}`: Coordinate values if known
- `units::String`: Units (e.g., "m", "s", "K", "" for dimensionless)
- `standard_name::Union{String, Nothing}`: Optional standard name
- `dim_type::DimensionType`: Semantic type of dimension
- `calendar::Union{String, Nothing}`: Calendar for time dimensions
- `is_unlimited::Bool`: Whether this is an unlimited (record) dimension
"""
@kwdef struct Dimension
    name::Symbol
    size::Union{Int, Nothing} = nothing
    coords::Union{AbstractVector, Nothing} = nothing
    units::String = ""
    standard_name::Union{String, Nothing} = nothing
    dim_type::DimensionType = DIM_GENERIC
    calendar::Union{String, Nothing} = nothing
    is_unlimited::Bool = false
end

"""
    CoordinateSystem

Describes the coordinate reference system of a dataset.

Defaults to cartesian (empty CRS). For geographic data, use CRS like "EPSG:4326".

# Fields
- `crs::String`: PROJ string, WKT, EPSG code, or "" for cartesian
- `proj4::Union{String, Nothing}`: PROJ.4 string if available
- `wkt::Union{String, Nothing}`: WKT representation if available
"""
@kwdef struct CoordinateSystem
    crs::String = ""  # empty = cartesian
    proj4::Union{String, Nothing} = nothing
    wkt::Union{String, Nothing} = nothing
end

"""
    is_cartesian(crs::CoordinateSystem) -> Bool

True if coordinate system is cartesian (no geographic CRS).
"""
is_cartesian(crs::CoordinateSystem) = isempty(crs.crs)

"""
    GeoArray{T, N, A <: AbstractArray{T, N}}

An n-dimensional array with coordinate metadata.

Wraps an `AbstractArray` (may be lazy/chunked) and attaches
dimension and coordinate information.

# Type Parameters
- `T`: Element type
- `N`: Number of dimensions
- `A`: Wrapped array type

# Fields
- `data::A`: The underlying array (may be lazy)
- `dims::NTuple{N, Dimension}`: Dimension metadata
- `crs::CoordinateSystem`: Coordinate reference system (cartesian by default)
- `attrs::Dict{String, Any}`: Variable attributes (units, long_name, etc.)
"""
struct GeoArray{T, N, A <: AbstractArray{T, N}}
    data::A
    dims::NTuple{N, Dimension}
    crs::CoordinateSystem
    attrs::Dict{String, Any}

    function GeoArray{T, N, A}(data::A, dims::NTuple{N, Dimension},
                                crs::CoordinateSystem, attrs::Dict{String, Any}) where {T, N, A}
        length(dims) == N || error("Number of dimensions ($(length(dims))) must match array ndims ($N)")
        new{T, N, A}(data, dims, crs, attrs)
    end
end

# Constructor with inference
function GeoArray(data::A; dims::Vector{Dimension} = Dimension[],
                  crs::CoordinateSystem = CoordinateSystem(),
                  attrs::Dict{String, Any} = Dict{String, Any}()) where {T, N, A <: AbstractArray{T, N}}
    if isempty(dims)
        # Infer dimensions from array axes
        dims = [Dimension(name=Symbol("dim$i"), size=size(data, i)) for i in 1:N]
    end
    GeoArray{T, N, A}(data, Tuple(dims), crs, attrs)
end

# Positional constructor for compatibility
function GeoArray(data::AbstractArray, dims::NTuple{N, Dimension}, crs::CoordinateSystem, attrs::Dict{String, Any}) where N
    T = eltype(data)
    GeoArray{T, N, typeof(data)}(data, dims, crs, attrs)
end

Base.size(A::GeoArray) = size(A.data)
Base.ndims(A::GeoArray) = ndims(A.data)
Base.eltype(A::GeoArray) = eltype(A.data)
Base.getindex(A::GeoArray, i::Vararg{Int, N}) where N = A.data[i...]
Base.IndexStyle(::Type{<:GeoArray}) = IndexStyle(A.data)

"""
    GeoDataset

A collection of named `GeoArray` variables sharing coordinates.

Primary container for multi-variable datasets (NetCDF, Zarr, Parquet, etc.).

# Fields
- `variables::Dict{String, GeoArray}`: Named data variables
- `coords::Dict{Symbol, GeoArray}`: Coordinate variables (dimension coordinates)
- `dims::Dict{Symbol, Dimension}`: Shared dimension metadata
- `crs::CoordinateSystem`: Coordinate reference system (cartesian by default)
- `attrs::Dict{String, Any}`: Global attributes
- `backend::GeoBackend`: The backend that loaded this dataset
- `source::String`: Source URI/path
"""
struct GeoDataset
    variables::Dict{String, GeoArray}
    coords::Dict{Symbol, GeoArray}
    dims::Dict{Symbol, Dimension}
    crs::CoordinateSystem
    attrs::Dict{String, Any}
    backend::Any  # GeoBackend (avoiding circular dep)
    source::String
end

# Convenience constructors
function GeoDataset(; variables::Dict{String, GeoArray} = Dict{String, GeoArray}(),
                     coords::Dict{Symbol, GeoArray} = Dict{Symbol, GeoArray}(),
                     dims::Dict{Symbol, Dimension} = Dict{Symbol, Dimension}(),
                     crs::CoordinateSystem = CoordinateSystem(),
                     attrs::Dict{String, Any} = Dict{String, Any}(),
                     backend::Any = nothing,
                     source::String = "")
    GeoDataset(variables, coords, dims, crs, attrs, backend, source)
end

Base.getproperty(ds::GeoDataset, name::Symbol) = begin
    if hasfield(GeoDataset, name)
        getfield(ds, name)
    elseif haskey(ds.variables, string(name))
        ds.variables[string(name)]
    elseif haskey(ds.coords, name)
        ds.coords[name]
    else
        error("Property $name not found in GeoDataset")
    end
end

Base.propertynames(ds::GeoDataset, private::Bool = false) = begin
    names = fieldnames(GeoDataset)
    vars = collect(keys(ds.variables))
    coords = collect(keys(ds.coords))
    if private
        return (names..., vars..., coords...)
    end
    return (names..., vars..., coords...)
end

Base.hasproperty(ds::GeoDataset, name::Symbol) = hasfield(GeoDataset, name) ||
    haskey(ds.variables, string(name)) || haskey(ds.coords, name)

Base.keys(ds::GeoDataset) = (keys(ds.variables)..., keys(ds.coords)...)
Base.values(ds::GeoDataset) = (values(ds.variables)..., values(ds.coords)...)
Base.iterate(ds::GeoDataset, state = nothing) = iterate(pairs(merge(ds.variables, ds.coords)), state)

"""
    GeoBackend

Abstract type for backend implementations.

Each backend must implement the `GeoBackendInterface` methods.
"""
abstract type GeoBackend end

"""
    BackendCapabilities

Describes what a backend supports.
"""
@kwdef struct BackendCapabilities
    read::Bool = true
    write::Bool = true
    lazy::Bool = false      # Supports lazy/chunked operations
    chunked::Bool = false   # Supports chunked storage
    compression::Bool = false
    parallel::Bool = false  # Supports parallel I/O
    remote::Bool = false    # Can read from remote (HTTP/S3/etc.)
    append::Bool = false    # Supports appending to existing
    transactions::Bool = false
end

export DimensionType, DIM_SPATIAL, DIM_TEMPORAL, DIM_PARAMETRIC, DIM_GENERIC
export Dimension, CoordinateSystem, GeoArray, GeoDataset, GeoBackend, BackendCapabilities
export is_cartesian

end # module GeoDataCoreTypes