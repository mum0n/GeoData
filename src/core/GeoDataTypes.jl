"""
Core type definitions for GeoData.

Generic n-dimensional array containers with coordinate metadata.
No geographic assumptions in the core - CRS defaults to cartesian.
"""
module GeoDataTypes

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
- `units::Union{String, Nothing}`: Units, or `nothing` when unknown (distinct from
  dimensionless, which is `"1"`)
- `standard_name::Union{String, Nothing}`: Optional standard name
- `dim_type::DimensionType`: Semantic type of dimension
- `calendar::Union{String, Nothing}`: Calendar for time dimensions
- `is_unlimited::Bool`: Whether this is an unlimited (record) dimension
"""
@kwdef struct Dimension
    name::Symbol
    size::Union{Int, Nothing} = nothing
    coords::Union{AbstractVector, Nothing} = nothing
    units::Union{String, Nothing} = nothing
    standard_name::Union{String, Nothing} = nothing
    dim_type::DimensionType = DIM_GENERIC
    calendar::Union{String, Nothing} = nothing
    is_unlimited::Bool = false
end

"""
    Dimension(existing::Dimension; kwargs...)

Copy `existing` with the given fields replaced. Needed by every operation that narrows
an axis: the units, calendar, and semantic type must survive.
"""
function Dimension(d::Dimension; kwargs...)
    base = NamedTuple(f => getfield(d, f) for f in fieldnames(Dimension))
    return Dimension(; base..., kwargs...)
end

"""
    Dimension(name::Symbol, size::Int; kwargs...)

Convenience constructor for tests and ad-hoc datasets; the remaining fields take their
defaults.
"""
Dimension(name::Symbol, size::Int; kwargs...) =
    Dimension(; name = name, size = size, kwargs...)

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
#
# `attrs` accepts any dictionary of attributes, whose keys may be strings or symbols:
# every backend spells attributes differently, and the stored form is Dict{String, Any}.
function GeoArray(data::A; dims::Vector{Dimension} = Dimension[],
                  crs::CoordinateSystem = CoordinateSystem(),
                  attrs = Dict{String, Any}()) where {T, N, A <: AbstractArray{T, N}}
    if isempty(dims)
        # Infer dimensions from array axes
        dims = [Dimension(name=Symbol("dim$i"), size=size(data, i)) for i in 1:N]
    end
    GeoArray{T, N, A}(data, Tuple(dims), crs, _string_attrs(attrs))
end

"""
    _string_attrs(attrs) -> Dict{String, Any}

Normalise an attribute dictionary to `Dict{String, Any}`, coercing symbol keys to
strings. Attributes arrive as `Dict{String, Any}` from backends, `Dict{Symbol, Any}`
from configs, and `Dict{Any, Any}` from stores.
"""
function _string_attrs(attrs)
    attrs isa Dict{String, Any} && return attrs
    return Dict{String, Any}(string(k) => v for (k, v) in pairs(attrs))
end

# Positional constructor for compatibility
function GeoArray(data::AbstractArray, dims::NTuple{N, Dimension}, crs::CoordinateSystem,
                  attrs = Dict{String, Any}()) where N
    T = eltype(data)
    GeoArray{T, N, typeof(data)}(data, dims, crs, _string_attrs(attrs))
end

Base.size(A::GeoArray) = size(A.data)
Base.ndims(A::GeoArray) = ndims(A.data)
Base.eltype(A::GeoArray) = eltype(A.data)
Base.length(A::GeoArray) = length(A.data)
Base.axes(A::GeoArray) = axes(A.data)
Base.getindex(A::GeoArray{T, N, A_}, i::Vararg{Int, N}) where {T, N, A_} = A.data[i...]
Base.setindex!(A::GeoArray{T, N, A_}, v, i::Vararg{Int, N}) where {T, N, A_} = A.data[i...] = v
Base.IndexStyle(::Type{<:GeoArray{T, N, A_}}) where {T, N, A_} = IndexStyle(A_)
Base.similar(A::GeoArray, ::Type{T}, dims::Dims) where {T} = similar(A.data, T, dims)

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
function GeoDataset(; variables::Dict = Dict{String, GeoArray}(),
                     coords::Dict = Dict{Symbol, GeoArray}(),
                     dims::Dict{Symbol, Dimension} = Dict{Symbol, Dimension}(),
                     crs::CoordinateSystem = CoordinateSystem(),
                     attrs::Dict{String, Any} = Dict{String, Any}(),
                     backend::Any = nothing,
                     source::String = "")
    vars_typed = Dict{String, GeoArray}(string(k) => v for (k, v) in variables)
    coords_typed = Dict{Symbol, GeoArray}(Symbol(k) => v for (k, v) in coords)
    GeoDataset(vars_typed, coords_typed, dims, crs, attrs, backend, source)
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

end # module GeoDataTypes