"""
Generic coordinate system utilities for GeoData.

Cartesian by default. Geographic specialization in GeoDataGeoCoordinates.
"""
module GeoDataCoordinates

using ..GeoDataCoreTypes
using LinearAlgebra
using Dates

# ============================================================
# Generic dimension handling (no CF assumptions)
# ============================================================

"""
    standardize_dimension_name(name::Symbol; mapping::Dict{Symbol,Symbol}=Dict()) -> Symbol

Map dimension name variants to canonical names using optional mapping.
Defaults to identity (no standardization).
"""
function standardize_dimension_name(name::Symbol; mapping::Dict{Symbol,Symbol}=Dict{Symbol,Symbol}())
    get(mapping, name, name)
end

"""
    infer_dimension(name::Symbol, array::AbstractArray; attrs::Dict = Dict(), mapping::Dict{Symbol,Symbol}=Dict()) -> Dimension

Infer a `Dimension` from a coordinate variable. Generic version.
"""
function infer_dimension(name::Symbol, array::AbstractArray; attrs::Dict = Dict(), mapping::Dict{Symbol,Symbol}=Dict())
    canon = standardize_dimension_name(name; mapping=mapping)
    coords = vec(Float64.(array))
    units = get(attrs, "units", "")
    standard_name = get(attrs, "standard_name", nothing)
    calendar = get(attrs, "calendar", nothing)
    dim_type = _infer_dim_type(canon, units)

    Dimension(
        name = canon,
        size = length(coords),
        coords = coords,
        units = units,
        standard_name = standard_name,
        dim_type = dim_type,
        calendar = calendar,
    )
end

function _infer_dim_type(canon::Symbol, units::String)
    canon in (:x, :y, :z, :lon, :lat, :depth, :height) && return DIM_SPATIAL
    canon in (:t, :time, :date, :datetime) && return DIM_TEMPORAL
    canon in (:chain, :draw, :sample, :param, :parameter) && return DIM_PARAMETRIC
    occursin("degree", lowercase(units)) && return DIM_SPATIAL
    occursin("since", lowercase(units)) && return DIM_TEMPORAL
    return DIM_GENERIC
end

"""
    wrap_longitude(lon::AbstractVector; convention::Symbol = :pm180) -> Vector{Float64}

Wrap longitude values (geographic specialization).
"""
function wrap_longitude(lon::AbstractVector; convention::Symbol = :pm180)
    if convention === :pm180
        return [mod(v + 180.0, 360.0) - 180.0 for v in Float64.(lon)]
    elseif convention === :0_360
        return [mod(v, 360.0) for v in Float64.(lon)]
    else
        error("Unknown longitude convention: $convention")
    end
end

"""
    sort_coordinates(coords::AbstractVector) -> (Vector{Float64}, Vector{Int})

Sort coordinates and return both sorted values and permutation indices.
"""
function sort_coordinates(coords::AbstractVector)
    vals = Float64.(coords)
    perm = sortperm(vals)
    return vals[perm], perm
end

"""
    normalize_depth(depth::AbstractVector) -> (Vector{Float64}, Bool)

Normalize depth to positive-down convention (geographic specialization).
"""
function normalize_depth(depth::AbstractVector)
    dep = Float64.(depth)
    dep_max = maximum(abs, dep)
    if minimum(dep) >= -0.5 * dep_max
        return dep, false
    else
        return -dep, true
    end
end

# ============================================================
# Generic coordinate matching and indexing
# ============================================================

"""
    find_coord_indices(coords::AbstractVector, target::Real; mode::Symbol = :nearest) -> Int

Find index of target value in coordinate array.

# Modes
- `:nearest` - nearest neighbor
- `:exact` - exact match (throws if not found)
- `:floor` - largest index with coord <= target
- `:ceil` - smallest index with coord >= target
"""
function find_coord_indices(coords::AbstractVector, target::Real; mode::Symbol = :nearest)
    c = Float64.(coords)
    target = Float64(target)

    if mode === :exact
        idx = findfirst(isequal(target), c)
        isnothing(idx) && error("Exact coordinate $target not found in $(c[1])..$(c[end])")
        return idx
    elseif mode === :nearest
        return argmin(abs.(c .- target))
    elseif mode === :floor
        idx = findlast(x -> x <= target, c)
        isnothing(idx) && error("No coordinate <= $target in $(c[1])..$(c[end])")
        return idx
    elseif mode === :ceil
        idx = findfirst(x -> x >= target, c)
        isnothing(idx) && error("No coordinate >= $target in $(c[1])..$(c[end])")
        return idx
    else
        error("Unknown mode: $mode")
    end
end

"""
    slice_indices(coords::AbstractVector, range::Tuple{Real, Real}; inclusive::Bool = true) -> UnitRange{Int}

Get index range for a coordinate slice.

# Arguments
- `coords`: Coordinate vector (must be sorted)
- `range`: (min, max) in coordinate units
- `inclusive`: Whether bounds are inclusive

# Returns
`UnitRange` of indices.
"""
function slice_indices(coords::AbstractVector, range::Tuple{Real, Real}; inclusive::Bool = true)
    c = Float64.(coords)
    lo, hi = Float64.(range)
    lo <= hi || error("Invalid range: $lo > $hi")

    i1 = find_coord_indices(c, lo; mode = inclusive ? :ceil : :floor)
    i2 = find_coord_indices(c, hi; mode = inclusive ? :floor : :ceil)

    i1 <= i2 || error("Empty slice: no coordinates in [$lo, $hi]")
    return i1:i2
end

"""
    slice_indices(coords::AbstractVector, values::AbstractVector) -> Vector{Int}

Get indices for exact coordinate values.
"""
function slice_indices(coords::AbstractVector, values::AbstractVector)
    c = Float64.(coords)
    return [find_coord_indices(c, Float64(v); mode = :nearest) for v in values]
end

# ============================================================
# Grid utilities
# ============================================================

"""
    is_regular_grid(coords::AbstractVector; rtol::Real = 1e-10) -> Bool

Check if coordinates form a regularly spaced grid.
"""
function is_regular_grid(coords::AbstractVector; rtol::Real = 1e-10)
    c = Float64.(coords)
    length(c) < 3 && return true
    d = diff(c)
    return all(isapprox.(d, d[1]; rtol = rtol))
end

"""
    grid_spacing(coords::AbstractVector) -> Float64

Get grid spacing (assumes regular grid).
"""
function grid_spacing(coords::AbstractVector)
    c = Float64.(coords)
    length(c) < 2 && return NaN
    return c[2] - c[1]
end

"""
    bounding_box(ds::GeoDataset) -> NamedTuple

Compute the bounding box of a dataset (all dimensions).
"""
function bounding_box(ds::GeoDataset)
    bb = NamedTuple()
    for (dim, d) in ds.dims
        if haskey(ds.coords, dim)
            c = ds.coords[dim].data
            bb = merge(bb, (dim => (Float64(minimum(c)), Float64(maximum(c))),))
        end
    end
    return bb
end

# ============================================================
# Time handling (generic, supports CF-like units)
# ============================================================

"""
    parse_time_units(units::String) -> (epoch::DateTime, calendar::String)

Parse CF time units string like "days since 1990-01-01".
"""
function parse_time_units(units::String)
    m = match(r"(\w+)\s+since\s+([\d\-T:Z]+)", units)
    isnothing(m) && error("Cannot parse time units: $units")
    unit, epoch_str = m.captures
    epoch = DateTime(epoch_str)
    calendar = "standard"
    return epoch, calendar
end

"""
    time_to_datetime(time_vals::AbstractVector, units::String) -> Vector{DateTime}

Convert numeric time values to DateTime.
"""
function time_to_datetime(time_vals::AbstractVector, units::String)
    epoch, _ = parse_time_units(units)
    unit = lowercase(first(split(units, ' ')))
    return [epoch + _time_delta(t, unit) for t in Float64.(time_vals)]
end

function _time_delta(val::Float64, unit::String)
    if unit in ("second", "seconds", "s")
        return Second(round(Int, val))
    elseif unit in ("minute", "minutes", "min")
        return Minute(round(Int, val))
    elseif unit in ("hour", "hours", "h")
        return Hour(round(Int, val))
    elseif unit in ("day", "days", "d")
        return Day(round(Int, val))
    elseif unit in ("month", "months")
        return Month(round(Int, val))
    elseif unit in ("year", "years", "y")
        return Year(round(Int, val))
    else
        error("Unknown time unit: $unit")
    end
end

"""
    datetime_to_time(dt::DateTime, units::String) -> Float64

Convert DateTime to numeric time value.
"""
function datetime_to_time(dt::DateTime, units::String)
    epoch, _ = parse_time_units(units)
    unit = lowercase(first(split(units, ' ')))
    delta = dt - epoch
    return _delta_to_unit(delta, unit)
end

function _delta_to_unit(delta::DateTime, unit::String)
    if unit in ("second", "seconds", "s")
        return Float64(Millisecond(delta)) / 1000
    elseif unit in ("minute", "minutes", "min")
        return Float64(Millisecond(delta)) / 60000
    elseif unit in ("hour", "hours", "h")
        return Float64(Millisecond(delta)) / 3600000
    elseif unit in ("day", "days", "d")
        return Float64(Millisecond(delta)) / 86400000
    else
        error("Unsupported time unit for conversion: $unit")
    end
end

# Default dimension mapping for geographic data (exported for opt-in use)
const GEO_DIMENSION_ALIASES = Dict{Symbol, Vector{Symbol}}(
    :lon => [:lon, :longitude, :x, :long],
    :lat => [:lat, :latitude, :y, :latit],
    :depth => [:depth, :z, :lev, :level, :depths, :height],
    :time => [:time, :t, :date, :datetime],
    :chain => [:chain, :chains],
    :draw => [:draw, :sample, :draws, :samples],
)

const GEO_UNITS = Dict{Symbol, String}(
    :lon => "degrees_east",
    :lat => "degrees_north",
    :depth => "meters",
    :time => "seconds since 1970-01-01 00:00:00",
)

export standardize_dimension_name, infer_dimension, wrap_longitude,
    sort_coordinates, normalize_depth, find_coord_indices, slice_indices,
    is_regular_grid, grid_spacing, bounding_box,
    parse_time_units, time_to_datetime, datetime_to_time,
    GEO_DIMENSION_ALIASES, GEO_UNITS

end # module GeoDataCoordinates