"""
High-level API: querying values at points.
"""

export geovalues, geovalue, geoprofile, geotimeseries

"""
    geovalues(ds::GeoDataset, varnames::Vector{String}; kwargs...) -> Dict{String, Any}

Query values at specific coordinate points.

# Arguments
- `ds`: Source dataset
- `varnames`: Variables to query
- `kwargs`: Coordinate points as `lon = -65.0, lat = 45.0, depth = 50.0, time = ...`

# Returns
Dict mapping variable names to queried values (scalars or arrays).

# Examples
```julia
vals = geovalues(ds, ["temperature", "salinity"], lon=-65.0, lat=45.0, depth=50.0)
vals = geovalues(ds, ["temperature"], lon=[-65.0, -66.0], lat=[45.0, 46.0])
```
"""
function geovalues(ds::GeoDataset, varnames::Vector{String}; kwargs...)
    if ds.backend !== nothing && hasmethod(backend_values, Tuple{typeof(ds.backend), GeoDataset, Vector{String}})
        return backend_values(ds.backend, ds, varnames; kwargs...)
    else
        return values_at(ds, varnames; kwargs...)
    end
end

"""
    geovalue(ds::GeoDataset, varname::String; kwargs...) -> Any

Query a single value for a single variable.
"""
function geovalue(ds::GeoDataset, varname::String; kwargs...)
    vals = geovalues(ds, [varname]; kwargs...)
    return vals[varname]
end

# Convenience: vertical profile
"""
    geoprofile(ds::GeoDataset, varnames::Vector{String}; lon, lat, time=0.0) -> Dict

Get vertical profile at a location.
"""
function geoprofile(ds::GeoDataset, varnames::Vector{String}; lon, lat, time=0.0)
    depths = ds.coords[:depth].data
    vals = geovalues(ds, varnames; lon=lon, lat=lat, depth=depths, time=time)
    return merge(vals, Dict(:depth => depths))
end

# Convenience: time series
"""
    geotimeseries(ds::GeoDataset, varnames::Vector{String}; lon, lat, depth=0.0) -> Dict

Get time series at a location.
"""
function geotimeseries(ds::GeoDataset, varnames::Vector{String}; lon, lat, depth=0.0)
    times = ds.coords[:time].data
    vals = geovalues(ds, varnames; lon=lon, lat=lat, depth=depth, time=times)
    return merge(vals, Dict(:time => times))
end