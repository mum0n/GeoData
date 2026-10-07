"""
High-level API: slicing datasets.
"""

export geoslice, geoslice_bbox

"""
    geoslice(ds::GeoDataset; kwargs...) -> GeoDataset

Slice a dataset by coordinate ranges.

# Keywords (per dimension)
- `lon = (min, max)` — longitude range
- `lat = (min, max)` — latitude range
- `depth = (min, max)` — depth range
- `time = (min, max)` — time range

A keyword may also be a single value (nearest point, drops the dimension)
or a vector of values (select those points).

# Examples
```julia
ds_sub = geoslice(ds, lon=(-70, -55), lat=(40, 50))
ds_sub = geoslice(ds, depth=(0, 100), time=0.0)  # single time, drops time dim
```
"""
function geoslice(ds::GeoDataset; kwargs...)
    if ds.backend !== nothing && hasmethod(backend_slice, Tuple{typeof(ds.backend), GeoDataset})
        return backend_slice(ds.backend, ds; kwargs...)
    else
        return slice(ds; kwargs...)
    end
end

# Convenience: slice by bounding box
"""
    geoslice_bbox(ds::GeoDataset; lon::Tuple, lat::Tuple, depth::Tuple=(0, Inf), time::Tuple=(0, Inf)) -> GeoDataset

Slice by bounding box.
"""
function geoslice_bbox(ds::GeoDataset; lon::Tuple, lat::Tuple, depth::Tuple=(0, Inf), time::Tuple=(0, Inf))
    return geoslice(ds; lon=lon, lat=lat, depth=depth, time=time)
end