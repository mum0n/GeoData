"""
High-level API: selecting specific coordinate values.
"""

export geoselect, geonearest

"""
    geoselect(ds::GeoDataset; kwargs...) -> GeoDataset

Select specific coordinate values (exact match or nearest).

# Keywords (per dimension)
- `lon = -65.0` — single longitude (nearest)
- `lon = [-65.0, -66.0]` — multiple longitudes
- `lat = 45.0` — single latitude
- `depth = 50.0` — single depth
- `time = 0.0` — single time

# Examples
```julia
ds_sub = geoselect(ds, lon=-65.0, lat=45.0)
ds_sub = geoselect(ds, lon=[-65.0, -66.0], lat=[45.0, 46.0])
```
"""
function geoselect(ds::GeoDataset; kwargs...)
    if ds.backend !== nothing && hasmethod(backend_select, Tuple{typeof(ds.backend), GeoDataset})
        return backend_select(ds.backend, ds; kwargs...)
    else
        return select_values(ds; kwargs...)
    end
end

# Convenience: select nearest point
"""
    geonearest(ds::GeoDataset; lon, lat, depth=0.0, time=0.0) -> GeoDataset

Select the nearest point to the given coordinates.
"""
function geonearest(ds::GeoDataset; lon, lat, depth=0.0, time=0.0)
    return geoselect(ds; lon=lon, lat=lat, depth=depth, time=time)
end