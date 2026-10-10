"""
High-level API: saving datasets.
"""

export geosave

"""
    geosave(uri::AbstractString, data::GeoDataset; backend=nothing, overwrite=false, kwargs...) -> GeoDataset

Save a GeoDataset to a URI, auto-detecting the backend from the file extension or scheme
when `backend` is not given.

Writes into an existing store when one is present (variables and dimensions the
destination already has are updated in place); pass `overwrite=true` to replace it
instead.

# Arguments
- `uri`: Destination path or URI
- `data`: GeoDataset to save
- `backend`: Backend name or instance
- `overwrite`: Replace an existing destination instead of updating it
- `kwargs`: Backend-specific options (`format`, ...)

# Examples
```julia
geosave("output/temperature.zarr", ds)
geosave("output/temperature.nc", ds; backend=:ncdatasets)
geosave("output/temperature.zarr", ds; overwrite=true)
```
"""
function geosave(uri::AbstractString, data::GeoDataset;
                 backend = nothing, overwrite::Bool = false, kwargs...)
    be = resolve_backend(backend, uri)
    path = strip_scheme(uri)
    backend_write(be, path, data; overwrite = overwrite, kwargs...)
    return data
end

"""
    _save_with_backend(be, uri, data; kwargs...)

Escape hatch for callers holding a backend already; `geosave` resolves one for you.
"""
_save_with_backend(be::GeoBackend, uri::AbstractString, data::GeoDataset; kwargs...) =
    backend_write(be, strip_scheme(uri), data; kwargs...)
