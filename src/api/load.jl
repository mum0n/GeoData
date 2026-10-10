"""
High-level API: loading datasets.
"""

export geoload, open_dataset, close_dataset, infer_backend

"""
    geoload(uri::String; backend::Union{Symbol, GeoBackend, Nothing}=nothing, kwargs...) -> GeoDataset

Load a dataset from a URI, auto-detecting the backend from the file extension or scheme.

# Arguments
- `uri`: Source URI (file path, `zarr://`, `netcdf://`, `nczarr://`, `geoparquet://`)
- `backend`: Backend name (`:zarr`, `:ncdatasets`, `:nczarr`, `:geoparquet`) or instance
- `kwargs`: Backend-specific options

# Examples
```julia
ds = geoload("data/temperature.zarr")
ds = geoload("data/temperature.nc", backend=:ncdatasets)
ds = geoload("data/temperature.parquet")
```
"""
function geoload(uri::String; backend::Union{Symbol, GeoBackend, Nothing}=nothing, kwargs...)
    be = resolve_backend(backend, uri)
    return backend_open(be, uri; kwargs...)
end

"""
    open_dataset(uri::String; backend::Union{Symbol, GeoBackend, Nothing}=nothing, kwargs...) -> GeoDataset

Alias for `geoload`.
"""
open_dataset(uri::String; backend=nothing, kwargs...) = geoload(uri; backend=backend, kwargs...)

"""
    close_dataset(dataset::GeoDataset) -> Nothing

Close a dataset and release resources.
"""
function close_dataset(dataset::GeoDataset)
    if dataset.backend !== nothing
        backend_close(dataset.backend, dataset)
    end
    return nothing
end