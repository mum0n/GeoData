"""
High-level API: saving datasets.
"""

export geosave

"""
    geosave(uri::String, data::GeoDataset; backend::Union{Symbol, GeoBackend, Nothing}=nothing, kwargs...) -> GeoDataset

Save a GeoDataset to a URI, auto-detecting the backend from the file extension or scheme.

# Arguments
- `uri`: Destination URI
- `data`: GeoDataset to save
- `backend`: Backend name or instance
- `kwargs`: Backend-specific options (chunks, compression, etc.)

# Examples
```julia
geosave("output/temperature.zarr", ds)
geosave("output/temperature.nc", ds; backend=:ncdatasets)
```
"""
function geosave(uri::String, data::GeoDataset; backend::Union{Symbol, GeoBackend, Nothing}=nothing, kwargs...)
    be = _resolve_backend(backend, uri)
    return _save_with_backend(be, uri, data; kwargs...)
end

function _save_with_backend(be::GeoBackend, uri::String, data::GeoDataset; kwargs...)
    # Use backend-specific URI path functions
    path = _get_backend_path(be, uri)
    
    # Prepare data for writing
    vars = Dict{String, AbstractArray}()
    for (name, ga) in data.variables
        vars[name] = ga.data
    end
    
    coords = Dict{Symbol, AbstractArray}()
    for (dim, ga) in data.coords
        coords[dim] = ga.data
    end
    
    dims = data.dims
    crs = data.crs
    attrs = data.attrs
    
    if isfile(path) || isdir(path)
        # Update existing
        backend_write(be, data; uri=path, variables=vars, coords=coords, attrs=attrs, mode="update", kwargs...)
    else
        # Create new
        backend_create(be, path; dims=dims, variables=vars, coords=coords, crs=crs, attrs=attrs, kwargs...)
    end
    
    return data
end

# Backend-specific path functions (call backend methods)
function _get_backend_path(be::NCDatasetsBackend, uri::String)
    return _ncdatasets_uri_to_path(uri)
end

function _get_backend_path(be::ZarrBackend, uri::String)
    return _zarr_uri_to_path(uri)
end

function _get_backend_path(be::YAXArraysBackend, uri::String)
    return _yax_uri_to_path(uri)
end

function _get_backend_path(be::GeoParquetBackend, uri::String)
    return _geoparquet_uri_to_path(uri)
end

# Fallback
_get_backend_path(be::GeoBackend, uri::String) = uri