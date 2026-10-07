"""
Backend registry for GeoData.
"""
module GeoDataRegistry

using ..GeoDataCoreTypes

const BACKEND_REGISTRY = Dict{Symbol, Any}()

"""
    register_backend(name::Symbol, backend)

Register a backend under a name.
"""
function register_backend(name::Symbol, backend)
    BACKEND_REGISTRY[name] = backend
    return backend
end

"""
    get_backend(name::Symbol) -> GeoBackend

Get a registered backend by name.
"""
function get_backend(name::Symbol)
    haskey(BACKEND_REGISTRY, name) || error("Unknown backend '$name'. Available: $(sort(collect(keys(BACKEND_REGISTRY))))")
    return BACKEND_REGISTRY[name]
end

"""
    list_backends() -> Vector{Symbol}

List all registered backend names.
"""
function list_backends()
    return sort(collect(keys(BACKEND_REGISTRY)))
end

"""
    infer_backend(uri::String) -> Symbol

Infer backend from URI scheme or file extension.
"""
function infer_backend(uri::String)
    u = lowercase(uri)
    if endswith(u, ".zarr") || startswith(u, "zarr://") || occursin(".zarr/", u)
        return :zarr
    elseif endswith(u, ".nc") || startswith(u, "nczarr://") || startswith(u, "netcdf://")
        return :ncdatasets
    elseif endswith(u, ".parquet") || startswith(u, "geoparquet://")
        return :geoparquet
    elseif endswith(u, ".yax") || occursin(".yax/", u) || startswith(u, "yaxarray://")
        return :yaxarray
    else
        return :ncdatasets
    end
end

export register_backend, get_backend, list_backends, infer_backend

end # module GeoDataRegistry