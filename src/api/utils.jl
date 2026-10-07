"""
Shared utilities for GeoData API.
"""

export _resolve_backend

"""
    _resolve_backend(backend, uri::String) -> GeoBackend

Resolve a backend specification to a backend instance.
"""
function _resolve_backend(backend, uri::String)
    if backend === nothing
        return get_backend(infer_backend(uri))
    elseif backend isa Symbol
        return get_backend(backend)
    else
        return backend
    end
end