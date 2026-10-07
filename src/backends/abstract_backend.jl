"""
Abstract backend type.
"""

"""
    GeoBackend

Abstract type for all backends. Concrete backends should subtype this.
"""
abstract type GeoBackend end

# Default backend capabilities (overridden by concrete backends)
function backend_capabilities(backend::GeoBackend)
    return BackendCapabilities()
end