"""
    GeoDataBackends

The backend interface: the four operations a backend must implement, and nothing else.

A backend is a *storage driver*. It opens, creates, writes, closes, and reports what it
can do. It does not plan queries: every coordinate operation (slice, select, values, join,
index, regrid, aggregate) is implemented once, generically, over
`GeoDataOperations` and `GeoDataCoordinates`, in `Core`. That is the whole design of this
file: keep the extension surface minimal enough that adding a backend is a day's work and
that no caller has code paths that can never be taken.

What was removed and why: `backend_slice`, `backend_select`, `backend_values`,
`backend_index`, `backend_join`, `backend_regrid`, `backend_aggregate`, `backend_subset`,
`backend_profile`, `backend_series`, and `backend_resample` were declared here and
implemented nowhere, while eight call sites in `API` guarded on them with
`backend_supports`. Every one of those guards was a branch that could never fire. Query
planning was never the bottleneck; it was complexity with no benefit. If a backend ever
does need its own query path, add exactly one operation, with a test proving the fast path
and the generic path agree before the guard ships.
"""

module GeoDataBackends

using ..GeoDataTypes:
    GeoDataset,
    GeoBackend,
    BackendCapabilities

export GeoBackendInterface,
    backend_open, backend_close, backend_create, backend_write, backend_capabilities

# ============================================================
# Lifecycle
# ============================================================

"""
    backend_open(backend::GeoBackend, uri::AbstractString; kwargs...) -> GeoDataset

Open a source and return it as a `GeoDataset`.

`mode` selects `"r"` (read-only, the default) or `"r+"` for a dataset that will be
written. A backend that cannot open `uri` errors naming the path and the reason; a missing
file is an error, not an empty dataset.
"""
function backend_open end

"""
    backend_create(backend::GeoBackend, uri::AbstractString, dims::Dict{Symbol, Dimension}; kwargs...)

Create an empty destination for `dims` and return a handle to write into. `overwrite =
true` replaces an existing destination. Backends that store filesystems-as-objects (Zarr)
overwrite by removing and recreating; backends that support in-place appends (NetCDF) do
not.
"""
function backend_create end

"""
    backend_write(backend::GeoBackend, uri::AbstractString, dataset::GeoDataset; kwargs...)

Write `dataset` to `uri`, creating the destination if needed and updating it in place if it
exists and the backend supports that. Backends that cannot update in place (Zarr) error
rather than failing part-way through a write.
"""
function backend_write end

"""
    backend_close(backend::GeoBackend, uri::AbstractString, dataset::GeoDataset; kwargs...)

Release any handle the backend holds for `dataset`. The default is a no-op: readers that
manage their own lifetime need nothing, and forcing every backend to implement this was
pure ceremony.
"""
function backend_close(backend::GeoBackend, uri::AbstractString, dataset::GeoDataset; kwargs...) end

"""
    backend_capabilities(backend::GeoBackend) -> BackendCapabilities

Return what this backend supports.

*Capabilities are reported, never used for dispatch*: a backend reporting `lazy = false`
must still produce correct results, it simply materialises on open. A caller branching on
capabilities is optimising, never deciding correctness.

The default is defined in `GeoData.jl`, alongside the four backend implementations it
covers, so it cannot be separated from the module that registers them.
"""
function backend_capabilities end

# ============================================================
# URI handling
# ============================================================

const SUPPORTED_URI_SCHEMES = ("zarr://", "nczarr://", "netcdf://", "geoparquet://", "file://")

"""
    strip_scheme(uri::AbstractString) -> String

Remove a recognised `scheme://` prefix from `uri`, returning the local path portion.

Recognised schemes: `zarr://`, `nczarr://`, `netcdf://`, `geoparquet://`, `file://`.
A URI without a recognised scheme is returned unchanged. Use this instead of slicing a
prefix by hand: every hand-written offset in this file's history was wrong by one.
"""
function strip_scheme(uri::AbstractString)
    s = string(uri)
    lower = lowercase(s)
    for prefix in SUPPORTED_URI_SCHEMES
        if startswith(lower, prefix)
            return s[length(prefix) + 1 : end]
        end
    end
    return s
end

export strip_scheme, SUPPORTED_URI_SCHEMES

end # module GeoDataBackends
