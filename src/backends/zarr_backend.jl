"""
    zarr_backend.jl

Zarr v2 storage backend.

Zarr v2 stores arrays with no dimension metadata at all, so axes are inferred:
1-D arrays whose canonical name is a standard dimension name (`lon`, `lat`, `depth`,
`time`, `x`, `y`, `z`) become coordinate variables; every other array is a data variable.
That heuristic is wrong for an unlabelled 1-D time series, so a self-describing backend
(the NetCDF one) is preferable when the file format can be chosen.

`backend_open` returns the chunked Zarr array itself, so opening a store reads nothing
from disk; indexing then materialises only the elements it touches. Coordinates are read
in full because every selection is resolved by searching them.
"""

using Zarr

struct ZarrBackend <: GeoBackend end

"""
    _zarr_dims_to_tuple(var_shape, dims, coords) -> Tuple{Vararg{Dimension}}

Dimension records for a variable of shape `var_shape`, matched against the declared
`dims` (or a 1-D coordinate array with the same name).
"""
function _zarr_dims_to_tuple(var_shape, dims::Dict{Symbol, Dimension},
                             coords::Dict{Symbol, GeoArray})
    names = infer_dim_names_from_size(var_shape, dims)
    refs = Dimension[]
    for i in 1:length(var_shape)
        name = names[i]
        coord_ga = get(coords, name, nothing)
        if coord_ga === nothing && !haskey(dims, name)
            # An axis no declared dimension matches: name it, values unknown.
            push!(refs, Dimension(name = name, size = var_shape[i], coords = nothing))
            continue
        end
        source = coord_ga === nothing ? dims[name] : coord_ga.dims[1]
        push!(refs, Dimension(name = name, size = var_shape[i],
                              coords = source.coords,
                              units = source.units,
                              standard_name = source.standard_name,
                              calendar = source.calendar))
    end
    return tuple(refs...)
end

function backend_open(backend::ZarrBackend, uri::AbstractString; kwargs...)
    path = strip_scheme(uri)
    isdir(path) || error("Cannot open Zarr store '$(uri)': directory does not exist.")
    mode = haskey(kwargs, :mode) ? kwargs[:mode] : "r"
    g = Zarr.zopen(path, mode)

    array_names = String[]
    for k in keys(g.arrays)
        obj = g.arrays[k]
        obj isa Zarr.ZArray || continue
        push!(array_names, k)
    end

    coord_names, data_names = String[], String[]
    for k in array_names
        canon = standardize_dimension_name(Symbol(k))
        if ndims(g.arrays[k]) == 1 && canon in (:lon, :lat, :depth, :time, :x, :y, :z)
            push!(coord_names, k)
        else
            push!(data_names, k)
        end
    end

    crs = CoordinateSystem(crs = something(get(g.attrs, "crs", nothing), "EPSG:4326"))
    coords = Dict{Symbol, GeoArray}()
    dims = Dict{Symbol, Dimension}()
    for k in coord_names
        obj = g.arrays[k]
        canon = standardize_dimension_name(Symbol(k))
        c = collect(Float64, vec(obj[:]))
        d = Dimension(name = canon, size = length(c), coords = c,
                      units = get(obj.attrs, "units", nothing),
                      standard_name = get(obj.attrs, "standard_name", nothing),
                      calendar = get(obj.attrs, "calendar", nothing))
        dims[canon] = d
        coords[canon] = GeoArray(c, (d,), crs, Dict(obj.attrs))
    end

    vars = Dict{String, GeoArray}()
    for k in data_names
        obj = g.arrays[k]
        vars[k] = GeoArray(obj, _zarr_dims_to_tuple(size(obj), dims, coords), crs,
                           Dict(obj.attrs))
    end

    return GeoDataset(vars, coords, dims, crs, Dict(g.attrs), ZarrBackend(), path)
end

function backend_create(backend::ZarrBackend, uri::AbstractString,
                        dims::Dict{Symbol, Dimension}; kwargs...)
    path = strip_scheme(uri)
    if isdir(path)
        get(kwargs, :overwrite, false) || error(
            "Cannot create Zarr store '$(uri)': it already exists. Pass overwrite=true " *
            "to replace it.")
        rm(path, recursive = true)
    end
    return zgroup(path; attrs = get(kwargs, :attrs, Dict{String, Any}()))
end
function backend_write(backend::ZarrBackend, uri::AbstractString,
                       ds::GeoDataset; kwargs...)
    path = strip_scheme(uri)
    if isdir(path)
        # Zarr.jl reopens existing arrays read-only and refuses to create a new array in a
        # reopened group, so an existing store cannot be updated in place. Say so instead
        # of failing half-way through a write.
        get(kwargs, :overwrite, false) || error(
            "Zarr store '$(uri)' already exists. An existing store cannot be updated in " *
            "place by this backend: pass overwrite=true to rewrite it, or choose " *
            "another path.")
        rm(path; recursive = true)
    end
    g = backend_create(backend, path, ds.dims; attrs = ds.attrs, overwrite = true)

    # Dimensions first: coordinate arrays define the axes the variables are chunked on. A
    # dimension with no coordinate array is not an error: a Zarr array's `.zarray` carries its
    # own shape, so the axis survives a round trip without values, and the read path already
    # reconstructs such an axis (`_zarr_dims_to_tuple` names it and leaves values unknown).
    # NetCDF's CF bounds convention creates exactly this - `climatology_bounds` over
    # `nbounds` - and refusing the whole save over it left WOA23 writing a store that held
    # nothing but `lat`.
    for (dim, d) in ds.dims
        d.coords === nothing && continue
        zcreate(Float64, g, string(dim), length(d.coords);
                chunks = default_chunks((length(d.coords),)),
                attrs = Dict{String, Any}("units" => something(d.units, "")))
        g.arrays[string(dim)][:] = collect(Float64, d.coords)
    end

    for (name, ga) in ds.variables
        sz = size(ga.data)
        names = Symbol[]
        for i in 1:length(sz)
            nm = dimension_name(ga, i)
            nm === nothing && error("Axis $i of variable $name carries no declared " *
                                     "dimension name. Declare every dimension a variable " *
                                     "spans before saving.")
            push!(names, nm)
        end
        check_declared_dims(sz, names, ds.dims)
        zcreate(eltype(ga.data), g, name, sz...;
                chunks = default_chunks(sz), attrs = ga.attrs)
        g.arrays[name][:] = ga.data
    end

    return ds
end

"""
    backend_capabilities(backend::ZarrBackend)

Zarr stores are chunked and compressible, and `backend_open` hands out the chunked
arrays themselves, so reads stay lazy and slicing a loaded dataset returns views.
"""
function backend_capabilities(backend::ZarrBackend)
    BackendCapabilities(read = true, write = true, lazy = true, chunked = true)
end

backend_close(backend::ZarrBackend, ds::GeoDataset; kwargs...) = nothing
