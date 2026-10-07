"""
Zarr backend using Zarr.jl.
"""

using Zarr

"""
    ZarrBackend

Backend for Zarr stores (local or remote).
"""
struct ZarrBackend <: GeoBackend end

function backend_capabilities(backend::ZarrBackend)
    return BackendCapabilities(read=true, write=true, lazy=true, chunked=true, compression=true, remote=true)
end

function backend_open(backend::ZarrBackend, uri::String; mode::String="r", kwargs...)
    path = _zarr_uri_to_path(uri)
    if mode == "r"
        g = zopen(path; kwargs...)
    else
        g = zopen(path; mode=mode, kwargs...)
    end
    
    vars = Dict{String, GeoArray}()
    coords = Dict{Symbol, GeoArray}()
    dims = Dict{Symbol, Dimension}()
    
    # First pass: identify coordinate variables (1D arrays)
    coord_names = Set{Symbol}()
    for name in keys(g.arrays)
        obj = g.arrays[name]
        obj isa ZArray || continue
        canon = standardize_dimension_name(Symbol(name))
        if length(size(obj)) == 1
            data = reshape(obj[:], size(obj)...)
            dim = Dimension(name=canon, size=length(data), coords=vec(data),
                          units=get(obj.attrs, "units", ""),
                          standard_name=string(canon))
            coords[canon] = GeoArray(data, (dim,), CoordinateSystem(), Dict(obj.attrs))
            dims[canon] = dim
            push!(coord_names, canon)
        end
    end
    
    # Second pass: data variables
    for name in keys(g.arrays)
        obj = g.arrays[name]
        obj isa ZArray || continue
        canon = standardize_dimension_name(Symbol(name))
        canon in coord_names && continue
        
        data = reshape(obj[:], size(obj)...)
        dim_names = _infer_dim_names(obj, size(obj), dims, coord_names)
        dim_objs = [_nc_dim(dn, dims, coords) for dn in dim_names]
        vars[string(name)] = GeoArray(data, tuple(dim_objs...), CoordinateSystem(), Dict(obj.attrs))
    end
    
    return GeoDataset(vars, coords, dims, CoordinateSystem(), Dict(g.attrs), backend, uri)
end

function backend_create(backend::ZarrBackend, uri::String; dims::Dict{Symbol, Dimension},
                       variables::Dict{String, <:AbstractArray}, coords::Dict{Symbol, <:AbstractArray},
                       crs::CoordinateSystem, attrs::Dict{String, Any}, kwargs...)
    path = _zarr_uri_to_path(uri)
    g = zgroup(path; attrs=attrs, kwargs...)
    
    # Write coordinate arrays first
    for (dim, coord) in coords
        d = dims[dim]
        arr = zcreate(eltype(coord), g, string(dim), length(coord);
                      chunks=(min(1000, length(coord)),), attrs=Dict("units" => d.units))
        arr[:] = coord
    end
    
    # Write data variables
    for (name, data) in variables
dim_names = _zarr_infer_dim_names_from_size(size(data), dims)
        chunk_sizes = _default_chunks(size(data))
        arr = zcreate(eltype(data), g, name, size(data)...;
                      chunks=chunk_sizes, attrs=Dict())
        arr[:] = data
    end
    
    return GeoDataset(
        Dict(k => GeoArray(v, _zarr_dims_to_tuple(k, v, dims, coords), crs, Dict()) for (k, v) in variables),
        Dict(k => GeoArray(v, (dims[k],), crs, Dict()) for (k, v) in coords),
        dims, crs, attrs, backend, uri
    )
end

function backend_write(backend::ZarrBackend, dataset::GeoDataset; variables::Dict{String, <:AbstractArray}=Dict(),
                      coords::Dict{Symbol, <:AbstractArray}=Dict(), attrs::Dict{String, Any}=Dict(),
                      mode::String="update", kwargs...)
    path = _zarr_uri_to_path(dataset.source)
    g = zopen(path; mode="r+")
    
    for (name, data) in variables
        if haskey(g, name)
            g[name][:] = data
        else
            dim_names = _zarr_infer_dim_names_from_size(size(data), dataset.dims)
            chunk_sizes = _default_chunks(size(data))
            arr = zcreate(eltype(data), g, name, size(data)...; chunks=chunk_sizes)
            arr[:] = data
        end
    end
    
    for (dim, data) in coords
        if haskey(g, string(dim))
            g[string(dim)][:] = data
        end
    end
    
    for (k, v) in attrs
        g.attrs[k] = v
    end
    
    return nothing
end

function backend_close(backend::ZarrBackend, dataset::GeoDataset)
    return nothing
end

function _zarr_uri_to_path(uri::String)
    startswith(uri, "zarr://") ? uri[8:end] : uri
end

function _infer_dim_names(obj::ZArray, shape, dims, coord_names)
    if haskey(obj.attrs, "_ARRAY_DIMENSIONS")
        return Symbol.(split(obj.attrs["_ARRAY_DIMENSIONS"], ","))
    end
    return _zarr_infer_dim_names_from_size(shape, dims)
end

function _zarr_infer_dim_names_from_size(shape, dims)
    names = Symbol[]
    used = Set{Symbol}()
    # Standard order for 3D: lon, lat, depth (or x, y, z)
    standard_order = [:lon, :lat, :depth, :time, :x, :y, :z]
    for s in shape
        found = nothing
        # First try standard order
        for dim in standard_order
            if haskey(dims, dim) && dims[dim].size == s && !(dim in used)
                found = dim
                break
            end
        end
        # Fallback: any unused dimension with matching size
        if isnothing(found)
            for (dim, d) in dims
                if d.size == s && !(dim in used)
                    found = dim
                    break
                end
            end
        end
        push!(names, isnothing(found) ? Symbol("dim$(length(names)+1)") : found)
        found !== nothing && push!(used, found)
    end
    return names
end

function _nc_dim(dim_name, dims, coords)
    d = dims[dim_name]
    coords_arr = get(coords, dim_name, nothing)
    coords_vec = coords_arr isa GeoArray ? coords_arr.data : coords_arr
    return Dimension(name=dim_name, size=d.size, coords=coords_vec, units=d.units,
                     standard_name=d.standard_name, calendar=d.calendar)
end

function _default_chunks(shape)
    n = length(shape)
    if n == 0 return ()
    elseif n == 1 return (min(1000, shape[1]),)
    elseif n == 2 return (min(100, shape[1]), min(100, shape[2]))
    else return (min(50, shape[1]), min(50, shape[2]), min(10, shape[3]), ntuple(_ -> 1, n-3)...)
    end
end

function _zarr_dims_to_tuple(name, data, dims, coords)
    dim_names = _zarr_infer_dim_names_from_size(size(data), dims)
    return tuple([_nc_dim(dn, dims, coords) for dn in dim_names]...)
end