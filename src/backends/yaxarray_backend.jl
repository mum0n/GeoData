"""
YAXArrays backend for lazy, labeled array operations.
"""

using YAXArrays
using DimensionalData

"""
    YAXArraysBackend

Backend using YAXArrays for lazy cube operations.
"""
struct YAXArraysBackend <: GeoBackend end

function backend_capabilities(backend::YAXArraysBackend)
    return BackendCapabilities(read=true, write=true, lazy=true, chunked=true, compression=true, parallel=true)
end

function backend_open(backend::YAXArraysBackend, uri::String; mode::String="r", kwargs...)
    mode == "r" || error("YAXArrays backend only supports read mode")
    path = _yax_uri_to_path(uri)
    
    ds = open_dataset(path; kwargs...)
    cube = Cube(ds; kwargs...)
    
    vars = Dict{String, GeoArray}()
    coords = Dict{Symbol, GeoArray}()
    dims = Dict{Symbol, Dimension}()
    
    # Extract coordinates from cube axes
    for (dim_name, axis) in axes(cube)
        canon = standardize_dimension_name(Symbol(dim_name))
        vals = collect(Float64, axis.val)
        dim = Dimension(name=canon, size=length(vals), coords=vals,
                      units=get(axis.properties, "units", ""),
                      standard_name=string(canon),
                      calendar=get(axis.properties, "calendar", nothing))
        coords[canon] = GeoArray(vals, (dim,), CoordinateSystem(), Dict(axis.properties))
        dims[canon] = dim
    end
    
    # Extract variables
    for name in keys(cube)
        yax = cube[name]
        data = yax.data  # lazy
        dim_names = [standardize_dimension_name(Symbol(d.name)) for d in dims(yax)]
        dim_objs = [dims[dn] for dn in dim_names]
        vars[string(name)] = GeoArray(data, tuple(dim_objs...), CoordinateSystem(), Dict(yax.properties))
    end
    
    return GeoDataset(vars, coords, dims, CoordinateSystem(), Dict(), backend, uri)
end

function backend_create(backend::YAXArraysBackend, uri::String; dims::Dict{Symbol, Dimension},
                       variables::Dict{String, <:AbstractArray}, coords::Dict{Symbol, <:AbstractArray},
                       crs::CoordinateSystem, attrs::Dict{String, Any}, kwargs...)
    path = _yax_uri_to_path(uri)
    
    # Build YAXArrays cube
    axes_list = [_yax_dim(dims[d], coords[d]) for d in keys(dims)]
    cubes = Dict{String, Any}()
    
    for (name, data) in variables
        axes_for_var = [axes_list[findfirst(d -> d.name == dn, axes_list)] for dn in _yax_infer_dim_names_from_size(size(data), dims)]
        cubes[name] = YAXArray(axes_for_var, data)
    end
    
    # Save using YAXArrays
    cube = Cube(cubes)
    savecube(cube, path; backend=:zarr, kwargs...)
    
    return GeoDataset(
        Dict(k => GeoArray(v, _yax_dims_to_tuple(k, v, dims, coords), crs, Dict()) for (k, v) in variables),
        Dict(k => GeoArray(v, (dims[k],), crs, Dict()) for (k, v) in coords),
        dims, crs, attrs, backend, uri
    )
end

function backend_write(backend::YAXArraysBackend, dataset::GeoDataset; variables::Dict{String, <:AbstractArray}=Dict(),
                      coords::Dict{Symbol, <:AbstractArray}=Dict(), attrs::Dict{String, Any}=Dict(),
                      mode::String="update", kwargs...)
    path = _yax_uri_to_path(dataset.source)
    cube = open_dataset(path)
    
    for (name, data) in variables
        cube[name][:] = data
    end
    for (dim, data) in coords
        cube[dim][:] = data
    end
    
    savecube(cube, path; backend=:zarr, overwrite=true, kwargs...)
    return nothing
end

function backend_close(backend::YAXArraysBackend, dataset::GeoDataset)
    return nothing
end

function _yax_uri_to_path(uri::String)
    startswith(uri, "yaxarray://") ? uri[11:end] : uri
end

function _yax_dim(dim::Dimension, coord::AbstractArray)
    Dim{:x}(coord)  # Simplified - in reality would map to proper dimension type
end

function _yax_infer_dim_names_from_size(shape, dims)
    names = Symbol[]
    for s in shape
        found = nothing
        for (dim, d) in dims
            if d.size == s
                found = dim
                break
            end
        end
        push!(names, isnothing(found) ? Symbol("dim$(length(names)+1)") : found)
    end
    return names
end

function _yax_dims_to_tuple(name, data, dims, coords)
    dim_names = _yax_infer_dim_names_from_size(size(data), dims)
    return tuple([dims[dn] for dn in dim_names]...)
end