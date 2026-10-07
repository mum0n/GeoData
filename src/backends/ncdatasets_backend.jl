"""
NCDatasets backend for NetCDF files (classic and Zarr).
"""

using NCDatasets
using CFTime

"""
    NCDatasetsBackend

Backend for NetCDF files using NCDatasets.jl. Supports classic NetCDF and NCZarr.
"""
struct NCDatasetsBackend <: GeoBackend
    format::String  # "netcdf4", "nczarr", etc.
end

NCDatasetsBackend(; format="netcdf4") = NCDatasetsBackend(format)

function backend_capabilities(backend::NCDatasetsBackend)
    return BackendCapabilities(read=true, write=true, lazy=false, chunked=true, compression=true, remote=false)
end

function backend_open(backend::NCDatasetsBackend, uri::String; mode::String="r", kwargs...)
    path = _ncdatasets_uri_to_path(uri)
    format = _infer_format(uri, backend.format)
    
    return NCDataset(path, mode; format=format, kwargs...) do ds
        vars = Dict{String, GeoArray}()
        coords = Dict{Symbol, GeoArray}()
        dims = Dict{Symbol, Dimension}()
        
        # First pass: identify coordinate variables (1D vars matching dimension names)
        coord_names = Set{Symbol}()
        for name in keys(ds)
            v = ds[name]
            if ndims(v) == 1
                canon = standardize_dimension_name(Symbol(name))
                if canon in (:lon, :lat, :depth, :time)
                    data = _read_coord(v)
                    dim = Dimension(name=canon, size=length(data), coords=data,
                                  units=get(v.attrib, "units", ""),
                                  standard_name=string(canon),
                                  calendar=get(v.attrib, "calendar", nothing))
                    coords[canon] = GeoArray(data, (dim,), CoordinateSystem(), Dict{String, Any}(v.attrib))
                    dims[canon] = dim
                    push!(coord_names, canon)
                end
            end
        end
        
        # Second pass: data variables
        for name in keys(ds)
            v = ds[name]
            canon = standardize_dimension_name(Symbol(name))
            canon in coord_names && continue
            
            data = reshape(Array{Float64}(v[:]), size(v)...)
            dn = String.(dimnames(v))
            dim_objs = [_nc_dim_from_var(dn[i], ds, dims, coords) for i in 1:length(dn)]
            vars[name] = GeoArray(data, tuple(dim_objs...), CoordinateSystem(), Dict{String, Any}(v.attrib))
        end
        
        GeoDataset(vars, coords, dims, CoordinateSystem(), Dict{String, Any}(ds.attrib), backend, uri)
    end
end

function backend_create(backend::NCDatasetsBackend, uri::String; dims::Dict{Symbol, Dimension},
                       variables::Dict{String, <:AbstractArray}, coords::Dict{Symbol, <:AbstractArray},
                       crs::CoordinateSystem, attrs::Dict{String, Any}, kwargs...)
    path = _ncdatasets_uri_to_path(uri)
    format = _infer_format(uri, backend.format)
    
    NCDataset(path, "c"; format=format, kwargs...) do ds
        # Define dimensions
        for (dim, d) in dims
            defDim(ds, string(dim), d.size)
        end
        
        # Write coordinate variables
        for (dim, coord) in coords
            v = defVar(ds, string(dim), eltype(coord), (string(dim),))
            v.attrib["units"] = dims[dim].units
            v[:] = coord
        end
        
        # Write data variables
        for (name, data) in variables
dim_names = _ncdatasets_infer_dim_names_from_size(size(data), dims)
            v = defVar(ds, name, eltype(data), tuple(String.(dim_names)...))
            v[:] = data
        end
        
        for (k, v) in attrs
            ds.attrib[k] = v
        end
    end
    
    return GeoDataset(
        Dict(k => GeoArray(v, _ncdatasets_dims_to_tuple(k, v, dims, coords), crs, Dict{String, Any}()) for (k, v) in variables),
        Dict(k => GeoArray(v, (dims[k],), crs, Dict{String, Any}()) for (k, v) in coords),
        dims, crs, attrs, backend, uri
    )
end

function backend_write(backend::NCDatasetsBackend, dataset::GeoDataset; variables::Dict{String, <:AbstractArray}=Dict(),
                      coords::Dict{Symbol, <:AbstractArray}=Dict(), attrs::Dict{String, Any}=Dict(),
                      mode::String="update", kwargs...)
    path = _ncdatasets_uri_to_path(dataset.source)
    NCDataset(path, "r+") do ds
        for (name, data) in variables
            if haskey(ds, name)
                ds[name][:] = data
            else
                dim_names = _ncdatasets_infer_dim_names_from_size(size(data), dataset.dims)
                v = defVar(ds, name, eltype(data), tuple(String.(dim_names)...))
                v[:] = data
            end
        end
        
        for (dim, data) in coords
            if haskey(ds, string(dim))
                ds[string(dim)][:] = data
            end
        end
        
        for (k, v) in attrs
            ds.attrib[k] = v
        end
    end
    
    return nothing
end

function backend_close(backend::NCDatasetsBackend, dataset::GeoDataset)
    return nothing
end

function _ncdatasets_uri_to_path(uri::String)
    for prefix in ("netcdf://", "nczarr://", "file://")
        startswith(uri, prefix) && return uri[length(prefix)+1:end]
    end
    return uri
end

function _infer_format(uri::String, default::String)
    lowercase(uri) |> u -> startswith(u, "nczarr://") ? :nczarr : Symbol(default)
end

function _read_coord(v)
    data = collect(Float64, v[:])
    if haskey(v.attrib, "units")
        units = v.attrib["units"]
        if occursin("since", units)  # time coordinate
            try
                epoch, calendar = parse_time_units(units)
                data = [datetime_to_time(DateTime(t), units) for t in CFTime.num2date(data, units)]
            catch
            end
        end
    end
    return data
end

function _nc_dim_from_var(dim_name, ds, dims, coords)
    if haskey(dims, Symbol(dim_name))
        d = dims[Symbol(dim_name)]
        coords_arr = get(coords, Symbol(dim_name), nothing)
        coords_vec = coords_arr isa GeoArray ? coords_arr.data : coords_arr
        return Dimension(name=Symbol(dim_name), size=d.size, coords=coords_vec, units=d.units,
                         standard_name=d.standard_name, calendar=d.calendar)
    end
    # Fallback
    v = ds[dim_name]
    return Dimension(name=Symbol(dim_name), size=length(v), coords=collect(Float64, v[:]),
                     units=get(v.attrib, "units", ""))
end

function _ncdatasets_infer_dim_names_from_size(shape, dims)
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

function _ncdatasets_dims_to_tuple(name, data, dims, coords)
    dim_names = _ncdatasets_infer_dim_names_from_size(size(data), dims)
    return tuple([_nc_dim_from_var(String(dn), nothing, dims, coords) for dn in dim_names]...)
end