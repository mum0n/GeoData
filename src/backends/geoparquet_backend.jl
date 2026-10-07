"""
GeoParquet backend for vector/tabular geospatial data.
"""

using GeoParquet
using DataFrames
import GeoInterface as GI
import GeoFormatTypes as GFT
using WellKnownGeometry

"""
    GeoParquetBackend

Backend for GeoParquet files (vector/tabular data with geometry column).
"""
struct GeoParquetBackend <: GeoBackend end

function backend_capabilities(backend::GeoParquetBackend)
    return BackendCapabilities(read=true, write=true, lazy=false, chunked=false, compression=true, remote=false)
end

function backend_open(backend::GeoParquetBackend, uri::String; mode::String="r", kwargs...)
    mode == "r" || error("GeoParquet backend only supports read mode")
    path = _geoparquet_uri_to_path(uri)
    
    df = GeoParquet.read(path; kwargs...)
    
    # Identify geometry column
    geom_cols = GI.geometrycolumns(df)
    isempty(geom_cols) && error("No geometry column found in GeoParquet file")
    geom_col = first(geom_cols)
    
    # Extract coordinate arrays from geometries
    geoms = df[!, geom_col]
    points = _extract_points(geoms)
    
    vars = Dict{String, GeoArray}()
    coords = Dict{Symbol, GeoArray}()
    dims = Dict{Symbol, Dimension}()
    
    # Non-geometry columns become variables
    for col in names(df)
        col == geom_col && continue
        data = collect(df[!, col])
        dim = Dimension(name=:points, size=length(data), coords=collect(1:length(data)))
        vars[string(col)] = GeoArray(data, (dim,), CoordinateSystem(), Dict{String, Any}())
    end
    
    # Geometry as coordinates
    if !isempty(points)
        lons = [p[1] for p in points]
        lats = [p[2] for p in points]
        
        for (dim, vals, units) in ((:lon, lons, "degrees_east"), (:lat, lats, "degrees_north"))
            dim_obj = Dimension(name=dim, size=length(vals), coords=vals, units=units)
            coords[dim] = GeoArray(vals, (dim_obj,), CoordinateSystem(), Dict{String, Any}())
            dims[dim] = dim_obj
        end
    end
    
    return GeoDataset(vars, coords, dims, CoordinateSystem(), Dict{String, Any}(), backend, uri)
end

function backend_create(backend::GeoParquetBackend, uri::String; dims::Dict{Symbol, Dimension},
                       variables::Dict{String, <:AbstractArray}, coords::Dict{Symbol, <:AbstractArray},
                       crs::CoordinateSystem, attrs::Dict{String, Any}, kwargs...)
    path = _geoparquet_uri_to_path(uri)
    
    # Build DataFrame with geometry column
    n_points = length(first(values(variables)))
    
    # Get coordinates
    lons = haskey(coords, :lon) ? collect(coords[:lon]) : collect(1:n_points)
    lats = haskey(coords, :lat) ? collect(coords[:lat]) : zeros(n_points)
    
    geom = [_make_point(lons[i], lats[i]) for i in 1:n_points]
    
    df = DataFrame()
    df[!, :geometry] = geom
    for (name, data) in variables
        df[!, name] = data
    end
    
    GeoParquet.write(path, df, (:geometry,); kwargs...)
    
    return GeoDataset(
        Dict(k => GeoArray(v, _geoparquet_dims_to_tuple(k, v, dims, coords), crs, Dict{String, Any}()) for (k, v) in variables),
        Dict(k => GeoArray(v, (dims[k],), crs, Dict{String, Any}()) for (k, v) in coords),
        dims, crs, attrs, backend, uri
    )
end

function backend_write(backend::GeoParquetBackend, dataset::GeoDataset; variables::Dict{String, <:AbstractArray}=Dict(),
                      coords::Dict{Symbol, <:AbstractArray}=Dict(), attrs::Dict{String, Any}=Dict(),
                      mode::String="update", kwargs...)
    path = _geoparquet_uri_to_path(dataset.source)
    df = GeoParquet.read(path)
    
    for (name, data) in variables
        df[!, name] = data
    end
    for (dim, data) in coords
        df[!, string(dim)] = data
    end
    
    GeoParquet.write(path, df, (:geometry,); overwrite=true, kwargs...)
    return nothing
end

function backend_close(backend::GeoParquetBackend, dataset::GeoDataset)
    return nothing
end

function _geoparquet_uri_to_path(uri::String)
    startswith(uri, "geoparquet://") ? uri[13:end] : uri
end

function _extract_points(geoms)
    points = []
    for g in geoms
        if GI.geomtrait(g) isa GI.PointTrait
            push!(points, (GI.x(g), GI.y(g)))
        elseif GI.geomtrait(g) isa GI.MultiPointTrait
            for pt in GI.getgeom(g)
                push!(points, (GI.x(pt), GI.y(pt)))
            end
        end
    end
    return points
end

function _make_point(lon, lat)
    return WellKnownGeometry.Point(lon, lat)
end

function _geoparquet_infer_dim_names_from_size(shape, dims)
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

function _geoparquet_dims_to_tuple(name, data, dims, coords)
    dim_names = _geoparquet_infer_dim_names_from_size(size(data), dims)
    return tuple([dims[dn] for dn in dim_names]...)
end