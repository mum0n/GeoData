"""
Shared test helpers and synthetic dataset factories for GeoData test suites.
"""

using Test
using Dates
using UUIDs
using Statistics
using GeoData
using Zarr, NCDatasets, DataFrames, GeoParquet, ArchGDAL, GeoInterface
using GeoData.GeoDataCoordinates: standardize_dimension_name, find_coord_indices,
    slice_indices, is_regular_grid, grid_spacing
using GeoData.GeoDataTypes: DIM_SPATIAL, DIM_TEMPORAL, DIM_PARAMETRIC, DIM_GENERIC,
    Dimension, GeoArray, GeoDataset
using GeoData: ZarrBackend, NCDatasetsBackend, GeoParquetBackend
using GeoData: GeoDataOperations

include("helpers.jl")

function make_dataset(; descending_lat = false)
    lon = collect(-65.0:0.5:-63.0)
    lat = descending_lat ? collect(45.0:-0.5:43.0) : collect(43.0:0.5:45.0)
    depth = collect(0.0:10.0:20.0)
    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    dd = Dimension(name = :depth, size = length(depth), coords = depth, units = "m",
                   dim_type = DIM_PARAMETRIC)
    data = rand(Float32, length(lon), length(lat), length(depth))
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict("temp" => GeoArray(data, (dl, dla, dd), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}()),
             :depth => GeoArray(depth, (dd,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla, :depth => dd), crs, Dict{String, Any}(), nothing,
        "synthetic")
end

function make_2d_dataset(data::AbstractMatrix, lon, lat, varname)
    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict(varname => GeoArray(data, (dl, dla,), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla), crs, Dict{String, Any}(), nothing, "synthetic")
end

function make_2d_dataset_ordered(data::AbstractMatrix, lon, lat, varname, order)
    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    dims = order == :lon_lat ? (dl, dla) : (dla, dl)
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict(varname => GeoArray(data, dims, crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla), crs, Dict{String, Any}(), nothing, "synthetic")
end

function make_3d_dataset(data::AbstractArray{<:Any, 3}, lon, lat, dep, varname)
    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    dd = Dimension(name = :depth, size = length(dep), coords = dep, units = "m",
                   dim_type = DIM_PARAMETRIC)
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict(varname => GeoArray(data, (dl, dla, dd), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}()),
             :depth => GeoArray(dep, (dd,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla, :depth => dd), crs, Dict{String, Any}(), nothing,
        "synthetic")
end
