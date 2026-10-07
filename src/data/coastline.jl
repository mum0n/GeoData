"""
    coastline.jl

Natural Earth coastline fetch, load, and polygon operations via GeoData.
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using Downloads
using JSON3
using GeoInterface
using GeoParquet
using DataFrames

"""
    fetch_natural_earth_coastline(;
        lon_range::Tuple{Real, Real},
        lat_range::Tuple{Real, Real},
        resolution::AbstractString = "10m",
        output_path::AbstractString = joinpath("inputs", "coastline.parquet"),
        margin_deg::Real = 1.0,
        backend::Symbol = :geoparquet,
        verbose::Bool = true
    ) -> GeoDataset

Fetch Natural Earth coastline as GeoParquet via GeoData.

# Arguments
- `lon_range`: Longitude bounds (min, max)
- `lat_range`: Latitude bounds (min, max)
- `resolution`: Natural Earth resolution ("10m", "50m", "110m")
- `output_path`: Output file path (GeoParquet)
- `margin_deg`: Margin in degrees around bounds
- `backend`: Storage backend (:geoparquet, :zarr)
- `verbose`: Print progress messages

# Returns
GeoDataset with coastline LineString geometries.
"""
function fetch_natural_earth_coastline(;
    lon_range::Tuple{Real, Real},
    lat_range::Tuple{Real, Real},
    resolution::AbstractString = "10m",
    output_path::AbstractString = joinpath("inputs", "coastline.parquet"),
    margin_deg::Real = 1.0,
    backend::Symbol = :geoparquet,
    verbose::Bool = true
)
    res = lowercase(strip(resolution))
    res in ("10m", "50m", "110m") || error(
        "Natural Earth resolution \"$(resolution)\" must be 10m, 50m, or 110m")

    url = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/" *
          "ne_$(res)_coastline.geojson"

    mkpath(dirname(output_path))
    lo_lon, hi_lon = minmax(Float64(lon_range[1]), Float64(lon_range[2]))
    lo_lat, hi_lat = minmax(Float64(lat_range[1]), Float64(lat_range[2]))
    lo_lon -= Float64(margin_deg); hi_lon += Float64(margin_deg)
    lo_lat -= Float64(margin_deg); hi_lat += Float64(margin_deg)

    if verbose
        println("Retrieving Natural Earth $(res) coastline...")
        println("  bbox: lon [$(lo_lon), $(hi_lon)], lat [$(lo_lat), $(hi_lat)]")
    end

    # Download and parse GeoJSON
    text = Downloads.download(url)
    gj = JSON3.read(text)

    # Build DataFrame for GeoParquet
    rows = []
    feature_id = 0
    for feat in gj.features
        geom = feat.geometry
        if geom.type == "MultiLineString"
            for (ring_idx, ring) in enumerate(geom.coordinates)
                feature_id += 1
                lons = Float64[pt[1] for pt in ring]
                lats = Float64[pt[2] for pt in ring]
                # Check if any point is in bbox
                any(lo_lon <= lons .<= hi_lon) && any(lo_lat <= lats .<= hi_lat) || continue
                
                # Create LineString geometry
                coords = [(lons[i], lats[i]) for i in 1:length(lons)]
                line = GeoInterface.LineString(coords)
                
                props = hasproperty(feat, :properties) ? feat.properties : nothing
                nm = (isnothing(props) || !haskey(props, :name)) ? "coastline" : String(props[:name])
                code = Symbol(lowercase(replace(nm, r"[^A-Za-z0-9]+" => "_")))
                
                push!(rows, (id=feature_id, name=nm, code=string(code), geometry=line))
            end
        elseif geom.type == "LineString"
            feature_id += 1
            lons = Float64[pt[1] for pt in geom.coordinates]
            lats = Float64[pt[2] for pt in geom.coordinates]
            any(lo_lon <= lons .<= hi_lon) && any(lo_lat <= lats .<= hi_lat) || continue
            
            coords = [(lons[i], lats[i]) for i in 1:length(lons)]
            line = GeoInterface.LineString(coords)
            
            props = hasproperty(feat, :properties) ? feat.properties : nothing
            nm = (isnothing(props) || !haskey(props, :name)) ? "coastline" : String(props[:name])
            code = Symbol(lowercase(replace(nm, r"[^A-Za-z0-9]+" => "_")))
            
            push!(rows, (id=feature_id, name=nm, code=string(code), geometry=line))
        end
    end

    if isempty(rows)
        error("No coastline features found in bounding box. Check domain bounds.")
    end

    df = DataFrame(rows)
    GeoParquet.write(output_path, df, :geometry)
    
    # Also load as GeoDataset for immediate use
    ds = geoload(output_path; backend=backend)
    
    verbose && println("  kept $(length(rows)) coastline segments -> $(output_path)")
    return ds
end

"""
    load_coastline_geodata(filepath::AbstractString; backend=:geoparquet) -> GeoDataset

Load coastline from GeoParquet file.
"""
function load_coastline_geodata(filepath::AbstractString; backend::Symbol = :geoparquet)
    return geoload(filepath; backend=backend)
end

"""
    load_coastline_polygons_geodata(filepath::AbstractString) -> Vector{NamedTuple}

Load coastline polygons from GeoParquet file.
Returns vector of NamedTuples with (name, code, lons, lats).
"""
function load_coastline_polygons_geodata(filepath::AbstractString = "inputs/coastline.parquet")
    if endswith(filepath, ".parquet")
        ds = geoload(filepath; backend=:geoparquet)
        polys = NamedTuple[]
        for row in eachrow(ds.variables[:geometry].data)
            geom = row.geometry
            if geom isa GeoInterface.LineString
                coords = collect(geom)
                lons = Float64[pt[1] for pt in coords]
                lats = Float64[pt[2] for pt in coords]
                push!(polys, (
                    name = row.name,
                    code = Symbol(row.code),
                    lons = lons,
                    lats = lats
                ))
            end
        end
        return polys
    else
        # Fall back to legacy .dat format
        return load_coastline_polygons(filepath)
    end
end

"""
    is_point_on_land_geodata(lon, lat; coastline, coastline_file=nothing) -> Bool

Check if a point is on land using coastline polygons from GeoData.
"""
function is_point_on_land_geodata(
    lon::Real, lat::Real;
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    coastline_file::Union{Nothing, AbstractString} = nothing
)
    polys = if !isnothing(coastline)
        coastline
    elseif !isnothing(coastline_file)
        load_coastline_polygons_geodata(coastline_file)
    else
        load_coastline_polygons_geodata()
    end

    x = Float64(lon)
    y = Float64(lat)

    for poly in polys
        min_x, max_x = extrema(poly.lons)
        min_y, max_y = extrema(poly.lats)
        if x >= min_x && x <= max_x && y >= min_y && y <= max_y
            if point_in_polygon(x, y, poly.lons, poly.lats)
                return true
            end
        end
    end

    return false
end

"""
    is_marine_water_geodata(lon, lat; bathymetry, min_seabed_depth=0.0, coastline_file) -> Bool

Check if coordinates are in marine water using GeoData bathymetry and coastline.
"""
function is_marine_water_geodata(
    lon::Real,
    lat::Real;
    bathymetry::Union{Function, GeoDataset, AbstractString, Nothing} = nothing,
    min_seabed_depth::Real = 0.0,
    coastline_file::Union{Nothing, AbstractString} = nothing
)
    # 1. Coastline land exclusion
    if is_point_on_land_geodata(lon, lat; coastline_file = coastline_file)
        return false
    end

    # 2. Bathymetric depth constraint
    bathy_fn = if isnothing(bathymetry)
        def_path = "inputs/bathymetry.zarr"
        isfile(def_path) ? get_bathymetry_interpolator(geoload(def_path)) : nothing
    elseif bathymetry isa GeoDataset
        get_bathymetry_interpolator(bathymetry)
    elseif bathymetry isa AbstractString
        get_bathymetry_interpolator(geoload(bathymetry))
    else
        bathymetry
    end

    if isnothing(bathy_fn)
        return true
    end

    z_bed = bathy_fn(Float64(lon), Float64(lat))
    h_threshold = -max(0.0, Float64(min_seabed_depth))
    return z_bed <= h_threshold
end