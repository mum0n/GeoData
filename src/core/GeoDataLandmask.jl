"""
    landmask.jl

Point-in-polygon and land/ocean tests over geometry the caller already holds.

This is generic geometry, not data acquisition: it decides whether a coordinate falls
inside a ring. It lives in `GeoData` because a workflow that has *already* downloaded a
coastline — or has coastline rings from any source at all — should not have to import a
package whose job is downloading. The boundary that matters is not "coastlines" but
"network access": `GeoData` contains no network code, and the fetcher that produces a
coastline file stays in `GeoDataSources`.
"""

using GeoData: GeoDataset

export point_in_polygon, is_point_on_land_geodata, is_marine_water_geodata,
    load_coastline_polygons

"""
    point_in_polygon(x::Real, y::Real, poly_lons::AbstractVector, poly_lats::AbstractVector) -> Bool

Determine whether 2D coordinates (x, y) fall within a polygon defined by vertex vectors
(poly_lons, poly_lats) using the standard Jordan curve theorem (ray-casting algorithm).
"""
function point_in_polygon(x::Real, y::Real, poly_lons::AbstractVector, poly_lats::AbstractVector)::Bool
    n = min(length(poly_lons), length(poly_lats))
    n < 3 && return false
    inside = false
    j = n
    px = Float64(x)
    py = Float64(y)
    @inbounds for i in 1:n
        xi = Float64(poly_lons[i])
        yi = Float64(poly_lats[i])
        xj = Float64(poly_lons[j])
        yj = Float64(poly_lats[j])

        if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

"""
    load_coastline_polygons(filepath::AbstractString) -> Vector{NamedTuple}

Load coastline rings from a GeoParquet file. Returns a vector of `NamedTuple`s with
fields `(name, code, lons, lats)`, keeping every `LineString` feature.

Empty when the file is absent, so a caller without a coastline file gets "no rings"
rather than a crash — which `is_point_on_land_geodata` treats as "unknown, not land".
"""
function load_coastline_polygons(filepath::AbstractString)
    is_cached(filepath) || return NamedTuple[]
    endswith(filepath, ".parquet") || error(
        "Unsupported coastline file format: $(filepath). Only .parquet is read here; " *
        "the legacy .dat reader was never implemented and has been removed rather than " *
        "left to fail at an undefined call.")

    ds = geoload(filepath; backend = :geoparquet)
    # Variables are keyed by String: :geometry is a KeyError, and the geometry column
    # is the one variable the GeoParquet backend keeps, so this is the column that
    # carries the shoreline geometry.
    haskey(ds.variables, "geometry") ||
        error("$(filepath) has no 'geometry' variable; it is not a coastline file. " *
              "Variables present: $(sort(collect(keys(ds.variables)))).")

    polys = NamedTuple[]
    for i in eachindex(ds.variables["geometry"].data)
        geom = ds.variables["geometry"].data[i]
        if geom isa GeoInterface.LineString
            coords = collect(geom)
            push!(polys, (
                name = string(i in eachindex(ds.variables["name"].data) ?
                              ds.variables["name"].data[i] : ""),
                code = Symbol(i in eachindex(ds.variables["code"].data) ?
                              ds.variables["code"].data[i] : "unknown"),
                lons = Float64[pt[1] for pt in coords],
                lats = Float64[pt[2] for pt in coords],
            ))
        end
    end
    return polys
end

"""
    is_point_on_land_geodata(
        lon::Real, lat::Real;
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        coastline_file::Union{Nothing, AbstractString} = nothing
    ) -> Bool

Check if a point is on land using coastline rings.

Natural Earth's coastline is a set of open LineStrings, and the ray-casting test used
below needs a ring: applied to an open line it reports land for the half-plane the
shoreline faces, which is exactly the wrong answer for a land mask. Rings are therefore
required to be closed, and an open input returns `false` (unknown, not land) rather than
answering the wrong question confidently.
"""
function is_point_on_land_geodata(
    lon::Real, lat::Real;
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    coastline_file::Union{Nothing, AbstractString} = nothing
)
    polys = if !isnothing(coastline)
        coastline
    elseif !isnothing(coastline_file)
        load_coastline_polygons(coastline_file)
    else
        NamedTuple[]
    end

    isempty(polys) && return false

    # A ring is closed only when its first and last vertex match; anything else is an open
    # line, and ray-casting it as a polygon is a silent wrong answer.
    any(p -> p.lons[begin] != p.lons[end] || p.lats[begin] != p.lats[end], polys) && return false

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
    is_marine_water_geodata(
        lon::Real, lat::Real;
        bathymetry::Union{Function, GeoDataset, AbstractString, Nothing} = nothing,
        coastlines::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        min_seabed_depth::Real = 0.0
    ) -> Bool

Check whether a coordinate is marine water: not land, and (when a bathymetry source is
given) at least `min_seabed_depth` of water below it.

With no bathymetry source the depth check is skipped rather than assumed, so the answer
reflects only what was actually asked.
"""
function is_marine_water_geodata(
    lon::Real,
    lat::Real;
    bathymetry::Union{Function, GeoDataset, AbstractString, Nothing} = nothing,
    coastlines::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    min_seabed_depth::Real = 0.0,
)
    # 1. Coastline land exclusion
    if is_point_on_land_geodata(lon, lat; coastline = coastlines)
        return false
    end

    isnothing(bathymetry) && return true    # depth not requested: answer the land test only

    bathy_fn = bathymetry isa GeoDataset ? get_bathymetry_interpolator(bathymetry) :
               bathymetry isa AbstractString ? get_bathymetry_interpolator(geoload(bathymetry)) :
               bathymetry

    z_bed = bathy_fn(Float64(lon), Float64(lat))
    h_threshold = -max(0.0, Float64(min_seabed_depth))
    return z_bed <= h_threshold
end
