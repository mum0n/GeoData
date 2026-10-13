"""
    landmask.jl

Point-in-polygon and land/ocean tests over geometry the caller already holds.

This is generic geometry, not data acquisition: it decides whether a coordinate falls
inside a ring. It lives in `GeoData` because a workflow that has *already* downloaded a
coastline — or has coastline rings from any source at all — should not have to import a
package whose job is downloading. The boundary that matters is not "coastlines" but
"network access": `GeoData` contains no network code, and the fetcher that produces a
coastline file stays in `GeoDataSources`.

It also holds no bundled geography. A ring set is data, and a data set that covers one
region rather than every region would make every consumer of this generic layer carry that
region; such sets live where regional shipped data lives.
"""

using GeoData: GeoDataset

export point_in_polygon, is_point_on_land_geodata, is_marine_water_geodata,
    load_coastline_polygons, is_point_on_land, DEFAULT_COASTLINE_PATH

"""
    DEFAULT_COASTLINE_PATH

Where a coastline file is looked for when a caller does not name one: the run's own
containerised `inputs/` directory, relative to the working directory.

It is a default rather than a compiled-in path so a run that stores its coastline elsewhere
says so; and it is *absent* rather than a bundled fallback, because a hand-drawn regional
outline masks most of a domain's land and looks like data. An absent file means the land test
answers "not land", which is what `is_point_on_land_geodata` documents for an empty ring set:
unknown, not assumed.

This layer owns no bundled geography. A run that wants a regional ring set without fetching
one asks the package that ships it -- `GeoDataSources.SCOTIAN_SHELF_COASTLINE` -- and passes
it, which keeps the generic layer free of any project's study region.
"""
const DEFAULT_COASTLINE_PATH = joinpath("inputs", "coastline.parquet")

"""
    is_point_on_land(lon::Real, lat::Real;
                     coastline = nothing,
                     coastline_file = DEFAULT_COASTLINE_PATH) -> Bool

Whether a coordinate is on land, reading the coastline file when the caller does not name
one and otherwise using the rings it passes.

This is `is_point_on_land_geodata` with the default file location applied, and it exists
because every caller that does its own land test would otherwise repeat that lookup: the
location of a run's coastline file is one decision, and repeating it is how two callers end
up disagreeing about whether a point is land.

Precedence, and why each step is what it is:

1. the rings the caller passed in -- the caller knows best what coastline it is using;
2. the rings read from `coastline_file`, when that file exists;
3. nothing -- and the answer is "not land", which is `is_point_on_land_geodata`'s documented
   answer for no rings: unknown, not assumed.

There is deliberately no step 4 with a bundled regional ring set. This layer is generic, and a
regional geography in it would make every consumer of the generic layer carry this project's
study region. A regional default is the *application's* decision: pass
`GeoDataSources.SCOTIAN_SHELF_COASTLINE` if that is what a run wants, or fetch a real
shoreline with `GeoDataSources.fetch_natural_earth_coastline`.
"""
function is_point_on_land(lon::Real, lat::Real;
                          coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
                          coastline_file::Union{Nothing, AbstractString} = DEFAULT_COASTLINE_PATH)
    if !isnothing(coastline)
        return is_point_on_land_geodata(lon, lat; coastline = coastline)
    end
    rings = isnothing(coastline_file) ? NamedTuple[] : load_coastline_polygons(String(coastline_file))
    return is_point_on_land_geodata(lon, lat; coastline = rings)
end

"""
    load_coastline_polygons(; coastline_path = DEFAULT_COASTLINE_PATH) -> Vector{NamedTuple}

Coastline rings from the default location when the caller does not name a file.

Absent file means no rings, and no rings means "unknown, not land" to
`is_point_on_land_geodata` -- which is the honest answer, and better than a coarse regional
outline that silently rejects land a run never asked about. A run that wants a regional ring
set without fetching one asks the package that ships it: `GeoDataSources`.
`SCOTIAN_SHELF_COASTLINE`.

The positional form `load_coastline_polygons(path)` is unchanged: naming a file means the file
was meant to be read.
"""
function load_coastline_polygons(; coastline_path::Union{Nothing, AbstractString} = DEFAULT_COASTLINE_PATH)
    isnothing(coastline_path) && return NamedTuple[]
    return load_coastline_polygons(String(coastline_path))
end

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
    _ring_named_tuple(name, code, coords) -> NamedTuple

One coastline ring as `(name, code, lons, lats)`, from GeoInterface's coordinate
representation. A ring is the first two coordinates of every point: a third is height, not
position along a shoreline.
"""
function _ring_named_tuple(name, code, coords)
    lons = Float64[Float64(pt[1]) for pt in coords]
    lats = Float64[Float64(pt[2]) for pt in coords]
    return (name = name, code = code, lons = lons, lats = lats)
end

"""
    load_coastline_polygons(filepath::AbstractString) -> Vector{NamedTuple}

Load coastline rings from a GeoParquet file. Returns a vector of `NamedTuple`s with
fields `(name, code, lons, lats)`, keeping every `LineString` feature, and every ring of a
`MultiLineString`.

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

    # The geometry column comes back as well-known binary, and GeoInterface's traits read it
    # directly -- the same idiom GeoData's own GeoParquet backend documents. The trait is what
    # must be tested: `GeoInterface.LineString` is a *trait type*, not the concrete geometry,
    # so `geom isa GeoInterface.LineString` is false for every feature and this function used
    # to return no rings at all, for every file, without failing. A land mask that silently
    # masks nothing is worse than one that errors.
    geoms = ds.variables["geometry"].data
    names = haskey(ds.variables, "name") ? ds.variables["name"].data : fill("", length(geoms))
    codes = haskey(ds.variables, "code") ? ds.variables["code"].data : fill("unknown", length(geoms))

    polys = NamedTuple[]
    for i in eachindex(geoms)
        geom = geoms[i]
        name = string(i in eachindex(names) ? names[i] : "")
        code = Symbol(i in eachindex(codes) ? codes[i] : "unknown")

        trait = GeoInterface.trait(geom)
        if trait isa GeoInterface.LineStringTrait
            push!(polys, _ring_named_tuple(name, code, GeoInterface.coordinates(geom)))
        elseif trait isa GeoInterface.MultiLineStringTrait
            # A MultiLineString is several shorelines in one feature, and `coordinates` nests
            # one level for the parts, so each becomes its own ring: the land test needs one
            # ring per curve. (GeoInterface 1.6 has no `geometries` accessor for this.)
            for part in GeoInterface.coordinates(geom)
                push!(polys, _ring_named_tuple(name, code, part))
            end
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
