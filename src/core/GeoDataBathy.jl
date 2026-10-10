"""
    bathymetry.jl (analysis)

Bathymetry analysis: continuous interpolators, smoothing, and marine-cell extraction.

These are all pure array math over a `GeoDataset` the caller already has in hand. They
belong in `GeoData` so a workflow can interrogate bathymetry without importing a package
whose job is *downloading* it. The fetchers that produce the files stay in
`GeoDataSources`.
"""

using Interpolations
using Statistics
using Random
using GeoData: variables_like

export smooth_bathymetry, get_bathymetry_interpolator, extract_marine_cells,
    sample_marine_cells

"""
    smooth_bathymetry(field; passes=4, alpha=0.5, h_min=0.0) -> Matrix{Float64}

Smooth a bathymetry field and floor the water depth at `h_min`.

Each pass blends every interior cell toward the weighted mean of its 3x3 neighbourhood
(centre weight 4, edge-adjacent 1, corners 1) by the factor `alpha`. Boundaries are left
untouched, so a land boundary does not bleed into the domain. A final clamp raises every
cell to at least `h_min`, which is what makes a field traversable by a model whose
minimum water depth is not zero.

# Arguments
- `field`: bathymetry shaped `(lon, lat)`, in metres, positive up or down as given
- `passes`: smoothing iterations
- `alpha`: blend factor per pass, 0 = unchanged, 1 = full replacement
- `h_min`: minimum water depth in metres (0 disables the floor)
"""
function smooth_bathymetry(field::AbstractMatrix; passes::Integer = 4, alpha::Real = 0.5,
                           h_min::Real = 0.0)
    smoothed = Float64.(collect(field))
    n_lon, n_lat = size(smoothed)
    (n_lon < 3 || n_lat < 3) && return smoothed

    for pass in 1:passes
        next = copy(smoothed)
        for j in 2:(n_lat - 1), i in 2:(n_lon - 1)
            block = smoothed[(i - 1):(i + 1), (j - 1):(j + 1)]
            mean = (4.0 * smoothed[i, j] + sum(block) - smoothed[i, j]) / 12.0
            next[i, j] = (1.0 - alpha) * smoothed[i, j] + alpha * mean
        end
        smoothed = next
    end

    if h_min > 0.0
        smoothed = max.(smoothed, Float64(h_min))
    end
    return smoothed
end

"""
    get_bathymetry_interpolator(ds::GeoDataset; varname::AbstractString = "elevation") -> Function

Create a continuous `(lon, lat) -> elevation` interpolator from a bathymetry GeoDataset.

Values outside the source domain are held at the nearest edge value (`Flat`
extrapolation) rather than failing, so a query slightly outside the survey box still
returns a real depth.

# Arguments
- `ds`: bathymetry GeoDataset
- `varname`: the elevation variable to use, or a keyword prefix to search for
"""
function get_bathymetry_interpolator(ds::GeoDataset; varname::AbstractString = "elevation")
    elev_var = haskey(ds.variables, varname) ? varname :
               variables_like(ds, ["elevation", "altitude", "z", "topo", "bedrock_altitude"])

    elev_ga = ds.variables[elev_var]
    (lons, _) = axis(ds, :lon)
    (lats, _) = axis(ds, :lat)

    elev = elev_ga.data

    # Ensure (lon, lat) ordering
    if size(elev) == (length(lats), length(lons))
        elev = permutedims(elev, (2, 1))
    end

    itp = interpolate((lons, lats), elev, Gridded(Linear()))
    itp_flat = extrapolate(itp, Flat())

    return (lon::Real, lat::Real) -> Float64(itp_flat(Float64(lon), Float64(lat)))
end

"""
    extract_marine_cells(
        ds::GeoDataset;
        elevation_var::AbstractString = "elevation",
        min_seabed_depth::Real = 0.0,
        coastlines::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
    ) -> NamedTuple

Extract marine grid cells from a bathymetry dataset.

Returns `(mask, indices, elevation, lons, lats)` where `mask` is a `BitMatrix` (true =
water) of shape `(lat, lon)`, `indices` the linear cell indices of water cells, and
`elevation`/`lons`/`lats` their values, row-major over latitude.

Two exclusions are applied, in order:

1. **Depth** - a cell is water only if `z <= -max(0, min_seabed_depth)`, so a positive
   `min_seabed_depth` requires at least that much water and zero admits sea level.
2. **Land** - when `coastlines` is given, cells inside any ring are removed. Skip this
   only if you have no land mask: it is what keeps a particle from spawning on a headland
   that the gridded elevation rounds to zero.

Orientation is read from the variable's own dimension records, not from `size`: on a square
grid both spellings have the same shape, and guessing from size is how a transposed read
stays silent. Descending axes are sorted, with the data permuted alongside them so a
coordinate and the values it describes stay tied.

# Arguments
- `ds`: bathymetry GeoDataset
- `elevation_var`: variable holding elevation (metres, positive down or up)
- `min_seabed_depth`: minimum water depth to accept, in metres
- `coastlines`: optional closed coastline rings, as returned by
  `load_coastline_polygons`. Cells inside any ring are excluded by `_exclude_land!`
"""
function extract_marine_cells(
    ds::GeoDataset;
    elevation_var::AbstractString = "elevation",
    min_seabed_depth::Real = 0.0,
    coastlines::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
)
    ga = ds.variables[elevation_var]
    elev = Float64.(ga.data)

    # Orientation is resolved from the variable's own dimension records, not inferred
    # from `size`: on a square grid both spellings have the same shape, and guessing
    # from size is how a transposed read stays silent. A 2-D field is the contract;
    # anything else is refused rather than flattened.
    (pos_lat, pos_lon) = dim_permutation(ga.dims, :lat, :lon)
    if ndims(elev) != 2 || length(unique([pos_lat, pos_lon])) != 2 ||
       sort([pos_lat, pos_lon]) != collect(1:ndims(elev))
        error("Variable $(repr(elevation_var)) is $(ndims(elev))-D with dims " *
              "$([standardize_dimension_name(d.name) for d in ga.dims]). A 2-D field " *
              "spanning :lon and :lat is what this function reads; anything else is " *
              "refused rather than flattened.")
    end

    (lons, _) = axis(ds, :lon)
    (lats, _) = axis(ds, :lat)
    if length(lons) != size(elev, pos_lon) || length(lats) != size(elev, pos_lat)
        error("Variable $(repr(elevation_var)) is $(size(elev)) with dims " *
              "$([standardize_dimension_name(d.name) for d in ga.dims]), but :lon has " *
              "$(length(lons)) values and :lat has $(length(lats)). A dimension record and " *
              "its coordinate vector must agree.")
    end

    # One canonical orientation for everything downstream: (lat, lon).
    elev_t = permutedims(elev, (pos_lat, pos_lon))

    # Sort both axes, permuting the data with them so a coordinate and the values it
    # describes stay tied. Sorting the axes alone would leave `elev_t` describing the
    # old order, which is the wrong answer delivered quietly.
    p_lon = sortperm(lons)
    p_lat = sortperm(lats)
    lons = lons[p_lon]
    lats = lats[p_lat]
    elev_t = elev_t[p_lat, p_lon]

    n_lon, n_lat = length(lons), length(lats)

    h_threshold = -max(0.0, Float64(min_seabed_depth))
    mask = BitMatrix(elev_t .<= h_threshold)

    if coastlines !== nothing
        _exclude_land!(mask, lons, lats, coastlines)
    end

    marine_indices = CartesianIndex{2}[]
    marine_elev = Float64[]
    marine_lons = Float64[]
    marine_lats = Float64[]
    for j in 1:n_lat, i in 1:n_lon
        if mask[j, i]
            push!(marine_indices, CartesianIndex(j, i))
            push!(marine_elev, elev_t[j, i])
            push!(marine_lons, lons[i])
            push!(marine_lats, lats[j])
        end
    end

    return (
        mask = mask,
        indices = marine_indices,
        elevation = marine_elev,
        lons = marine_lons,
        lats = marine_lats,
    )
end

"""
    _exclude_land!(mask, lons, lats, coastlines) -> BitMatrix

Clear `mask` inside each coastline ring, in place.

`mask` is `(lat, lon)` and `lons`/`lats` are ascending, both call-sorted together with the
data they index. Each ring's bounding box becomes an index range by binary search, so only
the cells that could be inside a ring get the ray-cast, rather than every cell against every
ring: the cost is the number of cells a ring could contain, not cells x rings x vertices.

Rings that are not closed are skipped. Ray-casting an open LineString reports land for the
half-plane it faces, which is a wrong answer delivered confidently.
"""
function _exclude_land!(mask::BitMatrix, lons::Vector{Float64}, lats::Vector{Float64},
                        coastlines::AbstractVector{<:NamedTuple})
    for poly in coastlines
        (poly.lons[begin] == poly.lons[end] && poly.lats[begin] == poly.lats[end]) || continue
        i0 = max(1, searchsortedfirst(lons, minimum(poly.lons)))
        i1 = min(length(lons), searchsortedlast(lons, maximum(poly.lons)))
        j0 = max(1, searchsortedfirst(lats, minimum(poly.lats)))
        j1 = min(length(lats), searchsortedlast(lats, maximum(poly.lats)))
        for j in j0:j1, i in i0:i1
            if mask[j, i] && point_in_polygon(lons[i], lats[j], poly.lons, poly.lats)
                mask[j, i] = false
            end
        end
    end
    return mask
end

"""
    sample_marine_cells(
        n_samples::Int;
        bathymetry_ds::GeoDataset,
        elevation_var::AbstractString = "elevation",
        min_seabed_depth::Real = 0.0,
        jitter_scale::Real = 0.5,
        rng::AbstractRNG = Random.default_rng()
    ) -> NamedTuple{(:lons, :lats, :depths), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}

Sample `n_samples` coordinates from the marine cells of a bathymetry dataset, weighted
by seabed depth magnitude, with a uniform positional jitter inside each cell.

The weighted sampling is not rejection sampling: cells are drawn directly from the
normalised depth profile via inverse-CDF sampling, so the result is exactly `n_samples`
points with no retry loop and no dependence on the acceptance rate.

# Arguments
- `n_samples`: number of coordinates to draw
- `bathymetry_ds`: bathymetry GeoDataset to sample from
- `elevation_var`: variable holding elevation
- `min_seabed_depth`: minimum water depth to accept, in metres
- `jitter_scale`: positional jitter as a fraction of local grid spacing, in `[0, 1]`
- `rng`: random generator
- `coastlines`: optional coastline rings, forwarded to `extract_marine_cells`; pass these
  for the same reason `extract_marine_cells` does, so a particle is not drawn onto a
  headland the gridded elevation rounds to sea level
"""
function sample_marine_cells(
    n_samples::Int;
    bathymetry_ds::GeoDataset,
    elevation_var::AbstractString = "elevation",
    min_seabed_depth::Real = 0.0,
    jitter_scale::Real = 0.5,
    rng::AbstractRNG = Random.default_rng(),
    coastlines::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
)
    marine = extract_marine_cells(bathymetry_ds; elevation_var = elevation_var,
                                  min_seabed_depth = min_seabed_depth, coastlines = coastlines)
    isempty(marine.indices) && error("No marine cells found in the provided bathymetry. " *
                                     "With no water cells the depth threshold, the `$(elevation_var)` " *
                                     "variable, or a coastline ring excluding everything are the likely causes.")

    weights = abs.(marine.elevation)
    weights_sum = sum(weights)
    if weights_sum == 0.0
        weights = ones(length(marine.elevation))
        weights_sum = length(marine.elevation)
    end
    cum_weights = cumsum(weights ./ weights_sum)

    # Local grid spacing, for a jitter that stays inside the cell rather than jumping
    # into a neighbour when the grid is irregular.
    lons = marine.lons
    lats = marine.lats
    d_lon = _axis_spacing(lons)
    d_lat = _axis_spacing(lats)

    lons_out = Float64[]
    lats_out = Float64[]
    depths_out = Float64[]

    for _ in 1:n_samples
        u = rand(rng)
        idx = searchsortedfirst(cum_weights, u)
        idx > length(cum_weights) && (idx = length(cum_weights))
        idx = clamp(idx, 1, length(marine.elevation))

        lon = marine.lons[idx] + (2 * rand(rng) - 1) * 0.5 * jitter_scale * d_lon
        lat = marine.lats[idx] + (2 * rand(rng) - 1) * 0.5 * jitter_scale * d_lat
        push!(lons_out, lon)
        push!(lats_out, lat)
        push!(depths_out, marine.elevation[idx])
    end

    return (lons = lons_out, lats = lats_out, depths = depths_out)
end

"""
    _axis_spacing(coords::Vector{Float64}) -> Float64

The local grid spacing of a coordinate axis: the smallest gap between consecutive
values. Falls back to 1.0 for an axis with fewer than two values, which is the
dimensionless case and keeps jitter harmless.
"""
function _axis_spacing(coords::Vector{Float64})
    length(coords) < 2 && return 1.0
    gaps = abs.(diff(sort(coords)))
    gaps = gaps[gaps .> 0]
    return isempty(gaps) ? 1.0 : minimum(gaps)
end
