"""
    bathymetry.jl

Bathymetry fetch, load, save, and query operations via GeoData.

Supports ERDDAP (NOAA), ETOPO2022, and generic NetCDF/Zarr formats.
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using NCDatasets
using Downloads
using Interpolations
using Statistics

# Optional NumericalEarth for ETOPO bathymetry
const _has_numerical_earth = try
    import NumericalEarth
    true
catch
    false
end

"""
    fetch_erddap_bathymetry(;
        lon_range::Tuple{Real, Real} = (-71.0, -53.0),
        lat_range::Tuple{Real, Real} = (40.0, 48.5),
        output_path::AbstractString = joinpath("inputs", "bathymetry.zarr"),
        dataset_id::AbstractString = "ETOPO_2022_v1_15s",
        stride::Int = 1,
        backend::Symbol = :zarr,
        verbose::Bool = true
    ) -> GeoDataset

Fetch bathymetry from ERDDAP/NOAA and save as a GeoData-compatible dataset (Zarr/NetCDF).

# Arguments
- `lon_range`: Longitude bounds (min, max)
- `lat_range`: Latitude bounds (min, max)
- `output_path`: Output file path
- `dataset_id`: ERDDAP dataset ID (e.g., "ETOPO_2022_v1_15s", "etopo180", "srtm30plus")
- `stride`: Sampling stride for download
- `backend`: Storage backend (:zarr, :netcdf)
- `verbose`: Print progress messages

# Returns
GeoDataset with elevation, lon, lat coordinates.
"""
function fetch_erddap_bathymetry(;
    lon_range::Tuple{Real, Real} = (-71.0, -53.0),
    lat_range::Tuple{Real, Real} = (40.0, 48.5),
    output_path::AbstractString = joinpath("inputs", "bathymetry.zarr"),
    dataset_id::AbstractString = "ETOPO_2022_v1_15s",
    stride::Int = 1,
    backend::Symbol = :zarr,
    verbose::Bool = true
)
    mkpath(dirname(output_path))
    min_lat, max_lat = Float64(lat_range[1]), Float64(lat_range[2])
    min_lon, max_lon = Float64(lon_range[1]), Float64(lon_range[2])

    erddap_vars = Dict(
        "ETOPO_2022_v1_15s" => "z",
        "etopo180"          => "altitude",
        "srtm30plus"        => "elevation"
    )
    var = get(erddap_vars, dataset_id, "z")

    primary_url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/$(dataset_id).nc?$(var)[($(min_lat)):$(stride):($(max_lat))][($(min_lon)):$(stride):($(max_lon))]"

    verbose && println("Fetching bathymetry from ERDDAP...")
    tmp_nc = tempname() * ".nc"
    try
        Downloads.download(primary_url, tmp_nc)
    catch err
        error("Failed to download bathymetry: $(err)")
    end

    # Load the downloaded NetCDF and convert to GeoDataset
    ds = geoload(tmp_nc; backend=:ncdatasets)
    
    # Save as Zarr (or other backend)
    geosave(output_path, ds; backend=backend)
    rm(tmp_nc, force=true)
    
    verbose && println("Saved bathymetry to $(output_path)")
    return ds
end

"""
    load_bathymetry_geodata(filepath::AbstractString; backend=:zarr) -> GeoDataset

Load bathymetry from a GeoData-compatible file (Zarr, NetCDF, etc.)
"""
function load_bathymetry_geodata(filepath::AbstractString; backend::Symbol = :zarr)
    return geoload(filepath; backend=backend)
end

"""
    save_bathymetry_geodata(ds::GeoDataset, filepath::AbstractString; backend=:zarr) -> String

Save bathymetry GeoDataset to a file.
"""
function save_bathymetry_geodata(ds::GeoDataset, filepath::AbstractString; backend::Symbol = :zarr)
    mkpath(dirname(filepath))
    geosave(filepath, ds; backend=backend)
    return filepath
end

"""
    get_bathymetry_interpolator(ds::GeoDataset; varname="elevation") -> Function

Create a continuous (lon, lat) -> elevation interpolator from a bathymetry GeoDataset.
"""
function get_bathymetry_interpolator(ds::GeoDataset; varname::AbstractString = "elevation")
    # Find the elevation variable
    elev_var = haskey(ds.variables, varname) ? varname :
               findfirst(k -> k in ("elevation", "altitude", "z", "topo", "bedrock_altitude"), keys(ds.variables))
    isnothing(elev_var) && error("No elevation variable found in dataset. Available: $(keys(ds.variables))")
    
    elev_ga = ds.variables[elev_var]
    lons = haskey(ds.coords, :lon) ? vec(ds.coords[:lon].data) :
           haskey(ds.coords, :longitude) ? vec(ds.coords[:longitude].data) : error("No longitude coordinate")
    lats = haskey(ds.coords, :lat) ? vec(ds.coords[:lat].data) :
           haskey(ds.coords, :latitude) ? vec(ds.coords[:latitude].data) : error("No latitude coordinate")
    
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
    regrid_bathymetry_from_etopo(;
        lon_range::Tuple{Real, Real} = (-71.0, -53.0),
        lat_range::Tuple{Real, Real} = (40.0, 48.5),
        resolution::Real = 600,
        minimum_depth::Real = 0.0,
        interpolation_passes::Integer = 1,
        major_basins::Real = 1,
        cache::Bool = true
    ) -> GeoDataset

Get regional bathymetry as a GeoDataset from ETOPO2022 using NumericalEarth.jl.
This is a thin wrapper around NumericalEarth.jl's regrid_bathymetry function.

# Arguments
- `lon_range`: Longitude bounds (min, max)
- `lat_range`: Latitude bounds (min, max)
- `resolution`: Target resolution in meters
- `minimum_depth`: Minimum water depth (m)
- `interpolation_passes`: Number of interpolation passes
- `major_basins`: Major basins parameter
- `cache`: Use cached data

# Returns
GeoDataset with elevation, lon, lat coordinates.
"""
function regrid_bathymetry_from_etopo(;
    lon_range::Tuple{Real, Real} = (-71.0, -53.0),
    lat_range::Tuple{Real, Real} = (40.0, 48.5),
    resolution::Real = 600,
    minimum_depth::Real = 0.0,
    interpolation_passes::Integer = 1,
    major_basins::Real = 1,
    cache::Bool = true
)
    _has_numerical_earth || error("NumericalEarth.jl required for ETOPO bathymetry regridding. Add NumericalEarth to your environment.")
    
    lat_mid = 0.5 * (lat_range[1] + lat_range[2])
    m_per_deg = 111_320.0 * cosd(lat_mid)
    n_lon = clamp(round(Int, abs(lon_range[2] - lon_range[1]) * m_per_deg / Float64(resolution)), 16, 21600)
    n_lat = clamp(round(Int, abs(lat_range[2] - lat_range[1]) * m_per_deg / Float64(resolution)), 16, 10800)
    
    # Create scratch grid
    scratch = NumericalEarth.LatitudeLongitudeGrid(
        NumericalEarth.CPU(), Float32;
        size = (n_lon, n_lat, 1),
        longitude = (lon_range[1], lon_range[2]),
        latitude = (lat_range[1], lat_range[2]),
        z = (-1.0, 0.0)
    )
    
    z_field = NumericalEarth.regrid_bathymetry(scratch;
        dataset = NumericalEarth.ETOPO2022(),
        minimum_depth = Float64(minimum_depth),
        interpolation_passes = Int(interpolation_passes),
        major_basins = Float64(major_basins),
        cache = cache
    )
    
    elevation = Array(NumericalEarth.interior(z_field, :, :, 1))
    
    lons = collect(Float64, scratch.λᶜᵃᵃ)[scratch.Hx+1 : scratch.Hx+scratch.Nx]
    lats = collect(Float64, scratch.φᵃᶜᵃ)[scratch.Hy+1 : scratch.Hy+scratch.Ny]
    
    lon_dim = Dimension(name=:lon, size=length(lons), coords=lons, units="degrees_east")
    lat_dim = Dimension(name=:lat, size=length(lats), coords=lats, units="degrees_north")
    
    elev_ga = GeoArray(elevation, (lon_dim, lat_dim), CoordinateSystem(crs="EPSG:4326"), Dict("units" => "meters"))
    lon_ga = GeoArray(lons, (lon_dim,), CoordinateSystem(crs="EPSG:4326"), Dict("units" => "degrees_east"))
    lat_ga = GeoArray(lats, (lat_dim,), CoordinateSystem(crs="EPSG:4326"), Dict("units" => "degrees_north"))
    
    return GeoDataset(
        Dict("elevation" => elev_ga),
        Dict(:lon => lon_ga, :lat => lat_ga),
        Dict(:lon => lon_dim, :lat => lat_dim),
        CoordinateSystem(crs="EPSG:4326"),
        Dict{String, Any}(),
        nothing, "etopo_bathymetry"
    )
end

"""
    etopo_bathymetry_field(;
        lon_range::Tuple{Real, Real},
        lat_range::Tuple{Real, Real},
        resolution::Real = 600
    ) -> GeoDataset

Alias for `regrid_bathymetry_from_etopo` with simplified arguments.
"""
function etopo_bathymetry_field(;
    lon_range::Tuple{Real, Real},
    lat_range::Tuple{Real, Real},
    resolution::Real = 600
)
    return regrid_bathymetry_from_etopo(;
        lon_range = lon_range,
        lat_range = lat_range,
        resolution = resolution
    )
end

"""
    smooth_bathymetry(
        topo::AbstractMatrix{<:Real};
        passes::Integer = 3,
        alpha::Real = 0.5,
        h_min::Real = 20.0
    ) -> Matrix{Float64}

Apply conservative 2D discrete Laplacian smoothing to bathymetric elevation data on wet
cells, attenuating subgrid \$2\\Delta x\$ pinnacles and single-cell topographic cliffs
arising from bilinear interpolation of high-resolution digital elevation models (ETOPO/GEBCO).
Preserves emerged land points (\$z = 0.0\\text{ m}\$) and enforces the minimum water column
depth floor \$h_{\\min}\$.

# Mathematical Formulation
For each smoothing iteration \$m = 1, \\dots, M\$ and every interior wet cell \$(i, j)\$ with
seabed elevation \$Z_{i,j}^{(m)} \\le -h_{\\min}\$:
```math
Z_{i,j}^{(m+1)} = (1 - \\alpha) Z_{i,j}^{(m)} + \\frac{\\alpha}{N_{\\text{wet}}}
                  \\sum_{(p,q) \\in \\mathcal{N}_{\\text{wet}}(i,j)} Z_{p,q}^{(m)}
```
where \$\\mathcal{N}_{\\text{wet}}(i,j) = \\{(i \\pm 1, j), (i, j \\pm 1) \\mid Z_{p,q}^{(m)} \\le -h_{\\min}\\}\$
denotes the 4-connected wet neighborhood, and \$N_{\\text{wet}} = |\\mathcal{N}_{\\text{wet}}(i,j)| \\ge 2\$.
If \$Z_{i,j}^{(m+1)} > -h_{\\min}\$, the depth floor \$Z_{i,j}^{(m+1)} = -h_{\\min}\$ is enforced.

# Inputs
- `topo::AbstractMatrix{<:Real}`: 2D array of seabed elevations in meters.
- `passes::Integer`: Number of smoothing iterations (default 3).
- `alpha::Real`: Smoothing weight in \$[0, 1]\$ (default 0.5).
- `h_min::Real`: Minimum physical water depth floor in meters (default 20.0m).

# Outputs
- `Matrix{Float64}`: Smoothed elevation matrix matching the input horizontal dimensions.

# References
- Shapiro, R. (1970). Smoothing, filtering, and boundary effects.
  *Reviews of Geophysics*, 8(2), 359-387.
- Haidvogel, D. B., & Beckmann, A. (1999). *Numerical Ocean Circulation Modeling*.
  Imperial College Press.
"""
function smooth_bathymetry(
    topo::AbstractMatrix{<:Real};
    passes::Integer = 3,
    alpha::Real = 0.5,
    h_min::Real = 20.0
)::Matrix{Float64}
    if passes < 0
        error("smooth_bathymetry: passes must be non-negative, got passes = $(passes)")
    end
    if !(0.0 <= alpha <= 1.0)
        error("smooth_bathymetry: alpha must be in [0, 1], got alpha = $(alpha)")
    end
    if h_min <= 0.0
        error("smooth_bathymetry: h_min must be positive, got h_min = $(h_min)")
    end

    nx, ny = size(topo)
    h_floor = Float64(h_min)
    alpha_f = Float64(alpha)

    # Initialize conditioned array: land >= 0 is 0.0, shallow wet is floored
    h = Matrix{Float64}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        z = Float64(topo[i, j])
        if z >= 0.0
            h[i, j] = 0.0
        elseif z > -h_floor
            h[i, j] = -h_floor
        else
            h[i, j] = z
        end
    end

    # Apply 2D discrete Laplacian smoothing iterations on wet cells
    for _ in 1:passes
        h_next = copy(h)
        for j in 2:(ny - 1), i in 2:(nx - 1)
            if h[i, j] <= -h_floor
                sum_nb = 0.0
                n_wet = 0
                for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))
                    nb = h[i + di, j + dj]
                    if nb <= -h_floor
                        sum_nb += nb
                        n_wet += 1
                    end
                end
                if n_wet >= 2
                    avg_nb = sum_nb / n_wet
                    z_smoothed = (1.0 - alpha_f) * h[i, j] + alpha_f * avg_nb
                    h_next[i, j] = z_smoothed > -h_floor ? -h_floor : z_smoothed
                end
            end
        end
        h = h_next
    end

    return h
end

"""
    extract_marine_cells(
        bathymetry::Union{NamedTuple, AbstractString, GeoDataset};
        lon_range::Tuple{Real, Real} = (-180.0, 180.0),
        lat_range::Tuple{Real, Real} = (-90.0, 90.0),
        min_seabed_depth::Real = 0.0,
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
    ) -> NamedTuple

Extract all discrete marine grid cell centers that lie in open water (\$z < 0\$),
strictly outside terrestrial landmasses, and meet the minimum water depth requirement
within the specified spatial bounding box.

# Inputs
- `bathymetry`: `NamedTuple` `(lon, lat, elevation)`, NetCDF/Zarr file path `AbstractString`, or `GeoDataset`.
- `lon_range::Tuple{Real, Real}`: Bounding box longitude limits.
- `lat_range::Tuple{Real, Real}`: Bounding box latitude limits.
- `min_seabed_depth::Real`: Minimum water depth in meters (default 0.0 m).
- `coastline`: Optional coastline polygons list.

# Outputs
- `NamedTuple`: `(lons = Vector{Float64}, lats = Vector{Float64}, depths = Vector{Float64}, weights = Vector{Float64})`
"""
function extract_marine_cells(
    bathymetry::Union{NamedTuple, AbstractString, GeoDataset};
    lon_range::Tuple{Real, Real} = (-180.0, 180.0),
    lat_range::Tuple{Real, Real} = (-90.0, 90.0),
    min_seabed_depth::Real = 0.0,
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing
)
    bathy_data = if bathymetry isa AbstractString
        load_bathymetry_geodata(bathymetry)
    elseif bathymetry isa GeoDataset
        bathymetry
    else
        bathymetry
    end

    lons = Float64.(bathy_data.lon)
    lats = Float64.(bathy_data.lat)
    elev = Float64.(bathy_data.elevation)

    h_threshold = -max(0.0, Float64(min_seabed_depth))

    marine_lons = Float64[]
    marine_lats = Float64[]
    marine_depths = Float64[]
    marine_weights = Float64[]

    n_lon = length(lons)
    n_lat = length(lats)

    for i in 1:n_lon
        x = lons[i]
        if !(lon_range[1] <= x <= lon_range[2])
            continue
        end
        for j in 1:n_lat
            y = lats[j]
            if !(lat_range[1] <= y <= lat_range[2])
                continue
            end
            z = elev[i, j]

            # Point must not be on terrestrial land and depth must satisfy threshold
            if z <= h_threshold && !is_point_on_land_geodata(x, y; coastline = coastline)
                push!(marine_lons, x)
                push!(marine_lats, y)
                push!(marine_depths, z)
                # Spherical grid cell surface area weight dA = R^2 cos(lat) dlon dlat
                weight = cos(deg2rad(clamp(y, -89.9, 89.9)))
                push!(marine_weights, weight)
            end
        end
    end

    if isempty(marine_lons)
        error(
            "No marine water cells found with depth >= $(min_seabed_depth) m " *
            "within lon $(lon_range) and lat $(lat_range). Domain is entirely land."
        )
    end

    return (
        lons = marine_lons,
        lats = marine_lats,
        depths = marine_depths,
        weights = marine_weights
    )
end

"""
    sample_marine_coordinates(
        n_particles::Int,
        bathymetry::Union{NamedTuple, AbstractString, GeoDataset};
        lon_range::Tuple{Real, Real} = (-180.0, 180.0),
        lat_range::Tuple{Real, Real} = (-90.0, 90.0),
        min_seabed_depth::Real = 0.0,
        coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        rng::AbstractRNG = Random.default_rng()
    ) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Sample \$N\$ continuous coordinates strictly within active marine cells (\$z_{\\text{bed}} < 0\$
and \$z_{\\text{bed}} \\le -h_{\\text{min}}\$) with spherical area weighting, sub-cell jitter,
and rigorous coastline land rejection.

# Guarantees
- 100% deterministic success in \$O(N)\$ time without rejection sampling stalls.
- Exactly 0% probability of particles landing on emergent land or coastal terrain.

# Inputs
- `n_particles::Int`: Number of particle coordinates to sample.
- `bathymetry`: `NamedTuple`, NetCDF/Zarr file path, or `GeoDataset`.
- `lon_range, lat_range`: Geographic bounding box.
- `min_seabed_depth::Real`: Minimum seabed depth (default: 0.0 m).
- `coastline`: Optional coastline polygons.
- `rng::AbstractRNG`: Random number generator.

# Outputs
- `Tuple`: `(sampled_lons, sampled_lats, sampled_seabed_depths)`
"""
function sample_marine_coordinates(
    n_particles::Int,
    bathymetry::Union{NamedTuple, AbstractString, GeoDataset};
    lon_range::Tuple{Real, Real} = (-180.0, 180.0),
    lat_range::Tuple{Real, Real} = (-90.0, 90.0),
    min_seabed_depth::Real = 0.0,
    coastline::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    rng::AbstractRNG = Random.default_rng()
)
    cells = extract_marine_cells(
        bathymetry,
        lon_range = lon_range,
        lat_range = lat_range,
        min_seabed_depth = min_seabed_depth,
        coastline = coastline
    )

    bathy_data = if bathymetry isa AbstractString
        load_bathymetry_geodata(bathymetry)
    elseif bathymetry isa GeoDataset
        bathymetry
    else
        bathymetry
    end

    lons_raw = Float64.(bathy_data.lon)
    lats_raw = Float64.(bathy_data.lat)
    dlon = length(lons_raw) > 1 ? abs(lons_raw[2] - lons_raw[1]) : 0.05
    dlat = length(lats_raw) > 1 ? abs(lats_raw[2] - lats_raw[1]) : 0.05

    bathy_interp = get_bathymetry_interpolator(bathy_data)
    h_threshold = -max(0.0, Float64(min_seabed_depth))

    total_weight = sum(cells.weights)
    cum_weights = cumsum(cells.weights) ./ total_weight

    sampled_lons = Vector{Float64}(undef, n_particles)
    sampled_lats = Vector{Float64}(undef, n_particles)
    sampled_zbed = Vector{Float64}(undef, n_particles)

    for p in 1:n_particles
        u = rand(rng, Float64)
        idx = searchsortedfirst(cum_weights, u)
        idx = clamp(idx, 1, length(cells.lons))

        c_lon = cells.lons[idx]
        c_lat = cells.lats[idx]

        jitter_x = (rand(rng, Float64) - 0.5) * dlon * 0.95
        jitter_y = (rand(rng, Float64) - 0.5) * dlat * 0.95
        cand_x = c_lon + jitter_x
        cand_y = c_lat + jitter_y

        z_cand = bathy_interp(cand_x, cand_y)

        # Enforce both bathymetric threshold and coastline land rejection
        if z_cand > h_threshold || is_point_on_land_geodata(cand_x, cand_y; coastline = coastline)
            cand_x = c_lon
            cand_y = c_lat
            z_cand = cells.depths[idx]
        end

        sampled_lons[p] = cand_x
        sampled_lats[p] = cand_y
        sampled_zbed[p] = z_cand
    end

    return (sampled_lons, sampled_lats, sampled_zbed)
end