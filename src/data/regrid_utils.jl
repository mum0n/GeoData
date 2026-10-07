"""
    regrid_utils.jl

Generic regridding and coordinate extraction utilities for GeoData.
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using Interpolations
using Statistics
using LinearAlgebra

"""
    regrid_2d_field(
        src_ds::GeoDataset,
        src_var::String,
        target_lon::AbstractVector,
        target_lat::AbstractVector
    ) -> Matrix{Float64}

Regrid a 2D field from a source GeoDataset to target coordinates using bilinear interpolation.
"""
function regrid_2d_field(
    src_ds::GeoDataset,
    src_var::String,
    target_lon::AbstractVector,
    target_lat::AbstractVector
)
    haskey(src_ds.variables, src_var) || error("Variable $src_var not found in dataset")
    
    src_ga = src_ds.variables[src_var]
    src_data = src_ga.data
    
    # Get source coordinates
    src_lons = haskey(src_ds.coords, :lon) ? vec(src_ds.coords[:lon].data) :
               haskey(src_ds.coords, :longitude) ? vec(src_ds.coords[:longitude].data) : error("No source lon")
    src_lats = haskey(src_ds.coords, :lat) ? vec(src_ds.coords[:lat].data) :
               haskey(src_ds.coords, :latitude) ? vec(src_ds.coords[:latitude].data) : error("No source lat")
    
    # Ensure (lon, lat) ordering of data
    if size(src_data) == (length(src_lats), length(src_lons))
        src_data = permutedims(src_data, (2, 1))
    end
    
    # Ensure coordinates are sorted
    if !issorted(src_lons)
        p = sortperm(src_lons); src_lons = src_lons[p]; src_data = src_data[p, :]
    end
    if !issorted(src_lats)
        p = sortperm(src_lats); src_lats = src_lats[p]; src_data = src_data[:, p]
    end
    
    # Bilinear interpolation
    itp = interpolate((src_lons, src_lats), src_data, Gridded(Linear()))
    itp_flat = extrapolate(itp, Flat())
    
    n_tgt_x = length(target_lon)
    n_tgt_y = length(target_lat)
    interpolated = Matrix{Float64}(undef, n_tgt_x, n_tgt_y)
    
    for (i_idx, x_val) in enumerate(target_lon)
        for (j_idx, y_val) in enumerate(target_lat)
            interpolated[i_idx, j_idx] = itp_flat(x_val, y_val)
        end
    end
    
    return interpolated
end

"""
    regrid_3d_field(
        src_ds::GeoDataset,
        src_var::String,
        target_lon::AbstractVector,
        target_lat::AbstractVector,
        target_depth::AbstractVector
    ) -> Array{Float64, 3}

Regrid a 3D field (lon, lat, depth) to target coordinates.
"""
function regrid_3d_field(
    src_ds::GeoDataset,
    src_var::String,
    target_lon::AbstractVector,
    target_lat::AbstractVector,
    target_depth::AbstractVector
)
    haskey(src_ds.variables, src_var) || error("Variable $src_var not found in dataset")
    
    src_ga = src_ds.variables[src_var]
    src_data = src_ga.data
    
    # Get source coordinates
    src_lons = haskey(src_ds.coords, :lon) ? vec(src_ds.coords[:lon].data) :
               haskey(src_ds.coords, :longitude) ? vec(src_ds.coords[:longitude].data) : error("No source lon")
    src_lats = haskey(src_ds.coords, :lat) ? vec(src_ds.coords[:lat].data) :
               haskey(src_ds.coords, :latitude) ? vec(src_ds.coords[:latitude].data) : error("No source lat")
    src_deps = haskey(src_ds.coords, :depth) ? vec(src_ds.coords[:depth].data) :
               haskey(src_ds.coords, :lev) ? vec(src_ds.coords[:lev].data) : error("No source depth")
    
    # Ensure coordinate ordering matches data
    if ndims(src_data) == 3
        if size(src_data, 1) == length(src_lats) && size(src_data, 2) == length(src_lons)
            src_data = permutedims(src_data, (2, 1, 3))
        end
    end
    
    # Ensure coordinates are sorted
    if !issorted(src_lons)
        p = sortperm(src_lons); src_lons = src_lons[p]; src_data = src_data[p, :, :]
    end
    if !issorted(src_lats)
        p = sortperm(src_lats); src_lats = src_lats[p]; src_data = src_data[:, p, :]
    end
    if !issorted(src_deps)
        p = sortperm(src_deps); src_deps = src_deps[p]; src_data = src_data[:, :, p]
    end
    
    # Trilinear interpolation
    itp = interpolate((src_lons, src_lats, src_deps), src_data, Gridded(Linear()))
    itp_flat = extrapolate(itp, Flat())
    
    n_tgt_x = length(target_lon)
    n_tgt_y = length(target_lat)
    n_tgt_z = length(target_depth)
    interpolated = Array{Float64}(undef, n_tgt_x, n_tgt_y, n_tgt_z)
    
    for (i_idx, x_val) in enumerate(target_lon)
        for (j_idx, y_val) in enumerate(target_lat)
            for (k_idx, z_val) in enumerate(target_depth)
                interpolated[i_idx, j_idx, k_idx] = itp_flat(x_val, y_val, z_val)
            end
        end
    end
    
    return interpolated
end

"""
    slice_bathymetry_geodata(ds::GeoDataset; lon_range, lat_range) -> GeoDataset

Slice a bathymetry GeoDataset to a geographic bounding box.
"""
function slice_bathymetry_geodata(ds::GeoDataset; 
    lon_range::Tuple{Real, Real} = (-180.0, 180.0),
    lat_range::Tuple{Real, Real} = (-90.0, 90.0))
    return geoslice(ds; lon=lon_range, lat=lat_range)
end

"""
    slice_wind_geodata(ds::GeoDataset; lon_range, lat_range, time_range=nothing) -> GeoDataset

Slice a wind GeoDataset spatially and optionally temporally.
"""
function slice_wind_geodata(ds::GeoDataset;
    lon_range::Tuple{Real, Real} = (-180.0, 180.0),
    lat_range::Tuple{Real, Real} = (-90.0, 90.0),
    time_range::Union{Nothing, Tuple{Real, Real}} = nothing)
    if isnothing(time_range)
        return geoslice(ds; lon=lon_range, lat=lat_range)
    else
        return geoslice(ds; lon=lon_range, lat=lat_range, time=time_range)
    end
end

"""
    extract_grid_coordinates_geodata(ds::GeoDataset) -> NamedTuple

Extract 1D coordinate vectors from a GeoDataset.
"""
function extract_grid_coordinates_geodata(ds::GeoDataset)
    lons = haskey(ds.coords, :lon) ? vec(ds.coords[:lon].data) :
           haskey(ds.coords, :longitude) ? vec(ds.coords[:longitude].data) : Float64[]
    lats = haskey(ds.coords, :lat) ? vec(ds.coords[:lat].data) :
           haskey(ds.coords, :latitude) ? vec(ds.coords[:latitude].data) : Float64[]
    depths = haskey(ds.coords, :depth) ? vec(ds.coords[:depth].data) :
             haskey(ds.coords, :lev) ? vec(ds.coords[:lev].data) : Float64[]
    
    return (lons = lons, lats = lats, depths = depths)
end

"""
    regrid_2d_field(
        src_lon::AbstractVector,
        src_lat::AbstractVector,
        src_data::AbstractMatrix,
        target_lon::AbstractVector,
        target_lat::AbstractVector
    ) -> Matrix{Float64}

Regrid a 2D field from source coordinates to target coordinates using bilinear interpolation.

This is a standalone version that works with raw arrays rather than GeoDataset objects.
"""
function regrid_2d_field(
    src_lon::AbstractVector,
    src_lat::AbstractVector,
    src_data::AbstractMatrix,
    target_lon::AbstractVector,
    target_lat::AbstractVector
)
    # Ensure (lon, lat) ordering of data
    if size(src_data) == (length(src_lat), length(src_lon))
        src_data = permutedims(src_data, (2, 1))
    end
    
    # Ensure coordinates are sorted
    if !issorted(src_lon)
        p = sortperm(src_lon); src_lon = src_lon[p]; src_data = src_data[p, :]
    end
    if !issorted(src_lat)
        p = sortperm(src_lat); src_lat = src_lat[p]; src_data = src_data[:, p]
    end
    
    # Bilinear interpolation
    itp = interpolate((src_lon, src_lat), src_data, Gridded(Linear()))
    itp_flat = extrapolate(itp, Flat())
    
    n_tgt_x = length(target_lon)
    n_tgt_y = length(target_lat)
    interpolated = Matrix{Float64}(undef, n_tgt_x, n_tgt_y)
    
    for (i_idx, x_val) in enumerate(target_lon)
        for (j_idx, y_val) in enumerate(target_lat)
            interpolated[i_idx, j_idx] = itp_flat(x_val, y_val)
        end
    end
    
    return interpolated
end

"""
    buffer_distance_to_degrees(
        buffer_km::Real,
        ref_lat::Real = 44.5
    ) -> Tuple{Float64, Float64}

Convert a linear physical buffer distance in kilometers \$d_{\\text{buf}}\$ into
equivalent geographic longitude and latitude degree increments \$(\\Delta\\lambda, \\Delta\\phi)\$
at a specified reference latitude \$\\phi_0\$.

# Mathematical Formulation
Using the spherical Earth model with mean radius \$R_{\\text{earth}} = 6371.0088\\text{ km}\$:
```math
\\Delta\\phi = \\frac{d_{\\text{buf}}}{R_{\\text{earth}}} \\times \\left(\\frac{180^\\circ}{\\pi}\\right)
```
```math
\\Delta\\lambda = \\frac{d_{\\text{buf}}}{R_{\\text{earth}} \\cos(\\deg2rad(\\phi_0))} \\times \\left(\\frac{180^\\circ}{\\pi}\\right)
```

# Inputs
- `buffer_km::Real`: Buffer distance in kilometers (e.g. 100.0 km).
- `ref_lat::Real`: Reference latitude in degrees North (default 44.5°N).

# Outputs
- `Tuple{Float64, Float64}`: `(dlon, dlat)` degree increments.
"""
function buffer_distance_to_degrees(
    buffer_km::Real,
    ref_lat::Real = 44.5
)
    if buffer_km < 0.0
        error("Buffer distance must be non-negative: $(buffer_km) km")
    end

    r_earth_km = 6371.0088
    km_per_deg_lat = (π * r_earth_km) / 180.0 # ~111.195 km/deg
    dlat = Float64(buffer_km) / km_per_deg_lat

    # Bounded cosine scaling with 85° Web Mercator cutoff to prevent division by near-zero at high latitudes
    cos_lat = cos(deg2rad(clamp(Float64(ref_lat), -85.0, 85.0)))
    km_per_deg_lon = km_per_deg_lat * max(cosd(85.0), cos_lat)
    dlon = Float64(buffer_km) / km_per_deg_lon

    return (dlon, dlat)
end

"""
    expand_domain_with_buffer(
        lon_range::Tuple{Real, Real},
        lat_range::Tuple{Real, Real};
        buffer_km::Real = 100.0
    ) -> Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}

Expand a geographic bounding box \$(\\lambda_{\\min}, \\lambda_{\\max}) \\times (\\phi_{\\min}, \\phi_{\\max})\$
outward by a user-defined physical buffer distance (default 100.0 km).

# Inputs
- `lon_range::Tuple{Real, Real}`: Input longitude bounds.
- `lat_range::Tuple{Real, Real}`: Input latitude bounds.
- `buffer_km::Real`: Buffer distance in kilometers (default 100.0 km).

# Outputs
- `Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}`: `(buffered_lon_range, buffered_lat_range)`
"""
function expand_domain_with_buffer(
    lon_range::Tuple{Real, Real},
    lat_range::Tuple{Real, Real};
    buffer_km::Real = 100.0
)
    ref_lat = 0.5 * (Float64(lat_range[1]) + Float64(lat_range[2]))
    dlon, dlat = buffer_distance_to_degrees(buffer_km, ref_lat)

    buf_lon = (
        max(-180.0, Float64(lon_range[1]) - dlon),
        min(180.0, Float64(lon_range[2]) + dlon)
    )
    buf_lat = (
        max(-90.0, Float64(lat_range[1]) - dlat),
        min(90.0, Float64(lat_range[2]) + dlat)
    )

    return (buf_lon, buf_lat)
end