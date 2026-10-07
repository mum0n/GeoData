"""
    boundary.jl

Boundary hydrography fetch and interpolator construction via GeoData.

Supports WOA23, GLORYS12v1 (Copernicus), HYCOM for sponge boundary conditions.
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using NCDatasets
using Interpolations
using Statistics

"""
    fetch_boundary_hydrography_geodata(source::Symbol;
        lon_range, lat_range, input_dir,
        month = 0, verbose = true) -> NamedTuple

Acquire and load boundary temperature and salinity as sponge targets using GeoData.

# Arguments
- `source`: Data source (:woa23, :glorys12v1, :glorys_climatology, :hycom)
- `lon_range`: Longitude bounds (min, max)
- `lat_range`: Latitude bounds (min, max)
- `input_dir`: Input/output directory
- `month`: Month (0 = annual mean)
- `verbose`: Print progress messages

# Returns
NamedTuple with interpolator functions (T, S).
"""
function fetch_boundary_hydrography_geodata(source::Symbol;
    lon_range, lat_range, input_dir,
    month::Int = 0, verbose::Bool = true)
    mkpath(input_dir)
    lon = (Float64(lon_range[1]), Float64(lon_range[2]))
    lat = (Float64(lat_range[1]), Float64(lat_range[2]))

    if source === :woa23
        mm = lpad(string(month), 2, '0')
        tp = joinpath(input_dir, "woa23_temperature_$(mm)_0.25deg.zarr")
        sp = joinpath(input_dir, "woa23_salinity_$(mm)_0.25deg.zarr")
        if !isfile(tp) || !isfile(sp)
            verbose && println("Fetching WOA23 climatology (keyless, NOAA)...")
            fetch_woa23(lon_range = lon_range, lat_range = lat_range,
                        month = month, output_dir = input_dir, backend=:zarr)
        end
        (isfile(tp) && isfile(sp)) || error(
            "WOA23 boundary needs $(basename(tp)) and $(basename(sp)) in $(input_dir)")
        return build_boundary_tracer_interpolators_geodata(tp, sp; T_name = "t_an", S_name = "s_an")
    end

    if source in (:glorys12v1, :glorys_climatology)
        p = joinpath(input_dir, "boundary_glorys_1993_$(lon[1])_$(lon[2])_$(lat[1])_$(lat[2]).zarr")
        if !isfile(p)
            verbose && println("Fetching GLORYS12V1 boundary subset (Copernicus Marine)...")
            # This would need the copernicus fetch to be adapted to GeoData
            Copernicus.fetch_copernicus_physics_subset(
                lon_range = lon_range, lat_range = lat_range,
                start_date = "1993-01-01", end_date = "1993-12-31",
                output_path = replace(p, ".zarr" => ".nc"),
                dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1M-m"
            )
            # Convert to zarr
            ds = geoload(replace(p, ".zarr" => ".nc"); backend=:ncdatasets)
            geosave(p, ds; backend=:zarr)
        end
        m = month == 0 ? collect(1:12) : [max(1, min(12, month))]
        return build_boundary_tracer_interpolators_geodata(p, p; T_name = "thetao", S_name = "so",
                                                            months = m)
    end

    if source === :hycom
        p = joinpath(input_dir, "boundary_hycom.zarr")
        if !isfile(p)
            error(
                "HYCOM GOFS OPeNDAP server (tds.hycom.org) is currently unreachable or returns " *
                "HTTP 400. Use source = :woa23 or source = :glorys12v1 instead, or place a " *
                "pre-downloaded 'boundary_hycom.zarr' in '$(input_dir)'."
            )
        end
        return build_boundary_tracer_interpolators_geodata(p, p; T_name = "water_temp",
                                                            S_name = "salinity", months = [1])
    end

    error(
        "Unknown boundary source \"$(source)\". Choose :woa23, :glorys12v1, :hycom, or :synthetic."
    )
end

"""
    build_boundary_tracer_interpolators_geodata(temp_file, salt_file; T_name, S_name, months=[1]) -> NamedTuple

Build (lon, lat, z, t) interpolators for boundary T and S from GeoData files.

# Arguments
- `temp_file`: Path to temperature GeoDataset
- `salt_file`: Path to salinity GeoDataset
- `T_name`: Temperature variable name
- `S_name`: Salinity variable name
- `months`: Months to average (if time dimension present)

# Returns
NamedTuple with (T, S) interpolator functions.
"""
function build_boundary_tracer_interpolators_geodata(
    temp_file::AbstractString, salt_file::AbstractString;
    T_name::AbstractString = "t_an", S_name::AbstractString = "s_an",
    months::Vector{Int} = [1]
)
    # Load temperature GeoDataset
    temp_ds = geoload(temp_file)
    salt_ds = geoload(salt_file)

    # Get coordinate arrays
    lons = haskey(temp_ds.coords, :lon) ? vec(temp_ds.coords[:lon].data) :
           haskey(temp_ds.coords, :longitude) ? vec(temp_ds.coords[:longitude].data) : error("No lon in temp")
    lats = haskey(temp_ds.coords, :lat) ? vec(temp_ds.coords[:lat].data) :
           haskey(temp_ds.coords, :latitude) ? vec(temp_ds.coords[:latitude].data) : error("No lat in temp")
    deps = haskey(temp_ds.coords, :depth) ? vec(temp_ds.coords[:depth].data) :
           haskey(temp_ds.coords, :lev) ? vec(temp_ds.coords[:lev].data) : error("No depth in temp")
    times = haskey(temp_ds.coords, :time) ? vec(temp_ds.coords[:time].data) : Float64[0.0]

    # Get data arrays
    T_data = temp_ds.variables[T_name].data
    S_data = salt_ds.variables[S_name].data

    # Handle time dimension - average over specified months if multiple
    if ndims(T_data) == 4 && length(months) > 1
        T_data = mean(T_data[:, :, :, months], dims=4)[:, :, :, 1]
        S_data = mean(S_data[:, :, :, months], dims=4)[:, :, :, 1]
    elseif ndims(T_data) == 4
        T_data = T_data[:, :, :, months[1]]
        S_data = S_data[:, :, :, months[1]]
    end

    # Ensure (lon, lat, depth) ordering
    if size(T_data, 1) == length(lats) && size(T_data, 2) == length(lons)
        T_data = permutedims(T_data, (2, 1, 3))
        S_data = permutedims(S_data, (2, 1, 3))
    end

    # Convert depth to negative-up if positive-down
    if !isempty(deps) && minimum(deps) > -0.5 * maximum(abs, deps)
        deps = -deps
    end

    # Sort coordinates
    if !issorted(lons)
        p = sortperm(lons); lons = lons[p]; T_data = T_data[p, :, :]; S_data = S_data[p, :, :]
    end
    if !issorted(lats)
        p = sortperm(lats); lats = lats[p]; T_data = T_data[:, p, :]; S_data = S_data[:, p, :]
    end
    if !issorted(deps)
        p = sortperm(deps); deps = deps[p]; T_data = T_data[:, :, p]; S_data = S_data[:, :, p]
    end

    # Create 3D interpolators (lon, lat, depth)
    itp_T = interpolate((lons, lats, deps), T_data, Gridded(Linear()))
    itp_T_ext = extrapolate(itp_T, Flat())
    itp_S = interpolate((lons, lats, deps), S_data, Gridded(Linear()))
    itp_S_ext = extrapolate(itp_S, Flat())

    # Return callable functions (lon, lat, z, t) -> value
    # Time dimension is handled by the sponge relaxation, not interpolation
    T_fn = (lon, lat, z, t) -> Float64(itp_T_ext(Float64(lon), Float64(lat), Float64(z)))
    S_fn = (lon, lat, z, t) -> Float64(itp_S_ext(Float64(lon), Float64(lat), Float64(z)))

    return (T = T_fn, S = S_fn)
end