"""
    winds.jl

Surface wind fetch, load, and processing operations via GeoData.

Supports Open-Meteo ERA5 reanalysis, with bulk heat flux computation.
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem, DIM_TEMPORAL
using Downloads
using HTTP
using JSON3
using Statistics
using LinearAlgebra
using Dates

"""
    fetch_open_meteo_winds(;
        lon_range::Tuple{Real, Real},
        lat_range::Tuple{Real, Real},
        time_iso::AbstractString = "2023-06-01T00:00:00Z",
        output_path::AbstractString = joinpath("inputs", "wind.zarr"),
        backend::Symbol = :zarr,
        verbose::Bool = true
    ) -> GeoDataset

Fetch ERA5 surface winds from Open-Meteo and save as GeoDataset.

# Arguments
- `lon_range`: Longitude bounds (min, max)
- `lat_range`: Latitude bounds (min, max)
- `time_iso`: ISO time string (date part used)
- `output_path`: Output file path
- `backend`: Storage backend (:zarr, :netcdf)
- `verbose`: Print progress messages

# Returns
GeoDataset with tau_x, tau_y, and optional bulk flux variables.
"""
function fetch_open_meteo_winds(;
    lon_range::Tuple{Real, Real},
    lat_range::Tuple{Real, Real},
    time_iso::AbstractString = "2023-06-01T00:00:00Z",
    output_path::AbstractString = joinpath("inputs", "wind.zarr"),
    backend::Symbol = :zarr,
    verbose::Bool = true
)
    mkpath(dirname(output_path))
    date_str = split(time_iso, "T")[1]
    clat = 0.5 * (lat_range[1] + lat_range[2])
    clon = 0.5 * (lon_range[1] + lon_range[2])

    # Build Open-Meteo API URL
    url = "https://archive-api.open-meteo.com/v1/archive?" *
          "latitude=$(clat)&longitude=$(clon)&start_date=$(date_str)&" *
          "end_date=$(date_str)&" *
          "hourly=wind_speed_10m,wind_direction_10m,temperature_2m,surface_pressure," *
          "relative_humidity_2m,cloud_cover,shortwave_radiation&" *
          "wind_speed_unit=ms"

    verbose && println("Requesting Open-Meteo ERA5 reanalysis winds for $(date_str)...")
    tmp_json = tempname() * ".json"
    try
        Downloads.download(url, tmp_json)
    catch err
        rm(tmp_json, force=true)
        throw(err)
    end

    json_str = read(tmp_json, String)
    rm(tmp_json, force=true)

    time_m = match(r"\"time\"\s*:\s*\[([^\]]+)\]", json_str)
    spd_m  = match(r"\"wind_speed_10m\"\s*:\s*\[([^\]]+)\]", json_str)
    dir_m  = match(r"\"wind_direction_10m\"\s*:\s*\[([^\]]+)\]", json_str)

    if isnothing(time_m) || isnothing(spd_m) || isnothing(dir_m)
        error("Malformed JSON received from Open-Meteo API.")
    end

    times_raw = [replace(strip(s), "\"" => "") for s in split(time_m.captures[1], ",")]
    speeds    = [parse(Float64, strip(s)) for s in split(spd_m.captures[1], ",")]
    dirs      = [parse(Float64, strip(s)) for s in split(dir_m.captures[1], ",")]
    nt = length(times_raw)

    # Surface state for bulk heat flux
    hourly_of(name) = begin
        m = match(Regex("\\\"$(name)\\\"\\s*:\\s*\\[([^\\]]+)\\]"), json_str)
        isnothing(m) ? nothing : [parse(Float64, strip(s)) for s in split(m.captures[1], ",")]
    end
    air_t   = hourly_of("temperature_2m")
    sfc_p   = hourly_of("surface_pressure")
    rel_hum = hourly_of("relative_humidity_2m")
    cloud   = hourly_of("cloud_cover")
    sw_down = hourly_of("shortwave_radiation")
    have_flux = !(isnothing(air_t) || isnothing(sfc_p) || isnothing(rel_hum) ||
                  isnothing(cloud) || isnothing(sw_down))

    # Build spatial grid
    n_lon, n_lat = 50, 50
    lon_coords = range(lon_range[1], lon_range[2], length = n_lon)
    lat_coords = range(lat_range[1], lat_range[2], length = n_lat)
    time_secs  = collect(range(0.0, step = 3600.0, length = nt))

    # Convert wind to stress
    u10_hourly = [-speeds[t] * sind(dirs[t]) for t in 1:nt]
    v10_hourly = [-speeds[t] * cosd(dirs[t]) for t in 1:nt]

    tau_x_3d = Array{Float64}(undef, n_lon, n_lat, nt)
    tau_y_3d = Array{Float64}(undef, n_lon, n_lat, nt)

    for t in 1:nt
        tx, ty = wind_speed_to_kinematic_stress(u10_hourly[t], v10_hourly[t])
        for i in 1:n_lon, j in 1:n_lat
            tau_x_3d[i, j, t] = tx
            tau_y_3d[i, j, t] = ty
        end
    end

    # Create GeoDataset
    lon_dim = Dimension(name=:lon, size=n_lon, coords=collect(lon_coords), units="degrees_east")
    lat_dim = Dimension(name=:lat, size=n_lat, coords=collect(lat_coords), units="degrees_north")
    time_dim = Dimension(name=:time, size=nt, coords=time_secs, units="seconds since $(date_str)T00:00:00Z", dim_type=DIM_TEMPORAL)
    
    crs = CoordinateSystem(crs="EPSG:4326")
    
    tau_x_ga = GeoArray(tau_x_3d, (lon_dim, lat_dim, time_dim), crs, Dict("units" => "m2 s-2", "standard_name" => "surface_downward_x_stress"))
    tau_y_ga = GeoArray(tau_y_3d, (lon_dim, lat_dim, time_dim), crs, Dict("units" => "m2 s-2", "standard_name" => "surface_downward_y_stress"))
    
    coords = Dict(
        :lon => GeoArray(collect(lon_coords), (lon_dim,), crs, Dict("units" => "degrees_east")),
        :lat => GeoArray(collect(lat_coords), (lat_dim,), crs, Dict("units" => "degrees_north")),
        :time => GeoArray(time_secs, (time_dim,), crs, Dict("units" => "seconds since $(date_str)T00:00:00Z"))
    )
    dims = Dict(:lon => lon_dim, :lat => lat_dim, :time => time_dim)
    variables = Dict("tau_x" => tau_x_ga, "tau_y" => tau_y_ga)

    if have_flux
        for (nm, vals, scale, offset, units, sname) in (
                ("air_temperature", air_t, 1.0, 273.15, "K", "air_temperature"),
                ("surface_pressure", sfc_p, 100.0, 0.0, "Pa", "surface_air_pressure"),
                ("relative_humidity", rel_hum, 1.0, 0.0, "percent", "relative_humidity"),
                ("cloud_fraction", cloud, 1.0, 0.0, "1", "cloud_area_fraction"),
                ("sw_down", sw_down, 1.0, 0.0, "W m-2", "surface_downwelling_shortwave_flux"),
            )
            f3 = Array{Float64}(undef, n_lon, n_lat, nt)
            for t in 1:nt, i in 1:n_lon, j in 1:n_lat
                f3[i, j, t] = vals[t] * scale + offset
            end
            ga = GeoArray(f3, (lon_dim, lat_dim, time_dim), crs, Dict("units" => units, "standard_name" => sname))
            variables[nm] = ga
        end
    end

    ds = GeoDataset(variables, coords, dims, crs, 
        Dict("title" => "Open-Meteo ERA5 Reanalysis Surface Wind Forcing",
             "source" => "Open-Meteo Historical Weather API (ERA5/ECMWF)"),
        nothing, output_path)

    geosave(output_path, ds; backend=backend)
    verbose && println("Successfully saved winds to: $(output_path)")
    return ds
end

"""
    load_wind_stress_geodata(filepath::AbstractString) -> NamedTuple

Read wind stress from a GeoData wind file.
Returns (tau_x, tau_y, tau_max, speed10_rms) where stresses are domain means.
"""
function load_wind_stress_geodata(filepath::AbstractString)
    ds = geoload(filepath)
    haskey(ds.variables, "tau_x") || error("Wind file has no tau_x variable")
    
    tx = Float64.(ds.variables["tau_x"].data)
    ty = Float64.(ds.variables["tau_y"].data)
    finite = isfinite.(tx) .& isfinite.(ty)
    any(finite) || error("Wind file contains no finite stress values")
    
    tau_x = sum(tx[finite]) / count(finite)
    tau_y = sum(ty[finite]) / count(finite)
    tau_max = maximum(hypot.(tx, ty)[finite])
    speed10_rms = wind_speed_from_stress(hypot(tau_x, tau_y))
    
    return (tau_x = tau_x, tau_y = tau_y, tau_max = tau_max, speed10_rms = speed10_rms)
end

"""
    wind_speed_to_kinematic_stress(u10, v10; ρ_air=1.225, ρ_water=1025.0)

Convert 10m wind velocity to kinematic surface wind stress.
"""
function wind_speed_to_kinematic_stress(
    u10::Real,
    v10::Real;
    ρ_air::Real = 1.225,
    ρ_water::Real = 1025.0
)
    speed = sqrt(u10^2 + v10^2)
    if speed == 0.0
        return (0.0, 0.0)
    end

    cd = if speed <= 11.0
        1.2e-3
    else
        (0.49 + 0.065 * speed) * 1e-3
    end

    factor = (ρ_air / ρ_water) * cd * speed
    tau_x = factor * u10
    tau_y = factor * v10

    return (Float64(tau_x), Float64(tau_y))
end

"""
    wind_speed_from_stress(tau; ρ_air=1.225, ρ_water=1025.0)

Recover 10m wind speed from kinematic surface stress (inverse of wind_speed_to_kinematic_stress).
"""
function wind_speed_from_stress(tau::Real; ρ_air::Real = 1.225, ρ_water::Real = 1025.0)
    t = abs(Float64(tau))
    t <= 0 && return 0.0
    lo, hi = 0.0, 120.0
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        s = wind_speed_to_kinematic_stress(mid, 0.0; ρ_air = ρ_air, ρ_water = ρ_water)[1]
        s < t ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

"""
    build_bulk_surface_flux_geodata(wind_file::AbstractString; albedo=0.06) -> Tuple

Read a surface-wind GeoDataset and return (tau_x_array, tau_y_array, heat_flux_function).

Returns raw arrays and a heat_flux function that computes flux at (x, y, t, T_surf).
The heat_flux function signature: (x, y, t, T_surf) -> flux [W/m²]
"""
function build_bulk_surface_flux_geodata(wind_file::AbstractString; albedo::Real = 0.06)
    ds = geoload(wind_file)
    needed = ("tau_x", "tau_y", "air_temperature", "surface_pressure",
              "relative_humidity", "cloud_fraction", "sw_down", "time")
    missing_vars = [v for v in needed if !haskey(ds.variables, v)]
    if !isempty(missing_vars)
        @warn "Wind file lacks $(join(missing_vars, ", ")); no bulk heat flux can be computed."
        return (nothing, nothing, nothing)
    end

    tx = Float64.(ds.variables["tau_x"].data)
    ty = Float64.(ds.variables["tau_y"].data)
    ta = Float64.(ds.variables["air_temperature"].data)
    ps = Float64.(ds.variables["surface_pressure"].data)
    rh = Float64.(ds.variables["relative_humidity"].data)
    cc = Float64.(ds.variables["cloud_fraction"].data)
    sw = Float64.(ds.variables["sw_down"].data)
    tvec = Float64.(ds.coords[:time].data)

    nt = size(tx, 3)
    at(k) = Float64(ta[1, 1, k])
    sp(k) = Float64(ps[1, 1, k])
    rhu(k) = Float64(rh[1, 1, k]) / 100
    swd(k) = Float64(sw[1, 1, k])

    function cf(k)
        c = Float64(cc[1, 1, k]) / 100.0
        if c < 0.0 || c > 1.0
            error("Cloud fraction $(c) at time index $(k) out of physical range [0, 1].")
        end
        return c
    end

    tstep = nt > 1 ? (tvec[end] - tvec[1]) / (nt - 1) : 3600.0
    horizon = nt > 1 ? tvec[end] - tvec[1] : 3600.0

    function kof(t)
        if horizon <= 0.0 || nt <= 1
            return 1
        end
        k = Int(floor(mod(Float64(t), horizon) / tstep)) + 1
        return k > nt ? 1 : (k < 1 ? 1 : k)
    end

    sigma = 5.670374419e-8
    rho_a, cp_a, Lv = 1.225, 1005.0, 2.501e6
    Cd, Ce = 1.3e-3, 1.5e-3

    esat(T) = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
    qsat(p, T) = 0.622 * esat(T) / (p - 0.378 * esat(T))

    function heat_flux(x, y, t, T_surf)
        k = kof(t)
        Ta = at(k)
        e_a = rhu(k) * esat(Ta)
        U = wind_speed_from_stress(hypot(Float64(tx[1, 1, k]), Float64(ty[1, 1, k])))
        sw_net = (1 - Float64(albedo)) * (1 - 0.65 * cf(k)^2) * swd(k)
        lw_net = sigma * T_surf^4 * (0.34 - 0.14 * e_a / 1000) * (1 - 0.8 * cf(k)) -
                 sigma * Ta^4 * (0.34 - 0.14 * e_a / 1000) * (1 - 0.8 * cf(k))
        sens = rho_a * cp_a * Cd * U * (T_surf - Ta)
        qa = 0.622 * e_a / (sp(k) - 0.378 * e_a)
        lat = rho_a * Lv * Ce * U * (qsat(sp(k), T_surf) - qa)
        return sw_net - lw_net - sens - lat
    end

    return (tx, ty, heat_flux)
end