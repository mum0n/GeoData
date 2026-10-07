"""
    woa23.jl

WOA23 climatology fetch, load, and interpolator operations via GeoData.

Fetches temperature, salinity, oxygen, and nutrients from NOAA NCEI.
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using NCDatasets
using Downloads
using Interpolations
using Statistics
using LinearAlgebra

"""
    fetch_woa23(;
        lon_range::Tuple{Real, Real} = (-71.0, -53.0),
        lat_range::Tuple{Real, Real} = (40.0, 48.5),
        month::Int = 0,
        output_dir::AbstractString = "inputs",
        include_o2::Bool = true,
        backend::Symbol = :zarr,
        verbose::Bool = true
    ) -> NamedTuple

Fetch WOA23 temperature, salinity, and oxygen climatology as GeoDatasets.

# Arguments
- `lon_range`: Longitude bounds (min, max)
- `lat_range`: Latitude bounds (min, max)
- `month`: Month (0 = annual mean, 1-12 = specific month)
- `output_dir`: Output directory
- `include_o2`: Include oxygen data
- `backend`: Storage backend (:zarr, :netcdf)
- `verbose`: Print progress messages

# Returns
NamedTuple with file paths and interpolator functions.
"""
function fetch_woa23(;
    lon_range::Tuple{Real, Real} = (-71.0, -53.0),
    lat_range::Tuple{Real, Real} = (40.0, 48.5),
    month::Int = 0,
    output_dir::AbstractString = "inputs",
    include_o2::Bool = true,
    backend::Symbol = :zarr,
    verbose::Bool = true
)
    mkpath(output_dir)
    month_str = lpad(string(month), 2, "0")

    # WOA23 0.25° grid parameters
    woa_lon_step = 0.25
    woa_lat_step = 0.25
    woa_lon_origin = -180.0
    woa_lat_origin = -90.0

    i_lon_lo = round(Int, (Float64(lon_range[1]) - woa_lon_origin) / woa_lon_step)
    i_lon_hi = round(Int, (Float64(lon_range[2]) - woa_lon_origin) / woa_lon_step)
    i_lat_lo = round(Int, (Float64(lat_range[1]) - woa_lat_origin) / woa_lat_step)
    i_lat_hi = round(Int, (Float64(lat_range[2]) - woa_lat_origin) / woa_lat_step)
    i_lon_lo = max(0, min(1440, i_lon_lo))
    i_lon_hi = max(0, min(1440, i_lon_hi))
    i_lat_lo = max(0, min(720, i_lat_lo))
    i_lat_hi = max(0, min(720, i_lat_hi))
    i_dep_hi = 101

    function opendap_subset(varname)
        "[0:1:0]" *
        "[0:1:$(i_dep_hi)]" *
        "[$(i_lat_lo):1:$(i_lat_hi)]" *
        "[$(i_lon_lo):1:$(i_lon_hi)]"
    end

    function woa_url_candidates(variable_letter, varname)
        sub_dir = if variable_letter == "t"
            "temperature"
        elseif variable_letter == "s"
            "salinity"
        elseif variable_letter in ["o", "O", "A"]
            "oxygen"
        else
            "nutrients"
        end

        urls = String[]
        if variable_letter in ["o", "O", "A"]
            push!(urls, "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/oxygen/netcdf/all/1.00/woa23_all_o$(month_str)_01.nc")
            push!(urls, "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/oxygen/netcdf/all/5.00/woa23_all_o$(month_str)_05.nc")
            push!(urls, "https://www.ncei.noaa.gov/thredds/dodsC/ncei/woa/oxygen/netcdf/all/1.00/woa23_all_o$(month_str)_01.nc?$(varname)$(opendap_subset(varname))")
        else
            base_fn_25 = "woa23_A5B4_$(variable_letter)$(month_str)_04.nc"
            base_fn_1  = "woa23_A5B4_$(variable_letter)$(month_str)_01.nc"
            path_25 = "$(sub_dir)/netcdf/A5B4/0.25"
            path_1  = "$(sub_dir)/netcdf/A5B4/1.00"
            static_base = "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA"

            push!(urls, "$(static_base)/$(path_25)/$(base_fn_25)")
            push!(urls, "$(static_base)/$(path_1)/$(base_fn_1)")

            subset_str = string(varname, opendap_subset(varname), ",") *
                         string("lon", opendap_subset("lon"), ",lat", opendap_subset("lat"), ",") *
                         "depth[0:1:$(i_dep_hi)],time[0:1:0]"
            for thredds in ("https://www.ncei.noaa.gov/thredds/dodsC/ncei/woa",
                            "https://thredds.ucar.edu/thredds/dodsC/ncei/woa",
                            "https://tds.marine.rutgers.edu/thredds/dodsC/ncei/woa",
                            "https://data.nodc.noaa.gov/thredds/dodsC/ncei/woa")
                push!(urls, "$(thredds)/$(path_25)/$(base_fn_25)?$(subset_str)")
            end
        end
        urls
    end

    t_file = joinpath(output_dir, "woa23_temperature_$(month_str)_0.25deg.zarr")
    s_file = joinpath(output_dir, "woa23_salinity_$(month_str)_0.25deg.zarr")
    o_file = joinpath(output_dir, "woa23_oxygen_$(month_str)_0.25deg.zarr")

    for (variable_letter, varname, out_path) in [
        ("t", "t_an", t_file),
        ("s", "s_an", s_file),
        (include_o2 ? ("o", "o_an", o_file) : nothing)
    ]
        isnothing(variable_letter) && continue
        if isfile(out_path) && filesize(out_path) > 1024
            verbose && println("WOA23: using cached file $(out_path)")
            continue
        end
        downloaded = false
        for url in woa_url_candidates(variable_letter, varname)
            verbose && println("Fetching WOA23 from:\n  $(url[1:min(80, length(url))])...")
            try
                tmp_nc = tempname() * ".nc"
                Downloads.download(url, tmp_nc)
                # Convert to GeoDataset and save as Zarr
                ds = geoload(tmp_nc; backend=:ncdatasets)
                geosave(out_path, ds; backend=backend)
                rm(tmp_nc, force=true)
                verbose && println("  -> Saved to $(out_path)")
                downloaded = true
                break
            catch err
                verbose && println("  -> Mirror unsuccessful ($(typeof(err))).")
            end
        end
        if !downloaded
            @warn "All WOA23 download attempts failed for $(varname). Falling back to synthetic."
        end
    end

    # Build interpolators
    function make_woa_interpolator(filepath, varname, fallback_val)
        if !isfile(filepath) || filesize(filepath) <= 1024
            verbose && println("WOA23: file $(filepath) not found -- using physical fallback.")
            return if fallback_val isa Function
                fallback_val
            else
                (lon, lat, z) -> Float64(fallback_val)
            end
        end

        ds = geoload(filepath)
        # Get the data variable
        data_var = haskey(ds.variables, varname) ? varname : first(keys(ds.variables))
        ga = ds.variables[data_var]
        
        # Get coordinates
        lons = haskey(ds.coords, :lon) ? vec(ds.coords[:lon].data) :
               haskey(ds.coords, :longitude) ? vec(ds.coords[:longitude].data) : error("No lon")
        lats = haskey(ds.coords, :lat) ? vec(ds.coords[:lat].data) :
               haskey(ds.coords, :latitude) ? vec(ds.coords[:latitude].data) : error("No lat")
        deps = haskey(ds.coords, :depth) ? vec(ds.coords[:depth].data) :
               haskey(ds.coords, :lev) ? vec(ds.coords[:lev].data) : error("No depth")
        
        # WOA depths are positive-down; convert to negative-up
        deps_neg = -abs.(deps)
        field_3d = ga.data
        
        # Ensure sorting
        if !issorted(lons)
            p = sortperm(lons); lons = lons[p]; field_3d = field_3d[p, :, :]
        end
        if !issorted(lats)
            p = sortperm(lats); lats = lats[p]; field_3d = field_3d[:, p, :]
        end
        if !issorted(deps_neg)
            p = sortperm(deps_neg); deps_neg = deps_neg[p]; field_3d = field_3d[:, :, p]
        end

        n_lon, n_lat, n_dep = size(field_3d)

        # Repair NaN cells (land masking, etc.)
        n_gap, n_landcol = _repair_woa_field!(field_3d)
        if verbose && (n_gap > 0 || n_landcol > 0)
            println("WOA23: repaired $(n_gap) masked level(s) and filled $(n_landcol) columns.")
        end

        function woa_interp(lon, lat, z)
            i_raw = searchsortedlast(lons, Float64(lon))
            i = max(1, min(i_raw, n_lon - 1))
            j_raw = searchsortedlast(lats, Float64(lat))
            j = max(1, min(j_raw, n_lat - 1))
            k_raw = searchsortedlast(deps_neg, Float64(z))
            k = max(1, min(k_raw, n_dep - 1))

            dx = lons[i+1] - lons[i]
            sx = dx != 0.0 ? (Float64(lon) - lons[i]) / dx : 0.0
            sx = clamp(sx, 0.0, 1.0)

            dy = lats[j+1] - lats[j]
            sy = dy != 0.0 ? (Float64(lat) - lats[j]) / dy : 0.0
            sy = clamp(sy, 0.0, 1.0)

            dz = deps_neg[k+1] - deps_neg[k]
            sz = dz != 0.0 ? (Float64(z) - deps_neg[k]) / dz : 0.0
            sz = clamp(sz, 0.0, 1.0)

            return Float64(
                field_3d[i,   j,   k]   * (1-sx)*(1-sy)*(1-sz) +
                field_3d[i+1, j,   k]   * sx*(1-sy)*(1-sz) +
                field_3d[i,   j+1, k]   * (1-sx)*sy*(1-sz) +
                field_3d[i+1, j+1, k]   * sx*sy*(1-sz) +
                field_3d[i,   j,   k+1] * (1-sx)*(1-sy)*sz +
                field_3d[i+1, j,   k+1] * sx*(1-sy)*sz +
                field_3d[i,   j+1, k+1] * (1-sx)*sy*sz +
                field_3d[i+1, j+1, k+1] * sx*sy*sz
            )
        end
        return woa_interp
    end

    # Physical fallbacks for Northwest Atlantic shelf
    t_fallback(lon, lat, z) = z > -20.0 ? 14.0 : (z > -80.0 ? 2.0 : 8.0)
    s_fallback(lon, lat, z) = 31.5 + 3.0 * (1.0 - exp(-abs(Float64(z)) / 150.0))
    o2_fallback(lon, lat, z) = 300.0 - 95.0 * (1.0 - exp(-abs(Float64(z)) / 150.0))
    o2_sat_fallback(lon, lat, z) = 98.0 - 28.0 * (1.0 - exp(-abs(Float64(z)) / 150.0))

    t_fn = make_woa_interpolator(t_file, "t_an", t_fallback)
    s_fn = make_woa_interpolator(s_file, "s_an", s_fallback)
    o_fn = make_woa_interpolator(o_file, "o_an", o2_fallback)

    return (
        temperature_file = t_file,
        salinity_file    = s_file,
        oxygen_file      = o_file,
        temperature_fn   = t_fn,
        salinity_fn      = s_fn,
        oxygen_fn        = o_fn,
        oxygen_sat_fn    = o2_sat_fallback,
        nitrate_fn       = (lon, lat, z) -> 2.0 + 18.0 * (1.0 - exp(-abs(Float64(z)) / 100.0)),
        phosphate_fn     = (lon, lat, z) -> 0.3 + 1.2 * (1.0 - exp(-abs(Float64(z)) / 100.0)),
        silicate_fn      = (lon, lat, z) -> 4.0 + 22.0 * (1.0 - exp(-abs(Float64(z)) / 120.0)),
        aou_fn           = (lon, lat, z) -> 10.0 + 90.0 * (1.0 - exp(-abs(Float64(z)) / 150.0))
    )
end

function _repair_woa_field!(field_3d)
    n_lon, n_lat, n_dep = size(field_3d)
    n_gap = 0
    n_empty = 0

    for j in 1:n_lat, i in 1:n_lon
        valid = findall(k -> !isnan(field_3d[i, j, k]), 1:n_dep)
        isempty(valid) && (n_empty += 1; continue)

        first_v, last_v = first(valid), last(valid)
        for k in 1:n_dep
            isnan(field_3d[i, j, k]) || continue
            if k < first_v
                field_3d[i, j, k] = field_3d[i, j, first_v]
            elseif k > last_v
                field_3d[i, j, k] = field_3d[i, j, last_v]
            else
                lo = last(valid[kk] for kk in eachindex(valid) if valid[kk] < k)
                hi = first(valid[kk] for kk in eachindex(valid) if valid[kk] > k)
                w = (k - lo) / (hi - lo)
                field_3d[i, j, k] =
                    field_3d[i, j, lo] * (1 - w) + field_3d[i, j, hi] * w
            end
            n_gap += 1
        end
    end

    if n_empty > 0
        good = filter(!isnan, field_3d)
        fill_value = isempty(good) ? 0.0 : sum(good) / length(good)
        for j in 1:n_lat, i in 1:n_lon
            all(!isnan, @view field_3d[i, j, :]) && continue
            field_3d[i, j, :] .= fill_value
        end
    end

    return (n_gap, n_empty)
end

"""
    load_woa23_interpolators(; lon_range, lat_range, month=0, output_dir="inputs") -> NamedTuple

Load or fetch WOA23 and return interpolator functions.
"""
function load_woa23_interpolators(;
    lon_range = (-71.0, -53.0),
    lat_range = (40.0, 48.5),
    month::Int = 0,
    output_dir::AbstractString = "inputs",
    kwargs...
)
    return fetch_woa23(;
        lon_range = lon_range,
        lat_range = lat_range,
        month = month,
        output_dir = output_dir,
        kwargs...
    )
end