"""
    fetch_baseline_geodata.jl

Baseline geospatial data cache initialization script for GeoData.

Fetches key datasets required by ParticleTracking simulations (matching the
specifications in `configs/snowcrab.toml`) and registers them in the persistent
`GeoDataCatalog` (`~/.geodata/catalog.json` by default).

# Datasets Fetched:
1. Coastline: Natural Earth 10m high-resolution vector coastline
   clipped to the Scotian Shelf / NW Atlantic domain (`coastline.parquet`).
2. Bathymetry: NOAA ETOPO 2022 15 arc-second gridded relief (`bathymetry.zarr` / `.nc`).
3. Surface Winds / Atmosphere: Open-Meteo ERA5 surface winds and radiative fluxes
   for the baseline period (`wind.zarr`).
4. Hydrography: World Ocean Atlas 2023 (WOA23) climatological temperature, salinity,
   and dissolved oxygen (`woa23_*_0.25deg.zarr`).
5. Tidal Solution: TPXO9-Atlas global harmonic velocity solution (`u_tpxo9.v1.nc`).

# Usage:
```sh
julia --project=. scripts/fetch_baseline_geodata.jl [options]
```

# Options:
- `--cache-dir <dir>` : Target directory for downloaded artefacts (default: `~/.geodata/cache`).
- `--catalog <path>`   : Path to GeoData catalog JSON file (default: `~/.geodata/catalog.json`).
- `--skip-tides`       : Skip fetching the large TPXO9 tidal solution archive.
- `--force`            : Force re-downloading even if cached files exist.
- `--dry-run`          : Display planned actions and domain coordinates without downloading.
"""

using Dates
using UUIDs

# Ensure GeoData is loaded from local repository
using GeoData
using GeoData.Data:
    fetch_natural_earth_coastline,
    fetch_erddap_bathymetry,
    fetch_open_meteo_winds,
    fetch_woa23,
    DATA_SOURCES,
    fetch_input
using GeoData:
    GeoDataCatalog,
    DatasetEntry,
    register_dataset!,
    save_catalog,
    default_catalog_path


# ==============================================================================
# Configuration Defaults (Derived from configs/snowcrab.toml)
# ==============================================================================

const DEFAULT_EMBEDDING_LON = (-71.0, -53.0)
const DEFAULT_EMBEDDING_LAT = (40.0, 48.5)
const DEFAULT_STUDY_LON     = (-68.0, -57.0)
const DEFAULT_STUDY_LAT     = (42.0, 47.5)
const DEFAULT_WIND_TIME_ISO = "2020-06-01T00:00:00Z"

"""
    parse_cli_args(args::Vector{String}) -> Dict{Symbol, Any}

Parse command-line arguments for cache initialization.
"""
function parse_cli_args(args::Vector{String})
    cfg = Dict{Symbol, Any}(
        :cache_dir     => joinpath(homedir(), ".geodata", "cache"),
        :catalog_path  => default_catalog_path(),
        :skip_tides    => false,
        :skip_woa      => false,
        :skip_bathy    => false,
        :skip_coast    => false,
        :skip_winds    => false,
        :force         => false,
        :dry_run       => false,
        :embedding_lon => DEFAULT_EMBEDDING_LON,
        :embedding_lat => DEFAULT_EMBEDDING_LAT,
        :wind_time_iso => DEFAULT_WIND_TIME_ISO
    )

    idx = 1
    while idx <= length(args)
        arg = args[idx]
        if arg == "--cache-dir" && idx < length(args)
            idx += 1
            cfg[:cache_dir] = normpath(abspath(args[idx]))
        elseif arg == "--catalog" && idx < length(args)
            idx += 1
            cfg[:catalog_path] = normpath(abspath(args[idx]))
        elseif arg == "--skip-tides"
            cfg[:skip_tides] = true
        elseif arg == "--skip-woa"
            cfg[:skip_woa] = true
        elseif arg == "--skip-bathy"
            cfg[:skip_bathy] = true
        elseif arg == "--skip-coast"
            cfg[:skip_coast] = true
        elseif arg == "--skip-winds"
            cfg[:skip_winds] = true
        elseif arg == "--force"
            cfg[:force] = true
        elseif arg == "--dry-run"
            cfg[:dry_run] = true
        elseif arg in ("-h", "--help")
            println(@doc @__MODULE__)
            exit(0)
        else
            @warn "Unrecognized option ignored: $(arg)"
        end
        idx += 1
    end

    return cfg
end

"""
    fetch_and_register_baseline(cfg::Dict{Symbol, Any}) -> GeoDataCatalog

Execute data acquisition routines and register each artefact in the catalog.
"""
function fetch_and_register_baseline(cfg::Dict{Symbol, Any})
    cache_dir = cfg[:cache_dir]
    catalog_path = cfg[:catalog_path]
    force = cfg[:force]
    dry_run = cfg[:dry_run]

    mkpath(cache_dir)
    mkpath(dirname(catalog_path))

    println("="^80)
    println("GeoData Baseline Cache Initialization")
    println("  Cache Directory : ", cache_dir)
    println("  Catalog File    : ", catalog_path)
    println("  Embedding Domain: lon $(cfg[:embedding_lon]), lat $(cfg[:embedding_lat])")
    println("  Force Re-fetch  : ", force)
    println("  Dry Run         : ", dry_run)
    println("="^80)

    catalog = isfile(catalog_path) ? GeoDataCatalog(catalog_path) :
                                     GeoDataCatalog(catalog_path; description="GeoData Baseline Catalog")

    # Helper to register or update entry
    function record_entry!(entry::DatasetEntry)
        if haskey(catalog.key_index, entry.key)
            old_uuid = catalog.key_index[entry.key]
            println("  Updating existing catalog entry for key: :$(entry.key) (UUID: $(old_uuid))")
            catalog.entries[old_uuid] = entry
            catalog.updated = Dates.now()
        else
            register_dataset!(catalog, entry)
            println("  Registered new dataset: :$(entry.key) -> $(entry.name)")
        end
    end

    # --------------------------------------------------------------------------
    # 1. Natural Earth 10m Coastline
    # --------------------------------------------------------------------------
    if !cfg[:skip_coast]
        coast_path = joinpath(cache_dir, "coastline.parquet")
        println("\n[1/5] Natural Earth 10m Coastline...")
        if dry_run
            println("  [DRY RUN] Would fetch to: $(coast_path)")
        else
            try
                fetch_natural_earth_coastline(
                    lon_range   = cfg[:embedding_lon],
                    lat_range   = cfg[:embedding_lat],
                    resolution  = "10m",
                    output_path = coast_path,
                    margin_deg  = 1.0,
                    verbose     = true
                )
                record_entry!(DatasetEntry(
                    name = "Natural Earth 10m Coastline (NW Atlantic)",
                    key = :natural_earth_coastline,
                    variables = ["geometry"],
                    spatial_bounds = (
                        Float64(cfg[:embedding_lon][1]),
                        Float64(cfg[:embedding_lon][2]),
                        Float64(cfg[:embedding_lat][1]),
                        Float64(cfg[:embedding_lat][2])
                    ),
                    source = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson",
                    location = coast_path,
                    access_method = :file,
                    format = :geoparquet,
                    credits = "Natural Earth, Public Domain",
                    license = "Public Domain",
                    notes = "Clipped 10m vector coastline for Scotian Shelf / NW Atlantic domain."
                ))
            catch err
                @error "Failed to fetch coastline: $(err)"
            end
        end
    end

    # --------------------------------------------------------------------------
    # 2. ETOPO 2022 Bathymetry
    # --------------------------------------------------------------------------
    if !cfg[:skip_bathy]
        bathy_path = joinpath(cache_dir, "bathymetry.zarr")
        println("\n[2/5] ETOPO 2022 15-arcsec Bathymetry...")
        if dry_run
            println("  [DRY RUN] Would fetch to: $(bathy_path)")
        else
            try
                fetch_erddap_bathymetry(
                    lon_range   = cfg[:embedding_lon],
                    lat_range   = cfg[:embedding_lat],
                    output_path = bathy_path,
                    dataset_id  = "ETOPO_2022_v1_15s",
                    stride      = 1,
                    backend     = :zarr,
                    verbose     = true
                )
                record_entry!(DatasetEntry(
                    name = "ETOPO 2022 15-arcsec Relief (NW Atlantic)",
                    key = :bathymetry_etopo2022,
                    variables = ["elevation"],
                    spatial_bounds = (
                        Float64(cfg[:embedding_lon][1]),
                        Float64(cfg[:embedding_lon][2]),
                        Float64(cfg[:embedding_lat][1]),
                        Float64(cfg[:embedding_lat][2])
                    ),
                    vertical_range = (-10000.0, 4000.0),
                    source = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ETOPO_2022_v1_15s.nc",
                    location = bathy_path,
                    access_method = :file,
                    format = :zarr,
                    credits = "NOAA NCEI",
                    license = "Public Domain",
                    notes = "High-resolution 15-arcsecond sea floor relief from NOAA ERDDAP."
                ))
            catch err
                @error "Failed to fetch bathymetry: $(err)"
            end
        end
    end

    # --------------------------------------------------------------------------
    # 3. Open-Meteo ERA5 Surface Winds & Flux
    # --------------------------------------------------------------------------
    if !cfg[:skip_winds]
        wind_path = joinpath(cache_dir, "wind.zarr")
        println("\n[3/5] Surface Winds & Atmospheric State (Open-Meteo ERA5)...")
        if dry_run
            println("  [DRY RUN] Would fetch to: $(wind_path)")
        else
            try
                fetch_open_meteo_winds(
                    lon_range   = cfg[:embedding_lon],
                    lat_range   = cfg[:embedding_lat],
                    time_iso    = cfg[:wind_time_iso],
                    output_path = wind_path,
                    backend     = :zarr,
                    verbose     = true
                )
                record_entry!(DatasetEntry(
                    name = "Open-Meteo ERA5 Atmospheric Forcing",
                    key = :surface_winds_open_meteo,
                    variables = ["tau_x", "tau_y", "t2m", "surface_pressure", "relative_humidity_2m"],
                    temporal_coverage = (
                        DateTime(2020, 6, 1, 0, 0),
                        DateTime(2020, 6, 1, 23, 0)
                    ),
                    spatial_bounds = (
                        Float64(cfg[:embedding_lon][1]),
                        Float64(cfg[:embedding_lon][2]),
                        Float64(cfg[:embedding_lat][1]),
                        Float64(cfg[:embedding_lat][2])
                    ),
                    source = "https://archive-api.open-meteo.com/v1/archive",
                    location = wind_path,
                    access_method = :file,
                    format = :zarr,
                    credits = "Open-Meteo / ECMWF ERA5",
                    license = "CC-BY-4.0",
                    notes = "Surface kinematic stress and bulk thermodynamic variables."
                ))
            catch err
                @error "Failed to fetch surface winds: $(err)"
            end
        end
    end

    # --------------------------------------------------------------------------
    # 4. World Ocean Atlas 2023 Climatology (T, S, O2)
    # --------------------------------------------------------------------------
    if !cfg[:skip_woa]
        println("\n[4/5] World Ocean Atlas 2023 Climatology...")
        if dry_run
            println("  [DRY RUN] Would fetch WOA23 files into: $(cache_dir)")
        else
            try
                res = fetch_woa23(
                    lon_range   = cfg[:embedding_lon],
                    lat_range   = cfg[:embedding_lat],
                    month       = 0,
                    output_dir  = cache_dir,
                    include_o2  = true,
                    backend     = :zarr,
                    verbose     = true
                )
                record_entry!(DatasetEntry(
                    name = "WOA23 Annual Temperature Climatology",
                    key = :woa23_temperature,
                    variables = ["t_an"],
                    spatial_bounds = (
                        Float64(cfg[:embedding_lon][1]),
                        Float64(cfg[:embedding_lon][2]),
                        Float64(cfg[:embedding_lat][1]),
                        Float64(cfg[:embedding_lat][2])
                    ),
                    vertical_range = (0.0, 5500.0),
                    source = "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA",
                    location = res.temperature_file,
                    access_method = :file,
                    format = :zarr,
                    credits = "NOAA NCEI",
                    license = "Public Domain",
                    notes = "WOA23 0.25-degree annual temperature climatology."
                ))
                record_entry!(DatasetEntry(
                    name = "WOA23 Annual Salinity Climatology",
                    key = :woa23_salinity,
                    variables = ["s_an"],
                    spatial_bounds = (
                        Float64(cfg[:embedding_lon][1]),
                        Float64(cfg[:embedding_lon][2]),
                        Float64(cfg[:embedding_lat][1]),
                        Float64(cfg[:embedding_lat][2])
                    ),
                    vertical_range = (0.0, 5500.0),
                    source = "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA",
                    location = res.salinity_file,
                    access_method = :file,
                    format = :zarr,
                    credits = "NOAA NCEI",
                    license = "Public Domain",
                    notes = "WOA23 0.25-degree annual salinity climatology."
                ))
                record_entry!(DatasetEntry(
                    name = "WOA23 Annual Dissolved Oxygen Climatology",
                    key = :woa23_oxygen,
                    variables = ["o_an"],
                    spatial_bounds = (
                        Float64(cfg[:embedding_lon][1]),
                        Float64(cfg[:embedding_lon][2]),
                        Float64(cfg[:embedding_lat][1]),
                        Float64(cfg[:embedding_lat][2])
                    ),
                    vertical_range = (0.0, 5500.0),
                    source = "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA",
                    location = res.oxygen_file,
                    access_method = :file,
                    format = :zarr,
                    credits = "NOAA NCEI",
                    license = "Public Domain",
                    notes = "WOA23 0.25-degree annual dissolved oxygen climatology."
                ))
            catch err
                @error "Failed to fetch WOA23 hydrography: $(err)"
            end
        end
    end

    # --------------------------------------------------------------------------
    # 5. TPXO9-Atlas Tidal Solution
    # --------------------------------------------------------------------------
    if !cfg[:skip_tides]
        tides_target = joinpath(cache_dir, "tpxo9", "u_tpxo9.v1.nc")
        println("\n[5/5] TPXO9-Atlas Global Harmonic Tidal Solution...")
        if dry_run
            println("  [DRY RUN] Would fetch to: $(tides_target)")
        else
            try
                fetch_input(:tides, cache_dir; force = force)
                record_entry!(DatasetEntry(
                    name = "TPXO9-Atlas Global Harmonic Tidal Solution",
                    key = :tpxo9_velocity,
                    variables = ["uRe", "uIm", "vRe", "vIm", "con"],
                    spatial_bounds = (-180.0, 180.0, -90.0, 90.0),
                    source = "https://zenodo.org/records/8074917",
                    location = tides_target,
                    access_method = :file,
                    format = :netcdf,
                    credits = "Egbert and Erofeeva / OSU",
                    license = "CC-BY-4.0",
                    notes = "Global harmonic tidal velocity solution (u_tpxo9.v1.nc)."
                ))
            catch err
                @error "Failed to fetch TPXO9 tidal solution: $(err)"
            end
        end
    end

    # --------------------------------------------------------------------------
    # Persist Catalog
    # --------------------------------------------------------------------------
    if !dry_run
        save_catalog(catalog, catalog_path)
        println("\nSaved updated catalog to $(catalog_path) with $(length(catalog.entries)) entries.")
    end

    return catalog
end

# CLI Entrypoint when invoked as a script
if abspath(PROGRAM_FILE) == @__FILE__
    cfg = parse_cli_args(ARGS)
    fetch_and_register_baseline(cfg)
end
