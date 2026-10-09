# GeoData API Reference

`GeoData.jl` provides backend-agnostic geospatial data access, manipulation, provenance tracking, and metadata cataloging across multiple scientific file formats (Zarr, NetCDF/NCZarr, GeoParquet).

---

## 1. Core Architecture and Types

### Core Types (`GeoData.GeoDataCoreTypes`)

- **`GeoDataset`**: Multi-variable geospatial dataset container:
  - `variables::Dict{String, GeoArray}`: Data variables indexed by name.
  - `coords::Dict{Symbol, GeoArray}`: Coordinate arrays (e.g., `:lon`, `:lat`, `:depth`, `:time`).
  - `dims::Dict{Symbol, Dimension}`: Coordinate dimensions and metadata.
  - `crs::CoordinateSystem`: Coordinate Reference System specifications (e.g. EPSG:4326).
  - `attrs::Dict{String, Any}`: Global dataset attributes.
  - `backend::GeoBackend`: Storage backend handle.
  - `uri::String`: URI or filesystem path to dataset.
- **`GeoArray{T, N}`**: N-dimensional array with dimension tags, CRS, and attributes.
- **`Dimension`**: Dimension specifications (`name`, `size`, `coords`, `units`, `dim_type`).
- **`CoordinateSystem`**: Geodetic or Cartesian CRS descriptor (`crs::String`).

---

## 2. High-Level Dataset Operations

All operations are backend-agnostic and operate uniformly across supported formats.

```julia
using GeoData

# Load from Zarr, NetCDF, or GeoParquet (backend auto-inferred from path or URI)
ds = geoload("data/bathymetry.zarr")

# Coordinate-aware slicing
subset = geoslice(ds, lon=(-68.0, -57.0), lat=(42.0, 47.5))

# Point selection or extraction
val = geovalues(subset, lon=-63.5, lat=44.0)

# Save to target format/backend
geosave("output/subset.zarr", subset; backend=:zarr)
```

### Core Functions

- **`geoload(uri::String; backend=nothing, kwargs...) -> GeoDataset`**: Load dataset from path or URI.
- **`geosave(uri::String, data::GeoDataset; backend=nothing, kwargs...) -> String`**: Write dataset to target backend.
- **`geoslice(data::GeoDataset; dims...) -> GeoDataset`**: Slice dataset along bounding coordinate intervals.
- **`geoselect(data::GeoDataset; dims...) -> GeoDataset`**: Extract nearest or exact coordinate slices.
- **`geovalues(data::GeoDataset, vars=nothing; dims...) -> Dict{String, Any}`**: Interpolate and extract values at specified spatial-temporal coordinates.
- **`georegrid(data::GeoDataset, target_grid; method=:bilinear) -> GeoDataset`**: Regrid fields onto target spatial coordinates.
- **`geojoin(datasets...; on=:time, how=:inner) -> GeoDataset`**: Join datasets along shared dimensions.
- **`bounding_box(data::GeoDataset) -> Dict{Symbol, Tuple}`**: Compute spatial and temporal bounding extents.
- **`variable_stats(data::GeoDataset, var::String) -> NamedTuple`**: Compute summary statistics (`min`, `max`, `mean`, `count`).

---

## 3. Storage Backends

| Backend | Identifier | URI Schemes / Extensions | Format Description |
|:---|:---|:---|:---|
| **Zarr** | `:zarr` | `zarr://`, `.zarr` | Chunked, compressed, cloud-optimized N-D arrays. |
| **NetCDF** | `:ncdatasets` | `netcdf://`, `.nc` | Classical CF-compliant NetCDF3/NetCDF4 data. |
| **NCZarr** | `:nczarr` | `nczarr://` | NetCDF-Zarr cloud storage mapping via NCDatasets. |
| **GeoParquet** | `:geoparquet` | `geoparquet://`, `.parquet` | Columnar vector geometries and tabular attributes. |
| **YAXArrays** | `:yaxarray` | `yaxarray://` | Optional lazy labeled multidimensional arrays. |

---

## 4. Oceanographic Data Acquisition (`GeoData.Data`)

`GeoData.Data` provides deterministic downloaders and loaders for key marine products:

### Coastlines
- **`fetch_natural_earth_coastline(; lon_range, lat_range, resolution="10m", output_path, margin_deg=1.0)`**:
  Downloads and clips Natural Earth vector coastline to specified bounds, writing to GeoParquet.
- **`load_coastline_geodata(path)`**: Reads vector coastline as a `GeoDataset`.
- **`is_marine_water_geodata(lon, lat, coastline_ds)`**: Evaluates whether coordinates fall in marine waters.

### Bathymetry & Topography
- **`fetch_erddap_bathymetry(; lon_range, lat_range, output_path, dataset_id="ETOPO_2022_v1_15s", stride=1)`**:
  Fetches NOAA ETOPO 2022 (15 arc-second relief) from ERDDAP and saves as Zarr or NetCDF.
- **`load_bathymetry_geodata(path)`**: Ingests bathymetry dataset returning `GeoDataset` with `:elevation`.
- **`get_bathymetry_interpolator(ds)`**: Builds a fast continuous spatial interpolator `f(lon, lat)`.

### Atmosphere & Winds
- **`fetch_open_meteo_winds(; lon_range, lat_range, time_iso, output_path)`**:
  Downloads ERA5 surface winds (`u10`, `v10`), air temperature (`t2m`), surface pressure, and radiation fluxes from Open-Meteo archive API.
- **`load_wind_stress_geodata(path)`**: Loads wind stress fields and kinematic stresses.
- **`build_bulk_surface_flux_geodata(path)`**: Computes radiative and turbulent bulk heat fluxes.

### Hydrography & Tracers
- **`fetch_woa23(; lon_range, lat_range, month=0, output_dir, include_o2=true)`**:
  Retrieves NOAA World Ocean Atlas 2023 climatological temperature, salinity, and oxygen.
- **`load_woa23_interpolators(res)`**: Returns continuous 3D tracer interpolators `f(lon, lat, depth)`.
- **`fetch_boundary_hydrography_geodata(...)`**: Prepares boundary conditions for ocean simulations.
- **`fetch_copernicus_physics_subset(; lon_range, lat_range, start_date, end_date, output_path, dataset_id, ...)`**:
  Downloads regional 3D ocean physics (temperature, salinity, currents) from Copernicus Marine (CMEMS).
- **`copernicus_credentials()`**: Discovers and decodes local Copernicus Marine authentication tokens.

### Lakehouse Architecture & Medallion Storage Tiers (`GeoData.LakehouseModule`)
`GeoData.jl` serves as a scientific lakehouse repository organized into three standard medallion tiers:
- **Bronze (`:raw`)**: Raw, immutable observations and external input files directly downloaded from data providers (NOAA ETOPO, Open-Meteo ERA5, WOA23, Natural Earth).
- **Silver (`:processed`)**: Harmonized, quality-controlled, regularized Zarr and GeoParquet grids suitable for direct numerical modeling.
- **Gold (`:derived`)**: Post-processed model predictions, dispersal kernels, GAM interpolations, and connectivity matrices produced by external consumer packages (such as `ParticleTracking.jl`) and registered back into the repository for consumption across other projects.

#### Lakehouse Infrastructure Functions
- **`lakehouse_root_dir() -> String`**: Root filesystem directory for lakehouse storage (default: `~/.geodata` or `ENV["GEODATA_LAKEHOUSE_ROOT"]`).
- **`lakehouse_path(key::Symbol; tier::Symbol = :processed, ext::String = "zarr") -> String`**: Computes canonical tiered filesystem storage path `<lakehouse_root>/<tier>/<key>.<ext>`.
- **`compute_lakehouse_checksum(path::String) -> String`**: Computes deterministic SHA-256 digest of file or Zarr directory store.
- **`geopublish_dataset!(catalog, location; key, name, tier, producer, derived_from, ...)`**: High-level publication and catalog registration.
- **`geofetch_dataset(catalog, key) -> DatasetEntry`**: Retrieve registered dataset with access tracking.

---

## 5. Metadata Catalog & Lineage (`GeoData.CatalogModule`)

The `GeoDataCatalog` provides database-like registration, discovery, querying, and persistent lineage tracking for local and remote geospatial datasets.

### Struct: `DatasetEntry`

```julia
@kwdef mutable struct DatasetEntry
    uuid::UUID                                   # Unique entry identifier
    name::String                                 # Descriptive title
    key::Symbol                                  # Stable lookup symbol (e.g. :bathymetry_etopo2022)
    variables::Vector{String}                    # Available field variables
    dimensions::Dict{String, Dimension}          # Dimension specifications
    temporal_coverage::Tuple{DateTime, DateTime} # Temporal validity span
    spatial_bounds::Tuple{Float64, Float64, Float64, Float64} # (min_lon, max_lon, min_lat, max_lat)
    vertical_range::Tuple{Float64, Float64}      # (min_depth, max_depth) in positive meters
    source::String                               # Origin URL or DOI
    location::String                             # Local storage path or remote URI
    access_method::Symbol                        # :file, :http, :opendap, :s3
    format::Symbol                               # :zarr, :netcdf, :geoparquet
    tier::Symbol                                 # Lakehouse tier: :raw, :processed, :derived
    producer::String                             # Producing project/organization
    derived_from::Vector{Symbol}                 # Lineage: upstream catalog dataset keys
    checksum::String                             # SHA-256 digest of payload
    credits::String                              # Provenance attribution
    license::String                              # Distribution license
    notes::String                                # Notes and specifications
end
```

### Catalog Functions

#### Construction and Persistence
- **`GeoDataCatalog(catalog_file::String=""; kwargs...)`**: Initialize or load catalog from disk.
- **`load_catalog!(catalog::GeoDataCatalog, file::String)`**: Ingest dataset entries from JSON file.
- **`save_catalog(catalog::GeoDataCatalog, file::String="")`**: Persist catalog to JSON file.
- **`default_catalog_path() -> String`**: Returns `~/.geodata/catalog.json`.
- **`get_global_catalog() -> GeoDataCatalog`**: Retrieve singleton global catalog instance.

#### Dataset Registration, Publishing, and Queries
- **`register_dataset!(catalog::GeoDataCatalog, entry::DatasetEntry)::UUID`**: Register dataset.
- **`unregister_dataset!(catalog::GeoDataCatalog, uuid::UUID)::Bool`**: Remove dataset from catalog.
- **`update_dataset!(catalog::GeoDataCatalog, uuid::UUID, updates...) -> DatasetEntry`**: Modify entry metadata.
- **`geopublish!(catalog::GeoDataCatalog, entry_or_path; key, tier, ...)`**: Publish dataset into catalog.
- **`geofetch(catalog::GeoDataCatalog, key::Symbol) -> DatasetEntry`**: Lookup dataset and track usage.
- **`get_dataset_by_key(catalog::GeoDataCatalog, key::Symbol) -> DatasetEntry`**: Fast lookup by key.
- **`get_dataset(catalog::GeoDataCatalog, uuid::UUID) -> DatasetEntry`**: Lookup by UUID.
- **`get_datasets_by_name(catalog::GeoDataCatalog, name::String) -> Vector{DatasetEntry}`**: Search by name.
- **`find_datasets(catalog::GeoDataCatalog; key, tier, producer, derived_from, ...)`**:
  Multi-criteria search filtering by tier, upstream lineage, coordinates, variables, date bounds, format, and priority.
- **`get_best_dataset(catalog::GeoDataCatalog; kwargs...) -> DatasetEntry`**: Retrieve highest-priority match.

#### Usage Tracking & Export
- **`access_dataset!(catalog::GeoDataCatalog, uuid::UUID)`**: Increment access counter and update timestamp.
- **`record_download!(catalog::GeoDataCatalog, uuid::UUID, bytes::UInt64)`**: Log transfer metrics.
- **`export_catalog_to_csv(catalog::GeoDataCatalog, file::String)`**: Export catalog table to CSV.
- **`get_stats(catalog::GeoDataCatalog) -> CatalogStats`**: Retrieve aggregate catalog statistics.

---

## 6. Provenance & Reproducibility (`GeoData.Data.GeoDataManifest`)

Every physical input consumed by simulation workflows is declared in `DATA_SOURCES`:

- **`DATA_SOURCES::Dict{Symbol, DataSource}`**: Declarative registry of verified physical inputs (`:bathymetry`, `:surface_winds`, `:hydrography`, `:coastline`, `:tides`, etc.).
- **`fetch_input(key::Symbol, output_dir::String; force=false) -> String`**:
  Downloads and caches the declared data input into `<output_dir>/inputs/`.
- **`data_provenance(output_dir; config_path, data_sources, extra) -> String`**:
  Writes `data_provenance.json` detailing exact files read, file sizes, and SHA-256 digests.
- **`describe_data_sources(io=stdout)`**: Prints human-readable data source registry.

---

## 7. Baseline Initialization Script

The script `scripts/fetch_baseline_geodata.jl` automates populating the local database cache:

```bash
julia --project=. scripts/fetch_baseline_geodata.jl [options]
```

### Options:
- `--cache-dir <dir>`: Target directory for cache files (default: `~/.geodata/cache`).
- `--catalog <path>`: Path to catalog JSON (default: `~/.geodata/catalog.json`).
- `--skip-tides`: Skip TPXO9 tidal velocity download.
- `--force`: Force re-download even if cache files exist.
- `--dry-run`: Display planned actions without executing network transfers.

---

## 8. Hierarchical Analytical Storage (`GeoStorage`)

`GeoStorage` provides an abstraction over multi-group scientific archives (such as simulation runs, Lagrangian particle trajectories, metrics, and demographic connectivity matrices) across chunked tensor stores (Zarr) and partitioned tabular collections (GeoParquet).

```julia
using GeoData

# Open or initialize storage container
st = open_geostorage("outputs/larval_trajectories.zarr"; backend=:zarr, read_only=false)

# Group navigation and variable serialization
write_storage_variable!(st, "trajectories/run1/lons", [-64.2, -64.1, -64.0])
lons = read_storage_variable(st, "trajectories/run1/lons")

# Release locks and flush buffers
close_geostorage(st)
```

### Core API Functions:
- **`open_geostorage(uri::AbstractString; backend=:auto, read_only=false, create=true) -> GeoStorage`**: Open or initialize storage at `uri`.
- **`close_geostorage(storage::GeoStorage) -> Nothing`**: Flush buffers and release file locks.
- **`create_storage_group(storage::GeoStorage, group_path::AbstractString)`**: Create nested hierarchical groups.
- **`has_storage_group(storage::GeoStorage, group_path::AbstractString) -> Bool`**: Query group presence.
- **`write_storage_variable!(storage::GeoStorage, var_path::AbstractString, data; chunks=nothing, compressor=nothing)`**: Write array, vector, or table.
- **`read_storage_variable(storage::GeoStorage, var_path::AbstractString)`**: Read array or table at `var_path`.

---

## 9. Unified Environmental Data Cube (`GeoData.Data.RegionalCube`)

The regional cube engine consolidates static geophysics (ETOPO bathymetry), physical oceanography (GLORYS12V1 / WOA23 $T, S, u, v$), and marine biogeochemistry (WOA23 / CMEMS $[\mathrm{O}_2], \mathrm{pH}, \mathrm{NO}_3, \mathrm{PO}_4, \mathrm{Si}, \text{Chl-}a$) into unified, self-describing chunked Zarr cubes.

### Configuration Schema (`RegionalCubeConfig`)

```julia
cfg = RegionalCubeConfig(
    region_name = "scotian_shelf",
    lon_range = (-71.0, -53.0),
    lat_range = (40.0, 48.5),
    resolution_deg = 1.0 / 12.0,
    depth_levels = standard_ocean_depths(),
    time_mode = :climatology, # or :timeseries, :hybrid
    time_step = :climatology_monthly,
    climatology_months = collect(1:12),
    variables = [:elevation, :bottom_depth, :temperature, :salinity, :u, :v, :dissolved_oxygen, :ph],
    output_path = "data/cubes/scotian_shelf_environmental.zarr"
)
```

### Ingestion & Bioenergetic Functions:
- **`load_cube_config(path::AbstractString) -> RegionalCubeConfig`**: Load declarative TOML specification.
- **`standard_ocean_depths() -> Vector{Float64}`**: 47 oceanographic depth levels (surface-refined down to 4000 m).
- **`assimilate_regional_cube(config::RegionalCubeConfig; verbose=true) -> String`**: Ingest, harmonize, and serialize data cube to Zarr store.
- **`build_cube_hsi_evaluator(cube_path::AbstractString) -> Function`**: Build high-performance callable `(lon, lat, depth, month) -> Float64` incorporating bathymetric and thermal habitat suitability constraints.
- **`evaluate_bioenergetic_scope(temp, dissolved_o2, ph; t_opt=3.0, t_max=10.0, k_o2=60.0, ph_ref=8.1) -> Float64`**: Compute physiological metabolic scope $\mu \in [0, 1]$:
  $$\mu = \mu_{\max}(T) \cdot \frac{[\mathrm{O}_2]}{K_{O2} + [\mathrm{O}_2]} \cdot f(\mathrm{pH})$$

### Standalone CLI Execution:
```bash
julia --project=projects/GeoData scripts/assimilate_environmental_cube.jl --config=configs/regions/scotian_shelf.toml
```

