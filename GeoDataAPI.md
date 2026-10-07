# GeoData API Documentation

## GeoDataCatalog

The `GeoDataCatalog` is a central registry for managing and querying geospatial datasets.

### Structure

The `GeoDataCatalog` mutable struct contains:

- `entries`: Dictionary of dataset entries keyed by UUID.
- `key_index`: Dictionary for quick lookup by dataset key (Symbol).
- `name_index`: Dictionary for quick lookup by dataset name (String).
- `temporal_index`: Interval Tree for efficient temporal range queries.
- `spatial_index`: Rect Tree for efficient spatial bounding box queries.
- `variable_index`: Inverted Index for efficient variable-based lookups.
- ... (other fields for configuration, statistics, metadata)

### Functions

#### Construction and Persistence

- `GeoDataCatalog(catalog_file::String=""; kwargs...)`: Create a new catalog, optionally loading from a file.
- `load_catalog!(catalog::GeoDataCatalog, file::String)`: Load catalog entries from a JSON file, initializing and populating indexes.
- `save_catalog(catalog::GeoDataCatalog, file::String="")`: Save catalog entries to a JSON file.

#### Catalog Modification

- `register_dataset!(catalog::GeoDataCatalog, entry::DatasetEntry)::UUID`: Register a new dataset, updating all indexes.
- `unregister_dataset!(catalog::GeoDataCatalog, uuid::UUID)::Bool`: Remove a dataset, updating all indexes.
- `update_dataset!(catalog::GeoDataCatalog, uuid::UUID, updates::Pair{Symbol, Any}...)::DatasetEntry`: Update a dataset, updating indexes for changed fields.

#### Query and Discovery

- `get_dataset(catalog::GeoDataCatalog, uuid::UUID)::DatasetEntry`: Get a dataset by UUID.
- `get_dataset_by_key(catalog::GeoDataCatalog, key::Symbol)::DatasetEntry`: Get a dataset by key.
- `get_datasets_by_name(catalog::GeoDataCatalog, name::String)::Vector{DatasetEntry}`: Get datasets by name.
- `find_datasets(; name::Union{Nothing, String}=nothing, key::Union{Nothing, Symbol}=nothing, variables::Union{Nothing, Vector{String}}=nothing, temporal_range::Union{Nothing, Tuple{DateTime, DateTime}}=nothing, spatial_bounds::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}}=nothing, vertical_range::Union{Nothing, Tuple{Float64, Float64}}=nothing, format::Union{Nothing, Symbol}=nothing, access_method::Union{Nothing, Symbol}=nothing, active_only::Bool=true, limit::Union{Nothing, Int}=nothing)::Vector{DatasetEntry}`: Find datasets matching various criteria.
- `get_best_dataset(; kwargs...)::DatasetEntry`: Get the single best matching dataset.

#### Usage Tracking

- `access_dataset!(catalog::GeoDataCatalog, uuid::UUID)::DatasetEntry`: Record dataset access.
- `record_download!(catalog::GeoDataCatalog, uuid::UUID, bytes_downloaded::UInt64)::DatasetEntry`: Record dataset download.

#### Statistics

- `get_stats(catalog::GeoDataCatalog)::CatalogStats`: Get catalog statistics.

#### Import/Export

- `export_catalog_to_csv(catalog::GeoDataCatalog, file::String)`: Export catalog to CSV.

#### Global Catalog

- `get_global_catalog()::GeoDataCatalog`: Get the global catalog instance.
- `init_global_catalog!(catalog_file::String=""; kwargs...)`: Initialize the global catalog.

### Indexing for Efficient Queries

The catalog maintains three indexes for efficient querying:

1. **Temporal Index** (`IntervalTree{DateTime, UUID}`): Allows fast lookup of datasets overlapping a given time range.
2. **Spatial Index** (`RectTree{Float64, UUID}`): Allows fast lookup of datasets overlapping a given spatial bounding box.
3. **Variable Index** (`InvertedIndex{String, UUID}`): Allows fast lookup of datasets containing a given variable.

These indexes are automatically maintained when datasets are registered, unregistered, or updated, and are rebuilt when the catalog is loaded from file.
