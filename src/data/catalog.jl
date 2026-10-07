"""
    catalog.jl

    Dataset catalog for GeoData - a centralized registry for tracking and managing
    geospatial datasets across multiple sources and backends.

    Provides database-like functionality for discovering, accessing, and managing
    oceanographic and geographic datasets with metadata, location tracking,
    versioning, and querying capabilities.
    """
module GeoDataCatalog

using ..GeoDataCoreTypes
using ..GeoDataRegistry
using Dates
using UUIDs
using SHA
using JSON3
using Interpolations
using DataFrames
using DelimitedFiles
using Downloads
using TranscodingStreams
using HTTP

# ============================================================================
# Catalog Entry Types
# ============================================================================

"""
    DatasetEntry

    A single entry in the dataset catalog representing a discoverable dataset.
    
    This is the core metadata record that tracks what datasets are available,
    where they are located, and how to access them.
"""
@kwdef struct DatasetEntry
    # Core identification
    uuid::UUID = UUIDs.uuid4()                    # Unique identifier for this entry
    name::String = ""                             # Human-readable name
    key::Symbol = Symbol("")                      # Short key for easy reference
    
    # Dataset characteristics
    variables::Vector{String} = String[]          # List of available variables
    dimensions::Dict{String, Dimension} = Dict{String, Dimension}()  # Dimension info
    temporal_coverage::Tuple{DateTime, DateTime} = (DateTime(1970,1,1), DateTime(1970,1,1))  # Time range
    spatial_bounds::Tuple{Float64, Float64, Float64, Float64} = (0.0, 0.0, 0.0, 0.0)  # (min_lon, max_lon, min_lat, max_lat)
    vertical_range::Tuple{Float64, Float64} = (0.0, 0.0)  # (min_depth, max_depth) in meters (positive down)
    
    # Source and access information
    source::String = ""                           # Original source (URL, DOI, etc.)
    location::String = ""                         # Local or remote location where dataset is stored/cached
    access_method::Symbol = :file                 # :file, :http, :opendap, :s3, etc.
    format::Symbol = :unknown                     # :zarr, :netcdf, :geoparquet, :csv, etc.
    
    # Metadata and provenance
    attributes::Dict{String, Any} = Dict{String, Any}()  # Global attributes from source
    credits::String = ""                          # Attribution/citation requirements
    license::String = ""                          # License information
    notes::String = ""                            # Additional notes
    
    # Catalog management
    registered::DateTime = Dates.now()            # When this entry was added to catalog
    updated::DateTime = Dates.now()               # When this entry was last updated
    version::String = "1.0"                       # Version of this catalog entry
    active::Bool = true                           # Whether this entry is currently active
    priority::Int = 0                             # Priority for selection (higher = preferred)
    
    # Usage tracking
    access_count::Int = 0                         # Number of times this dataset has been accessed
    last_accessed::DateTime = Dates.now()         # When it was last accessed
    total_download_size::UInt64 = 0               # Total bytes downloaded for this dataset
end

"""
    CatalogStats

    Statistics about the catalog itself.
"""
@kwdef struct CatalogStats
    total_entries::Int = 0
    active_entries::Int = 0
    total_size_estimate::UInt64 = 0
    last_updated::DateTime = Dates.now()
    unique_formats::Set{Symbol} = Set{Symbol}()
    unique_access_methods::Set{Symbol} = Set{Symbol}()
end

# ============================================================================
# Catalog Storage and Management
# ============================================================================

"""
    GeoDataCatalog

    Main catalog structure that holds all dataset entries and provides
    querying and management capabilities.
"""
mutable struct GeoDataCatalog
    # Storage
    entries::Dict{UUID, DatasetEntry} = Dict{UUID, DatasetEntry}()
    key_index::Dict{Symbol, UUID} = Dict{Symbol, UUID}()
    name_index::Dict{String, Vector{UUID}} = Dict{String, Vector{UUID}}()
    
    # Configuration
    catalog_file::String = ""                     # Path to persistent catalog storage
    auto_save::Bool = true                        # Whether to automatically save to disk
    readonly::Bool = false                        # Whether catalog is read-only
    
    # Statistics (cached)
    stats::CatalogStats = CatalogStats()
    stats_dirty::Bool = true                      # Whether stats need recalculation
    
    # Metadata
    description::String = ""                      # Description of this catalog
    version::String = "1.0"                       # Catalog schema version
    created::DateTime = Dates.now()               # When catalog was created
    updated::DateTime = Dates.now()               # When catalog was last modified
end

# ============================================================================
# Catalog Construction and Persistence
# ============================================================================

"""
    GeoDataCatalog(catalog_file::String=""; kwargs...)

    Create a new dataset catalog, optionally loading from a file.
"""
function GeoDataCatalog(catalog_file::String=""; kwargs...)
    catalog = GeoDataCatalog(catalog_file = catalog_file)
    
    # Apply keyword arguments
    for (k, v) in kwargs
        setfield!(catalog, k, v)
    end
    
    # Load from file if specified and exists
    if !isempty(catalog_file) && isfile(catalog_file)
        load_catalog!(catalog, catalog_file)
    end
    
    return catalog
end

"""
    load_catalog!(catalog::GeoDataCatalog, file::String)

    Load catalog entries from a JSON file.
"""
function load_catalog!(catalog::GeoDataCatalog, file::String)
    if !isfile(file)
        throw(ArgumentError("Catalog file does not exist: $file"))
    end
    
    try
        # Read JSON data
        data = JSON3.read(read(file, String))
        
        # Restore catalog metadata
        catalog.description = get(data, "description", "")
        catalog.version = get(data, "version", "1.0")
        catalog.created = Dates.DateTime(get(data, "created", string(Dates.now())))
        catalog.updated = Dates.DateTime(get(data, "updated", string(Dates.now())))
        catalog.readonly = get(data, "readonly", false)
        catalog.auto_save = get(data, "auto_save", true)
        
        # Restore entries
        entries_data = get(data, "entries", Dict{String, Any}())
        for (uuid_str, entry_data) in entries_data
            uuid = UUIDs.UUID(uuid_str)
            entry = _deserialize_entry(entry_data)
            catalog.entries[uuid] = entry
            
            # Update basic indices
            if !isempty(entry.key)
                catalog.key_index[entry.key] = uuid
            end
            if !isempty(entry.name)
                if !haskey(catalog.name_index, entry.name)
                    catalog.name_index[entry.name] = UUID[]
                end
                push!(catalog.name_index[entry.name], uuid)
            end
        end
        
        # Initialize and populate indexes for efficient querying
        catalog.temporal_index = IntervalTree{DateTime, UUID}()
        catalog.spatial_index = RectTree{Float64, UUID}()
        catalog.variable_index = InvertedIndex{String, UUID}()
        for (uuid, entry) in catalog.entries
            insert!(catalog.temporal_index, entry.temporal_coverage[1], entry.temporal_coverage[2], uuid)
            insert!(catalog.spatial_index, entry.spatial_bounds[1], entry.spatial_bounds[3], entry.spatial_bounds[2], entry.spatial_bounds[4], uuid)
            for variable in entry.variables
                insert!(catalog.variable_index, variable, uuid)
            end
        end
        
        # Rebuild statistics
        _update_stats!(catalog)
        
        catalog.updated = Dates.now()
        println("Loaded $(length(catalog.entries)) dataset entries from $file")
        
    catch e
        throw(RuntimeError("Failed to load catalog from $file: $e"))
    end
    
    return catalog
function save_catalog(catalog::GeoDataCatalog, file::String="")
    target_file = isempty(file) ? catalog.catalog_file : file
    
    if isempty(target_file)
        throw(ArgumentError("No catalog file specified for saving"))
    end
    
    if catalog.readonly
        throw(ArgumentError("Catalog is read-only and cannot be saved"))
    end
    
    try
        # Prepare data for serialization
        data = Dict{String, Any}(
            "description" => catalog.description,
            "version" => catalog.version,
            "created" => string(catalog.created),
            "updated" => string(catalog.updated),
            "readonly" => catalog.readonly,
            "auto_save" => catalog.auto_save,
            "entries" => Dict{String, Any}()
        )
        
        # Serialize entries
        for (uuid, entry) in catalog.entries
            data["entries"][string(uuid)] = _serialize_entry(entry)
        end
        
        # Write to file
        mkpath(dirname(target_file))  # Ensure directory exists
        write(target_file, JSON3.write(data, 4))  # Pretty print with 4-space indent
        
        catalog.updated = Dates.now()
        println("Saved $(length(catalog.entries)) dataset entries to $target_file")
        
    catch e
        throw(RuntimeError("Failed to save catalog to $target_file: $e"))
    end
    
    return nothing
end

# ============================================================================
# Entry Serialization/Deserialization
# ============================================================================

function _serialize_entry(entry::DatasetEntry)::Dict{String, Any}
    return Dict{String, Any}(
        "uuid" => string(entry.uuid),
        "name" => entry.name,
        "key" => string(entry.key),
        "variables" => entry.variables,
        "dimensions" => Dict{String, Any}([string(k) => _serialize_dimension(v) for (k, v) in entry.dimensions]),
        "temporal_coverage" => [string(entry.temporal_coverage[1]), string(entry.temporal_coverage[2])],
        "spatial_bounds" => collect(entry.spatial_bounds),
        "vertical_range" => collect(entry.vertical_range),
        "source" => entry.source,
        "location" => entry.location,
        "access_method" => string(entry.access_method),
        "format" => string(entry.format),
        "attributes" => entry.attributes,
        "credits" => entry.credits,
        "license" => entry.license,
        "notes" => entry.notes,
        "registered" => string(entry.registered),
        "updated" => string(entry.updated),
        "version" => entry.version,
        "active" => entry.active,
        "priority" => entry.priority,
        "access_count" => entry.access_count,
        "last_accessed" => string(entry.last_accessed),
        "total_download_size" => entry.total_download_size
    )
end

function _deserialize_entry(data::Dict{String, Any})::DatasetEntry
    # Handle dimensions specially
    dimensions_dict = Dict{String, Dimension}()
    if haskey(data, "dimensions")
        for (k, v) in data["dimensions"]
            dimensions_dict[k] = _deserialize_dimension(v)
        end
    end
    
    return DatasetEntry(
        uuid = UUIDs.UUID(get(data, "uuid", string(UUIDs.uuid4()))),
        name = get(data, "name", ""),
        key = Symbol(get(data, "key", "")),
        variables = get(data, "variables", String[]),
        dimensions = dimensions_dict,
        temporal_coverage = (
            Dates.DateTime(get(data, "temporal_coverage", ["1970-01-01T00:00:00", "1970-01-01T00:00:00"])[1]),
            Dates.DateTime(get(data, "temporal_coverage", ["1970-01-01T00:00:00", "1970-01-01T00:00:00"])[2])
        ),
        spatial_bounds = Tuple(get(data, "spatial_bounds", [0.0, 0.0, 0.0, 0.0])...),
        vertical_range = Tuple(get(data, "vertical_range", [0.0, 0.0])...),
        source = get(data, "source", ""),
        location = get(data, "location", ""),
        access_method = Symbol(get(data, "access_method", "file")),
        format = Symbol(get(data, "format", "unknown")),
        attributes = get(data, "attributes", Dict{String, Any}()),
        credits = get(data, "credits", ""),
        license = get(data, "license", ""),
        notes = get(data, "notes", ""),
        registered = Dates.DateTime(get(data, "registered", string(Dates.now()))),
        updated = Dates.DateTime(get(data, "updated", string(Dates.now()))),
        version = get(data, "version", "1.0"),
        active = get(data, "active", true),
        priority = get(data, "priority", 0),
        access_count = get(data, "access_count", 0),
        last_accessed = Dates.DateTime(get(data, "last_accessed", string(Dates.now()))),
        total_download_size = get(data, "total_download_size", UInt64(0))
    )
end

function _serialize_dimension(dim::Dimension)::Dict{String, Any}
    return Dict{String, Any}(
        "name" => string(dim.name),
        "size" => dim.size,
        "coords" => dim.coords !== nothing ? collect(dim.coords) : nothing,
        "units" => dim.units,
        "standard_name" => dim.standard_name,
        "dim_type" => string(dim.dim_type),
        "calendar" => dim.calendar,
        "is_unlimited" => dim.is_unlimited
    )
end

function _deserialize_dimension(data::Dict{String, Any})::Dimension
    return Dimension(
        name = Symbol(get(data, "name", "")),
        size = get(data, "size", nothing),
        coords = get(data, "coords") !== nothing ? collect(get(data, "coords")) : nothing,
        units = get(data, "units", ""),
        standard_name = get(data, "standard_name", nothing),
        dim_type = DimensionType(get(data, "dim_type", "DIM_GENERIC")),
        calendar = get(data, "calendar", nothing),
        is_unlimited = get(data, "is_unlimited", false)
    )
end

# ============================================================================
# Catalog Statistics
# ============================================================================

"""
    _update_stats!(catalog::GeoDataCatalog)

    Recalculate cached statistics for the catalog.
"""
function _update_stats!(catalog::GeoDataCatalog)
    active_count = 0
    total_size = UInt64(0)
    formats = Set{Symbol}()
    methods = Set{Symbol}()
    
    for entry in values(catalog.entries)
        if entry.active
            active_count += 1
        end
        # Rough size estimate (could be improved)
        total_size += entry.total_download_size
        push!(formats, entry.format)
        push!(methods, entry.access_method)
    end
    
    catalog.stats = CatalogStats(
        total_entries = length(catalog.entries),
        active_entries = active_count,
        total_size_estimate = total_size,
        last_updated = Dates.now(),
        unique_formats = formats,
        unique_access_methods = methods
    )
    
    catalog.stats_dirty = false
    return nothing
end

"""
    get_stats(catalog::GeoDataCatalog) -> CatalogStats

    Get current catalog statistics (recalculates if dirty).
"""
function get_stats(catalog::GeoDataCatalog)
    if catalog.stats_dirty
        _update_stats!(catalog)
    end
    return catalog.stats
end

# ============================================================================
# Catalog Modification Functions
# ============================================================================

"""
    register_dataset!(catalog::GeoDataCatalog, entry::DatasetEntry) -> UUID

    Register a new dataset entry in the catalog.
    Returns the UUID of the registered entry.
"""
function register_dataset!(catalog::GeoDataCatalog, entry::DatasetEntry)::UUID
    if catalog.readonly
        throw(ArgumentError("Catalog is read-only and cannot be modified"))
    end
    
    # Generate UUID if not provided
    if entry.uuid == UUIDs.uuid4()
        entry.uuid = UUIDs.uuid4()
    end
    
    # Check for conflicts with key index
    if !isempty(entry.key) && haskey(catalog.key_index, entry.key)
        existing_uuid = catalog.key_index[entry.key]
        if haskey(catalog.entries, existing_uuid)
            existing_entry = catalog.entries[existing_uuid]
            throw(ArgumentError(
                "Dataset key '$(entry.key)' already exists for entry '$(existing_entry.name)' " *
                "(UUID: $(existing_entry.uuid)). Use a different key or unregister first."))
        end
    end
    
    # Store the entry
    catalog.entries[entry.uuid] = entry
    
    # Update basic indices
    if !isempty(entry.key)
        catalog.key_index[entry.key] = entry.uuid
    end
    if !isempty(entry.name)
        if !haskey(catalog.name_index, entry.name)
            catalog.name_index[entry.name] = UUID[]
        end
        push!(catalog.name_index[entry.name], entry.uuid)
    end
    
    # Update indexes for efficient querying
    insert!(catalog.temporal_index, entry.temporal_coverage[1], entry.temporal_coverage[2], entry.uuid)
    insert!(catalog.spatial_index, entry.spatial_bounds[1], entry.spatial_bounds[3], entry.spatial_bounds[2], entry.spatial_bounds[4], entry.uuid)
    for variable in entry.variables
        insert!(catalog.variable_index, variable, entry.uuid)
    end
    
    # Update timestamps
    entry.registered = Dates.now()
    entry.updated = Dates.now()
    catalog.updated = Dates.now()
    catalog.stats_dirty = true
    
    # Auto-save if enabled
    if catalog.auto_save && !isempty(catalog.catalog_file)
        save_catalog(catalog)
    end
    
    println("Registered dataset: '$(entry.name)' (key: :$(entry.key), UUID: $(entry.uuid))")
    return entry.uuid
end
function unregister_dataset!(catalog::GeoDataCatalog, uuid::UUID)::Bool
    if catalog.readonly
        throw(ArgumentError("Catalog is read-only and cannot be modified"))
    end
    
    if !haskey(catalog.entries, uuid)
        return false
    end
    
    entry = catalog.entries[uuid]
    
    # Remove from basic indices
    if !isempty(entry.key) && haskey(catalog.key_index, entry.key) && catalog.key_index[entry.key] == uuid
        delete!(catalog.key_index, entry.key)
    end
    if !isempty(entry.name) && haskey(catalog.name_index, entry.name)
        idxs = catalog.name_index[entry.name]
        filter!(x -> x != uuid, idxs)
        if isempty(idxs)
            delete!(catalog.name_index, entry.name)
        end
    end
    
    # Remove from indexes for efficient querying
    delete!(catalog.temporal_index, entry.temporal_coverage[1], entry.temporal_coverage[2], entry.uuid)
    delete!(catalog.spatial_index, entry.spatial_bounds[1], entry.spatial_bounds[3], entry.spatial_bounds[2], entry.spatial_bounds[4], entry.uuid)
    for variable in entry.variables
        delete!(catalog.variable_index, variable, entry.uuid)
    end
    
    # Remove entry
    delete!(catalog.entries, uuid)
    
    # Update timestamps
    catalog.updated = Dates.now()
    catalog.stats_dirty = true
    
    # Auto-save if enabled
    if catalog.auto_save && !isempty(catalog.catalog_file)
        save_catalog(catalog)
    end
    
    println("Unregistered dataset: '$(entry.name)' (UUID: $(uuid))")
    return true
end
function update_dataset!(catalog::GeoDataCatalog, uuid::UUID, updates::Pair{Symbol, Any}...)::DatasetEntry
    if catalog.readonly
        throw(ArgumentError("Catalog is read-only and cannot be modified"))
    end
    
    if !haskey(catalog.entries, uuid)
        throw(ArgumentError("Dataset entry not found: $uuid"))
    end
    
    entry = catalog.entries[uuid]
    
    # Apply updates
    for (key, value) in updates
        if hasfield(DatasetEntry, key)
            setfield!(entry, key, value)
        else
            @warn "Ignoring unknown field: $key"
        end
    end
    
    # Update timestamp
    entry.updated = Dates.now()
    catalog.updated = Dates.now()
    catalog.stats_dirty = true
    
    # Handle index updates if key, name, temporal_coverage, spatial_bounds, or variables changed
    if haskey(dict(updates), :key) || haskey(dict(updates), :name) || haskey(dict(updates), :temporal_coverage) || haskey(dict(updates), :spatial_bounds) || haskey(dict(updates), :variables)
       
        # Create a copy of the entry as it was before updates for index removal
        old_entry = DatasetEntry(
            uuid = entry.uuid,
            name = entry.name,
            key = entry.key,
            variables = copy(entry.variables),
            dimensions = copy(entry.dimensions),
            temporal_coverage = entry.temporal_coverage,
            spatial_bounds = entry.spatial_bounds,
            vertical_range = entry.vertical_range,
            source = entry.source,
            location = entry.location,
            access_method = entry.access_method,
            format = entry.format,
            attributes = copy(entry.attributes),
            credits = entry.credits,
            license = entry.license,
            notes = entry.notes,
            registered = entry.registered,
            updated = entry.updated,
            version = entry.version,
            active = entry.active,
            priority = entry.priority,
            access_count = entry.access_count,
            last_accessed = entry.last_accessed,
            total_download_size = entry.total_download_size
        )
       
        # Apply updates has already been done above
       
        # Remove old values from basic indices
        if !isempty(old_entry.key)
            delete!(catalog.key_index, old_entry.key)
        end
        if !isempty(old_entry.name)
            idxs = catalog.name_index[old_entry.name]
            filter!(x -> x != entry.uuid, idxs)  # Note: we want to remove the UUID, not the old_entry.uuid
            if isempty(idxs)
                delete!(catalog.name_index, old_entry.name)
            end
        end
       
        # Remove old values from indexes
        delete!(catalog.temporal_index, old_entry.temporal_coverage[1], old_entry.temporal_coverage[2], entry.uuid)
        delete!(catalog.spatial_index, old_entry.spatial_bounds[1], old_entry.spatial_bounds[3], old_entry.spatial_bounds[2], old_entry.spatial_bounds[4], entry.uuid)
        for variable in old_entry.variables
            delete!(catalog.variable_index, variable, entry.uuid)
        end
       
        # Add new values to basic indices
        if !isempty(entry.key)
            catalog.key_index[entry.key] = entry.uuid
        end
        if !isempty(entry.name)
            if !haskey(catalog.name_index, entry.name)
                catalog.name_index[entry.name] = UUID[]
            end
            push!(catalog.name_index[entry.name], entry.uuid)
        end
       
        # Add new values to indexes
        insert!(catalog.temporal_index, entry.temporal_coverage[1], entry.temporal_coverage[2], entry.uuid)
        insert!(catalog.spatial_index, entry.spatial_bounds[1], entry.spatial_bounds[3], entry.spatial_bounds[2], entry.spatial_bounds[4], entry.uuid)
        for variable in entry.variables
            insert!(catalog.variable_index, variable, entry.uuid)
        end
    end
    
    # Auto-save if enabled
    if catalog.auto_save && !isempty(catalog.catalog_file)
        save_catalog(catalog)
    end
    
    println("Updated dataset: $(entry.name) (UUID: $(uuid))")
    return entry
end
function get_dataset(catalog::GeoDataCatalog, uuid::UUID)::DatasetEntry
    if !haskey(catalog.entries, uuid)
        throw(ArgumentError("Dataset entry not found: $uuid"))
    end
    return catalog.entries[uuid]
end

"""
    get_dataset_by_key(catalog::GeoDataCatalog, key::Symbol) -> DatasetEntry

    Get a dataset entry by its key.
    Throws ArgumentError if not found or multiple matches.
"""
function get_dataset_by_key(catalog::GeoDataCatalog, key::Symbol)::DatasetEntry
    if !haskey(catalog.key_index, key)
        throw(ArgumentError("No dataset found with key: :$key"))
    end
    
    uuid = catalog.key_index[key]
    if !haskey(catalog.entries, uuid)
        throw(ArgumentError("Inconsistent catalog: key index points to non-existent entry"))
    end
    
    return catalog.entries[uuid]
end

"""
    get_datasets_by_name(catalog::GeoDataCatalog, name::String)::Vector{DatasetEntry}

    Get all dataset entries matching a name (exact match).
    Returns empty vector if none found.
"""
function get_datasets_by_name(catalog::GeoDataCatalog, name::String)::Vector{DatasetEntry}
    if !haskey(catalog.name_index, name)
        return DatasetEntry[]
    end
    
    uuids = catalog.name_index[name]
    entries = DatasetEntry[]
    for uuid in uuids
        if haskey(catalog.entries, uuid)
            push!(entries, catalog.entries[uuid])
        end
    end
    return entries
end

"""
    find_datasets(catalog::GeoDataCatalog; 
                  name::Union{Nothing, String} = nothing,
                  key::Union{Nothing, Symbol} = nothing,
                  variables::Union{Nothing, Vector{String}} = nothing,
                  temporal_range::Union{Nothing, Tuple{DateTime, DateTime}} = nothing,
                  spatial_bounds::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}} = nothing,
                  vertical_range::Union{Nothing, Tuple{Float64, Float64}} = nothing,
                  format::Union{Nothing, Symbol} = nothing,
                  access_method::Union{Nothing, Symbol} = nothing,
                  active_only::Bool = true,
                  limit::Union{Nothing, Int} = nothing)::Vector{DatasetEntry}

    Find datasets matching various criteria.
    Returns vector of matching entries, sorted by priority (descending) then name.
"""
function find_datasets(catalog::GeoDataCatalog; 
                       name::Union{Nothing, String} = nothing,
                       key::Union{Nothing, Symbol} = nothing,
                       variables::Union{Nothing, Vector{String}} = nothing,
                       temporal_range::Union{Nothing, Tuple{DateTime, DateTime}} = nothing,
                       spatial_bounds::Union{Nothing, Tuple{Float64, Float64, Float64, Float64}} = nothing,
                       vertical_range::Union{Nothing, Tuple{Float64, Float64}} = nothing,
                       format::Union{Nothing, Symbol} = nothing,
                       access_method::Union{Nothing, Symbol} = nothing,
                       active_only::Bool = true,
                       limit::Union{Nothing, Int} = nothing)::Vector{DatasetEntry}
    
    matches = DatasetEntry[]
    
    for entry in values(catalog.entries)
        # Skip inactive entries if requested
        if active_only && !entry.active
            continue
        end
        
        # Check name
        if !(name === nothing) && entry.name != name
            continue
        end
        
        # Check key
        if !(key === nothing) && entry.key != key
            continue
        end
        
        # Check variables (if specified, all must be present)
        if !(variables === nothing)
            missing_vars = setdiff(variables, entry.variables)
            if !isempty(missing_vars)
                continue
            end
        end
        
        # Check temporal range (overlap required)
        if !(temporal_range === nothing)
            req_start, req_end = temporal_range
            entry_start, entry_end = entry.temporal_coverage
            if entry_end < req_start || req_end < entry_start
                continue  # No overlap
            end
        end
        
        # Check spatial bounds (overlap required)
        if !(spatial_bounds === nothing)
            req_min_lon, req_max_lon, req_min_lat, req_max_lat = spatial_bounds
            entry_min_lon, entry_max_lon, entry_min_lat, entry_max_lat = entry.spatial_bounds
            if entry_max_lon < req_min_lon || req_max_lon < entry_min_lon ||
               entry_max_lat < req_min_lat || req_max_lat < entry_min_lat
                continue  # No overlap
            end
        end
        
        # Check vertical range (overlap required)
        if !(vertical_range === nothing)
            req_min_depth, req_max_depth = vertical_range
            entry_min_depth, entry_max_depth = entry.vertical_range
            if entry_max_depth < req_min_depth || req_max_depth < entry_min_depth
                continue  # No overlap
            end
        end
        
        # Check format
        if !(format === nothing) && entry.format != format
            continue
        end
        
        # Check access method
        if !(access_method === nothing) && entry.access_method != access_method
            continue
        end
        
        # All checks passed - add to matches
        push!(matches, entry)
    end
    
    # Sort by priority (descending), then by name
    sort!(matches, by = x -> (-x.priority, lowercase(x.name)))
    
    # Apply limit if specified
    if !(limit === nothing)
        if length(matches) > limit
            resize!(matches, limit)
        end
    end
    
    return matches
end

"""
    get_best_dataset(catalog::GeoDataCatalog; kwargs...) -> DatasetEntry

    Get the single best matching dataset based on criteria.
    Throws ArgumentError if none found or multiple equally good matches.
    Returns the highest priority match.
"""
function get_best_dataset(catalog::GeoDataCatalog; kwargs...)::DatasetEntry
    matches = find_datasets(catalog; kwargs...)
    
    if isempty(matches)
        throw(ArgumentError("No datasets found matching criteria"))
    end
    
    # If there's a clear winner (highest priority), return it
    if length(matches) == 1 || matches[1].priority > matches[2].priority
        return matches[1]
    else
        # Multiple entries with same highest priority
        names = join([entry.name for entry in matches if entry.priority == matches[1].priority], ", ")
        throw(ArgumentError(
            "Multiple datasets found with equal highest priority: $names. " *
            "Please refine your search criteria or specify a key/name directly."))
    end
end

# ============================================================================
# Dataset Access and Usage Tracking
# ============================================================================

"""
    access_dataset!(catalog::GeoDataCatalog, uuid::UUID) -> DatasetEntry

    Record that a dataset has been accessed, updating usage statistics.
    Returns the dataset entry.
"""
function access_dataset!(catalog::GeoDataCatalog, uuid::UUID)::DatasetEntry
    if !haskey(catalog.entries, uuid)
        throw(ArgumentError("Dataset entry not found: $uuid"))
    end
    
    entry = catalog.entries[uuid]
    entry.access_count += 1
    entry.last_accessed = Dates.now()
    catalog.updated = Dates.now()
    catalog.stats_dirty = true
    
    # Auto-save if enabled
    if catalog.auto_save && !isempty(catalog.catalog_file)
        save_catalog(catalog)
    end
    
    return entry
end

"""
    record_download!(catalog::GeoDataCatalog, uuid::UUID, bytes_downloaded::UInt64) -> DatasetEntry

    Record that data was downloaded for a dataset, updating storage statistics.
    Returns the dataset entry.
"""
function record_download!(catalog::GeoDataCatalog, uuid::UUID, bytes_downloaded::UInt64)::DatasetEntry
    if !haskey(catalog.entries, uuid)
        throw(ArgumentError("Dataset entry not found: $uuid"))
    end
    
    entry = catalog.entries[uuid]
    entry.total_download_size += bytes_downloaded
    catalog.updated = Dates.now()
    catalog.stats_dirty = true
    
    # Auto-save if enabled
    if catalog.auto_save && !isempty(catalog.catalog_file)
        save_catalog(catalog)
    end
    
    return entry
end

# ============================================================================
# Catalog Import/Export Utilities
# ============================================================================

"""
    export_catalog_to_csv(catalog::GeoDataCatalog, file::String)

    Export catalog entries to a CSV file for external analysis.
"""
function export_catalog_to_csv(catalog::GeoDataCatalog, file::String)
    open(file, "w") do io
        # Write header
        header = [
            "uuid", "name", "key", "variables", "temporal_start", "temporal_end",
            "min_lon", "max_lon", "min_lat", "max_lat", "min_depth", "max_depth",
            "source", "location", "access_method", "format", "credits", "license",
            "registered", "updated", "version", "active", "priority", 
            "access_count", "last_accessed", "total_download_size"
        ]
        writedlm(io, [header], ',')
        
        # Write data rows
        for entry in values(catalog.entries)
            row = [
                string(entry.uuid),
                entry.name,
                string(entry.key),
                join(entry.variables, "|"),  # Pipe-separated list of variables
                string(entry.temporal_coverage[1]),
                string(entry.temporal_coverage[2]),
                entry.spatial_bounds[1], entry.spatial_bounds[2],  # min_lon, max_lon
                entry.spatial_bounds[3], entry.spatial_bounds[4],  # min_lat, max_lat
                entry.vertical_range[1], entry.vertical_range[2],  # min_depth, max_depth
                entry.source,
                entry.location,
                string(entry.access_method),
                string(entry.format),
                entry.credits,
                entry.license,
                string(entry.registered),
                string(entry.updated),
                entry.version,
                entry.active ? "true" : "false",
                entry.priority,
                entry.access_count,
                string(entry.last_accessed),
                entry.total_download_size
            ]
            writedlm(io, [row], ',')
        end
    end
    
    println("Exported $(length(catalog.entries)) dataset entries to CSV: $file")
    return nothing
end

# ============================================================================
# Convenience Functions and Default Catalog
# ============================================================================

"""
    default_catalog_path() -> String

    Get the default path for the GeoData catalog file.
"""
function default_catalog_path()
    homedir() |> joinpath -> ".geodata" |> joinpath -> "catalog.json"
end

"""
    GlobalCatalog

    A singleton global catalog instance for convenient access.
"""
const GlobalCatalog = Ref{GeoDataCatalog}(GeoDataCatalog(default_catalog_path()))

"""
    get_global_catalog() -> GeoDataCatalog

    Get the global catalog instance.
"""
function get_global_catalog()
    return GlobalCatalog[]
end

"""
    init_global_catalog!(catalog_file::String=""; kwargs...)

    Initialize the global catalog with optional file and parameters.
"""
function init_global_catalog!(catalog_file::String=""; kwargs...)
    GlobalCatalog[] = GeoDataCatalog(isempty(catalog_file) ? default_catalog_path() : catalog_file; kwargs...)
    return GlobalCatalog[]
end

# ============================================================================
# Export Public Interface
# ============================================================================

export GeoDataCatalog,
       DatasetEntry,
       CatalogStats,
       
       # Construction and persistence
       load_catalog!,
       save_catalog,
       
       # Catalog modification
       register_dataset!,
       unregister_dataset!,
       update_dataset!,
       
       # Query and discovery
       get_dataset,
       get_dataset_by_key,
       get_datasets_by_name,
       find_datasets,
       get_best_dataset,
       
       # Usage tracking
       access_dataset!,
       record_download!,
       
       # Statistics
       get_stats,
       
       # Import/export
       export_catalog_to_csv,
       
       # Global catalog
       get_global_catalog,
       init_global_catalog!

end # module GeoDataCatalog
