"""
    storage.jl

Hierarchical analytical and trajectory storage abstraction for GeoData.

Provides a unified interface (`GeoStorage`) for managing multi-group, multi-dataset
scientific archives (such as simulation runs, particle trajectories, metrics,
and demographic connectivity matrices) across both chunked tensor arrays (Zarr)
and tabular partitioned collections (GeoParquet).
"""

using ..GeoDataCoreTypes
using DataFrames
using Dates
import Zarr
import GeoParquet

export GeoStorage,
       open_geostorage,
       close_geostorage,
       create_storage_group,
       has_storage_group,
       write_storage_variable!,
       read_storage_variable,
       list_storage_variables,
       list_storage_groups

"""
    GeoStorage

Abstract container representing a hierarchical analytical storage archive.
Supports both Zarr multi-group stores and GeoParquet directory collections.
"""
mutable struct GeoStorage
    uri::String
    backend::Symbol
    read_only::Bool
    handle::Any

    function GeoStorage(
        uri::AbstractString;
        backend::Union{Symbol, Nothing} = nothing,
        read_only::Bool = false,
        create::Bool = true
    )
        resolved_be = if !isnothing(backend)
            backend
        elseif endswith(lowercase(uri), ".parquet")
            :geoparquet
        elseif endswith(lowercase(uri), ".zarr")
            :zarr
        else
            :zarr
        end

        norm_path = if resolved_be == :zarr
            endswith(uri, ".zarr") ? String(uri) : String(uri) * ".zarr"
        elseif resolved_be == :geoparquet
            endswith(uri, ".parquet") ? String(uri) : String(uri) * ".parquet"
        else
            String(uri)
        end

        h = if resolved_be == :zarr
            _init_zarr_storage(norm_path, read_only, create)
        elseif resolved_be == :geoparquet
            _init_parquet_storage(norm_path, read_only, create)
        else
            error("Unsupported storage backend: $(resolved_be). Supported: :zarr, :geoparquet")
        end

        return new(norm_path, resolved_be, read_only, h)
    end
end

function _init_zarr_storage(path::String, read_only::Bool, create::Bool)
    if read_only
        return Zarr.zopen(path, "r")
    end
    if isdir(path) && isfile(joinpath(path, ".zgroup"))
        return Zarr.zopen(path, "w")
    else
        create && mkpath(path)
        return Zarr.zgroup(path)
    end
end

mutable struct ParquetDirectoryHandle
    path::String
    read_only::Bool
    cached_tables::Dict{String, Any}
end

function _init_parquet_storage(path::String, read_only::Bool, create::Bool)
    base_dir = if isdir(path)
        path
    elseif endswith(path, ".parquet")
        d = dirname(path)
        isempty(d) ? "." : d
    else
        path
    end
    if create && !read_only
        mkpath(base_dir)
    end
    return ParquetDirectoryHandle(path, read_only, Dict{String, Any}())
end

"""
    open_geostorage(uri::AbstractString; backend=:auto, read_only=false, create=true) -> GeoStorage

Open or initialize a hierarchical storage container at `uri`.
"""
function open_geostorage(
    uri::AbstractString;
    backend::Union{Symbol, Nothing} = nothing,
    read_only::Bool = false,
    create::Bool = true
)::GeoStorage
    return GeoStorage(uri; backend = backend, read_only = read_only, create = create)
end

"""
    close_geostorage(storage::GeoStorage) -> Nothing

Flush buffered data and release file locks.
"""
function close_geostorage(storage::GeoStorage)
    GC.gc()
    return nothing
end

"""
    create_storage_group(storage::GeoStorage, group_path::AbstractString)

Create a nested group hierarchy in the storage archive.
"""
function create_storage_group(storage::GeoStorage, group_path::AbstractString)
    if storage.backend == :zarr
        grp = storage.handle
        parts = split(replace(group_path, "\\" => "/"), "/")
        cur = grp
        for p in parts
            s = String(p)
            isempty(s) && continue
            cur = haskey(cur, s) ? cur[s] : Zarr.zgroup(cur, s)
        end
        return cur
    elseif storage.backend == :geoparquet
        target_dir = joinpath(dirname(storage.uri), group_path)
        mkpath(target_dir)
        return target_dir
    end
end

"""
    has_storage_group(storage::GeoStorage, group_path::AbstractString) -> Bool
"""
function has_storage_group(storage::GeoStorage, group_path::AbstractString)::Bool
    if storage.backend == :zarr
        parts = split(replace(group_path, "\\" => "/"), "/")
        cur = storage.handle
        for p in parts
            s = String(p)
            isempty(s) && continue
            !haskey(cur, s) && return false
            cur = cur[s]
        end
        return true
    elseif storage.backend == :geoparquet
        target_dir = joinpath(dirname(storage.uri), group_path)
        return isdir(target_dir)
    end
end

"""
    write_storage_variable!(storage::GeoStorage, var_path::AbstractString, data; kwargs...)

Write an array, vector, or table into `storage` at `var_path`.
"""
function write_storage_variable!(
    storage::GeoStorage,
    var_path::AbstractString,
    data;
    chunks = nothing,
    compressor = nothing
)
    storage.read_only && error("Cannot write to read-only GeoStorage at $(storage.uri)")

    if storage.backend == :zarr
        norm = replace(var_path, "\\" => "/")
        parts = split(norm, "/")
        var_name = String(last(parts))
        parent_group = storage.handle
        if length(parts) > 1
            for p in parts[1:end-1]
                s = String(p)
                isempty(s) && continue
                parent_group = haskey(parent_group, s) ?
                    parent_group[s] : Zarr.zgroup(parent_group, s)
            end
        end

        if haskey(parent_group, var_name)
            parent_group[var_name][:] = data
        else
            c_size = isnothing(chunks) ? size(data) : chunks
            kw = Dict{Symbol, Any}(:chunks => c_size)
            if !isnothing(compressor)
                kw[:compressor] = compressor
            end
            arr = Zarr.zcreate(
                eltype(data), parent_group, var_name, size(data)...;
                kw...
            )
            arr[:] = data
        end
        return nothing
    elseif storage.backend == :geoparquet
        target_file = if endswith(var_path, ".parquet")
            joinpath(dirname(storage.uri), var_path)
        else
            joinpath(dirname(storage.uri), var_path * ".parquet")
        end
        mkpath(dirname(target_file))
        if data isa DataFrame
            GeoParquet.write(target_file, data)
        else
            df = DataFrame(value = vec(data))
            GeoParquet.write(target_file, df)
        end
        return nothing
    end
end

"""
    read_storage_variable(storage::GeoStorage, var_path::AbstractString)

Read an array or tabular dataset from `storage` at `var_path`.
"""
function read_storage_variable(storage::GeoStorage, var_path::AbstractString)
    if storage.backend == :zarr
        parts = split(replace(var_path, "\\" => "/"), "/")
        cur = storage.handle
        for p in parts
            s = String(p)
            isempty(s) && continue
            cur = cur[s]
        end
        return cur[:]
    elseif storage.backend == :geoparquet
        target_file = if endswith(var_path, ".parquet")
            joinpath(dirname(storage.uri), var_path)
        else
            joinpath(dirname(storage.uri), var_path * ".parquet")
        end
        isfile(target_file) || error("Variable file not found at $(target_file)")
        return GeoParquet.read(target_file)
    end
end
