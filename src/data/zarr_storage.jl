"""
    zarr_storage.jl

General Zarr storage utilities for geospatial data.
"""
module ZarrStorage

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using Zarr
using JSON3
using Dates
using Statistics
using LinearAlgebra

# General Zarr storage exports
export
    open_zarr_storage,
    close_zarr_storage,
    create_zarr_array,
    write_zarr_array,
    read_zarr_array,
    zarr_create_group,
    zarr_open_group,
    zarr_append

"""
    open_zarr_storage(
        store_path::AbstractString;
        mode::Symbol = :rw
    ) -> Zarr.ZGroup

Open or create a Zarr storage group.
"""
function open_zarr_storage(
    store_path::AbstractString;
    mode::Symbol = :rw
)
    mkpath(store_path)
    store = Zarr.DirectoryStore(store_path)
    return Zarr.open_group(store; mode = mode)
end

"""
    close_zarr_storage(group::Zarr.ZGroup)

Close the Zarr storage group.
"""
function close_zarr_storage(group::Zarr.ZGroup)
    Zarr.close(group.store)
    nothing
end

"""
    create_zarr_array(
        group::Zarr.ZGroup,
        name::AbstractString;
        shape,
        chunks,
        dtype,
        compressor = Zarr.Blosc(cname = "zstd", clevel = 3)
    )

Create a Zarr array with compression.
"""
function create_zarr_array(
    group::Zarr.ZGroup,
    name::AbstractString;
    shape,
    chunks,
    dtype,
    compressor = Zarr.Blosc(cname = "zstd", clevel = 3)
)
    if !haskey(group, name)
        Zarr.create_array(
            group, name;
            shape = shape,
            chunks = chunks,
            dtype = dtype,
            compressor = compressor
        )
    end
    return group[name]
end

"""
    write_zarr_array(array::Zarr.ZArray, data)

Write data to a Zarr array.
"""
function write_zarr_array(array::Zarr.ZArray, data)
    Zarr.write(array, data)
end

"""
    read_zarr_array(array::Zarr.ZArray)

Read data from a Zarr array.
"""
function read_zarr_array(array::Zarr.ZArray)
    Zarr.read(array)
end

"""
    zarr_create_group(parent::Zarr.ZGroup, name::AbstractString)

Create a Zarr group.
"""
function zarr_create_group(parent::Zarr.ZGroup, name::AbstractString)
    if !haskey(parent, name)
        Zarr.create_group(parent, name)
    end
    return parent[name]
end

"""
    zarr_open_group(store_path::AbstractString; mode::Symbol = :rw)

Open a Zarr group directly from path.
"""
function zarr_open_group(store_path::AbstractString; mode::Symbol = :rw)
    mkpath(store_path)
    store = Zarr.DirectoryStore(store_path)
    Zarr.open_group(store; mode = mode)
end

"""
    zarr_append(array::Zarr.ZArray, new_data)

Append data to a Zarr array along the first dimension.
Note: Zarr doesn't support true append, so this reads, concatenates, and rewrites.
"""
function zarr_append(array::Zarr.ZArray, new_data)
    existing = Zarr.read(array)
    combined = vcat(existing, new_data)
    # Resize array
    Zarr.write(array, combined)
end

end # module ZarrStorage