# GeoData Reference

`GeoData.jl` provides backend-agnostic geospatial data access, manipulation, and
storage abstraction across multiple scientific file formats (Zarr, NetCDF/NCZarr,
GeoParquet).

It holds **containers, backends, coordinate operations, and the API**. It does not
contain the fetching layer, the source registry, the catalog, or the data-cube engine:
those are in [`GeoDataSources.jl`](GeoDataSources.md), which depends on this package and
never the other way round.

---

## 1. Core Architecture and Types

### Core types (`GeoData.GeoDataCoreTypes`)

- **`GeoDataset`** — a multi-variable geospatial container:
  - `variables::Dict{String, GeoArray}` — data variables, keyed by **name**.
  - `coords::Dict{Symbol, GeoArray}` — coordinate arrays (`:lon`, `:lat`, `:depth`, `:time`).
  - `dims::Dict{Symbol, Dimension}` — dimension records.
  - `crs::CoordinateSystem` — coordinate reference system (e.g. `"EPSG:4326"`).
  - `attrs::Dict{String, Any}` — global attributes. Any `Dict` is accepted on
    construction and normalised to `Dict{String, Any}`.
  - `backend::GeoBackend` — the backend the dataset was loaded from, or `nothing`.
  - `source::String` — URI or path the dataset came from. (This field is `source`, not
    `uri`; a common source of confusion.)
- **`GeoArray{T, N}`** — an N-dimensional array with dimension tags, CRS, and attributes.
  It forwards `size`, `ndims`, `eltype`, `axes`, `getindex`, `setindex!`, `IndexStyle`,
  and `similar`.
- **`Dimension`** — `name`, `size`, `coords`, `units`, `standard_name`, `calendar`,
  `dim_type`, `is_unlimited`. `units` is `Union{String, Nothing}`: `nothing` means
  unknown, which is different from dimensionless (`"1"`). `Dimension(existing; kwargs...)`
  copies a dimension with fields replaced, which every narrowing operation needs.
- **`CoordinateSystem`** — geodetic or Cartesian CRS descriptor (`crs::String`).

### The dimension selection  

Every keyword selection on `dim` follows one rule, enforced by
`assert_dataset_invariants` (in the test helpers):

| Keyword form | Meaning |
|:---|:---|
| `dim = (lo, hi)` | the index range spanning the two coordinates; bounds snap inward |
| `dim = [v1, v2, ...]` | the nearest index to each value, in that order (duplicates kept) |
| `dim = v` | the nearest index, and the dimension is **dropped** from the result's variables, `dims`, *and* `coords`; the value is recorded in `attrs["selections"]` |

A dimension not named is returned whole. A scalar selection reduces the rank everywhere
consistently, so a sliced dataset never claims a size its variables do not have.

---

## 2. High-level dataset operations

```julia
using GeoData

# Load: the backend is inferred from the path or URI scheme
ds = geoload("data/bathymetry.zarr")

# Coordinate-aware slicing
subset = geoslice(ds, lon = (-68.0, -57.0), lat = (42.0, 47.5))

# Point value (scalar when every selected dimension is dropped)
val = geovalues(ds, ["temp"], lon = -63.5, lat = 44.0, depth = 10.0)["temp"]

# Reduce along a dimension
monthly = geoaggregate(ds, dim = :time, func = mean)

# Join along shared coordinates: missing values are NaN, never nearest-neighbour
merged = geojoin([ds_a, ds_b], on = [:lon, :lat], how = :left)
```

The `geo*` entry points are in `GeoData` and dispatch through `backend_supports`, which
checks whether the loaded backend implements the corresponding `backend_*` method. No
backend currently implements one, so every call uses the generic path; adding a backend
method can make an operation faster but can never change its result.

---

## 3. Storage backends

Each backend implements `GeoDataInterfaces.backend_open`, `backend_create`,
`backend_write`, and `backend_close`, and declares its own
`BackendCapabilities`.

| Backend | Identifier | URI schemes / extensions | Notes |
|:---|:---|:---|:---|
| Zarr | `:zarr` | `zarr://`, `.zarr` | Chunked, compressed. `backend_open` returns the chunked arrays themselves, so opening reads nothing and indexing materialises only the elements it touches |
| NetCDF | `:ncdatasets` | `netcdf://`, `.nc` | Axes are self-describing. Variables are materialised on read |
| NCZarr | `:nczarr` | `nczarr://` | NetCDF-Zarr mapping through NCDatasets |
| GeoParquet | `:geoparquet` | `geoparquet://`, `.parquet` | One `:points` dimension of length `nrow`; `:lon`/`:lat` coordinates when the geometry column is point-like; the geometry column is kept as a variable so a round trip returns the same file |

Backend detection is explicit:

```julia
infer_backend("a.zarr")        # :zarr
infer_backend("nczarr://store") # :nczarr
infer_backend("a.tif")          # ERROR: names the supported schemes and the `backend=` keyword
```

A source URI's `scheme://` prefix is removed by `strip_scheme` using URI parsing.

### Writing

`geosave(uri, ds; backend, overwrite = false)` resolves a backend and delegates to
`backend_write`.

- NetCDF updates an existing file in place (mode `"a"`), adding or replacing variables.
- Zarr datasets cannot be modified in place when reopened as read-only groups:
  pass `overwrite = true` to recreate the store.
- `backend_create` validates that `dims` covers every axis of every variable;
  omitted dimensions raise an `ArgumentError` or `ErrorException`.

---

## 4. Hierarchical analytical storage (`GeoStorage`)

For simulation archives and multi-group experiments, `GeoStorage` manages hierarchical
Zarr or GeoParquet group directories:

```julia
storage = open_geostorage("work/run/archive.zarr")
create_storage_group(storage, "outputs")
write_storage_variable!(storage, "outputs", "temperature", arr)
arr2 = read_storage_variable(storage, "outputs", "temperature")
close_geostorage(storage)
```

`open_geostorage` on an existing store opens it for inspection without truncation;
pass `create = true` or `overwrite = true` to initialize or replace storage.

---

## 5. Dataset Invariants

Every dataset produced in testing is validated with `assert_dataset_invariants`:

1. Every variable axis length matches its corresponding `Dimension` record.
2. Every dimension spanned by a variable is registered in `ds.dims`.
3. Every declared dimension is referenced by a variable or defined as a coordinate.
4. Every coordinate array is 1D with a matching dimension record.

## Related

- [`GeoDataSources.md`](../GeoDataSources/README.md) — fetching, provenance, the source registry,
  the catalog, and the data-cube engine. 