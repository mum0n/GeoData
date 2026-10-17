### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ 6b4da31c-2e9f-4a05-9c43-9d2f7a0e3c66
# ╠═╡
using Pkg

# ╔═╡ 8d6fc53e-401b-4c27-9b65-bef4e5f10a88
# ╠═╡
Pkg.activate(@__DIR__)

# ╔═╡ 4f2b8a19-0c7d-4e63-9a21-7b0d5e8c1a44
# ╠═╡
md"""
# GeoData — the backend-agnostic geospatial layer

How to read, write, slice, query, join, regrid, and land-mask geospatial data without caring
which file format it lives in — and what each option does, because in this package an unclear
option is a silent wrong answer.

Everything here is offline: datasets are built in memory, written to a scratch directory, and
read back, so a cell can be re-run while you inspect what it returned.

### The contract in one paragraph

`GeoDataset` is a container of `GeoArray`s, each carrying its `Dimension` records, a
`CoordinateSystem` (CRS), and units. Four invariants define a *self-consistent* dataset: every
variable's shape matches its dimension records; every dimension a variable spans is declared in
`ds.dims`; every declared dimension is used by a variable or exists as a coordinate; and every
coordinate is 1-D with a matching size. A dataset that fails them is a bug, not a style problem —
which is why every cell below ends by printing a summary rather than trusting that a call worked.
"""

# ╔═╡ 5a3c920b-1d8e-4f74-8b32-8c1e6f9d2b55
# ╠═╡
md"""
## Running this notebook

Launch Pluto from an environment that already has it, and open this file:

```
julia -e 'using Pluto; Pluto.run()'
```

The first three cells pin the working directory to the package root and activate **this
package's environment**, so `scratch/` paths stay inside the project and `using GeoData` loads
the working tree.

One consequence of the format is worth knowing, because it is what makes this notebook useful: a
cell defines something and every later cell sees it, and re-running a cell re-runs everything
downstream of it. That is what makes "change the slice and watch the answer move" cheap.
"""

# ╔═╡ 7c5eb42d-3f0a-4b16-8a54-ae3a8b1f4d77
# ╠═╡
cd(joinpath(@__DIR__, ".."))

# ╔═╡ d213b7d0-ca2e-11f1-84d9-e7b99ec2f5a2
# ╠═╡
using GeoData
using GeoData: GeoDataset, GeoArray, Dimension, CoordinateSystem
using GeoData.GeoDataTypes: DIM_SPATIAL, DIM_PARAMETRIC, DIM_TEMPORAL

# ╔═╡ af81e750-623d-4e49-9a87-d0b6b7130cfa
# ╠═╡
md"""
## The types, and what they are for

- **`GeoDataset`** — variables + coordinates + dimension records + CRS + attributes. It is the
  unit everything else takes and returns.
- **`GeoArray`** — one variable's data with the dimensions it spans.
- **`Dimension`** — a name, a size, and optionally a coordinate array (`coords=nothing` for a
  dimension that has no coordinate values, which is a statement, not an omission).
- **``CoordinateSystem(crs = "EPSG:4326")`** — longitude/latitude degrees. One CRS is in use
  here; the field exists so that a host that says otherwise cannot be assumed away.

The constructor takes explicit dimension records because a dataset built from a bare array has
no idea what its axes mean, and every later operation is a coordinate operation.

The cell below builds a 3-D field (lon × lat × depth) with a separable analytic value field, so
every later cell has something real to work on.
"""

# ╔═╡ bf92f86a-734e-4f5a-9b98-e1c7c8241d0b
# ╠═╡
# The one constructor helper this notebook needs. It is written out in full rather than imported
# from the test suite, so the construction of a dataset is legible where it is used.
function make_2d(varname::AbstractString, data::AbstractMatrix, lon, lat)
    crs = CoordinateSystem(crs = "EPSG:4326")
    dlon = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                     dim_type = DIM_SPATIAL)
    dlat = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                     dim_type = DIM_SPATIAL)
    return GeoDataset(
        Dict(varname => GeoArray(data, (dlon, dlat), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dlon,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dlat,), crs, Dict{String, Any}())),
        Dict(:lon => dlon, :lat => dlat), crs, Dict{String, Any}(), nothing, "synthetic")
end

function make_3d(varname::AbstractString, data::AbstractArray, lon, lat, depth)
    crs = CoordinateSystem(crs = "EPSG:4326")
    dlon = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                     dim_type = DIM_SPATIAL)
    dlat = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                     dim_type = DIM_SPATIAL)
    ddep = Dimension(name = :depth, size = length(depth), coords = depth, units = "m",
                     dim_type = DIM_PARAMETRIC)
    return GeoDataset(
        Dict(varname => GeoArray(data, (dlon, dlat, ddep), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dlon,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dlat,), crs, Dict{String, Any}()),
             :depth => GeoArray(depth, (ddep,), crs, Dict{String, Any}())),
        Dict(:lon => dlon, :lat => dlat, :depth => ddep), crs, Dict{String, Any}(),
        nothing, "synthetic")
end

const LON = collect(-65.0:0.5:-62.0)
const LAT = collect(43.0:0.5:45.0)
const DEPTH = collect(0.0:25.0:100.0)

# Separable in each axis, so a wrong axis order shows up in the value rather than in the shape.
function _surface(lon, lat)
    return 10.0 + 0.2 * (lon - minimum(LON)) - 0.1 * (lat - maximum(LAT))
end

const ds = make_3d("temperature",
                   [_surface(lon, lat) - 0.01 * dep for lon in LON, lat in LAT, dep in DEPTH],
                   LON, LAT, DEPTH)

println(geosummary(ds))
println("bounding box: ", bounding_box(ds))
println("temperature stats: ", variable_stats(ds, "temperature"))

# ╔═╡ c0a39b14-846c-4f6b-9cb9-f2d8d9352e1c
# ╠═╡
md"""
## Reading and writing: the backend is inferred, and it can be named

`geoload(uri; backend=nothing)` picks a backend from the URI when you do not name one, and
`geosave(uri, ds; overwrite=false)` does the same in the other direction. Two options matter:

- **`backend`** — `:zarr`, `:ncdatasets`, `:nczarr`, `:geoparquet`. Name it when the extension
  is ambiguous or absent. `infer_backend(uri)` shows what would be picked.
- **`overwrite`** — the one place the layer is *not* uniform, so it is worth reading closely.
  `overwrite=false` on a destination that is already there **replaces the variables it contains
  in place** for NetCDF and GeoParquet, and is an **error** for Zarr — Zarr.jl reopens an
  existing store read-only and refuses to create an array in a reopened group, so the backend
  refuses rather than leaving a half-written store. The cell below shows both behaviours, because
  a notebook that hid the second one would be teaching a contract the Zarr backend does not keep.

What differs between backends is *when* the data are read, not what you can do with them:

| backend | format | strength | materialised |
|---|---|---|---|
| `:zarr` | chunked store | cloud-native, lazy | on index |
| `:ncdatasets` | NetCDF | CF conventions, interop | on open |
| `:geoparquet` | GeoParquet | vector/tabular, Arrow ecosystem | on open |

The GeoParquet backend is for *tables*: one dimension, `:points`, of length `nrow(df)`. It
extracts `:lon`/`:lat` coordinate arrays when the geometry column is point-like, and for line
and polygon geometries leaves them undefined rather than faking a grid — a coastline is not a
grid, and pretending otherwise is how a figure ends up with an invented axis.
"""

# ╔═╡ d1b4c225-957d-4b6a-9dcb-a13e9f462f0d
# ╠═╡
const SCRATCH = normpath(joinpath("scratch", "notebook_geodata"))
isdir(SCRATCH) || mkpath(SCRATCH)

println("inferred backend for .zarr   : ", infer_backend(joinpath(SCRATCH, "field.zarr")))
println("inferred backend for .nc     : ", infer_backend(joinpath(SCRATCH, "field.nc")))
println("inferred backend for .parquet: ", infer_backend(joinpath(SCRATCH, "field.parquet")))

zarr_path = joinpath(SCRATCH, "temperature.zarr")
geosave(zarr_path, ds; overwrite = true)
back = geoload(zarr_path)
println("\nround trip through Zarr   : ", size(back.variables["temperature"].data),
        "  max diff ",
        maximum(abs, back.variables["temperature"].data .- ds.variables["temperature"].data))

nc_path = joinpath(SCRATCH, "temperature.nc")
geosave(nc_path, ds; backend = :ncdatasets, overwrite = true)
back_nc = geoload(nc_path; backend = :ncdatasets)
println("round trip through NetCDF : ", size(back_nc.variables["temperature"].data),
        "  max diff ",
        maximum(abs, back_nc.variables["temperature"].data .- ds.variables["temperature"].data))

# The same save twice, at the default `overwrite = false`. NetCDF updates the file in place;
# Zarr refuses, and says why. Both are the documented behaviour, not an accident.
geosave(nc_path, ds; backend = :ncdatasets)
println("NetCDF written twice      : still ", size(geoload(nc_path).variables["temperature"].data))

try
    geosave(zarr_path, ds)          # no overwrite
catch err
    println("Zarr written twice       : ", first(split(sprint(showerror, err), "\n")))
end

# ╔═╡ e2c5d336-068e-4d7a-9aec-b2fa057f3a1b
# ╠═╡
md"""
## Slicing: one contract, three spellings

`geoslice` and `geoselect` share a single contract (both end at the same operation; the two
names exist because callers read differently):

| value passed | meaning |
|---|---|
| `lon = (-66.0, -63.0)` | a **range**, given as a `(min, max)` pair |
| `lon = [-65.0, -64.0]` | a **vector** of coordinate values, in that order, duplicates allowed |
| `lon = -64.5` | a **single value**, which selects the point *and drops the dimension* |

That last row is the one that surprises people: selecting one length along an axis removes it,
so a 3-D query collapses to 2-D. It is the behaviour you want for a map, and the reason a
"missing" dimension error later usually means a single value was passed earlier.

Three helpers sit on top:

- **`geoslice_bbox(ds; lon_range, lat_range, ...)`** — a bounding box in degrees. A range is
  applied only when the dataset has that dimension, so a 2-D bathymetry set and a 4-D
  temperature set take the same keyword list. Ranges outside the extent **snap inward**: the
  result is never larger than the data.
- **`geosubset_ids(ds, var, coord, ids)`** — explicit integer indices in an order you choose,
  translated to coordinate values so later operations keep working.
- **`geosubset_cube`** — the multi-variable box.

`backend` is accepted and ignored by the coordinate operations: slicing is generic and runs the
same way for every backend. It is kept so callers do not have to special-case a dataset they
happened to load from Zarr.
"""

# ╔═╡ f3d6e447-179f-4b8b-9afd-c3ab168c4b2c
# ╠═╡
function _dims(ds_local)
    return join([string(k, "=", v.size) for (k, v) in
                 sort(collect(ds_local.dims), by = first)], "  ")
end

ranged = geoslice(ds; lon = (-65.0, -63.5), lat = (43.5, 44.5))
println("range slice        : ", _dims(ranged))

picked = geoslice(ds; lon = [-65.0, -64.5, -63.0])
println("vector slice lons  : ", picked.coords[:lon].data)

point = geoslice(ds; lon = -64.5, lat = 44.0)
println("single value       : dims ", _dims(point),
        "  value at depth 0 m ", geovalue(point, "temperature"; depth = 0.0))
# Note the query above can no longer name `lon` or `lat`: that is the collapse from row 3 of the
# table, seen from the consumer side. A later "no coordinate for :lon" error usually means a
# single value was selected earlier, not that the dataset is broken.

bbox = geoslice_bbox(ds; lon_range = (-70.0, -60.0), lat_range = (44.0, 50.0),
                     depth_range = (0.0, 10.0))
println("bbox (snapped)     : lon ", extrema(bbox.coords[:lon].data),
        "  lat ", extrema(bbox.coords[:lat].data), "  depth ", extrema(bbox.coords[:depth].data))

reordered = geosubset_ids(ds, "temperature", :lat, [3, 1])
println("subset ids [3, 1]  : ", reordered.coords[:lat].data)

# ╔═╡ 04a7f558-2a90-4c9b-8b0f-d4bc279d5c3d
# ╠═╡
md"""
## Querying values at points

- **`geovalues(ds, varnames; lon, lat, depth, time)`** — returns a *dict* of variable name to
  value. Pass scalars for one point, or vectors for many; vectors are paired positionally, so
  the first lon is matched with the first lat. A mismatched pair is an error, not a ragged
  result.
- **`geovalue(ds, var; ...)`** — one variable, one value.
- **`geoprofile(ds, vars; lon, lat, time)`** — a vertical profile, returning the depths it used
  alongside the values.
- **`geotimeseries(ds, vars; lon, lat, depth)`** — the time axis of the dataset at one location.

The keys you pass must name a dimension the dataset has. An axis a dataset does not carry is a
fact about the dataset, reported rather than invented — no default depth or time is assumed.
"""

# ╔═╡ 159b8604-3ca1-4dac-9d0f-e5cd390ef2d9
# ╠═╡
one_point = geovalues(ds, ["temperature"]; lon = -64.5, lat = 44.0, depth = 25.0)
println("one point  : ", one_point["temperature"])

many = geovalues(ds, ["temperature"]; lon = [-64.5, -64.0, -63.5], lat = [44.0, 44.5, 43.5],
                 depth = 0.0)
println("three pts  : ", many["temperature"])

prof = geoprofile(ds, ["temperature"]; lon = -64.5, lat = 44.0)
println("profile    : depths ", prof[:depth])
println("             values ", prof["temperature"])

# ╔═╡ 26bc9715-4b2e-4adb-8e10-f6de4b1fa2b3
# ╠═╡
md"""
## Reducing and joining

**`geoaggregate(ds; dim, func)`** reduces every variable along `dim` and drops that dimension.
`func` must accept the array-positional form `func(data; dims = i)`, which is why `mean`, `sum`,
`maximum` and `extrema` all work while a closure taking a matrix does not. A variable that does
not span `dim` comes back unchanged rather than being dropped, and the reduction is recorded in
`attrs["reduced"]` so the result says what happened to it.

**`geojoin(datasets, on = [...], how = :inner)`** joins along shared dimensions. The `how` choices
are the relational ones — `:inner` (intersection), `:outer` (union), `:left`, `:right` — and the
join key is a dimension name list, because that is the axis two datasets can genuinely be said
to share. A variable name present in **both** datasets is an error, not an automatic rename: two
products contributing the same name is a decision the caller has to make.

**`geomerged(datasets)`** is the narrow case: the same grid, and you want the variables together.
It renames collisions by suffixing the dataset index, and it refuses when the grid counts differ
rather than resampling for you.
"""

# ╔═╡ 37cad826-5c3e-4ebb-9f21-070f5c2ab3c4
# ╠═╡
using Statistics

surf = geoaggregate(ds; dim = :depth, func = mean)
println("mean over depth : dims ", sort(collect(keys(surf.dims))),
        "  attrs ", get(surf.attrs, "reduced", "(none)"))

zmax = geoaggregate(ds; dim = :depth, func = maximum)
println("max  over depth : ", extrema(zmax.variables["temperature"].data))

# Two products on the same grid, which is what a join is for.
salinity = make_2d("salinity",
                   [30.0 + 0.05 * (lat - minimum(LAT)) for lon in LON, lat in LAT], LON, LAT)
joined = geojoin([geoslice(ds; depth = 0.0), salinity], on = [:lon, :lat])
println("join on [:lon, :lat] : vars ", sort(collect(keys(joined.variables))),
        "  dims ", sort(collect(keys(joined.dims))))

# The same variable name in both datasets is refused, not silently merged.
try
    geojoin([geoslice(ds; depth = 0.0), geoslice(ds; depth = 25.0)], on = [:lon, :lat])
catch err
    println("shared name refused   : ", first(split(sprint(showerror, err), "\n")))
end

merged = geomerged([geoslice(ds; depth = 0.0), salinity])
println("geomerged            : vars ", sort(collect(keys(merged.variables))))

# ╔═╡ 48dbe937-6d4c-4fcc-9a32-181f6d3bc4d5
# ╠═╡
md"""
## Regridding: source grid → target grid

`regrid_2d_field` and `regrid_3d_field` interpolate a field from the grid it is on onto another.
They accept the pieces separately (`lon, lat, field, tgt_lon, tgt_lat`) or a dataset with a
variable name, and they return the array only — a regridded field is not yet a dataset, because
what it is depends on the grid it landed on.

Two properties worth knowing before you trust the numbers:

- the target grid may be given in any order; a **descending source latitude is sorted, not
  mirrored**. A silently flipped axis produces a plausible-looking field with north and south
  exchanged, which is the hardest wrong answer to notice.
- the output shape is `(length(tgt_lon), length(tgt_lat))` (and `..., length(tgt_depth)` for
  the 3-D form), so the axis order follows the arguments you passed.
"""

# ╔═╡ 59ecf048-7e5d-4fdd-9b3c-2924088cd4d6
# ╠═╡
tgt_lon = collect(-64.8:0.2:-63.0)
tgt_lat = collect(44.0:0.25:44.75)

coarse = regrid_2d_field(LON, LAT, _spatial(LON, LAT), tgt_lon, tgt_lat)
println("regrid_2d_field, lon/lat/field -> size ", size(coarse))

fine_lat = regrid_2d_field(LON, reverse(LAT), [_spatial(lon, lat) for lon in LON, lat in reverse(LAT)],
                           tgt_lon, tgt_lat)
println("descending lat handled        -> size ", size(fine_lat),
        "  (not transposed: ", size(fine_lat) == size(coarse), ")")

from_ds = regrid_2d_field(geoslice(ds; depth = 0.0), "temperature", tgt_lon, tgt_lat)
println("regrid_2d_field, ds + var     -> size ", size(from_ds))

# ╔═╡ 6afd1597-8a0e-4fde-9c54-3a35179ce5e7
# ╠═╡
md"""
## Land and water: rings, points, and marine cells

The land test is deliberately generic. `GeoData` holds **no bundled geography**: a ring set is
data, and a set that covers one region would make every consumer of this generic layer carry
that region.

- **`point_in_polygon(lon, lat, lons, lats)`** — ray casting over vertex vectors, determined for
  points on the boundary. Edge cases are the test's job, not an epsilon's: a point on the edge
  is answered by a rule, not by a tolerance.
- **`load_coastline_polygons(file)`** — a GeoParquet coastline to rings. Each ring is
  `(name, code, lons, lats)`, and every `LineString` feature and every part of a
  `MultiLineString` becomes one ring. An absent file returns no rings, which the land test
  treats as "unknown, not land".
- **`is_point_on_land_geodata(lon, lat; coastline, coastline_file)`** — rings win over file;
  no rings means not land. An **open** ring (first point ≠ last) is skipped, because an open
  line cannot enclose anything.
- **`is_marine_water_geodata(lon, lat; coastlines, bathymetry, min_seabed_depth)`** — land test
  *and* a depth test, which is what a workflow that places something in the water actually wants.
- **`extract_marine_cells(ds; coastlines)`** and **`sample_marine_cells(n; bathymetry_ds, ...)`**
  — the mask over a grid, and random points inside the water cells.

The geometry column of a GeoParquet file comes back as well-known binary, so the reader dispatches
on `GeoInterface.trait`, not on `isa GeoInterface.LineString` — that is a trait type, and
testing it as a class silently returns no rings for every file.
"""

# ╔═╡ 7b02a0a8-9b1f-4ef0-9a65-4b4618adf6e8
# ╠═╡
# A closed ring around one area, plus a shallow band along the northern edge that is land in
# its own right, so the ring and the seabed disagree and each can be checked separately.
const RING = (name = "island", code = :island,
              lons = [-64.3, -64.3, -63.8, -63.8, -64.3],
              lats = [43.8, 44.3, 44.3, 43.8, 43.8])

# Elevation, positive up: negative is water, zero is the shoreline.
_elevation = Float64[(lat >= 45.0) ? 5.0 : -80.0 for lon in LON, lat in LAT]
const bathy = make_2d("elevation", _elevation, LON, LAT)

println("point inside  the ring : ", is_point_on_land_geodata(-64.0, 44.0; coastline = [RING]))
println("point outside the ring : ", is_point_on_land_geodata(-62.2, 44.0; coastline = [RING]))
println("no rings at all        : ", is_point_on_land_geodata(-64.0, 44.0),
        "   (unknown, not assumed land)")

masked = extract_marine_cells(bathy; coastlines = [RING])
println("\nmarine cells (ring only)  : ", sum(masked.mask), " of ", length(masked.mask))
masked_bathy = extract_marine_cells(bathy)
println("marine cells (seabed only): ", sum(masked_bathy.mask), " of ", length(masked_bathy.mask),
        "   (the northern land band is excluded by depth, not by the ring)")

pts = sample_marine_cells(6; bathymetry_ds = bathy, coastlines = [RING], jitter_scale = 0.0)
println("sampled water points     : ", length(pts.lons), " -> ", pts.lons)

println("\ndepth floor honoured      : ",
        is_marine_water_geodata(-62.5, 44.0; bathymetry = bathy, min_seabed_depth = 50.0),
        "   (80 m of water clears a 50 m floor)")

# ╔═╡ 8c1b2c09-0a2b-4f31-9b76-5c5729be6d8e
# ╠═╡
md"""
## Grouped storage: the small-file problem, solved once

`GeoStorage` is for a run that writes several artefacts belonging together: a grid, a bathymetry,
and a per-interval field. It exposes groups and variables without the caller building a path
convention:

- **`open_geostorage(path; backend)`** / **`close_geostorage(st)`**
- **`create_storage_group(st, "inputs")`**, **`has_storage_group(st, name)`**
- **`write_storage_variable!(st, path, data; attrs)`**, **`read_storage_variable(st, path)`**
- **`list_storage_groups(st)`**, **`list_storage_variables(st; group)`**

This is what keeps a run's `inputs/` addressable after the fact — the group is a path convention
that a reader can walk, rather than a naming rule each caller has to remember.
"""

# ╔═╡ 9d2c3d1a-1b3c-4a42-9c87-6d683a0f7e9f
# ╠═╡
store = open_geostorage(joinpath(SCRATCH, "run_artefacts.zarr"))
create_storage_group(store, "inputs")
has_storage_group(store, "inputs") || error("group was not created")

write_storage_variable!(store, "inputs/elevation", Float32.(bathy.variables["elevation"].data);
                        attrs = Dict("units" => "m", "convention" => "positive_up"))
write_storage_variable!(store, "inputs/lon", Float64.(LON); attrs = Dict("units" => "degrees_east"))

println("groups     : ", list_storage_groups(store))
println("inputs     : ", list_storage_variables(store; group = "inputs"))
println("read back  : ", size(read_storage_variable(store, "inputs/elevation")))
close_geostorage(store)

# ╔═╡ a1ef4f42-2c5d-4b53-9d98-7e794b1f0a1b
# ╠═╡
md"""
## What this package deliberately does not do

- **No network code.** `GeoData` contains no downloader, no URL, no credentials. Acquisition
  lives in `GeoDataSources`; the split is what lets a consumer depend on the data layer alone.
- **No bundled geography.** Regional ring sets are passed in, not defaulted. An absent file
  means "unknown", which is an answer the land test documents rather than hides.
- **No silent defaults.** An unknown backend, URI, dimension, or missing coordinate is an error.
  Nothing is filled in to make a call succeed.
- **No guessed depth or time.** A dimension either has coordinates or it does not.

### Global Shared Lakehouse Location

`GeoData` provides the canonical location definition across all `Geo*` projects:
- `lakehouse_dir(subpaths...)` resolves paths under `ENV["GEODATA_LAKEHOUSE"]`, `LocalPreferences.toml`, `~/.geodata/storage.toml`, or the workspace fallback.
- `set_lakehouse_dir!(path)` saves the default lakehouse path into `LocalPreferences.toml`.
"""

# ╔═╡ d213dee0-ca2e-11f1-ab07-3dea5077f600
# Cell order:
# ╠═4f2b8a19-0c7d-4e63-9a21-7b0d5e8c1a44
# ╠═5a3c920b-1d8e-4f74-8b32-8c1e6f9d2b55
# ╠═6b4da31c-2e9f-4a05-9c43-9d2f7a0e3c66
# ╠═7c5eb42d-3f0a-4b16-8a54-ae3a8b1f4d77
# ╠═8d6fc53e-401b-4c27-9b65-bef4e5f10a88
# ╠═9e70d64f-512c-4d38-9a76-cfa5a602b9e
# ╠═af81e750-623d-4e49-9a87-d0b6b7130cfa
# ╠═bf92f86a-734e-4f5a-9b98-e1c7c8241d0b
# ╠═c0a39b14-846c-4f6b-9cb9-f2d8d9352e1c
# ╠═d1b4c225-957d-4b6a-9dcb-a13e9f462f0d
# ╠═e2c5d336-068e-4d7a-9aec-b2fa057f3a1b
# ╠═f3d6e447-179f-4b8b-9afd-c3ab168c4b2c
# ╠═04a7f558-2a90-4c9b-8b0f-d4bc279d5c3d
# ╠═159b8604-3ca1-4dac-9d0f-e5cd390ef2d9
# ╠═26bc9715-4b2e-4adb-8e10-f6de4b1fa2b3
# ╠═37cad826-5c3e-4ebb-9f21-070f5c2ab3c4
# ╠═48dbe937-6d4c-4fcc-9a32-181f6d3bc4d5
# ╠═59ecf048-7e5d-4fdd-9b3c-2924088cd4d6
# ╠═6afd1597-8a0e-4fde-9c54-3a35179ce5e7
# ╠═7b02a0a8-9b1f-4ef0-9a65-4b4618adf6e8
# ╠═8c1b2c09-0a2b-4f31-9b76-5c5729be6d8e
# ╠═9d2c3d1a-1b3c-4a42-9c87-6d683a0f7e9f
# ╠═a1ef4f42-2c5d-4b53-9d98-7e794b1f0a1b
# ╠═b1a20f9b-24bc-4d34-9e23-59d08652a9

# ╔═╡ Cell order:
# ╠═4f2b8a19-0c7d-4e63-9a21-7b0d5e8c1a44
# ╠═5a3c920b-1d8e-4f74-8b32-8c1e6f9d2b55
# ╠═6b4da31c-2e9f-4a05-9c43-9d2f7a0e3c66
# ╠═7c5eb42d-3f0a-4b16-8a54-ae3a8b1f4d77
# ╠═8d6fc53e-401b-4c27-9b65-bef4e5f10a88
# ╠═d213b7d0-ca2e-11f1-84d9-e7b99ec2f5a2
# ╠═af81e750-623d-4e49-9a87-d0b6b7130cfa
# ╠═bf92f86a-734e-4f5a-9b98-e1c7c8241d0b
# ╠═c0a39b14-846c-4f6b-9cb9-f2d8d9352e1c
# ╠═d1b4c225-957d-4b6a-9dcb-a13e9f462f0d
# ╠═e2c5d336-068e-4d7a-9aec-b2fa057f3a1b
# ╠═f3d6e447-179f-4b8b-9afd-c3ab168c4b2c
# ╠═04a7f558-2a90-4c9b-8b0f-d4bc279d5c3d
# ╠═159b8604-3ca1-4dac-9d0f-e5cd390ef2d9
# ╠═26bc9715-4b2e-4adb-8e10-f6de4b1fa2b3
# ╠═37cad826-5c3e-4ebb-9f21-070f5c2ab3c4
# ╠═48dbe937-6d4c-4fcc-9a32-181f6d3bc4d5
# ╠═59ecf048-7e5d-4fdd-9b3c-2924088cd4d6
# ╠═6afd1597-8a0e-4fde-9c54-3a35179ce5e7
# ╠═7b02a0a8-9b1f-4ef0-9a65-4b4618adf6e8
# ╠═8c1b2c09-0a2b-4f31-9b76-5c5729be6d8e
# ╠═9d2c3d1a-1b3c-4a42-9c87-6d683a0f7e9f
# ╠═a1ef4f42-2c5d-4b53-9d98-7e794b1f0a1b
# ╠═d213dee0-ca2e-11f1-ab07-3dea5077f600
