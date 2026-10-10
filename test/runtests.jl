"""
Test suite for `GeoData.jl`.

Every operation test asserts the new contract (see `PLAN.md`, Phase 1), and every
dataset produced is checked with `assert_dataset_invariants`. Network access is
never used: sources tests live in `GeoDataSources` against fixtures.
"""

using Test
using Dates
using UUIDs
using Statistics
using GeoData
using Zarr, NCDatasets, DataFrames, GeoParquet, ArchGDAL, GeoInterface
using GeoData.GeoDataCoordinates: standardize_dimension_name, find_coord_indices,
    slice_indices, is_regular_grid, grid_spacing
using GeoData.GeoDataTypes: DIM_SPATIAL, DIM_TEMPORAL, DIM_PARAMETRIC, DIM_GENERIC,
    Dimension, GeoArray, GeoDataset
using GeoData: ZarrBackend, NCDatasetsBackend, GeoParquetBackend
using GeoData: GeoDataOperations

include("helpers.jl")

function make_dataset(; descending_lat = false)
    lon = collect(-65.0:0.5:-63.0)
    lat = descending_lat ? collect(45.0:-0.5:43.0) : collect(43.0:0.5:45.0)
    depth = collect(0.0:10.0:20.0)
    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    dd = Dimension(name = :depth, size = length(depth), coords = depth, units = "m",
                   dim_type = DIM_PARAMETRIC)
    data = rand(Float32, length(lon), length(lat), length(depth))
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict("temp" => GeoArray(data, (dl, dla, dd), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}()),
              :depth => GeoArray(depth, (dd,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla, :depth => dd), crs, Dict{String, Any}(), nothing,
        "synthetic")
end

function make_2d_dataset(data::AbstractMatrix, lon, lat, varname)
    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict(varname => GeoArray(data, (dl, dla,), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla), crs, Dict{String, Any}(), nothing, "synthetic")
end

# Like `make_2d_dataset` but with the dimension order stated explicitly, so a test can
# build the same logical grid in either storage orientation and check that the reader
# resolves orientation rather than inferring it from `size`.
function make_2d_dataset_ordered(data::AbstractMatrix, lon, lat, varname, order)
    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    dims = order == :lon_lat ? (dl, dla) : (dla, dl)
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict(varname => GeoArray(data, dims, crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla), crs, Dict{String, Any}(), nothing, "synthetic")
end

function make_3d_dataset(data::AbstractArray{<:Any, 3}, lon, lat, dep, varname)    dl = Dimension(name = :lon, size = length(lon), coords = lon, units = "degrees_east",
                   dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = length(lat), coords = lat, units = "degrees_north",
                    dim_type = DIM_SPATIAL)
    dd = Dimension(name = :depth, size = length(dep), coords = dep, units = "m",
                   dim_type = DIM_PARAMETRIC)
    crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
    return GeoDataset(
        Dict(varname => GeoArray(data, (dl, dla, dd), crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(lon, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat, (dla,), crs, Dict{String, Any}()),
             :depth => GeoArray(dep, (dd,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla, :depth => dd), crs, Dict{String, Any}(), nothing,
        "synthetic")
end

@testset "GeoData Lakehouse Suite" begin
    test_run_dir = test_scratch()
    try
        @testset "Core Types" begin
            dim = Dimension(name = :x, size = 5, coords = collect(1.0:5.0), units = "m",
                            dim_type = DIM_SPATIAL)
            @test dim.name == :x && dim.size == 5 && dim.dim_type == DIM_SPATIAL
            @test dim.units == "m"

            # Dimension(existing) replaced fields: every other field survives.
            narrowed = Dimension(dim; size = 3)
            @test narrowed.units == dim.units
            @test narrowed.dim_type == dim.dim_type
            @test narrowed.size == 3

            crs = GeoData.CoordinateSystem(crs = "EPSG:4326")
            @test !GeoData.GeoDataTypes.is_cartesian(crs)
            @test GeoData.GeoDataTypes.is_cartesian(GeoData.CoordinateSystem())

            ga = GeoArray(rand(5, 4); dims = [Dimension(name = :x, size = 5, coords = 1.0:5.0),
                                              Dimension(name = :y, size = 4, coords = 1.0:4.0)])
            @test size(ga) == (5, 4) && ndims(ga) == 2
            # IndexStyle is looked up by array utilities, not only by indexing.
            @test Base.IndexStyle(typeof(ga)) == Base.IndexStyle(Array)
        end

        @testset "Dataset invariants" begin
            ds = make_dataset()
            assert_dataset_invariants(ds)

            dl = ds.temp.dims[1]
            # a variable spanning a dimension the dataset does not declare
            undeclared = GeoDataset(
                Dict("temp" => GeoArray(rand(5, 3), (dl, Dimension(:time, 3)), ds.crs,
                                        Dict{String, Any}())),
                Dict{Symbol, GeoArray}(), Dict(:lon => dl), ds.crs, Dict{String, Any}(),
                nothing, "undeclared")
            @test_throws Exception assert_dataset_invariants(undeclared)

            # a coordinate whose length does not match its dimension
            mismatched = GeoDataset(
                Dict("temp" => GeoArray(rand(5, 5, 3), ds.temp.dims, ds.crs, Dict{String, Any}())),
                Dict(:lon => GeoArray(collect(1.0:4.0), (dl,), ds.crs, Dict{String, Any}())),
                ds.dims, ds.crs, Dict{String, Any}(), nothing, "mismatched")
            @test_throws Exception assert_dataset_invariants(mismatched)
        end

        @testset "Slicing and selection contract" begin
            ds = make_dataset()

            # A range keeps the dimension.
            s = GeoData.geoslice(ds, lon = (-64.5, -63.5))
            assert_dataset_invariants(s)
            @test ndims(s.temp) == 3
            @test s.dims[:lon].size == 3

            # A single value drops the dimension everywhere, consistently.
            sel = GeoData.geoselect(ds, lon = -64.0, lat = 44.0)
            assert_dataset_invariants(sel)
            @test ndims(sel.temp) == 1
            @test !haskey(sel.dims, :lon) && !haskey(sel.coords, :lon)
            @test sel.temp.dims[1].name == :depth
            @test sel.attrs["selections"][:lon] == -64.0

            # A vector selects those points, in order, with duplicates kept.
            vec = GeoData.geoslice(ds, lon = [-64.0, -63.5, -64.0])
            assert_dataset_invariants(vec)
            @test vec.dims[:lon].size == 3
            @test vec.coords[:lon].data == [-64.0, -63.5, -64.0]

            # An unknown dimension is an error, not a no-op.
            @test_throws Exception GeoData.geoslice(ds, chlor = 1.0)

            # An empty range snaps to nothing rather than returning the whole axis.
            @test_throws Exception GeoData.geoslice(ds, lon = (0.0, 10.0))

            # Bounding-box helper ignores dimensions the dataset does not have.
            b = GeoData.geoslice_bbox(ds, lon_range = (-64.5, -63.5), lat_range = (43.0, 44.0),
                                      depth_range = (0.0, Inf), time_range = (0.0, 1.0))
            assert_dataset_invariants(b)
            @test b.dims[:lon].size == 3
        end

        @testset "Point value queries" begin
            ds = make_dataset()

            v = GeoData.geovalues(ds, ["temp"], lon = -64.0, lat = 44.0, depth = 10.0)
            @test v["temp"] isa Real

            v2 = GeoData.geovalues(ds, ["temp"], lon = [-64.0, -63.5], lat = [44.0, 44.5])
            @test v2["temp"] isa AbstractArray
            # A dimension that is not selected is kept whole: four points, all depths.
            @test size(v2["temp"]) == (2, 2, 3)

            v3 = GeoData.geovalues(ds, ["temp"], lon = -64.0, lat = 44.0)
            @test v3["temp"] isa AbstractVector
            @test length(v3["temp"]) == 3

            @test_throws Exception GeoData.geovalues(ds, ["temp"])
            @test_throws Exception GeoData.geovalues(ds, ["temp"], lon = [-64.0, -63.5, -64.0],
                                                     lat = [44.0, 44.5])
            @test_throws Exception GeoData.geovalues(ds, ["nope"], lon = -64.0)
        end

        @testset "Joining" begin
            ds1 = make_dataset()
            ds2 = GeoDataset(
                Dict("sal" => GeoArray(rand(Float32, 5, 5, 3), ds1.temp.dims,
                                       ds1.temp.crs, Dict{String, Any}())),
                ds1.coords, ds1.dims, ds1.crs, Dict{String, Any}(), nothing, "synthetic2")

            inner = GeoData.geojoin([ds1, ds2], on = [:lon, :lat, :depth], how = :inner)
            assert_dataset_invariants(inner)
            @test haskey(inner.variables, "temp") && haskey(inner.variables, "sal")
            @test inner.temp.dims[1].size == 5

            # Left join over a partially overlapping dataset fills NaN, never nearest.
            other_lon = collect(-64.0:0.5:-63.0)
            lon_dim = Dimension(name = :lon, size = 3, coords = other_lon)
            ds3 = GeoDataset(
                Dict("sal" => GeoArray(rand(Float32, 3, 5, 3),
                                       (lon_dim, ds1.temp.dims[2], ds1.temp.dims[3]),
                                       ds1.crs, Dict{String, Any}())),
                merge(ds1.coords, Dict(:lon => GeoArray(other_lon, (lon_dim,), ds1.crs,
                                                         Dict{String, Any}()))),
                merge(ds1.dims, Dict(:lon => lon_dim)), ds1.crs, Dict{String, Any}(), nothing,
                "synthetic3")
            # Left join over a partially overlapping dataset fills NaN, never nearest:
            # `sal` exists only at the three longitudes ds3 carries.
            left = GeoData.geojoin([ds1, ds3], on = [:lon, :lat, :depth], how = :left)
            assert_dataset_invariants(left)
            @test left.dims[:lon].size == 5
            @test any(isnan, left.sal.data)
            @test !any(isnan, left.temp.data)

            @test_throws Exception GeoData.geojoin([ds1, ds2], on = [])
            @test_throws Exception GeoData.geojoin([ds1, ds2], on = [:lon], how = :cross)
            # Name collisions must not be silently suffixed.
            ds_dup = GeoDataset(Dict("temp" => ds1.temp), ds1.coords, ds1.dims, ds1.crs,
                                Dict{String, Any}(), nothing, "dup")
            @test_throws Exception GeoData.geojoin([ds1, ds_dup], on = [:lon, :lat, :depth])
        end

        @testset "Aggregation" begin
            ds = make_dataset()
            agg = GeoData.geoaggregate(ds, dim = :lat, func = mean)
            assert_dataset_invariants(agg)
            @test ndims(agg.temp) == 2
            @test agg.temp.dims[1].name == :lon
            @test !haskey(agg.dims, :lat)
            @test agg.dims[:lon].size == 5
            @test size(agg.temp.data) == (5, 3)
            @test haskey(agg.attrs, "reduced")

            # A variable that does not span the dimension is kept, not dropped.
            surf = GeoDataset(Dict("temp" => ds.temp, "mask" => GeoArray(trues(5), (ds.temp.dims[1],),
                                                                        ds.crs, Dict{String, Any}())),
                              ds.coords, ds.dims, ds.crs, Dict{String, Any}(), nothing, "s")
            agg2 = GeoData.geoaggregate(surf, dim = :lat, func = sum)
            @test haskey(agg2.variables, "mask")
            @test_throws Exception GeoData.geoaggregate(ds, dim = :time, func = mean)
        end

        @testset "Backend capabilities" begin
            # Every backend must report capabilities that exist as fields on the struct.
            # `create` was passed by all three and silently broke the accessors.
            for name in GeoData.list_backends()
                be = GeoData.get_backend(name)
                caps = GeoData.backend_capabilities(be)
                @test caps isa GeoData.BackendCapabilities
                @test caps.read
                @test :lazy in propertynames(caps)
                @test :chunked in propertynames(caps)
            end
            # Zarr is the lazy backend; NetCDF and GeoParquet materialise.
            @test GeoData.backend_capabilities(GeoData.get_backend(:zarr)).lazy
            @test !GeoData.backend_capabilities(GeoData.get_backend(:ncdatasets)).lazy
            # The default advertises read and write only. Tested through a minimal
            # backend because `GeoBackend` is abstract and cannot be instantiated.
            struct ProbeBackend <: GeoData.GeoBackend end
            @test GeoData.backend_capabilities(ProbeBackend()) ==
                  GeoData.BackendCapabilities(read = true, write = true)
        end

        @testset "Zarr & lazy reads" begin
            ds = make_dataset()
            zarr_path = joinpath(test_run_dir, "grid.zarr")
            GeoData.geosave(zarr_path, ds; backend = :zarr)
            @test isdir(zarr_path)

            loaded = GeoData.geoload(zarr_path)
            assert_dataset_invariants(loaded)
            @test hasproperty(loaded, :temp)
            @test size(loaded.temp.data) == size(ds.temp.data)
            @test haskey(loaded.coords, :lon)
            @test loaded.temp.data isa Zarr.ZArray     # nothing was read at open time
            # The Zarr round trip preserves values and keeps the store lazy.
            @test GeoData.geoslice(loaded, lon = (-64.5, -63.5)).temp.data[1, 1, 1] ==
                  ds.temp.data[2, 1, 1]

            nc_path = joinpath(test_run_dir, "grid.nc")
            GeoData.geosave(nc_path, ds; backend = :ncdatasets)
            ds_nc = GeoData.geoload(nc_path; backend = :ncdatasets)
            assert_dataset_invariants(ds_nc)
            @test hasproperty(ds_nc, :temp)
            @test size(ds_nc.temp.data) == size(ds.temp.data)

            # A store that already exists cannot be updated in place; replacing is explicit.
            @test_throws Exception GeoData.geosave(zarr_path, ds; backend = :zarr)
            GeoData.geosave(zarr_path, ds; backend = :zarr, overwrite = true)
            @test GeoData.geoload(zarr_path).temp.data[2, 1, 1] == ds.temp.data[2, 1, 1]
        end

        @testset "GeoParquet Backend" begin
            df = DataFrame(lon = [-64.0, -63.5], lat = [43.5, 44.0], val = [10.5, 12.3])
            df[!, :geometry] = [GeoInterface.Point(-64.0, 43.5), GeoInterface.Point(-63.5, 44.0)]
            pq_path = joinpath(test_run_dir, "points.parquet")
            GeoParquet.write(pq_path, df, (:geometry,))

            ds_pq = GeoData.geoload(pq_path)
            assert_dataset_invariants(ds_pq)
            @test haskey(ds_pq.variables, "val")
            @test haskey(ds_pq.variables, "geometry")   # geometry survives the round trip
            @test haskey(ds_pq.coords, :lon)
            @test ds_pq.coords[:lon].data == [-64.0, -63.5]
            @test ds_pq.dims[:points].size == 2
            @test haskey(ds_pq.coords, :lat)

            out = joinpath(test_run_dir, "points_out.parquet")
            GeoData.geosave(out, ds_pq)
            rt = GeoData.geoload(out)
            @test rt.coords[:lon].data == ds_pq.coords[:lon].data
        end

        @testset "Registry, schemes, infer_backend" begin
            @test sort(GeoData.list_backends()) == [:geoparquet, :ncdatasets, :nczarr, :zarr]
            @test GeoData.infer_backend("a.zarr") == :zarr
            @test GeoData.infer_backend("zarr://a.zarr") == :zarr
            @test GeoData.infer_backend("a.nc") == :ncdatasets
            @test GeoData.infer_backend("netcdf://a.nc") == :ncdatasets
            @test GeoData.infer_backend("nczarr://store") == :nczarr
            @test GeoData.infer_backend("a.parquet") == :geoparquet
            @test GeoData.infer_backend("geoparquet://a.parquet") == :geoparquet
            @test_throws Exception GeoData.infer_backend("a.tif")

            @test GeoData.strip_scheme("zarr://data/grid.zarr") == "data/grid.zarr"
            @test GeoData.strip_scheme("geoparquet://data/p.parquet") == "data/p.parquet"
            @test GeoData.strip_scheme("nczarr://store/group") == "store/group"
            @test GeoData.strip_scheme("file://data/x.nc") == "data/x.nc"
            @test GeoData.strip_scheme("/plain/path.zarr") == "/plain/path.zarr"
        end

        @testset "Coordinate Helpers" begin
            coords = collect(-65.0:0.5:-63.0)
            idx = find_coord_indices(coords, -64.0; mode = :nearest)
            @test coords[idx] ≈ -64.0
            @test is_regular_grid(coords)
            @test grid_spacing(coords) ≈ 0.5
            @test standardize_dimension_name(:x) == :lon
            @test GeoData.normalize_depth([-10.0, -20.0])[1] == [10.0, 20.0]

            # Descending axes: the mirroring bug is a real one, so check both orders.
            for desc in (false, true)
                ds = make_dataset(; descending_lat = desc)
                s = GeoData.geoslice(ds, lat = (44.0, 44.5))
                assert_dataset_invariants(s)
                @test s.coords[:lat].data == (desc ? [44.5, 44.0] : [44.0, 44.5])
            end
        end

        @testset "Bathymetry analysis and land exclusion" begin
            lons = collect(1.0:10.0)
            lats = collect(1.0:10.0)
            elev = fill(-50.0, 10, 10)          # uniformly deep
            ds_b = make_2d_dataset(elev, lons, lats, "elevation")

            water_only = extract_marine_cells(ds_b)
            @test all(water_only.mask)            # depth alone says water everywhere

            # A closed square ring over lon/lat [3.5, 6.5]
            ring = [(lon = 3.5, lat = 3.5), (lon = 6.5, lat = 3.5),
                    (lon = 6.5, lat = 6.5), (lon = 3.5, lat = 6.5),
                    (lon = 3.5, lat = 3.5)]
            rings = [(name = "island", code = :island, lons = [p.lon for p in ring],
                      lats = [p.lat for p in ring])]

            with_land = extract_marine_cells(ds_b; coastlines = rings)
            @test with_land.mask[4, 4] == false      # inside the ring
            @test with_land.mask[5, 5] == false
            @test with_land.mask[4, 6] == false
            @test with_land.mask[3, 3] == true      # outside
            @test count(with_land.mask) == count(water_only.mask) - 9
            @test length(with_land.indices) == 100 - 9

            # An *open* ring is skipped, not ray-cast: the mask is unchanged rather than
            # land being reported over half the plane.
            open_ring = [(lons = [3.5, 6.5, 6.5, 3.5], lats = [3.5, 3.5, 6.5, 6.5],
                          name = "open", code = Symbol("open"))]
            open_test = extract_marine_cells(ds_b; coastlines = open_ring)
            @test open_test.mask == water_only.mask

            # The depth threshold applies too.
            @test isempty(extract_marine_cells(ds_b; min_seabed_depth = 60.0).indices)

            # Sampling returns exactly n points, inside the water cells. Jitter is off, so
            # every point lands on a cell centre, and none lands inside the land ring.
            pts = sample_marine_cells(25; bathymetry_ds = ds_b, coastlines = rings,
                                       jitter_scale = 0.0)
            @test length(pts.lons) == 25
            @test all(1.0 .<= pts.lons .<= 10.0) && all(1.0 .<= pts.lats .<= 10.0)
            in_land = [(lon in (4.0, 5.0, 6.0)) && (lat in (4.0, 5.0, 6.0))
                       for (lon, lat) in zip(pts.lons, pts.lats)]
            @test !any(in_land)
            @test all(pts.depths .<= -50.0)
        end

        @testset "extract_marine_cells orientation and depth pairing" begin
            # A square, uniform grid hides two orientation errors: a shape check against
            # the wrong axis order, and a transposed `elev[i, j]` read. This testset uses a
            # non-square grid whose elevation identifies the cell it belongs to, so both
            # fail loudly. Orientation is taken from the variable's own dimension records.
            lon2 = collect(-65.0:0.25:-61.0)          # 17 lon
            lat2 = collect(43.0:0.25:45.0)            # 9 lat
            n_lon, n_lat = length(lon2), length(lat2)

            # Cell-identifying depth: -(100 + lat + lon/100)
            function answer(lon, lat)
                -(100.0 + lat + lon / 100)
            end
            for order in (:lon_lat, :lat_lon), desc in (false, true)
                lo = desc ? reverse(lon2) : lon2
                la = desc ? reverse(lat2) : lat2
                raw = [answer(lo[i], la[j]) for i in 1:n_lon, j in 1:n_lat]
                data = order == :lon_lat ? raw : permutedims(raw, (2, 1))
                ds = make_2d_dataset_ordered(data, lo, la, "elevation", order)

                mc = extract_marine_cells(ds; min_seabed_depth = 10.0)
                tag = "$(order), $(desc ? "descending" : "ascending")"
                @test length(mc.indices) == n_lon * n_lat
                @test size(mc.mask) == (n_lat, n_lon)
                # Every returned cell's depth must be the one at its own (lon, lat).
                @test (all(answer.(mc.lons, mc.lats) .≈ mc.elevation))
                # Both axes come back ascending, whatever order the source stored them in.
                @test all(-65.0 .<= mc.lons .<= -61.0) && all(43.0 .<= mc.lats .<= 45.0)
                @test length(unique(mc.lons)) == n_lon && length(unique(mc.lats)) == n_lat

                # Sampling inherits the same pairing, and its depth weights come from it.
                s = sample_marine_cells(40; bathymetry_ds = ds, jitter_scale = 0.0)
                @test (all(answer.(s.lons, s.lats) .≈ s.depths))
            end

            # A variable that does not span both axes is refused, not flattened.
            bad = make_3d_dataset(rand(n_lon, n_lat, 3), lon2, lat2,
                                  collect(0.0:1.0:2.0), "elevation")
            @test_throws ErrorException extract_marine_cells(bad)

            # Data whose extent disagrees with the declared axes is named, not guessed.
            dl = Dimension(name = :lon, size = n_lon, coords = lon2,
                           dim_type = DIM_SPATIAL)
            dla = Dimension(name = :lat, size = n_lat, coords = lat2,
                            dim_type = DIM_SPATIAL)
            crs = CoordinateSystem(crs = "EPSG:4326")
            mismatched = GeoDataset(
                Dict("elevation" => GeoArray(rand(n_lon + 2, n_lat), (dl, dla), crs,
                                             Dict{String, Any}())),
                Dict(:lon => GeoArray(lon2, (dl,), crs, Dict{String, Any}()),
                     :lat => GeoArray(lat2, (dla,), crs, Dict{String, Any}())),
                Dict(:lon => dl, :lat => dla), crs, Dict{String, Any}(), nothing, "bad")
            @test_throws ErrorException extract_marine_cells(mismatched)
        end

        @testset "Migrated pure-computation functions" begin
            # Every function moved out of GeoDataSources into core/. They were unexercised
            # after the split and several of them were broken when first run here, so the
            # split is now verified rather than assumed.
            rings = [(name = "island", code = :island,
                       lons = [0.0, 1.0, 1.0, 0.0, 0.0], lats = [0.0, 0.0, 1.0, 1.0, 0.0])]

            @test buffer_distance_to_degrees(111_000.0, 60.0)[1] >
                  buffer_distance_to_degrees(111_000.0, 0.0)[1]

            (blon, blat) = expand_domain_with_buffer((-65.0, -63.0), (43.0, 45.0);
                                                     buffer_km = 20.0)
            @test blon[1] < -65.0 && blat[2] > 45.0

            sq_lons = [0.0, 1.0, 1.0, 0.0]; sq_lats = [0.0, 0.0, 1.0, 1.0]
            @test point_in_polygon(0.5, 0.5, sq_lons, sq_lats)
            @test !point_in_polygon(2.0, 0.5, sq_lons, sq_lats)
            # A vertex is a boundary case: even-odd ray casting says false, and that is the
            # contract, not a gap to close with an epsilon.
            @test !point_in_polygon(1.0, 1.0, sq_lons, sq_lats)

            lon = collect(-65.0:0.25:-63.0); lat = collect(43.0:0.25:45.0)
            dep = collect(0.0:10.0:50.0)
            tgt_lon = collect(-65.0:0.05:-63.8)
            tgt_lat = collect(43.0:0.05:43.8)

            src2 = make_2d_dataset(rand(9, 9), lon, lat, "ssh")
            @test size(regrid_2d_field(lon, lat, src2.variables["ssh"].data, tgt_lon, tgt_lat)) ==
                  (25, 17)
            @test size(regrid_2d_field(src2, "ssh", tgt_lon, tgt_lat)) == (25, 17)

            src3 = make_3d_dataset(rand(9, 9, 6), lon, lat, dep, "t")
            @test size(regrid_3d_field(src3, "t", tgt_lon, tgt_lat, dep)) == (25, 17, 6)
            # Descending source latitude must be sorted, not silently mirrored.
            src3r = make_3d_dataset(rand(9, 9, 6), lon, reverse(lat), dep, "t")
            @test size(regrid_3d_field(src3r, "t", tgt_lon, tgt_lat, dep)) == (25, 17, 6)

            @test size(smooth_bathymetry(sign.(rand(12, 12) .- 0.5) .* 2.0; passes = 2)) ==
                  (12, 12)

            it = get_bathymetry_interpolator(make_2d_dataset(fill(-100.0, 4, 4),
                                                             [1.0, 2.0, 3.0, 4.0],
                                                             [1.0, 2.0, 3.0, 4.0], "elevation"))
            @test it(2.0, 2.0) == -100.0

            # An absent file yields no rings ("unknown"), documented, not a crash.
            @test isempty(load_coastline_polygons("C:/nope/nowhere.parquet"))

            @test is_point_on_land_geodata(0.5, 0.5; coastline = rings)
            @test !is_point_on_land_geodata(9.9, 9.9; coastline = rings)
            @test !is_point_on_land_geodata(9.9, 9.9)      # no input: unknown, not land
            open = [(name = "o", code = :o, lons = [0.0, 1.0, 1.0, 0.0],
                     lats = [0.0, 0.0, 1.0, 1.0])]
            @test !is_point_on_land_geodata(0.5, 0.5; coastline = open)   # open line skipped

            bds = make_2d_dataset(fill(-100.0, 6, 6), collect(0.0:1.0:5.0),
                                  collect(0.0:1.0:5.0), "elevation")
            @test is_marine_water_geodata(9.9, 9.9; coastlines = rings)
            @test !is_marine_water_geodata(0.5, 0.5; coastlines = rings)
            @test is_marine_water_geodata(2.5, 2.5; bathymetry = bds, min_seabed_depth = 50.0)
            @test !is_marine_water_geodata(2.5, 2.5; bathymetry = bds, min_seabed_depth = 200.0)

            nm, ga = variables_like(make_2d_dataset(rand(3, 3), [1.0, 2.0, 3.0],
                                                    [1.0, 2.0, 3.0], "temp_C"), ["temp_C"])
            @test nm == "temp_C" && ga isa GeoArray

            v, flipped = normalize_depth([-5.0, -10.0])
            @test v == [5.0, 10.0] && flipped
            v2, flipped2 = normalize_depth([5.0, 10.0])
            @test v2 == [5.0, 10.0] && !flipped2

            # axis resolves aliases, and errors on an absent axis rather than inventing one
            ali = make_2d_dataset(rand(3, 3), [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], "t")
            vl, dl = axis(ali, :longitude)
            @test vl == [1.0, 2.0, 3.0] && dl.name == :lon
            @test axis(ali, :x)[1] == [1.0, 2.0, 3.0]
            @test_throws ErrorException axis(ali, :depth)

            bb = bounding_box(ali)
            @test bb.lon == (1.0, 3.0) && bb.lat == (1.0, 3.0)
            @test is_regular_grid([1.0, 2.0, 3.0]) && !is_regular_grid([1.0, 2.0, 4.0])
            @test grid_spacing([1.0, 2.0, 3.0]) == 1.0
            @test find_coord_indices(collect(0.0:1.0:10.0), 4.4; mode = :nearest) == 5
            @test slice_indices(collect(0.0:1.0:10.0), (3.0, 6.0)) == 4:7
        end
    finally
        GC.gc()
        sleep(0.1)
        try
            rm(test_run_dir; recursive = true, force = true)
        catch
            sleep(0.5)
            try rm(test_run_dir; recursive = true, force = true) catch end
        end
    end
end
