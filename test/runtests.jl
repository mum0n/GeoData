using Test
using Dates
using UUIDs
using GeoData
using Zarr, NCDatasets, DataFrames, GeoParquet, Statistics, ArchGDAL, GeoInterface
using GeoData.GeoDataCoordinates: standardize_dimension_name, find_coord_indices, slice_indices, is_regular_grid, grid_spacing
using GeoData.GeoDataCoreTypes: DIM_SPATIAL, DIM_TEMPORAL, DIM_PARAMETRIC, DIM_GENERIC
using GeoData: ZarrBackend, NCDatasetsBackend, YAXArraysBackend, GeoParquetBackend

@testset "GeoData Lakehouse Suite" begin
    # Isolate test files strictly inside scratch directory
    scratch_base = joinpath(@__DIR__, "..", "scratch")
    test_run_dir = joinpath(scratch_base, "test_$(UUIDs.uuid4())")
    mkpath(test_run_dir)

    try
        @testset "Core Types" begin
            dim = Dimension(
                name = :x,
                size = 5,
                coords = collect(1.0:5.0),
                units = "m",
                dim_type = DIM_SPATIAL
            )
            @test dim.name == :x
            @test dim.size == 5
            @test dim.dim_type == DIM_SPATIAL

            crs = CoordinateSystem(crs="EPSG:4326")
            @test !is_cartesian(crs)
            crs_def = CoordinateSystem()
            @test is_cartesian(crs_def)

            data = rand(5, 4)
            dims = [
                Dimension(name=:x, size=5, coords=1.0:5.0),
                Dimension(name=:y, size=4, coords=1.0:4.0)
            ]
            ga = GeoArray(data; dims=dims, crs=crs_def, attrs=Dict{String, Any}("units" => "m"))
            @test size(ga) == (5, 4)
            @test ga.dims[1].name == :x

            ds = GeoDataset(
                variables = Dict("elev" => ga),
                coords = Dict(:x => GeoArray(collect(1.0:5.0); dims=[dims[1]])),
                dims = Dict(:x => dims[1], :y => dims[2])
            )
            @test hasproperty(ds, :elev)
            @test ds.elev === ga
        end

        @testset "Zarr & NetCDF Storage" begin
            lon = collect(-65.0:0.5:-63.0)
            lat = collect(43.0:0.5:45.0)
            depth = collect(0.0:10.0:20.0)
            temp = rand(Float32, length(lon), length(lat), length(depth))

            zarr_path = joinpath(test_run_dir, "grid.zarr")
            g = zgroup(zarr_path; attrs=Dict("title" => "Lakehouse Test"))
            zcreate(Float32, g, "temp", length(lon), length(lat), length(depth); chunks=(length(lon), length(lat), length(depth)))[:] = temp
            zcreate(Float64, g, "lon", length(lon))[:] = lon
            zcreate(Float64, g, "lat", length(lat))[:] = lat
            zcreate(Float64, g, "depth", length(depth))[:] = depth

            ds = geoload(zarr_path)
            @test hasproperty(ds, :temp)
            @test length(ds.dims) == 3
            @test haskey(ds.coords, :lon)

            # Test slice and select
            ds_slice = geoslice(ds, lon=(-64.5, -63.5))
            @test ds_slice.dims[:lon].size <= ds.dims[:lon].size

            ds_sel = geoselect(ds, lon=-64.0, lat=44.0)
            @test ds_sel.dims[:lon].size == 1

            # Point values
            vals = geovalues(ds, ["temp"], lon=-64.0, lat=44.0, depth=10.0)
            @test haskey(vals, "temp")

            # NetCDF export and reload
            nc_path = joinpath(test_run_dir, "grid.nc")
            geosave(nc_path, ds; backend=:ncdatasets)
            ds_nc = geoload(nc_path; backend=:ncdatasets)
            @test hasproperty(ds_nc, :temp)
            @test size(ds_nc.temp.data) == size(ds.temp.data)
        end

        @testset "GeoParquet Backend" begin
            df = DataFrame(lon=[-64.0, -63.5], lat=[43.5, 44.0], val=[10.5, 12.3])
            geoms = [GeoInterface.Point(-64.0, 43.5), GeoInterface.Point(-63.5, 44.0)]
            df[!, :geometry] = geoms
            pq_path = joinpath(test_run_dir, "points.parquet")
            GeoParquet.write(pq_path, df, (:geometry,))

            ds_pq = geoload(pq_path)
            @test hasproperty(ds_pq, :val)
            @test haskey(ds_pq.coords, :lon)
            @test haskey(ds_pq.coords, :lat)
        end

        @testset "Lakehouse Catalog & Medallion Tiering" begin
            cat_file = joinpath(test_run_dir, "catalog.json")
            catalog = GeoDataCatalog(cat_file; auto_save=true)

            # Register Bronze (:raw) entry
            uuid_raw = geopublish!(
                catalog,
                "https://example.com/raw_obs.nc";
                key = :noaa_obs,
                name = "NOAA In-Situ Observations",
                tier = :raw,
                producer = "NOAA",
                variables = ["temp", "salinity"],
                format = :netcdf
            )
            @test haskey(catalog.entries, uuid_raw)

            # Register Silver (:processed) entry
            uuid_proc = geopublish!(
                catalog,
                joinpath(test_run_dir, "grid.zarr");
                key = :temp_harmonized,
                name = "Harmonized Temperature Grid",
                tier = :processed,
                producer = "GeoData.jl",
                derived_from = [:noaa_obs],
                variables = ["temp"],
                format = :zarr
            )

            # Register Gold (:derived) entry
            uuid_derived = geopublish!(
                catalog,
                joinpath(test_run_dir, "kernel.zarr");
                key = :larval_dispersal,
                name = "Larval Dispersal Kernels",
                tier = :derived,
                producer = "ParticleTracking.jl",
                derived_from = [:temp_harmonized],
                variables = ["kernel"],
                format = :zarr
            )

            # Query by tier
            raw_entries = find_datasets(catalog; tier=:raw)
            @test length(raw_entries) == 1
            @test raw_entries[1].key == :noaa_obs

            derived_entries = find_datasets(catalog; tier=:derived)
            @test length(derived_entries) == 1
            @test derived_entries[1].key == :larval_dispersal

            # Lineage tracking
            lineage_proc = find_datasets(catalog; derived_from=:noaa_obs)
            @test length(lineage_proc) == 1
            @test lineage_proc[1].key == :temp_harmonized

            # Fetch dataset
            entry = geofetch(catalog, :larval_dispersal)
            @test entry.tier == :derived
            @test entry.access_count == 1

            # Persistence reload test
            catalog2 = GeoDataCatalog(cat_file)
            @test length(catalog2.entries) == 3
            entry_reloaded = get_dataset_by_key(catalog2, :temp_harmonized)
            @test entry_reloaded.derived_from == [:noaa_obs]
            @test entry_reloaded.tier == :processed
        end

        @testset "Coordinate Helpers & Backend Registry" begin
            @test :zarr in list_backends()
            @test :ncdatasets in list_backends()
            @test :geoparquet in list_backends()

            coords = collect(-65.0:0.5:-63.0)
            idx = find_coord_indices(coords, -64.0; mode=:nearest)
            @test coords[idx] ≈ -64.0
            @test is_regular_grid(coords)
            @test grid_spacing(coords) ≈ 0.5
        end

        include("cube_test.jl")

    finally
        # Target test clean up with retry for Windows file handles
        GC.gc()
        sleep(0.1)
        try
            rm(test_run_dir, recursive=true, force=true)
        catch
            # Fallback if Windows file lock delays deletion
            sleep(0.5)
            try rm(test_run_dir, recursive=true, force=true) catch end
        end
    end
end