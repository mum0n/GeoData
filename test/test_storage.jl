include("test_fixtures.jl")

@testset "Storage, Backends & Formats" begin
    test_run_dir = test_scratch()
    try
        # Backend capabilities
        for name in GeoData.list_backends()
            be = GeoData.get_backend(name)
            caps = GeoData.backend_capabilities(be)
            @test caps isa GeoData.BackendCapabilities
            @test caps.read
            @test :lazy in propertynames(caps)
            @test :chunked in propertynames(caps)
        end
        @test GeoData.backend_capabilities(GeoData.get_backend(:zarr)).lazy
        @test !GeoData.backend_capabilities(GeoData.get_backend(:ncdatasets)).lazy
        struct ProbeBackend <: GeoData.GeoBackend end
        @test GeoData.backend_capabilities(ProbeBackend()) ==
              GeoData.BackendCapabilities(read = true, write = true)

        # Zarr & lazy reads
        ds = make_dataset()
        zarr_path = joinpath(test_run_dir, "grid.zarr")
        GeoData.geosave(zarr_path, ds; backend = :zarr)
        @test isdir(zarr_path)

        loaded = GeoData.geoload(zarr_path)
        assert_dataset_invariants(loaded)
        @test hasproperty(loaded, :temp)
        @test size(loaded.temp.data) == size(ds.temp.data)
        @test haskey(loaded.coords, :lon)
        @test loaded.temp.data isa Zarr.ZArray
        @test GeoData.geoslice(loaded, lon = (-64.5, -63.5)).temp.data[1, 1, 1] ==
              ds.temp.data[2, 1, 1]

        nc_path = joinpath(test_run_dir, "grid.nc")
        GeoData.geosave(nc_path, ds; backend = :ncdatasets)
        ds_nc = GeoData.geoload(nc_path; backend = :ncdatasets)
        assert_dataset_invariants(ds_nc)
        @test hasproperty(ds_nc, :temp)
        @test size(ds_nc.temp.data) == size(ds.temp.data)

        @test_throws Exception GeoData.geosave(zarr_path, ds; backend = :zarr)
        GeoData.geosave(zarr_path, ds; backend = :zarr, overwrite = true)
        @test GeoData.geoload(zarr_path).temp.data[2, 1, 1] == ds.temp.data[2, 1, 1]

        # GeoParquet Backend
        df = DataFrame(lon = [-64.0, -63.5], lat = [43.5, 44.0], val = [10.5, 12.3])
        df[!, :geometry] = [GeoInterface.Point(-64.0, 43.5), GeoInterface.Point(-63.5, 44.0)]
        pq_path = joinpath(test_run_dir, "points.parquet")
        GeoParquet.write(pq_path, df, (:geometry,))

        ds_pq = GeoData.geoload(pq_path)
        assert_dataset_invariants(ds_pq)
        @test haskey(ds_pq.variables, "val")
        @test haskey(ds_pq.variables, "geometry")
        @test haskey(ds_pq.coords, :lon)
        @test ds_pq.coords[:lon].data == [-64.0, -63.5]
        @test ds_pq.dims[:points].size == 2
        @test haskey(ds_pq.coords, :lat)

        out = joinpath(test_run_dir, "points_out.parquet")
        GeoData.geosave(out, ds_pq)
        rt = GeoData.geoload(out)
        @test rt.coords[:lon].data == ds_pq.coords[:lon].data

        # Registry, schemes, infer_backend
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

        # Appending NetCDF
        scratch_nc = test_scratch()
        ds_nc3 = make_3d_dataset([Float32(z + d) for x in 1:4, y in 1:3, (d, z) in enumerate(-2.0:-2.0:-6.0)],
                             collect(1.0:4.0), collect(1.0:3.0), collect(-2.0:-2.0:-6.0), "temp")
        app_path = joinpath(scratch_nc, "append.nc")
        geosave(app_path, ds_nc3; backend = :ncdatasets)
        geosave(app_path, ds_nc3; backend = :ncdatasets)
        @test size(geoload(app_path; backend = :ncdatasets).variables["temp"].data) == size(ds_nc3.variables["temp"].data)

        ds2_nc = make_2d_dataset([Float32(1.0 * (x + y)) for x in 1:4, y in 1:3],
                              collect(1.0:4.0), collect(1.0:3.0), "salt")
        geosave(app_path, ds2_nc; backend = :ncdatasets)
        both = geoload(app_path; backend = :ncdatasets)
        @test sort(collect(keys(both.variables))) == ["salt", "temp"]

        # Profile and timeseries
        ds_prof = make_3d_dataset([Float32(x + y + d) for x in 1:5, y in 1:4, d in 1:3],
                             collect(1.0:5.0), collect(1.0:4.0), collect(-5.0:-5.0:-15.0), "temp")
        prof = geoprofile(ds_prof, ["temp"]; lon = 3.0, lat = 2.0)
        @test length(prof[:depth]) == 3
        @test length(prof["temp"]) == 3
        @test_throws ErrorException geoprofile(ds_prof, ["temp"]; lon = 3.0, lat = 2.0, time = 0.0)

        crs_time = GeoData.CoordinateSystem(crs = "EPSG:4326")
        dl_t = Dimension(name = :lon, size = 5, coords = collect(1.0:5.0), units = "degrees_east",
                       dim_type = DIM_SPATIAL)
        dla_t = Dimension(name = :lat, size = 4, coords = collect(1.0:4.0), units = "degrees_north",
                        dim_type = DIM_SPATIAL)
        dd_t = Dimension(name = :depth, size = 3, coords = collect(-5.0:-5.0:-15.0), units = "m",
                       dim_type = DIM_PARAMETRIC)
        dt_t = Dimension(name = :time, size = 2, coords = [0.0, 30.0], units = "days",
                       dim_type = DIM_TEMPORAL)
        tds = GeoDataset(
            Dict("temp" => GeoArray(
                Float32[x + y + d + t for x in 1:5, y in 1:4, d in 1:3, t in 1:2],
                (dl_t, dla_t, dd_t, dt_t), crs_time, Dict{String, Any}())),
            Dict(:lon => GeoArray(collect(1.0:5.0), (dl_t,), crs_time, Dict{String, Any}()),
                 :lat => GeoArray(collect(1.0:4.0), (dla_t,), crs_time, Dict{String, Any}()),
                 :depth => GeoArray(collect(-5.0:-5.0:-15.0), (dd_t,), crs_time, Dict{String, Any}()),
                 :time => GeoArray([0.0, 30.0], (dt_t,), crs_time, Dict{String, Any}())),
            Dict(:lon => dl_t, :lat => dla_t, :depth => dd_t, :time => dt_t), crs_time, Dict{String, Any}(),
            nothing, "synthetic")

        ts = geotimeseries(tds, ["temp"]; lon = 3.0, lat = 2.0)
        @test ts[:time] == [0.0, 30.0]
        @test length(ts["temp"]) == 2
        ts_d = geotimeseries(tds, ["temp"]; lon = 3.0, lat = 2.0, depth = -10.0)
        @test ts_d["temp"] ≠ ts["temp"]

        # GeoStorage listing
        scratch_store = test_scratch()
        store = open_geostorage(joinpath(scratch_store, "store.zarr"))
        create_storage_group(store, "inputs")
        create_storage_group(store, "inputs/nested")
        m_mat = [Float32(1.0 * (x + y)) for x in 1:7, y in 1:5]
        write_storage_variable!(store, "inputs/elevation", m_mat)
        write_storage_variable!(store, "inputs/lon", Float64.(1.0:7.0))

        @test list_storage_groups(store) == ["inputs"]
        @test list_storage_groups(store; group_path = "inputs") == ["nested"]
        @test list_storage_groups(store; group_path = "absent") == String[]
        @test list_storage_variables(store) == String[]
        @test list_storage_variables(store; group = "inputs") == ["elevation", "lon"]

        back = read_storage_variable(store, "inputs/elevation")
        @test size(back) == size(m_mat)
        @test vec(back) == vec(m_mat)
        close_geostorage(store)
    finally
        GC.gc()
        sleep(0.1)
        try rm(test_run_dir; recursive = true, force = true) catch end
    end
end
