using Test
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using GeoData
using Zarr, NCDatasets, DataFrames, GeoParquet, Statistics, ArchGDAL, GeoInterface
using GeoData.GeoDataCoordinates: standardize_dimension_name, find_coord_indices, slice_indices, is_regular_grid, grid_spacing
using GeoData.GeoDataCoreTypes: DIM_SPATIAL, DIM_TEMPORAL, DIM_PARAMETRIC, DIM_GENERIC
using GeoData: ZarrBackend, NCDatasetsBackend, YAXArraysBackend, GeoParquetBackend

@testset "GeoData.jl" begin
    # Clean up test directory
    test_dir = joinpath(@__DIR__, "..", "test_data")
    isdir(test_dir) && rm(test_dir, recursive=true, force=true)
    mkpath(test_dir)

    @testset "Core Types" begin
        # Test Dimension
        dim = Dimension(name=:x, size=10, coords=collect(1.0:10.0), units="m", dim_type=DIM_SPATIAL)
        @test dim.name == :x
        @test dim.size == 10
        @test dim.dim_type == DIM_SPATIAL

        # Test CoordinateSystem
        crs = CoordinateSystem(crs="EPSG:4326")
        @test !is_cartesian(crs)
        crs2 = CoordinateSystem()
        @test is_cartesian(crs2)

        # Test GeoArray
        data = rand(10, 5)
        dims = [Dimension(name=:x, size=10, coords=1:10), Dimension(name=:y, size=5, coords=1:5)]
        ga = GeoArray(data; dims=dims, crs=CoordinateSystem(), attrs=Dict{String, Any}("units" => "K"))
        @test size(ga) == (10, 5)
        @test ga.dims[1].name == :x

        # Test GeoDataset (simplified)
        ds = GeoDataset(
            variables=Dict("temp" => ga),
            coords=Dict(:x => GeoArray(collect(1.0:10.0); dims=[dims[1]])),
            dims=Dict(:x => dims[1], :y => dims[2])
        )
        @test hasproperty(ds, :temp)
        @test ds.temp === ga
    end

    @testset "Zarr Backend" begin
        # Create test data
        lon = collect(-70.0:1.0:-60.0)
        lat = collect(40.0:1.0:50.0)
        depth = collect(0.0:10.0:100.0)
        temp = rand(Float32, length(lon), length(lat), length(depth))

        g = zgroup(joinpath(test_dir, "test.zarr"); attrs=Dict("title" => "Test"))
        zcreate(Float32, g, "temperature", length(lon), length(lat), length(depth); chunks=(10, 10, 5), attrs=Dict("units" => "degC"))[:] = temp
        zcreate(Float64, g, "lon", length(lon); attrs=Dict("units" => "degrees_east"))[:] = lon
        zcreate(Float64, g, "lat", length(lat); attrs=Dict("units" => "degrees_north"))[:] = lat
        zcreate(Float64, g, "depth", length(depth); attrs=Dict("units" => "meters"))[:] = depth

        # Load
        ds = geoload(joinpath(test_dir, "test.zarr"))
        @test hasproperty(ds, :temperature)
        @test length(ds.dims) == 3
        @test haskey(ds.coords, :lon)
        @test haskey(ds.coords, :lat)
        @test haskey(ds.coords, :depth)

        # Slice
        ds_slice = geoslice(ds, lon=(-68, -65), lat=(42, 48))
        @test ds_slice.dims[:lon].size < ds.dims[:lon].size
        @test ds_slice.dims[:lat].size < ds.dims[:lat].size

        # Select
        ds_sel = geoselect(ds, lon=-65.0, lat=45.0)
        @test ds_sel.dims[:lon].size == 1
        @test ds_sel.dims[:lat].size == 1
        @test ds_sel.dims[:depth].size == ds.dims[:depth].size

        # Values at point
        vals = geovalues(ds, ["temperature"], lon=-65.0, lat=45.0, depth=50.0)
        @test haskey(vals, "temperature")

        # Stats
        st = variable_stats(ds, "temperature")
        @test st.count == length(temp)
        @test st.min <= st.max

        # Summary
        s = geosummary(ds)
        @test occursin("GeoDataset", s)

        # Bounding box
        bb = bounding_box(ds)
        @test haskey(bb, :lon)
        @test haskey(bb, :lat)

        # Cartesian check
        @test is_cartesian(ds.crs)
    end

    @testset "NetCDF Backend" begin
        ds = geoload(joinpath(test_dir, "test.zarr"))
        geosave(joinpath(test_dir, "test.nc"), ds; backend=:ncdatasets)
        ds2 = geoload(joinpath(test_dir, "test.nc"), backend=:ncdatasets)

        @test hasproperty(ds2, :temperature)
        @test size(ds2.temperature.data) == size(ds.temperature.data)
        @test haskey(ds2.coords, :lon)
        @test haskey(ds2.coords, :lat)
        @test haskey(ds2.coords, :depth)
    end

    @testset "Generic Dimensions (MCMC)" begin
        chain = collect(1:4)
        draw = collect(1:100)
        param = [:alpha, :beta, :sigma]
        mcmc_data = rand(Float64, length(chain), length(draw), length(param))

        g2 = zgroup(joinpath(test_dir, "mcmc.zarr"))
        zcreate(Float64, g2, "samples", length(chain), length(draw), length(param); chunks=(1, 20, 1))[:] = mcmc_data
        zcreate(Int, g2, "chain", length(chain))[:] = chain
        zcreate(Int, g2, "draw", length(draw))[:] = draw
        zcreate(Int, g2, "param", length(param); attrs=Dict("names" => ["alpha", "beta", "sigma"]))[:] = collect(1:length(param))

        ds_mcmc = geoload(joinpath(test_dir, "mcmc.zarr"))
        @test hasproperty(ds_mcmc, :samples)
        @test length(ds_mcmc.dims) == 3
        @test haskey(ds_mcmc.dims, :chain)
        @test haskey(ds_mcmc.dims, :draw)
        @test haskey(ds_mcmc.dims, :param)

        # Slice MCMC
        ds_slice = geoslice(ds_mcmc, chain=1:2, draw=1:50)
        @test ds_slice.dims[:chain].size == 2
        @test ds_slice.dims[:draw].size == 50
    end

    @testset "GeoParquet Backend" begin
        df = DataFrame(lon=[-65.0, -66.0], lat=[45.0, 46.0], temp=[15.0, 16.0])
        geoms = [GeoInterface.Point(-65.0, 45.0), GeoInterface.Point(-66.0, 46.0)]
        df[!, :geometry] = geoms
        GeoParquet.write(joinpath(test_dir, "test.parquet"), df, (:geometry,))

        ds_gp = geoload(joinpath(test_dir, "test.parquet"))
        @test hasproperty(ds_gp, :temp)
        @test haskey(ds_gp.coords, :lon)
        @test haskey(ds_gp.coords, :lat)
    end

    @testset "Join Operations" begin
        ds1 = geoload(joinpath(test_dir, "test.zarr"))
        ds2 = geoload(joinpath(test_dir, "test.zarr"))

        # geomerged
        ds_merged = geomerged([ds1, ds2])
        @test length(ds_merged.variables) == 2  # renamed

        # geojoin on time (not applicable here since no time dim, skip)
    end

    @testset "Index and Subset" begin
        ds = geoload(joinpath(test_dir, "test.zarr"))
        idx = geoindex(ds)
        @test hasproperty(idx, :bbox)

        ds_sub = geosubset(ds, bbox=(-68, -62, 42, 48))
        @test ds_sub.dims[:lon].size < ds.dims[:lon].size
    end

    @testset "Backend Registry" begin
        @test :zarr in list_backends()
        @test :ncdatasets in list_backends()
        @test :yaxarray in list_backends()
        @test :geoparquet in list_backends()
        @test :nczarr in list_backends()

        be = get_backend(:zarr)
        @test be isa ZarrBackend
    end

    @testset "Coordinate Utilities" begin
        # Identity mapping by default
        @test standardize_dimension_name(:longitude) == :longitude
        @test standardize_dimension_name(:x) == :x

        # With custom mapping
        my_mapping = Dict(:longitude => :lon, :latitude => :lat)
        @test standardize_dimension_name(:longitude; mapping=my_mapping) == :lon

        coords = collect(-70.0:1.0:-60.0)
        idx = find_coord_indices(coords, -65.0; mode=:nearest)
        @test coords[idx] ≈ -65.0 atol=1.0

        rng = slice_indices(coords, (-68.0, -62.0))
        @test rng[1] <= rng[end]

        @test is_regular_grid(coords)
        @test grid_spacing(coords) ≈ 1.0
    end

    # Clean up
    rm(test_dir, recursive=true, force=true)
end

println("All tests passed!")