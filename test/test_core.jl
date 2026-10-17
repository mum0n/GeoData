include("test_fixtures.jl")

@testset "Core Types & Dataset Invariants" begin
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
    @test Base.IndexStyle(typeof(ga)) == Base.IndexStyle(Array)

    ds = make_dataset()
    assert_dataset_invariants(ds)

    dl = ds.temp.dims[1]
    undeclared = GeoDataset(
        Dict("temp" => GeoArray(rand(5, 3), (dl, Dimension(:time, 3)), ds.crs,
                                Dict{String, Any}())),
        Dict{Symbol, GeoArray}(), Dict(:lon => dl), ds.crs, Dict{String, Any}(),
        nothing, "undeclared")
    @test_throws Exception assert_dataset_invariants(undeclared)

    mismatched = GeoDataset(
        Dict("temp" => GeoArray(rand(5, 5, 3), ds.temp.dims, ds.crs, Dict{String, Any}())),
        Dict(:lon => GeoArray(collect(1.0:4.0), (dl,), ds.crs, Dict{String, Any}())),
        ds.dims, ds.crs, Dict{String, Any}(), nothing, "mismatched")
    @test_throws Exception assert_dataset_invariants(mismatched)
end
