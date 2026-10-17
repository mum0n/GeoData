include("test_fixtures.jl")

@testset "Slicing, Selection & Operations" begin
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

    # Point value queries
    v = GeoData.geovalues(ds, ["temp"], lon = -64.0, lat = 44.0, depth = 10.0)
    @test v["temp"] isa Real

    v2 = GeoData.geovalues(ds, ["temp"], lon = [-64.0, -63.5], lat = [44.0, 44.5])
    @test v2["temp"] isa AbstractArray
    @test size(v2["temp"]) == (2, 2, 3)

    v3 = GeoData.geovalues(ds, ["temp"], lon = -64.0, lat = 44.0)
    @test v3["temp"] isa AbstractVector
    @test length(v3["temp"]) == 3

    @test_throws Exception GeoData.geovalues(ds, ["temp"])
    @test_throws Exception GeoData.geovalues(ds, ["temp"], lon = [-64.0, -63.5, -64.0],
                                             lat = [44.0, 44.5])
    @test_throws Exception GeoData.geovalues(ds, ["nope"], lon = -64.0)

    # Joining
    ds1 = make_dataset()
    ds2 = GeoDataset(
        Dict("sal" => GeoArray(rand(Float32, 5, 5, 3), ds1.temp.dims,
                               ds1.temp.crs, Dict{String, Any}())),
        ds1.coords, ds1.dims, ds1.crs, Dict{String, Any}(), nothing, "synthetic2")

    inner = GeoData.geojoin([ds1, ds2], on = [:lon, :lat, :depth], how = :inner)
    assert_dataset_invariants(inner)
    @test haskey(inner.variables, "temp") && haskey(inner.variables, "sal")
    @test inner.temp.dims[1].size == 5

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
    left = GeoData.geojoin([ds1, ds3], on = [:lon, :lat, :depth], how = :left)
    assert_dataset_invariants(left)
    @test left.dims[:lon].size == 5
    @test any(isnan, left.sal.data)
    @test !any(isnan, left.temp.data)

    @test_throws Exception GeoData.geojoin([ds1, ds2], on = [])
    @test_throws Exception GeoData.geojoin([ds1, ds2], on = [:lon], how = :cross)
    ds_dup = GeoDataset(Dict("temp" => ds1.temp), ds1.coords, ds1.dims, ds1.crs,
                        Dict{String, Any}(), nothing, "dup")
    @test_throws Exception GeoData.geojoin([ds1, ds_dup], on = [:lon, :lat, :depth])

    # Aggregation
    agg = GeoData.geoaggregate(ds, dim = :lat, func = mean)
    assert_dataset_invariants(agg)
    @test ndims(agg.temp) == 2
    @test agg.temp.dims[1].name == :lon
    @test !haskey(agg.dims, :lat)
    @test agg.dims[:lon].size == 5
    @test size(agg.temp.data) == (5, 3)
    @test haskey(agg.attrs, "reduced")

    surf = GeoDataset(Dict("temp" => ds.temp, "mask" => GeoArray(trues(5), (ds.temp.dims[1],),
                                                                ds.crs, Dict{String, Any}())),
                      ds.coords, ds.dims, ds.crs, Dict{String, Any}(), nothing, "s")
    agg2 = GeoData.geoaggregate(surf, dim = :lat, func = sum)
    @test haskey(agg2.variables, "mask")
    @test_throws Exception GeoData.geoaggregate(ds, dim = :time, func = mean)
end
