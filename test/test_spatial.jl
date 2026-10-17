include("test_fixtures.jl")

@testset "Spatial Utilities & Geometry" begin
    coords = collect(-65.0:0.5:-63.0)
    idx = find_coord_indices(coords, -64.0; mode = :nearest)
    @test coords[idx] ≈ -64.0
    @test is_regular_grid(coords)
    @test grid_spacing(coords) ≈ 0.5
    @test standardize_dimension_name(:x) == :lon
    @test GeoData.normalize_depth([-10.0, -20.0])[1] == [10.0, 20.0]

    for desc in (false, true)
        ds = make_dataset(; descending_lat = desc)
        s = GeoData.geoslice(ds, lat = (44.0, 44.5))
        assert_dataset_invariants(s)
        @test s.coords[:lat].data == (desc ? [44.5, 44.0] : [44.0, 44.5])
    end

    # Bathymetry analysis and land exclusion
    lons = collect(1.0:10.0)
    lats = collect(1.0:10.0)
    elev = fill(-50.0, 10, 10)
    ds_b = make_2d_dataset(elev, lons, lats, "elevation")

    water_only = extract_marine_cells(ds_b)
    @test all(water_only.mask)

    ring = [(lon = 3.5, lat = 3.5), (lon = 6.5, lat = 3.5),
            (lon = 6.5, lat = 6.5), (lon = 3.5, lat = 6.5),
            (lon = 3.5, lat = 3.5)]
    rings = [(name = "island", code = :island, lons = [p.lon for p in ring],
              lats = [p.lat for p in ring])]

    with_land = extract_marine_cells(ds_b; coastlines = rings)
    @test with_land.mask[4, 4] == false
    @test with_land.mask[5, 5] == false
    @test with_land.mask[4, 6] == false
    @test with_land.mask[3, 3] == true
    @test count(with_land.mask) == count(water_only.mask) - 9
    @test length(with_land.indices) == 100 - 9

    open_ring = [(lons = [3.5, 6.5, 6.5, 3.5], lats = [3.5, 3.5, 6.5, 6.5],
                  name = "open", code = Symbol("open"))]
    open_test = extract_marine_cells(ds_b; coastlines = open_ring)
    @test open_test.mask == water_only.mask

    @test isempty(extract_marine_cells(ds_b; min_seabed_depth = 60.0).indices)

    pts = sample_marine_cells(25; bathymetry_ds = ds_b, coastlines = rings,
                               jitter_scale = 0.0)
    @test length(pts.lons) == 25
    @test all(1.0 .<= pts.lons .<= 10.0) && all(1.0 .<= pts.lats .<= 10.0)
    in_land = [(lon in (4.0, 5.0, 6.0)) && (lat in (4.0, 5.0, 6.0))
               for (lon, lat) in zip(pts.lons, pts.lats)]
    @test !any(in_land)
    @test all(pts.depths .<= -50.0)

    # Orientation and depth pairing
    lon2 = collect(-65.0:0.25:-61.0)
    lat2 = collect(43.0:0.25:45.0)
    n_lon, n_lat = length(lon2), length(lat2)

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
        @test length(mc.indices) == n_lon * n_lat
        @test size(mc.mask) == (n_lat, n_lon)
        @test all(answer.(mc.lons, mc.lats) .≈ mc.elevation)
        @test all(-65.0 .<= mc.lons .<= -61.0) && all(43.0 .<= mc.lats .<= 45.0)
        @test length(unique(mc.lons)) == n_lon && length(unique(mc.lats)) == n_lat

        s_sample = sample_marine_cells(40; bathymetry_ds = ds, jitter_scale = 0.0)
        @test all(answer.(s_sample.lons, s_sample.lats) .≈ s_sample.depths)
    end

    bad = make_3d_dataset(rand(n_lon, n_lat, 3), lon2, lat2,
                          collect(0.0:1.0:2.0), "elevation")
    @test_throws ErrorException extract_marine_cells(bad)

    dl = Dimension(name = :lon, size = n_lon, coords = lon2, dim_type = DIM_SPATIAL)
    dla = Dimension(name = :lat, size = n_lat, coords = lat2, dim_type = DIM_SPATIAL)
    crs = CoordinateSystem(crs = "EPSG:4326")
    mismatched = GeoDataset(
        Dict("elevation" => GeoArray(rand(n_lon + 2, n_lat), (dl, dla), crs,
                                     Dict{String, Any}())),
        Dict(:lon => GeoArray(lon2, (dl,), crs, Dict{String, Any}()),
             :lat => GeoArray(lat2, (dla,), crs, Dict{String, Any}())),
        Dict(:lon => dl, :lat => dla), crs, Dict{String, Any}(), nothing, "bad")
    @test_throws ErrorException extract_marine_cells(mismatched)

    # Migrated pure-computation functions
    rings_comp = [(name = "island", code = :island,
               lons = [0.0, 1.0, 1.0, 0.0, 0.0], lats = [0.0, 0.0, 1.0, 1.0, 0.0])]

    @test buffer_distance_to_degrees(111_000.0, 60.0)[1] >
          buffer_distance_to_degrees(111_000.0, 0.0)[1]

    (blon, blat) = expand_domain_with_buffer((-65.0, -63.0), (43.0, 45.0);
                                             buffer_km = 20.0)
    @test blon[1] < -65.0 && blat[2] > 45.0

    sq_lons = [0.0, 1.0, 1.0, 0.0]; sq_lats = [0.0, 0.0, 1.0, 1.0]
    @test point_in_polygon(0.5, 0.5, sq_lons, sq_lats)
    @test !point_in_polygon(2.0, 0.5, sq_lons, sq_lats)
    @test !point_in_polygon(1.0, 1.0, sq_lons, sq_lats)

    lon_reg = collect(-65.0:0.25:-63.0); lat_reg = collect(43.0:0.25:45.0)
    dep_reg = collect(0.0:10.0:50.0)
    tgt_lon = collect(-65.0:0.05:-63.8)
    tgt_lat = collect(43.0:0.05:43.8)

    src2 = make_2d_dataset(rand(9, 9), lon_reg, lat_reg, "ssh")
    @test size(regrid_2d_field(lon_reg, lat_reg, src2.variables["ssh"].data, tgt_lon, tgt_lat)) ==
          (25, 17)
    @test size(regrid_2d_field(src2, "ssh", tgt_lon, tgt_lat)) == (25, 17)

    src3 = make_3d_dataset(rand(9, 9, 6), lon_reg, lat_reg, dep_reg, "t")
    @test size(regrid_3d_field(src3, "t", tgt_lon, tgt_lat, dep_reg)) == (25, 17, 6)
    src3r = make_3d_dataset(rand(9, 9, 6), lon_reg, reverse(lat_reg), dep_reg, "t")
    @test size(regrid_3d_field(src3r, "t", tgt_lon, tgt_lat, dep_reg)) == (25, 17, 6)

    @test size(smooth_bathymetry(sign.(rand(12, 12) .- 0.5) .* 2.0; passes = 2)) ==
          (12, 12)

    it = get_bathymetry_interpolator(make_2d_dataset(fill(-100.0, 4, 4),
                                                     [1.0, 2.0, 3.0, 4.0],
                                                     [1.0, 2.0, 3.0, 4.0], "elevation"))
    @test it(2.0, 2.0) == -100.0

    @test isempty(load_coastline_polygons("C:/nope/nowhere.parquet"))

    scratch_geo = test_scratch()
    ring_lons_c = [-60.0, -60.0, -59.0, -59.0, -60.0]
    ring_lats_c = [44.0, 45.0, 45.0, 44.0, 44.0]
    line_c = GeoInterface.LineString(
        [(ring_lons_c[i], ring_lats_c[i]) for i in eachindex(ring_lons_c)])
    mls_c = GeoInterface.MultiLineString([
        [(-60.0, 44.0), (-59.0, 45.0)],
        [(-58.0, 44.0), (-57.0, 45.0)]])
    coast_path_c = joinpath(scratch_geo, "coastline.parquet")
    GeoParquet.write(coast_path_c,
                     DataFrame(id = [1, 2], name = ["test ring", "two parts"],
                               code = ["test_ring", "two_parts"],
                               geometry = [line_c, mls_c]), (:geometry,))

    read_rings_c = load_coastline_polygons(coast_path_c)
    @test length(read_rings_c) == 3
    @test read_rings_c[1].name == "test ring" && read_rings_c[1].code == :test_ring
    @test read_rings_c[1].lons == ring_lons_c
    @test read_rings_c[1].lats == ring_lats_c
    @test read_rings_c[2].lons == [-60.0, -59.0] && read_rings_c[2].lats == [44.0, 45.0]
    @test read_rings_c[3].lons == [-58.0, -57.0] && read_rings_c[3].lats == [44.0, 45.0]

    @test is_point_on_land_geodata(-59.5, 44.5; coastline = [read_rings_c[1]])
    @test !is_point_on_land_geodata(-58.5, 44.5; coastline = [read_rings_c[1]])
    @test !is_point_on_land_geodata(-58.5, 44.5; coastline_file = coast_path_c)

    @test is_point_on_land_geodata(0.5, 0.5; coastline = rings_comp)
    @test !is_point_on_land_geodata(9.9, 9.9; coastline = rings_comp)
    @test !is_point_on_land_geodata(9.9, 9.9)
    open_c = [(name = "o", code = :o, lons = [0.0, 1.0, 1.0, 0.0],
             lats = [0.0, 0.0, 1.0, 1.0])]
    @test !is_point_on_land_geodata(0.5, 0.5; coastline = open_c)

    bds_c = make_2d_dataset(fill(-100.0, 6, 6), collect(0.0:1.0:5.0),
                          collect(0.0:1.0:5.0), "elevation")
    @test is_marine_water_geodata(9.9, 9.9; coastlines = rings_comp)
    @test !is_marine_water_geodata(0.5, 0.5; coastlines = rings_comp)
    @test is_marine_water_geodata(2.5, 2.5; bathymetry = bds_c, min_seabed_depth = 50.0)
    @test !is_marine_water_geodata(2.5, 2.5; bathymetry = bds_c, min_seabed_depth = 200.0)

    nm, ga = variables_like(make_2d_dataset(rand(3, 3), [1.0, 2.0, 3.0],
                                            [1.0, 2.0, 3.0], "temp_C"), ["temp_C"])
    @test nm == "temp_C" && ga isa GeoArray

    v_norm, flipped = normalize_depth([-5.0, -10.0])
    @test v_norm == [5.0, 10.0] && flipped
    v2_norm, flipped2 = normalize_depth([5.0, 10.0])
    @test v2_norm == [5.0, 10.0] && !flipped2

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
