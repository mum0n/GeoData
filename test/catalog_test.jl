using GeoData
using GeoData.GeoDataCatalog
using Dates

println("Testing GeoData Catalog...")

# Create a test catalog
catalog_path = joinpath(pwd(), "test_catalog.json")
if isfile(catalog_path)
    rm(catalog_path)
end

catalog = GeoDataCatalog(catalog_path)

# Test creating a dataset entry
entry1 = DatasetEntry(
    name = "WOA23 Temperature Climatology",
    key = :woa23_temp,
    variables = ["temperature"],
    temporal_coverage = (DateTime(1955,1,1), DateTime(2018,12,31)),
    spatial_bounds = (-180.0, 180.0, -90.0, 90.0),
    vertical_range = (0.0, 5500.0),
    source = "https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/",
    location = "data/woa23/temperature.zarr",
    access_method = :file,
    format = :zarr,
    credits = "NOAA National Centers for Environmental Information",
    license = "Public Domain",
    notes = "World Ocean Atlas 2023 temperature climatology at 0.25 degree resolution"
)

# Register the dataset
uuid1 = register_dataset!(catalog, entry1)
println("Registered dataset with UUID: $uuid1")

# Test getting dataset by key
retrieved = get_dataset_by_key(catalog, :woa23_temp)
println("Retrieved dataset by key: $(retrieved.name)")

# Test getting dataset by UUID
retrieved2 = get_dataset(catalog, uuid1)
println("Retrieved dataset by UUID: $(retrieved2.name)")

# Test finding datasets
matches = find_datasets(catalog, variables = ["temperature"])
println("Found $(length(matches)) datasets with temperature variable")

# Test getting best dataset
best = get_best_dataset(catalog, variables = ["temperature"])
println("Best dataset for temperature: $(best.name)")

# Test access tracking
accessed = access_dataset!(catalog, uuid1)
println("Access count after access: $(accessed.access_count)")

# Test download recording
downloaded = record_download!(catalog, uuid1, 1024*1024*50)  # 50 MB
println("Total download size: $(downloaded.total_download_size) bytes")

# Test statistics
stats = get_stats(catalog)
println("Catalog stats: $(stats.total_entries) entries, $(stats.active_entries) active")

# Test export to CSV
csv_path = joinpath(pwd(), "test_catalog_export.csv")
export_catalog_to_csv(catalog, csv_path)
println("Exported catalog to CSV: $csv_path")

# Test saving and reloading
save_catalog(catalog)
println("Saved catalog to: $(catalog.catalog_file)")

# Create new catalog instance and load
catalog2 = GeoDataCatalog(catalog.catalog_file)
println("Reloaded catalog has $(length(catalog2.entries)) entries")

# Clean up
if isfile(catalog_path)
    rm(catalog_path)
end
if isfile(csv_path)
    rm(csv_path)
end

println("All tests passed!")