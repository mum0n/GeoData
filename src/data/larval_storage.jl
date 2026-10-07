"""
    larval_storage.jl

Zarr-based analytical storage backend for larval particle tracking (*Chionoecetes opilio*).

Provides relational persistence for simulation runs, multi-million particle
trajectory steps, cohort recruitment outcomes, demographic connectivity
matrices, and spatial dispersal fields with fast array I/O, multi-scenario
benchmarking, and ensemble model averaging.

Uses Zarr directly for efficiency (native chunked/compressed arrays) and JSON
for metadata (simulation runs, recruitment metrics, connectivity).
"""
module LarvalStorage

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using Zarr
using JSON3
using Dates
using Statistics
using LinearAlgebra
using DataFrames

# Storage-specific exports
export
    # Core storage operations
    open_larval_storage,
    close_larval_storage,
    initialize_larval_storage_schema!,
    save_larval_simulation_run!,
    list_larval_simulation_runs,
    load_larval_trajectories,
    load_larval_recruitment_metrics,
    load_larval_connectivity,
    load_larval_gridded_dispersal,
    load_larval_hydrodynamic_fields,
    compare_larval_scenarios,
    compute_larval_ensemble_model_average

include("larval_storage_impl.jl")

end # module LarvalStorage