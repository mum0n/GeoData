"""
    larval_storage_impl.jl

Zarr-based storage implementation for larval particle tracking data.
"""

using GeoData
using GeoData.GeoDataCoreTypes: GeoDataset, GeoArray, Dimension, CoordinateSystem
using Zarr
using GeoData.Data.ZarrStorage
using JSON3
using Dates
using Statistics
using LinearAlgebra
using DataFrames

"""
    open_larval_storage(
        store_path::AbstractString = joinpath("outputs", "particle_tracking.zarr");
        mode::Symbol = :rw
    ) -> Zarr.ZGroup

Open or create a Zarr storage group for larval particle tracking data.

# Arguments
- `store_path`: Path to the Zarr directory store
- `mode`: :ro (read-only) or :rw (read-write)

# Returns
Zarr.ZGroup handle to the storage
"""
function open_larval_storage(
    store_path::AbstractString = joinpath("outputs", "particle_tracking.zarr");
    mode::Symbol = :rw
)
    mkpath(store_path)
    store = Zarr.DirectoryStore(store_path)
    return Zarr.open_group(store; mode = mode)
end

"""
    close_larval_storage(group::Zarr.ZGroup)

Close the Zarr storage group (flushes all data).
"""
function close_larval_storage(group::Zarr.ZGroup)
    Zarr.close(group.store)
    nothing
end

"""
    initialize_larval_storage_schema!(group::Zarr.ZGroup)

Initialize the Zarr storage schema with arrays for trajectories, metadata, etc.
"""
function initialize_larval_storage_schema!(group::Zarr.ZGroup)
    # Create metadata group
    if !haskey(group, "metadata")
        ZarrStorage.zarr_create_group(group, "metadata")
    end
    
    # Create trajectories array (chunked for efficiency)
    # Shape: (n_particles, n_timesteps) for each variable
    if !haskey(group, "trajectories")
        traj_group = ZarrStorage.zarr_create_group(group, "trajectories")
        # We'll create arrays on first write since we don't know dimensions upfront
    end
    
    # Create gridded fields group
    if !haskey(group, "gridded_fields")
        ZarrStorage.zarr_create_group(group, "gridded_fields")
    end
    
    # Create hydrodynamic fields group
    if !haskey(group, "hydrodynamic_fields")
        ZarrStorage.zarr_create_group(group, "hydrodynamic_fields")
    end
    
    # Create runs metadata as JSON
    if !haskey(group, "runs")
        ZarrStorage.create_zarr_array(
            group, "runs";
            shape = (0,),
            chunks = (1000,),
            dtype = String,
            compressor = Zarr.Blosc(cname = "zstd", clevel = 3)
        )
    end
    
    nothing
end

"""
    save_larval_simulation_run!(
        group::Zarr.ZGroup,
        run_id::AbstractString,
        opts;
        trajectories::NamedTuple,
        metrics::Union{Nothing, NamedTuple} = nothing,
        connectivity::Union{Nothing, NamedTuple} = nothing,
        gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
        config::Union{Nothing, AbstractDict, AbstractString} = nothing,
        notes::AbstractString = ""
    )::String

Save a complete simulation run to Zarr storage.
"""
function save_larval_simulation_run!(
    group::Zarr.ZGroup,
    run_id::AbstractString,
    opts;
    trajectories::NamedTuple,
    metrics::Union{Nothing, NamedTuple} = nothing,
    connectivity::Union{Nothing, NamedTuple} = nothing,
    gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
    config::Union{Nothing, AbstractDict, AbstractString} = nothing,
    notes::AbstractString = ""
)
    initialize_larval_storage_schema!(group)
    
    # Canonicalize trajectories
    canonical_trajs = canonicalize_trajectories(trajectories)
    
    # Extract metadata
    scenario = string(hasproperty(opts, :scenario) ? opts.scenario : :baseline)
    proj_year = Int(hasproperty(opts, :projection_year) ? opts.projection_year : 2050)
    n_parts = Int(hasproperty(opts, :n_particles) ? opts.n_particles : size(canonical_trajs.lons, 1))
    duration = Float64(hasproperty(opts, :track_duration) ? opts.track_duration : (canonical_trajs.times[end] - canonical_trajs.times[1]))
    dt_val = Float64(hasproperty(opts, :track_dt) ? opts.track_dt : (length(canonical_trajs.times) > 1 ? canonical_trajs.times[2] - canonical_trajs.times[1] : 300.0))
    tides = Bool(hasproperty(opts, :enable_tides) ? opts.enable_tides : false)
    dvm = Bool(hasproperty(opts, :enable_dvm) ? opts.enable_dvm : true)
    molt = Bool(hasproperty(opts, :enable_molting) ? opts.enable_molting : true)
    diff_h = Float64(hasproperty(opts, :diffusivity_h) ? opts.diffusivity_h : 10.0)
    diff_v = Float64(hasproperty(opts, :diffusivity_v) ? opts.diffusivity_v : 1e-4)
    min_depth = Float64(hasproperty(opts, :min_seabed_depth) ? opts.min_seabed_depth : 100.0)
    seed_val = Int(hasproperty(opts, :seed) ? opts.seed : 42)
    created_at = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    
    # Serialize configuration
    config_toml_str = if !isnothing(config)
        if config isa AbstractString
            config
        else
            s_io = IOBuffer()
            TOML.print(s_io, config; sorted = true)
            String(take!(s_io))
        end
    elseif isdefined(Main, :options_to_configuration) && opts isa HydrodynamicOptions
        s_io = IOBuffer()
        TOML.print(s_io, options_to_configuration(opts); sorted = true)
        String(take!(s_io))
    else
        ""
    end
    
    # 1. Save run metadata as JSON
    run_metadata = Dict(
        "run_id" => run_id,
        "scenario" => scenario,
        "projection_year" => proj_year,
        "created_at" => created_at,
        "n_particles" => n_parts,
        "duration_seconds" => duration,
        "dt_seconds" => dt_val,
        "enable_tides" => tides,
        "enable_dvm" => dvm,
        "enable_molting" => molt,
        "diffusivity_h" => diff_h,
        "diffusivity_v" => diff_v,
        "min_seabed_depth" => min_depth,
        "seed" => seed_val,
        "config_toml" => config_toml_str,
        "notes" => notes
    )
    
    # Save run metadata
    runs_array = group["runs"]
    # Append to runs array (Zarr doesn't have native append, so we rewrite)
    existing_runs = Zarr.read(runs_array)
    push!(existing_runs, JSON3.write(run_metadata))
    Zarr.write(runs_array, existing_runs)
    
    # 2. Save trajectories as Zarr arrays
    traj_group = group["trajectories"]
    
    n_p, n_t = size(canonical_trajs.lons)
    
    # Create arrays if they don't exist
    for (name, data) in (
        "lons" => canonical_trajs.lons,
        "lats" => canonical_trajs.lats,
        "depths" => canonical_trajs.depths,
        "temperatures" => canonical_trajs.temperatures,
        "degree_days" => canonical_trajs.degree_days,
        "degree_days_timeseries" => canonical_trajs.degree_days_timeseries,
        "survival_probability" => canonical_trajs.survival_probability,
        "stages" => canonical_trajs.stages,
        "alive" => canonical_trajs.alive,
        "settlement_status" => canonical_trajs.settlement_status,
        "settlement_age" => get(canonical_trajs, :settlement_age, nothing),
        "ids" => canonical_trajs.ids,
        "times" => canonical_trajs.times
    )
        if !isnothing(data)
            if !haskey(traj_group, name)
                # Create array with compression
                Zarr.create_array(
                    traj_group, name;
                    shape = size(data),
                    chunks = (min(1000, size(data, 1)), size(data, 2)),
                    dtype = eltype(data),
                    compressor = Zarr.Blosc(cname = "zstd", clevel = 3)
                )
            end
            # Write data - for variable-length arrays we store as JSON
            if name in ["ids", "times", "stages", "settlement_status"]
                Zarr.write(traj_group[name], JSON3.write(data))
            else
                Zarr.write(traj_group[name], data)
            end
        end
    end
    
    # 3. Save recruitment metrics
    if !isnothing(metrics)
        metrics_group = Zarr.create_group(group, "metrics")
        if !haskey(metrics_group, run_id)
            Zarr.create_group(metrics_group, run_id)
        end
        for (k, v) in pairs(metrics)
            Zarr.write(metrics_group[run_id][k], v)
        end
    end
    
    # 4. Save connectivity
    if !isnothing(connectivity)
        conn_group = Zarr.create_group(group, "connectivity")
        if !haskey(conn_group, run_id)
            Zarr.create_group(conn_group, run_id)
        end
        Zarr.write(conn_group[run_id]["matrix"], connectivity.matrix)
        Zarr.write(conn_group[run_id]["strata_names"], connectivity.strata_names)
        if hasproperty(connectivity, :counts_unweighted)
            Zarr.write(conn_group[run_id]["counts_unweighted"], connectivity.counts_unweighted)
        end
        if hasproperty(connectivity, :counts_matrix)
            Zarr.write(conn_group[run_id]["counts_matrix"], connectivity.counts_matrix)
        end
    end
    
    # 5. Save gridded dispersal
    if !isnothing(gridded_dispersal)
        grid_group = Zarr.create_group(group, "gridded_dispersal")
        if !haskey(grid_group, run_id)
            Zarr.create_group(grid_group, run_id)
        end
        for (k, v) in pairs(gridded_dispersal)
            if v isa Matrix || v isa Vector
                Zarr.write(grid_group[run_id][k], v)
            end
        end
    end
    
    return run_id
end

"""
    list_larval_simulation_runs(
        group::Zarr.ZGroup;
        scenario::Union{Nothing, AbstractString, Symbol} = nothing,
        projection_year::Union{Nothing, Int} = nothing
    ) -> DataFrame

Query and return a summary DataFrame of all simulation runs.
"""
function list_larval_simulation_runs(
    group::Zarr.ZGroup;
    scenario::Union{Nothing, AbstractString, Symbol} = nothing,
    projection_year::Union{Nothing, Int} = nothing
)::DataFrame
    runs_array = group["runs"]
    runs_json = Zarr.read(runs_array)
    
    runs = DataFrame()
    for r_json in runs_json
        r = JSON3.read(r_json)
        if !isnothing(scenario) && r["scenario"] != string(scenario)
            continue
        end
        if !isnothing(projection_year) && r["projection_year"] != projection_year
            continue
        end
        push!(runs, r)
    end
    
    sort!(runs, :created_at, rev = true)
    return runs
end

"""
    load_larval_trajectories(
        group::Zarr.ZGroup,
        run_id::AbstractString;
        particle_ids::Union{Nothing, AbstractVector{Int}} = nothing,
        stage::Union{Nothing, Symbol, AbstractString} = nothing,
        time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
        max_particles::Union{Nothing, Int} = nothing
    ) -> NamedTuple

Load particle trajectories from Zarr storage.
"""
function load_larval_trajectories(
    group::Zarr.ZGroup,
    run_id::AbstractString;
    particle_ids::Union{Nothing, AbstractVector{Int}} = nothing,
    stage::Union{Nothing, Symbol, AbstractString} = nothing,
    time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
    max_particles::Union{Nothing, Int} = nothing
)::NamedTuple
    
    traj_group = group["trajectories"]
    
    # Read all trajectory arrays
    lons = Zarr.read(traj_group["lons"])
    lats = Zarr.read(traj_group["lats"])
    depths = Zarr.read(traj_group["depths"])
    temperatures = Zarr.read(traj_group["temperatures"])
    degree_days = Zarr.read(traj_group["degree_days"])
    degree_days_ts = Zarr.read(traj_group["degree_days_timeseries"])
    survival_prob = Zarr.read(traj_group["survival_probability"])
    stages_json = JSON3.read(Zarr.read(traj_group["stages"]))
    alive = Zarr.read(traj_group["alive"])
    settle_status_json = JSON3.read(Zarr.read(traj_group["settlement_status"]))
    ids = JSON3.read(Zarr.read(traj_group["ids"]))
    times = JSON3.read(Zarr.read(traj_group["times"]))
    
    # Apply filters
    n_p = size(lons, 1)
    p_ids = collect(1:n_p)
    
    if !isnothing(particle_ids) && !isempty(particle_ids)
        p_ids = particle_ids
    elseif !isnothing(max_particles) && max_particles > 0
        p_ids = p_ids[1:min(max_particles, n_p)]
    end
    
    # Return filtered trajectories
    return (
        lons = lons[p_ids, :],
        lats = lats[p_ids, :],
        depths = depths[p_ids, :],
        temperatures = temperatures[p_ids, :],
        degree_days = degree_days[p_ids, :],
        degree_days_timeseries = degree_days_ts[p_ids, :],
        survival_probability = survival_prob[p_ids, :],
        stages = stages_json[p_ids, :],
        alive = alive[p_ids, :],
        settlement_status = settle_status_json[p_ids, :],
        settlement_age = haskey(traj_group, "settlement_age") ? JSON3.read(Zarr.read(traj_group["settlement_age"]))[p_ids, :] : nothing,
        times = times,
        ids = ids[p_ids]
    )
end

"""
    load_larval_recruitment_metrics(group::Zarr.ZGroup, run_id::AbstractString) -> NamedTuple

Load recruitment metrics for a run.
"""
function load_larval_recruitment_metrics(group::Zarr.ZGroup, run_id::AbstractString)
    metrics_group = group["metrics"]
    if haskey(metrics_group, run_id)
        run_group = metrics_group[run_id]
        return (
            total_released = Zarr.read(run_group["total_released"]),
            total_settled_successful = Zarr.read(run_group["total_settled_successful"]),
            total_settled_unsuitable = Zarr.read(run_group["total_settled_unsuitable"]),
            total_dead_thermal = Zarr.read(run_group["total_dead_thermal"]),
            total_pelagic_remaining = Zarr.read(run_group["total_pelagic_remaining"]),
            settlement_success_rate = Zarr.read(run_group["settlement_success_rate"]),
            mean_pld_days = Zarr.read(run_group["mean_pld_days"]),
            mean_degree_days = Zarr.read(run_group["mean_degree_days"]),
            mean_exposure_temperature = Zarr.read(run_group["mean_exposure_temperature"]),
            mean_dispersal_distance_km = Zarr.read(run_group["mean_dispersal_distance_km"])
        )
    end
    return nothing
end

"""
    load_larval_connectivity(group::Zarr.ZGroup, run_id::AbstractString) -> NamedTuple

Load connectivity matrix for a run.
"""
function load_larval_connectivity(group::Zarr.ZGroup, run_id::AbstractString)
    conn_group = group["connectivity"]
    if haskey(conn_group, run_id)
        run_group = conn_group[run_id]
        return (
            matrix = Zarr.read(run_group["matrix"]),
            strata_names = Zarr.read(run_group["strata_names"]),
            counts_unweighted = haskey(run_group, "counts_unweighted") ? Zarr.read(run_group["counts_unweighted"]) : nothing,
            counts_matrix = haskey(run_group, "counts_matrix") ? Zarr.read(run_group["counts_matrix"]) : nothing
        )
    end
    return nothing
end

"""
    load_larval_gridded_dispersal(group::Zarr.ZGroup, run_id::AbstractString) -> NamedTuple

Load gridded dispersal fields for a run.
"""
function load_larval_gridded_dispersal(group::Zarr.ZGroup, run_id::AbstractString)
    grid_group = group["gridded_dispersal"]
    if haskey(grid_group, run_id)
        run_group = grid_group[run_id]
        fields = NamedTuple()
        for k in keys(run_group)
            fields = merge(fields, (Symbol(k) => Zarr.read(run_group[k]),))
        end
        return fields
    end
    return nothing
end

"""
    load_larval_hydrodynamic_fields(group::Zarr.ZGroup, run_id::AbstractString) -> NamedTuple

Load hydrodynamic fields for a run.
"""
function load_larval_hydrodynamic_fields(group::Zarr.ZGroup, run_id::AbstractString)
    hydro_group = group["hydrodynamic_fields"]
    if haskey(hydro_group, run_id)
        run_group = hydro_group[run_id]
        fields = NamedTuple()
        for k in keys(run_group)
            fields = merge(fields, (Symbol(k) => Zarr.read(run_group[k]),))
        end
        return fields
    end
    return nothing
end

"""
    compare_larval_scenarios(
        group::Zarr.ZGroup;
        scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
        projection_year::Int = 2050
    ) -> DataFrame

Compare scenarios across runs.
"""
function compare_larval_scenarios(
    group::Zarr.ZGroup;
    scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
    projection_year::Int = 2050
)::DataFrame
    runs = list_larval_simulation_runs(group; scenario = nothing, projection_year = projection_year)
    
    # Filter by scenarios
    filter!(r -> r.scenario in [string(s) for s in scenarios], runs)
    
    return runs
end

"""
    compute_ensemble_model_average(
        group::Zarr.ZGroup;
        variable::Symbol = :mean_pld_days,
        scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
        projection_year::Int = 2050
    ) -> NamedTuple

Compute ensemble model average across scenarios.
"""
function compute_larval_ensemble_model_average(
    group::Zarr.ZGroup;
    variable::Symbol = :mean_pld_days,
    scenarios::Vector{Symbol} = [:baseline, :ssp245, :ssp585],
    projection_year::Int = 2050
)
    runs = list_larval_simulation_runs(group; scenario = nothing, projection_year = projection_year)
    filter!(r -> r.scenario in [string(s) for s in scenarios], runs)
    filter!(r -> r.scenario in [string(s) for s in scenarios], runs)
    
    if isempty(runs)
        return (mean = NaN, std = NaN, n = 0)
    end
    
    values = runs[!, variable]
    valid_values = values[.!isnan.(values)]
    
    return (
        mean = mean(valid_values),
        std = std(valid_values),
        n = length(valid_values)
    )
end

"""
    canonicalize_trajectories(trajectories) -> NamedTuple

Ensure trajectories have consistent structure.
"""
function canonicalize_trajectories(trajectories::NamedTuple)
    # Ensure all required fields exist with proper types
    required_fields = (:lons, :lats, :depths, :temperatures, :degree_days, 
                       :degree_days_timeseries, :survival_probability, :stages, 
                       :alive, :settlement_status, :times, :ids)
    
    result = Dict{Symbol, Any}()
    for field in required_fields
        if hasproperty(trajectories, field)
            result[field] = getproperty(trajectories, field)
        else
            # Provide defaults
            if field == :stages
                result[field] = fill(:unknown, size(trajectories.lons))
            elseif field == :alive
                result[field] = fill(true, size(trajectories.lons))
            elseif field == :settlement_status
                result[field] = fill(:pelagic, size(trajectories.lons))
            elseif field == :degree_days_timeseries
                result[field] = zeros(size(trajectories.lons))
            elseif field == :ids
                n_p = size(trajectories.lons, 1)
                result[field] = collect(1:n_p)
            end
        end
    end
    
    return NamedTuple(result)
end