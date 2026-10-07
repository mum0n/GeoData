"""
Generic backend-agnostic operations on GeoDataset.

These functions work on the abstract GeoDataset container regardless of
which backend produced it. Backends may override individual operations
(`backend_slice`, `backend_select`, ...) for lazy/push-down evaluation;
the default implementation lives here.
"""
module GeoDataOperations

using ..GeoDataCoreTypes
using ..GeoDataCoordinates
using ..GeoDataInterfaces
using LinearAlgebra
using Statistics
using SparseArrays

export
    slice,
    select_values,
    values_at,
    build_index,
    join_datasets,
    regrid,
    aggregate,
    geosummary,
    data_stats,
    variable_stats,
    bounding_box,
    _apply_slice,
    _gather_dim,
    _intersect_sorted,
    _union_sorted,
    _infer_dim_names

# ============================================================
# Slicing
# ============================================================

"""
    slice(ds::GeoDataset; kwargs...) -> GeoDataset

Return a sliced view by coordinate ranges.

# Keywords (per dimension)
- `lon = (min, max)` — longitude range
- `lat = (min, max)` — latitude range
- `depth = (min, max)` — depth range
- `time = (min, max)` — time range

A keyword may also be a single value (nearest point, drops the dimension)
or a vector of values (select those points).
"""
function slice(ds::GeoDataset; kwargs...)
    ranges = _parse_dim_specs(ds, kwargs)

    new_vars = Dict{String, GeoArray}()
    for (name, ga) in ds.variables
        nd, new_dims = _apply_slice(ga, ranges)
        new_vars[name] = GeoArray(nd, new_dims, ga.crs, ga.attrs)
    end

    new_coords = Dict{Symbol, GeoArray}()
    new_dims = Dict{Symbol, Dimension}()
    for (dim, ga) in ds.coords
        if haskey(ranges, dim)
            spec = ranges[dim]
            c = vec(ga.data)
            if spec isa UnitRange
                nc = c[spec]
                new_coords[dim] = _slice_coord(dim, ga, nc)
                new_dims[dim] = Dimension(name=dim, size=length(nc), coords=nc,
                                          units=ga.dims[1].units,
                                          standard_name=ga.dims[1].standard_name,
                                          calendar=ga.dims[1].calendar)
            elseif spec isa Int
                new_coords[dim] = _slice_coord(dim, ga, c[spec:spec])
                new_dims[dim] = Dimension(name=dim, size=1, coords=c[spec:spec],
                                          units=ga.dims[1].units,
                                          standard_name=ga.dims[1].standard_name,
                                          calendar=ga.dims[1].calendar)
            elseif spec isa AbstractVector
                nc = c[spec]
                new_coords[dim] = _slice_coord(dim, ga, nc)
                new_dims[dim] = Dimension(name=dim, size=length(nc), coords=nc,
                                          units=ga.dims[1].units,
                                          standard_name=ga.dims[1].standard_name,
                                          calendar=ga.dims[1].calendar)
            end
        else
            new_coords[dim] = ga
            new_dims[dim] = ds.dims[dim]
        end
    end

    return GeoDataset(new_vars, new_coords, new_dims, ds.crs, ds.attrs, ds.backend, ds.source)
end

function _parse_dim_specs(ds, kwargs)
    ranges = Dict{Symbol, Any}()
    for (k, v) in kwargs
        dim = standardize_dimension_name(Symbol(k))
        if !haskey(ds.coords, dim)
            ranges[dim] = v
            continue
        end
        c = vec(ds.coords[dim].data)
        if v isa Tuple && length(v) == 2
            ranges[dim] = slice_indices(c, Float64.(v))
        elseif v isa AbstractVector
            ranges[dim] = slice_indices(c, Float64.(v))
        else
            ranges[dim] = find_coord_indices(c, Float64(v))
        end
    end
    return ranges
end

function _slice_coord(dim, ga, nc)
    return GeoArray(nc,
        (Dimension(name=dim, size=length(nc), coords=nc,
                   units=ga.dims[1].units,
                   standard_name=ga.dims[1].standard_name,
                   calendar=ga.dims[1].calendar),),
        ga.crs, ga.attrs)
end

function _apply_slice(ga::GeoArray, ranges::Dict{Symbol, Any})
    inds = Any[]
    new_dims = Dimension[]
    for dim in ga.dims
        if haskey(ranges, dim.name)
            r = ranges[dim.name]
            if r isa Integer
                push!(inds, Int(r))
            elseif r isa AbstractVector
                push!(inds, r)
                coords = isnothing(dim.coords) ? nothing : dim.coords[r]
                push!(new_dims, Dimension(name=dim.name, size=length(r), coords=coords,
                                          units=dim.units, standard_name=dim.standard_name,
                                          calendar=dim.calendar, is_unlimited=dim.is_unlimited))
            else
                push!(inds, r)
                coords = isnothing(dim.coords) ? nothing : dim.coords[r]
                push!(new_dims, Dimension(name=dim.name, size=length(r), coords=coords,
                                          units=dim.units, standard_name=dim.standard_name,
                                          calendar=dim.calendar, is_unlimited=dim.is_unlimited))
            end
        else
            push!(inds, :)
            push!(new_dims, dim)
        end
    end
    return ga.data[inds...], tuple(new_dims...)
end

# ============================================================
# Selection (exact/nearest values)
# ============================================================

"""
    select_values(ds::GeoDataset; kwargs...) -> GeoDataset

Select specific coordinate values (nearest or exact match).
Like `slice` but intended for discrete value selection.
"""
function select_values(ds::GeoDataset; kwargs...)
    slice(ds; kwargs...)
end

# ============================================================
# Point value queries
# ============================================================

"""
    values_at(ds::GeoDataset, varnames::Vector{String}; kwargs...) -> Dict

Query values at specific coordinate points.

# Arguments
- `ds`: Source dataset
- `varnames`: Variables to query
- `kwargs`: Coordinate points as `lon = -65.0, lat = 45.0, depth = 50.0, time = ...`

# Returns
Dict mapping variable names to queried values (scalars or arrays).
"""
function values_at(ds::GeoDataset, varnames::Vector{String}; kwargs...)
    idx = Dict{Symbol, Vector{Int}}()
    for (k, v) in kwargs
        dim = standardize_dimension_name(Symbol(k))
        haskey(ds.coords, dim) || error("Unknown dimension: $dim")
        c = vec(ds.coords[dim].data)
        if v isa AbstractVector
            idx[dim] = [find_coord_indices(c, Float64(x)) for x in v]
        else
            idx[dim] = [find_coord_indices(c, Float64(v))]
        end
    end

    result = Dict{String, Any}()
    for name in varnames
        haskey(ds.variables, name) || error("Unknown variable: $name")
        ga = ds.variables[name]
        inds = Any[]
        for dim in ga.dims
            if haskey(idx, dim.name)
                push!(inds, idx[dim.name])
            else
                push!(inds, :)
            end
        end
        result[name] = ga.data[inds...]
    end
    return result
end

# ============================================================
# Indexing
# ============================================================

"""
    GeoIndex

Lightweight spatial/temporal index for fast subsetting.
"""
struct GeoIndex
    coords::Dict{Symbol, Vector{Float64}}
    bbox::NamedTuple
end

"""
    build_index(ds::GeoDataset; spatial::Bool=true, temporal::Bool=true) -> GeoIndex

Build a coordinate index for fast spatial/temporal queries.
"""
function build_index(ds::GeoDataset; spatial::Bool=true, temporal::Bool=true)
    coords = Dict{Symbol, Vector{Float64}}()
    for dim in (:lon, :lat, :depth, :time)
        if haskey(ds.coords, dim)
            coords[dim] = vec(Float64.(ds.coords[dim].data))
        end
    end
    bbox = bounding_box(ds)
    return GeoIndex(coords, bbox)
end

# ============================================================
# Joining
# ============================================================

function _intersect_sorted(vs::Vector{Float64}...)
    isempty(vs) && return Float64[]
    result = vs[1]
    for v in vs[2:end]
        result = intersect(result, v)
        isempty(result) && break
    end
    return sort(result)
end

function _union_sorted(vs::Vector{Float64}...)
    isempty(vs) && return Float64[]
    result = vcat(vs...)
    return sort(unique(result))
end

"""
    join_datasets(dss::Vector{GeoDataset}; on::Vector{Symbol}, how::Symbol=:inner) -> GeoDataset

Join multiple datasets along shared dimensions.

# Arguments
- `dss`: Datasets to join
- `on`: Dimension names to join on (e.g., `[:time]`, `[:lon, :lat]`)
- `how`: Join type - `:inner`, `:outer`, `:left`, `:right`

# Returns
Joined `GeoDataset`.
"""
function join_datasets(dss::Vector{GeoDataset}; on::Vector{Symbol}, how::Symbol=:inner)
    length(dss) >= 2 || error("join requires at least 2 datasets")

    # Compute target coords on `on` dims
    target = Dict{Symbol, Vector{Float64}}()
    for dim in on
        all(haskey(ds.coords, dim) for ds in dss) || error("All datasets must have dimension $dim to join on it")
        allvecs = [vec(Float64.(ds.coords[dim].data)) for ds in dss]
        if how === :inner
            target[dim] = _intersect_sorted(allvecs...)
        elseif how === :outer
            target[dim] = _union_sorted(allvecs...)
        elseif how === :left
            target[dim] = allvecs[1]
        elseif how === :right
            target[dim] = allvecs[end]
        else
            error("Unknown join type: $how")
        end
    end

    # Build new dataset: for each ds, reindex variables to target coords
    all_vars = Dict{String, GeoArray}()
    all_coords = Dict{Symbol, GeoArray}()
    all_dims = Dict{Symbol, Dimension}()

    # Set up coords for `on` dims
    for dim in on
        tc = target[dim]
        all_coords[dim] = _make_coord(dim, tc, dss[1])
        all_dims[dim] = Dimension(name=dim, size=length(tc), coords=tc)
    end

    # Carry over non-joined dims from first ds
    for (dim, d) in dss[1].dims
        if !(dim in on)
            all_dims[dim] = d
            all_coords[dim] = dss[1].coords[dim]
        end
    end

    # For each ds, for each variable, gather to target coords
    for (i, ds) in enumerate(dss)
        for (name, ga) in ds.variables
            new_data = _gather_to_target(ga, target, on, how)
            new_name = length(dss) == 1 ? name : "$(name)_ds$i"
            all_vars[new_name] = GeoArray(new_data, ga.dims, ga.crs, ga.attrs)
        end
    end

    return GeoDataset(all_vars, all_coords, all_dims, dss[1].crs, dss[1].attrs, dss[1].backend, join([ds.source for ds in dss], "+"))
end

function _make_coord(dim, vals, ds)
    ref = ds.coords[dim]
    return GeoArray(vals,
        (Dimension(name=dim, size=length(vals), coords=vals,
                   units=ref.dims[1].units,
                   standard_name=ref.dims[1].standard_name,
                   calendar=ref.dims[1].calendar),),
        ref.crs, ref.attrs)
end

function _gather_to_target(ga::GeoArray, target::Dict{Symbol, Vector{Float64}}, on::Vector{Symbol}, how::Symbol)
    inds = Any[]
    for dim in ga.dims
        if dim.name in on && haskey(target, dim.name)
            c = isnothing(dim.coords) ? Float64[] : dim.coords
            tc = target[dim.name]
            if isempty(c)
                push!(inds, 1:length(tc))
            else
                idxs = [find_coord_indices(c, Float64(t)) for t in tc]
                push!(inds, idxs)
            end
        else
            push!(inds, :)
        end
    end
    return ga.data[inds...]
end

# ============================================================
# Regridding
# ============================================================

"""
    regrid(ds::GeoDataset, target::GeoDataset; method::Symbol=:bilinear, vars::Vector{String}=String[]) -> GeoDataset

Regrid dataset to a target grid.

# Arguments
- `ds`: Source dataset
- `target`: Dataset defining target coordinates
- `method`: Interpolation method (`:bilinear`, `:nearest`)
- `vars`: Variables to regrid (empty = all data variables)

# Returns
Regridded `GeoDataset`.
"""
function regrid(ds::GeoDataset, target::GeoDataset; method::Symbol=:bilinear, vars::Vector{String}=String[])
    isempty(vars) && (vars = collect(keys(ds.variables)))

    target_coords = Dict{Symbol, Vector{Float64}}()
    for dim in (:lon, :lat, :depth, :time)
        if haskey(target.coords, dim)
            target_coords[dim] = vec(Float64.(target.coords[dim].data))
        end
    end

    new_vars = Dict{String, GeoArray}()
    for name in vars
        ga = ds.variables[name]
        new_data = _regrid_array(ga, target_coords, method)
        new_dims = _build_dims_from_target(ga, target_coords)
        new_vars[name] = GeoArray(new_data, new_dims, target.crs, ga.attrs)
    end

    return GeoDataset(new_vars, target.coords, target.dims, target.crs, ds.attrs, ds.backend, ds.source)
end

function _build_dims_from_target(ga::GeoArray, target_coords::Dict{Symbol, Vector{Float64}})
    new_dims = Dimension[]
    for dim in ga.dims
        if haskey(target_coords, dim.name)
            tc = target_coords[dim.name]
            push!(new_dims, Dimension(name=dim.name, size=length(tc), coords=tc,
                                      units=dim.units, standard_name=dim.standard_name,
                                      calendar=dim.calendar, is_unlimited=dim.is_unlimited))
        else
            push!(new_dims, dim)
        end
    end
    return tuple(new_dims...)
end

function _regrid_array(ga::GeoArray, target_coords::Dict{Symbol, Vector{Float64}}, method::Symbol)
    data = ga.data
    for (i, dim) in enumerate(ga.dims)
        if haskey(target_coords, dim.name) && !isnothing(dim.coords) && !isempty(dim.coords)
            src = vec(Float64.(dim.coords))
            tgt = target_coords[dim.name]
            if method === :nearest
                data = _resample_nearest(data, i, src, tgt)
            else
                data = _resample_linear(data, i, src, tgt)
            end
        end
    end
    return data
end

function _resample_linear(data, dim::Int, src::Vector{Float64}, tgt::Vector{Float64})
    d = ndims(data)
    perm = (dim, setdiff(1:d, dim)...)
    dp = permutedims(data, perm)
    nout = length(tgt)
    tail_size = size(dp)[2:end]
    out = Array{eltype(data)}(undef, nout, tail_size...)
    for I in Iterators.product(ntuple(_ -> 1, d-1)...)
        idxs = Tuple([1] + [I...])
        line = view(dp, :, ntuple(k -> idxs[k+1], d-1)...)
        out[:, ntuple(k -> idxs[k+1], d-1)...] = _interp1d(line, src, tgt)
    end
    invperm = zeros(Int, d)
    for (i, p) in enumerate(perm)
        invperm[p] = i
    end
    return permutedims(out, invperm)
end

function _resample_nearest(data, dim::Int, src::Vector{Float64}, tgt::Vector{Float64})
    idxs = [find_coord_indices(src, t) for t in tgt]
    return _gather_dim(data, dim, idxs)
end

function _interp1d(line::AbstractVector, src::Vector{Float64}, tgt::Vector{Float64})
    n = length(src)
    out = similar(line, length(tgt))
    for (j, t) in enumerate(tgt)
        if t <= src[1]
            out[j] = line[1]
        elseif t >= src[end]
            out[j] = line[end]
        else
            k = searchsortedlast(src, t)
            k = clamp(k, 1, n-1)
            dx = src[k+1] - src[k]
            sx = dx == 0 ? 0.0 : (t - src[k]) / dx
            out[j] = line[k] * (1-sx) + line[k+1] * sx
        end
    end
    return out
end

function _gather_dim(data, dim::Int, idxs::Vector{Int})
    d = ndims(data)
    perm = (dim, setdiff(1:d, dim)...)
    dp = permutedims(data, perm)
    sub = dp[idxs, ntuple(_ -> (:), d-1)...]
    invperm = zeros(Int, d)
    for (i, p) in enumerate(perm)
        invperm[p] = i
    end
    return permutedims(sub, invperm)
end

# ============================================================
# Aggregation
# ============================================================

"""
    aggregate(ds::GeoDataset; dim::Symbol, func::Function) -> GeoDataset

Aggregate along a dimension.
"""
function aggregate(ds::GeoDataset; dim::Symbol, func::Function)
    new_vars = Dict{String, GeoArray}()
    for (name, ga) in ds.variables
        dim_idx = findfirst(d -> d.name == dim, ga.dims)
        isnothing(dim_idx) && continue
        new_data = dropdims(func(ga.data, dims=dim_idx); dims=dim_idx)
        new_dims = Tuple(d for (i, d) in enumerate(ga.dims) if i != dim_idx)
        new_vars[name] = GeoArray(new_data, new_dims, ga.crs, ga.attrs)
    end
    new_coords = Dict{Symbol, GeoArray}(dim => ga for (dim, ga) in ds.coords if dim != dim)
    new_dims = Dict{Symbol, Dimension}(dim => d for (dim, d) in ds.dims if dim != dim)
    return GeoDataset(new_vars, new_coords, new_dims, ds.crs, ds.attrs, ds.backend, ds.source)
end

# ============================================================
# Summary statistics
# ============================================================

"""
    data_stats(arr::AbstractArray) -> NamedTuple

Compute basic statistics for an array (ignoring NaN/missing).
"""
function data_stats(arr::AbstractArray)
    vals = collect(vec(Float64.(arr)))
    vals = vals[isfinite.(vals)]
    isempty(vals) && return (min=NaN, max=NaN, mean=NaN, std=NaN, count=0)
    return (min=minimum(vals), max=maximum(vals), mean=mean(vals), std=std(vals), count=length(vals))
end

"""
    variable_stats(ds::GeoDataset, varname::String) -> NamedTuple

Stats for a single variable.
"""
function variable_stats(ds::GeoDataset, varname::String)
    haskey(ds.variables, varname) || error("Unknown variable: $varname")
    return data_stats(ds.variables[varname].data)
end

"""
    geosummary(ds::GeoDataset) -> String

Human-readable summary of dataset.
"""
function geosummary(ds::GeoDataset)
    lines = String[]
    push!(lines, "GeoDataset: $(ds.source)")
    push!(lines, "Backend: $(ds.backend)")
    push!(lines, "Dimensions:")
    for (dim, d) in ds.dims
        push!(lines, "  $dim: $(d.size) points, units=$(d.units)")
    end
    push!(lines, "Variables:")
    for (name, ga) in ds.variables
        st = data_stats(ga.data)
        push!(lines, "  $name: $(size(ga.data)), min=$(st.min), max=$(st.max), mean=$(st.mean), count=$(st.count)")
    end
    return join(lines, "\n")
end

end # module GeoDataOperations