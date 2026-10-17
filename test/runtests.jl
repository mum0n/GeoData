# ==============================================================================
# GeoData Comprehensive Test Suite Orchestrator
# ==============================================================================

import Pkg

const RUNNING_UNDER_PKG_TEST = basename(dirname(Base.active_project())) == "test"
if !RUNNING_UNDER_PKG_TEST
    Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
    Pkg.instantiate(; io = devnull)
end

# Usage:
#   julia --project=. test/runtests.jl                # Run every segment
#   julia --project=. test/runtests.jl core           # Run only that segment
#   julia --project=. test/runtests.jl operations     # Run operations
#   julia --project=. test/runtests.jl storage        # Run storage & backends
#   julia --project=. test/runtests.jl spatial        # Run spatial utilities
#
# Under `Pkg.test()` the whole suite runs; pass segments as test_args:
#   Pkg.test(test_args = ["core", "storage"])
# ==============================================================================

using Test
include(joinpath(@__DIR__, "test_fixtures.jl"))

const AVAILABLE_SEGMENTS = Dict{Symbol, @NamedTuple{file::String, desc::String}}(
    :core => (
        file = "test_core.jl",
        desc = "Core Types & Dataset Invariants"
    ),
    :operations => (
        file = "test_operations.jl",
        desc = "Slicing, Selection, Joining & Aggregation"
    ),
    :storage => (
        file = "test_storage.jl",
        desc = "Storage Engines, Backends (Zarr, NC, Parquet) & Time Series"
    ),
    :spatial => (
        file = "test_spatial.jl",
        desc = "Spatial Utilities, Bathymetry & Marine Cell Masking"
    )
)

"""
    parse_requested_segments(args::Vector{String}) -> Vector{Symbol}

Resolves CLI segment arguments to segment keys.
"""
function parse_requested_segments(args::Vector{String})::Vector{Symbol}
    if isempty(args)
        return collect(keys(AVAILABLE_SEGMENTS))
    end

    selected = Symbol[]
    for a in args
        key = Symbol(lowercase(strip(a)))
        if haskey(AVAILABLE_SEGMENTS, key)
            push!(selected, key)
        else
            avail = sort(string.(collect(keys(AVAILABLE_SEGMENTS))))
            @warn "Unrecognized test segment: '$a'. Available segments: $(join(avail, ", "))"
        end
    end

    return isempty(selected) ? sort!(collect(keys(AVAILABLE_SEGMENTS))) : unique(selected)
end

requested_segments = parse_requested_segments(ARGS)

@testset "GeoData Comprehensive Suite" begin
    for seg_key in requested_segments
        info = AVAILABLE_SEGMENTS[seg_key]
        println("\n" * "="^80)
        println(">>> Running Segment: $(info.desc) ($(info.file))")
        println("="^80)
        include(joinpath(@__DIR__, info.file))
    end
end
