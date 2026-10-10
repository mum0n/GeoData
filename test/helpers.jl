"""
Test helpers: the invariant checker and the scratch directory conventions.
"""

using GeoData
using GeoData: GeoDataset
using Dates, UUIDs

"""
    assert_dataset_invariants(ds::GeoDataset) -> nothing

Check the invariants that make a `GeoDataset` self-consistent, raising with the
offending name instead of leaving the mistake to surface as a wrong result:

1. Every variable's axis sizes match its dimension records.
2. Every dimension a variable spans is declared in `ds.dims`.
3. Every declared dimension is used by a variable or exists as a coordinate.
4. Every coordinate is 1-D and has a dimension record of matching size.

Correctness in this package is defined by these; a dataset that fails them is a bug, not
a style problem.
"""
function assert_dataset_invariants(ds::GeoDataset)
    for (name, ga) in ds.variables
        size(ga.data) == Tuple(d.size for d in ga.dims) || error(
            "Variable '$name' has shape $(size(ga.data)) but dimensions of size " *
            "$(Tuple(d.size for d in ga.dims)).")
        for d in ga.dims
            haskey(ds.dims, d.name) || error(
                "Variable '$name' spans :$(d.name), which ds.dims does not declare.")
        end
    end
    for (dim, d) in ds.dims
        used_by_var = any(any(dd.name == dim for dd in ga.dims) for (_, ga) in ds.variables)
        (used_by_var || haskey(ds.coords, dim)) || error(
            "Dimension :$dim is declared but no variable spans it and it has no " *
            "coordinate array.")
    end
    for (dim, ga) in ds.coords
        ndims(ga) == 1 || error("Coordinate :$dim must be 1-D, got $(ndims(ga)) dimensions.")
        haskey(ds.dims, dim) || error("Coordinate :$dim has no dimension record.")
        length(ga.data) == ds.dims[dim].size || error(
            "Coordinate :$dim has $(length(ga.data)) values but its dimension says " *
            "$(ds.dims[dim].size).")
    end
    return nothing
end

"""
    test_scratch() -> String

A fresh directory inside `test/scratch` - never `pwd()` and never the OS temp directory,
so a test run leaves nothing behind outside the project.
"""
function test_scratch()
    base = joinpath(@__DIR__, "..", "scratch", "test_$(UUIDs.uuid4())")
    mkpath(base)
    return base
end
