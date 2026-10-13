# Notebooks

Notebooks for the generic geospatial layer: the containers, the three backends, and the
coordinate operations that work on all of them.

| Notebook | What it covers |
|---|---|
| [`geodata_tour.jl`](geodata_tour.jl) | `GeoDataset` and its invariants; reading and writing through Zarr, NetCDF and GeoParquet; slicing and the one contract behind it; point queries; aggregation and joins; regridding; the land mask; grouped storage. |

Pluto notebooks: `.jl` files in Pluto's own format, with `md"""..."""` cells between the code.
Each step is an ordinary package call, so a cell can be re-run while you inspect what it
returned — and the markdown between the cells documents what every option means, which
assumptions the layer bakes in, and the failure modes it chooses instead of guessing.

## Running them

Pluto is deliberately **not** a dependency of this package: a notebook is how you look at the
work, not something the run needs. Launched from the package directory, `Pkg.activate(@__DIR__)`
inside the first cells pins the notebook to this package's environment:

```
julia --project=. -e 'using Pluto; Pluto.run()'
```

then open the file from the notebook list. The first cells also `cd` to the package root, so
`scratch/` paths stay inside the project.

The notebook builds its own data and writes it under `scratch/`, which is git-ignored, so it runs
with no downloads and leaves nothing outside the project. It does not touch the shared lakehouse.

## What the notebook reads

Nothing. Every input is synthesised in a cell: a 3-D temperature field, a 2-D seabed, and a small
coastline written as GeoParquet. That is deliberate — the layer's contract is that a dataset says
what its axes mean, so the clearest test of it is to watch a call refuse to guess. Where the
notebook shows an error (a Zarr store written twice, a joint dataset sharing a variable name), it
is the documented behaviour, not a defect.
