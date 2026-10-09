"""
    assimilate_environmental_cube.jl

CLI driver for assimilating unified regional environmental and biogeochemical data cubes.

# Usage
    julia --project=. scripts/assimilate_environmental_cube.jl --config=configs/regions/scotian_shelf.toml
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using GeoData

function parse_args(args::Vector{String})
    cfg_file = nothing
    for a in args
        if startswith(a, "--config=")
            cfg_file = split(a, "=", limit = 2)[2]
        end
    end
    return cfg_file
end

function main()
    cfg_path = parse_args(ARGS)
    if isnothing(cfg_path)
        println("Error: --config=<path> argument is required.")
        println("Example:")
        println("  julia --project=. scripts/assimilate_environmental_cube.jl --config=configs/regions/scotian_shelf.toml")
        return 1
    end

    println("Loading configuration: $(cfg_path)")
    cfg = load_cube_config(cfg_path)
    output_path = assimilate_regional_cube(cfg; verbose = true)
    println("Successfully generated data cube at: $(output_path)")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
