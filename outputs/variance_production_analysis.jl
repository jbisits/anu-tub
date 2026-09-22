# This script computes and saves the time and spatial mean, and the global integral of the
# variance prodiction diagnostics for all experiments in a `catalogue` defined in `plotting.jl`.
include("plotting.jl")

nm_file = "variance_production_025.jld2"
if isfile(nm_file)
    @info "output file for variance production and numerical mixing already exists: $(nm_file)"
else
    output_grid = "monthly" # this uses the saved data on the native layer grid
    nm = variance_production(catalogue, output_grid)
    save(nm_file, nm)
    @info "variance production and numerical mixing saved to $(nm_file)"
end
