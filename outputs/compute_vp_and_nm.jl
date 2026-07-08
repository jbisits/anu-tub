# This script computes and saves the time and spatial mean of the
# variance prodiction and numerical mixing for all experiments in a `catalogue`
# which is defined in `plotting.jl`.
include("plotting.jl")

nm_file = "vp_and_nm.jld2"
if isfile(nm_file)
    @info "output file for variance production and numerical mixing already exists: $(nm_file)"
else
    output_grid = "monthly" # this uses the saved data on the native layer grid
    nm = variance_production_and_numerical_diffusivity(catalogue, output_grid)
    save(nm_file, nm)
    @info "variance production and numerical mixing saved to $(nm_file)"
end
