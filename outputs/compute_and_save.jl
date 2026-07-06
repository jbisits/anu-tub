# This script computes and saves the vertical coordinates, model state,
# variance prodiction and numerical mixing for all experiments in a `catalogue`
# which is defined in `plotting.jl`.
include("plotting.jl")

vg_file = "vertical_grids.jld2"
if isfile(vg_file)
    @info "output file for vertical grids already exists: $(vg_file)"
else
    vg = vertical_grid(catalogue)
    jldsave(vg_file; data=vg)
    @info "vertical grids saved to $(vg_file)"
end

ms_file = "model_states.jld2"
if isfile(ms_file)
    @info "output file for model states already exists: $(ms_file)"
else
    ms = model_state(catalogue)
    jldsave(ms_file; data=ms)
    @info "model states saved to $(ms_file)"
end

nm_file = "vp_and_nm.jld2"
if isfile(nm_file)
    @info "output file for variance production and numerical mixing already exists: $(nm_file)"
else
    nm = variance_production_and_numerical_diffusivity(catalogue, output_grid)
    jldsave(nm_file; data=nm)
    @info "variance production and numerical mixing saved to $(nm_file)"
end
