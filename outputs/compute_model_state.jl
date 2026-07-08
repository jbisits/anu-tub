# This script computes and saves the time and spatial mean of the
# model state diagnostics for all experiments in a `catalogue`
# which is defined in `plotting.jl`.
include("plotting.jl")

vg_file = "vertical_grids.jld2"
if isfile(vg_file)
    @info "output file for vertical grids already exists: $(vg_file)"
else
    vg = vertical_grid(catalogue, time_idx = 24)
    save(vg_file, vg)
    @info "vertical grids saved to $(vg_file)"
end

ms_file = "model_states.jld2"
if isfile(ms_file)
    @info "output file for model states already exists: $(ms_file)"
else
    ms = model_state(catalogue)
    save(ms_file, ms)
    @info "model states saved to $(ms_file)"
end
