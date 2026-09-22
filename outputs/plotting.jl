# Single script to setup the plotting and analysis environemt
using Revise

# This is needed for the jupyter notebooks
using Pkg
Pkg.activate(".")

# Packages required
using NCDatasets, CairoMakie, NaNStatistics, Interpolations, Statistics
using GibbsSeaWater, Dates, Printf, Glob, JLD2

include("plotting_config.jl")
Revise.includet("plotting_utils.jl")
Revise.includet("generate_tub_ics.jl")

@info "Plotting and analysis environment setup!"

# the directories with the output need to be set manully because of saving
expt_output = Dict(
                   "zstar"              => ["output01"*i for i ∈ string.(6:9)],
                   "hycom1"             => ["output01"*i for i ∈ string.(6:9)],
                   "alt-hycom"          => ["output01"*i for i ∈ string.(4:7)],
                   "adapt"              => ["output01"*i for i ∈ string.(4:7)],
                   "zstar-weaker-KD"    => ["output01"*i for i ∈ string.(2:5)],
                   "hycom1-weaker-KD"   => ["output01"*i for i ∈ string.(2:5)],
                   "alt-hycom-weaker-KD"=> ["output01"*i for i ∈ string.(2:5)],
                   "adapt-weaker-KD"    => ["output01"*i for i ∈ string.(2:5)]
)

catalogue = Dict{String, Any}()

for (i, expt) ∈ enumerate(keys(expt_output))

   expt_path = joinpath(pwd(), "anu-tub-nm-" * expt)
   odir = expt_path .* "/" .* expt_output[expt]
   monthly = [glob("ocean_month.nc", d)[1] for d ∈ odir]
   monthlyz = [glob("ocean_month_z.nc", d)[1] for d ∈ odir]
   monthlyrho2 = [glob("ocean_month_rho2.nc", d)[1] for d ∈ odir]
   static = joinpath(odir[end], "ocean_static.nc")
   vc = joinpath(odir[end], "Vertical_coordinate.nc")

   catalogue[expt] = Dict("path" => expt_path,
                          "monthly" => monthly,
                          "monthlyrho2" => monthlyrho2,
                          "monthlyz" => monthlyz,
                          "static" => static, 
                          "vc" => vc)
end

@info "Experiments in the catalogue are: $(keys(expt_output))
      Available data is:
      - monthly
      - monthlyz
      - monthlyrho2
      - static
      - vertical coordinate"

vgrids = joinpath(pwd(), "vertical_grids.jld2")
model_states = joinpath(pwd(), "model_states.jld2")
vp025 = joinpath(pwd(), "variance_production_025.jld2")

@info "Analysis files available:
      - vgrids: $(vgrids)
      - model_states: $(model_states)
      - variance_production_025: $(vp025)"
