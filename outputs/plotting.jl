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

expts = ("zstar", "AG", "hycom", "alt-hycom")
# needs to be set manually because of the way the experiments have been run
output = (
   # vcat(["output00"*i for i ∈ string.(2:9)], ["output01"*i for i ∈ string.(0:5)]),
   # vcat(["output00"*i for i ∈ string.(2:9)], ["output01"*i for i ∈ string.(0:5)]),
   # vcat(["output00"*i for i ∈ string.(2:9)], ["output01"*i for i ∈ string.(0:3)]),
   # vcat(["output00"*i for i ∈ string.(2:9)], ["output01"*i for i ∈ string.(0:3)])
   # susbset for testing
   vcat(["output00"*i for i ∈ string.(2:3)]),
   vcat(["output00"*i for i ∈ string.(2:3)]),
   vcat(["output00"*i for i ∈ string.(2:3)]),
   vcat(["output00"*i for i ∈ string.(2:3)]),
)

catalogue = Dict{String, Any}()

for (i, expt) ∈ enumerate(expts)

   expt_path = joinpath(pwd(), "anu-tub-nm-" * expt)
   odir = expt_path .* "/" .* output[i]
   daily = [glob("ocean_daily.nc", d)[1] for d ∈ odir]
   monthly = [glob("ocean_month.nc", d)[1] for d ∈ odir]
   monthlyz = [glob("ocean_month_z.nc", d)[1] for d ∈ odir]
   monthlyrho2 = [glob("ocean_month_rho2.nc", d)[1] for d ∈ odir]
   static = joinpath(odir[end], "ocean_static.nc")
   vc = joinpath(odir[end], "Vertical_coordinate.nc")

   catalogue[expt] = Dict("path" => expt_path,
                          "daily" => daily,
                          "monthly" => monthly,
                          "monthlyrho2" => monthlyrho2,
                          "monthlyz" => monthlyz,
                          "static" => static, 
                          "vc" => vc)
end

@info "Experiments in the catalogue are:
       - $(expts[1])
       - $(expts[2])
       - $(expts[3])
       - $(expts[4])
       Available data is:
          - daily
          - monthly
          - monthlyz
          - monthlyrho2
          - static
          - vertical coordinate"

vgrids = joinpath(pwd(), "vertical_grids.jld2")
model_states = joinpath(pwd(), "model_states.jld2")
vp_and_nm = joinpath(pwd(), "vp_and_nm.jld2")

@info "Analysis files available:
      - vgrids: $(vgrids)
      - model_states: $(model_states)
      - vp_and_nm: $(vp_and_nm)"