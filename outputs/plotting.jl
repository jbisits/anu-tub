# Single script to setup the plotting and analysis environemt
using Revise

# This is needed for the jupyter notebooks
using Pkg
Pkg.activate(".")

# Packages required
using NCDatasets, CairoMakie, NaNStatistics, Interpolations, Statistics
using GibbsSeaWater, Dates, Printf, Glob

include("plotting_config.jl")
Revise.includet("plotting_utils.jl")
Revise.includet("generate_tub_ics.jl")

@info "Plotting and analysis environment setup!"

expts = ("zstar", "hycom", "alt-hycom", "AG")
output = (["output007"], ["output007"], ["output005"], ["output005"])

catalogue = Dict{String, Any}()

for (i, expt) ∈ enumerate(expts)

   expt_path = joinpath(pwd(), "anu-tub-nm-" * expt)
   # odir = joinpath(expt_path, output[i])
   # daily = joinpath(odir, "ocean_daily.nc")
   # monthly = joinpath(odir, "ocean_month.nc")
   # monthlyz = joinpath(odir, "ocean_month_z.nc")
   # monthlyrho2 = joinpath(odir, "ocean_month_rho2.nc")
   # static = joinpath(odir, "ocean_static.nc")
   # vc = joinpath(odir, "Vertical_coordinate.nc")
   odir = expt_path .* "/" .* output[i]
   println(odir)
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
