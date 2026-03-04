# Single script to setup the plotting and analysis environemt
using Revise

# This is needed for the jupyter notebooks
using Pkg
Pkg.activate(".")

# Packages required
using NCDatasets, CairoMakie, NaNStatistics, Interpolations, Statistics
using GibbsSeaWater, Dates, Printf

include("plotting_config.jl")
Revise.includet("plotting_utils.jl")
Revise.includet("generate_tub_ics.jl")

@info "Plotting and analysis environment setup!"

expts = ("zstar-PPMH3", "hycom-PPMH3")

catalogue = Dict{String, Any}()

for expt ∈ expts

    expt_path = joinpath(pwd(), "anu-tub-nm-" * expt)
    odir = expt == "zstar-PPMH3" ? joinpath(expt_path, "output442") : joinpath(expt_path, "output000")
    daily = joinpath(odir, "ocean_daily.nc")
    monthly = joinpath(odir, "ocean_month.nc")
    monthlyz = joinpath(odir, "ocean_month_z.nc")
    monthlyrho2 = joinpath(odir, "ocean_month_rho2.nc")
    static = joinpath(odir, "ocean_static.nc")
    vc = joinpath(odir, "Vertical_coordinate.nc")

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
       Available data is:
          - daily
          - monthly
          - monthlyz
          - monthlyrho2
          - static
          - vertical coordinate"
