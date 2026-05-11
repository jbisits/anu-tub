# Utilities for plotting ANU-TUB output
"""
    struct GeoCoords{V, M}
Geographical coordinates for MOM6 output. Returned is `struct` with:
- centre point lat/lon
- corner point lat/lon
- seafloor depth
"""
struct GeoCoords{V, M}
    lonh :: V
    lath :: V
    lonq :: V
    latq :: V
    zl :: V
    zi :: V
    sf :: M
end

"""
    function get_geo_coords(static_file::AbstractString)
Get the geographical coordinates of the grid in `static_file`.
"""
function get_geo_coords(static_file::AbstractString, vc_file::AbstractString)

    ds = NCDataset(static_file, maskingvalue = NaN)

    lonh = ds["geolon"][:, 1]
    lath = ds["geolat"][1, :]
    lonq = ds["geolon_c"][:, 1]
    latq = ds["geolat_c"][1, :]
    depth = -ds["deptho"][:, :]

    close(ds)

    ds = NCDataset(vc_file, maskingvalue = NaN)

    zl = -ds["Layer"][:]
    zi = -ds["Interface"][:]

    close(ds)
    return GeoCoords(lonh, lath, lonq, latq, zl, zi, depth)
end

"""
    function barotropic_streamfunction(output_file; timestamps = Colon(), ρ₀ = 1035)
Compute the barotropic streamfunction from data saved in output_file.
The required output is:
- `umo_2d`, the vertically integrated mass transport. 
By default, the `mean` over all timesteps in the file will be returned.
Otherwise pass a range to take `mean` over or a single timestamp.
"""
function barotropic_streamfunction(output_file::AbstractString;
                                   timestamps = Colon(),
                                   ρ₀ = 1035)

    ds = NCDataset(output_file, maskingvalue = NaN)
    umo_2d = ds["umo_2d"][:, :, timestamps] # kgs⁻¹
    close(ds)

    umo_2d ./= (ρ₀ * 1e6)                      # Sv

    # Cumulatively integrate over meridional dimension
    ψ = nancumsum(umo_2d, dims = 2)

    return mean(ψ, dims = 3)
end
"""
    function overturning_circulation(output_file::AbstractString, layer::AbstractString; timestamps = Colon(), ρ₀ = 1035)
Compute the overturning from data in `output_file`. The required output is:
- `vmo` saved on the `rho2` (potential density referenced to 2000dbar) grid
By default, the `mean` over all timesteps in the file will be returned.
Otherwise pass a range to take `mean` over or a single timestamp.
"""
function overturning_circulation(output_file::AbstractString, layer::AbstractString;
                                 timestamps = Colon(),
                                 ρ₀ = 1035)

    ds = NCDataset(output_file, maskingvalue = NaN)
    vmo = ds["vmo"][:, :, :, timestamps] # kgs⁻¹
    l = ds[layer][:]
    close(ds)

    vmo ./= (ρ₀ * 1e6)                   # Sv

    # dims = 1 <=> sum over xh
    ψ = nansum(vmo, dims = 1)
    # Cumulatively sum over vertical dimension
    ψ = nancumsum(ψ, dims = 3)

    return l, mean(ψ, dims = 4)
end
"""
    function dθ_dz(output_file, timestamps = Colon())
Compute temperature gradient from data in `output_file`. The required output is:
- `thetao`, potential temperature
- `thkcello`, cell thickness
By default, the `mean` over all timesteps in the file will be returned.
Otherwise pass a range to take `mean` over or a single timestamp.
"""
function dθ_dz(output_file::AbstractString; timestamps = Colon())

    ds = NCDataset(output_file, maskingvalue = NaN)
    θ = ds["thetao"][:, :, :, timestamps]
    h = ds["thkcello"][:, :, :, timestamps]
    close(ds)

    Δθ = θ[:, :, 1:end-1, :] .- θ[:, :, 2:end, :]
    # Distance between cell centres from cell thicknesses
    Δz = 0.5 * (h[:, :, 1:end-1, :] .+ h[:, :, 2:end, :])

    Δθ ./= Δz # avoid more allocations

    return nanmean(Δθ, dims = (1, 4)) # return zonal and time mean
end
"""
    function zonal_variance_dissipation(output_file; timestamps = Colon())
Return the zonal variance dissipation calculated as:

        Zonal sum of numerical mixing diagnostic
       -----------------------------------------
              Zonal sum of layer thickness
"""
function zonal_variance_dissipation(output_file::AbstractString; 
                                    timestamps = Colon())

    # take zonal sum and drop dimension when reading in
    ds = NCDataset(output_file, maskingvalue = NaN)
    vd = nansum(ds["T_advection_scheme_variance_production"][:, :, :, timestamps], dims = 1)  # °C²ms⁻¹
    vd ./= nansum(ds["thkcello"][:, :, :, timestamps], dims = 1)                              # °C²s⁻¹
    close(ds)

    return mean(vd, dims = 4)
end

"""
    function depth_variance_dissipation(output_file::AbstractString, timestamps = Colon())
"""
function depth_variance_dissipation(output_file::AbstractString; timestamps = Colon())
    
    ds = NCDataset(output_file, maskingvalue = NaN)
    vd = nansum(ds["T_advection_scheme_variance_production"][:, :, :, timestamps], dims = 3)  # °C²ms⁻¹
    vd ./= nansum(ds["thkcello"][:, :, :, timestamps], dims = 3)                              # °C²s⁻¹
    close(ds)
    
    return mean(vd, dims = 4)
end
"""
    function zonal_numerical_mixing_diffusivity(output_file::AbstractString; timestamps = Colon())
Calculate the zonal numerical mixing as a diffusivity:

            zonal mean variance dissipation   # °C²s⁻¹
            -------------------------------
                    zonal mean |∇θ|²          # °C²m⁻²
"""
function zonal_numerical_mixing_diffusivity(output_file::AbstractString, static_file::AbstractString; 
                                            timestamps = Colon())

    ds = NCDataset(output_file, maskingvalue = NaN)
    vd = nansum(ds["T_advection_scheme_variance_production"][:, :, :, timestamps], dims = 1)  # °C²ms⁻¹
    vd ./= nansum(ds["thkcello"][:, :, :, timestamps], dims = 1)                              # °C²s⁻¹
    close(ds)
    
    interp_vd = nanmean(vd, dim = 1)
    interp_vd = 0.5 * (interp_vd[:, 1:end-1, :] .+ interp_vd[:, 2:end, :])
    interp_vd = 0.5 * (interp_vd[1:end-1, :, :] .+ interp_vd[2:end, :, :])
    
    ∇θ² = abs_temperature_gradient(output_file, static_file; timestamps)
    ∇θ² = nanmean(∇θ², dim = 1)

    κ_nm = interp_vd ./ ∇θ²
    replace!(κ_nm, -Inf => NaN)
    replace!(κ_nm, Inf => NaN)

    return mean(κ_nm, dims = 3)[:, :, 1]
end
"""
    function depth_numerical_mixing(output_file::AbstractString, static_file::AbstractString; 
                                   timestamps = Colon())
Caclulate the depth integrated numerical mixing as diffusivity:

            depth mean variance dissipation   # °C²s⁻¹
            -------------------------------
                    depth mean |∇θ|²          # °C²m⁻²

"""
function depth_numerical_mixing_diffusivity(output_file::AbstractString, static_file::AbstractString; 
                                            timestamps = Colon())
    
    ds = NCDataset(output_file, maskingvalue = NaN)
    vd = nansum(ds["T_advection_scheme_variance_production"][:, :, :, timestamps], dims = 3)  # °C²ms⁻¹
    vd ./= nansum(ds["thkcello"][:, :, :, timestamps], dims = 3)                              # °C²s⁻¹
    close(ds)
    interp_vd = 0.5 * (vd[1:end-1, :, :, :] .+ vd[2:end, :, :, :])
    interp_vd = 0.5 * (interp_vd[:, 1:end-1, :, :] .+ interp_vd[:, 2:end, :, :])
    interp_vd = nanmean(interp_vd, dim = 3)
    
    ∇θ² = abs_temperature_gradient(output_file, static_file; timestamps)
    ∇θ² = nanmean(∇θ², dim = 3)
    
    κ_nm = interp_vd ./ ∇θ²
    replace!(κ_nm, -Inf => NaN)
    replace!(κ_nm, Inf => NaN)

    return mean(κ_nm, dims = 3)[:, :, 1]
end
"""
    function abs_temperature_gradient(output_file::AbstractString, static_file::AbstractString; 
                                     timestamps = Colon())
Return the squared norm of the temperature gradient at each grid cell.
"""
function abs_temperature_gradient(output_file::AbstractString, static_file::AbstractString; 
                                 timestamps = Colon())
    
    ds = NCDataset(output_file, maskingvalue = NaN)
    θ = ds["thetao"][:, :, :, timestamps]
    h = ds["thkcello"][:, :, :, timestamps]
    Δθx = θ[1:end-1, :, :, :] .- θ[2:end, :, :, :]
    Δθy = θ[:, 1:end-1, :, :] .- θ[:, 2:end, :, :]
    Δθz = θ[:, :, 1:end-1, :] .- θ[:, :, 2:end, :]
    Δθz ./= 0.5 * (h[:, :, 1:end-1, :] .+ h[:, :, 2:end, :])     # dθ/dz
    close(ds)
    
    ds = NCDataset(static_file, maskingvalue = NaN)
    dx = 0.5 * (ds["dxt"][1:end-1, :] .+ ds["dxt"][2:end, :])
    dy = 0.5 * (ds["dyt"][:, 1:end-1] .+ ds["dxt"][:, 2:end])
    close(ds)
    
    for t in axes(Δθx, 4)
        for k in axes(Δθx, 3)
            Δθx[:, :, k, t] ./= dx  # dθ/dx           
        end
    end
    
    for t in axes(Δθy, 4)
        for k in axes(Δθy, 3)
            Δθy[:, :, k, t] ./= dy  # dθ/dy           
        end
    end
    
    # interpolation so all arrays are the same size
    Δθx = 0.5*(Δθx[:, 1:end-1, :, :] .+ Δθx[:, 2:end, :, :])
    Δθx = 0.5*(Δθx[:, :, 1:end-1, :] .+ Δθx[:, :, 2:end, :])
    
    Δθy = 0.5*(Δθy[1:end-1, :, :, :] .+ Δθy[2:end, :, :, :])
    Δθy = 0.5*(Δθy[:, :, 1:end-1, :] .+ Δθy[:, :, 2:end, :])
    
    Δθz = 0.5*(Δθz[1:end-1, :, :, :] .+ Δθz[2:end, :, :, :])
    Δθz = 0.5*(Δθz[:, 1:end-1, :, :] .+ Δθz[:, 2:end, :, :])
    
    return Δθx.^2 .+ Δθy.^2 .+ Δθz.^2
end
"""
    function vertical_sum(output_file::AbstractString; timestamps = Colon())
Sum the thickness weighted variance dissipation over the vertical dimension.
Because of the thickness weighting this is a vertical integral.
"""
function vertical_sum(output_file::AbstractString; timestamps = Colon())

    ds = NCDataset(output_file, maskingvalue = NaN)
    ∫nm = nansum(ds["T_advection_scheme_variance_production"][:, :, :, timestamps], dims = 3)
    close(ds)

    return mean(∫nm, dims = 4)
end
"""
    function global_integral(output_file::AbstractString; timestamps = Colon())
Global integral of numerical mixing at each saved timestep.
"""
function global_integral(output_file::AbstractString; timestamps = Colon())
    
    ds = NCDataset(output_file, maskingvalue = NaN)
    ∫nm = nansum(ds["T_advection_scheme_variance_production"][:, :, :, timestamps], dim = (1, 2, 3))
    close(ds)

    return ∫nm
end
"""
    function interface_depth(output_file::AbstractString; lonslice = 80, timestamp = 1)
Find the height of the model interfaces by cumulatively summing the thickness field.
By default, the slice from the middle of the longitude domain, and the initial timestamp,
are used. **NOTE:** the returned interfaces are as depths so they are positive.
"""
function inteface_depth(output_file::AbstractString; lonslice = 80, timestamp = 1)

    ds = NCDataset(output_file, maskingvalue = NaN)
    h = ds["thkcello"][80, :, :, timestamp]
    int_depth = cumsum(h, dims = 2)
    close(ds)
    
    return int_depth
end
"""
    function interface_height(output_file::AbstractString; lonslice = 80, timestamp = 1
Return the interfeace_height which is the negative depth. This is just a convenience function.
"""
interface_height(output_file::AbstractString; lonslice = 80, timestamp = 1) = 
    -inteface_depth(output_file; lonslice, timestamp)