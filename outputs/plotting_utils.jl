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
    nm = nansum(ds["T_numerical_mixing"][:, :, :, timestamps], dims = 1) # °c²ms⁻¹
    nm ./= nansum(ds["thkcello"][:, :, :, timestamps], dims = 1) # °C²s⁻¹
    close(ds)

    return mean(nm, dims = 4)
end
"""
    function zonal_numerical_mixing_diffusivity(output_file::AbstractString; timestamps = Colon())
Calculate the zonal numerical mixing as a diffusivity:

            variance dissipation        # °C²s⁻¹
            --------------------
                 (dθ/dz)²               # °C²m⁻²
"""
function zonal_numerical_mixing_diffusivity(output_file::AbstractString; 
                                            timestamps = Colon())

    variance_dissipation = zonal_variance_dissipation(output_file; timestamps)
    interp_vd = 0.5 * (variance_dissipation[:, :, 1:end-1, :] .+ variance_dissipation[:, :, 2:end, :])
    dθdz = dθ_dz(output_file; timestamps)
    dθdz .^= 2

    κ_nm = interp_vd ./ dθdz
    replace!(κ_nm, -Inf => NaN)
    replace!(κ_nm, Inf => NaN)

    return κ_nm
end
"""
    function vertical_sum(output_file::AbstractString; timestamps = Colon())
Sum the thickness weighted variance dissipation over the vertical dimension.
Because of the thickness weighting this is a vertical integral.
"""
function vertical_sum(output_file::AbstractString; timestamps = Colon())

    ds = NCDataset(output_file, maskingvalue = NaN)
    ∫nm = nansum(ds["T_numerical_mixing"], dims = 3)
    close(ds)

    return mean(∫nm, dims = 4)
end
