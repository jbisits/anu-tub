"""
    function generate_tub_ics!(ic_filename, monthly_output, static_output)
Generate a netcdf file that containes salinity, temperature, thickness and interface height
from mean sea level for use with an ANU-TUB run. **NOTE:** `ic_filename` must end in `.nc`.
The snapshot used will be the last saved monthly output.
"""
function generate_tub_ics!(ic_filename::AbstractString, monthly_output::AbstractString, static_output::AbstractString)
   
    # Get the ocean depth
    ds = NCDataset(static_output, maskingvalue = Float32(NaN))
    depth = ds["deptho"][:, :]
    close(ds)
    depth_nan = .!isnan.(depth) .== true

    ds = NCDataset(monthly_output, maskingvalue = Float32(NaN))

    # Read thickness and find the `NaN` values
    h = ds["thkcello"][:, :, :, end]
    find_nan = .!isnan.(h) .== true

    # Generate e from thickness and ocean depth
    nx, ny, nz = size(h)
    e = Array{eltype(h)}(undef, nx, ny, nz + 1)
    e[:, :, nz + 1] .= -depth
    for k ∈ nz:-1:1
        e[:, :, k] .= e[:, :, k + 1] + h[:, :, k]
    end

    # Correctly set the missing value for saving
    h .*= find_nan
    replace!(h, 0 => Float32(1.0e20))

    e[:, :, 1:nz] .*= find_nan
    e[:, :, nz] .*= depth_nan
    replace!(e, 0 => Float32(1.0e20))

    # Get temperature
    T = ds["thetao"][:, :, :, end]
    T .*= find_nan
    replace!(T, 0 => Float32(1.0e20))

    # Create a salinity array that is constant value of 35
    S = similar(ds["thetao"][:, :, :, end])
    S .= 35
    S .*= find_nan
    replace!(S, 0 => Float32(1.0e20))

    # Save the new initial conditions file   
    _ds = NCDataset(ic_filename, "c")
    defVar(_ds, "h", h, ("lonh", "lath", "Layer"), fillvalue = Float32(1.0e20), attrib = ds["thkcello"].attrib)
    defVar(_ds, "eta", e, ("lonh", "lath", "Interface"), fillvalue = Float32(1.0e20), 
            attrib = Dict("long_name" => "Interface Height Relative to Mean Sea Level",
                          "units" => "m", "missing_value" => 1.0e20)
    )
    defVar(_ds, "Temp", T, ("lonh", "lath", "Layer"), fillvalue = Float32(1.0e20), attrib = ds["thetao"].attrib)
    defVar(_ds, "Salt", S, ("lonh", "lath", "Layer"), fillvalue = Float32(1.0e20), attrib = Dict("long_name" => "Salinity",
                                                                    "units" => "PPT",
                                                                    "missing_value" => 1.0e20)
    )
    defVar(_ds, "time", [ds["time"][end]], ("time",), attrib = ds["time"].attrib)

    close(_ds)
    @info "Initial conditions saved to $(ic_filename)"
    close(ds)

    return nothing
end
"""
    function generate_tubTS_ics!(ic_filename::AbstractString, monthly_output::AbstractString)
Generate salinity and temperature conditions from the last temperature snapshot in`monthly_output`.
Salinity is also set but this is constant (35psu).
These initial conditions are then written to `ic_filename` which **must** have a netcdf file extension.
"""
function generate_tubTS_ics!(ic_filename::AbstractString, monthly_output::AbstractString)

    ds = NCDataset(monthly_output, maskingvalue = Float32(NaN))
    # Get temperature
    T = ds["thetao"][:, :, :, end]
    T .*= find_nan
    replace!(T, 0 => Float32(1.0e20))

    # Create a salinity array that is constant value of 35
    S = similar(ds["thetao"][:, :, :, end])
    S .= 35
    S .*= find_nan
    replace!(S, 0 => Float32(1.0e20))

    # Save the new initial conditions file   
    _ds = NCDataset(ic_filename, "c")
    defVar(_ds, "Temp", T, ("lonh", "lath", "Layer"), fillvalue = Float32(1.0e20), attrib = ds["thetao"].attrib)
    defVar(_ds, "Salt", S, ("lonh", "lath", "Layer"), fillvalue = Float32(1.0e20), attrib = Dict("long_name" => "Salinity",
                                                                    "units" => "PPT",
                                                                    "missing_value" => 1.0e20)
    )
    defVar(_ds, "time", [ds["time"][end]], ("time",), attrib = ds["time"].attrib)

    close(_ds)
    @info "Initial conditions saved to $(ic_filename)"
    close(ds)

    return nothing
end
