# anu-tub
A bathtub-sector ocean configuration for MOM6.
This configuration is designed primarily as a test case for implementing new vertical coordinate systems into MOM6, and to compare with existing vertical coordinate.
Thus, the design criteria is to produce the simplest configuration that incorporates buoyancy and wind forcing, gyre circulation, overturning circulation and topographically-constrained overflows.
The configuration includes the following features:
* The sector is 40° wide, and goes from ~70.3°S to ~70.3°N.
* It has a 1° nominal resolution and a Mercator grid refinement.
* the default vertical coordinate is ZSTAR with 75 vertical levels.
* The bathymetry is simple, with a vertical wall to the north, an "Antarctic shelf and slope" in the south and sloping sidewalls on the east and west.
* The domain is periodic in the east-west direction, allowing zonal flow in a narrow "Drake Passage" between ~65°S and ~52°S.
* Surface momentum forcing is via a prescribed zonal wind field that is constant, but varies with latitude.
* Thermal forcing is through relation to a latitude-dependent SST profile, and there is (currently) no freshwater forcing so that salinity is constant. 
* We use the WRIGHT equation of state, without frazil formation or sea ice.
* The model includes the PBL surface boundary layer scheme with a contant background diffusivity of $2 \times 10^{-5}$.
* The MEKE eddy parameterisation scheme is used with the default configuration being what is currently set in the [ACCESS-OM3 100km](https://github.com/ACCESS-NRI/access-om3-configs/tree/dev-MC_100km_jra_ryf) configuration.

The model has $40 \times 200$ grid points, with a tile layout of $4 \times 13$ to run efficiently on 52 cores.
With a 1200-second timestep, the standard ZSTAR case runs at ~40 simulated years/day and consumes ~90 SU per model year.
