# Synopsis
* State: rho, rho\*u, rho\*v, rho\*w, rho\*theta
* Conserves: rho, rho\*u, rho\*v, rho\*w, rho\*theta
* Acoustics: explicit, no sub-cycling
* Staggering: A-grid

# Method Description
Cell-centered Finite-Volume method, dimensionally split, upwind fluxes, ADER-DT and Runge-Kutta time discretizations, WENO limiting, positive-definite transport.

# Definitions
rho:total density
rho_d: dry density
u,v,w: air velocity
theta: potential temperature
