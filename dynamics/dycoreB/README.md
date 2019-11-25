# Synopsis
* State: rho, u, v, w, rho\*theta
* Conserves: rho, u, v, w, rho\*theta
* Acoustics: HEVI, infinite and implicit options
* Staggering: A-grid

# Method Description
Cell-centered Finite-Volume method, dimensionally split, upwind fluxes, ADER-DT and Runge-Kutta time discretizations, WENO limiting, positive-definite transport.

# Definitions
rho:total density
rho_d: dry density
u,v,w: air velocity
theta: potential temperature
