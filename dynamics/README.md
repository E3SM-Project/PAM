# Dynamical Core Overviews

## `dycoreA`
* State: rho, rho\*u, rho\*v, rho\*w, rho\*theta
* Conserves: rho, rho\*u, rho\*v, rho\*w, rho\*theta
* Acoustics: explicit, no sub-cycling
* Staggering: A-grid

Cell-centered Finite-Volume method, dimensionally split, upwind fluxes, ADER-DT and Runge-Kutta time discretizations, WENO limiting, positive-definite transport.

## `dycoreB`
* State: rho, u, v, w, rho\*theta
* Conserves: rho, u, v, w, rho\*theta
* Acoustics: HEVI, infinite and implicit options
* Staggering: A-grid

Cell-centered Finite-Volume method, dimensionally split, upwind fluxes, ADER-DT and Runge-Kutta time discretizations, WENO limiting, positive-definite transport.

