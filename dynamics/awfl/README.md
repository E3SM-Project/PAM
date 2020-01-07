# Synopsis
* State: rho, u, v, w, theta
* Conserves: non-conservative form, so basically it conserves nothing :)
* Acoustics: explicit, no sub-cycling
* Order of accuracy: any odd order greater than one in time or space (time can be different than space)
* Staggering: A-grid
* Transports: Total mass, not mixing ratios
* Other: None

# Method Description
Cell-centered Finite-Volume method, dimensionally split, high-order upwind flux difference splitting, ADER-DT time stepping with temporally linearized acoustic dynamics within a time step, WENO limiting, and positive-definite transport. Hydrostasis is removed from the vertical momentum equation by an additive term.

# Definitions
* rho: total density
* u,v,w: air velocity
* theta: potential temperature
