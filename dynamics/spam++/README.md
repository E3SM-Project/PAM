# Synopsis
* SPAM++ is based on a Hamiltonian formulation of the equations of motion, using a discrete exterior calculus based on mimetic finite differences and (centered) finite-volume or WENO reconstructions for structured grids, with explicit RK time stepping
* Written in C++ using YAKL
* Prognostic variables are velocity (horizontal for layer model, horizontal and vertical for extruded model) and an arbitrary set of densities. The specification of a Hamiltonian determines the dynamics of the system.
* Currently implemented equation sets ie Hamiltonians are advection, (thermal) shallow water, compressible Euler and multicomponent compressible Euler
* For CE and MCE, the thermodynamics (ie equation of state) are arbitrary. Currently only ideal gas with potential temperature and constant kappa thermo with virtual potential temperature are implemented.
* There are two basic models: a layer model (2D) and an extruded model (1D+1D or 2D+1D). The former assumes doubly periodic topology with uniform rectangular geometry for an arbitrary number of layers. The latter (doubly) periodic horizontal topology with uniform rectangular geometry, along with an "extrusion" into the vertical with boundaries at the top and bottom. The geometry is currently uniform in the vertical, but this will soon be modified to support arbitrary vertical geometry.
* Targeted primarily at GPUs; with basic, unoptimized support for MPI and multicore CPUs (no vectorization, etc.)

