# Synopsis
* SPAM++ is based on Hamiltonian, using mimetic finite differences and finite-volume or WENO reconstructions, RK or ADER-DT time stepping, etc.
* Written in C++ using YAKL
* Advection, (thermal) shallow water, compressible Euler, multicomponent compressible Euler
* Arbitrary equations of state
* Prognostic variables are velocity, component densities and the density of an arbitrary function of specific entropy and wet mixing ratios (thermodynamic scalar density)
* Targeted primarily at GPUs; with basic, unoptimized support for MPI and multicore CPUs (no vectorization, etc.)

