# The AWFL Cloud Dynamical Core

## Note that this is a work in progress

A stratified, non-hydrostatic, fully compressible 3-D Cartesian grid cloud model dynamical core using [A]DER [W]ENO [F]inite-Vo[L]ume (AWFL) Numerics. The numerical scheme is designed to perform well on accelerators by performing significant amounts of local work in between DRAM memory accesses. It is written in C++ using the Kokkos library to accelerate the computations on a variety of hardware backends.

Some attributes of the numerical scheme:
* __ADER__: ADER time stepping using Differential Transforms for single-step, single-stage high-order time stepping
* __WENO__: A new function-based Weighted Essentially Non-Oscillatory (WENO) limiting for crisp resolution of discontinuities without oscillations. 
* __Finite-Volume__: Flux-vector Finite-Volume spatial integration for local conservation
* __Upwind__ fluxes using the upwind Godunov linear state via flux-vector splitting for slight dissipation in the numerical scheme
* __High-Order__: Up to ninth-order accuracy with WENO limiting and 15th-order accuracy without WENO limiting
* __Compressible__: Solves the fully compressible Euler equations for applicability to the full range of atmospheric motions
* __Non-Hydrostatic__: Resolves non-hydrostatic motions, thus admitting vertically propagating acoustic waves
* __Dimensionally-Split__ for reduced cost
* Uses __Perturbed__ quantities with the hydrostatic background state removed for better resolution of hydrostatic balance
* Includes a __2-D__ option to simulate the traditional NH atmospheric test cases
* Uses __MPI__ (nearest neighbor communication) for off-node parallelism
* Uses __Kokkos__ for on-node parallelism

## Future Work
The most pressing developments needed are:
1) A more efficient handling of fast acoustics, particularly in the vertical direction. This could include implicit, multi-level, and Lagrangian treatments.
2) A variable vertical grid spacing (simply a matter of altering the SageMath scripts to compute Vandermonde matrices that respond to grid spacing and inverting those at model initialization).
3) Add tracer transport routines
4) Create a plan for incorporating moisture, possibly through an equivalent potential temperature.

## Software Dependencies
* Nvidia GPUs require CUB: https://github.com/NVlabs/cub
* AMD GPUs require hipCub and rocPrim:
  * https://github.com/ROCmSoftwarePlatform/hipcub
  * https://github.com/ROCmSoftwarePlatform/rocPRIM

## How to Build
```
ln -s build/[whatever_machine].inc mach.inc
make
```

* Note that you'll need to have the kokkos repo (https://github.com/kokkos/kokkos) cloned somewhere and pointed to by your mach.inc file. See the existing files for examples.
* Requires an MPI library (that provides mpic++) as well as parallel-netcdf (https://trac.mcs.anl.gov/projects/parallel-netcdf)

