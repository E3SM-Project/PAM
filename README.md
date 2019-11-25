# SPAM: Scalable and Portable Atmospheric Model

PAM is a portable atmospheric model written in C++ with performance portability in mind using a kernel launching approach. It works on CPUs, Nvidia GPUs, and AMD GPUs. The goal for PAM is to allow easily exchangeable dynamical core and physics options through a simple and clear interface layer that's common to each. An emphasis is placed on algorithms that give good hardware uitilization on accelerated architectures and MPI patterns that give good scalability.

## Directory Structure

* `build`: Machine files for compiling
* `common`: Source & header files with functions and variables used in multiple places
* `driver`: Code to drive the standalone model
* `dynamics`: Dynamical core code, each version in its own sub-directory. Each dycore must provide:
  * `math_desc.pdf`: documenting the numerical discretizations used
  * `README.md`: giving a brief synopsis of the key numerical features
* `externals`: Submodules of external code (`cub`, `hipcub`, `rocPRIM`, `kokkos`, `YAKL`)
* `physics`: Physics parameterizations (microphysics, radiation, etc.)
* `sage`: SageMath code for generating C++ code
* `utils`: Miscillaneous utilities

## Software Dependencies
* MPI library
* parallel-netcdf (https://trac.mcs.anl.gov/projects/parallel-netcdf)
  * See `utils/build_parallel_netcdf.sh` for a sample build script.
