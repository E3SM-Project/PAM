#pragma once

#include <cmath>
#include <iostream>
//#include <cstring>
#include "mpi.h"
#include "yaml-cpp/yaml.h"
#include <array>
#include <complex>
#include <extensions/YAKL_fft.h>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>

using yakl::c::Bounds;
using yakl::c::parallel_for;
using yakl::c::SimpleBounds;

typedef unsigned long ulong;
typedef unsigned int uint;

////////////// These control the settings for SPAM++    //////////////

// Declaring the precision for the model
typedef double real;
typedef std::complex<real> complex;
typedef yakl::Array<complex, 5, yakl::memDevice, yakl::styleC> complex5d;
typedef yakl::Array<complex, 4, yakl::memDevice, yakl::styleC> complex4d;
typedef yakl::Array<complex, 3, yakl::memDevice, yakl::styleC> complex3d;
typedef yakl::Array<complex, 2, yakl::memDevice, yakl::styleC> complex2d;
typedef yakl::Array<complex, 1, yakl::memDevice, yakl::styleC> complex1d;
#define REAL_MPI MPI_DOUBLE
//#define REAL_NC NC_DOUBLE

// Spatial derivatives order of accuracy ie Hodge stars [2,4,6] (vert only
// supports 2 for now)
uint constexpr diff_ord = 2;
uint constexpr vert_diff_ord = 2;

// Hodge stars order for diffusion, for now only 2nd order is supported
uint constexpr diffusion_diff_ord = 2;
uint constexpr vert_diffusion_diff_ord = 2;

// Reconstruction types and order
enum class RECONSTRUCTION_TYPE { CFV, WENO, WENOFUNC };

RECONSTRUCTION_TYPE constexpr reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr dual_reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr dual_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr vert_reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr vert_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr dual_vert_reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr dual_vert_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr coriolis_reconstruction_type =
    RECONSTRUCTION_TYPE::CFV;
uint constexpr coriolis_reconstruction_order = 3;

RECONSTRUCTION_TYPE constexpr coriolis_vert_reconstruction_type =
    RECONSTRUCTION_TYPE::CFV;
uint constexpr coriolis_vert_reconstruction_order = 1;

uint constexpr max_reconstruction_order =
    std::max({reconstruction_order, dual_reconstruction_order,
              coriolis_reconstruction_order});

uint constexpr max_vert_reconstruction_order =
    std::max({vert_reconstruction_order, dual_vert_reconstruction_order,
              coriolis_vert_reconstruction_order});

enum class UPWIND_TYPE { HEAVISIDE, TANH };
constexpr real tanh_upwind_coeff = 250;
UPWIND_TYPE constexpr upwind_type = UPWIND_TYPE::TANH;
UPWIND_TYPE constexpr dual_upwind_type = UPWIND_TYPE::TANH;
UPWIND_TYPE constexpr vert_upwind_type = UPWIND_TYPE::TANH;
UPWIND_TYPE constexpr dual_vert_upwind_type = UPWIND_TYPE::TANH;

// How to handle PV flux term
// ADD AL81-TYPE SCHEME HERE EVENTUALLY AS WELL
enum class QF_MODE { EC, NOEC };
QF_MODE constexpr qf_choice = QF_MODE::EC;

// initial condition quadrature pts
uint constexpr ic_quad_pts_x = 5;
uint constexpr ic_quad_pts_y = 5;
uint constexpr ic_quad_pts_z = 5;

// Halo sizes
uint constexpr maxhalosize =
    std::max({(max_reconstruction_order - 1) / 2, diff_ord / 2});
uint constexpr mirroringhalo =
    std::max({(max_vert_reconstruction_order - 1) / 2, vert_diff_ord / 2});

int constexpr si_monitor_convergence = 2;
// 0 = do not monitor (does si_max_iters iterations)
// 1 = computes initial and final residual but still does si_max_iter iterations
// 2 = iterates until convergence or si_max_iter is reached

int constexpr si_verbosity_level = si_monitor_convergence;
// 0 = do not print
// 1 = print initial and final
// 2 = print every iteration

int constexpr si_max_iters = si_monitor_convergence > 1 ? 50 : 5;

#if defined _EXTRUDED && !defined _AN && !defined _MAN &&                      \
    (defined _IDEAL_GAS_POTTEMP || defined _CONST_KAPPA_VIRPOTTEMP)
bool constexpr si_compute_functional_derivatives_quadrature = false;
#else
bool constexpr si_compute_functional_derivatives_quadrature = true;
#endif
uint constexpr si_quad_pts = 4;

//////////////////////////////////////////////////////////////////////////////////////////////

real constexpr pi =
    3.141592653589793238462643383279502884197169399375105820974944_fp;

//#define PNETCDF_PUT_VAR ncmpi_put_vara_double
//#define PNETCDF_PUT_VAR_ALL ncmpi_put_vara_double_all

// GPU compilers sometimes have issues with zero-sized arrays than can occur
// for some parameter choices. For this reason we sometimes have to pad arrays
#if defined YAKL_ARCH_CUDA || defined YAKL_ARCH_HIP || defined YAKL_ARCH_SYCL
#define GPU_PAD 1
#else
#define GPU_PAD 0
#endif

// Specifying templated min and max functions
template <class T> YAKL_INLINE T mymin(T const v1, T const v2) {
  if (v1 < v2) {
    return v1;
  } else {
    return v2;
  }
}
template <class T> YAKL_INLINE T mymax(T const v1, T const v2) {
  if (v1 > v2) {
    return v1;
  } else {
    return v2;
  }
}

// Utility functions for serial output
void inline debug_print(std::string out, int masterproc) {
#ifdef PAM_DEBUG
  if (masterproc) {
    std::cout << out << "\n";
  }
#endif
}
void inline serial_print(std::string out, int masterproc) {
  if (masterproc) {
    std::cout << out << "\n";
  }
}

// Add mode for various operators
enum class ADD_MODE { REPLACE, ADD };

// Boundary types
enum class BND_TYPE { PERIODIC, NONE };

#if defined _HAMILTONIAN && defined _LAYER
#include "layermodel-common.h"
#elif defined _HAMILTONIAN && defined _EXTRUDED
#include "extrudedmodel-common.h"
#elif defined _ADVECTION && defined _LAYER
#include "layeradvection-common.h"
#elif defined _ADVECTION && defined _EXTRUDED
#include "extrudedadvection-common.h"
#endif

void finalize_parallel(Parameters &params, Parallel &par) {
  int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &par.actualrank);
  par.masterproc = par.actualrank == 0;
  params.masterproc = par.masterproc;
  if (params.inner_mpi) {
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &par.myrank);
  } else {
    par.nranks = 1;
    par.myrank = 0;
  }

  // if (!(par.nprocx * par.nprocy == par.nranks)) {endrun("Error: nranks !=
  // nprocx * nprocy");}

  // Get my process grid IDs
  par.py = floor(par.myrank / par.nprocx);
  par.px = par.myrank - par.nprocx * par.py;

  // Get my beginning and ending global indices; and domain sizes
  double nper;
  nper = ((double)params.nx_glob) / par.nprocx;
  par.i_beg = (int)round(nper * par.px);
  par.i_end = (int)round(nper * (par.px + 1)) - 1;
  par.nx = par.i_end - par.i_beg + 1;
  par.nx_glob = params.nx_glob;

  if (ndims >= 2) {
    nper = ((double)params.ny_glob) / par.nprocy;
    par.j_beg = (int)round(nper * par.py);
    par.j_end = (int)round(nper * (par.py + 1)) - 1;
    par.ny = par.j_end - par.j_beg + 1;
    par.ny_glob = params.ny_glob;
  }

  par.nz = params.nz_dual;
  par.nens = params.nens;

  // Determine my neighbors
  // x-dir
  par.x_neigh(0) = ij_to_l(wrap(par.px - 1, par.nprocx), par.py, par.nprocx);
  par.x_neigh(1) = ij_to_l(wrap(par.px + 1, par.nprocx), par.py, par.nprocx);

  // y-dir
  par.y_neigh(0) = -1;
  par.y_neigh(1) = -1;
  if (ndims == 2) {
    par.y_neigh(0) = ij_to_l(par.px, wrap(par.py - 1, par.nprocy), par.nprocx);
    par.y_neigh(1) = ij_to_l(par.px, wrap(par.py + 1, par.nprocy), par.nprocx);

    par.ll_neigh = ij_to_l(wrap(par.px - 1, par.nprocx),
                           wrap(par.py - 1, par.nprocy), par.nprocx);
    par.lr_neigh = ij_to_l(wrap(par.px + 1, par.nprocx),
                           wrap(par.py - 1, par.nprocy), par.nprocx);
    par.ur_neigh = ij_to_l(wrap(par.px + 1, par.nprocx),
                           wrap(par.py + 1, par.nprocy), par.nprocx);
    par.ul_neigh = ij_to_l(wrap(par.px - 1, par.nprocx),
                           wrap(par.py + 1, par.nprocy), par.nprocx);
  }

  // set boundaries
  par.xbnd = BND_TYPE::PERIODIC;
  par.ybnd = BND_TYPE::PERIODIC;

  // set halos
  par.halox = maxhalosize;
  par.haloy = maxhalosize;

// Debug output for the parallel decomposition
#ifdef PAM_DEBUG
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  for (int rr = 0; rr < par.nranks; rr++) {
    if (rr == par.myrank) {
      std::cout << "Hello! My rank is: " << par.myrank << "\n";
      std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
      std::cout << "I have: " << par.nx << " x " << par.ny << " grid points."
                << "\n";
      std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg
                << "\n";
      std::cout << "I end at index: " << par.i_end << " x " << par.j_end
                << "\n";
      std::cout << "My x neighbors are: " << par.x_neigh(0) << " "
                << par.x_neigh(1) << "\n";
      std::cout << "My y neighbors are: " << par.y_neigh(0) << " "
                << par.y_neigh(1) << "\n";
      std::cout << "My corner neighbors are: " << par.ll_neigh << " "
                << par.ul_neigh << " " << par.ur_neigh << " " << par.lr_neigh
                << "\n";
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }
  ierr = MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void check_and_print_parameters(const Parameters &params, const Parallel &par) {
  // Check time stepping params
  if (not(params.dtphys > 0.0_fp) or not(params.simSteps > 0) or
      not(params.crm_per_phys > 0) or not(params.Nout > 0) or
      not(params.Nstat > 0)) {
    endrun("spam++ must use step-based time control logic ie set simSteps >0, "
           "dtphys>0, crm_per_phys >0, outSteps >0, statSteps >0");
  }

  // Print out the values
  if (par.masterproc) {
    std::cout << "nx:         " << params.nx_glob << "\n";
    std::cout << "ny:         " << params.ny_glob << "\n";
    std::cout << "nl dual:         " << params.nz_dual << "\n";
    std::cout << "ni dual:         " << params.nz_dual + 1 << "\n";
    std::cout << "nens:         " << params.nens << "\n";

    std::cout << "halox:      " << par.halox << "\n";
    std::cout << "haloy:      " << par.haloy << "\n";

    std::cout << "dtcrm:         " << params.dtcrm << "\n";
    std::cout << "dtphys:         " << params.dtphys << "\n";
    std::cout << "Nsteps:     " << params.Nsteps << "\n";
    std::cout << "simSteps:     " << params.simSteps << "\n";
    std::cout << "crm per phys:     " << params.crm_per_phys << "\n";
    std::cout << "Nout:       " << params.Nout << "\n";
    std::cout << "outputName: " << params.outputName << "\n";

    std::cout << "nranks:     " << par.nranks << "\n";
    std::cout << "nprocx:     " << par.nprocx << "\n";
    std::cout << "nprocy:     " << par.nprocy << "\n";

    std::cout << "xlen:       " << params.xlen << "\n";
    std::cout << "ylen:       " << params.ylen << "\n";
    std::cout << "xc:         " << params.xc << "\n";
    std::cout << "yc:         " << params.yc << "\n";
  }
};
