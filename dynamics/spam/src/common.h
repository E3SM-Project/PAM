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
#include <optional>
#include <sstream>
#include <string>

using yakl::c::Bounds;
using yakl::c::parallel_for;
using yakl::c::SimpleBounds;

using ulong = unsigned long;
using uint = unsigned int;

////////////// These control the settings for SPAM++    //////////////

// Declaring the precision for the model
using real = double;

using complex = std::complex<real>;
using complex5d = yakl::Array<complex, 5, yakl::memDevice, yakl::styleC>;
using complex4d = yakl::Array<complex, 4, yakl::memDevice, yakl::styleC>;
using complex3d = yakl::Array<complex, 3, yakl::memDevice, yakl::styleC>;
using complex2d = yakl::Array<complex, 2, yakl::memDevice, yakl::styleC>;
using complex1d = yakl::Array<complex, 1, yakl::memDevice, yakl::styleC>;

using optional_real1d = std::optional<real1d>;
using optional_real2d = std::optional<real2d>;
using optional_real3d = std::optional<real3d>;
using optional_real4d = std::optional<real4d>;
using optional_real5d = std::optional<real5d>;

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

int constexpr si_monitor_convergence = 0;
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
