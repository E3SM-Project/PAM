
#ifndef __CONST_H__
#define __CONST_H__

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Kokkos_Core.hpp>

typedef double real;

#ifdef __USE_CUDA__
  typedef KokkosView<real*      , KokkosLayoutRight, KokkosCudaUVMSpace> real1d;
  typedef KokkosView<real**     , KokkosLayoutRight, KokkosCudaUVMSpace> real2d;
  typedef KokkosView<real***    , KokkosLayoutRight, KokkosCudaUVMSpace> real3d;
  typedef KokkosView<real****   , KokkosLayoutRight, KokkosCudaUVMSpace> real4d;
  typedef KokkosView<real*****  , KokkosLayoutRight, KokkosCudaUVMSpace> real5d;
  typedef KokkosView<real****** , KokkosLayoutRight, KokkosCudaUVMSpace> real6d;
  typedef KokkosView<real*******, KokkosLayoutRight, KokkosCudaUVMSpace> real7d;
#else
  typedef KokkosView<real*      , KokkosLayoutRight, KokkosHostSpace> real1d;
  typedef KokkosView<real**     , KokkosLayoutRight, KokkosHostSpace> real2d;
  typedef KokkosView<real***    , KokkosLayoutRight, KokkosHostSpace> real3d;
  typedef KokkosView<real****   , KokkosLayoutRight, KokkosHostSpace> real4d;
  typedef KokkosView<real*****  , KokkosLayoutRight, KokkosHostSpace> real5d;
  typedef KokkosView<real****** , KokkosLayoutRight, KokkosHostSpace> real6d;
  typedef KokkosView<real*******, KokkosLayoutRight, KokkosHostSpace> real7d;
#endif

typedef KokkosView<real*      , KokkosLayoutRight, KokkosHostSpace> hostReal1d;
typedef KokkosView<real**     , KokkosLayoutRight, KokkosHostSpace> hostReal2d;
typedef KokkosView<real***    , KokkosLayoutRight, KokkosHostSpace> hostReal3d;
typedef KokkosView<real****   , KokkosLayoutRight, KokkosHostSpace> hostReal4d;
typedef KokkosView<real*****  , KokkosLayoutRight, KokkosHostSpace> hostReal5d;
typedef KokkosView<real****** , KokkosLayoutRight, KokkosHostSpace> hostReal6d;
typedef KokkosView<real*******, KokkosLayoutRight, KokkosHostSpace> hostReal7d;

int constexpr YES3D = YES3DVAL;   // Domain dimensionality: 1 - 3D, 0 - 2D
int constexpr nx_gl = crm_nx;     // Number of grid points in X
int constexpr ny_gl = crm_ny;     // Number of grid points in Y
int constexpr nz_gl = crm_nz;     // Number of pressure (scalar) levels
int constexpr nsubdomains_x = 1;  // No of subdomains in x
int constexpr nsubdomains_y = 1;  // No of subdomains in y
int constexpr navgmom_x = -1;
int constexpr navgmom_y = -1;
int constexpr ntracers = 0;       // number of transported tracers (dotracers=.true.)

// #ifndef CRM
real constexpr cp    = 1004.;          // Specific heat of air, J/kg/K
real constexpr ggr   = 9.81;           // Gravity acceleration, m/s2
real constexpr lcond = 2.5104e+06;     // Latent heat of condensation, J/kg
real constexpr lfus  = 0.3336e+06;     // Latent heat of fusion, J/kg
real constexpr lsub  = 2.8440e+06;     // Latent heat of sublimation, J/kg
real constexpr rv    = 461.;           // Gas constant for water vapor, J/kg/K
real constexpr rgas  = 287.;           // Gas constant for dry air, J/kg/K
// #else
//   real constexpr  cp    = real( shr_const_cpdair ,crm_rknd)
//   real constexpr  ggr   = real( shr_const_g      ,crm_rknd)
//   real constexpr  lcond = real( shr_const_latvap ,crm_rknd)
//   real constexpr  lfus  = real( shr_const_latice ,crm_rknd)
//   real constexpr  lsub  = real( lcond + lfus     ,crm_rknd)
//   real constexpr  rgas  = real( shr_const_rdair  ,crm_rknd)
//   real constexpr  rv    = real( shr_const_rgas/shr_const_mwwv ,crm_rknd)
// #endif
real constexpr diffelq = 2.21e-05;     // Diffusivity of water vapor, m2/s
real constexpr therco = 2.40e-02;      // Thermal conductivity of air, J/m/s/K
real constexpr muelq = 1.717e-05;      // Dynamic viscosity of air

real constexpr fac_cond = lcond/cp;
real constexpr fac_fus  = lfus/cp;
real constexpr fac_sub  = lsub/cp;

real constexpr pi = 3.141592653589793;

real constexpr ug = 0.;        // Velocity of the Domain's drift in x direction
real constexpr vg = 0.;        // Velocity of the Domain's drift in y direction

int  constexpr les =0;         // flag for Large-Eddy Simulation

int  constexpr sfc_flx_fxd = 0; // surface sensible flux is fixed
int  constexpr sfc_tau_fxd = 0; // surface drag is fixed

int constexpr dodamping        = 1;
int constexpr doprecip         = 1;
int constexpr dosgs            = 1;
int constexpr dosurface        = 1;
int constexpr docloud          = 1;
int constexpr docam_sfc_fluxes = 0;   // Apply the surface fluxes within CAM
int constexpr docoriolis       = 0;
int constexpr dowallx          = 0;
int constexpr dowally          = 0;
int constexpr docolumn         = 0;
int constexpr dotracers        = 0;
int constexpr dosmoke          = 0;

