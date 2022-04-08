#pragma once

#include "common.h"
#include "topology.h"
#include "geometry.h"
#include "stats.h"

// #include "ext_deriv.h"
// #include "hodge_star.h"
// #include "fct.h"
// #include "recon.h"
// #include "wedge.h"
// #include "hamiltonian.h"
// #include "params.h"

// *******   Functionals/Hamiltonians   ***********//
// 
// Functional_PVPE PVPE;
// Hamiltonian_Hk Hk;
// 
// // EVENTUALLY THESE NEED VARIABLESETS, AND A MORE CLEVER UNIFIED HAMILTONIAN IDEALLY?
// 
// #ifdef _SWE
// Hamiltonian_SWE_Hs Hs;
// #elif _TSWE
// Hamiltonian_TSWE_Hs Hs;
// #elif _CE
// Hamiltonian_CE_Hs Hs;
// #elif _CEp
// Hamiltonian_CE_p_Hs Hs;
// #elif _MCErho
// Hamiltonian_MCE_rho_Hs Hs;
// #elif _MCErhop
// Hamiltonian_MCE_rho_p_Hs Hs;
// #elif _MCErhod
// Hamiltonian_MCE_rhod_Hs Hs;
// #elif _MCErhodp
// Hamiltonian_MCE_rhod_p_Hs Hs;
// #endif
// //ADD ANELASTIC + MOIST ANELASTIC
// //MIGHT NEED MORE CLEVER Hk/PVPE HERE? NOT SURE...
// 
// 
// #ifdef _THERMONONE
// ThermoPotential thermo;
// #elif _IDEAL_GAS_POTTEMP
// IdealGas_Pottemp thermo;
// #elif _IDEAL_GAS_ENTROPY
// IdealGas_Entropy thermo;
// #elif _CONST_KAPPA_VIRPOTTEMP
// ConstantKappa_VirtualPottemp thermo;
// #elif _UNAPPROX_POTTEMP
// Unapprox_Pottemp thermo;
// #elif _UNAPPROX_ENTROPY
// Unapprox_Entropy thermo;
// #endif

// *******   Statistics   ***********//

//THIS STUFF SHOULD BE CLEANED UP AND GENERALIZED LIKE VARIABLE SETS IF POSSIBLE...
//ONLY COMPUTE FUNCTION NEEDS TO CHANGE!



// 
// 
// template <uint nprog, uint nconst, uint nstats> class Stats
// {
// public:
//   std::array<Stat,nstats> stats_arr;
//   MPI_Request Req [nstats];
//   MPI_Status  Status[nstats];
//   int ierr;
//   int masterproc;
//   const Topology *primal_topology;
//   const Topology *dual_topology;
//   Geometry<1,1,1> *primal_geometry;
//   Geometry<1,1,1> *dual_geometry;
//   int nens;
// 
//   real3d TEarr, KEarr, PEarr, IEarr, PVarr, PENSarr, trimmed_density;
// 
//   void initialize(Parameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
//   {
//     this->primal_topology = &primal_topo;
//     this->dual_topology = &dual_topo;
//     this->primal_geometry = &primal_geom;
//     this->dual_geometry = &dual_geom;
//     this->nens = params.nens;
//     masterproc = par.masterproc;
// 
//     stats_arr[DENSSTAT].initialize("mass", ndensity, params, par);
//     stats_arr[DENSMAXSTAT].initialize("densmax", ndensity, params, par);
//     stats_arr[DENSMINSTAT].initialize("densmin", ndensity, params, par);
//     stats_arr[ESTAT].initialize("energy", 4, params, par);
//     stats_arr[PVSTAT].initialize("pv", 1, params, par);
//     stats_arr[PESTAT].initialize("pens", 1, params, par);
// 
//     TEarr = real3d("TE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
//     KEarr = real3d("KE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
//     IEarr = real3d("IE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
//     PEarr = real3d("PE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
//     PVarr = real3d("PV", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);
//     PENSarr = real3d("PENS", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);
//     trimmed_density = real3d("trimmed_density", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
// 
//   }
// 
// 
// 
//   void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int tind)
//   {
// 
//     for (int n=0;n<nens;n++)
//     {
// 
//       SArray<real,1,ndensity> masslocal, massglobal;
//       SArray<real,1,ndensity> densmaxlocal, densmaxglobal;
//       SArray<real,1,ndensity> densminlocal, densminglobal;
//       SArray<real,1,1> pvlocal, pvglobal;
//       SArray<real,1,4> elocal, eglobal;
//       SArray<real,1,1> pelocal, peglobal;
// 
// 
//       pvlocal(0) = 0.;
//       pvglobal(0) = 0.;
//       pelocal(0) = 0.;
//       peglobal(0) = 0.;
//       elocal(0) = 0.;
//       elocal(1) = 0.;
//       elocal(2) = 0.;
//       elocal(3) = 0.;
//       eglobal(0) = 0.;
//       eglobal(1) = 0.;
//       eglobal(2) = 0.;
//       eglobal(3) = 0.;
//       for (int l=0;l<ndensity;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
//       for (int l=0;l<ndensity;l++) {densmaxlocal(l) = 0.; densmaxglobal(l) = 0.;}
//       for (int l=0;l<ndensity;l++) {densminlocal(l) = 0.; densminglobal(l) = 0.;}
// 
// int dis = dual_topology->is;
// int djs = dual_topology->js;
// int dks = dual_topology->ks;
// 
//       parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
//          real KE, PE, IE;
// KE = Hk.compute_KE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k, n);
// PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, k, n);
// IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k, n);
// TEarr(k, j, i) = KE + PE + IE;
// KEarr(k, j, i) = KE;
// PEarr(k, j, i) = PE;
// IEarr(k, j, i) = IE;
// });
// 
// elocal(0) = yakl::intrinsics::sum(TEarr);
// elocal(1) = yakl::intrinsics::sum(KEarr);
// elocal(2) = yakl::intrinsics::sum(PEarr);
// elocal(3) = yakl::intrinsics::sum(IEarr);
// 
// int pis = primal_topology->is;
// int pjs = primal_topology->js;
// int pks = primal_topology->ks;
// 
// 
//   parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
//    pvpe vals_pvpe;
//    vals_pvpe = PVPE.compute_PVPE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISVAR].data, pis, pjs, pks, i, j, k, n);
//    PVarr(k, j, i) = vals_pvpe.pv;
//    PENSarr(k, j, i) = vals_pvpe.pe;
//     });
// 
//     pvlocal(0) = yakl::intrinsics::sum(PVarr);
//     pelocal(0) = yakl::intrinsics::sum(PENSarr);
// 
//     for (int l=0;l<ndensity;l++)
//     {
//       parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
//         trimmed_density(k,j,i) = progvars.fields_arr[DENSVAR].data(l,k+dks,j+djs,i+dis,n);
//       });
// 
// 
//     masslocal(l) = yakl::intrinsics::sum(trimmed_density);
//     densmaxlocal(l) = yakl::intrinsics::maxval(trimmed_density);
//     densminlocal(l) = yakl::intrinsics::minval(trimmed_density);
//   }
// 
//     //MPI sum/min/max
//     this->ierr = MPI_Ireduce( &masslocal, &massglobal, ndensity, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSSTAT]);
//     this->ierr = MPI_Ireduce( &densmaxlocal, &densmaxglobal, ndensity, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
//     this->ierr = MPI_Ireduce( &densminlocal, &densminglobal, ndensity, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
//     this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
//     this->ierr = MPI_Ireduce( &pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
//     this->ierr = MPI_Ireduce( &elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);
// 
//     this->ierr = MPI_Waitall(nstats, this->Req, this->Status);
// 
// 
//   if (masterproc)
//   {
//     for (int l=0;l<ndensity;l++)
//     {
//   this->stats_arr[DENSSTAT].data(l,tind,n) = massglobal(l);
//   this->stats_arr[DENSMAXSTAT].data(l,tind,n) = densmaxglobal(l);
//   this->stats_arr[DENSMINSTAT].data(l,tind,n) = densminglobal(l);
// }
// 
//   this->stats_arr[ESTAT].data(0,tind,n) = eglobal(0);
//   this->stats_arr[ESTAT].data(1,tind,n) = eglobal(1);
//   this->stats_arr[ESTAT].data(2,tind,n) = eglobal(2);
//   this->stats_arr[ESTAT].data(3,tind,n) = eglobal(3);
//   this->stats_arr[PVSTAT].data(0,tind,n) = pvglobal(0);
//   this->stats_arr[PESTAT].data(0,tind,n) = peglobal(0);
// 
// }
// 
// }
// }
//};


// *******   VariableSet Initialization   ***********//
template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int,2, nprog, 3> &prog_ndofs_arr, SArray<int,2, nconst, 3> &const_ndofs_arr, SArray<int,2, naux, 3> &aux_ndofs_arr, SArray<int,2, ndiag, 3> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{

  //primal grid represents straight quantities, dual grid twisted quantities
  
  // v, dens
  prog_topo_arr[VVAR] = &ptopo;
  prog_topo_arr[DENSVAR] = &dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[DENSVAR] = "dens";
  set_dofs_arr(prog_ndofs_arr, VVAR, 1, 1, 1); //v = straight 1-form
  set_dofs_arr(prog_ndofs_arr, DENSVAR, ndims, 1, ndensity); //dens = twisted n-form

  // hs, coriolis
  const_topo_arr[HSVAR] = &dtopo;
  const_topo_arr[CORIOLISVAR] = &ptopo;
  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";
  set_dofs_arr(const_ndofs_arr, HSVAR, ndims, 1, 1); //hs = twisted n-form
  set_dofs_arr(const_ndofs_arr, CORIOLISVAR, 2, 1, 1); //f = straight 2-form

  //functional derivatives = F, B, K, he, U
  aux_topo_arr[BVAR] = &ptopo;
  aux_topo_arr[FVAR] = &dtopo;
  aux_topo_arr[UVAR] = &dtopo;
  aux_topo_arr[HEVAR] = &dtopo;
  aux_topo_arr[KVAR] = &dtopo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  set_dofs_arr(aux_ndofs_arr, BVAR, 0, 1, ndensity); //B = straight 0-form
  set_dofs_arr(aux_ndofs_arr, KVAR, ndims, 1, 1);   //K = twisted n-form
  set_dofs_arr(aux_ndofs_arr, FVAR, ndims-1, 1, 1);  //F = twisted (n-1)-form
  set_dofs_arr(aux_ndofs_arr, UVAR, ndims-1, 1, 1); //U = twisted (n-1)-form
  set_dofs_arr(aux_ndofs_arr, HEVAR, ndims-1, 1, 1); //he lives on dual edges, associated with F

  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &dtopo;
  aux_topo_arr[DENSEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENS0VAR] = &ptopo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  set_dofs_arr(aux_ndofs_arr, DENSRECONVAR, ndims-1, 1, ndensity);  //densrecon lives on dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSEDGERECONVAR, ndims, 1, 2*ndims*ndensity); //densedgerecon lives on dual cells, associated with F
  set_dofs_arr(aux_ndofs_arr, DENS0VAR, 0, 1, ndensity); //dens0 = straight 0-form
  
  //dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon, coriolisedgercon, coriolisrecon
  aux_topo_arr[FTVAR] = &ptopo;
  aux_topo_arr[CORIOLISRECONVAR] = &ptopo;
  aux_topo_arr[CORIOLISEDGERECONVAR] = &ptopo; 
  aux_topo_arr[Q0VAR] = &dtopo;
  aux_topo_arr[F0VAR] = &dtopo;
  aux_topo_arr[QRECONVAR] = &ptopo; 
  aux_topo_arr[QEDGERECONVAR] = &ptopo; 
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[CORIOLISRECONVAR] = "coriolisrecon";
  aux_names_arr[CORIOLISEDGERECONVAR] = "coriolisedgerecon";
  aux_names_arr[Q0VAR] = "q";
  aux_names_arr[F0VAR] = "f";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QEDGERECONVAR] = "qedgerecon";
  set_dofs_arr(aux_ndofs_arr, FTVAR, 1, 1, 1); //FT = straight 1-form
  set_dofs_arr(aux_ndofs_arr, Q0VAR, 0, 1, 1);  //q0 = twisted 0-form
  set_dofs_arr(aux_ndofs_arr, F0VAR, 0, 1, 1);  //f0 = twisted 0-form
  set_dofs_arr(aux_ndofs_arr, QRECONVAR, 1, 1, 1);  //qrecon lives on primal edges, associated with FT
  set_dofs_arr(aux_ndofs_arr, QEDGERECONVAR, 2, 1, 4);  //qedgerecon lives on primal cells
  set_dofs_arr(aux_ndofs_arr, CORIOLISRECONVAR, 1, 1, 1);  //coriolisrecon lives on primal edges, associated with FT
  set_dofs_arr(aux_ndofs_arr, CORIOLISEDGERECONVAR, 2, 1, 4);  //coriolisedgerecon lives on primal cells

  // q, concentration 0-forms for dens
  diag_topo_arr[QDIAGVAR] = &dtopo;
  diag_topo_arr[DENSLDIAGVAR] = &ptopo;
  diag_names_arr[QDIAGVAR] = "q";
  diag_names_arr[DENSLDIAGVAR] = "densl";
  set_dofs_arr(diag_ndofs_arr, QDIAGVAR, 0, 1, 1);  //qdiag = twisted 0-form
  set_dofs_arr(diag_ndofs_arr, DENSLDIAGVAR, 0, 1, ndensity); //densldiag = straight 0-form

  //fct stuff- Phi, Mf, edgeflux
  aux_topo_arr[PHIVAR] = &dtopo;
  aux_topo_arr[MFVAR] = &dtopo;
  aux_topo_arr[EDGEFLUXVAR] = &dtopo;
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims-1, 1, ndensity); 
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ndensity);  
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims-1, 1, ndensity); 

  #if defined _AN || defined _MAN
  aux_topo_arr[PVAR] = &ptopo; //p = straight 0-form
  aux_names_arr[PVAR] = "p";
  set_dofs_arr(aux_ndofs_arr, PVAR, 0, 1, 1);  //p = straight 0-form
  #endif
  
}


//***************** Set Initial Conditions ***************************//

struct dbv_constants {
real const g = 9.80616;
real const Lx = 5000. * 1000.;
real const Ly = 5000. * 1000.;
real const coriolis = 0.00006147;
real const H0 = 750.0;
real const ox = 0.1;
real const oy = 0.1;
real const sigmax = 3./40.*Lx;
real const sigmay = 3./40.*Ly;
real const dh = 75.0;
real const xc1 = (0.5-ox) * Lx;
real const yc1 = (0.5-oy) * Ly;
real const xc2 = (0.5+ox) * Lx;
real const yc2 = (0.5+oy) * Ly;
real const xc = 0.5 * Lx;
real const yc = 0.5 * Ly;
real const c = 0.05;
real const a = 1.0/3.0;
real const D = 0.5 * Lx;
};
dbv_constants dbl_vortex_constants;


real YAKL_INLINE double_vortex_coriolis(real x, real y)
{
  return dbl_vortex_constants.coriolis;
}

real YAKL_INLINE double_vortex_h(real x, real y)
{
  real xprime1 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
  real yprime1 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
  real xprime2 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
  real yprime2 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));
  real xprimeprime1 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
  real yprimeprime1 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
  real xprimeprime2 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
  real yprimeprime2 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));

  return dbl_vortex_constants.H0 - dbl_vortex_constants.dh * (exp(-0.5 * (xprime1 * xprime1 + yprime1 * yprime1)) + exp(-0.5 * (xprime2 * xprime2 + yprime2 * yprime2)) - 4. * M_PI * dbl_vortex_constants.sigmax * dbl_vortex_constants.sigmay / dbl_vortex_constants.Lx / dbl_vortex_constants.Ly);
}

vec<2> YAKL_INLINE double_vortex_v(real x, real y) {
vec<2> vvec;

real xprime1 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
real yprime1 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
real xprime2 = dbl_vortex_constants.Lx / (M_PI * dbl_vortex_constants.sigmax) * sin(M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
real yprime2 = dbl_vortex_constants.Ly / (M_PI * dbl_vortex_constants.sigmay) * sin(M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));
real xprimeprime1 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc1));
real yprimeprime1 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc1));
real xprimeprime2 = dbl_vortex_constants.Lx / (2.0 * M_PI * dbl_vortex_constants.sigmax) * sin(2 * M_PI / dbl_vortex_constants.Lx * (x - dbl_vortex_constants.xc2));
real yprimeprime2 = dbl_vortex_constants.Ly / (2.0 * M_PI * dbl_vortex_constants.sigmay) * sin(2 * M_PI / dbl_vortex_constants.Ly * (y - dbl_vortex_constants.yc2));

vvec.u = - dbl_vortex_constants.g * dbl_vortex_constants.dh / dbl_vortex_constants.coriolis / dbl_vortex_constants.sigmay * (yprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + yprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
vvec.v = dbl_vortex_constants.g * dbl_vortex_constants.dh / dbl_vortex_constants.coriolis / dbl_vortex_constants.sigmax * (xprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + xprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
return vvec;
}

real YAKL_INLINE double_vortex_S(real x, real y)
{
  //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. * M_PI / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
  real sval = dbl_vortex_constants.g * (1. + dbl_vortex_constants.c * exp(-((x-dbl_vortex_constants.xc)*(x-dbl_vortex_constants.xc) + (y-dbl_vortex_constants.yc)*(y-dbl_vortex_constants.yc))/(dbl_vortex_constants.a*dbl_vortex_constants.a*dbl_vortex_constants.D*dbl_vortex_constants.D)));
  //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
  //real sval = g;
  //real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35 * Ly && y < 0.65 * Ly ) ? 1. : 0.));
  return sval * double_vortex_h(x,y);
}


// CAN WE GENERALIZE THESE? HOW? NEED TO CALL OUT TO THE CORRECT HEIGHT FUNCTION...
// ALSO NEED TO SET VARIOUS LX/LY SIZES
// MAYBE THIS LATTER IS DOABLE VIA PARAMS?
// SIMILAR WITH GAUSSIANS
// MAYBE WHAT WE DO IS HAVE A TRACER STRUCT, AND THEN SET THE VALUES FOR THE TRACER STRUCTURE ACCORDINGLY?
// AND THE TRACER STRUCT CAN HAVE A LINK TO DENSITY/HEIGHT FUNCTION!
// IN FACT, I THINK ALL TEST CASES STRUCTURES SHOULD INHERIT FROM A MASTER STRUCTURE THAT HAS LX/LY/ETC ALREADY SET UP?
// SOMETHING LIKE THIS...

real YAKL_INLINE tracer_square_cent(real x, real y)         {return (x > 0.35*dbl_vortex_constants.Lx && x < 0.65*dbl_vortex_constants.Lx && y > 0.35*dbl_vortex_constants.Ly && y < 0.65*dbl_vortex_constants.Ly                        ) ? 0.005 : 0.;}
real YAKL_INLINE tracer_square_ur(real x, real y)         {return (x > 0.6*dbl_vortex_constants.Lx && x < 0.9*dbl_vortex_constants.Lx && y > 0.6*dbl_vortex_constants.Ly && y < 0.9*dbl_vortex_constants.Ly                        ) ? 0.005 : 0.;}
real YAKL_INLINE tracer_square_ll(real x, real y)         {return (x > 0.1*dbl_vortex_constants.Lx && x < 0.4*dbl_vortex_constants.Lx && y > 0.1*dbl_vortex_constants.Ly && y < 0.4*dbl_vortex_constants.Ly                        ) ? 0.005 : 0.;}
real YAKL_INLINE tracer_square_urpll(real x, real y)         {return tracer_square_ur(x,y) + tracer_square_ll(x,y);}

real YAKL_INLINE double_vortex_tracer_square_cent (real x, real y) {return double_vortex_h(x,y) * tracer_square_cent(x, y);}
real YAKL_INLINE double_vortex_tracer_square_urpll (real x, real y) {return double_vortex_h(x,y) * tracer_square_urpll(x, y);}

real YAKL_INLINE double_vortex_tracer_gaussian(real x, real y)     {return double_vortex_h(x,y) * 0.005 * exp(-((x-dbl_vortex_constants.xc)*(x-dbl_vortex_constants.xc) + (y-dbl_vortex_constants.yc)*(y-dbl_vortex_constants.yc))/(dbl_vortex_constants.a*dbl_vortex_constants.a*dbl_vortex_constants.D*dbl_vortex_constants.D));}
//{ return 0.005 *double_vortex_h(x,y) * exp(-100. * pow((x-xc)/Lx,2.)) * exp(-100. * pow((y-yc)/Ly,2.)); }



template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (Parameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, 
Geometry<nquadx, nquady, nquadz> &primal_geom, Geometry<nquadx, nquady, nquadz> &dual_geom)
{

  if (params.initdataStr == "doublevortex")
  {
    std::cout << "IC: double vortex " << "\n";

      dual_geom.set_2form_values(double_vortex_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
      dual_geom.set_2form_values(double_vortex_S, progvars.fields_arr[DENSVAR], 1);
#endif
  primal_geom.set_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
  primal_geom.set_2form_values(double_vortex_coriolis, constvars.fields_arr[CORIOLISVAR], 0);

// HOW DO GENERALIZE THESE?
// WANT TO SCALE TRACER FIELDS BY ACTUAL HEIGHT FIELDS...
// SHOULD BE USABLE FOR ANY IC!
for (int i=0; i<ntracers; i++)
{
if (params.tracerdataStr[i] == "gaussian") {dual_geom.set_2form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracerdataStr[i] == "square") {dual_geom.set_2form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
if (params.tracerdataStr[i] == "doublesquare") {dual_geom.set_2form_values(double_vortex_tracer_square_urpll, progvars.fields_arr[DENSVAR], i+ndensity-ntracers);}
}

  }
  
}