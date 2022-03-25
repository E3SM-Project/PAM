#ifndef _MODEL_H_
#define _MODEL_H_

#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "ext_deriv.h"
#include "hodge_star.h"
#include "fct.h"
#include "recon.h"
#include "wedge.h"
#include "geometry.h"
#include "params.h"
//#include "string.h"
#include "hamiltonian.h"


// Number of variables
// v, w, dens, densfct
uint constexpr nprognostic = 3;
#define VVAR 0
#define WVAR 1
#define DENSVAR 2

// hs, coriolis
uint constexpr nconstant = 2;
#define HSVAR 0
#define CORIOLISXZVAR 1

//functional derivatives = F, FW, B, K, he, hew
//primal grid reconstruction stuff- U, W, dens0, edgerecon, recon, vertedgerecon, vertrecon
//fct stuff- Phi, Mf, edgeflux
//Q/W STUFF?

uint constexpr nauxiliary = 32;

#define FVAR 0
#define BVAR 1
#define KVAR 2
#define HEVAR 3
#define UVAR 4
#define FWVAR 5
#define HEWVAR 6
#define UWVAR 7

#define DENS0VAR 8
#define DENSRECONVAR 9
#define DENSEDGERECONVAR 10
#define DENSVERTRECONVAR 11
#define DENSVERTEDGERECONVAR 12

#define PHIVAR 13
#define PHIVERTVAR 14
#define EDGEFLUXVAR 15
#define VERTEDGEFLUXVAR 16
#define MFVAR 17

#define QXZ0VAR 18
#define QXZRECONVAR 19
#define QXZVERTRECONVAR 20
#define QXZEDGERECONVAR 21
#define QXZVERTEDGERECONVAR 22
#define QXZFLUXVAR 23
#define QXZVERTFLUXVAR 24

#define FTVAR 25
#define FTWVAR 26

#define FXZ0VAR 27
#define CORIOLISXZRECONVAR 28
#define CORIOLISXZEDGERECONVAR 29
#define CORIOLISXZVERTRECONVAR 30
#define CORIOLISXZVERTEDGERECONVAR 31

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 9;
#define DENSLDIAGVAR 0
#define QXZDIAGVAR 1
#define FDIAGVAR 2
#define KDIAGVAR 3
#define HEDIAGVAR 4
#define FWDIAGVAR 5
#define HEWDIAGVAR 6
#define UDIAGVAR 7
#define UWDIAGVAR 8

//track total densities, dens min/max, energy (total, K, P, I), PV, PE
uint constexpr nstats = 6;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define PVSTAT 4
#define PESTAT 5

// *******   Functionals/Hamiltonians   ***********//


Functional_PVPE_extruded PVPE;
Hamiltonian_Hk_extruded Hk;
// ANELASTIC MIGHT POSSIBLY NEED SEPARATE FUNCTIONALS HERE TO ACCOUNT FOR CONSTANT DENSITY PROFILE?
// OR REALLY THE NUMERICS BETTER CONSERVE IT TO MACHINE PRECISION, SO MAYBE IT IS OK AS A FIELD?

// NEED TO SET WHERE/HOW TOTAL DENSITY IS SET

#ifdef _SWE
Hamiltonian_SWE_Hs Hs;
ThermoPotential thermo;
#endif

#ifdef _TSWE
Hamiltonian_TSWE_Hs Hs;
ThermoPotential thermo;
#endif

#ifdef _CE
#ifdef _USE_P_VARIANT
Hamiltonian_CE_p_Hs Hs;
#else
Hamiltonian_CE_Hs Hs;
#endif
#ifdef _IDEAL_GAS_POTTEMP
IdealGas_Pottemp thermo;
#endif
#ifdef _IDEAL_GAS_ENTROPY
IdealGas_Entropy thermo;
#endif
#endif


#ifdef _AN
Hamiltonian_AN_Hs Hs;
#ifdef _IDEAL_GAS_POTTEMP
IdealGas_Pottemp thermo;
#endif
#ifdef _IDEAL_GAS_ENTROPY
IdealGas_Entropy thermo;
#endif
#endif


#ifdef _MCE
#ifdef _USE_RHOD_VARIANT
#ifdef _USE_P_VARIANT
Hamiltonian_MCE_rhod_p_Hs Hs;
#else
Hamiltonian_MCE_rhod_Hs Hs;
#endif
#else
#ifdef _USE_P_VARIANT
Hamiltonian_MCE_rho_p_Hs Hs;
#else
Hamiltonian_MCE_rho_Hs Hs;
#endif
#endif
#ifdef _CONST_KAPPA_VIRPOTTEMP
ConstantKappa_VirtualPottemp thermo;
#endif
#ifdef _CONST_KAPPA_ENTROPY
ConstantKappa_Entropy thermo;
#endif
#ifdef _UNAPPROX_POTTEMP
Unapprox_Pottemp thermo;
#endif
#ifdef _UNAPPROX_ENTROPY
Unapprox_Entropy thermo;
#endif
#endif


#ifdef _MAN
Hamiltonian_MAN_Hs Hs;
#ifdef _CONST_KAPPA_VIRPOTTEMP
ConstantKappa_VirtualPottemp thermo;
#endif
#ifdef _CONST_KAPPA_ENTROPY
ConstantKappa_Entropy thermo;
#endif
#ifdef _UNAPPROX_POTTEMP
Unapprox_Pottemp thermo;
#endif
#ifdef _UNAPPROX_ENTROPY
Unapprox_Entropy thermo;
#endif
#endif



// *******   Model Specific Parameters   ***********//

void set_model_specific_params(std::string inFile, ModelParameters &params)
{

  params.etime = 0.0;

  // Read in equals-separated key = value file line by line
  std::ifstream fInStream(inFile);
  std::string line;
  while (std::getline(fInStream, line)) {
    // Remove spaces and tabs from the line
    line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
    line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());

    // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
    if (!line.empty() && line.find("//",0) != 0) {
      // Find the colon
      uint splitloc = line.find('=',0);
      // Store the key and value strings
      std::string key   = line.substr(0,splitloc);
      std::string value = line.substr(splitloc+1,line.length()-splitloc);

      // Transform the value into a string stream for convenience
      std::stringstream ssVal(value);

      // Match the key, and store the value
      if      ( !strcmp( "acoustic_balance"         , key.c_str() ) ) { ssVal >> params.acoustic_balance    ; }
      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }
  
}

// ******* Diagnostics *************//

template <uint nprog, uint nconst, uint ndiag> class Diagnostics {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;
  ExchangeSet<ndiag> *diag_exchange;

  bool is_initialized;

   Diagnostics() {
     this->is_initialized = false;
     std::cout << "CREATED DIAGNOSTICS\n";
   }

   void initialize(const Topology &ptopo, const Topology &dtopo, Geometry<1,1,1> &pgeom, Geometry<1,1,1> &dgeom, ExchangeSet<ndiag> &diag_exchange)
   {
     this->primal_topology = &ptopo;
     this->dual_topology = &dtopo;
     this->primal_geometry = &pgeom;
     this->dual_geometry = &dgeom;
     this->diag_exchange = &diag_exchange;

     this->is_initialized = true;
   }



   void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
   {

     int dis = dual_topology->is;
     int djs = dual_topology->js;
     int dks = dual_topology->ks;

     int pis = primal_topology->is;
     int pjs = primal_topology->js;
     int pks = primal_topology->ks;
     
    // yakl::parallel_for("ComputeDiagIp", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   SArray<real,1> zeta;
    //   real hv;
    //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
// //LOOP CORRECT?


 // ARE THESE REALLY DUAL FORMS?
 // MAYBE THEY ARE STRAIGHT 0-FORMS?

 //THIS RELIES ON DENSITY BEING IN DENSVAR0
 //IS THERE ANOTHER/BETTER WAY?
 //YES- CAN USE THE STATIC INDEXING ARRAY THAT GIVES TOTAL DENSITY AS A SUM OF CERTAIN ELEMENTS OF DENSVAR...
 
 //OTHER DIAGNOSTICS?

 //THIS CAN REALLY BE DONE OFFLINE
 //MAYBE DO HODGE STARS HERE?
 // compute dens0var
 //for(int l=0; l<ndensity; l++)
 //{diagnostic_vars.fields_arr[DENSLDIAGVAR].data(l, k, j+djs, i+dis) = x.fields_arr[DENSVAR].data(l, k+dks, j+djs, i+dis) / x.fields_arr[DENSVAR].data(0, k+dks, j+djs, i+dis);}

// yakl::parallel_for("ComputeDens0", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
   parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
 compute_Iext<ndensity, diff_ord, vert_diff_ord>(diagnostic_vars.fields_arr[DENSLDIAGVAR].data, x.fields_arr[DENSVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
 });


 // yakl::parallel_for("ComputeUVAR", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
 //   int k, j, i;
 //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
 // compute_Hext<1, diff_ord>(diagnostic_vars.fields_arr[UDIAGVAR].data, x.fields_arr[VVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);
 // });
 
 // yakl::parallel_for("ComputeUWVAR", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
 //   int k, j, i;
 //   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
 // compute_Hv<1, vert_diff_ord>(diagnostic_vars.fields_arr[UWDIAGVAR].data, x.fields_arr[WVAR].data, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k+1);
 // });    
 // 
 // yakl::parallel_for("ComputeQxz0VAR", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
 //   int k, j, i;
 //   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
 // PVPE.compute_qxz0(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1);
 // });
 // 
 //  this->diag_exchange->exchanges_arr[UDIAGVAR].exchange_field(diagnostic_vars.fields_arr[UDIAGVAR]);
 //  this->diag_exchange->exchanges_arr[UWDIAGVAR].exchange_field(diagnostic_vars.fields_arr[UWDIAGVAR]);
 // this->diag_exchange->exchanges_arr[DENSLDIAGVAR].exchange_field(diagnostic_vars.fields_arr[DENSLDIAGVAR]);
 // 
 // yakl::parallel_for("ComputeDiagII", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
 //   int k, j, i;
 //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
 // //Hk.compute_F(diagnostic_vars.fields_arr[FDIAGVAR].data, diagnostic_vars.fields_arr[HEDIAGVAR].data, diagnostic_vars.fields_arr[UDIAGVAR].data, diagnostic_vars.fields_arr[DENSLDIAGVAR].data, dis, djs, dks, i, j, k);
 // Hk.compute_K(diagnostic_vars.fields_arr[KDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, diagnostic_vars.fields_arr[UDIAGVAR].data, diagnostic_vars.fields_arr[UWDIAGVAR].data, dis, djs, dks, i, j, k);
 // });
 // yakl::parallel_for("ComputeDiagII", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
 //  int k, j, i;
 //  yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
 //  Hk.compute_Fw(diagnostic_vars.fields_arr[FWDIAGVAR].data, diagnostic_vars.fields_arr[HEWDIAGVAR].data, diagnostic_vars.fields_arr[UWDIAGVAR].data, diagnostic_vars.fields_arr[DENSLDIAGVAR].data, dis, djs, dks, i, j, k+1);
 //  });
 // 


// yakl::parallel_for("ComputeQxz0VAR", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
   parallel_for( Bounds<3>( dual_topology->ni-3, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
 PVPE.compute_qxz0(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data, dis, djs, dks, i, j, k+2);
 });

//Bottom is k=1 and top is k=dual_topology->ni-2
 parallel_for( Bounds<2>(dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
PVPE.compute_qxz0_bottom(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data, dis, djs, dks, i, j, 1);
PVPE.compute_qxz0_top(diagnostic_vars.fields_arr[QXZDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data, dis, djs, dks, i, j, dual_topology->ni-2);
});

   diagnostic_vars.fields_arr[QXZDIAGVAR].set_bnd(0.0);
}

};
// *******   Tendencies   ***********//

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE
template <uint nprog, uint nconst, uint naux> class Tendencies {
public:

  const Topology *primal_topology;
  const Topology *dual_topology;
  ExchangeSet<naux> *aux_exchange;
  ExchangeSet<nconst> *const_exchange;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;

  //TransformMatrices<real> trans;

  SArray<real,2,reconstruction_order,2> primal_to_gll;
  SArray<real,3,reconstruction_order,reconstruction_order,reconstruction_order> primal_wenoRecon;
  SArray<real,1,(reconstruction_order-1)/2+2> primal_wenoIdl;
  real primal_wenoSigma;

  SArray<real,2,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,3,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dual_wenoRecon;
  SArray<real,1,(dual_reconstruction_order-1)/2+2> dual_wenoIdl;
  real dual_wenoSigma;

  SArray<real,2,coriolis_reconstruction_order,2> coriolis_to_gll;
  SArray<real,3,coriolis_reconstruction_order,coriolis_reconstruction_order,coriolis_reconstruction_order> coriolis_wenoRecon;
  SArray<real,1,(coriolis_reconstruction_order-1)/2+2> coriolis_wenoIdl;
  real coriolis_wenoSigma;

  SArray<real,2,vert_reconstruction_order,2> primal_vert_to_gll;
  SArray<real,3,vert_reconstruction_order,vert_reconstruction_order,vert_reconstruction_order> primal_vert_wenoRecon;
  SArray<real,1,(vert_reconstruction_order-1)/2+2> primal_vert_wenoIdl;
  real primal_vert_wenoSigma;

  SArray<real,2,dual_vert_reconstruction_order,2> dual_vert_to_gll;
  SArray<real,3,dual_vert_reconstruction_order,dual_vert_reconstruction_order,dual_vert_reconstruction_order> dual_vert_wenoRecon;
  SArray<real,1,(dual_vert_reconstruction_order-1)/2+2> dual_vert_wenoIdl;
  real dual_vert_wenoSigma;

  SArray<real,2,coriolis_vert_reconstruction_order,2> coriolis_vert_to_gll;
  SArray<real,3,coriolis_vert_reconstruction_order,coriolis_vert_reconstruction_order,coriolis_vert_reconstruction_order> coriolis_vert_wenoRecon;
  SArray<real,1,(coriolis_vert_reconstruction_order-1)/2+2> coriolis_vert_wenoIdl;
  real coriolis_vert_wenoSigma;
  
  bool is_initialized;
  
   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

   void initialize(ModelParameters &params, const Topology &primal_topo, const Topology &dual_topo, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom, ExchangeSet<naux> &aux_exchange, ExchangeSet<nconst> &const_exchange)
   {
     this->primal_topology = &primal_topo;
     this->dual_topology = &dual_topo;
     this->primal_geometry = &primal_geom;
     this->dual_geometry = &dual_geom;
     this->aux_exchange = &aux_exchange;
     this->const_exchange = &const_exchange;

     TransformMatrices::coefs_to_gll_lower( primal_to_gll );
     TransformMatrices::weno_sten_to_coefs(primal_wenoRecon);
     wenoSetIdealSigma<reconstruction_order>(primal_wenoIdl,primal_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( dual_to_gll );
     TransformMatrices::weno_sten_to_coefs(dual_wenoRecon);
     wenoSetIdealSigma<dual_reconstruction_order>(dual_wenoIdl,dual_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( coriolis_to_gll );
     TransformMatrices::weno_sten_to_coefs(coriolis_wenoRecon);
     wenoSetIdealSigma<coriolis_reconstruction_order>(coriolis_wenoIdl,coriolis_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( primal_vert_to_gll );
     TransformMatrices::weno_sten_to_coefs(primal_vert_wenoRecon);
     wenoSetIdealSigma<vert_reconstruction_order>(primal_vert_wenoIdl,primal_vert_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( dual_vert_to_gll );
     TransformMatrices::weno_sten_to_coefs(dual_vert_wenoRecon);
     wenoSetIdealSigma<dual_vert_reconstruction_order>(dual_vert_wenoIdl,dual_vert_wenoSigma);

     TransformMatrices::coefs_to_gll_lower( coriolis_vert_to_gll );
     TransformMatrices::weno_sten_to_coefs(coriolis_vert_wenoRecon);
     wenoSetIdealSigma<coriolis_vert_reconstruction_order>(coriolis_vert_wenoIdl,coriolis_vert_wenoSigma);
     
     PVPE.initialize(params);
     Hk.initialize(params, *this->primal_geometry, *this->dual_geometry);
     Hs.initialize(params, thermo, *this->primal_geometry, *this->dual_geometry);
    
    this->is_initialized = true;
  }


  void compute_constants(VariableSet<nconst> &const_vars, VariableSet<nprog> &x)
  {}
 
   void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
    real4d Uvar, real4d UWvar, real4d qxz0var, real4d fxz0var, real4d dens0var, 
    const real4d Vvar, const real4d Wvar, const real4d densvar, const real4d coriolisxzvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;


// Compute DENS0
//yakl::parallel_for("ComputeDens0VAR", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//  int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
   parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Iext<ndensity, diff_ord, vert_diff_ord>(dens0var, densvar, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
});

// yakl::parallel_for("ComputeUVAR", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Hext<1, diff_ord>(Uvar, Vvar, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);
});

// yakl::parallel_for("ComputeUWVAR", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Hv<1, vert_diff_ord>(UWvar, Wvar, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k+1);
});    

// yakl::parallel_for("ComputeQxz0VAR", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->ni-3, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
PVPE.compute_qxz0fxz0(qxz0var, fxz0var, Vvar, Wvar, densvar, coriolisxzvar, dis, djs, dks, i, j, k+2);
});
parallel_for( Bounds<2>( dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
PVPE.compute_qxz0fxz0_bottom(qxz0var, fxz0var, Vvar, Wvar, densvar, coriolisxzvar, dis, djs, dks, i, j, 1);
PVPE.compute_qxz0fxz0_top(qxz0var, fxz0var, Vvar, Wvar, densvar, coriolisxzvar, dis, djs, dks, i, j,  dual_topology->ni-2);
});
    }


    void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
      real4d Fvar, real4d FWvar, real4d Kvar, real4d HEvar, real4d HEWvar, const real4d Vvar, const real4d Uvar, const real4d Wvar, const real4d UWvar, const real4d dens0var) {

        int dis = dual_topology->is;
        int djs = dual_topology->js;
        int dks = dual_topology->ks;

        int pis = primal_topology->is;
        int pjs = primal_topology->js;
        int pks = primal_topology->ks;
        
  //THIS WILL NEED SOME SLIGHT MODIFICATIONS FOR CASE OF NON-ZERO UWVAR_B IE BOUNDARY FLUXES
  //BUT FOR NOW IT IS FINE SINCE UWVAR=0 on BND AND THEREFORE K COMPUTATIONS IGNORE IT
        // yakl::parallel_for("ComputeDiagII", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
        //   int k, j, i;
        //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
          parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        Hk.compute_F(Fvar, HEvar, Uvar, dens0var, dis, djs, dks, i, j, k);
        Hk.compute_K(Kvar, Vvar, Uvar, Wvar, UWvar, dis, djs, dks, i, j, k);
      });
      // yakl::parallel_for("ComputeDiagII", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
        parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
      Hk.compute_Fw(FWvar, HEWvar, UWvar, dens0var, dis, djs, dks, i, j, k+1);
    });
    
      }

  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    real4d Bvar, real4d FTvar, real4d FTWvar,
    const real4d Fvar, const real4d Uvar, const real4d FWvar, const real4d UWvar,
    const real4d Kvar, const real4d dens0var, const real4d HSvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

//compute FT = "Wxz" Fw ie Fw at v
//compute FTW = "Wxz" F ie F at w
//yakl::parallel_for("ComputeFTVAR", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//  int k, j, i;
//  yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( primal_topology->ni-2, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Wxz_u(FTvar, FWvar, pis, pjs, pks, i, j, k+1);
});
parallel_for( Bounds<2>( primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
compute_Wxz_u_bottom(FTvar, FWvar, pis, pjs, pks, i, j, 0);
compute_Wxz_u_top(FTvar, FWvar, pis, pjs, pks, i, j, primal_topology->ni-1);
});
//yakl::parallel_for("ComputeFTWVAR", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//  int k, j, i;
//  yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( primal_topology->nl-2, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Wxz_w(FTWvar, Fvar, pis, pjs, pks, i, j, k+1);
});
parallel_for( Bounds<2>( primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
compute_Wxz_w_bottom(FTWvar, Fvar, pis, pjs, pks, i, j, 0);
compute_Wxz_w_top(FTWvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1);
});
//      yakl::parallel_for("ComputeDiagIII", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//        int k, j, i;
//        yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
        parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
Hs.compute_dHsdx(Bvar, dens0var, HSvar, pis, pjs, pks, i, j, k);
Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j, k);
    });

    }


  void YAKL_INLINE compute_edge_reconstructions(real4d densedgereconvar, real4d densvertedgereconvar, 
    real4d qxzedgereconvar, real4d qxzvertedgereconvar, 
    real4d coriolisxzedgereconvar, real4d coriolisxzvertedgereconvar,
    const real4d dens0var, const real4d qxz0var, const real4d fxz0var) {

      int dis = dual_topology->is;
      int djs = dual_topology->js;
      int dks = dual_topology->ks;

      int pis = primal_topology->is;
      int pjs = primal_topology->js;
      int pks = primal_topology->ks;
      
      // yakl::parallel_for("ComputeDensEdgeRecon", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
       parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

        compute_twisted_edge_recon<ndensity, dual_reconstruction_type, dual_reconstruction_order>(
          densedgereconvar, dens0var, dis, djs, dks, i, j, k,
          dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
          
          compute_twisted_vert_edge_recon<ndensity, dual_vert_reconstruction_type, dual_vert_reconstruction_order>(
            densvertedgereconvar, dens0var, dis, djs, dks, i, j, k, 
            dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl, dual_vert_wenoSigma);
        });

      //  yakl::parallel_for("ComputeQEdgeRecon", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //    int k, j, i;
        //  yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
          parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

           compute_straight_xz_edge_recon<1, reconstruction_type, reconstruction_order>(
             qxzedgereconvar, qxz0var, pis, pjs, pks, i, j, k,
             primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);

           compute_straight_xz_vert_edge_recon<1, vert_reconstruction_type, vert_reconstruction_order>(
             qxzvertedgereconvar, qxz0var, pis, pjs, pks, i, j, k,
             primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl, primal_vert_wenoSigma);
  
             compute_straight_xz_edge_recon<1, coriolis_reconstruction_type, coriolis_reconstruction_order>(
               coriolisxzedgereconvar, fxz0var, pis, pjs, pks, i, j, k,
               coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl, coriolis_wenoSigma);

             compute_straight_xz_vert_edge_recon<1, coriolis_vert_reconstruction_type, coriolis_vert_reconstruction_order>(
               coriolisxzvertedgereconvar, fxz0var, pis, pjs, pks, i, j, k,
               coriolis_vert_wenoRecon, coriolis_vert_to_gll, coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
               
        });
        
  }
  
  void YAKL_INLINE compute_recons(
  real4d densreconvar,  real4d densvertreconvar,
  real4d qxzreconvar,  real4d qxzvertreconvar,
  real4d coriolisxzreconvar,  real4d coriolisxzvertreconvar,
  const real4d densedgereconvar, const real4d densvertedgereconvar, 
  const real4d qxzedgereconvar, const real4d qxzvertedgereconvar,
  const real4d coriolisxzedgereconvar, const real4d coriolisxzvertedgereconvar,
  const real4d HEvar, const real4d HEWvar,
  const real4d Uvar, const real4d UWvar,
  const real4d FTvar, const real4d FTWvar) {

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

// yakl::parallel_for("ComputeDensRECON", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

  compute_twisted_recon<ndensity, dual_reconstruction_type>(
    densreconvar, densedgereconvar, Uvar, dis, djs, dks, i, j, k);
    //scale twisted recons
    for (int d=0;d<ndims;d++) {
      //std::cout << "HE in recon " << i << " " << j << " " << k << " " << HEvar(d,k+dks,j+djs,i+dis) << "\n" << std::flush;
    for (int l=0;l<ndensity;l++) {
    densreconvar(l+d*ndensity,k+dks,j+djs,i+dis) = densreconvar(l+d*ndensity,k+dks,j+djs,i+dis) / HEvar(d,k+dks,j+djs,i+dis);
  }}
  });

    // yakl::parallel_for("ComputeDensVertRECON", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

    compute_twisted_vert_recon<ndensity, dual_vert_reconstruction_type>(
      densvertreconvar, densvertedgereconvar, UWvar, dis, djs, dks, i, j, k+1);
      //scale twisted recons
      for (int l=0;l<ndensity;l++) {
      //  std::cout << "HEW in recon " << i << " " << j << " " << k << " " << HEWvar(0,k+dks+1,j+djs,i+dis) << "\n" << std::flush;
      densvertreconvar(l,k+dks+1,j+djs,i+dis) = densvertreconvar(l,k+dks+1,j+djs,i+dis) / HEWvar(0,k+dks+1,j+djs,i+dis);
    }
});

// yakl::parallel_for("ComputeQRECON", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

  compute_straight_xz_recon<1, reconstruction_type>(
    qxzreconvar, qxzedgereconvar, FTWvar, pis, pjs, pks, i, j, k);
    
    compute_straight_xz_recon<1, coriolis_reconstruction_type>(
      coriolisxzreconvar, coriolisxzedgereconvar, FTWvar, pis, pjs, pks, i, j, k);
});
// yakl::parallel_for("ComputeQVERTRECON", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_straight_xz_vert_recon<1, vert_reconstruction_type>(
    qxzvertreconvar, qxzvertedgereconvar, FTvar, pis, pjs, pks, i, j, k);
    
    compute_straight_xz_vert_recon<1, coriolis_vert_reconstruction_type>(
      coriolisxzvertreconvar, coriolisxzvertedgereconvar, FTvar, pis, pjs, pks, i, j, k);
});

}

  void YAKL_INLINE compute_tendencies(
  real4d denstendvar, real4d Vtendvar, real4d Wtendvar,
  const real4d densreconvar, const real4d densvertreconvar,
  const real4d qxzreconvar, const real4d qxzvertreconvar,
  const real4d coriolisxzreconvar, const real4d coriolisxzvertreconvar,
  const real4d Bvar, const real4d Fvar, const real4d FWvar, const real4d Phivar, const real4d Phivertvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

// yakl::parallel_for("ComputePrimalTendencies", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( primal_topology->nl-2, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_wDv_fct<ndensity> (Wtendvar, densvertreconvar, Phivertvar, Bvar, pis, pjs, pks, i, j, k+1);
   if (qf_choice == QF_MODE::EC)
  { compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i, j, k+1);}
   if (qf_choice == QF_MODE::NOEC)
  { compute_Qxz_w_nonEC<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, k+1);}
  compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis, pjs, pks, i, j, k+1);
});
parallel_for( Bounds<2>(primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
compute_wDv_fct<ndensity> (Wtendvar, densvertreconvar, Phivertvar, Bvar, pis, pjs, pks, i, j, 0);
compute_wDv_fct<ndensity> (Wtendvar, densvertreconvar, Phivertvar, Bvar, pis, pjs, pks, i, j, primal_topology->nl-1);
 if (qf_choice == QF_MODE::EC)
{ compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i, j, 0);
  compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1);}
 if (qf_choice == QF_MODE::NOEC)
{ compute_Qxz_w_nonEC_bottom<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, 0);
  compute_Qxz_w_nonEC_top<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1);}
compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis, pjs, pks, i, j, 0);
compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis, pjs, pks, i, j, primal_topology->nl-1);
});
// yakl::parallel_for("ComputePrimalTendencies", primal_topology->n_cells_interfaces, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->ni, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( primal_topology->ni-2, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_wD1_fct<ndensity> (Vtendvar, densreconvar, Phivar, Bvar, pis, pjs, pks, i, j, k+1);
  if (qf_choice == QF_MODE::EC)
  { compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar,qxzreconvar, qxzvertreconvar,FWvar, pis, pjs, pks, i, j, k+1);}
  if (qf_choice == QF_MODE::NOEC)
  { compute_Qxz_u_nonEC<1, ADD_MODE::ADD>(Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j, k+1);}
  compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar,coriolisxzreconvar, coriolisxzvertreconvar,FWvar, pis, pjs, pks, i, j, k+1);
});
parallel_for( Bounds<2>( primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
compute_wD1_fct<ndensity> (Vtendvar, densreconvar, Phivar, Bvar, pis, pjs, pks, i, j, 0);
compute_wD1_fct<ndensity> (Vtendvar, densreconvar, Phivar, Bvar, pis, pjs, pks, i, j, primal_topology->ni-1);
if (qf_choice == QF_MODE::EC)
{ compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(Vtendvar,qxzreconvar, qxzvertreconvar,FWvar, pis, pjs, pks, i, j, 0);
  compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(Vtendvar,qxzreconvar, qxzvertreconvar,FWvar, pis, pjs, pks, i, j, primal_topology->ni-1);}
if (qf_choice == QF_MODE::NOEC)
{ compute_Qxz_u_nonEC_bottom<1, ADD_MODE::ADD>(Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j, 0);
  compute_Qxz_u_nonEC_top<1, ADD_MODE::ADD>(Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j, primal_topology->ni-1);}
compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(Vtendvar,coriolisxzreconvar, coriolisxzvertreconvar,FWvar, pis, pjs, pks, i, j, 0);
compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(Vtendvar,coriolisxzreconvar, coriolisxzvertreconvar,FWvar, pis, pjs, pks, i, j, primal_topology->ni-1);
});

// yakl::parallel_for("ComputeDualTendencies", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_wDbar2_fct<ndensity> (denstendvar, densreconvar, Phivar, Fvar, dis, djs, dks, i, j, k);
  compute_wDvbar_fct<ndensity, ADD_MODE::ADD> (denstendvar, densvertreconvar, Phivertvar, FWvar, dis, djs, dks, i, j, k);
  });

  }



  
  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {
       
      //Compute U, W, q0, dens0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data, auxiliary_vars.fields_arr[QXZ0VAR].data, auxiliary_vars.fields_arr[FXZ0VAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISXZVAR].data);
  
      auxiliary_vars.fields_arr[QXZ0VAR].set_bnd(0.0);
      auxiliary_vars.fields_arr[FXZ0VAR].set_bnd(0.0);
      auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
      
      this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
      this->aux_exchange->exchanges_arr[UWVAR].exchange_field(auxiliary_vars.fields_arr[UWVAR]);
      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);
      this->aux_exchange->exchanges_arr[QXZ0VAR].exchange_field(auxiliary_vars.fields_arr[QXZ0VAR]);
      this->aux_exchange->exchanges_arr[FXZ0VAR].exchange_field(auxiliary_vars.fields_arr[FXZ0VAR]);

      //Compute K, F, FW, he, hew
      compute_functional_derivatives_and_diagnostic_quantities_II(
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[HEWVAR].data,
      x.fields_arr[VVAR].data, auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[WVAR].data, auxiliary_vars.fields_arr[UWVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data);

      //std::cout << "HE min/max " << auxiliary_vars.fields_arr[HEVAR].min() << " " << auxiliary_vars.fields_arr[HEVAR].max() << " " << auxiliary_vars.fields_arr[HEVAR].sum() << "\n" <<std::flush;
      //std::cout << "HEw min/max " << auxiliary_vars.fields_arr[HEWVAR].min() << " " << auxiliary_vars.fields_arr[HEWVAR].max() << " " << auxiliary_vars.fields_arr[HEWVAR].sum() << "\n" <<std::flush;

      auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);

      this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
      this->aux_exchange->exchanges_arr[FWVAR].exchange_field(auxiliary_vars.fields_arr[FWVAR]);
      this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
      this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);
      this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(auxiliary_vars.fields_arr[HEWVAR]);
      
      //std::cout << "HE min/max " << auxiliary_vars.fields_arr[HEVAR].min() << " " << auxiliary_vars.fields_arr[HEVAR].max() << " " << auxiliary_vars.fields_arr[HEVAR].sum() << "\n" <<std::flush;
      //std::cout << "HEw min/max " << auxiliary_vars.fields_arr[HEWVAR].min() << " " << auxiliary_vars.fields_arr[HEWVAR].max() << " " << auxiliary_vars.fields_arr[HEWVAR].sum() << "\n" <<std::flush;
      
      //Compute FT, B
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[BVAR].data,  auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[FTWVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars.fields_arr[UWVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, const_vars.fields_arr[HSVAR].data);

      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);
      this->aux_exchange->exchanges_arr[FTVAR].exchange_field(auxiliary_vars.fields_arr[FTVAR]);
      this->aux_exchange->exchanges_arr[FTWVAR].exchange_field(auxiliary_vars.fields_arr[FTWVAR]);
          
      // Compute densrecon, densvertrecon, qrecon and frecon
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data, 
      auxiliary_vars.fields_arr[QXZEDGERECONVAR].data, auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[QXZ0VAR].data, auxiliary_vars.fields_arr[FXZ0VAR].data);

      this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[DENSVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[QXZEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[QXZEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISXZVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISXZVERTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR]);

      compute_recons(
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data, 
      auxiliary_vars.fields_arr[QXZRECONVAR].data, auxiliary_vars.fields_arr[QXZVERTRECONVAR].data, 
      auxiliary_vars.fields_arr[CORIOLISXZRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR].data, 
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[QXZEDGERECONVAR].data, auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[HEWVAR].data,
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[UWVAR].data,
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[FTWVAR].data);
      
      this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);
      this->aux_exchange->exchanges_arr[DENSVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSVERTRECONVAR]);
      this->aux_exchange->exchanges_arr[QXZRECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZRECONVAR]);
      this->aux_exchange->exchanges_arr[QXZVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[QXZVERTRECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISXZRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZRECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISXZVERTRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR]);


//Compute fct
int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

// yakl::parallel_for("ComputeEdgeFlux", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

compute_edgefluxes<ndensity> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, dis, djs, dks, i, j, k);
});
// yakl::parallel_for("ComputeVertEdgeFlux", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_vertedgefluxes<ndensity> (auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data, auxiliary_vars.fields_arr[FWVAR].data, dis, djs, dks, i, j, k);
});

this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);
this->aux_exchange->exchanges_arr[VERTEDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[VERTEDGEFLUXVAR]);


// yakl::parallel_for("ComputeMf", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
// int k, j, i;
// yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Mfext<ndensity> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, dt, dis, djs, dks, i, j, k);
});

this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

// yakl::parallel_for("ComputePhiVert", dual_topology->n_cells_interfaces_internal, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->ni-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_Phivert<ndensity> (auxiliary_vars.fields_arr[PHIVERTVAR].data, auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1);
});

// yakl::parallel_for("ComputePhi", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_Phi<ndensity> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k);
});

//Don't do FCT for non-FCT vars
    for (int l=0; l<ndensity_nofct; l++)
    {
    auxiliary_vars.fields_arr[PHIVAR].set(l, 1.0);
    auxiliary_vars.fields_arr[PHIVERTVAR].set(l, 1.0);
    }
    
this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);
this->aux_exchange->exchanges_arr[PHIVERTVAR].exchange_field(auxiliary_vars.fields_arr[PHIVERTVAR]);

      // Compute tendencies
      compute_tendencies(
      xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data, xtend.fields_arr[WVAR].data,
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
      auxiliary_vars.fields_arr[QXZRECONVAR].data, auxiliary_vars.fields_arr[QXZVERTRECONVAR].data,
      auxiliary_vars.fields_arr[CORIOLISXZRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR].data,
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[FWVAR].data, 
      auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[PHIVERTVAR].data);
}

};



// *******   Statistics Calculations   ***********//

// THIS STUFF SHOULD BE CLEANED UP AND GENERALIZED LIKE VARIABLE SETS IF POSSIBLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE!
class Stat
{
public:
  real2d data;
  std::string name;
  int ndofs;


void initialize(std::string statName, int ndof, ModelParameters &params, Parallel &par)
{
  name = statName;
  ndofs = ndof;

  if (par.masterproc)
  {
    data = real2d(name.c_str(), ndofs, params.Nsteps/params.Nstat + 1);
  }
}

};



template <uint nprog, uint nconst, uint nstats> class Stats
{
public:
  std::array<Stat,nstats> stats_arr;
  MPI_Request Req [nstats];
  MPI_Status  Status[nstats];
  int ierr;
  int statsize;
  int masterproc;
  const Topology *primal_topology;
  const Topology *dual_topology;
  Geometry<1,1,1> *primal_geometry;
  Geometry<1,1,1> *dual_geometry;

  real3d TEarr, KEarr, PEarr, IEarr, PVarr, PENSarr, trimmed_density;

  void initialize(ModelParameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
  {
    this->primal_topology = &primal_topo;
    this->dual_topology = &dual_topo;
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;

    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[DENSSTAT].initialize("mass", ndensity, params, par);
    stats_arr[DENSMAXSTAT].initialize("densmax", ndensity, params, par);
    stats_arr[DENSMINSTAT].initialize("densmin", ndensity, params, par);
    stats_arr[ESTAT].initialize("energy", 4, params, par);
    stats_arr[PVSTAT].initialize("pv", 1, params, par);
    stats_arr[PESTAT].initialize("pens", 1, params, par);
    masterproc = par.masterproc;
  
    TEarr = real3d("TE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
    KEarr = real3d("KE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
    IEarr = real3d("IE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
    PEarr = real3d("PE", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);
    PVarr = real3d("PV", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);
    PENSarr = real3d("PENS", this->primal_topology->nl, this->primal_topology->n_cells_y, this->primal_topology->n_cells_x);
    trimmed_density = real3d("trimmed_density", this->dual_topology->nl, this->dual_topology->n_cells_y, this->dual_topology->n_cells_x);

  }



  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int tind)
  {


      SArray<real,1,ndensity> masslocal, massglobal;
      SArray<real,1,ndensity> densmaxlocal, densmaxglobal;
      SArray<real,1,ndensity> densminlocal, densminglobal;
      SArray<real,1,4> elocal, eglobal;
      SArray<real,1,1> pvlocal, pvglobal;
      SArray<real,1,1> pelocal, peglobal;

      elocal(0) = 0.;
      elocal(1) = 0.;
      elocal(2) = 0.;
      elocal(3) = 0.;
      eglobal(0) = 0.;
      eglobal(1) = 0.;
      eglobal(2) = 0.;
      eglobal(3) = 0.;
      pvlocal(0) = 0.;
      pvglobal(0) = 0.;
      pelocal(0) = 0.;
      peglobal(0) = 0.;
      for (int l=0;l<ndensity;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densmaxlocal(l) = 0.; densmaxglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densminlocal(l) = 0.; densminglobal(l) = 0.;}

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

      // yakl::parallel_for("ComputeDualStats", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

      parallel_for( Bounds<3>( dual_topology->nl-2, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
         real KE, PE, IE;
KE = Hk.compute_KE(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, k+1);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k+1);
TEarr(k+1, j, i) = KE + PE + IE;
KEarr(k+1, j, i) = KE;
PEarr(k+1, j, i) = PE;
IEarr(k+1, j, i) = IE;
});
parallel_for( Bounds<2>( dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
   real KE, PE, IE;
KE = Hk.compute_KE_bottom(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, 0);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, 0);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, 0);
TEarr(0, j, i) = KE + PE + IE;
KEarr(0, j, i) = KE;
PEarr(0, j, i) = PE;
IEarr(0, j, i) = IE;
KE = Hk.compute_KE_top(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, dual_topology->nl-1);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, dual_topology->nl-1);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, dual_topology->nl-1);
TEarr(dual_topology->nl-1, j, i) = KE + PE + IE;
KEarr(dual_topology->nl-1, j, i) = KE;
PEarr(dual_topology->nl-1, j, i) = PE;
IEarr(dual_topology->nl-1, j, i) = IE;
});

elocal(0) = yakl::intrinsics::sum(TEarr);
elocal(1) = yakl::intrinsics::sum(KEarr);
elocal(2) = yakl::intrinsics::sum(PEarr);
elocal(3) = yakl::intrinsics::sum(IEarr);

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

// yakl::parallel_for("ComputePrimalStats", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//  int k, j, i;
//  yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( primal_topology->nl-2, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
   pvpe vals_pvpe;
   vals_pvpe = PVPE.compute_PVPE(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j, k+1);
   PVarr(k+1, j, i) = vals_pvpe.pv;
   PENSarr(k+1, j, i) = vals_pvpe.pe;
    });
    parallel_for( Bounds<2>( primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int j, int i) { 
     pvpe vals_pvpe;
     vals_pvpe = PVPE.compute_PVPE_bottom(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j, 0);
     PVarr(0, j, i) = vals_pvpe.pv;
     PENSarr(0, j, i) = vals_pvpe.pe;
     vals_pvpe = PVPE.compute_PVPE_top(progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j, primal_topology->nl-1);
     PVarr(primal_topology->nl-1, j, i) = vals_pvpe.pv;
     PENSarr(primal_topology->nl-1, j, i) = vals_pvpe.pe;
      });
      
    pvlocal(0) = yakl::intrinsics::sum(PVarr);
    pelocal(0) = yakl::intrinsics::sum(PENSarr);

    for (int l=0;l<ndensity;l++)
    {
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        trimmed_density(k,j,i) = progvars.fields_arr[DENSVAR].data(l,k+dks,j+djs,i+dis);
      });

    masslocal(l) = yakl::intrinsics::sum(trimmed_density);
    densmaxlocal(l) = yakl::intrinsics::maxval(trimmed_density);
    densminlocal(l) = yakl::intrinsics::minval(trimmed_density);
    }

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, ndensity, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSSTAT]);
    this->ierr = MPI_Ireduce( &densmaxlocal, &densmaxglobal, ndensity, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
    this->ierr = MPI_Ireduce( &densminlocal, &densminglobal, ndensity, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
    this->ierr = MPI_Ireduce( &elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);
    this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
    for (int l=0;l<ndensity;l++)
    {
  this->stats_arr[DENSSTAT].data(l,tind) = massglobal(l);
  this->stats_arr[DENSMAXSTAT].data(l,tind) = densmaxglobal(l);
  this->stats_arr[DENSMINSTAT].data(l,tind) = densminglobal(l);
}

  this->stats_arr[ESTAT].data(0,tind) = eglobal(0);
  this->stats_arr[ESTAT].data(1,tind) = eglobal(1);
  this->stats_arr[ESTAT].data(2,tind) = eglobal(2);
  this->stats_arr[ESTAT].data(3,tind) = eglobal(3);
  this->stats_arr[PVSTAT].data(0,tind) = pvglobal(0);
  this->stats_arr[PESTAT].data(0,tind) = peglobal(0);
  }
  }
};


template <uint nvars> void set_dofs_arr(SArray<int, 2, nvars, 3> &dofs_arr, int var, int basedof, int extdof, int ndofs)
{
  dofs_arr(var, 0) = basedof;
  dofs_arr(var, 1) = extdof;
  dofs_arr(var, 2) = ndofs;
}

// *******   VariableSet Initialization   ***********//
template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int, 2, nprog, 3> &prog_ndofs_arr, SArray<int, 2, nconst, 3> &const_ndofs_arr, SArray<int, 2, naux, 3> &aux_ndofs_arr, SArray<int, 2, ndiag, 3> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{

  //primal grid represents straight quantities, dual grid twisted quantities
  // "edges" are 0-forms
// ndims is the BASEDIM size!


  // ADD Q RELATED stuff


  // v, w, dens
  prog_topo_arr[VVAR] = &ptopo;
  prog_topo_arr[WVAR] = &ptopo;
  prog_topo_arr[DENSVAR] = &dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[WVAR] = "w";
  prog_names_arr[DENSVAR] = "dens";
  set_dofs_arr(prog_ndofs_arr, VVAR, 1, 0, 1); //v = straight (1,0)-form
  set_dofs_arr(prog_ndofs_arr, WVAR, 0, 1, 1); //w = straight (0,1)-form
  set_dofs_arr(prog_ndofs_arr, DENSVAR, ndims, 1, ndensity); //dens = twisted (n,1)-form

  // hs  
  const_topo_arr[HSVAR] = &dtopo;
  const_names_arr[HSVAR] = "hs";
  set_dofs_arr(const_ndofs_arr, HSVAR, ndims, 1, 1); //hs = twisted (n,1)-form
  const_topo_arr[CORIOLISXZVAR] = &ptopo;
  const_names_arr[CORIOLISXZVAR] = "coriolisxz";
  set_dofs_arr(const_ndofs_arr, CORIOLISXZVAR, 1, 1, 1); //f = straight (1,1)-form
  
  //functional derivatives = F, B, K, he, U
  aux_topo_arr[BVAR] = &ptopo;
  aux_topo_arr[FVAR] = &dtopo;
  aux_topo_arr[UVAR] = &dtopo;
  aux_topo_arr[HEVAR] = &dtopo;
  aux_topo_arr[FWVAR] = &dtopo;
  aux_topo_arr[UWVAR] = &dtopo;
  aux_topo_arr[HEWVAR] = &dtopo;
  aux_topo_arr[KVAR] = &dtopo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  aux_names_arr[FWVAR] = "Fw";
  aux_names_arr[UWVAR] = "Uw";
  aux_names_arr[HEWVAR] = "hew";
  set_dofs_arr(aux_ndofs_arr, BVAR, 0, 0, ndensity); //B = straight (0,0)-form
  set_dofs_arr(aux_ndofs_arr, KVAR, ndims, 1, 1);   //K = twisted (n,1)-form
  set_dofs_arr(aux_ndofs_arr, FVAR, ndims-1, 1, 1);  //F = twisted (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, UVAR, ndims-1, 1, 1); //U = twisted (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, HEVAR, ndims-1, 1, 1); //he lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, FWVAR, ndims, 0, 1);  //Fw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, UWVAR, ndims, 0, 1); //Uw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, HEWVAR, ndims, 0, 1); //hew lives on vert dual edges, associated with Fw
  
  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &dtopo;
  aux_topo_arr[DENSEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENSVERTRECONVAR] = &dtopo;
  aux_topo_arr[DENSVERTEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENS0VAR] = &ptopo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSVERTRECONVAR] = "densvertrecon";
  aux_names_arr[DENSVERTEDGERECONVAR] = "densvertedgerecon";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  set_dofs_arr(aux_ndofs_arr, DENSRECONVAR, ndims-1, 1, ndensity);  //densrecon lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSEDGERECONVAR, ndims, 1, 2*ndims*ndensity); //densedgerecon lives on dual cells, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSVERTRECONVAR, ndims, 0, ndensity);  //densvertrecon lives on vert dual edges, associated with Fw
  set_dofs_arr(aux_ndofs_arr, DENSVERTEDGERECONVAR, ndims, 1, 2*ndensity); //densedgerecon lives on dual cells, associated with Fw
  set_dofs_arr(aux_ndofs_arr, DENS0VAR, 0, 0, ndensity); //dens0 = straight (0,0)-form

  //fct stuff- Phi, Mf, edgeflux
  aux_topo_arr[PHIVAR] = &dtopo;
  aux_topo_arr[PHIVERTVAR] = &dtopo;
  aux_topo_arr[MFVAR] = &dtopo;
  aux_topo_arr[EDGEFLUXVAR] = &dtopo;
  aux_topo_arr[VERTEDGEFLUXVAR] = &dtopo;
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[PHIVERTVAR] = "PhiVert";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_names_arr[VERTEDGEFLUXVAR] = "vertedgeflux";
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims-1, 1, ndensity); 
  set_dofs_arr(aux_ndofs_arr, PHIVERTVAR, ndims, 0, ndensity); 
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ndensity);  
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims-1, 1, ndensity); 
  set_dofs_arr(aux_ndofs_arr, VERTEDGEFLUXVAR, ndims, 0, ndensity); 
  
  
  //Q stuff
  aux_topo_arr[QXZ0VAR] = &dtopo;
  aux_names_arr[QXZ0VAR] = "QXZ0";
  set_dofs_arr(aux_ndofs_arr, QXZ0VAR, 0, 0, 1); //Q0 = twisted (0,0)-form
  aux_topo_arr[QXZRECONVAR] = &ptopo;
  aux_topo_arr[QXZEDGERECONVAR] = &ptopo;
  aux_topo_arr[QXZVERTRECONVAR] = &ptopo;
  aux_topo_arr[QXZVERTEDGERECONVAR] = &ptopo;
  aux_topo_arr[QXZFLUXVAR] = &ptopo;
  aux_topo_arr[QXZVERTFLUXVAR] = &ptopo;
  aux_names_arr[QXZRECONVAR] = "qxzrecon";
  aux_names_arr[QXZEDGERECONVAR] = "qxzedgerecon";
  aux_names_arr[QXZVERTRECONVAR] = "qxzvertrecon";
  aux_names_arr[QXZVERTEDGERECONVAR] = "qxzvertedgerecon";
  aux_names_arr[QXZFLUXVAR] = "qxzflux";
  aux_names_arr[QXZVERTFLUXVAR] = "qxzvertflux";
  set_dofs_arr(aux_ndofs_arr, QXZRECONVAR, 0, 1, 1);  //qxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, QXZEDGERECONVAR, ndims, 1, 2*1); //qxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(aux_ndofs_arr, QXZVERTRECONVAR, 1, 0, 1);  //qxzsvertrecon lives on horiz primal edges, associated with v
  set_dofs_arr(aux_ndofs_arr, QXZVERTEDGERECONVAR, ndims, 1, 2*1); //qxzvertedgerecon lives on primal cells, associated with F/v
  set_dofs_arr(aux_ndofs_arr, QXZFLUXVAR, 0, 1, 1); //qxzflux lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, QXZVERTFLUXVAR, 1, 0, 1); //qxzvertflux lives on horiz primal edges, associated with v

  aux_topo_arr[FTVAR] = &ptopo;
  aux_topo_arr[FTWVAR] = &ptopo;
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[FTWVAR] = "FTW";
  set_dofs_arr(aux_ndofs_arr, FTVAR, 1, 0, 1); //FT = straight (1,0)-form ie Fw at v pts
  set_dofs_arr(aux_ndofs_arr, FTWVAR, 0, 1, 1); //FTW = straight (0,1)-form ie F at w pts
  
  aux_topo_arr[FXZ0VAR] = &dtopo;
  aux_topo_arr[CORIOLISXZRECONVAR] = &ptopo;
  aux_topo_arr[CORIOLISXZEDGERECONVAR] = &ptopo; 
  aux_topo_arr[CORIOLISXZVERTRECONVAR] = &ptopo;
  aux_topo_arr[CORIOLISXZVERTEDGERECONVAR] = &ptopo; 
  aux_names_arr[FXZ0VAR] = "fxz0";
  aux_names_arr[CORIOLISXZRECONVAR] = "coriolisxzrecon";
  aux_names_arr[CORIOLISXZEDGERECONVAR] = "coriolisxzedgerecon";
  aux_names_arr[CORIOLISXZVERTRECONVAR] = "coriolisxzvertrecon";
  aux_names_arr[CORIOLISXZVERTEDGERECONVAR] = "coriolisxzvertedgerecon";
  set_dofs_arr(aux_ndofs_arr, FXZ0VAR, 0, 0, 1);  //fxz0 is a twisted 0,0 form
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZRECONVAR, 0, 1, 1);  //coriolisxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZEDGERECONVAR, ndims, 1, 2*1); //coriolisxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZVERTRECONVAR, 1, 0, 1);  //coriolisxzsvertrecon lives on horiz primal edges, associated with v
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZVERTEDGERECONVAR, ndims, 1, 2*1); //coriolisxzvertedgerecon lives on primal cells, associated with F/v  
  
  // q, concentration 0-forms for dens
  diag_topo_arr[DENSLDIAGVAR] = &ptopo;
  diag_names_arr[DENSLDIAGVAR] = "densl";
  set_dofs_arr(diag_ndofs_arr, DENSLDIAGVAR, 0, 0, ndensity); //densldiag = straight (0,0)-form
  diag_topo_arr[QXZDIAGVAR] = &dtopo;
  diag_names_arr[QXZDIAGVAR] = "QXZl";
  set_dofs_arr(diag_ndofs_arr, QXZDIAGVAR, 0, 0, 1); //Qldiag = twisted (0,0)-form

  diag_topo_arr[FDIAGVAR] = &dtopo;
  diag_topo_arr[HEDIAGVAR] = &dtopo;
  diag_topo_arr[FWDIAGVAR] = &dtopo;
  diag_topo_arr[HEWDIAGVAR] = &dtopo;
  diag_topo_arr[KDIAGVAR] = &dtopo;
  diag_names_arr[KDIAGVAR] = "K";
  diag_names_arr[FDIAGVAR] = "F";
  diag_names_arr[HEDIAGVAR] = "he";
  diag_names_arr[FWDIAGVAR] = "Fw";
  diag_names_arr[HEWDIAGVAR] = "hew";
  set_dofs_arr(diag_ndofs_arr, KDIAGVAR, ndims, 1, 1);   //K = twisted (n,1)-form
  set_dofs_arr(diag_ndofs_arr, FDIAGVAR, ndims-1, 1, 1);  //F = twisted (n-1,1)-form
  set_dofs_arr(diag_ndofs_arr, HEDIAGVAR, ndims-1, 1, 1); //he lives on horiz dual edges, associated with F
  set_dofs_arr(diag_ndofs_arr, FWDIAGVAR, ndims, 0, 1);  //Fw = twisted (n,0)-form
  set_dofs_arr(diag_ndofs_arr, HEWDIAGVAR, ndims, 0, 1); //hew lives on vert dual edges, associated with Fw
  diag_topo_arr[UDIAGVAR] = &dtopo;
  diag_topo_arr[UWDIAGVAR] = &dtopo;  
  diag_names_arr[UDIAGVAR] = "U";
  diag_names_arr[UWDIAGVAR] = "UW";
  set_dofs_arr(diag_ndofs_arr, UDIAGVAR, ndims-1, 1, 1); //U = twisted (n-1,1)-form
  set_dofs_arr(diag_ndofs_arr, UWDIAGVAR, ndims, 0, 1); //Uw = twisted (n,0)-form  
  
}

#endif
