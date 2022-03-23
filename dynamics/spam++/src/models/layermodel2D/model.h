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
// v, dens
uint constexpr nprognostic = 2;
#define VVAR 0
#define DENSVAR 1

// hs, coriolis
uint constexpr nconstant = 2;
#define HSVAR 0
#define CORIOLISVAR 1

//functional derivatives = F, B, K, he
//dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon, coriolisedgercon, coriolisrecon
//primal grid reconstruction stuff- U, dens0, edgerecon, recon
//fct stuff- Phi, Mf, edgeflux

#if defined _AN || defined _MAN
uint constexpr nauxiliary = 19;
#else
uint constexpr nauxiliary = 18;
#endif

#define FVAR 0
#define BVAR 1
#define KVAR 2
#define HEVAR 3
#define UVAR 4

#define DENS0VAR 5
#define DENSRECONVAR 6
#define DENSEDGERECONVAR 7

#define FTVAR 8
#define Q0VAR 9
#define F0VAR 10
#define QRECONVAR 11
#define QEDGERECONVAR 12
#define CORIOLISRECONVAR 13
#define CORIOLISEDGERECONVAR 14

//RE-ARRANGE THIS STUFF?
#define PHIVAR 15
#define EDGEFLUXVAR 16
#define MFVAR 17

#if defined _AN || defined _MAN
#define PVAR 18
#endif

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 2;
#define QDIAGVAR 0
#define DENSLDIAGVAR 1

//track total densities, dens min/max, densfct min/max, energy (total, K, P, I), PV, PE,
uint constexpr nstats = 6;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define PVSTAT 4
#define PESTAT 5


// *******   Functionals/Hamiltonians   ***********//

//THERE HAS TO BE A BETTER WAY TO DO THIS...maybe once Hk/Hs/functionals/thermo are all one big class we can substantially reduce switches?
// ie have 1 line per "eqnset"?

Functional_PVPE PVPE;
Hamiltonian_Hk Hk;


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
Functional_PVPE_AN PVPE;
Hamiltonian_Hk_AN Hk;
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
Functional_PVPE_rhod PVPE;
Hamiltonian_Hk_rhod Hk;

#ifdef _USE_P_VARIANT
Hamiltonian_MCE_rhod_p_Hs Hs;
#else
Hamiltonian_MCE_rhod_Hs Hs;
#endif

#else
Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;

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

Functional_PVPE_AN PVPE;
Hamiltonian_Hk_AN Hk;
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


// JUST A NEW VARIABLE SET!
// IE EVERYTHING FCT!
// SUPER SIMPLE!
// THESE NEED FIXING- BASICALLY STUFF JUST MOVED FROM DENS TO DENSFCT...
#ifdef _CEFCT

Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;

#ifdef _USE_P_VARIANT
//std::cout << "using p" << "\n";
Hamiltonian_CE_p_Hs Hs;
#else
//std::cout << "not using p" << "\n";
Hamiltonian_CE_Hs Hs;
#endif

#ifdef _IDEAL_GAS_POTTEMP
//std::cout << "using ideal gas potential temperature" << "\n";
IdealGas_Pottemp thermo;
#endif
#ifdef _IDEAL_GAS_ENTROPY
//std::cout << "using ideal gas entropy" << "\n";
IdealGas_Entropy thermo;
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

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE

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
     
     // yakl::parallel_for("ComputeDens0", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
     //   int k, j, i;
     //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
       parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
     compute_I<ndensity, diff_ord>(diagnostic_vars.fields_arr[DENSLDIAGVAR].data, x.fields_arr[DENSVAR].data, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
     });

    // yakl::parallel_for("ComputeQ0VAR", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
       parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
     PVPE.compute_q0(diagnostic_vars.fields_arr[QDIAGVAR].data, x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k);
     });
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

    PVPE.initialize(params);
    Hk.initialize(params, *this->primal_geometry, *this->dual_geometry);
    Hs.initialize(params, thermo, *this->primal_geometry, *this->dual_geometry);
    
    this->is_initialized = true;
  }


  void compute_constants(VariableSet<nconst> &const_vars, VariableSet<nprog> &x)
  {}
 
   void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
    real4d Uvar, real4d Q0var, real4d f0var, real4d dens0var,
    const real4d Vvar, const real4d densvar, const real4d coriolisvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

      // yakl::parallel_for("ComputeDiagIp", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
        parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        // compute dens0var = I densvar
        compute_I<ndensity, diff_ord>(dens0var, densvar, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
});

// yakl::parallel_for("ComputeDiagId", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

  // compute U = H v, q0, f0
  compute_H<1, diff_ord>(Uvar, Vvar, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);
  PVPE.compute_q0f0(Q0var, f0var, Vvar, densvar, coriolisvar, dis, djs, dks, i, j, k);
      });
    }

    void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
      real4d Fvar, real4d Kvar, real4d HEvar, const real4d Vvar, const real4d Uvar, const real4d dens0var) {

        int dis = dual_topology->is;
        int djs = dual_topology->js;
        int dks = dual_topology->ks;

        // yakl::parallel_for("ComputeDiagII", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
        //   int k, j, i;
        //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
          parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
        Hk.compute_dKdv(Fvar, Kvar, HEvar, Vvar, Uvar, dens0var, dis, djs, dks, i, j, k);
      });

      }


  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    real4d FTvar, real4d Bvar,
    const real4d Fvar, const real4d Uvar,
    const real4d Kvar, const real4d dens0var, const real4d HSvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

      // yakl::parallel_for("ComputeDiagIII", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
        parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_W(FTvar, Fvar, pis, pjs, pks, i, j, k);
Hs.compute_dHsdx(Bvar, dens0var, HSvar, pis, pjs, pks, i, j, k);
Hk.compute_dKddens(Bvar, Kvar, pis, pjs, pks, i, j, k);
    });

    }




  void YAKL_INLINE compute_edge_reconstructions(real4d densedgereconvar, real4d Qedgereconvar, real4d fedgereconvar,
    const real4d dens0var, const real4d Q0var, const real4d f0var) {

      int pis = primal_topology->is;
      int pjs = primal_topology->js;
      int pks = primal_topology->ks;

      int dis = dual_topology->is;
      int djs = dual_topology->js;
      int dks = dual_topology->ks;

    // yakl::parallel_for("ComputePrimalEdgeRecon", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
      compute_straight_edge_recon<1, reconstruction_type, reconstruction_order>(Qedgereconvar, Q0var, pis, pjs, pks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
      compute_straight_edge_recon<1, coriolis_reconstruction_type, coriolis_reconstruction_order>(fedgereconvar, f0var, pis, pjs, pks, i, j, k, coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl, coriolis_wenoSigma);
    });

    // yakl::parallel_for("ComputeDualEdgeRecon", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
      compute_twisted_edge_recon<ndensity, dual_reconstruction_type, dual_reconstruction_order>(densedgereconvar, dens0var, dis, djs, dks, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
    });



  }

  void YAKL_INLINE compute_recons(
  real4d densreconvar, real4d Qreconvar, real4d Coriolisreconvar,
  const real4d densedgereconvar, const real4d Qedgereconvar, const real4d fedgereconvar, const real4d HEvar,
  const real4d FTvar, const real4d Uvar) {

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

    // yakl::parallel_for("ComputePrimalRecon", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
      compute_straight_recon<1, reconstruction_type>(Qreconvar, Qedgereconvar, FTvar, pis, pjs, pks, i, j, k);
      compute_straight_recon<1, coriolis_reconstruction_type>(Coriolisreconvar, fedgereconvar, FTvar, pis, pjs, pks, i, j, k);
    });

    // yakl::parallel_for("ComputeDualRecon", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
    //   int k, j, i;
    //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
      compute_twisted_recon<ndensity, dual_reconstruction_type>(densreconvar, densedgereconvar, Uvar, dis, djs, dks, i, j, k);
    //scale primal recons
    for (int d=0;d<ndims;d++) {
    for (int l=0;l<ndensity;l++) {
    densreconvar(l+d*ndensity,k+dks,j+djs,i+dis) = densreconvar(l+d*ndensity,k+dks,j+djs,i+dis) / HEvar(d,k+dks,j+djs,i+dis);
  }}
    });

}


  void YAKL_INLINE compute_tendencies(
  real4d denstendvar, real4d Vtendvar,
  const real4d densreconvar, const real4d Qreconvar, const real4d Coriolisreconvar,
  const real4d Bvar, const real4d Fvar, const real4d Phivar) {

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

      // yakl::parallel_for("ComputePrimalTendencies", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
        parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

    compute_wD1_fct<ndensity> (Vtendvar, densreconvar, Phivar, Bvar, pis, pjs, pks, i, j, k);
    if (qf_choice == QF_MODE::EC)
    { compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, pis, pjs, pks, i, j, k);}
    if (qf_choice == QF_MODE::NOEC)
    { compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, pis, pjs, pks, i, j, k);}
    compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Coriolisreconvar, Fvar, pis, pjs, pks, i, j, k);
});


// yakl::parallel_for("ComputeDualTendencies", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
  compute_wDbar2_fct<ndensity> (denstendvar, densreconvar, Phivar, Fvar, dis, djs, dks, i, j, k);
  });

  }




  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

      //Compute U, q0, hf, dens0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, 
      x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISVAR].data);

      this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);
      this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(auxiliary_vars.fields_arr[Q0VAR]);
      this->aux_exchange->exchanges_arr[F0VAR].exchange_field(auxiliary_vars.fields_arr[F0VAR]);


      //Compute K, F, he
      compute_functional_derivatives_and_diagnostic_quantities_II(
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      x.fields_arr[VVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data);

      this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
      this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
      this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);

      //Compute FT, B
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[BVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, const_vars.fields_arr[HSVAR].data);

      this->aux_exchange->exchanges_arr[FTVAR].exchange_field(auxiliary_vars.fields_arr[FTVAR]);
      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);

      // Compute densrecon, qrecon and frecon
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data);

      this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR]);

      compute_recons(
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
      auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[UVAR].data);
      
      this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);
      this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISRECONVAR]);


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
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);

// yakl::parallel_for("ComputeMf", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
// int k, j, i;
// yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 

compute_Mf<ndensity> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, dis, djs, dks, i, j, k);
});

this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

// yakl::parallel_for("ComputePhi", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
// int k, j, i;
// yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
compute_Phi<ndensity> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k);
});


//Don't do FCT for non-FCT vars
    for (int l=0; l<ndensity_nofct; l++)
    {
    auxiliary_vars.fields_arr[PHIVAR].set(l, 1.0);
    }

    this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);

      // Compute tendencies
      compute_tendencies(
      xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[PHIVAR].data);
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
      SArray<real,1,1> pvlocal, pvglobal;
      SArray<real,1,4> elocal, eglobal;
      SArray<real,1,1> pelocal, peglobal;

      pvlocal(0) = 0.;
      pvglobal(0) = 0.;
      pelocal(0) = 0.;
      peglobal(0) = 0.;
      elocal(0) = 0.;
      elocal(1) = 0.;
      elocal(2) = 0.;
      elocal(3) = 0.;
      eglobal(0) = 0.;
      eglobal(1) = 0.;
      eglobal(2) = 0.;
      eglobal(3) = 0.;
      for (int l=0;l<ndensity;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densmaxlocal(l) = 0.; densmaxglobal(l) = 0.;}
      for (int l=0;l<ndensity;l++) {densminlocal(l) = 0.; densminglobal(l) = 0.;}

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

      // yakl::parallel_for("ComputeDualStats", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
      parallel_for( Bounds<3>( dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
         real KE, PE, IE;
KE = Hk.compute_KE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, k);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, k);
TEarr(k, j, i) = KE + PE + IE;
KEarr(k, j, i) = KE;
PEarr(k, j, i) = PE;
IEarr(k, j, i) = IE;
});

elocal(0) = yakl::intrinsics::sum(TEarr);
elocal(1) = yakl::intrinsics::sum(KEarr);
elocal(2) = yakl::intrinsics::sum(PEarr);
elocal(3) = yakl::intrinsics::sum(IEarr);

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

// yakl::parallel_for("ComputePrimalStats", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  parallel_for( Bounds<3>( primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x) , YAKL_LAMBDA(int k, int j, int i) { 
   pvpe vals_pvpe;
   vals_pvpe = PVPE.compute_PVPE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, constvars.fields_arr[CORIOLISVAR].data, pis, pjs, 0, i, j, k);
   PVarr(k, j, i) = vals_pvpe.pv;
   PENSarr(k, j, i) = vals_pvpe.pe;
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
    this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Ireduce( &elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);

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



template <uint nvars> void set_dofs_arr(SArray<int,2, nvars, 3> &dofs_arr, int var, int basedof, int extdof, int ndofs)
{
  dofs_arr(var, 0) = basedof;
  dofs_arr(var, 1) = extdof;
  dofs_arr(var, 2) = ndofs;
}

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

#endif
