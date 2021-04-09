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
uint constexpr nprognostic = 4;
#define VVAR 0
#define WVAR 1
#define DENSVAR 2
#define DENSFCTVAR 3

// hs, coriolis
uint constexpr nconstant = 2;
#define HSVAR 0
#define CORIOLISVAR 1

//functional derivatives = F, B, BFCT, K, he
//primal grid reconstruction stuff- U, dens0, densfct0, edgerecon, recon, edgereconfct, reconfct
//fct stuff- Phi, Mf, edgeflux

uint constexpr nauxiliary = 24;

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

#define BFCTVAR 13
#define DENSFCT0VAR 14
#define DENSFCTRECONVAR 15
#define DENSFCTEDGERECONVAR 16
#define DENSFCTVERTRECONVAR 17
#define DENSFCTVERTEDGERECONVAR 18
#define PHIVAR 19
#define VERTPHIVAR 20
#define EDGEFLUXVAR 21
#define VERTEDGEFLUXVAR 22
#define MFVAR 23

// ADD Q STUFF BACK HERE!!!
// FIGURE OUT EXACTLY WHAT WE NEED...
// #define FTVAR 8
// #define Q0VAR 9
// #define F0VAR 10
// #define QRECONVAR 11
// #define QEDGERECONVAR 12
// #define CORIOLISRECONVAR 13
// #define CORIOLISEDGERECONVAR 14

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 3;
#define QDIAGVAR 0
#define DENSLDIAGVAR 1
#define DENSFCTLDIAGVAR 2

//track total densities, dens min/max, densfct min/max, energy (total, K, P, I), PV, PE
uint constexpr nstats = 9;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define DENSFCTSTAT 4
#define DENSFCTMAXSTAT 5
#define DENSFCTMINSTAT 6
#define PVSTAT 7
#define PESTAT 8

// *******   Functionals/Hamiltonians   ***********//

//THERE HAS TO BE A BETTER WAY TO DO THIS...maybe once Hk/Hs/functionals/thermo are all one big class we can substantially reduce switches?
// ie have 1 line per "eqnset"?

// THESE KINETIC ENERGIES ARE MISSING W....

#ifdef _SWE
//std::cout << "SWE Hamiltonian" << "\n";
Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;
Hamiltonian_SWE_Hs<ntracers, ntracers_fct> Hs;
ThermoPotential thermo;
#endif

#ifdef _TSWE
//std::cout << "TSWE Hamiltonian" << "\n";
Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;
Hamiltonian_TSWE_Hs<ntracers, ntracers_fct> Hs;
ThermoPotential thermo;
#endif

#ifdef _CE
//std::cout << "CE Hamiltonian" << "\n";

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


#ifdef _AN
//std::cout << "AN Hamiltonian" << "\n";

Functional_PVPE_AN PVPE;
Hamiltonian_Hk_AN Hk;
Hamiltonian_AN_Hs Hs;

#ifdef _IDEAL_GAS_POTTEMP
//std::cout << "using ideal gas potential temperature" << "\n";
IdealGas_Pottemp thermo;
#endif
#ifdef _IDEAL_GAS_ENTROPY
//std::cout << "using ideal gas entropy" << "\n";
IdealGas_Entropy thermo;
#endif

#endif


#ifdef _MCE
//std::cout << "MCE Hamiltonian" << "\n";

#ifdef _USE_RHOD_VARIANT
//std::cout << "using rho d" << "\n";
Functional_PVPE_rhod PVPE;
Hamiltonian_Hk_rhod Hk;

#ifdef _USE_P_VARIANT
//std::cout << "using p" << "\n";
Hamiltonian_MCE_rhod_p_Hs Hs;
#else
//std::cout << "not using p" << "\n";
Hamiltonian_MCE_rhod_Hs Hs;
#endif

#else
//std::cout << "using rho" << "\n";
Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;

#ifdef _USE_P_VARIANT
//std::cout << "using p" << "\n";
Hamiltonian_MCE_rho_p_Hs Hs;
#else
//std::cout << "not using p" << "\n";
Hamiltonian_MCE_rho_Hs Hs;
#endif

#endif

#ifdef _CONST_KAPPA_VIRPOTTEMP
//std::cout << "using constant kappa potential temperature" << "\n";
ConstantKappa_VirtualPottemp thermo;
#endif
#ifdef _CONST_KAPPA_ENTROPY
//std::cout << "using constant kappa entropy" << "\n";
ConstantKappa_Entropy thermo;
#endif
#ifdef _UNAPPROX_POTTEMP
//std::cout << "using unapprox potential temperature" << "\n";
Unapprox_Pottemp thermo;
#endif
#ifdef _UNAPPROX_ENTROPY
std::cout << "using unapprox entropy" << "\n";
Unapprox_Entropy thermo;
#endif

#endif


#ifdef _MAN
//std::cout << "MAN Hamiltonian" << "\n";

Functional_PVPE_AN PVPE;
Hamiltonian_Hk_AN Hk;
Hamiltonian_MAN_Hs Hs;

#ifdef _CONST_KAPPA_VIRPOTTEMP
//std::cout << "using constant kappa potential temperature" << "\n";
ConstantKappa_VirtualPottemp thermo;
#endif
#ifdef _CONST_KAPPA_ENTROPY
//std::cout << "using constant kappa entropy" << "\n";
ConstantKappa_Entropy thermo;
#endif
#ifdef _UNAPPROX_POTTEMP
//std::cout << "using unapprox potential temperature" << "\n";
Unapprox_Pottemp thermo;
#endif
#ifdef _UNAPPROX_ENTROPY
//std::cout << "using unapprox entropy" << "\n";
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
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;

  bool is_initialized;

   Diagnostics() {
     this->is_initialized = false;
     std::cout << "CREATED DIAGNOSTICS\n";
   }

   void initialize(const Topology &ptopo, const Topology &dtopo, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom)
   {
     this->primal_topology = &ptopo;
     this->dual_topology = &dtopo;
     this->primal_geometry = &pgeom;
     this->dual_geometry = &dgeom;
     this->is_initialized = true;
   }


// ADD Q STUFF BACK HERE!!!

   void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
   {

     int dis = dual_topology->is;
     int djs = dual_topology->js;
        
    yakl::parallel_for("ComputeDiagIp", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      SArray<real,1> zeta;
      real hv;
      yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

 // ARE THESE REALLY DUAL FORMS?
 // MAYBE THEY ARE STRAIGHT 0-FORMS?

 //THIS RELIES ON DENSITY BEING IN DENSVAR0
 //IS THERE ANOTHER/BETTER WAY?
 //OTHER DIAGNOSTICS?

 //THIS CAN REALLY BE DONE OFFLINE
 //MAYBE DO HODGE STARS HERE?
 // compute dens0var
 for(int l=0; l<ndensity; l++)
 {diagnostic_vars.fields_arr[DENSLDIAGVAR].data(l, k, j+djs, i+dis) = x.fields_arr[DENSVAR].data(l, k, j+djs, i+dis) / x.fields_arr[DENSVAR].data(0, k, j+djs, i+dis);}

 // compute densfct0var
 for(int l=0; l<ndensityfct; l++)
 {diagnostic_vars.fields_arr[DENSFCTLDIAGVAR].data(l, k, j+djs, i+dis) = x.fields_arr[DENSFCTVAR].data(l, k, j+djs, i+dis) / x.fields_arr[DENSVAR].data(0, k, j+djs, i+dis);}
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
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;

  TransformMatrices<real> trans;

  SArray<real,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dual_wenoRecon;
  SArray<real,(dual_reconstruction_order-1)/2+2> dual_wenoIdl;
  real dual_wenoSigma;

  bool is_initialized;
  
   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

   void initialize(ModelParameters &params, const Topology &primal_topo, const Topology &dual_topo, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom, ExchangeSet<naux> &aux_exchange, ExchangeSet<nconst> &const_exchange)
   {
     this->primal_topology = &primal_topo;
     this->dual_topology = &dual_topo;
     this->primal_geometry = &primal_geom;
     this->dual_geometry = &dual_geom;
     this->aux_exchange = &aux_exchange;
     this->const_exchange = &const_exchange;

    trans.coefs_to_gll_lower( dual_to_gll );
    trans.weno_sten_to_coefs(dual_wenoRecon);
    wenoSetIdealSigma<dual_reconstruction_order>(dual_wenoIdl,dual_wenoSigma);

    Hk.initialize(params, *this->primal_geometry, *this->dual_geometry);
    Hs.initialize(params, thermo, *this->primal_geometry, *this->dual_geometry);
    
    this->is_initialized = true;
  }


  void compute_constants(VariableSet<nconst> &const_vars, VariableSet<nprog> &x)
  {}
 
   void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
    realArr Uvar, realArr dens0var, realArr densfct0var,
    const realArr Vvar, const realArr densvar, const realArr densfctvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;

int dis = dual_topology->is;
int djs = dual_topology->js;

      yakl::parallel_for("ComputeDiagIp", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        real hv;
        yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

        // compute dens0var = I densvar, densfct0var = I densfctvar
        compute_I<ndensity, diff_ord>(dens0var, densvar, *this->primal_geometry, *this->dual_geometry, pis, pjs, 0, i, j, k);
        compute_I<ndensityfct, diff_ord>(densfct0var, densfctvar, *this->primal_geometry, *this->dual_geometry, pis, pjs, 0, i, j, k);

});

yakl::parallel_for("ComputeDiagId", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  real hv;
  yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

  // compute U = H v, q0, f0
  compute_H<1, diff_ord>(Uvar, Vvar, *this->primal_geometry, *this->dual_geometry, dis, djs, 0, i, j, k);

  // ADD Q/W RELATED stuff

      });

    }

    void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
      realArr Fvar, realArr Kvar, realArr HEvar, const realArr Vvar, const realArr Uvar, const realArr dens0var, const realArr densfct0var) {


        int dis = dual_topology->is;
        int djs = dual_topology->js;

        yakl::parallel_for("ComputeDiagII", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
        Hk.compute_dKdv(Fvar, Kvar, HEvar, Vvar, Uvar, dens0var, densfct0var, dis, djs, 0, i, j, k);

// ADD Q/W RELATED stuff


      });

      }


  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    realArr Bvar, realArr Bfctvar,
    const realArr Fvar, const realArr Uvar,
    const realArr Kvar, const realArr dens0var, const realArr densfct0var, const realArr HSvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;

      yakl::parallel_for("ComputeDiagIII", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

Hs.compute_dHsdx(Bvar, Bfctvar, dens0var, densfct0var, HSvar, pis, pjs, 0, i, j, k);
Hk.compute_dKddens(Bvar, Bfctvar, Kvar, pis, pjs, 0, i, j, k);

// ADD Q/W RELATED stuff

    });

    }




  void YAKL_INLINE compute_edge_reconstructions(realArr densedgereconvar, realArr densfctedgereconvar,
    const realArr dens0var, const realArr densfct0var) {


int dis = dual_topology->is;
int djs = dual_topology->js;

    yakl::parallel_for("ComputeDualEdgeRecon", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

// ADD Q/W RELATED stuff

      compute_twisted_edge_recon<ndensity, dual_reconstruction_type, dual_reconstruction_order>(densedgereconvar, dens0var, dis, djs, 0, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
      compute_twisted_edge_recon<ndensityfct, dual_reconstruction_type, dual_reconstruction_order>(densfctedgereconvar, densfct0var, dis, djs, 0, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
    });



  }

  void YAKL_INLINE compute_recons(
  realArr densreconvar, realArr densfctreconvar,
  const realArr densedgereconvar, const realArr densfctedgereconvar, const realArr HEvar,
  const realArr Uvar) {

int dis = dual_topology->is;
int djs = dual_topology->js;

    yakl::parallel_for("ComputeDualRecon", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

      // ADD Q/W RELATED stuff

      compute_twisted_recon<ndensity, dual_reconstruction_type>(densreconvar, densedgereconvar, Uvar, dis, djs, 0, i, j, k);
      compute_twisted_recon<ndensityfct, dual_reconstruction_type>(densfctreconvar, densfctedgereconvar, Uvar, dis, djs, 0, i, j, k);

    //scale twisted recons
    for (int d=0;d<ndims;d++) {
    for (int l=0;l<ndensity;l++) {
    densreconvar(l+d*ndensity,k,j+djs,i+dis) = densreconvar(l+d*ndensity,k,j+djs,i+dis) / HEvar(d,k,j+djs,i+dis);
  }}
  for (int d=0;d<ndims;d++) {
  for (int l=0;l<ndensityfct;l++) {
  densfctreconvar(l+d*ndensityfct,k,j+djs,i+dis) = densfctreconvar(l+d*ndensityfct,k,j+djs,i+dis) / HEvar(d,k,j+djs,i+dis);
}}
    });

}


  void YAKL_INLINE compute_tendencies(
  realArr denstendvar, realArr densfcttendvar, realArr Vtendvar,
  const realArr densreconvar, const realArr densfctreconvar,
  const realArr Bvar, const realArr Bfctvar, const realArr Fvar, const realArr Phivar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;

int dis = dual_topology->is;
int djs = dual_topology->js;

// ADD Q/W RELATED stuff

      yakl::parallel_for("ComputePrimalTendencies", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

    compute_wD1<ndensity> (Vtendvar, densreconvar, Bvar, pis, pjs, 0, i, j, k);
    compute_wD1_fct<ndensityfct, ADD_MODE::ADD> (Vtendvar, densfctreconvar, Phivar, Bfctvar, pis, pjs, 0, i, j, k);

});

// ADD Q/W RELATED stuff

yakl::parallel_for("ComputeDualTendencies", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

  compute_wDbar2<ndensity> (denstendvar, densreconvar, Fvar, dis, djs, 0, i, j, k);
  compute_wDbar2_fct<ndensityfct> (densfcttendvar, densfctreconvar, Phivar, Fvar, dis, djs, 0, i, j, k);

  });

  }




  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

      //Compute U, q0, hf, dens0, densfct0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data);

      this->aux_exchange->exchanges_arr[DENSFCT0VAR].exchange_field(auxiliary_vars.fields_arr[DENSFCT0VAR]);
      this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);


      //Compute K, F, he
      compute_functional_derivatives_and_diagnostic_quantities_II(
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      x.fields_arr[VVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data);

      this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
      this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
      this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);

      //Compute FT, B, Bfct
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data, const_vars.fields_arr[HSVAR].data);

      this->aux_exchange->exchanges_arr[BFCTVAR].exchange_field(auxiliary_vars.fields_arr[BFCTVAR]);
      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);

      // Compute densrecon, densfctrecon, qrecon and frecon
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data);

      this->aux_exchange->exchanges_arr[DENSFCTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);

      compute_recons(
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data,
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[HEVAR].data,
      auxiliary_vars.fields_arr[UVAR].data);
      
      this->aux_exchange->exchanges_arr[DENSFCTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTRECONVAR]);
      this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);


//Compute fct

int dis = dual_topology->is;
int djs = dual_topology->js;

yakl::parallel_for("ComputeEdgeFlux", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
compute_edgefluxes<ndensityfct> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, dis, djs, 0, i, j, k);
});
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


yakl::parallel_for("ComputeMf", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
compute_Mf<ndensityfct> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, dis, djs, 0, i, j, k);
});

this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

yakl::parallel_for("ComputePhi", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
compute_Phi<ndensityfct> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSFCTVAR].data, dis, djs, 0, i, j, k);
});

this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);

      // Compute tendencies
      compute_tendencies(
      xtend.fields_arr[DENSVAR].data, xtend.fields_arr[DENSFCTVAR].data, xtend.fields_arr[VVAR].data,
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data,
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[PHIVAR].data);
}

};



// *******   Statistics Calculations   ***********//

// THIS STUFF SHOULD BE CLEANED UP AND GENERALIZED LIKE VARIABLE SETS IF POSSIBLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE!
class Stat
{
public:
  realArr data;
  std::string name;
  int ndofs;


void initialize(std::string statName, int ndof, ModelParameters &params, Parallel &par)
{
  name = statName;
  ndofs = ndof;

  if (par.masterproc)
  {
    data = realArr(name.c_str(), ndofs, params.Nsteps/params.Nstat + 1);
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
  Geometry<ndims,1,1,1> *primal_geometry;
  Geometry<ndims,1,1,1> *dual_geometry;
  
  void initialize(ModelParameters &params, Parallel &par, const Topology &primal_topo, const Topology &dual_topo, Geometry<ndims,1,1,1> &primal_geom, Geometry<ndims,1,1,1> &dual_geom)
  {
    this->primal_topology = &primal_topo;
    this->dual_topology = &dual_topo;
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;

    // ADD Q/W RELATED stuff

    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[DENSSTAT].initialize("mass", ndensity, params, par);
    stats_arr[DENSMAXSTAT].initialize("densmax", ndensity, params, par);
    stats_arr[DENSMINSTAT].initialize("densmin", ndensity, params, par);
    stats_arr[DENSFCTSTAT].initialize("massfct", ndensityfct, params, par);
    stats_arr[DENSFCTMAXSTAT].initialize("densfctmax", ndensityfct, params, par);
    stats_arr[DENSFCTMINSTAT].initialize("densfctmin", ndensityfct, params, par);
    stats_arr[ESTAT].initialize("energy", 4, params, par);
    stats_arr[PVSTAT].initialize("pv", 1, params, par);
    stats_arr[PESTAT].initialize("pens", 1, params, par);
    masterproc = par.masterproc;
    
  }



  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int i)
  {


      SArray<real,ndensity> masslocal, massglobal;
      SArray<real,ndensity> densmaxlocal, densmaxglobal;
      SArray<real,ndensity> densminlocal, densminglobal;
      SArray<real,ndensityfct> massfctlocal, massfctglobal;
      SArray<real,ndensityfct> densfctmaxlocal, densfctmaxglobal;
      SArray<real,ndensityfct> densfctminlocal, densfctminglobal;
      SArray<real,4> elocal, eglobal;

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
      for (int l=0;l<ndensityfct;l++) {massfctlocal(l) = 0.; massfctglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {densfctmaxlocal(l) = 0.; densfctmaxglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {densfctminlocal(l) = 0.; densfctminglobal(l) = 0.;}

int dis = dual_topology->is;
int djs = dual_topology->js;

      yakl::parallel_for("ComputeDualStats", dual_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
         real KE, PE, IE;

KE = Hk.compute_KE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, dis, djs, 0, i, j, k);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, 0, i, j, k);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, dis, djs, 0, i, j, k);

elocal(3) += IE;
elocal(2) += PE;
elocal(1) += KE;
elocal(0) += KE + PE + IE;

});

// FIX THIS

// int pis = primal_topology->is;
// int pjs = primal_topology->js;
// 
// yakl::parallel_for("ComputePrimalStats", primal_topology->n_cells_layers, YAKL_LAMBDA (int iGlob) {
//   int k, j, i;
//   yakl::unpackIndices(iGlob, primal_topology->nl, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
// 
//    pvpe vals_pvpe;
//    vals_pvpe = PVPE.compute_PVPE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, constvars.fields_arr[CORIOLISVAR].data, pis, pjs, 0, i, j, k);
//    pvlocal(0) += vals_pvpe.pv;
//     pelocal(0) += vals_pvpe.pe;
// 
//     });
    
    
    for (int l=0;l<ndensity;l++)
    {
    masslocal(l) = progvars.fields_arr[DENSVAR].sum(l);
    densmaxlocal(l) = progvars.fields_arr[DENSVAR].max(l);
    densminlocal(l) = progvars.fields_arr[DENSVAR].min(l);
    }

    for (int l=0;l<ndensityfct;l++)
    {
      massfctlocal(l) = progvars.fields_arr[DENSFCTVAR].sum(l);
      densfctmaxlocal(l) = progvars.fields_arr[DENSFCTVAR].max(l);
      densfctminlocal(l) = progvars.fields_arr[DENSFCTVAR].min(l);
    }

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, ndensity, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSSTAT]);
    this->ierr = MPI_Ireduce( &densmaxlocal, &densmaxglobal, ndensity, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
    this->ierr = MPI_Ireduce( &densminlocal, &densminglobal, ndensity, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
    this->ierr = MPI_Ireduce( &massfctlocal, &massfctglobal, ndensityfct, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSFCTSTAT]);
    this->ierr = MPI_Ireduce( &densfctmaxlocal, &densfctmaxglobal, ndensityfct, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSFCTMAXSTAT]);
    this->ierr = MPI_Ireduce( &densfctminlocal, &densfctminglobal, ndensityfct, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSFCTMINSTAT]);
    this->ierr = MPI_Ireduce( &elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);
    this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
    for (int l=0;l<ndensity;l++)
    {
  this->stats_arr[DENSSTAT].data(l,i) = massglobal(l);
  this->stats_arr[DENSMAXSTAT].data(l,i) = densmaxglobal(l);
  this->stats_arr[DENSMINSTAT].data(l,i) = densminglobal(l);
}

  for (int l=0;l<ndensityfct;l++)
  {
  this->stats_arr[DENSFCTSTAT].data(l,i) = massfctglobal(l);
  this->stats_arr[DENSFCTMAXSTAT].data(l,i) = densfctmaxglobal(l);
  this->stats_arr[DENSFCTMINSTAT].data(l,i) = densfctminglobal(l);
}

  this->stats_arr[ESTAT].data(0,i) = eglobal(0);
  this->stats_arr[ESTAT].data(1,i) = eglobal(1);
  this->stats_arr[ESTAT].data(2,i) = eglobal(2);
  this->stats_arr[ESTAT].data(3,i) = eglobal(3);
  this->stats_arr[PVSTAT].data(0,i) = pvglobal(0);
  this->stats_arr[PESTAT].data(0,i) = peglobal(0);
  }
  }
};


template <uint nvars> void set_dofs_arr(SArray<int, nvars, 3> &dofs_arr, int var, int basedof, int extdof, int ndofs)
{
  dofs_arr(var, 0) = basedof;
  dofs_arr(var, 1) = extdof;
  dofs_arr(var, 2) = ndofs;
}

// *******   VariableSet Initialization   ***********//
template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int, nprog, 3> &prog_ndofs_arr, SArray<int, nconst, 3> &const_ndofs_arr, SArray<int, naux, 3> &aux_ndofs_arr, SArray<int, ndiag, 3> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{

  //primal grid represents straight quantities, dual grid twisted quantities
  // "edges" are 0-forms
// ndims is the BASEDIM size!


  // ADD Q RELATED stuff


  // v, dens
// THIS IS DISTINCT FROM LAYERMODEL VVAR
// WHY?
// LAYERMODEL DOESN'T REALLY ASSUME ANY STAGGERING IN THE VERTICAL!
// IT HAS COLLOCATED PRIMAL AND DUAL GRID IN THE VERTICAL!  
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

// CAREFUL HERE WITH CORIOLIS- THIS IS IN THE "VERTICAL SLICE"
// HOW DO WE GENERALIZE TO 3D?
// we should really have qxy, qyz, qxz + associated coriolis components...
  const_topo_arr[HSVAR] = &dtopo;
  const_topo_arr[CORIOLISVAR] = &ptopo;
  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";
  set_dofs_arr(const_ndofs_arr, HSVAR, ndims, 1, 1); //hs = twisted (n,1)-form
  set_dofs_arr(const_ndofs_arr, CORIOLISVAR, 1, 1, 1); //f = straight (1,1)-form
  
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

  // q, concentration 0-forms for dens
  diag_topo_arr[QDIAGVAR] = &dtopo;
  diag_topo_arr[DENSLDIAGVAR] = &ptopo;
  diag_names_arr[QDIAGVAR] = "q";
  diag_names_arr[DENSLDIAGVAR] = "densl";
  set_dofs_arr(diag_ndofs_arr, QDIAGVAR, 0, 0, 1);  //qdiag = twisted (0,0)-form
  set_dofs_arr(diag_ndofs_arr, DENSLDIAGVAR, 0, 0, ndensity); //densldiag = straight (0,0)-form
  
  //densfct stuff- densfct, BFCT, densfct0, edgereconfct, reconfct, Phi, Mf, edgeflux, concentration 0-forms for densfct
  prog_topo_arr[DENSFCTVAR] = &dtopo;
  prog_names_arr[DENSFCTVAR] = "densfct";
  set_dofs_arr(prog_ndofs_arr, DENSFCTVAR, ndims, 1, ndensityfct);  //densfct = twisted (n,1)-form
  
  aux_topo_arr[BFCTVAR] = &ptopo;
  aux_topo_arr[DENSFCTRECONVAR] = &dtopo;
  aux_topo_arr[DENSFCTEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENSFCTVERTRECONVAR] = &dtopo;
  aux_topo_arr[DENSFCTVERTEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENSFCT0VAR] = &ptopo;
  aux_topo_arr[PHIVAR] = &dtopo;
  aux_topo_arr[VERTPHIVAR] = &dtopo;
  aux_topo_arr[MFVAR] = &dtopo;
  aux_topo_arr[EDGEFLUXVAR] = &dtopo;
  aux_topo_arr[VERTEDGEFLUXVAR] = &dtopo;
  aux_names_arr[BFCTVAR] = "Bfct";
  aux_names_arr[DENSFCT0VAR] = "densfct0";
  aux_names_arr[DENSFCTRECONVAR] = "densfctrecon";
  aux_names_arr[DENSFCTEDGERECONVAR] = "densfctedgerecon";
  aux_names_arr[DENSFCTVERTRECONVAR] = "densfctvertrecon";
  aux_names_arr[DENSFCTVERTEDGERECONVAR] = "densfctvertedgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[VERTPHIVAR] = "VertPhi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_names_arr[VERTEDGEFLUXVAR] = "vertedgeflux";
  set_dofs_arr(aux_ndofs_arr, BFCTVAR, 0, 0, ndensityfct);  //Bfct = straight (0,0)-form
  set_dofs_arr(aux_ndofs_arr, DENSFCTRECONVAR, ndims-1, 1, ndensityfct);  //densfctrecon lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, DENSFCTEDGERECONVAR, ndims, 1, 2*ndims*ndensityfct);  //densfctedgerecon lives on dual cells
  set_dofs_arr(aux_ndofs_arr, DENSFCTVERTRECONVAR, ndims, 0, ndensityfct);  //densfctvertrecon lives on vert dual edges, associated with Fw
  set_dofs_arr(aux_ndofs_arr, DENSFCTVERTEDGERECONVAR, ndims, 1, 2*ndensityfct);  //densfctvertedgerecon lives on dual cells
  set_dofs_arr(aux_ndofs_arr, DENSFCT0VAR, 0, 0, ndensityfct);  //densfct0 = straight (0,0)-form
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims-1, 1, ndensityfct); 
  set_dofs_arr(aux_ndofs_arr, VERTPHIVAR, ndims, 0, ndensityfct); 
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ndensityfct);  
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims-1, 1, ndensityfct); 
  set_dofs_arr(aux_ndofs_arr, VERTEDGEFLUXVAR, ndims, 0, ndensityfct); 


  diag_topo_arr[DENSFCTLDIAGVAR] = &ptopo;
  diag_names_arr[DENSFCTLDIAGVAR] = "densfctl";
  set_dofs_arr(diag_ndofs_arr, DENSFCTLDIAGVAR, 0, 0, ndensityfct);  //densfctldiag = straight (0,0)-form
  
}

#endif
