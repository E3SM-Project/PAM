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
#include "string.h"
#include "hamiltonian.h"

// Number of variables
// v, dens, densfct
uint constexpr nprognostic = 3;
#define VVAR 0
#define DENSVAR 1
#define DENSFCTVAR 2

// hs, coriolis
uint constexpr nconstant = 2;
#define HSVAR 0
#define CORIOLISVAR 1

//functional derivatives = F, B, BFCT, K, he
//dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon, coriolisedgercon, coriolisrecon
//primal grid reconstruction stuff- U, dens0, densfct0, edgerecon, recon, edgereconfct, reconfct
//fct stuff- Phi, Mf, edgeflux

uint constexpr nauxiliary = 22;

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
#define BFCTVAR 15
#define DENSFCT0VAR 16
#define DENSFCTRECONVAR 17
#define DENSFCTEDGERECONVAR 18
#define PHIVAR 19
#define EDGEFLUXVAR 20
#define MFVAR 21

// q, associated concentration 0-forms for den

uint constexpr ndiagnostic = 3;
#define QDIAGVAR 0
#define DENSLDIAGVAR 1
#define DENSFCTLDIAGVAR 2

//track total densities, dens min/max, densfct min/max, energy (total, K, P, I), PV, PE,
uint constexpr nstats = 9;

#define DENSSTAT 0
#define DENSMINSTAT 1
#define DENSMAXSTAT 2
#define ESTAT 3
#define PVSTAT 4
#define PESTAT 5
#define DENSFCTSTAT 6
#define DENSFCTMAXSTAT 7
#define DENSFCTMINSTAT 8


// *******   Functionals/Hamiltonians   ***********//


#ifdef _SWE
Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;
Hamiltonian_SWE_Hs<ntracers, ntracers_fct> Hs;
ThermoPotential thermo;
#endif

#ifdef _TSWE
Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;
Hamiltonian_TSWE_Hs<ntracers, ntracers_fct> Hs;
ThermoPotential thermo;
#endif

#ifdef _CE

Functional_PVPE_rho PVPE;
Hamiltonian_Hk_rho Hk;

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


   void YAKL_INLINE compute_diagnostic_quantities(
     realArr Q0var, realArr dens0var, realArr densfct0var,
     const realArr Vvar, const realArr densvar, const realArr densfctvar, const realArr coriolisvar) {

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;
       
   yakl::parallel_for("ComputeDiagIp", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
     int k, j, i;
     SArray<real,1> zeta;
     real hv;
     yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

       PVPE.compute_q0(Q0var, Vvar, densvar, densfctvar, dis, djs, dks, i, j, k);

// ARE THESE REALLY DUAL FORMS?
// MAYBE THEY ARE STRAIGHT 0-FORMS?

//THIS RELIES ON DENSITY BEING IN DENSVAR0
//IS THERE ANOTHER/BETTER WAY?
//OTHER DIAGNOSTICS?

//THIS CAN REALLY BE DONE OFFLINE
//MAYBE DO HODGE STARS HERE?
// compute dens0var
for(int l=0; l<ndensity; l++)
{dens0var(l, k+dks, j+djs, i+dis) = densvar(l, k+dks, j+djs, i+dis) / densvar(0, k+dks, j+djs, i+dis);}

// compute densfct0var
for(int l=0; l<ndensityfct; l++)
{densfct0var(l, k+dks, j+djs, i+dis) = densfctvar(l, k+dks, j+djs, i+dis) / densvar(0, k+dks, j+djs, i+dis);}
       });

     }

   void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
   {

   compute_diagnostic_quantities(
   diagnostic_vars.fields_arr[QDIAGVAR].data, diagnostic_vars.fields_arr[DENSLDIAGVAR].data, diagnostic_vars.fields_arr[DENSFCTLDIAGVAR].data,
   x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data, const_vars.fields_arr[CORIOLISVAR].data);
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

  SArray<real,reconstruction_order,2> primal_to_gll;
  SArray<real,reconstruction_order,reconstruction_order,reconstruction_order> primal_wenoRecon;
  SArray<real,(reconstruction_order-1)/2+2> primal_wenoIdl;
  real primal_wenoSigma;

  SArray<real,dual_reconstruction_order,2> dual_to_gll;
  SArray<real,dual_reconstruction_order,dual_reconstruction_order,dual_reconstruction_order> dual_wenoRecon;
  SArray<real,(dual_reconstruction_order-1)/2+2> dual_wenoIdl;
  real dual_wenoSigma;

  SArray<real,coriolis_reconstruction_order,2> coriolis_to_gll;
  SArray<real,coriolis_reconstruction_order,coriolis_reconstruction_order,coriolis_reconstruction_order> coriolis_wenoRecon;
  SArray<real,(coriolis_reconstruction_order-1)/2+2> coriolis_wenoIdl;
  real coriolis_wenoSigma;

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

    trans.coefs_to_gll_lower( primal_to_gll );
    trans.weno_sten_to_coefs(primal_wenoRecon);
    wenoSetIdealSigma<reconstruction_order>(primal_wenoIdl,primal_wenoSigma);

    trans.coefs_to_gll_lower( dual_to_gll );
    trans.weno_sten_to_coefs(dual_wenoRecon);
    wenoSetIdealSigma<dual_reconstruction_order>(dual_wenoIdl,dual_wenoSigma);

    trans.coefs_to_gll_lower( coriolis_to_gll );
    trans.weno_sten_to_coefs(coriolis_wenoRecon);
    wenoSetIdealSigma<coriolis_reconstruction_order>(coriolis_wenoIdl,coriolis_wenoSigma);

    PVPE.initialize(params);
    Hk.initialize(params, *this->primal_geometry, *this->dual_geometry);
    Hs.initialize(params, thermo, *this->primal_geometry, *this->dual_geometry);
    
    this->is_initialized = true;
  }


  void compute_constants(VariableSet<nconst> &const_vars, VariableSet<nprog> &x)
  {}
 
   void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
    realArr Uvar, realArr Q0var, realArr f0var, realArr dens0var, realArr densfct0var,
    const realArr Vvar, const realArr densvar, const realArr densfctvar, const realArr coriolisvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

      yakl::parallel_for("ComputeDiagIp", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        real hv;
        yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

        // compute dens0var = I densvar, densfct0var = I densfctvar
        compute_I<ndensity, diff_ord>(dens0var, densvar, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);
        compute_I<ndensityfct, diff_ord>(densfct0var, densfctvar, *this->primal_geometry, *this->dual_geometry, pis, pjs, pks, i, j, k);

});

yakl::parallel_for("ComputeDiagId", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  real hv;
  yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

  // compute U = H v, q0, f0
  compute_H<1, diff_ord>(Uvar, Vvar, *this->primal_geometry, *this->dual_geometry, dis, djs, dks, i, j, k);
  PVPE.compute_q0f0(Q0var, f0var, Vvar, densvar, densfctvar, coriolisvar, dis, djs, dks, i, j, k);

      });

    }

    void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
      realArr Fvar, realArr Kvar, realArr HEvar, const realArr Vvar, const realArr Uvar, const realArr dens0var, const realArr densfct0var) {


        int dis = dual_topology->is;
        int djs = dual_topology->js;
        int dks = dual_topology->ks;

        yakl::parallel_for("ComputeDiagII", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
        Hk.compute_dKdv(Fvar, Kvar, HEvar, Vvar, Uvar, dens0var, densfct0var, dis, djs, dks, i, j, k);
      });

      }


  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    realArr FTvar, realArr Bvar, realArr Bfctvar,
    const realArr Fvar, const realArr Uvar,
    const realArr Kvar, const realArr dens0var, const realArr densfct0var, const realArr HSvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

      yakl::parallel_for("ComputeDiagIII", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

      compute_W(FTvar, Fvar, pis, pjs, pks, i, j, k);
Hs.compute_dHsdx(Bvar, Bfctvar, dens0var, densfct0var, HSvar, pis, pjs, pks, i, j, k);
Hk.compute_dKddens(Bvar, Bfctvar, Kvar, pis, pjs, pks, i, j, k);

    });

    }




  void YAKL_INLINE compute_edge_reconstructions(realArr densedgereconvar, realArr densfctedgereconvar, realArr Qedgereconvar, realArr fedgereconvar,
    const realArr dens0var, const realArr densfct0var, const realArr Q0var, const realArr f0var) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

    yakl::parallel_for("ComputePrimalEdgeRecon", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
      
      compute_straight_edge_recon<1, reconstruction_type, reconstruction_order>(Qedgereconvar, Q0var, pis, pjs, pks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
      compute_straight_edge_recon<1, coriolis_reconstruction_type, coriolis_reconstruction_order>(fedgereconvar, f0var, pis, pjs, pks, i, j, k, coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl, coriolis_wenoSigma);

    });

    yakl::parallel_for("ComputeDualEdgeRecon", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

      compute_twisted_edge_recon<ndensity, dual_reconstruction_type, dual_reconstruction_order>(densedgereconvar, dens0var, dis, djs, dks, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
      compute_twisted_edge_recon<ndensityfct, dual_reconstruction_type, dual_reconstruction_order>(densfctedgereconvar, densfct0var, dis, djs, dks, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
    });



  }

  void YAKL_INLINE compute_recons(
  realArr densreconvar, realArr densfctreconvar, realArr Qreconvar, realArr Coriolisreconvar,
  const realArr densedgereconvar, const realArr densfctedgereconvar, const realArr Qedgereconvar, const realArr fedgereconvar, const realArr HEvar,
  const realArr FTvar, const realArr Uvar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

    yakl::parallel_for("ComputePrimalRecon", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

      compute_straight_recon<1, reconstruction_type>(Qreconvar, Qedgereconvar, FTvar, pis, pjs, pks, i, j, k);
      compute_straight_recon<1, coriolis_reconstruction_type>(Coriolisreconvar, fedgereconvar, FTvar, pis, pjs, pks, i, j, k);
      
    });

    yakl::parallel_for("ComputeDualRecon", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

      compute_twisted_recon<ndensity, dual_reconstruction_type>(densreconvar, densedgereconvar, Uvar, dis, djs, dks, i, j, k);
      compute_twisted_recon<ndensityfct, dual_reconstruction_type>(densfctreconvar, densfctedgereconvar, Uvar, dis, djs, dks, i, j, k);

    //scale primal recons
    for (int d=0;d<ndims;d++) {
    for (int l=0;l<ndensity;l++) {
    densreconvar(l+d*ndensity,k+dks,j+djs,i+dis) = densreconvar(l+d*ndensity,k+dks,j+djs,i+dis) / HEvar(d,k+dks,j+djs,i+dis);
  }}
  for (int d=0;d<ndims;d++) {
  for (int l=0;l<ndensityfct;l++) {
  densfctreconvar(l+d*ndensityfct,k+dks,j+djs,i+dis) = densfctreconvar(l+d*ndensityfct,k+dks,j+djs,i+dis) / HEvar(d,k+dks,j+djs,i+dis);
}}
    });

}


  void YAKL_INLINE compute_tendencies(
  realArr denstendvar, realArr densfcttendvar, realArr Vtendvar,
  const realArr densreconvar, const realArr densfctreconvar, const realArr Qreconvar, const realArr Coriolisreconvar,
  const realArr Bvar, const realArr Bfctvar, const realArr Fvar, const realArr Phivar) {

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

      yakl::parallel_for("ComputePrimalTendencies", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);

    compute_wD1<ndensity> (Vtendvar, densreconvar, Bvar, pis, pjs, pks, i, j, k);
    compute_wD1_fct<ndensityfct, ADD_MODE::ADD> (Vtendvar, densfctreconvar, Phivar, Bfctvar, pis, pjs, pks, i, j, k);

    if (qf_choice == QF_MODE::EC)
    { compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, pis, pjs, pks, i, j, k);}
    if (qf_choice == QF_MODE::NOEC)
    { compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, pis, pjs, pks, i, j, k);}
    compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Coriolisreconvar, Fvar, pis, pjs, pks, i, j, k);


});


yakl::parallel_for("ComputeDualTendencies", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);

  compute_wDbar2<ndensity> (denstendvar, densreconvar, Fvar, dis, djs, dks, i, j, k);
  compute_wDbar2_fct<ndensityfct> (densfcttendvar, densfctreconvar, Phivar, Fvar, dis, djs, dks, i, j, k);

  });

  }




  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

      //Compute U, q0, hf, dens0, densfct0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data, const_vars.fields_arr[CORIOLISVAR].data);

      this->aux_exchange->exchanges_arr[DENSFCT0VAR].exchange_field(auxiliary_vars.fields_arr[DENSFCT0VAR]);


      this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);
      this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(auxiliary_vars.fields_arr[Q0VAR]);
      this->aux_exchange->exchanges_arr[F0VAR].exchange_field(auxiliary_vars.fields_arr[F0VAR]);


      //Compute K, F, he
      compute_functional_derivatives_and_diagnostic_quantities_II(
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      x.fields_arr[VVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data);

      this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
      this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
      this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);

      //Compute FT, B, Bfct
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data, const_vars.fields_arr[HSVAR].data);

      this->aux_exchange->exchanges_arr[BFCTVAR].exchange_field(auxiliary_vars.fields_arr[BFCTVAR]);
      this->aux_exchange->exchanges_arr[FTVAR].exchange_field(auxiliary_vars.fields_arr[FTVAR]);
      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);

      // Compute densrecon, densfctrecon, qrecon and frecon
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data);

      this->aux_exchange->exchanges_arr[DENSFCTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR]);

      compute_recons(
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[UVAR].data);
      
      this->aux_exchange->exchanges_arr[DENSFCTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTRECONVAR]);
      this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);
      this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISRECONVAR]);


//Compute fct

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

yakl::parallel_for("ComputeEdgeFlux", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
compute_edgefluxes<ndensityfct> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, dis, djs, dks, i, j, k);
});
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


yakl::parallel_for("ComputeMf", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
compute_Mf<ndensityfct> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, dis, djs, dks, i, j, k);
});

this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

yakl::parallel_for("ComputePhi", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
compute_Phi<ndensityfct> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSFCTVAR].data, dis, djs, dks, i, j, k);
});

this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);

      // Compute tendencies
      compute_tendencies(
      xtend.fields_arr[DENSVAR].data, xtend.fields_arr[DENSFCTVAR].data, xtend.fields_arr[VVAR].data,
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
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
      SArray<real,1> pvlocal, pvglobal;
      SArray<real,1> pelocal, peglobal;

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
      for (int l=0;l<ndensityfct;l++) {massfctlocal(l) = 0.; massfctglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {densfctmaxlocal(l) = 0.; densfctmaxglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {densfctminlocal(l) = 0.; densfctminglobal(l) = 0.;}

int dis = dual_topology->is;
int djs = dual_topology->js;
int dks = dual_topology->ks;

      yakl::parallel_for("ComputeDualStats", dual_topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, dual_topology->n_cells_z, dual_topology->n_cells_y, dual_topology->n_cells_x, k, j, i);
         real KE, PE, IE;

KE = Hk.compute_KE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, dis, djs, dks, i, j, k);
PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, constvars.fields_arr[HSVAR].data, dis, djs, dks, i, j, k);
IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, dis, djs, dks, i, j, k);

elocal(3) += IE;
elocal(2) += PE;
elocal(1) += KE;
elocal(0) += KE + PE + IE;

});

int pis = primal_topology->is;
int pjs = primal_topology->js;
int pks = primal_topology->ks;

yakl::parallel_for("ComputePrimalStats", primal_topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, primal_topology->n_cells_z, primal_topology->n_cells_y, primal_topology->n_cells_x, k, j, i);
  
   pvpe vals_pvpe;
   vals_pvpe = PVPE.compute_PVPE(progvars.fields_arr[VVAR].data, progvars.fields_arr[DENSVAR].data, progvars.fields_arr[DENSFCTVAR].data, constvars.fields_arr[CORIOLISVAR].data, pis, pjs, pks, i, j, k);
   pvlocal(0) += vals_pvpe.pv;
    pelocal(0) += vals_pvpe.pe;
    
    });

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
    this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Ireduce( &elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);

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




// *******   VariableSet Initialization   ***********//
template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &ptopo, const Topology &dtopo,
SArray<int, nprog, 4> &prog_ndofs_arr, SArray<int, nconst, 4> &const_ndofs_arr, SArray<int, naux, 4> &aux_ndofs_arr, SArray<int, ndiag, 4> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{

  //primal grid represents straight quantities, dual grid twisted quantities
  
  // v, dens
  prog_topo_arr[VVAR] = &ptopo;
  prog_topo_arr[DENSVAR] = &dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[DENSVAR] = "dens";
  prog_ndofs_arr(VVAR,1) = 1; //v = straight 1-form
  prog_ndofs_arr(DENSVAR,2) = ndensity; //dens = twisted 2-form

  // hs, coriolis
  const_topo_arr[HSVAR] = &dtopo;
  const_topo_arr[CORIOLISVAR] = &ptopo;
  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";
  const_ndofs_arr(HSVAR,2) = 1; //hs = twisted 2-form
  const_ndofs_arr(CORIOLISVAR,2) = 1; //f = straight 2-form

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
  aux_ndofs_arr(BVAR,0) = ndensity; //B = straight 0-form
  aux_ndofs_arr(KVAR,2) = 1; //K = twisted 2-form
  aux_ndofs_arr(FVAR,1) = 1; //F = twisted 1-form
  aux_ndofs_arr(UVAR,1) = 1; //U = twisted 1-form
  aux_ndofs_arr(HEVAR,1) = 1; //he lives on dual edges, associated with F

  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &dtopo;
  aux_topo_arr[DENSEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENS0VAR] = &ptopo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  aux_ndofs_arr(DENSRECONVAR,1) = ndensity; //densrecon lives on dual edges, associated with F
  aux_ndofs_arr(DENSEDGERECONVAR,2) = 4*ndensity; //densedgerecon lives on dual cells, associated with F
  aux_ndofs_arr(DENS0VAR,0) = ndensity; //dens0 = straight 0-form

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
  aux_ndofs_arr(FTVAR,1) = 1; //FT = straight 1-form
  aux_ndofs_arr(Q0VAR,0) = 1; //q0 = twisted 0-form
  aux_ndofs_arr(F0VAR,0) = 1; //f0 = twisted 0-form
  aux_ndofs_arr(QRECONVAR,1) = 1; //qrecon lives on primal edges, associated with FT
  aux_ndofs_arr(QEDGERECONVAR,2) = 4; //qedgerecon lives on primal cells
  aux_ndofs_arr(CORIOLISRECONVAR,1) = 1; //coriolisrecon lives on primal edges, associated with FT
  aux_ndofs_arr(CORIOLISEDGERECONVAR,2) = 4; //coriolisedgerecon lives on primal cells

  // q, concentration 0-forms for dens
  diag_topo_arr[QDIAGVAR] = &dtopo;
  diag_topo_arr[DENSLDIAGVAR] = &ptopo;
  diag_names_arr[QDIAGVAR] = "q";
  diag_names_arr[DENSLDIAGVAR] = "densl";
  diag_ndofs_arr(QDIAGVAR,0) = 1; //qdiag = twisted 0-form
  diag_ndofs_arr(DENSLDIAGVAR,0) = ndensity; //densldiag = straight 0-form

  //densfct stuff- densfct, BFCT, densfct0, edgereconfct, reconfct, Phi, Mf, edgeflux, concentration 0-forms for densfct
  prog_topo_arr[DENSFCTVAR] = &dtopo;
  prog_names_arr[DENSFCTVAR] = "densfct";
  prog_ndofs_arr(DENSFCTVAR,2) = ndensityfct; //densfct = twisted 2-form

  aux_topo_arr[BFCTVAR] = &ptopo;
  aux_topo_arr[DENSFCTRECONVAR] = &dtopo;
  aux_topo_arr[DENSFCTEDGERECONVAR] = &dtopo;
  aux_topo_arr[DENSFCT0VAR] = &ptopo;
  aux_topo_arr[PHIVAR] = &dtopo;
  aux_topo_arr[MFVAR] = &dtopo;
  aux_topo_arr[EDGEFLUXVAR] = &dtopo;
  aux_names_arr[BFCTVAR] = "Bfct";
  aux_names_arr[DENSFCT0VAR] = "densfct0";
  aux_names_arr[DENSFCTRECONVAR] = "densfctrecon";
  aux_names_arr[DENSFCTEDGERECONVAR] = "densfctedgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_ndofs_arr(BFCTVAR,0) = ndensityfct; //Bfct = straight 0-form
  aux_ndofs_arr(DENSFCTRECONVAR,1) = ndensityfct; //densfctrecon lives on dual edges, associated with F
  aux_ndofs_arr(DENSFCTEDGERECONVAR,2) = 4*ndensityfct; //densfctedgerecon lives on dual cells
  aux_ndofs_arr(DENSFCT0VAR,0) = ndensityfct; //densfct0 = straight 0-form

  aux_ndofs_arr(PHIVAR,1) = ndensityfct;
  aux_ndofs_arr(MFVAR,2) = ndensityfct;
  aux_ndofs_arr(EDGEFLUXVAR,1) = ndensityfct;

  diag_topo_arr[DENSFCTLDIAGVAR] = &ptopo;
  diag_names_arr[DENSFCTLDIAGVAR] = "densfctl";
  diag_ndofs_arr(DENSFCTLDIAGVAR,0) = ndensityfct; //densfctldiag = straight 0-form

}

#endif
