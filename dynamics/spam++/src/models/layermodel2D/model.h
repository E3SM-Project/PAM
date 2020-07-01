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

#ifdef _SWE
uint constexpr ndensity = 1;
#endif
#ifdef _TSWE
uint constexpr ndensity = 2;
#endif

#if NTRACERS > 0
uint constexpr ndensityfct = NTRACERS;
#endif

// Number of variables
// v, dens, densfct
#if NTRACERS > 0
uint constexpr nprognostic = 3;
#define DENSFCTVAR 2
#else
uint constexpr nprognostic = 2;
#endif

#define VVAR 0
#define DENSVAR 1


// hs, coriolis
uint constexpr nconstant = 2;
#define HSVAR 0
#define CORIOLISVAR 1

//functional derivatives = F, B, BFCT, K, he
//dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon, coriolisedgercon, coriolisrecon
//primal grid reconstruction stuff- U, dens0, densfct0, edgerecon, recon, edgereconfct, reconfct
//fct stuff- Phi, Mf, edgeflux

#if NTRACERS > 0
uint constexpr nauxiliary = 22;
#else
uint constexpr nauxiliary = 15;
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

#if NTRACERS > 0
#define BFCTVAR 15
#define DENSFCT0VAR 16
#define DENSFCTRECONVAR 17
#define DENSFCTEDGERECONVAR 18
#define PHIVAR 19
#define EDGEFLUXVAR 20
#define MFVAR 21
#endif

// q, associated concentration 0-forms for den

#if NTRACERS > 0
uint constexpr ndiagnostic = 3;
#define DENSFCTLDIAGVAR 2
#else
uint constexpr ndiagnostic = 2;
#endif
#define QDIAGVAR 0
#define DENSLDIAGVAR 1


//track total densities, densfct mix/max, energy (P+K+I+T), PV, PE,
#if NTRACERS > 0
uint constexpr nstats = 7;
#else
uint constexpr nstats = 4;
#endif

#define DENSSTAT 0
#define ESTAT 1
#define PVSTAT 2
#define PESTAT 3
#if NTRACERS > 0
#define DENSFCTSTAT 4
#define DENSFCTMAXSTAT 5
#define DENSFCTMINSTAT 6
#endif





// *******   Model Specific Parameters   ***********//

// THESE ARE DOUBLE VORTEX SPECIFIC...
real constexpr Lx = 5000. * 1000.;
real constexpr Ly = 5000. * 1000.;
real constexpr H0 = 750.0;
real constexpr ox = 0.1;
real constexpr oy = 0.1;
real constexpr sigmax = 3./40.*Lx;
real constexpr sigmay = 3./40.*Lx;
real constexpr dh = 75.0;
real constexpr g = 9.80616;
real constexpr coriolis = 0.00006147;
real constexpr xc1 = (0.5-ox) * Lx;
real constexpr yc1 = (0.5-oy) * Ly;
real constexpr xc2 = (0.5+ox) * Lx;
real constexpr yc2 = (0.5+oy) * Ly;
real constexpr xc = 0.5 * Lx;
real constexpr yc = 0.5 * Ly;
real constexpr c = 0.05;
real constexpr a = 1.0/3.0;
real constexpr D = 0.5 * Lx;

void set_model_specific_params(std::string inFile, ModelParameters &params)
{


  std::string strDataInit = "";

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
           if ( !strcmp( "dataInit"   , key.c_str() ) ) { ssVal >> strDataInit       ;}
      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

  if (!strcmp("",strDataInit.c_str())) { std::cout << "Error: key " << "dataInit" << " not set.\n"; exit(-1); }

    size_t splitloc = strDataInit.find("//",0);
    std::string sub_str;
    if (splitloc != std::string::npos){
      sub_str = strDataInit.substr(0,splitloc);
    } else {
      sub_str = strDataInit;
    }
    if      ( !strcmp(sub_str.c_str(),"doublevortex" ) ) { params.data_init_cond = DATA_INIT::DOUBLEVORTEX  ; }
    else if ( !strcmp(sub_str.c_str(),"golo" ) ) { params.data_init_cond = DATA_INIT::GOLO  ; }
    else if ( !strcmp(sub_str.c_str(),"tc2"   ) ) { params.data_init_cond = DATA_INIT::TC2    ; }
    else if ( !strcmp(sub_str.c_str(),"tc5"   ) ) { params.data_init_cond = DATA_INIT::TC5    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataInit << "\n";
      exit(-1);
    }

  params.etime = 0.0;


  if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
  {
  params.xlen = Lx;
  params.xc = Lx/2.;
  params.ylen = Ly;
  params.yc = Ly/2.;
  }

  // if (params.data_init_cond == DATA_INIT::GOLO)
  // {
  // params.xlen = 1.0;
  // params.xc = 0.5;
  // params.ylen = 1.0;
  // params.yc = 0.5;
  // }
  //
  // if (params.data_init_cond == DATA_INIT::TC2 || params.data_init_cond == DATA_INIT::TC5)
  // {
  // params.xlen = 1.0;
  // params.xc = 0.5;
  // params.ylen = 1.0;
  // params.yc = 0.5;
  // }

}

// ******* Diagnostics *************//

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE

template <uint nprog, uint nconst, uint ndiag> class Diagnostics {
public:

  const Topology *topology;
  Geometry<ndims,1,1,1> *geom;

  bool is_initialized;

   Diagnostics() {
     this->is_initialized = false;
     std::cout << "CREATED DIAGNOSTICS\n";
   }

   void initialize(const Topology &topo, Geometry<ndims,1,1,1> &geom)
   {
     this->topology = &topo;
     this->geom = &geom;
     this->is_initialized = true;
   }


#if NTRACERS > 0
   void YAKL_INLINE compute_diagnostic_quantities(
     realArr Q0var, realArr dens0var, realArr densfct0var,
     const realArr Vvar, const realArr densvar, const realArr densfctvar, const realArr coriolisvar) {
#else
void YAKL_INLINE compute_diagnostic_quantities(
  realArr Q0var, realArr dens0var,
  const realArr Vvar, const realArr densvar, const realArr coriolisvar) {
#endif

       int is = topology->is;
       int js = topology->js;
       int ks = topology->ks;

       yakl::parallel_for("ComputeDiagI", topology->n_cells, YAKL_LAMBDA (int iGlob) {
         SArray<real,1> zeta;
         real hv;

         int k, j, i;
         yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

 // compute zeta = D2 v
 compute_D2<1>(zeta, Vvar, is, js, ks, i, j, k);

 // compute q0 = zeta / R h
   compute_R<0> (hv, densvar, is, js, ks, i, j, k);
   Q0var(0, k+ks, j+js, i+is) = (zeta(0) + coriolisvar(0, k+ks, j+js, i+is)) / hv;

// compute dens0var
for(int l=0; l<ndensity; l++)
{dens0var(l, k+ks, j+js, i+is) = densvar(l, k+ks, j+js, i+is) / densvar(0, k+ks, j+js, i+is);}

// compute densfct0var
#if NTRACERS > 0
for(int l=0; l<ndensityfct; l++)
{densfct0var(l, k+ks, j+js, i+is) = densfctvar(l, k+ks, j+js, i+is) / densvar(0, k+ks, j+js, i+is);}
#endif
       });

     }



   void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
   {

#if NTRACERS > 0
   compute_diagnostic_quantities(
   diagnostic_vars.fields_arr[QDIAGVAR].data, diagnostic_vars.fields_arr[DENSLDIAGVAR].data, diagnostic_vars.fields_arr[DENSFCTLDIAGVAR].data,
   x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data, const_vars.fields_arr[CORIOLISVAR].data);
#else
compute_diagnostic_quantities(
diagnostic_vars.fields_arr[QDIAGVAR].data, diagnostic_vars.fields_arr[DENSLDIAGVAR].data,
x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISVAR].data);
#endif
}

};
// *******   Tendencies   ***********//

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE
template <uint nprog, uint nconst, uint naux> class Tendencies {
public:

  const Topology *topology;
  ExchangeSet<naux> *aux_exchange;
  ExchangeSet<nconst> *const_exchange;
  Geometry<ndims,1,1,1> *geom;

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

  void initialize(const Topology &topo, Geometry<ndims,1,1,1> &geom, ExchangeSet<naux> &aux_exchange, ExchangeSet<nconst> &const_exchange)
  {
    this->topology = &topo;
    this->geom = &geom;
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

    this->is_initialized = true;
  }




#if NTRACERS > 0
  void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
    realArr Uvar, realArr Q0var, realArr f0var, realArr dens0var, realArr densfct0var,
    const realArr Vvar, const realArr densvar, const realArr densfctvar, const realArr coriolisvar) {
#else
void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
  realArr Uvar, realArr Q0var, realArr f0var, realArr dens0var,
  const realArr Vvar, const realArr densvar, const realArr coriolisvar) {
#endif
      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeDiagI", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        real hv;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

        // compute dens0var = I densvar, U = H v, densfct0var = I densfctvar
        compute_H<1, diff_ord>(Uvar, Vvar, *this->geom, is, js, ks, i, j, k);
        compute_I<ndensity, diff_ord>(dens0var, densvar, *this->geom, is, js, ks, i, j, k);
#if NTRACERS > 0
        compute_I<ndensityfct, diff_ord>(densfct0var, densfctvar, *this->geom, is, js, ks, i, j, k);
#endif

// compute zeta = D2 v
compute_D2<1>(Q0var, Vvar, is, js, ks, i, j, k);

// compute q0 = zeta / R h
  compute_R<0> (hv, densvar, is, js, ks, i, j, k);
  Q0var(0, k+ks, j+js, i+is) = Q0var(0, k+ks, j+js, i+is) / hv;
  f0var(0, k+ks, j+js, i+is) = coriolisvar(0, k+ks, j+js, i+is) / hv;

      });

    }

    void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
      realArr Fvar, realArr Kvar, realArr HEvar, const realArr Vvar, const realArr Uvar, const realArr dens0var) {

        int is = topology->is;
        int js = topology->js;
        int ks = topology->ks;

        yakl::parallel_for("ComputeDiagII", topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

        // compute he = phi * h0
        compute_phi<0>(HEvar, dens0var, is, js, ks, i, j, k);

        //compute F = he * U
        Fvar(0, k+ks, j+js, i+is) = Uvar(0, k+ks, j+js, i+is) * HEvar(0, k+ks, j+js, i+is);
        Fvar(1, k+ks, j+js, i+is) = Uvar(1, k+ks, j+js, i+is) * HEvar(1, k+ks, j+js, i+is);

        //compute KE = 0.5 * phiT(u,v)
        compute_phiT(Kvar, Uvar, Vvar, is, js, ks, i, j, k);
        Kvar(0, k+ks, j+js, i+is) *= 0.5;

      });

      }


#if NTRACERS > 0
  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    realArr FTvar, realArr Bvar, realArr Bfctvar,
    const realArr Fvar, const realArr Uvar,
    const realArr Kvar, const realArr dens0var, const realArr densfct0var, const realArr HSvar) {
#else
void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
  realArr FTvar, realArr Bvar,
  const realArr Fvar, const realArr Uvar,
  const realArr Kvar, const realArr dens0var, const realArr HSvar) {
#endif
      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeDiagIII", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        SArray<real,1> hs0;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

      compute_W(FTvar, Fvar, is, js, ks, i, j, k);

#ifdef _TSWE
      // Compute B = IK + S0/2.
      compute_I<1, diff_ord>(Bvar, Kvar, *this->geom, is, js, ks, i, j, k);
      Bvar(0, k+ks, j+js, i+is) += dens0var(1, k+ks, j+js, i+is)/2.;
      //Compute T = Ihs + h0/2;
      compute_I<1, diff_ord>(hs0, HSvar, *this->geom, is, js, ks, i, j, k);
      Bvar(1, k+ks, j+js, i+is) = dens0var(0, k+ks, j+js, i+is)/2. + hs0(0);
#endif

#ifdef _SWE
      // Compute B = IK + gh0 + ghs0
      compute_I<1, diff_ord>(Bvar, Kvar, *this->geom, is, js, ks, i, j, k);
      compute_I<1, diff_ord>(hs0, HSvar, *this->geom, is, js, ks, i, j, k);
      Bvar(0, k+ks, j+js, i+is) += g*dens0var(0, k+ks, j+js, i+is) + g*hs0(0);
#endif

#if NTRACERS > 0
// Compute Xil = h0/2; and add \sum trl/2. to B
for (int l=0; l<ndensityfct; l++)
{
 Bvar(0, k+ks, j+js, i+is) += densfct0var(l, k+ks, j+js, i+is)/2.;
  Bfctvar(l, k+ks, j+js, i+is) = dens0var(0, k+ks, j+js, i+is)/2.;
}
#endif

    });

    }




#if NTRACERS > 0
  void YAKL_INLINE compute_edge_reconstructions(realArr densedgereconvar, realArr densfctedgereconvar, realArr Qedgereconvar, realArr fedgereconvar,
    const realArr dens0var, const realArr densfct0var, const realArr Q0var, const realArr f0var) {
#else
  void YAKL_INLINE compute_edge_reconstructions(realArr densedgereconvar, realArr Qedgereconvar, realArr fedgereconvar,
    const realArr dens0var, const realArr Q0var, const realArr f0var) {
#endif
    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

    yakl::parallel_for("ComputeEdgeRecon", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);


      compute_primal_edge_recon<ndensity, reconstruction_type, reconstruction_order>(densedgereconvar, dens0var, is, js, ks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
#if NTRACERS > 0
      compute_primal_edge_recon<ndensityfct, reconstruction_type, reconstruction_order>(densfctedgereconvar, densfct0var, is, js, ks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
#endif
      compute_dual_edge_recon<1, dual_reconstruction_type, dual_reconstruction_order>(Qedgereconvar, Q0var, is, js, ks, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
      compute_dual_edge_recon<1, coriolis_reconstruction_type, coriolis_reconstruction_order>(fedgereconvar, f0var, is, js, ks, i, j, k, coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl, coriolis_wenoSigma);

    });

  }

#if NTRACERS > 0
  void YAKL_INLINE compute_recons(
  realArr densreconvar, realArr densfctreconvar, realArr Qreconvar, realArr Coriolisreconvar,
  const realArr densedgereconvar, const realArr densfctedgereconvar, const realArr Qedgereconvar, const realArr fedgereconvar, const realArr HEvar,
  const realArr FTvar, const realArr Uvar) {
#else
void YAKL_INLINE compute_recons(
realArr densreconvar, realArr Qreconvar, realArr Coriolisreconvar,
const realArr densedgereconvar, const realArr Qedgereconvar, const realArr fedgereconvar, const realArr HEvar,
const realArr FTvar, const realArr Uvar) {
#endif

    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

    yakl::parallel_for("ComputeRecon", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

      compute_primal_recon<ndensity, reconstruction_type>(densreconvar, densedgereconvar, Uvar, is, js, ks, i, j, k);
#if NTRACERS > 0
      compute_primal_recon<ndensityfct, reconstruction_type>(densfctreconvar, densfctedgereconvar, Uvar, is, js, ks, i, j, k);
#endif
      compute_dual_recon<1, dual_reconstruction_type>(Qreconvar, Qedgereconvar, FTvar, is, js, ks, i, j, k);
      compute_dual_recon<1, coriolis_reconstruction_type>(Coriolisreconvar, fedgereconvar, FTvar, is, js, ks, i, j, k);

    //scale primal recons
    for (int d=0;d<ndims;d++) {
    for (int l=0;l<ndensity;l++) {
    densreconvar(l+d*ndensity,k+ks,j+js,i+is) = densreconvar(l+d*ndensity,k+ks,j+js,i+is) / HEvar(d,k+ks,j+js,i+is);
  }}
#if NTRACERS > 0
  for (int d=0;d<ndims;d++) {
  for (int l=0;l<ndensityfct;l++) {
  densfctreconvar(l+d*ndensityfct,k+ks,j+js,i+is) = densfctreconvar(l+d*ndensityfct,k+ks,j+js,i+is) / HEvar(d,k+ks,j+js,i+is);
}}
#endif
    });

}


#if NTRACERS > 0
  void YAKL_INLINE compute_tendencies(
  realArr denstendvar, realArr densfcttendvar, realArr Vtendvar,
  const realArr densreconvar, const realArr densfctreconvar, const realArr Qreconvar, const realArr Coriolisreconvar,
  const realArr Bvar, const realArr Bfctvar, const realArr Fvar, const realArr Phivar) {
#else
void YAKL_INLINE compute_tendencies(
realArr denstendvar, realArr Vtendvar,
const realArr densreconvar, const realArr Qreconvar, const realArr Coriolisreconvar,
const realArr Bvar, const realArr Fvar) {
#endif
    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

      yakl::parallel_for("ComputeTendencies", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

    compute_wDbar2<ndensity> (denstendvar, densreconvar, Fvar, is, js, ks, i, j, k);
    compute_wD1<ndensity> (Vtendvar, densreconvar, Bvar, is, js, ks, i, j, k);

#if NTRACERS > 0
    compute_wDbar2_fct<ndensityfct> (densfcttendvar, densfctreconvar, Phivar, Fvar, is, js, ks, i, j, k);
    compute_wD1_fct<ndensityfct, ADD_MODE::ADD> (Vtendvar, densfctreconvar, Phivar, Bfctvar, is, js, ks, i, j, k);
#endif

    if (qf_choice == QF_MODE::EC)
    { compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, is, js, ks, i, j, k);}
    if (qf_choice == QF_MODE::NOEC)
    { compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, is, js, ks, i, j, k);}
    compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Coriolisreconvar, Fvar, is, js, ks, i, j, k);

  });

  }




  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

      //Compute U, q0, hf, dens0, densfct0
#if NTRACERS > 0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data, const_vars.fields_arr[CORIOLISVAR].data);

      this->aux_exchange->exchanges_arr[DENSFCT0VAR].exchange_field(auxiliary_vars.fields_arr[DENSFCT0VAR]);
#else
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, const_vars.fields_arr[CORIOLISVAR].data);
#endif

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

      //Compute FT, B, Bfct
#if NTRACERS > 0
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data, const_vars.fields_arr[HSVAR].data);

      this->aux_exchange->exchanges_arr[BFCTVAR].exchange_field(auxiliary_vars.fields_arr[BFCTVAR]);
#else
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[BVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, const_vars.fields_arr[HSVAR].data);
#endif

      this->aux_exchange->exchanges_arr[FTVAR].exchange_field(auxiliary_vars.fields_arr[FTVAR]);
      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);

      // Compute densrecon, densfctrecon, qrecon and frecon
#if NTRACERS > 0
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data);

      this->aux_exchange->exchanges_arr[DENSFCTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR]);
#else
compute_edge_reconstructions(
auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[F0VAR].data);
#endif

      this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR]);

#if NTRACERS > 0
      compute_recons(
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[UVAR].data);
      this->aux_exchange->exchanges_arr[DENSFCTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTRECONVAR]);
#else
    compute_recons(
    auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
    auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
    auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
    auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[UVAR].data);
#endif
      this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);
      this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISRECONVAR]);


//Compute fct
#if NTRACERS > 0

int is = topology->is;
int js = topology->js;
int ks = topology->ks;

yakl::parallel_for("ComputeEdgeFlux", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_edgefluxes<ndensityfct> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, is, js, ks, i, j, k);
});
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


yakl::parallel_for("ComputeMf", topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_Mf<ndensityfct> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, is, js, ks, i, j, k);
});

this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

yakl::parallel_for("ComputePhi", topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_Phi<ndensityfct> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSFCTVAR].data, is, js, ks, i, j, k);
});


this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);
#endif

      // Compute tendencies
#if NTRACERS > 0
      compute_tendencies(
      xtend.fields_arr[DENSVAR].data, xtend.fields_arr[DENSFCTVAR].data, xtend.fields_arr[VVAR].data,
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[PHIVAR].data);
#else
  compute_tendencies(
  xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
  auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
  auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[FVAR].data);
#endif
};

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
  Geometry<ndims,1,1,1> *geom;
  const Topology *topology;

  void initialize(ModelParameters &params, Parallel &par, const Topology &topo, Geometry<ndims,1,1,1> &geom)
  {
      this->topology = &topo;
      this->geom = &geom;
    statsize = params.Nsteps/params.Nstat + 1;

    stats_arr[DENSSTAT].initialize("mass", ndensity, params, par);
#if NTRACERS > 0
    stats_arr[DENSFCTSTAT].initialize("massfct", ndensityfct, params, par);
    stats_arr[DENSFCTMAXSTAT].initialize("fctmax", ndensityfct, params, par);
    stats_arr[DENSFCTMINSTAT].initialize("fctmin", ndensityfct, params, par);
#endif
    stats_arr[ESTAT].initialize("energy", 3, params, par);
    stats_arr[PVSTAT].initialize("pv", 1, params, par);
    stats_arr[PESTAT].initialize("pens", 1, params, par);
    masterproc = par.masterproc;
  }



  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int i)
  {

      SArray<real,ndensity> masslocal, massglobal;
#if NTRACERS > 0
      SArray<real,ndensityfct> massfctlocal, massfctglobal;
      SArray<real,ndensityfct> maxfctlocal, maxfctglobal;
      SArray<real,ndensityfct> minfctlocal, minfctglobal;
#endif
      SArray<real,3> elocal, eglobal;
      SArray<real,1> pvlocal, pvglobal;
      SArray<real,1> pelocal, peglobal;

      pvlocal(0) = 0.;
      pvglobal(0) = 0.;
      pelocal(0) = 0.;
      peglobal(0) = 0.;
      elocal(0) = 0.;
      elocal(1) = 0.;
      elocal(2) = 0.;
      eglobal(0) = 0.;
      eglobal(1) = 0.;
      eglobal(2) = 0.;
      for (int l=0;l<ndensity;l++) {masslocal(l) = 0.; massglobal(l) = 0.;}
#if NTRACERS > 0
      for (int l=0;l<ndensityfct;l++) {massfctlocal(l) = 0.; massfctglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {maxfctlocal(l) = 0.; maxfctglobal(l) = 0.;}
      for (int l=0;l<ndensityfct;l++) {minfctlocal(l) = 0.; minfctglobal(l) = 0.;}
#endif
      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeStats", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
        real eta, hv, q0, KE, PE;
        SArray<real,ndensity> dens0;
        SArray<real,1> zeta, h0im1, h0jm1;
        SArray<real,2> U, he;
        SArray<real,2,2> h0arr;

        //compute stats locally
        compute_I<ndensity,diff_ord> (dens0, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i, j, k);
        compute_H<1,diff_ord> (U, progvars.fields_arr[VVAR].data, *this->geom, is, js, ks, i, j, k);

        //Slight hack that only computes I for first dof
        compute_I<1,diff_ord> (h0im1, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i-1, j, k);
        compute_I<1,diff_ord> (h0jm1, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i, j-1, k);

h0arr(0,0) = dens0(0);
h0arr(1,0) = dens0(0);
h0arr(0,1) = h0im1(0);
h0arr(1,1) = h0jm1(0);
phi(he, h0arr);

        KE = 1./2. * (he(0) * ( U(0) * progvars.fields_arr[VVAR].data(0,k+ks,j+js,i+is)) +
                    + he(1) * ( U(1) * progvars.fields_arr[VVAR].data(1,k+ks,j+js,i+is)));

#ifdef _TSWE
PE = 0.5*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(1) + dens0(1)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif
#ifdef _SWE
PE = PE = 0.5*g*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(0) + g*dens0(0)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif

#if NTRACERS > 0
for (int l=0;l<ndensityfct;l++)
{PE += 0.5*dens0(0)*progvars.fields_arr[DENSFCTVAR].data(l,k+ks,j+js,i+is);}
#endif

        compute_D2<1>(zeta, progvars.fields_arr[VVAR].data, is, js, ks, i, j, k);
        eta = zeta(0) + constvars.fields_arr[CORIOLISVAR].data(0,k+ks,j+js,i+is);
        compute_R<0> (hv, progvars.fields_arr[DENSVAR].data, is, js, ks, i, j, k);
        q0 = eta / hv;

       elocal(0) += PE;
       elocal(1) += KE;
       elocal(2) += KE + PE;
       pvlocal(0) += eta;
       pelocal(0) += 0.5 * eta * q0;

    });

    for (int l=0;l<ndensity;l++)
    {
    masslocal(l) = progvars.fields_arr[DENSVAR].sum(l);
    }

#if NTRACERS > 0
    for (int l=0;l<ndensityfct;l++)
    {
      massfctlocal(l) = progvars.fields_arr[DENSFCTVAR].sum(l);
      maxfctlocal(l) = progvars.fields_arr[DENSFCTVAR].max(l);
      minfctlocal(l) = progvars.fields_arr[DENSFCTVAR].min(l);
    }
#endif

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &masslocal, &massglobal, ndensity, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSSTAT]);
#if NTRACERS > 0
    this->ierr = MPI_Ireduce( &massfctlocal, &massfctglobal, ndensityfct, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[DENSFCTSTAT]);
    this->ierr = MPI_Ireduce( &maxfctlocal, &maxfctglobal, ndensityfct, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSFCTMAXSTAT]);
    this->ierr = MPI_Ireduce( &minfctlocal, &minfctglobal, ndensityfct, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSFCTMINSTAT]);
#endif
    this->ierr = MPI_Ireduce( &pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Ireduce( &elocal, &eglobal, 3, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);

    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
    for (int l=0;l<ndensity;l++)
    {
  this->stats_arr[DENSSTAT].data(l,i) = massglobal(l);
}

#if NTRACERS > 0
  for (int l=0;l<ndensityfct;l++)
  {
  this->stats_arr[DENSFCTSTAT].data(l,i) = massfctglobal(l);
  this->stats_arr[DENSFCTMAXSTAT].data(l,i) = maxfctglobal(l);
  this->stats_arr[DENSFCTMINSTAT].data(l,i) = minfctglobal(l);
}
#endif
  this->stats_arr[ESTAT].data(0,i) = eglobal(0);
  this->stats_arr[ESTAT].data(1,i) = eglobal(1);
  this->stats_arr[ESTAT].data(2,i) = eglobal(2);
  this->stats_arr[PVSTAT].data(0,i) = pvglobal(0);
  this->stats_arr[PESTAT].data(0,i) = peglobal(0);
  }
  }
};




// *******   VariableSet Initialization   ***********//
template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &topo,
SArray<int, nprog, 4> &prog_ndofs_arr, SArray<int, nconst, 4> &const_ndofs_arr, SArray<int, naux, 4> &aux_ndofs_arr, SArray<int, ndiag, 4> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{

//primal grid represents twisted quantities, dual grid straight quantities

  // v, dens
  prog_topo_arr[VVAR] = &topo;
  prog_topo_arr[DENSVAR] = &topo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[DENSVAR] = "dens";
  prog_ndofs_arr(VVAR,1) = 1; //v = straight 1-form
  prog_ndofs_arr(DENSVAR,2) = ndensity; //dens = twisted 2-form

  // hs, coriolis
  const_topo_arr[HSVAR] = &topo;
  const_topo_arr[CORIOLISVAR] = &topo;
  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";
  const_ndofs_arr(HSVAR,2) = 1; //hs = twisted 2-form
  const_ndofs_arr(CORIOLISVAR,0) = 1; //f = straight 2-form

  //functional derivatives = F, B, K, he, U
  aux_topo_arr[BVAR] = &topo;
  aux_topo_arr[FVAR] = &topo;
  aux_topo_arr[UVAR] = &topo;
  aux_topo_arr[HEVAR] = &topo;
  aux_topo_arr[KVAR] = &topo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  aux_ndofs_arr(BVAR,2) = ndensity; //B = straight 0-form
  aux_ndofs_arr(KVAR,2) = 1; //K = twisted 2-form
  aux_ndofs_arr(FVAR,1) = 1; //F = twisted 1-form
  aux_ndofs_arr(UVAR,1) = 1; //U = twisted 1-form
  aux_ndofs_arr(HEVAR,1) = 1; //he lives on edges

  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &topo;
  aux_topo_arr[DENSEDGERECONVAR] = &topo;
  aux_topo_arr[DENS0VAR] = &topo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  aux_ndofs_arr(DENSRECONVAR,1) = ndensity; //densrecon lives on edges
  aux_ndofs_arr(DENSEDGERECONVAR,2) = 4*ndensity; //densedgerecon lives on cells
  aux_ndofs_arr(DENS0VAR,2) = ndensity; //dens0 = straight 0-form

  //dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon, coriolisedgercon, coriolisrecon
  aux_topo_arr[FTVAR] = &topo;
  aux_topo_arr[CORIOLISRECONVAR] = &topo;
  aux_topo_arr[CORIOLISEDGERECONVAR] = &topo;
  aux_topo_arr[Q0VAR] = &topo;
  aux_topo_arr[F0VAR] = &topo;
  aux_topo_arr[QRECONVAR] = &topo;
  aux_topo_arr[QEDGERECONVAR] = &topo;
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
  aux_ndofs_arr(QRECONVAR,1) = 1; //qrecon lives on edges
  aux_ndofs_arr(QEDGERECONVAR,0) = 4; //qedgerecon lives on dual cells
  aux_ndofs_arr(CORIOLISRECONVAR,1) = 1; //coriolisrecon lives on edges
  aux_ndofs_arr(CORIOLISEDGERECONVAR,0) = 4; //coriolisedgerecon lives on dual cells

  // q, concentration 0-forms for dens
  diag_topo_arr[QDIAGVAR] = &topo;
  diag_topo_arr[DENSLDIAGVAR] = &topo;
  diag_names_arr[QDIAGVAR] = "q";
  diag_names_arr[DENSLDIAGVAR] = "densl";
  diag_ndofs_arr(QDIAGVAR,0) = 1; //qdiag = twisted 0-form
  diag_ndofs_arr(DENSLDIAGVAR,2) = ndensity; //densldiag = straight 0-form

  //densfct stuff- densfct, BFCT, densfct0, edgereconfct, reconfct, Phi, Mf, edgeflux, concentration 0-forms for densfct
#if NTRACERS > 0
  prog_topo_arr[DENSFCTVAR] = &topo;
  prog_names_arr[DENSFCTVAR] = "densfct";
  prog_ndofs_arr(DENSFCTVAR,2) = ndensityfct; //densfct = twisted 2-form

  aux_topo_arr[BFCTVAR] = &topo;
  aux_topo_arr[DENSFCTRECONVAR] = &topo;
  aux_topo_arr[DENSFCTEDGERECONVAR] = &topo;
  aux_topo_arr[DENSFCT0VAR] = &topo;
  aux_topo_arr[PHIVAR] = &topo;
  aux_topo_arr[MFVAR] = &topo;
  aux_topo_arr[EDGEFLUXVAR] = &topo;
  aux_names_arr[BFCTVAR] = "Bfct";
  aux_names_arr[DENSFCT0VAR] = "densfct0";
  aux_names_arr[DENSFCTRECONVAR] = "densfctrecon";
  aux_names_arr[DENSFCTEDGERECONVAR] = "densfctedgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_ndofs_arr(BFCTVAR,2) = ndensityfct; //Bfct = straight 0-form
  aux_ndofs_arr(DENSFCTRECONVAR,1) = ndensityfct; //densfctrecon lives on edges
  aux_ndofs_arr(DENSFCTEDGERECONVAR,2) = 4*ndensityfct; //densfctedgerecon lives on cells
  aux_ndofs_arr(DENSFCT0VAR,2) = ndensityfct; //densfct0 = straight 0-form
  aux_ndofs_arr(PHIVAR,1) = ndensityfct;
  aux_ndofs_arr(MFVAR,2) = ndensityfct;
  aux_ndofs_arr(EDGEFLUXVAR,1) = ndensityfct;

  diag_topo_arr[DENSFCTLDIAGVAR] = &topo;
  diag_names_arr[DENSFCTLDIAGVAR] = "densfctl";
  diag_ndofs_arr(DENSFCTLDIAGVAR,2) = ndensityfct; //densfctldiag = straight 0-form
#endif


}

  // *******   Initial Conditions   ***********//

real YAKL_INLINE double_vortex_coriolis(real x, real y)
{
    return coriolis;
}

real YAKL_INLINE double_vortex_h(real x, real y)
{
    real xprime1 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc1));
    real yprime1 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc1));
    real xprime2 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc2));
    real yprime2 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc2));
    real xprimeprime1 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc1));
    real yprimeprime1 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc1));
    real xprimeprime2 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc2));
    real yprimeprime2 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc2));

    return H0 - dh * (exp(-0.5 * (xprime1 * xprime1 + yprime1 * yprime1)) + exp(-0.5 * (xprime2 * xprime2 + yprime2 * yprime2)) - 4. * M_PI * sigmax * sigmay / Lx / Ly);
}

vec<2> YAKL_INLINE double_vortex_v(real x, real y) {
  vec<2> vvec;

  real xprime1 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc1));
  real yprime1 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc1));
  real xprime2 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc2));
  real yprime2 = Ly / (M_PI * sigmay) * sin(M_PI / Ly * (y - yc2));
  real xprimeprime1 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc1));
  real yprimeprime1 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc1));
  real xprimeprime2 = Lx / (2.0 * M_PI * sigmax) * sin(2 * M_PI / Lx * (x - xc2));
  real yprimeprime2 = Ly / (2.0 * M_PI * sigmay) * sin(2 * M_PI / Ly * (y - yc2));

  vvec.u = - g * dh / coriolis / sigmay * (yprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + yprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
  vvec.v = g * dh / coriolis / sigmax * (xprimeprime1 * exp(-0.5*(xprime1 * xprime1 + yprime1 * yprime1)) + xprimeprime2 * exp(-0.5*(xprime2 * xprime2 + yprime2 * yprime2)));
  return vvec;
}

real YAKL_INLINE double_vortex_S(real x, real y)
{
    //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. * M_PI / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
    real sval = g * (1. + c * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
    //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
    //real sval = g;
    //real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35 * Ly && y < 0.65 * Ly ) ? 1. : 0.));
    return sval * double_vortex_h(x,y);
}

real YAKL_INLINE double_vortex_tracer_square_cent(real x, real y)         {return (x > 0.35*Lx && x < 0.65*Lx && y > 0.35*Ly && y < 0.65*Ly                        ) ? 0.005*double_vortex_h(x,y) : 0.;}
real YAKL_INLINE double_vortex_tracer_square_ur(real x, real y)         {return (x > 0.6*Lx && x < 0.9*Lx && y > 0.6*Ly && y < 0.9*Ly                        ) ? 0.005*double_vortex_h(x,y) : 0.;}
real YAKL_INLINE double_vortex_tracer_square_ll(real x, real y)         {return (x > 0.1*Lx && x < 0.4*Lx && y > 0.1*Ly && y < 0.4*Ly                        ) ? 0.005*double_vortex_h(x,y) : 0.;}
real YAKL_INLINE double_vortex_tracer_square_urpll(real x, real y)         {return double_vortex_tracer_square_ur(x,y) + double_vortex_tracer_square_ll(x,y);}


real YAKL_INLINE double_vortex_tracer_gaussian(real x, real y)     {return  0.005 *double_vortex_h(x,y) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D));}
//{ return 0.005 *double_vortex_h(x,y) * exp(-100. * pow((x-xc)/Lx,2.)) * exp(-100. * pow((y-yc)/Ly,2.)); }


#define vortex1_2d(x,y)   (0.005 * double_vortex_h(x,y) *  exp(-100. * pow((x-0.75*Lx)/Lx,2.)) * exp(-100. * pow((y-0.75*Ly)/Ly,2.)))
#define vortex2_2d(x,y)   (0.005 * double_vortex_h(x,y) * exp(-50.  * pow((x-0.25*Lx)/Lx,2.)) * exp(-75.  * pow((y-0.25*Ly)/Ly,2.)))
real YAKL_INLINE double_vortex_tracer_vortices(real x, real y)         { return vortex1_2d(x,y)   + vortex2_2d(x,y); }

// ADD MORE ICs HERE!!!

//wavespeed = sqrt(g * H0)
//dt = Constant(get_dt(wavespeed, cval, order, variant, Lx, Ly, nx, ny))


template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{

    if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
    {
        geom.set_primal_2form_values(double_vortex_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
        geom.set_primal_2form_values(double_vortex_S, progvars.fields_arr[DENSVAR], 1);
#endif
#if NTRACERS>0
        geom.set_primal_2form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSFCTVAR], 0);
#endif
#if NTRACERS>1
        geom.set_primal_2form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSFCTVAR], 1);
#endif
#if NTRACERS>2
        geom.set_primal_2form_values(double_vortex_tracer_square_urpll, progvars.fields_arr[DENSFCTVAR], 2);
#endif
        geom.set_dual_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
        geom.set_dual_2form_values(double_vortex_coriolis, constvars.fields_arr[CORIOLISVAR], 0);
    }

    // if (params.data_init_cond == DATA_INIT::GOLO)
    // {
    //     geom.set_primal_2form_values(gaussian, progvars.fields_arr[HVAR], 0);
    //     geom.set_dual_1form_values(gaussian, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    //     geom.set_primal_0form_values(gaussian, constvars.fields_arr[CORIOLISVAR], 0);
    // }
    //
    // if (params.data_init_cond == DATA_INIT::TC2)
    // {
    //     geom.set_primal_2form_values(gaussian, progvars.fields_arr[HVAR], 0);
    //     geom.set_dual_1form_values(gaussian, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    //     geom.set_primal_0form_values(gaussian, constvars.fields_arr[CORIOLISVAR], 0);
    // }
    //
    // if (params.data_init_cond == DATA_INIT::TC5)
    // {
    //     geom.set_primal_2form_values(gaussian, progvars.fields_arr[HVAR], 0);
    //     geom.set_dual_1form_values(gaussian, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    //     geom.set_primal_2form_values(gaussian, constvars.fields_arr[HSVAR], 0);
    //     geom.set_primal_0form_values(gaussian, constvars.fields_arr[CORIOLISVAR], 0);
    // }

}



#endif
