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
#ifdef _SWEVT
uint constexpr ndensity = 2;
#endif
#ifdef _TSWEVT
uint constexpr ndensity = 3;
#endif
#ifdef _TSWEVTC
uint constexpr ndensity = 3;
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


// hs
uint constexpr nconstant = 1;
#define HSVAR 0

//functional derivatives = F, B, BFCT, K, he
//primal grid reconstruction stuff- U, dens0, densfct0, edgerecon, recon, edgereconfct, reconfct
//fct stuff- Phi, Mf, edgeflux

#if NTRACERS > 0
uint constexpr nauxiliary = 15;
#else
uint constexpr nauxiliary = 8;
#endif


#define FVAR 0
#define BVAR 1
#define KVAR 2
#define HEVAR 3
#define UVAR 4

#define DENS0VAR 5
#define DENSRECONVAR 6
#define DENSEDGERECONVAR 7

#if NTRACERS > 0
#define BFCTVAR 8
#define DENSFCT0VAR 9
#define DENSFCTRECONVAR 10
#define DENSFCTEDGERECONVAR 11
#define PHIVAR 12
#define EDGEFLUXVAR 13
#define MFVAR 14
#endif

// concentration 0-forms for den

#if NTRACERS > 0
uint constexpr ndiagnostic = 2;
#define DENSFCTLDIAGVAR 1
#else
uint constexpr ndiagnostic = 1;
#endif
#define DENSLDIAGVAR 0


//track total densities, densfct mix/max, energy (P+K+I+T)
#if NTRACERS > 0
uint constexpr nstats = 5;
#else
uint constexpr nstats = 2;
#endif

#define DENSSTAT 0
#define ESTAT 1
#if NTRACERS > 0
#define DENSFCTSTAT 2
#define DENSFCTMAXSTAT 3
#define DENSFCTMINSTAT 4
#endif





// *******   Model Specific Parameters   ***********//


real constexpr Lx = 5000. * 1000.;
real constexpr H0 = 750.0;
real constexpr ox = 0.2;
real constexpr sigmax = 3./60.*Lx;
real constexpr dh = 75.0;
real constexpr g = 9.80616;
real constexpr xc1 = (0.5-ox) * Lx;
real constexpr xc2 = (0.5+ox) * Lx;
real constexpr xc = 0.5 * Lx;
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

// FIX THESE

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

// FIX THESE

  if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
  {
  params.xlen = Lx;
  params.xc = Lx/2.;
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
     realArr dens0var, realArr densfct0var,
     const realArr Vvar, const realArr densvar, const realArr densfctvar) {
#else
void YAKL_INLINE compute_diagnostic_quantities(
  realArr dens0var,
  const realArr Vvar, const realArr densvar) {
#endif

       int is = topology->is;
       int js = topology->js;
       int ks = topology->ks;

       yakl::parallel_for("ComputeDiagI", topology->n_cells, YAKL_LAMBDA (int iGlob) {
         SArray<real,1> zeta;
         real hv;

         int k, j, i;
         yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);


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
   diagnostic_vars.fields_arr[DENSLDIAGVAR].data, diagnostic_vars.fields_arr[DENSFCTLDIAGVAR].data,
   x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data);
#else
compute_diagnostic_quantities(
diagnostic_vars.fields_arr[DENSLDIAGVAR].data,
x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data);
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
    realArr Uvar, realArr dens0var, realArr densfct0var,
    const realArr Vvar, const realArr densvar, const realArr densfctvar) {
#else
void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
  realArr Uvar, realArr dens0var,
  const realArr Vvar, const realArr densvar) {
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

        //compute KE = 0.5 * phiT(u,v)
        compute_phiT(Kvar, Uvar, Vvar, is, js, ks, i, j, k);
        Kvar(0, k+ks, j+js, i+is) *= 0.5;

      });

      }


#if NTRACERS > 0
  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    realArr Bvar, realArr Bfctvar,
    const realArr Kvar, const realArr dens0var, const realArr densfct0var, const realArr HSvar) {
#else
void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
  realArr Bvar,
  const realArr Kvar, const realArr dens0var, const realArr HSvar) {
#endif
      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeDiagIII", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        SArray<real,1> hs0;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

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

//FIX THESE
#ifdef _SWEVT
#endif
#ifdef _TSWEVT
#endif
#ifdef _TSWEVTC
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
  void YAKL_INLINE compute_edge_reconstructions(realArr densedgereconvar, realArr densfctedgereconvar,
    const realArr dens0var, const realArr densfct0var) {
#else
  void YAKL_INLINE compute_edge_reconstructions(realArr densedgereconvar,
    const realArr dens0var) {
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

    });

  }

#if NTRACERS > 0
  void YAKL_INLINE compute_recons(
  realArr densreconvar, realArr densfctreconvar,
  const realArr densedgereconvar, const realArr densfctedgereconvar, const realArr HEvar, const realArr Uvar) {
#else
void YAKL_INLINE compute_recons(
realArr densreconvar,
const realArr densedgereconvar, const realArr HEvar, const realArr Uvar) {
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
  const realArr densreconvar, const realArr densfctreconvar,
  const realArr Bvar, const realArr Bfctvar, const realArr Fvar, const realArr Phivar) {
#else
void YAKL_INLINE compute_tendencies(
realArr denstendvar, realArr Vtendvar,
const realArr densreconvar,
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

//FIX THESE
#ifdef _TSWEVTC
#endif

  });

  }




  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

      //Compute U, dens0, densfct0
#if NTRACERS > 0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data, x.fields_arr[DENSFCTVAR].data);

      this->aux_exchange->exchanges_arr[DENSFCT0VAR].exchange_field(auxiliary_vars.fields_arr[DENSFCT0VAR]);
#else
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data);
#endif

      this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(auxiliary_vars.fields_arr[DENS0VAR]);


      //Compute K, F, he
      compute_functional_derivatives_and_diagnostic_quantities_II(
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      x.fields_arr[VVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data);

      this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
      this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
      this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);

      //Compute B, Bfct
#if NTRACERS > 0
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data, const_vars.fields_arr[HSVAR].data);

      this->aux_exchange->exchanges_arr[BFCTVAR].exchange_field(auxiliary_vars.fields_arr[BFCTVAR]);
#else
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[BVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data, const_vars.fields_arr[HSVAR].data);
#endif

      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);

      // Compute densrecon, densfctrecon
#if NTRACERS > 0
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[DENS0VAR].data, auxiliary_vars.fields_arr[DENSFCT0VAR].data);

      this->aux_exchange->exchanges_arr[DENSFCTEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR]);
#else
compute_edge_reconstructions(
auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
auxiliary_vars.fields_arr[DENS0VAR].data);
#endif

      this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSEDGERECONVAR]);

#if NTRACERS > 0
      compute_recons(
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data,
      auxiliary_vars.fields_arr[DENSEDGERECONVAR].data, auxiliary_vars.fields_arr[DENSFCTEDGERECONVAR].data,
      auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[UVAR].data);
      this->aux_exchange->exchanges_arr[DENSFCTRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSFCTRECONVAR]);
#else
    compute_recons(
    auxiliary_vars.fields_arr[DENSRECONVAR].data,
    auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
    auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[UVAR].data);
#endif
      this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(auxiliary_vars.fields_arr[DENSRECONVAR]);


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
      auxiliary_vars.fields_arr[DENSRECONVAR].data, auxiliary_vars.fields_arr[DENSFCTRECONVAR].data,
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[BFCTVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[PHIVAR].data);
#else
  compute_tendencies(
  xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
  auxiliary_vars.fields_arr[DENSRECONVAR].data,
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
        real KE, PE;
        SArray<real,ndensity> dens0;
        SArray<real,1> h0im1;
        SArray<real,1> U, he;
        SArray<real,1,2> h0arr;

        //compute stats locally
        compute_I<ndensity,diff_ord> (dens0, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i, j, k);
        compute_H<1,diff_ord> (U, progvars.fields_arr[VVAR].data, *this->geom, is, js, ks, i, j, k);

        //Slight hack that only computes I for first dof
        compute_I<1,diff_ord> (h0im1, progvars.fields_arr[DENSVAR].data, *this->geom, is, js, ks, i-1, j, k);

h0arr(0,0) = dens0(0);
h0arr(0,1) = h0im1(0);
phi(he, h0arr);

        KE = 1./2. * (he(0) * ( U(0) * progvars.fields_arr[VVAR].data(0,k+ks,j+js,i+is)));

#ifdef _TSWE
PE = 0.5*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(1) + dens0(1)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif
#ifdef _SWE
PE = PE = 0.5*g*progvars.fields_arr[DENSVAR].data(0,k+ks,j+js,i+is)*dens0(0) + g*dens0(0)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);
#endif

//FIX THESE
#ifdef _SWEVT
#endif
#ifdef _TSWEVT
#endif
#ifdef _TSWEVTC
#endif

#if NTRACERS > 0
for (int l=0;l<ndensityfct;l++)
{PE += 0.5*dens0(0)*progvars.fields_arr[DENSFCTVAR].data(l,k+ks,j+js,i+is);}
#endif

       elocal(0) += PE;
       elocal(1) += KE;
       elocal(2) += KE + PE;

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
  prog_ndofs_arr(VVAR,0) = 1; //v = straight 1-form
  prog_ndofs_arr(DENSVAR,1) = ndensity; //dens = twisted 1-form

  // hs
  const_topo_arr[HSVAR] = &topo;
  const_names_arr[HSVAR] = "hs";
  const_ndofs_arr(HSVAR,1) = 1; //hs = twisted 1-form

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
  aux_ndofs_arr(BVAR,1) = ndensity; //B = straight 0-form
  aux_ndofs_arr(KVAR,1) = 1; //K = twisted 1-form
  aux_ndofs_arr(FVAR,0) = 1; //F = twisted 0-form
  aux_ndofs_arr(UVAR,0) = 1; //U = twisted 0-form
  aux_ndofs_arr(HEVAR,0) = 1; //he lives on edges

  //dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = &topo;
  aux_topo_arr[DENSEDGERECONVAR] = &topo;
  aux_topo_arr[DENS0VAR] = &topo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  aux_ndofs_arr(DENSRECONVAR,0) = ndensity; //densrecon lives on edges
  aux_ndofs_arr(DENSEDGERECONVAR,1) = 2*ndensity; //densedgerecon lives on cells
  aux_ndofs_arr(DENS0VAR,1) = ndensity; //dens0 = straight 0-form

  // concentration 0-forms for dens
  diag_topo_arr[DENSLDIAGVAR] = &topo;
  diag_names_arr[DENSLDIAGVAR] = "densl";
  diag_ndofs_arr(DENSLDIAGVAR,1) = ndensity; //densldiag = straight 0-form

  //densfct stuff- densfct, BFCT, densfct0, edgereconfct, reconfct, Phi, Mf, edgeflux, concentration 0-forms for densfct
#if NTRACERS > 0
  prog_topo_arr[DENSFCTVAR] = &topo;
  prog_names_arr[DENSFCTVAR] = "densfct";
  prog_ndofs_arr(DENSFCTVAR,1) = ndensityfct; //densfct = twisted 1-form

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
  aux_ndofs_arr(BFCTVAR,1) = ndensityfct; //Bfct = straight 0-form
  aux_ndofs_arr(DENSFCTRECONVAR,0) = ndensityfct; //densfctrecon lives on edges
  aux_ndofs_arr(DENSFCTEDGERECONVAR,1) = 2*ndensityfct; //densfctedgerecon lives on cells
  aux_ndofs_arr(DENSFCT0VAR,1) = ndensityfct; //densfct0 = straight 0-form
  aux_ndofs_arr(PHIVAR,0) = ndensityfct; // edges
  aux_ndofs_arr(MFVAR,1) = ndensityfct; //cells
  aux_ndofs_arr(EDGEFLUXVAR,0) = ndensityfct; //edges

  diag_topo_arr[DENSFCTLDIAGVAR] = &topo;
  diag_names_arr[DENSFCTLDIAGVAR] = "densfctl";
  diag_ndofs_arr(DENSFCTLDIAGVAR,1) = ndensityfct; //densfctldiag = straight 0-form
#endif


}

  // *******   Initial Conditions   ***********//

// FIX THIS


real YAKL_INLINE double_vortex_h(real x)
{
  real xprime1 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc1));
  real xprime2 = Lx / (M_PI * sigmax) * sin(M_PI / Lx * (x - xc2));

  return H0 + dh * (exp(-0.5 * (xprime1 * xprime1)) + exp(-0.5 * (xprime2 * xprime2)) - 4. * M_PI * sigmax / Lx);

}

real YAKL_INLINE double_vortex_v(real x) {
  return 0.;
}

real YAKL_INLINE double_vortex_S(real x)
{
    //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. * M_PI / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
    real sval = g * (1. + c * exp(-((x-xc)*(x-xc))/(a*a*D*D)));
    //real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
    //real sval = g;
    //real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35 * Ly && y < 0.65 * Ly ) ? 1. : 0.));
    return sval * double_vortex_h(x);
}

real YAKL_INLINE double_vortex_tracer_square_cent(real x)         {return (x > 0.35*Lx && x < 0.65*Lx                       ) ? 0.005*double_vortex_h(x) : 0.;}
real YAKL_INLINE double_vortex_tracer_square_r(real x)         {return (x > 0.6*Lx && x < 0.9*Lx                        ) ? 0.005*double_vortex_h(x) : 0.;}
real YAKL_INLINE double_vortex_tracer_square_l(real x)         {return (x > 0.1*Lx && x < 0.4*Lx                        ) ? 0.005*double_vortex_h(x) : 0.;}
real YAKL_INLINE double_vortex_tracer_square_rpl(real x)         {return double_vortex_tracer_square_r(x) + double_vortex_tracer_square_l(x);}


real YAKL_INLINE double_vortex_tracer_gaussian(real x)     {return  0.005 *double_vortex_h(x) * exp(-((x-xc)*(x-xc))/(a*a*D*D));}
//{ return 0.005 *double_vortex_h(x,y) * exp(-100. * pow((x-xc)/Lx,2.)) * exp(-100. * pow((y-yc)/Ly,2.)); }


#define vortex1_1d(x)   (0.005 * double_vortex_h(x) *  exp(-100. * pow((x-0.75*Lx)/Lx,2.)))
#define vortex2_1d(x)   (0.005 * double_vortex_h(x) * exp(-50.  * pow((x-0.25*Lx)/Lx,2.)))
real YAKL_INLINE double_vortex_tracer_vortices(real x)         { return vortex1_1d(x)   + vortex2_1d(x); }

// ADD MORE ICs HERE!!!

//wavespeed = sqrt(g * H0)
//dt = Constant(get_dt(wavespeed, cval, order, variant, Lx, Ly, nx, ny))
template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, Geometry<1, nquadx, nquady, nquadz> &geom)
{

// FIX THIS
    if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
    {
        geom.set_primal_1form_values(double_vortex_h, progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
        geom.set_primal_1form_values(double_vortex_S, progvars.fields_arr[DENSVAR], 1);
#endif
#if NTRACERS>0
        geom.set_primal_1form_values(double_vortex_tracer_gaussian, progvars.fields_arr[DENSFCTVAR], 0);
#endif
#if NTRACERS>1
        geom.set_primal_1form_values(double_vortex_tracer_square_cent, progvars.fields_arr[DENSFCTVAR], 1);
#endif
#if NTRACERS>2
        geom.set_primal_1form_values(double_vortex_tracer_square_rpl, progvars.fields_arr[DENSFCTVAR], 2);
#endif
        geom.set_dual_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0);
    }

//FIX THESE
#ifdef _SWEVT
#endif
#ifdef _TSWEVT
#endif
#ifdef _TSWEVTC
#endif

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
