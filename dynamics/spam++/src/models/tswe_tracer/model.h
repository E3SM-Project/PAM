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



// Number of variables

uint constexpr nprognostic = 4; // h, v, S, tr
uint constexpr nconstant = 2;   // hs, coriolis
uint constexpr nauxiliary = 25; // B, F, hv, hrecon, qrecon, FT, U, coriolisrecon, hedegerecon, qedgerecon, h0, q0, K, he, T, srecon, sedgerecon, S0, trrecon, tr0, tredgerecon, Phi, edgeflux, Mf, Xi
uint constexpr nstats = 9;      // M, B, PE, KE, TE, PV, TR, TRmin, TRmax
uint constexpr ndiagnostic = 3;      // qdiag, sldiag, trdiag

#define HVAR 0
#define VVAR 1
#define SVAR 2
#define TRVAR 3

#define HSVAR 0
#define CORIOLISVAR 1

#define UVAR 0
#define H0VAR 1
#define HVVAR 2
#define KVAR 3
#define FVAR 4
#define HEVAR 5
#define BVAR 6
#define FTVAR 7
#define Q0VAR 8
#define HRECONVAR 9
#define HEDGERECONVAR 10
#define QRECONVAR 11
#define QEDGERECONVAR 12
#define CORIOLISRECONVAR 13
#define SRECONVAR 14
#define SEDGERECONVAR 15
#define S0VAR 16
#define TVAR 17
#define TR0VAR 18
#define TRRECONVAR 19
#define TREDGERECONVAR 20
#define PHIVAR 21
#define EDGEFLUXVAR 22
#define MFVAR 23
#define XIVAR 24

#define QDIAGVAR 0
#define SLDIAGVAR 1
#define TRDIAGVAR 2

#define MSTAT 0
#define PESTAT 1
#define KESTAT 2
#define TESTAT 3
#define PVSTAT 4
#define BSTAT 5
#define TRSTAT 6
#define TRMINSTAT 7
#define TRMAXSTAT 8

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


   void YAKL_INLINE compute_diagnostic_quantities(
     realArr Q0var, realArr SLvar, realArr TR0var,
     const realArr Vvar, const realArr Hvar, const realArr Svar, const realArr TRvar, const realArr coriolisvar) {

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
   compute_R<0> (hv, Hvar, is, js, ks, i, j, k);
   Q0var(0, k+ks, j+js, i+is) = (zeta(0) + coriolisvar(0, k+ks, j+js, i+is)) / hv;

// compute SL
SLvar(0, k+ks, j+js, i+is) = Svar(0, k+ks, j+js, i+is) / Hvar(0, k+ks, j+js, i+is);

// compute tr0
TR0var(0, k+ks, j+js, i+is) = TRvar(0, k+ks, j+js, i+is) / Hvar(0, k+ks, j+js, i+is);
       });

     }



   void compute_diag(const VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<ndiag> &diagnostic_vars)
   {


   compute_diagnostic_quantities(
   diagnostic_vars.fields_arr[QDIAGVAR].data, diagnostic_vars.fields_arr[SLDIAGVAR].data, diagnostic_vars.fields_arr[TRDIAGVAR].data,
   x.fields_arr[VVAR].data, x.fields_arr[HVAR].data, x.fields_arr[SVAR].data, x.fields_arr[TRVAR].data, const_vars.fields_arr[CORIOLISVAR].data);
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

    this->is_initialized = true;
  }





  void YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_I(
    realArr Uvar, realArr HVvar, realArr Q0var, realArr H0var, realArr S0var, realArr TR0var,
    const realArr Vvar, const realArr Hvar, const realArr Svar, const realArr TRvar) {

      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeDiagI", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

        // compute h0 = I h, U = H v, S0 = I S, tr0 = I TR
        compute_H<1, diff_ord>(Uvar, Vvar, *this->geom, is, js, ks, i, j, k);
        compute_I<1, diff_ord>(H0var, Hvar, *this->geom, is, js, ks, i, j, k);
        compute_I<1, diff_ord>(S0var, Svar, *this->geom, is, js, ks, i, j, k);
        compute_I<1, diff_ord>(TR0var, TRvar, *this->geom, is, js, ks, i, j, k);

// compute zeta = D2 v
compute_D2<1>(Q0var, Vvar, is, js, ks, i, j, k);

// compute q0 = zeta / R h
  compute_R<0> (HVvar, Hvar, is, js, ks, i, j, k);
  Q0var(0, k+ks, j+js, i+is) = Q0var(0, k+ks, j+js, i+is) / HVvar(0, k+ks, j+js, i+is);

      });

    }

    void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_II(
      realArr Fvar, realArr Kvar, realArr HEvar, const realArr Vvar, const realArr Uvar, const realArr H0var) {

        int is = topology->is;
        int js = topology->js;
        int ks = topology->ks;

        yakl::parallel_for("ComputeDiagII", topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

        // compute he = phi * h0
        compute_phi<0>(HEvar, H0var, is, js, ks, i, j, k);

        //compute F = he * U
        Fvar(0, k+ks, j+js, i+is) = Uvar(0, k+ks, j+js, i+is) * HEvar(0, k+ks, j+js, i+is);
        Fvar(1, k+ks, j+js, i+is) = Uvar(1, k+ks, j+js, i+is) * HEvar(1, k+ks, j+js, i+is);

        //compute KE = 0.5 * phiT(u,v)
        compute_phiT(Kvar, Uvar, Vvar, is, js, ks, i, j, k);
        Kvar(0, k+ks, j+js, i+is) *= 0.5;

      });

      }

  void  YAKL_INLINE compute_functional_derivatives_and_diagnostic_quantities_III(
    realArr FTvar, realArr Q0var, realArr Bvar, realArr Tvar, realArr Xivar,
    const realArr Fvar, const realArr Uvar, const realArr HVvar,
    const realArr Kvar, const realArr H0var, const realArr S0var, const realArr TR0var, const realArr HSvar) {

      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeDiagIII", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

      compute_W(FTvar, Fvar, is, js, ks, i, j, k);

      //Compute T = Ihs + h0/2;
      compute_I<1, diff_ord>(Tvar, HSvar, *this->geom, is, js, ks, i, j, k);
      Tvar(0, k+ks, j+js, i+is) += H0var(0, k+ks, j+js, i+is)/2.;

      // Compute B = IK + S0
      compute_I<1, diff_ord>(Bvar, Kvar, *this->geom, is, js, ks, i, j, k);
      Bvar(0, k+ks, j+js, i+is) += S0var(0, k+ks, j+js, i+is)/2.;
      Bvar(0, k+ks, j+js, i+is) += TR0var(0, k+ks, j+js, i+is)/2.;

      // Compute Xi = h0/2;
      Xivar(0, k+ks, j+js, i+is) = H0var(0, k+ks, j+js, i+is)/2.;

    });

    }

  void YAKL_INLINE compute_edge_reconstructions(realArr Hedgereconvar, realArr Sedgereconvar, realArr TRedgereconvar, realArr Qedgereconvar,
    const realArr H0var, const realArr S0var, const realArr TR0var, const realArr Q0var) {

    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

    yakl::parallel_for("ComputeEdgeRecon", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

      compute_primal_edge_recon<1, reconstruction_type, reconstruction_order>(Hedgereconvar, H0var, is, js, ks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
      compute_primal_edge_recon<1, reconstruction_type, reconstruction_order>(Sedgereconvar, S0var, is, js, ks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
      compute_primal_edge_recon<1, reconstruction_type, reconstruction_order>(TRedgereconvar, TR0var, is, js, ks, i, j, k, primal_wenoRecon, primal_to_gll, primal_wenoIdl, primal_wenoSigma);
      compute_dual_edge_recon<1, dual_reconstruction_type, dual_reconstruction_order>(Qedgereconvar, Q0var, is, js, ks, i, j, k, dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);

    });

  }

  void YAKL_INLINE compute_recons(
  realArr Hreconvar, realArr Sreconvar, realArr TRreconvar, realArr Qreconvar, realArr Coriolisreconvar,
  const realArr Hedgereconvar, const realArr Sedgereconvar, const realArr TRedgereconvar, const realArr Qedgereconvar, const realArr HEvar, const realArr HVvar, const realArr Coriolisvar,
  const realArr FTvar, const realArr Uvar) {

    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

    yakl::parallel_for("ComputeRecon", topology->n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

      compute_primal_recon<1, reconstruction_type>(Hreconvar, Hedgereconvar, Uvar, is, js, ks, i, j, k);
      compute_primal_recon<1, reconstruction_type>(Sreconvar, Sedgereconvar, Uvar, is, js, ks, i, j, k);
      compute_primal_recon<1, reconstruction_type>(TRreconvar, TRedgereconvar, Uvar, is, js, ks, i, j, k);
      compute_dual_recon<1, dual_reconstruction_type>(Qreconvar, Qedgereconvar, FTvar, is, js, ks, i, j, k);

    //scale primal recon
    Hreconvar(0,k+ks,j+js,i+is) = Hreconvar(0,k+ks,j+js,i+is) / HEvar(0,k+ks,j+js,i+is);
    Hreconvar(1,k+ks,j+js,i+is) = Hreconvar(1,k+ks,j+js,i+is) / HEvar(1,k+ks,j+js,i+is);
    Sreconvar(0,k+ks,j+js,i+is) = Sreconvar(0,k+ks,j+js,i+is) / HEvar(0,k+ks,j+js,i+is);
    Sreconvar(1,k+ks,j+js,i+is) = Sreconvar(1,k+ks,j+js,i+is) / HEvar(1,k+ks,j+js,i+is);
    TRreconvar(0,k+ks,j+js,i+is) = TRreconvar(0,k+ks,j+js,i+is) / HEvar(0,k+ks,j+js,i+is);
    TRreconvar(1,k+ks,j+js,i+is) = TRreconvar(1,k+ks,j+js,i+is) / HEvar(1,k+ks,j+js,i+is);

    //coriolis recon
    //x-dir
    Coriolisreconvar(1, k+ks, j+js, i+is) = 0.5*(Coriolisvar(0,k+ks,j+js,i+is)/HVvar(0,k+ks,j+js,i+is) + Coriolisvar(0,k+ks,j+js,i+is+1) / HVvar(0,k+ks,j+js,i+is+1));
    //y-dir
    Coriolisreconvar(0, k+ks, j+js, i+is) = 0.5*(Coriolisvar(0,k+ks,j+js,i+is)/HVvar(0,k+ks,j+js,i+is) + Coriolisvar(0,k+ks,j+js+1,i+is) / HVvar(0,k+ks,j+js+1,i+is));
    });

}



  void YAKL_INLINE compute_tendencies(
  realArr Htendvar, realArr Stendvar, realArr TRtendvar, realArr Vtendvar,
  const realArr Hreconvar, const realArr Sreconvar, const realArr TRreconvar, const realArr Qreconvar, const realArr Coriolisreconvar,
  const realArr Bvar, const realArr Tvar, const realArr Xivar, const realArr Fvar, const realArr Phivar) {

    int is = topology->is;
    int js = topology->js;
    int ks = topology->ks;

      yakl::parallel_for("ComputeTendencies", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

    compute_wDbar2<1> (Htendvar, Hreconvar, Fvar, is, js, ks, i, j, k);
    compute_wDbar2<1> (Stendvar, Sreconvar, Fvar, is, js, ks, i, j, k);
    compute_wDbar2_fct<1> (TRtendvar, TRreconvar, Phivar, Fvar, is, js, ks, i, j, k);

    compute_wD1<1> (Vtendvar, Hreconvar, Bvar, is, js, ks, i, j, k);
    compute_wD1<1, ADD_MODE::ADD> (Vtendvar, Sreconvar, Tvar, is, js, ks, i, j, k);
    compute_wD1_fct<1, ADD_MODE::ADD> (Vtendvar, TRreconvar, Phivar, Xivar, is, js, ks, i, j, k);

    if (qf_choice == QF_MODE::EC)
    { compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, is, js, ks, i, j, k);}
    if (qf_choice == QF_MODE::NOEC)
    { compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, is, js, ks, i, j, k);}
    compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Coriolisreconvar, Fvar, is, js, ks, i, j, k);

  });

  }




  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconst> &const_vars, VariableSet<nprog> &x, VariableSet<naux> &auxiliary_vars, VariableSet<nprog> &xtend)
  {

      //Compute U, q0, h0, S0, TR0
      compute_functional_derivatives_and_diagnostic_quantities_I(
      auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[HVVAR].data, auxiliary_vars.fields_arr[Q0VAR].data,
      auxiliary_vars.fields_arr[H0VAR].data, auxiliary_vars.fields_arr[S0VAR].data, auxiliary_vars.fields_arr[TR0VAR].data,
      x.fields_arr[VVAR].data, x.fields_arr[HVAR].data, x.fields_arr[SVAR].data, x.fields_arr[TRVAR].data);


      this->aux_exchange->exchanges_arr[UVAR].exchange_field(auxiliary_vars.fields_arr[UVAR]);
      this->aux_exchange->exchanges_arr[H0VAR].exchange_field(auxiliary_vars.fields_arr[H0VAR]);
      this->aux_exchange->exchanges_arr[S0VAR].exchange_field(auxiliary_vars.fields_arr[S0VAR]);
      this->aux_exchange->exchanges_arr[TR0VAR].exchange_field(auxiliary_vars.fields_arr[TR0VAR]);
      this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(auxiliary_vars.fields_arr[Q0VAR]);
      this->aux_exchange->exchanges_arr[HVVAR].exchange_field(auxiliary_vars.fields_arr[HVVAR]);


      //Compute K, F, he
      compute_functional_derivatives_and_diagnostic_quantities_II(
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[HEVAR].data,
      x.fields_arr[VVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[H0VAR].data);

      this->aux_exchange->exchanges_arr[FVAR].exchange_field(auxiliary_vars.fields_arr[FVAR]);
      this->aux_exchange->exchanges_arr[KVAR].exchange_field(auxiliary_vars.fields_arr[KVAR]);
      this->aux_exchange->exchanges_arr[HEVAR].exchange_field(auxiliary_vars.fields_arr[HEVAR]);

      //Compute FT, B, T
      compute_functional_derivatives_and_diagnostic_quantities_III(
      auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[Q0VAR].data, auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[TVAR].data, auxiliary_vars.fields_arr[XIVAR].data,
      auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[UVAR].data, auxiliary_vars.fields_arr[HVVAR].data,
      auxiliary_vars.fields_arr[KVAR].data, auxiliary_vars.fields_arr[H0VAR].data, auxiliary_vars.fields_arr[S0VAR].data, auxiliary_vars.fields_arr[TR0VAR].data, const_vars.fields_arr[HSVAR].data);

      this->aux_exchange->exchanges_arr[FTVAR].exchange_field(auxiliary_vars.fields_arr[FTVAR]);
      this->aux_exchange->exchanges_arr[BVAR].exchange_field(auxiliary_vars.fields_arr[BVAR]);
      this->aux_exchange->exchanges_arr[TVAR].exchange_field(auxiliary_vars.fields_arr[TVAR]);
      this->aux_exchange->exchanges_arr[XIVAR].exchange_field(auxiliary_vars.fields_arr[XIVAR]);

      // Compute hrecon, srecon, qrecon and frecon
      compute_edge_reconstructions(
      auxiliary_vars.fields_arr[HEDGERECONVAR].data, auxiliary_vars.fields_arr[SEDGERECONVAR].data, auxiliary_vars.fields_arr[TREDGERECONVAR].data, auxiliary_vars.fields_arr[QEDGERECONVAR].data,
      auxiliary_vars.fields_arr[H0VAR].data, auxiliary_vars.fields_arr[S0VAR].data, auxiliary_vars.fields_arr[TR0VAR].data, auxiliary_vars.fields_arr[Q0VAR].data);

      this->aux_exchange->exchanges_arr[HEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[HEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[SEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[SEDGERECONVAR]);
      this->aux_exchange->exchanges_arr[TREDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[TREDGERECONVAR]);
      this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(auxiliary_vars.fields_arr[QEDGERECONVAR]);

      compute_recons(
      auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[TRRECONVAR].data,
      auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[HEDGERECONVAR].data, auxiliary_vars.fields_arr[SEDGERECONVAR].data, auxiliary_vars.fields_arr[TREDGERECONVAR].data,
      auxiliary_vars.fields_arr[QEDGERECONVAR].data, auxiliary_vars.fields_arr[HEVAR].data, auxiliary_vars.fields_arr[HVVAR].data,
      const_vars.fields_arr[CORIOLISVAR].data, auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[UVAR].data);

      this->aux_exchange->exchanges_arr[HRECONVAR].exchange_field(auxiliary_vars.fields_arr[HRECONVAR]);
      this->aux_exchange->exchanges_arr[SRECONVAR].exchange_field(auxiliary_vars.fields_arr[SRECONVAR]);
      this->aux_exchange->exchanges_arr[TRRECONVAR].exchange_field(auxiliary_vars.fields_arr[TRRECONVAR]);
      this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(auxiliary_vars.fields_arr[QRECONVAR]);
      this->aux_exchange->exchanges_arr[CORIOLISRECONVAR].exchange_field(auxiliary_vars.fields_arr[CORIOLISRECONVAR]);


//Compute fct
int is = topology->is;
int js = topology->js;
int ks = topology->ks;

yakl::parallel_for("ComputeEdgeFlux", topology->n_cells, YAKL_LAMBDA (int iGlob) {
  int k, j, i;
  yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_edgefluxes<1> (auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[TRRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, is, js, ks, i, j, k);
});
this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(auxiliary_vars.fields_arr[EDGEFLUXVAR]);


yakl::parallel_for("ComputeMf", topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_Mf<1> (auxiliary_vars.fields_arr[MFVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, is, js, ks, i, j, k);
});

this->aux_exchange->exchanges_arr[MFVAR].exchange_field(auxiliary_vars.fields_arr[MFVAR]);

yakl::parallel_for("ComputePhi", topology->n_cells, YAKL_LAMBDA (int iGlob) {
int k, j, i;
yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
compute_Phi<1> (auxiliary_vars.fields_arr[PHIVAR].data, auxiliary_vars.fields_arr[EDGEFLUXVAR].data, auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[TRVAR].data, is, js, ks, i, j, k);
});


this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(auxiliary_vars.fields_arr[PHIVAR]);

      // Compute tendencies
      compute_tendencies(
      xtend.fields_arr[HVAR].data, xtend.fields_arr[SVAR].data, xtend.fields_arr[TRVAR].data, xtend.fields_arr[VVAR].data,
      auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[TRRECONVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
      auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[TVAR].data, auxiliary_vars.fields_arr[XIVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[PHIVAR].data);
};

};
// *******   Statistics Calculations   ***********//







// THIS STUFF SHOULD BE CLEANED UP AND GENERALIZED LIKE VARIABLE SETS IF POSSIBLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE!
class Stat
{
public:
  realArr data;
  real local_dat;
  real global_dat;
  std::string name;


void initialize(std::string statName, ModelParameters &params, Parallel &par)
{
  name = statName;

  if (par.masterproc)
  {
    data = realArr(name.c_str(), params.Nsteps/params.Nstat + 1);
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
    stats_arr[MSTAT].initialize("mass", params, par);
    stats_arr[PESTAT].initialize("potential_energy", params, par);
    stats_arr[KESTAT].initialize("kinetic_energy", params, par);
    stats_arr[TESTAT].initialize("total_energy", params, par);
    stats_arr[BSTAT].initialize("bouyancy", params, par);
    stats_arr[TRSTAT].initialize("tracer", params, par);
    stats_arr[TRMINSTAT].initialize("tracermin", params, par);
    stats_arr[TRMAXSTAT].initialize("tracermax", params, par);
    stats_arr[PVSTAT].initialize("potential_vorticity", params, par);
    masterproc = par.masterproc;
  }



  void compute( VariableSet<nprog> &progvars,  VariableSet<nconst> &constvars, int i)
  {

      this->stats_arr[MSTAT].local_dat = 0;
      this->stats_arr[PESTAT].local_dat = 0;
      this->stats_arr[KESTAT].local_dat = 0;
      this->stats_arr[TESTAT].local_dat = 0;
      this->stats_arr[BSTAT].local_dat = 0;
      this->stats_arr[TRSTAT].local_dat = 0;
      this->stats_arr[TRMINSTAT].local_dat = 0;
      this->stats_arr[TRMAXSTAT].local_dat = 0;
      this->stats_arr[PVSTAT].local_dat = 0;

      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      yakl::parallel_for("ComputeStats", topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);
        real eta, KE, PE;
        SArray<real,1> zeta, h0, S0, h0im1, h0jm1;
        SArray<real,2> U, he;
        SArray<real,2,2> h0arr;

        //compute stats locally

        compute_I<1,diff_ord> (h0, progvars.fields_arr[HVAR].data, *this->geom, is, js, ks, i, j, k);
        compute_I<1,diff_ord> (S0, progvars.fields_arr[SVAR].data, *this->geom, is, js, ks, i, j, k);
        compute_I<1,diff_ord> (h0im1, progvars.fields_arr[HVAR].data, *this->geom, is, js, ks, i-1, j, k);
        compute_I<1,diff_ord> (h0jm1, progvars.fields_arr[HVAR].data, *this->geom, is, js, ks, i, j-1, k);
        compute_H<1,diff_ord> (U, progvars.fields_arr[VVAR].data, *this->geom, is, js, ks, i, j, k);

h0arr(0,0) = h0(0);
h0arr(1,0) = h0(0);
h0arr(0,1) = h0im1(0);
h0arr(1,1) = h0jm1(0);
phi(he, h0arr);

        KE = 1./2. * (he(0) * ( U(0) * progvars.fields_arr[VVAR].data(0,k+ks,j+js,i+is)) +
                    + he(1) * ( U(1) * progvars.fields_arr[VVAR].data(1,k+ks,j+js,i+is)));

      compute_D2<1>(zeta, progvars.fields_arr[VVAR].data, is, js, ks, i, j, k);
      eta = zeta(0) + constvars.fields_arr[CORIOLISVAR].data(0,k+ks,j+js,i+is);

        PE = 0.5*progvars.fields_arr[HVAR].data(0,k+ks,j+js,i+is)*S0(0) + S0(0)*constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is) + 0.5*h0(0)*progvars.fields_arr[TRVAR].data(0,k+ks,j+js,i+is);

       this->stats_arr[MSTAT].local_dat += progvars.fields_arr[HVAR].data(0,k+ks,j+js,i+is);
       this->stats_arr[BSTAT].local_dat += progvars.fields_arr[SVAR].data(0,k+ks,j+js,i+is);
       this->stats_arr[TRSTAT].local_dat += progvars.fields_arr[TRVAR].data(0,k+ks,j+js,i+is);
       this->stats_arr[PESTAT].local_dat += PE;
       this->stats_arr[KESTAT].local_dat += KE;
       this->stats_arr[TESTAT].local_dat += KE + PE;
       this->stats_arr[PVSTAT].local_dat += eta;

    });

    this->stats_arr[TRMINSTAT].local_dat = progvars.fields_arr[TRVAR].min();
    this->stats_arr[TRMAXSTAT].local_dat = progvars.fields_arr[TRVAR].max();

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &this->stats_arr[MSTAT].local_dat, &this->stats_arr[MSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[MSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[PESTAT].local_dat, &this->stats_arr[PESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[KESTAT].local_dat, &this->stats_arr[KESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[KESTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[TESTAT].local_dat, &this->stats_arr[TESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[TESTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[PVSTAT].local_dat, &this->stats_arr[PVSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[BSTAT].local_dat, &this->stats_arr[BSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[BSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[TRSTAT].local_dat, &this->stats_arr[TRSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[TRSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[TRMINSTAT].local_dat, &this->stats_arr[TRMINSTAT].global_dat, 1, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[TRMINSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[TRMAXSTAT].local_dat, &this->stats_arr[TRMAXSTAT].global_dat, 1, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[TRMAXSTAT]);

    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
  this->stats_arr[MSTAT].data(i) = this->stats_arr[MSTAT].global_dat;
  this->stats_arr[PESTAT].data(i) = this->stats_arr[PESTAT].global_dat;
  this->stats_arr[KESTAT].data(i) = this->stats_arr[KESTAT].global_dat;
  this->stats_arr[TESTAT].data(i) = this->stats_arr[TESTAT].global_dat;
  this->stats_arr[BSTAT].data(i) = this->stats_arr[BSTAT].global_dat;
  this->stats_arr[TRSTAT].data(i) = this->stats_arr[TRSTAT].global_dat;
  this->stats_arr[TRMINSTAT].data(i) = this->stats_arr[TRMINSTAT].global_dat;
  this->stats_arr[TRMAXSTAT].data(i) = this->stats_arr[TRMAXSTAT].global_dat;
  this->stats_arr[PVSTAT].data(i) = this->stats_arr[PVSTAT].global_dat;
  }
  }
};




// *******   VariableSet Initialization   ***********//
// h, v
// hs
// B, F, q

template <uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology &topo,
SArray<int, nprog, 4> &prog_ndofs_arr, SArray<int, nconst, 4> &const_ndofs_arr, SArray<int, naux, 4> &aux_ndofs_arr, SArray<int, ndiag, 4> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology *, nprog> &prog_topo_arr, std::array<const Topology *, nconst> &const_topo_arr, std::array<const Topology *, naux> &aux_topo_arr, std::array<const Topology *, ndiag> &diag_topo_arr)
{
  prog_topo_arr[HVAR] = &topo;
  prog_topo_arr[VVAR] = &topo;
  prog_topo_arr[SVAR] = &topo;
  prog_topo_arr[TRVAR] = &topo;

  const_topo_arr[HSVAR] = &topo;
  const_topo_arr[CORIOLISVAR] = &topo;

  aux_topo_arr[BVAR] = &topo;
  aux_topo_arr[FVAR] = &topo;
  aux_topo_arr[UVAR] = &topo;
  aux_topo_arr[FTVAR] = &topo;
  aux_topo_arr[HRECONVAR] = &topo;
  aux_topo_arr[HEDGERECONVAR] = &topo;
  aux_topo_arr[CORIOLISRECONVAR] = &topo;
  aux_topo_arr[HVVAR] = &topo;
  aux_topo_arr[Q0VAR] = &topo;
  aux_topo_arr[H0VAR] = &topo;
  aux_topo_arr[QRECONVAR] = &topo;
  aux_topo_arr[QEDGERECONVAR] = &topo;
  aux_topo_arr[KVAR] = &topo;
  aux_topo_arr[HEVAR] = &topo;
  aux_topo_arr[TVAR] = &topo;
  aux_topo_arr[XIVAR] = &topo;
  aux_topo_arr[SRECONVAR] = &topo;
  aux_topo_arr[SEDGERECONVAR] = &topo;
  aux_topo_arr[S0VAR] = &topo;
  aux_topo_arr[TRRECONVAR] = &topo;
  aux_topo_arr[TREDGERECONVAR] = &topo;
  aux_topo_arr[TR0VAR] = &topo;
  aux_topo_arr[PHIVAR] = &topo;
  aux_topo_arr[MFVAR] = &topo;
  aux_topo_arr[EDGEFLUXVAR] = &topo;

  diag_topo_arr[QDIAGVAR] = &topo;
  diag_topo_arr[SLDIAGVAR] = &topo;
  diag_topo_arr[TRDIAGVAR] = &topo;


  prog_names_arr[HVAR] = "h";
  prog_names_arr[VVAR] = "v";
  prog_names_arr[SVAR] = "S";
  prog_names_arr[TRVAR] = "tr";

  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";

  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[HRECONVAR] = "hrecon";
  aux_names_arr[HEDGERECONVAR] = "hedgerecon";
  aux_names_arr[KVAR] = "K";
  aux_names_arr[HEVAR] = "he";
  aux_names_arr[CORIOLISRECONVAR] = "coriolisrecon";
  aux_names_arr[HVVAR] = "hv";
  aux_names_arr[Q0VAR] = "q";
  aux_names_arr[H0VAR] = "h0";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QEDGERECONVAR] = "qedgerecon";
  aux_names_arr[TVAR] = "T";
  aux_names_arr[XIVAR] = "Xi";
  aux_names_arr[SEDGERECONVAR] = "Sedgerecon";
  aux_names_arr[SRECONVAR] = "Srecon";
  aux_names_arr[S0VAR] = "S0";
  aux_names_arr[TR0VAR] = "tr0";
  aux_names_arr[TRRECONVAR] = "trrecon";
  aux_names_arr[TREDGERECONVAR] = "tredgerecon";
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";


  diag_names_arr[QDIAGVAR] = "q";
  diag_names_arr[SLDIAGVAR] = "sl";
  diag_names_arr[TRDIAGVAR] = "tr0";




//primal grid represents twisted quantities, dual grid straight quantities
    prog_ndofs_arr(HVAR,2) = 1; //h = twisted 2-form
    prog_ndofs_arr(VVAR,1) = 1; //v = straight 1-form
    prog_ndofs_arr(SVAR,2) = 1; //S = twisted 2-form
    prog_ndofs_arr(TRVAR,2) = 1; //tr = twisted 2-form

    const_ndofs_arr(HSVAR,2) = 1; //hs = twisted 2-form
    const_ndofs_arr(CORIOLISVAR,0) = 1; //f = straight 2-form

    aux_ndofs_arr(BVAR,2) = 1; //B = straight 0-form
    aux_ndofs_arr(KVAR,2) = 1; //K = twisted 2-form
    aux_ndofs_arr(FVAR,1) = 1; //F = twisted 1-form
    aux_ndofs_arr(UVAR,1) = 1; //U = twisted 1-form
    aux_ndofs_arr(FTVAR,1) = 1; //FT = straight 1-form
    aux_ndofs_arr(HRECONVAR,1) = 1; //hrecon lives on edges
    aux_ndofs_arr(HEDGERECONVAR,2) = 4; //hedgerecon lives on cells
    aux_ndofs_arr(HEVAR,1) = 1; //he lives on edges
    aux_ndofs_arr(CORIOLISRECONVAR,1) = 1; //coriolisrecon lives on edges
    aux_ndofs_arr(HVVAR,0) = 1; //hv = straight 2-form
    aux_ndofs_arr(H0VAR,2) = 1; //h0 = straight 0-form
    aux_ndofs_arr(Q0VAR,0) = 1; //q0 = twisted 0-form
    aux_ndofs_arr(QRECONVAR,1) = 1; //qrecon lives on edges
    aux_ndofs_arr(QEDGERECONVAR,0) = 4; //qedgerecon lives on dual cells
    aux_ndofs_arr(SEDGERECONVAR,2) = 4; //sedgerecon lives on cells
    aux_ndofs_arr(SRECONVAR,1) = 1; //srecon lives on edges
    aux_ndofs_arr(TVAR,2) = 1; //T = straight 0-form
    aux_ndofs_arr(XIVAR,2) = 1; //Xi = straight 0-form
    aux_ndofs_arr(S0VAR,2) = 1; //S0 = straight 0-form
    aux_ndofs_arr(TR0VAR,2) = 1; //tr0 = straight 0-form
    aux_ndofs_arr(TRRECONVAR,1) = 1; //trrecon lives on edges
    aux_ndofs_arr(TREDGERECONVAR,2) = 4; //tredgerecon lives on cells
    aux_ndofs_arr(PHIVAR,1) = 1;
    aux_ndofs_arr(MFVAR,2) = 1;
    aux_ndofs_arr(EDGEFLUXVAR,1) = 1;


    diag_ndofs_arr(QDIAGVAR,0) = 1; //qdiag = twisted 0-form
    diag_ndofs_arr(SLDIAGVAR,2) = 1; //sl = straight 0-form
    diag_ndofs_arr(TRDIAGVAR,2) = 1; //tr0 = straight 0-form


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

real YAKL_INLINE double_vortex_tracer(real x, real y)         {return (x > 0.35*Lx && x < 0.65*Lx && y > 0.35*Ly && y < 0.65*Ly                        ) ? 0.005*double_vortex_h(x,y) : 0.;}

//wavespeed = sqrt(g * H0)
//dt = Constant(get_dt(wavespeed, cval, order, variant, Lx, Ly, nx, ny))


template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<nprog> &progvars, VariableSet<nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{

    if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
    {
        geom.set_primal_2form_values(double_vortex_h, progvars.fields_arr[HVAR], 0);
        geom.set_primal_2form_values(double_vortex_S, progvars.fields_arr[SVAR], 0);
        geom.set_primal_2form_values(double_vortex_tracer, progvars.fields_arr[TRVAR], 0);
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
