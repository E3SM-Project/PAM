#ifndef _MODEL_H_
#define _MODEL_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "divergence.h"
#include "gradient.h"
#include "curl.h"
#include "weno.h"
#include "weno_dual.h"
#include "cfv.h"
#include "cfv_dual.h"
#include "Q2D.h"
#include "W2D.h"
#include "geometry.h"
#include "params.h"
#include "string.h"



// Number of variables
// THIS REALLY NEEDS TO CHANGE DEPENDING ON 1D/2D!

uint constexpr nprognostic = 3; // h, v, S
uint constexpr nconstant = 2;   // hs, coriolis
uint constexpr nauxiliary = 9; // B, F, T, q, sl, hrecon, qrecon, srecon, FT
uint constexpr nstats = 6;      // M, PE, KE, TE, B, PV
uint constexpr ndiagnostic = 2;      // q0, sl0

#define HVAR 0
#define VVAR 1
#define SVAR 2

// THIS REALLY NEEDS TO CHANGE DEPENDING ON 1D/2D!
#define HSVAR 0
#define CORIOLISVAR 1

// THIS REALLY NEEDS TO CHANGE DEPENDING ON 1D/2D!
#define BVAR 0
#define FVAR 1
#define TVAR 2
#define HRECONVAR 3
#define SRECONVAR 4
#define SLVAR 5
#define QVAR 6
#define QRECONVAR 7
#define FTVAR 8

// THIS REALLY NEEDS TO CHANGE DEPENDING ON 1D/2D!
#define Q0VAR 0
#define SL0VAR 1

#define MSTAT 0
#define PESTAT 1
#define KESTAT 2
#define TESTAT 3
#define BSTAT 4
// THIS REALLY NEEDS TO CHANGE DEPENDING ON 1D/2D!
#define PVSTAT 5

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

template<uint ndims> void set_model_specific_params(std::string inFile, ModelParameters &params)
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

template <uint ndims, uint nprog, uint nconst, uint ndiag> class Diagnostics {
public:

  const Topology<ndims> *topology;
  Geometry<ndims,1,1,1> *geom;

  bool is_initialized;

   Diagnostics() {
     this->is_initialized = false;
     std::cout << "CREATED DIAGNOSTICS\n";
   }

   void initialize(const Topology<ndims> &topo, Geometry<ndims,1,1,1> &geom)
   {
     this->topology = &topo;
     this->geom = &geom;
     this->is_initialized = true;
   }

   void YAKL_INLINE compute_diagnostic_quantities(realArr q, realArr sl, const realArr v, const realArr h, const realArr S, const realArr coriolis, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom)
   {

       int is = topology.is;
       int js = topology.js;
       int ks = topology.ks;

       real he0, he1, KE, zeta, eta, hv, U, V;

         yakl::parallel_for("ComputeDiagnostics", topology.n_cells, YAKL_LAMBDA (int iGlob) {
           int k, j, i;
           yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);

        sl(0,k+ks,j+js,i+is) = S(0,k+ks,j+js,i+is) / h(0,k+ks,j+js,i+is);

        if (ndims == 2) {
       zeta = (v(1,k+ks,j+js,i+is) - v(0,k+ks,j+js,i+is) - v(1,k+ks,j+js,i+is-1) + v(0,k+ks,j+js-1,i+is));
       hv = 1./4. * (h(0,k+ks,j+js,i+is) + h(0,k+ks,j+js-1,i+is) + h(0,k+ks,j+js,i+is-1) + h(0,k+ks,j+js-1,i+is-1));
       eta = zeta + coriolis(0,k+ks,j+js,i+is);
       // This q is the primal (twisted) point value (0-form)
       q(0,k+ks,j+js,i+is) = eta / hv ; //* geom.get_J_dual_cell(k+ks, j+js, i+is)
   }
           });

   }

   void compute_diag(const VariableSet<ndims, nconst> &const_vars, VariableSet<ndims, nprog> &x, VariableSet<ndims, ndiag> &diagnostic_vars)
   {

   compute_diagnostic_quantities(diagnostic_vars.fields_arr[Q0VAR].data, diagnostic_vars.fields_arr[SL0VAR].data, x.fields_arr[VVAR].data, x.fields_arr[HVAR].data, x.fields_arr[SVAR].data, const_vars.fields_arr[CORIOLISVAR].data, *this->topology, *this->geom);

   }

};
// *******   Tendencies   ***********//

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE
template <uint ndims, uint nprog, uint nconst, uint naux> class Tendencies {
public:

  const Topology<ndims> *topology;
  ExchangeSet<ndims, naux> *aux_exchange;
  Geometry<ndims,1,1,1> *geom;

  bool is_initialized;

   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

  void initialize(const Topology<ndims> &topo, Geometry<ndims,1,1,1> &geom, ExchangeSet<ndims, naux> &aux_exchange)
  {
    this->topology = &topo;
    this->geom = &geom;
    this->aux_exchange = &aux_exchange;
    this->is_initialized = true;
  }

void YAKL_INLINE compute_primal_reconstruction(realArr reconvar, realArr densityvar, realArr fluxvar)
{
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 2)
    { cfv2_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 4)
    { cfv4_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 6)
    { cfv6_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 8)
    { cfv8_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 10)
    { cfv10_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}

    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 1)
    { weno1_recon<ndims, 1>(reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 3)
    { weno3_recon<ndims, 1>(reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 5)
    { weno5_recon<ndims, 1>(reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 7)
    { weno7_recon<ndims, 1>(reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 9)
    { weno9_recon<ndims, 1>(reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 11)
    { weno11_recon<ndims, 1>(reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
}

void YAKL_INLINE compute_dual_reconstruction(realArr reconvar, realArr densityvar, realArr fluxvar)
{
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 2)
    { cfv2_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 4)
    { cfv4_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 6)
    { cfv6_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 8)
    { cfv8_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 10)
    { cfv10_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}

    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 1)
    { weno1_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 3)
    { weno3_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 5)
    { weno5_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 7)
    { weno7_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 9)
    { weno9_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 11)
    { weno11_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
}


void YAKL_INLINE compute_auxiliary_quantities(realArr B, realArr F, realArr T, realArr q, realArr hrecon, realArr srecon, const realArr v, const realArr h, const realArr S, const realArr hs, const realArr coriolis, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom)
{

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

    real he0, he1, KE, zeta, eta, hv, U, V;

      yakl::parallel_for("ComputeAuxiliary", topology.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);

        he0 = 0.5 * (h(0, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + h(0, k+ks, j+js, i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1));
        U = geom.get_H_edge(0, k+ks, j+js, i+is) * v(0,k+ks,j+js,i+is);
        F(0,k+ks,j+js,i+is) = he0 * U;
        hrecon(0,k+ks,j+js,i+is) = hrecon(0,k+ks,j+js,i+is) / he0;
        //srecon(0,k+ks,j+js,i+is) = srecon(0,k+ks,j+js,i+is) / he0;

        if (ndims == 1) {
        KE = 1./2. * ( v(0,k+ks,j+js,i+is) * U + v(0,k+ks,j+js,i+is+1) * geom.get_H_edge(0, k+ks, j+js, i+is+1) * v(0,k+ks,j+js,i+is+1));
        }

        if (ndims == 2) {

        he1 = 0.5 * (h(0, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + h(0, k+ks, j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is));;
        V = geom.get_H_edge(1, k+ks, j+js, i+is) * v(1,k+ks,j+js,i+is);
        F(1,k+ks,j+js,i+is) = he1 * V;
        hrecon(1,k+ks,j+js,i+is) = hrecon(1,k+ks,j+js,i+is) / he1;
        //srecon(1,k+ks,j+js,i+is) = srecon(1,k+ks,j+js,i+is) / he1;

        KE = 1./4. * ( v(0,k+ks,j+js,i+is) * U + v(0,k+ks,j+js,i+is+1) * geom.get_H_edge(0, k+ks, j+js, i+is+1) * v(0,k+ks,j+js,i+is+1) +
                       v(1,k+ks,j+js,i+is) * V + v(1,k+ks,j+js+1,i+is) * geom.get_H_edge(1, k+ks, j+js+1, i+is) * v(1,k+ks,j+js+1,i+is));

        zeta = (v(1,k+ks,j+js,i+is) - v(0,k+ks,j+js,i+is) - v(1,k+ks,j+js,i+is-1) + v(0,k+ks,j+js-1,i+is));
        hv = 1./4. * (h(0,k+ks,j+js,i+is) + h(0,k+ks,j+js-1,i+is) + h(0,k+ks,j+js,i+is-1) + h(0,k+ks,j+js-1,i+is-1));
        eta = zeta + coriolis(0,k+ks,j+js,i+is);
        q(0,k+ks,j+js,i+is) = eta / hv * geom.get_J_dual_cell(k+ks, j+js, i+is);
        }

        B(0,k+ks,j+js,i+is) = (S(0,k+ks,j+js,i+is)/2.0 + KE)/geom.get_J_cell(k+ks, j+js, i+is);
        T(0,k+ks,j+js,i+is) = (h(0,k+ks,j+js,i+is)/2.0 + hs(0,k+ks,j+js,i+is))/geom.get_J_cell(k+ks, j+js, i+is);

    });




}



  void compute_rhs(const VariableSet<ndims, nconst> &const_vars, VariableSet<ndims, nprog> &x, VariableSet<ndims, naux> &auxiliary_vars, VariableSet<ndims, nprog> &xtend)
  {

      //Compute h and S reconstructions
   compute_primal_reconstruction(auxiliary_vars.fields_arr[HRECONVAR].data, x.fields_arr[HVAR].data, x.fields_arr[VVAR].data);

   int is = topology->is;
   int js = topology->js;
   int ks = topology->ks;

     yakl::parallel_for("ComputeAuxiliary", topology->n_cells, YAKL_LAMBDA (int iGlob) {
       int k, j, i;
       yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

   auxiliary_vars.fields_arr[SLVAR].data(0,k+ks,j+js,i+is) = x.fields_arr[SVAR].data(0,k+ks,j+js,i+is) / x.fields_arr[HVAR].data(0,k+ks,j+js,i+is) * geom->get_J_cell(k+ks, j+js, i+is);
    });

    compute_primal_reconstruction(auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[SLVAR].data, x.fields_arr[VVAR].data);

   //compute_primal_reconstruction(auxiliary_vars.fields_arr[SRECONVAR].data, x.fields_arr[SVAR].data, auxiliary_vars.fields_arr[VVAR].data);


//Compute B, F and q; also scale hrecon
compute_auxiliary_quantities(
auxiliary_vars.fields_arr[BVAR].data, auxiliary_vars.fields_arr[FVAR].data, auxiliary_vars.fields_arr[TVAR].data,
auxiliary_vars.fields_arr[QVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data,
x.fields_arr[VVAR].data, x.fields_arr[HVAR].data, x.fields_arr[SVAR].data,
const_vars.fields_arr[HSVAR].data, const_vars.fields_arr[CORIOLISVAR].data,
*this->topology, *this->geom);


//Compute FT and q reconstruction
if (ndims == 2) {
    // STILL BROKEN- DUAL grid flux is W H v!
    // ONLY WORKS FOR UNIFORM GRIDS...
//W2D_2(auxiliary_vars.fields_arr[FTVAR].data, x.fields_arr[VVAR].data, *this->topology);

// This is correct- dual grid flux is W F
W2D_2(auxiliary_vars.fields_arr[FTVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology);
compute_dual_reconstruction(auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[QVAR].data, auxiliary_vars.fields_arr[FTVAR].data);
}


   this->aux_exchange->exchange_variable_set(auxiliary_vars);

   //compute h rhs = D (hrecon* U) = D (hrecon/he F) with F = he U
   if (differential_order == 2)
   { divergence2<ndims, 1>(xtend.fields_arr[HVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 4)
   { divergence4<ndims, 1>(xtend.fields_arr[HVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 6)
   { divergence6<ndims, 1>(xtend.fields_arr[HVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 8)
   { divergence8<ndims, 1>(xtend.fields_arr[HVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }

   //compute S rhs = D (srecon F)
   if (differential_order == 2)
   { divergence2<ndims, 1>(xtend.fields_arr[SVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 4)
   { divergence4<ndims, 1>(xtend.fields_arr[SVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 6)
   { divergence6<ndims, 1>(xtend.fields_arr[SVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 8)
   { divergence8<ndims, 1>(xtend.fields_arr[SVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology); }


   //compute u rhs = hrecon/he G B + srecon/he G T + Q F
if (differential_order == 2)
{ gradient2<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[BVAR].data, *this->topology); }
if (differential_order == 4)
{ gradient4<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[BVAR].data, *this->topology); }
if (differential_order == 6)
{ gradient6<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[BVAR].data, *this->topology); }
if (differential_order == 8)
{ gradient8<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[HRECONVAR].data, auxiliary_vars.fields_arr[BVAR].data, *this->topology); }

if (differential_order == 2)
{ gradient2_add<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[TVAR].data, *this->topology); }
if (differential_order == 4)
{ gradient4_add<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[TVAR].data, *this->topology); }
if (differential_order == 6)
{ gradient6_add<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[TVAR].data, *this->topology); }
if (differential_order == 8)
{ gradient8_add<ndims, 1>(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[SRECONVAR].data, auxiliary_vars.fields_arr[TVAR].data, *this->topology); }

if (ndims == 2) {
Q2D_2_add(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology);
//Q2D_nonEC_2_add(xtend.fields_arr[VVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, auxiliary_vars.fields_arr[FVAR].data, *this->topology);
}
 }

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



template <uint ndims, uint nprog, uint nconst, uint nstats> class Stats
{
public:
  std::array<Stat,nstats> stats_arr;
  MPI_Request Req [nstats];
  MPI_Status  Status[nstats];
  int ierr;
  int statsize;
  int masterproc;
  Geometry<ndims,1,1,1> *geom;
  const Topology<ndims> *topology;

  void initialize(ModelParameters &params, Parallel &par, const Topology<ndims> &topo, Geometry<ndims,1,1,1> &geom)
  {
      this->topology = &topo;
      this->geom = &geom;
    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[MSTAT].initialize("mass", params, par);
    stats_arr[PESTAT].initialize("potential_energy", params, par);
    stats_arr[KESTAT].initialize("kinetic_energy", params, par);
    stats_arr[TESTAT].initialize("total_energy", params, par);
    stats_arr[BSTAT].initialize("bouyancy", params, par);

    if (ndims == 2) {
    stats_arr[PVSTAT].initialize("potential_vorticity", params, par);
    }
    masterproc = par.masterproc;
  }



  void compute( VariableSet<ndims, nprog> &progvars,  VariableSet<ndims, nconst> &constvars, int i)
  {

      this->stats_arr[MSTAT].local_dat = 0;
      this->stats_arr[PESTAT].local_dat = 0;
      this->stats_arr[KESTAT].local_dat = 0;
      this->stats_arr[TESTAT].local_dat = 0;
      this->stats_arr[BSTAT].local_dat = 0;
      if (ndims == 2) {
      this->stats_arr[PVSTAT].local_dat = 0;
  }

      int is = topology->is;
      int js = topology->js;
      int ks = topology->ks;

      real h, S, h1, h2, h3, u0, u1, v0, v1, um1, vm1, hs, coriolis;
      real zeta, eta, hv, q0, he0, he1, KE, PE;

        yakl::parallel_for("ComputeStats", topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, topology->n_cells_z, topology->n_cells_y, topology->n_cells_x, k, j, i);

          //compute stats locally
          h = progvars.fields_arr[HVAR].data(0,k+ks,j+js,i+is);
          S = progvars.fields_arr[SVAR].data(0,k+ks,j+js,i+is);
          h2 = progvars.fields_arr[HVAR].data(0,k+ks,j+js,i+is-1);
          u0 = progvars.fields_arr[VVAR].data(0,k+ks,j+js,i+is);
          u1 = progvars.fields_arr[VVAR].data(0,k+ks,j+js,i+is+1);
          hs = constvars.fields_arr[HSVAR].data(0,k+ks,j+js,i+is);

          if (ndims == 2)
          {
          h1 = progvars.fields_arr[HVAR].data(0,k+ks,j+js-1,i+is);
          h3 = progvars.fields_arr[HVAR].data(0,k+ks,j+js-1,i+is-1);
          v0 = progvars.fields_arr[VVAR].data(1,k+ks,j+js,i+is);
          v1 = progvars.fields_arr[VVAR].data(1,k+ks,j+js+1,i+is);
          vm1 = progvars.fields_arr[VVAR].data(1,k+ks,j+js,i+is-1); //note i instead of j
          um1 = progvars.fields_arr[VVAR].data(0,k+ks,j+js-1,i+is); //note j instead of i
          coriolis = constvars.fields_arr[CORIOLISVAR].data(0,k+ks,j+js,i+is);
            }

            he0 = 0.5 * (h/geom->get_J_cell(k+ks, j+js, i+is) + h2/geom->get_J_cell(k+ks, j+js, i+is-1));
          if (ndims == 1) {
          KE = 1./2. * he0 * ( u0   * geom->get_H_edge(0, k+ks, j+js, i+is)   * u0);
          }

          if (ndims == 2) {
              he1 = 0.5 * (h/geom->get_J_cell(k+ks, j+js, i+is) + h1/geom->get_J_cell(k+ks, j+js-1, i+is));

          KE = 1./2. * (he0 * ( u0 * geom->get_H_edge(0, k+ks, j+js, i+is) * u0 ) +
                      + he1 * ( v0 * geom->get_H_edge(1, k+ks, j+js, i+is) * v0 ));

        zeta = (v0 - u0 - vm1 + um1);
        hv = 1./4. * (h + h1 + h2 + h3);
        eta = zeta + coriolis;
        q0 = eta / hv;
          }

          PE = 0.5*S*h/geom->get_J_cell(k+ks, j+js, i+is) + S*hs/geom->get_J_cell(k+ks, j+js, i+is);

         this->stats_arr[MSTAT].local_dat += h;
         this->stats_arr[BSTAT].local_dat += S;
         this->stats_arr[PESTAT].local_dat += PE;
         this->stats_arr[KESTAT].local_dat += KE;
         this->stats_arr[TESTAT].local_dat += KE + PE;
         if (ndims == 2) {
         this->stats_arr[PVSTAT].local_dat += hv * q0;
     }

      });


    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &this->stats_arr[MSTAT].local_dat, &this->stats_arr[MSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[MSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[PESTAT].local_dat, &this->stats_arr[PESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PESTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[KESTAT].local_dat, &this->stats_arr[KESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[KESTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[TESTAT].local_dat, &this->stats_arr[TESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[TESTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[PVSTAT].local_dat, &this->stats_arr[PVSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[BSTAT].local_dat, &this->stats_arr[BSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[BSTAT]);

    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


  if (masterproc)
  {
  this->stats_arr[MSTAT].data(i) = this->stats_arr[MSTAT].global_dat;
  this->stats_arr[PESTAT].data(i) = this->stats_arr[PESTAT].global_dat;
  this->stats_arr[KESTAT].data(i) = this->stats_arr[KESTAT].global_dat;
  this->stats_arr[TESTAT].data(i) = this->stats_arr[TESTAT].global_dat;
  this->stats_arr[BSTAT].data(i) = this->stats_arr[BSTAT].global_dat;
  if (ndims == 2) {
  this->stats_arr[PVSTAT].data(i) = this->stats_arr[PVSTAT].global_dat;
    }
  }
  }
};



// *******   VariableSet Initialization   ***********//
// h, v
// hs
// B, F, q

template <uint ndims, uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology<ndims> &topo,
SArray<int, nprog, 4> &prog_ndofs_arr, SArray<int, nconst, 4> &const_ndofs_arr, SArray<int, naux, 4> &aux_ndofs_arr, SArray<int, ndiag, 4> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology<ndims> *, nprog> &prog_topo_arr, std::array<const Topology<ndims> *, nconst> &const_topo_arr, std::array<const Topology<ndims> *, naux> &aux_topo_arr, std::array<const Topology<ndims> *, ndiag> &diag_topo_arr)
{
  prog_topo_arr[HVAR] = &topo;
  prog_topo_arr[VVAR] = &topo;
  prog_topo_arr[SVAR] = &topo;
  const_topo_arr[HSVAR] = &topo;
  aux_topo_arr[BVAR] = &topo;
  aux_topo_arr[FVAR] = &topo;
  aux_topo_arr[FTVAR] = &topo;
  aux_topo_arr[TVAR] = &topo;
  aux_topo_arr[HRECONVAR] = &topo;
  aux_topo_arr[SRECONVAR] = &topo;
  diag_topo_arr[SL0VAR] = &topo;
  aux_topo_arr[SLVAR] = &topo;
  if (ndims == 2) {
  const_topo_arr[CORIOLISVAR] = &topo;
  aux_topo_arr[QVAR] = &topo;
  diag_topo_arr[Q0VAR] = &topo;
  aux_topo_arr[QRECONVAR] = &topo;
}

  prog_names_arr[HVAR] = "h";
  prog_names_arr[VVAR] = "v";
  prog_names_arr[SVAR] = "S";
  const_names_arr[HSVAR] = "hs";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[TVAR] = "T";
  aux_names_arr[HRECONVAR] = "hrecon";
  aux_names_arr[SRECONVAR] = "srecon";
  diag_names_arr[SL0VAR] = "sl";
  aux_names_arr[SLVAR] = "sl";

  if (ndims == 2) {
  const_names_arr[CORIOLISVAR] = "coriolis";
  aux_names_arr[QVAR] = "q";
  diag_names_arr[Q0VAR] = "q";
  aux_names_arr[QRECONVAR] = "qrecon";
}

//primal grid represents twisted quantities, dual grid straight quantities
    prog_ndofs_arr(HVAR,2) = 1; //h = twisted 2-form
    prog_ndofs_arr(SVAR,2) = 1; //S = twisted 2-form
    prog_ndofs_arr(VVAR,1) = 1; //v = straight 1-form
    const_ndofs_arr(HSVAR,2) = 1; //hs = twisted 2-form
    aux_ndofs_arr(BVAR,2) = 1; //B = straight 0-form
    aux_ndofs_arr(FVAR,1) = 1; //F = twisted 1-form
    aux_ndofs_arr(FTVAR,1) = 1; //FT = straight 1-form
    aux_ndofs_arr(TVAR,2) = 1; //T = straight 0-form
    aux_ndofs_arr(HRECONVAR,1) = 1; //hrecon lives on edges
    aux_ndofs_arr(SRECONVAR,1) = 1; //srecon lives on edges
    aux_ndofs_arr(SLVAR,2) = 1; //sl = straight 0-form
    diag_ndofs_arr(SL0VAR,2) = 1; //sl0 = straight 0-form

    if (ndims == 2) {
    const_ndofs_arr(CORIOLISVAR,0) = 1; //f = straight 2-form
    aux_ndofs_arr(QVAR,0) = 1; //q = straight 2-form
    diag_ndofs_arr(Q0VAR,0) = 1; //q0 = twisted 0-form
    aux_ndofs_arr(QRECONVAR,1) = 1; //qrecon lives on edges
}
}

  // *******   Initial Conditions   ***********//


// FIX THESE FOR 1D/2D
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


//wavespeed = sqrt(g * H0)
//dt = Constant(get_dt(wavespeed, cval, order, variant, Lx, Ly, nx, ny))

// FIX THESE FOR 1D/2D

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<1, nprog> &progvars, VariableSet<1, nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{

}

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<2, nprog> &progvars, VariableSet<2, nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{

    if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
    {
        geom.set_primal_2form_values(double_vortex_h, progvars.fields_arr[HVAR], 0);
        geom.set_primal_2form_values(double_vortex_S, progvars.fields_arr[SVAR], 0);
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
