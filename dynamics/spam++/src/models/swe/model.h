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
#include "W2D.h"
#include "Q2D.h"
#include "geometry.h"
#include "params.h"
#include "string.h"



// Number of variables
uint constexpr nprognostic = 2; // h, v
uint constexpr nconstant = 2;   // hs, coriolis
uint constexpr ndiagnostic = 5; // B, F, q, hrecon, qrecon
uint constexpr nstats = 4;      // M, TE, PE, PV

#define HVAR 0
#define VVAR 1

#define HSVAR 0
#define CORIOLISVAR 1

#define BVAR 0
#define FVAR 1
#define QVAR 2
#define HRECONVAR 3
#define QRECONVAR 4

#define MSTAT 0
#define TESTAT 1
#define PVSTAT 2
#define PESTAT 3

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


// *******   Tendencies   ***********//

// THIS SHOULD BE GENERALIZABLE...
// ONLY COMPUTE FUNCTION NEEDS TO CHANGE
template <uint ndims, uint nprog, uint nconst, uint ndiag> class Tendencies {
public:

  const Topology<ndims> *topology;
  ExchangeSet<ndims, ndiag> *diag_exchange;
  Geometry<ndims,1,1,1> *geom;

  bool is_initialized;

   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

  void initialize(const Topology<ndims> &topo, Geometry<ndims,1,1,1> &geom, ExchangeSet<ndims, ndiag> &diag_exchange)
  {
    this->topology = &topo;
    this->geom = &geom;
    this->diag_exchange = &diag_exchange;
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
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 2)
    { cfv2_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 4)
    { cfv4_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 6)
    { cfv6_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 8)
    { cfv8_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}
    if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 10)
    { cfv10_dual_recon<ndims, 1>(reconvar, densityvar, *this->topology, *this->geom);}

    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 1)
    { weno1_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 3)
    { weno3_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 5)
    { weno5_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 7)
    { weno7_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 9)
    { weno9_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
    if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 11)
    { weno11_dual_recon<ndims, 1>(false, reconvar, densityvar, fluxvar, *this->topology, *this->geom); }
}

void YAKL_INLINE compute_B()
{

}

void YAKL_INLINE compute_F()
{

}

void YAKL_INLINE compute_q()
{

}

// THIS NEEDS MAJOR FIXING
// compute B, F, q
// compute rhs for h and v
  void compute_rhs(const VariableSet<ndims, nconst> &const_vars, VariableSet<ndims, nprog> &x, VariableSet<ndims, ndiag> &diagnostic_vars, VariableSet<ndims, nprog> &xtend)
  {

   //Compute B: requires v, h, hs

//Compute F: requires v, h

   //Compute q: requires v, f, h

   //Compute h reconstructions
compute_primal_reconstruction(diagnostic_vars.fields_arr[HRECONVAR].data, x.fields_arr[HVAR].data, diagnostic_vars.fields_arr[FVAR].data);

   //Compute q reconstructions
compute_dual_reconstruction(diagnostic_vars.fields_arr[QRECONVAR].data, diagnostic_vars.fields_arr[QVAR].data, diagnostic_vars.fields_arr[FVAR].data);

   this->diag_exchange->exchange_variable_set(diagnostic_vars);

   //compute h rhs = D (hrecon F)
   if (differential_order == 2)
   { divergence2<ndims, 1>(xtend.fields_arr[HVAR].data, diagnostic_vars.fields_arr[HRECONVAR].data, diagnostic_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 4)
   { divergence4<ndims, 1>(xtend.fields_arr[HVAR].data, diagnostic_vars.fields_arr[HRECONVAR].data, diagnostic_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 6)
   { divergence6<ndims, 1>(xtend.fields_arr[HVAR].data, diagnostic_vars.fields_arr[HRECONVAR].data, diagnostic_vars.fields_arr[FVAR].data, *this->topology); }
   if (differential_order == 8)
   { divergence8<ndims, 1>(xtend.fields_arr[HVAR].data, diagnostic_vars.fields_arr[HRECONVAR].data, diagnostic_vars.fields_arr[FVAR].data, *this->topology); }

// FIX THIS
   //compute u rhs = hrecon G B + Q F
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

  void initialize(ModelParameters &params, Parallel &par)
  {
    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[MSTAT].initialize("mass", params, par);
    stats_arr[TESTAT].initialize("energy", params, par);
    stats_arr[PVSTAT].initialize("potential vorticity", params, par);
    stats_arr[PESTAT].initialize("potential enstrophy", params, par);
    masterproc = par.masterproc;
  }



// FIX THESE COMPUTATIONS
  void compute( VariableSet<ndims, nprog> &progvars,  VariableSet<ndims, nprog> &constvars, int i)
  {
     // DO A LOOP HERE!

     //compute q locally

    //compute stats locally
    this->stats_arr[MSTAT].local_dat = progvars.fields_arr[HVAR].sum();
    this->stats_arr[TESTAT].local_dat = progvars.fields_arr[0].sum();
    this->stats_arr[PVSTAT].local_dat = progvars.fields_arr[0].sum();
    this->stats_arr[PESTAT].local_dat = progvars.fields_arr[0].sum();

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &this->stats_arr[MSTAT].local_dat, &this->stats_arr[MSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[0]);
    this->ierr = MPI_Ireduce( &this->stats_arr[TESTAT].local_dat, &this->stats_arr[TESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[1]);
    this->ierr = MPI_Ireduce( &this->stats_arr[PVSTAT].local_dat, &this->stats_arr[PVSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[2]);
    this->ierr = MPI_Ireduce( &this->stats_arr[PESTAT].local_dat, &this->stats_arr[PESTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[3]);
    this->ierr = MPI_Waitall(4, this->Req, this->Status);

  if (masterproc)
  {
  this->stats_arr[MSTAT].data(i) = this->stats_arr[0].global_dat;
  this->stats_arr[TESTAT].data(i) = this->stats_arr[1].global_dat;
  this->stats_arr[PVSTAT].data(i) = this->stats_arr[2].global_dat;
  this->stats_arr[PESTAT].data(i) = this->stats_arr[3].global_dat;
  }
  }
};



// *******   VariableSet Initialization   ***********//
// h, v
// hs
// B, F, q

template <uint ndims, uint nprog, uint nconst, uint ndiag> void initialize_variables(const Topology<ndims> &topo,
SArray<int, nprognostic, 4> &prog_ndofs_arr, SArray<int, nconstant, 4> &const_ndofs_arr, SArray<int, ndiagnostic, 4> &diag_ndofs_arr,
std::array<std::string, nprognostic> &prog_names_arr, std::array<std::string, nconstant> &const_names_arr, std::array<std::string, ndiagnostic> &diag_names_arr,
std::array<const Topology<ndims> *, nprognostic> &prog_topo_arr, std::array<const Topology<ndims> *, nconstant> &const_topo_arr, std::array<const Topology<ndims> *, ndiagnostic> &diag_topo_arr)
{
  prog_topo_arr[HVAR] = &topo;
  prog_topo_arr[VVAR] = &topo;
  const_topo_arr[HSVAR] = &topo;
  const_topo_arr[CORIOLISVAR] = &topo;
  diag_topo_arr[BVAR] = &topo;
  diag_topo_arr[FVAR] = &topo;
  diag_topo_arr[QVAR] = &topo;
  diag_topo_arr[HRECONVAR] = &topo;
  diag_topo_arr[QRECONVAR] = &topo;


  prog_names_arr[HVAR] = "h";
  prog_names_arr[VVAR] = "v";
  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";
  diag_names_arr[BVAR] = "B";
  diag_names_arr[FVAR] = "F";
  diag_names_arr[QVAR] = "q";
  diag_names_arr[HRECONVAR] = "hrecon";
  diag_names_arr[QRECONVAR] = "qrecon";

//primal grid represents twisted quantities, dual grid straight quantities
    prog_ndofs_arr(HVAR,2) = 1; //h = twisted 2-form
    prog_ndofs_arr(VVAR,1) = 1; //v = straight 1-form
    const_ndofs_arr(HSVAR,2) = 1; //hs = twisted 2-form
    const_ndofs_arr(CORIOLISVAR,0) = 1; //f = twisted 0-form
    diag_ndofs_arr(BVAR,2) = 1; //B = straight 0-form
    diag_ndofs_arr(FVAR,1) = 1; //F = twisted 1-form
    diag_ndofs_arr(QVAR,0) = 1; //q = twisted 0-form
    diag_ndofs_arr(HRECONVAR,1) = 1; //hrecon lives on edges
    diag_ndofs_arr(QRECONVAR,1) = 1; //qrecon lives on edges
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


//wavespeed = sqrt(g * H0)
//dt = Constant(get_dt(wavespeed, cval, order, variant, Lx, Ly, nx, ny))

template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<2, nprog> &progvars, VariableSet<2, nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{

    if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
    {
        geom.set_primal_2form_values(double_vortex_h, progvars.fields_arr[HVAR], 0);
        geom.set_dual_1form_values(double_vortex_v, progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
        geom.set_primal_0form_values(double_vortex_coriolis, constvars.fields_arr[CORIOLISVAR], 0);
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
