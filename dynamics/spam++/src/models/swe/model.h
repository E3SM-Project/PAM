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
uint constexpr nconstant = 1;   // hs
uint constexpr ndiagnostic = 3; // B, F, q
uint constexpr nstats = 4;      // M, TE, PE, PV

// ADD SUPPORT FOR 1D rswe version!
// WOULD INTRODUCE AN ADDITIONAL PROGNOSTIC VARIABLE THOUGH...UGGH
// Maybe we need actual prognostic varible set classes?
// Basically a bunch of base classes that have empty functions, and then overloaded versions for subclasses as needed?
// We know everything at compile time anyways!

// Initial conditions related variables and functions
// FIX THIS STUFF UP

#define gaussian_1d(x)     (1. * exp(-100. * pow(x-0.5,2.)))
#define gaussian_2d(x,y)   (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)))
#define gaussian_3d(x,y,z) (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)) * exp(-100. * pow(z-0.5,2.)))
real YAKL_INLINE gaussian(real x)                 { return gaussian_1d(x); }
real YAKL_INLINE gaussian(real x, real y)         { return gaussian_2d(x,y); }
real YAKL_INLINE gaussian(real x, real y, real z) { return gaussian_3d(x,y,z); }

#define C_UNIFORM_WIND 1.

real YAKL_INLINE uniform_x_wind(real x) {
  return C_UNIFORM_WIND;
}
vec<2> YAKL_INLINE uniform_x_wind(real x, real y) {
  vec<2> vvec;
  vvec.u = C_UNIFORM_WIND;
  return vvec;
}
vec<3> YAKL_INLINE uniform_x_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.u = C_UNIFORM_WIND;
  return vvec;
}

}

// *******   Model Specific Parameters   ***********//


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

// FIX THESE
  if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX || params.data_init_cond == DATA_INIT::GOLO || params.data_init_cond == DATA_INIT::TC2 || params.data_init_cond == DATA_INIT::TC5)
  {
  params.xlen = 1.0;
  params.xc = 0.5;
  params.ylen = 1.0;
  params.yc = 0.5;
  }

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

// THIS NEEDS MAJOR FIXING
// compute B, F, q
// compute rhs for h and v
  void compute_rhs(const VariableSet<ndims, nconst> &const_vars, VariableSet<ndims, nprog> &x, VariableSet<ndims, ndiag> &diagnostic_vars, VariableSet<ndims, nprog> &xtend)
  {

   //compute reconstructions


   if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 2)
   { cfv2_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, *this->topology, *this->geom);}
   if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 4)
   { cfv4_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, *this->topology, *this->geom);}
   if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 6)
   { cfv6_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, *this->topology, *this->geom);}
   if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 8)
   { cfv8_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, *this->topology, *this->geom);}
   if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 10)
   { cfv10_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, *this->topology, *this->geom);}

   if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 1)
   { weno1_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology, *this->geom); }
   if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 3)
   { weno3_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology, *this->geom); }
   if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 5)
   { weno5_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology, *this->geom); }
   if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 7)
   { weno7_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology, *this->geom); }
   if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 9)
   { weno9_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology, *this->geom); }
   if (reconstruction_type == RECONSTRUCTION_TYPE::WENO && reconstruction_order == 11)
   { weno11_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology, *this->geom); }

   this->diag_exchange->exchange_variable_set(diagnostic_vars);

   //compute D (qrecon U)
   if (differential_order == 2)
   { divergence2<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology); }
   if (differential_order == 4)
   { divergence4<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology); }
   if (differential_order == 6)
   { divergence6<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology); }
   if (differential_order == 8)
   { divergence8<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology); }
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
    stats_arr[0].initialize("mass", params, par);
    stats_arr[1].initialize("energy", params, par);
    stats_arr[2].initialize("potential vorticity", params, par);
    stats_arr[2].initialize("potential enstrophy", params, par);
    masterproc = par.masterproc;
  }



// FIX THESE COMPUTATIONS
  void compute( VariableSet<ndims, nprog> &progvars,  VariableSet<ndims, nprog> &constvars, int i)
  {

    //compute locally
    this->stats_arr[0].local_dat = progvars.fields_arr[0].sum();
    this->stats_arr[1].local_dat = progvars.fields_arr[0].sum();
    this->stats_arr[2].local_dat = progvars.fields_arr[0].sum();
    this->stats_arr[3].local_dat = progvars.fields_arr[0].sum();

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &this->stats_arr[0].local_dat, &this->stats_arr[0].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[0]);
    this->ierr = MPI_Ireduce( &this->stats_arr[1].local_dat, &this->stats_arr[1].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[1]);
    this->ierr = MPI_Ireduce( &this->stats_arr[2].local_dat, &this->stats_arr[2].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[2]);
    this->ierr = MPI_Ireduce( &this->stats_arr[2].local_dat, &this->stats_arr[2].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[2]);
    this->ierr = MPI_Waitall(4, this->Req, this->Status);

  if (masterproc)
  {
  this->stats_arr[0].data(i) = this->stats_arr[0].global_dat;
  this->stats_arr[1].data(i) = this->stats_arr[1].global_dat;
  this->stats_arr[2].data(i) = this->stats_arr[2].global_dat;
  this->stats_arr[3].data(i) = this->stats_arr[3].global_dat;
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
  prog_topo_arr[0] = &topo;
  prog_topo_arr[1] = &topo;
  const_topo_arr[0] = &topo;
  diag_topo_arr[0] = &topo;
  diag_topo_arr[1] = &topo;
  diag_topo_arr[2] = &topo;


  prog_names_arr[0] = "h";
  prog_names_arr[1] = "v";
  const_names_arr[0] = "hs";
  diag_names_arr[0] = "B";
  diag_names_arr[0] = "F";
  diag_names_arr[0] = "q";

    prog_ndofs_arr(0,2) = 1;
    prog_ndofs_arr(1,1) = 1;
    const_ndofs_arr(0,2) = 1;
    diag_ndofs_arr(0,2) = 1;
    diag_ndofs_arr(1,1) = 1;
    diag_ndofs_arr(2,0) = 1;
}

  // *******   Initial Conditions   ***********//

// FIX THIS
// Set h, v, hs, maybe coriolis?
template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<2, nprog> &progvars, VariableSet<2, nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{

    if (params.data_init_cond == DATA_INIT::DOUBLEVORTEX)
    {
        geom.set_primal_2form_values(gaussian, progvars.fields_arr[0], i);
    }

    if (params.data_init_cond == DATA_INIT::GOLO)
    {
        geom.set_primal_2form_values(gaussian, progvars.fields_arr[0], i);
    }

    if (params.data_init_cond == DATA_INIT::TC2)
    {
        geom.set_primal_2form_values(gaussian, progvars.fields_arr[0], i);
    }

    if (params.data_init_cond == DATA_INIT::TC5)
    {
        geom.set_primal_2form_values(gaussian, progvars.fields_arr[0], i);
    }

}



#endif
