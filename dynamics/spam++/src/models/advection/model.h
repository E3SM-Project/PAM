#ifndef _MODEL_H_
#define _MODEL_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "divergence.h"
#include "weno.h"
#include "cfv.h"
#include "geometry.h"
#include "params.h"
#include "string.h"



// Number of variables
uint constexpr nprognostic = 1;
uint constexpr nconstant = 1;
uint constexpr ndiagnostic = 1;
uint constexpr nstats = 3;

// Initial conditions related variables and functions
#define gaussian_1d(x)     (1. * exp(-100. * pow(x-0.5,2.)))
#define gaussian_2d(x,y)   (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)))
#define gaussian_3d(x,y,z) (1. * exp(-100. * pow(x-0.5,2.)) * exp(-100. * pow(y-0.5,2.)) * exp(-100. * pow(z-0.5,2.)))
real YAKL_INLINE gaussian(real x)                 { return gaussian_1d(x); }
real YAKL_INLINE gaussian(real x, real y)         { return gaussian_2d(x,y); }
real YAKL_INLINE gaussian(real x, real y, real z) { return gaussian_3d(x,y,z); }

#define vortex1_1d(x)     (1. *  exp(-100. * pow(x-0.75,2.)))
#define vortex1_2d(x,y)   (1. *  exp(-100. * pow(x-0.75,2.)) * exp(-100. * pow(y-0.75,2.)))
#define vortex1_3d(x,y,z) (1. *  exp(-100. * pow(x-0.75,2.)) * exp(-100. * pow(y-0.75,2.)) * exp(-100. * pow(z-0.75,2.)))
#define vortex2_1d(x)     (0.5 * exp(-50.  * pow(x-0.25,2.)))
#define vortex2_2d(x,y)   (0.5 * exp(-50.  * pow(x-0.25,2.)) * exp(-75.  * pow(y-0.25,2.)))
#define vortex2_3d(x,y,z) (0.5 * exp(-50.  * pow(x-0.25,2.)) * exp(-75.  * pow(y-0.25,2.)) * exp(-100. * pow(z-0.25,2.)))
real YAKL_INLINE vortices(real x)                 { return vortex1_1d(x)     + vortex2_1d(x); }
real YAKL_INLINE vortices(real x, real y)         { return vortex1_2d(x,y)   + vortex2_2d(x,y); }
real YAKL_INLINE vortices(real x, real y, real z) { return vortex1_3d(x,y,z) + vortex2_3d(x,y,z); }

real YAKL_INLINE square(real x)                 {return (x > 0.35 && x < 0.65                                                ) ? 1. : 0.;}
real YAKL_INLINE square(real x, real y)         {return (x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65                        ) ? 1. : 0.;}
real YAKL_INLINE square(real x, real y, real z) {return (x > 0.35 && x < 0.65 && y > 0.35 && y < 0.65 && z > 0.35 && z < 0.65) ? 1. : 0.;}





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



vec<2> YAKL_INLINE uniform_y_wind(real x, real y) {
  vec<2> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}
vec<3> YAKL_INLINE uniform_y_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}

vec<2> YAKL_INLINE uniform_xy_wind(real x, real y) {
  vec<2> vvec;
  vvec.u = sqrt(C_UNIFORM_WIND/2.);
  vvec.v = sqrt(C_UNIFORM_WIND/2.);
  return vvec;
}


vec<3> YAKL_INLINE uniform_z_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.w = C_UNIFORM_WIND;
  return vvec;
}

// FIX THIS
vec<2> YAKL_INLINE deformational_wind(real x, real y) {
  vec<2> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}
vec<3> YAKL_INLINE deformational_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.v = C_UNIFORM_WIND;
  return vvec;
}

// *******   Model Specific Parameters   ***********//


template<uint ndims> void set_model_specific_params(std::string inFile, ModelParameters &params)
{


  std::string strDataInit = "";
  std::string strWindInit = "";

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
      else if ( !strcmp( "windInit"   , key.c_str() ) ) { ssVal >> strWindInit       ;}
      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

  if (!strcmp("",strDataInit.c_str())) { std::cout << "Error: key " << "dataInit" << " not set.\n"; exit(-1); }
  if (!strcmp("",strWindInit.c_str())) { std::cout << "Error: key " << "windInit" << " not set.\n"; exit(-1); }

    size_t splitloc = strDataInit.find("//",0);
    std::string sub_str;
    if (splitloc != std::string::npos){
      sub_str = strDataInit.substr(0,splitloc);
    } else {
      sub_str = strDataInit;
    }
    if      ( !strcmp(sub_str.c_str(),"gaussian" ) ) { params.data_init_cond = DATA_INIT::GAUSSIAN  ; }
    else if ( !strcmp(sub_str.c_str(),"vortices" ) ) { params.data_init_cond = DATA_INIT::VORTICES  ; }
    else if ( !strcmp(sub_str.c_str(),"square"   ) ) { params.data_init_cond = DATA_INIT::SQUARE    ; }
    else  {
      std::cout << "Error: unrecognized dataInit " << strDataInit << "\n";
      exit(-1);
    }

    splitloc = strWindInit.find("//",0);
    if (splitloc != std::string::npos){
      sub_str = strWindInit.substr(0,splitloc);
    } else {
      sub_str = strWindInit;
    }
    if      ( !strcmp(sub_str.c_str(),"uniform_x"     ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_X     ; }
    else if ( !strcmp(sub_str.c_str(),"uniform_y"     ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_Y     ; }
    else if ( !strcmp(sub_str.c_str(),"uniform_z"     ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_Z     ; }
    else if ( !strcmp(sub_str.c_str(),"uniform_xy"    ) ) { params.wind_init_cond = WIND_INIT::UNIFORM_XY    ; }
    else if ( !strcmp(sub_str.c_str(),"deformational" ) ) { params.wind_init_cond = WIND_INIT::DEFORMATIONAL ; }
    else  {
      std::cout << "Error: unrecognized windInit " << strWindInit << "\n";
      exit(-1);
    }

  params.etime = 0.0;

  if (params.data_init_cond == DATA_INIT::GAUSSIAN || params.data_init_cond == DATA_INIT::VORTICES || params.data_init_cond == DATA_INIT::SQUARE)
  {
  params.xlen = 1.0;
  params.xc = 0.5;
  if (ndims>=2)
  {
    params.ylen = 1.0;
    params.yc = 0.5;
  }
  if (ndims ==3)
  {
  params.zlen = 1.0;
  params.zc = 0.5;
  }
  }

}


// *******   Tendencies   ***********//

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
  MPI_Request Req [3];
  MPI_Status  Status[3];
  int ierr;
  int statsize;
  int masterproc;

  void initialize(ModelParameters &params, Parallel &par)
  {
    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[0].initialize("qmass", params, par);
    stats_arr[1].initialize("qmin", params, par);
    stats_arr[2].initialize("qmax", params, par);
    masterproc = par.masterproc;
  }




  void compute( VariableSet<ndims, nprog> &progvars,  VariableSet<ndims, nprog> &constvars, int i)
  {

    //compute locally
    this->stats_arr[0].local_dat = progvars.fields_arr[0].sum();
    this->stats_arr[1].local_dat = progvars.fields_arr[0].min();
    this->stats_arr[2].local_dat = progvars.fields_arr[0].max();

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &this->stats_arr[0].local_dat, &this->stats_arr[0].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[0]);
    this->ierr = MPI_Ireduce( &this->stats_arr[1].local_dat, &this->stats_arr[1].global_dat, 1, REAL_MPI, MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[1]);
    this->ierr = MPI_Ireduce( &this->stats_arr[2].local_dat, &this->stats_arr[2].global_dat, 1, REAL_MPI, MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[2]);
    this->ierr = MPI_Waitall(3, this->Req, this->Status);

  if (masterproc)
  {
  this->stats_arr[0].data(i) = this->stats_arr[0].global_dat;
  this->stats_arr[1].data(i) = this->stats_arr[1].global_dat;
  this->stats_arr[2].data(i) = this->stats_arr[2].global_dat;
  }
  }
};



// *******   VariableSet Initialization   ***********//

template <uint ndims, uint nprog, uint nconst, uint ndiag> void initialize_variables(const Topology<ndims> &topo,
SArray<int, nprognostic, 4> &prog_ndofs_arr, SArray<int, nconstant, 4> &const_ndofs_arr, SArray<int, ndiagnostic, 4> &diag_ndofs_arr,
std::array<std::string, nprognostic> &prog_names_arr, std::array<std::string, nconstant> &const_names_arr, std::array<std::string, ndiagnostic> &diag_names_arr,
std::array<const Topology<ndims> *, nprognostic> &prog_topo_arr, std::array<const Topology<ndims> *, nconstant> &const_topo_arr, std::array<const Topology<ndims> *, ndiagnostic> &diag_topo_arr)
{
  prog_topo_arr[0] = &topo;
  const_topo_arr[0] = &topo;
  diag_topo_arr[0] = &topo;
  prog_names_arr[0] = "q";
  const_names_arr[0] = "u";
  diag_names_arr[0] = "qrecon";

  if (ndims == 1) {
    prog_ndofs_arr(0,1) = nqdofs;
    const_ndofs_arr(0,0) = 1;
    diag_ndofs_arr(0,0) = nqdofs;
  }

  if (ndims == 2) {
    prog_ndofs_arr(0,2) = nqdofs;
    const_ndofs_arr(0,1) = 1;
    diag_ndofs_arr(0,1) = nqdofs;
  }

  if (ndims == 3) {
    prog_ndofs_arr(0,3) = nqdofs;
    const_ndofs_arr(0,2) = 1;
    diag_ndofs_arr(0,2) = nqdofs;
  }

}

  // *******   Initial Conditions   ***********//

template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<1, nprog> &progvars, VariableSet<1, nconst> &constvars, Geometry<1, nquadx, nquady, nquadz> &geom)
{

    for (int i=0; i<nqdofs; i++)
    {
    if (params.data_init_cond == DATA_INIT::GAUSSIAN) {geom.set_1form_values(gaussian, progvars.fields_arr[0], i);}
    if (params.data_init_cond == DATA_INIT::SQUARE)   {geom.set_1form_values(square,   progvars.fields_arr[0], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X ) {geom.set_0form_values(uniform_x_wind,     constvars.fields_arr[0], 0);}
}

template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<2, nprog> &progvars, VariableSet<2, nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{
    for (int i=0; i<nqdofs; i++)
    {
    if (params.data_init_cond == DATA_INIT::GAUSSIAN) {geom.set_2form_values(gaussian, progvars.fields_arr[0], i);}
    if (params.data_init_cond == DATA_INIT::VORTICES) {geom.set_2form_values(vortices, progvars.fields_arr[0], i);}
    if (params.data_init_cond == DATA_INIT::SQUARE)   {geom.set_2form_values(square,   progvars.fields_arr[0], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_1form_values(uniform_x_wind,     constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Y    ) {geom.set_1form_values(uniform_y_wind,     constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_XY   ) {geom.set_1form_values(uniform_xy_wind,    constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (params.wind_init_cond == WIND_INIT::DEFORMATIONAL) {geom.set_1form_values(deformational_wind, constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}

}


template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<3, nprog> &progvars, VariableSet<3, nconst> &constvars, Geometry<3, nquadx, nquady, nquadz> &geom)
{
    for (int i=0; i<nqdofs; i++)
    {
    if (params.data_init_cond == DATA_INIT::GAUSSIAN) {geom.set_3form_values(gaussian, progvars.fields_arr[0], i);}
    if (params.data_init_cond == DATA_INIT::SQUARE)   {geom.set_3form_values(square,   progvars.fields_arr[0], i);}
    }
    if (params.wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_2form_values(uniform_x_wind,     constvars.fields_arr[0], 0);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Y    ) {geom.set_2form_values(uniform_y_wind,     constvars.fields_arr[0], 0);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Z    ) {geom.set_2form_values(uniform_z_wind,     constvars.fields_arr[0], 0);}
    if (params.wind_init_cond == WIND_INIT::DEFORMATIONAL) {geom.set_2form_values(deformational_wind, constvars.fields_arr[0], 0);}
}

#endif
