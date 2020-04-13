#ifndef _MODEL_H_
#define _MODEL_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "Q2D.h"
#include "W2D.h"
#include "curl.h"
#include "weno_dual.h"
#include "cfv_dual.h"
#include "geometry.h"
#include "params.h"
#include "string.h"



// Number of variables
uint constexpr nprognostic = 1;
uint constexpr nconstant = 2;
uint constexpr nauxiliary = 2;
uint constexpr ndiagnostic = 0;
uint constexpr nstats = 3;

#define QVAR 0

#define UVAR 0
#define UTVAR 1

#define QRECONVAR 0
#define QFLUXVAR 1

#define MSTAT 0
#define MINSTAT 1
#define MAXSTAT 2

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

vec<3> YAKL_INLINE uniform_xy_wind(real x, real y, real z) {
  vec<3> vvec;
  vvec.u = sqrt(C_UNIFORM_WIND/2.);
  vvec.v = sqrt(C_UNIFORM_WIND/2.);
  vvec.w = 0;
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


      void compute_diag(const VariableSet<ndims, nconst> &const_vars, VariableSet<ndims, nprog> &x, VariableSet<ndims, ndiag> &diagnostic_vars)
      {
      }

   };

// *******   Tendencies   ***********//

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

  void compute_rhs(const VariableSet<ndims, nconst> &const_vars, VariableSet<ndims, nprog> &x, VariableSet<ndims, naux> &auxiliary_vars, VariableSet<ndims, nprog> &xtend)
  {

      //compute W U = dual grid flux
      W2D_2(const_vars.fields_arr[UTVAR].data, const_vars.fields_arr[UVAR].data, *this->topology);


   //compute reconstructions

   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 2)
   { cfv2_dual_recon<ndims, nqdofs>(auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, *this->topology, *this->geom);}
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 4)
   { cfv4_dual_recon<ndims, nqdofs>(auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, *this->topology, *this->geom);}
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 6)
   { cfv6_dual_recon<ndims, nqdofs>(auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, *this->topology, *this->geom);}
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 8)
   { cfv8_dual_recon<ndims, nqdofs>(auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, *this->topology, *this->geom);}
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::CFV && dual_reconstruction_order == 10)
   { cfv10_dual_recon<ndims, nqdofs>(auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, *this->topology, *this->geom);}

   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 1)
   { weno1_dual_recon<ndims, nqdofs>(false, auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, const_vars.fields_arr[UTVAR].data, *this->topology, *this->geom); }
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 3)
   { weno3_dual_recon<ndims, nqdofs>(false, auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, const_vars.fields_arr[UTVAR].data, *this->topology, *this->geom); }
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 5)
   { weno5_dual_recon<ndims, nqdofs>(false, auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, const_vars.fields_arr[UTVAR].data, *this->topology, *this->geom); }
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 7)
   { weno7_dual_recon<ndims, nqdofs>(false, auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, const_vars.fields_arr[UTVAR].data, *this->topology, *this->geom); }
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 9)
   { weno9_dual_recon<ndims, nqdofs>(false, auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, const_vars.fields_arr[UTVAR].data, *this->topology, *this->geom); }
   if (dual_reconstruction_type == RECONSTRUCTION_TYPE::WENO && dual_reconstruction_order == 11)
   { weno11_dual_recon<ndims, nqdofs>(false, auxiliary_vars.fields_arr[QRECONVAR].data, x.fields_arr[QVAR].data, const_vars.fields_arr[UTVAR].data, *this->topology, *this->geom); }

   //compute Q(qrecon, U)
   Q2D_2( auxiliary_vars.fields_arr[QFLUXVAR].data, auxiliary_vars.fields_arr[QRECONVAR].data, const_vars.fields_arr[UVAR].data, *this->topology);


   this->aux_exchange->exchange_variable_set(auxiliary_vars);

   //compute C Q(qrecon, U)

    curl2D_2( xtend.fields_arr[QVAR].data, auxiliary_vars.fields_arr[QFLUXVAR].data,  *this->topology);

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
  MPI_Request Req [nstats];
  MPI_Status  Status[nstats];
  int ierr;
  int statsize;
  int masterproc;

  void initialize(ModelParameters &params, Parallel &par, const Topology<ndims> &topo, Geometry<ndims,1,1,1> &geom)
  {
    statsize = params.Nsteps/params.Nstat + 1;
    stats_arr[MSTAT].initialize("qmass", params, par);
    stats_arr[MINSTAT].initialize("qmin", params, par);
    stats_arr[MAXSTAT].initialize("qmax", params, par);
    masterproc = par.masterproc;
  }




  void compute( VariableSet<ndims, nprog> &progvars,  VariableSet<ndims, nconst> &constvars, int i)
  {

    //compute locally
    this->stats_arr[MSTAT].local_dat = progvars.fields_arr[QVAR].sum();
    this->stats_arr[MINSTAT].local_dat = progvars.fields_arr[QVAR].min();
    this->stats_arr[MAXSTAT].local_dat = progvars.fields_arr[QVAR].max();

    //MPI sum/min/max
    this->ierr = MPI_Ireduce( &this->stats_arr[MSTAT].local_dat, &this->stats_arr[MSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[MSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[MINSTAT].local_dat, &this->stats_arr[MINSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[MINSTAT]);
    this->ierr = MPI_Ireduce( &this->stats_arr[MAXSTAT].local_dat, &this->stats_arr[MAXSTAT].global_dat, 1, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[MAXSTAT]);
    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);

  if (masterproc)
  {
  this->stats_arr[MSTAT].data(i) = this->stats_arr[MSTAT].global_dat;
  this->stats_arr[MINSTAT].data(i) = this->stats_arr[MINSTAT].global_dat;
  this->stats_arr[MAXSTAT].data(i) = this->stats_arr[MAXSTAT].global_dat;
  }
  }
};



// *******   VariableSet Initialization   ***********//

template <uint ndims, uint nprog, uint nconst, uint naux, uint ndiag> void initialize_variables(const Topology<ndims> &topo,
SArray<int, nprog, 4> &prog_ndofs_arr, SArray<int, nconst, 4> &const_ndofs_arr, SArray<int, naux, 4> &aux_ndofs_arr, SArray<int, ndiag, 4> &diag_ndofs_arr,
std::array<std::string, nprog> &prog_names_arr, std::array<std::string, nconst> &const_names_arr, std::array<std::string, naux> &aux_names_arr, std::array<std::string, ndiag> &diag_names_arr,
std::array<const Topology<ndims> *, nprog> &prog_topo_arr, std::array<const Topology<ndims> *, nconst> &const_topo_arr, std::array<const Topology<ndims> *, naux> &aux_topo_arr, std::array<const Topology<ndims> *, ndiag> &diag_topo_arr)
{
  prog_topo_arr[QVAR] = &topo;
  const_topo_arr[UVAR] = &topo;
  const_topo_arr[UTVAR] = &topo;
  aux_topo_arr[QRECONVAR] = &topo;
  aux_topo_arr[QFLUXVAR] = &topo;
  prog_names_arr[QVAR] = "q";
  const_names_arr[UVAR] = "u";
  const_names_arr[UTVAR] = "ut";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QFLUXVAR] = "qflux";


  if (ndims == 2) {
    prog_ndofs_arr(QVAR,0) = nqdofs;
    const_ndofs_arr(UVAR,1) = 1;
    const_ndofs_arr(UTVAR,1) = 1;
    aux_ndofs_arr(QRECONVAR,1) = nqdofs;
    aux_ndofs_arr(QFLUXVAR,1) = nqdofs;
  }



}

  // *******   Initial Conditions   ***********//

template <int nprog, int nconst, int nquadx, int nquady, int nquadz> void set_initial_conditions (ModelParameters &params, VariableSet<2, nprog> &progvars, VariableSet<2, nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{
    for (int i=0; i<nqdofs; i++)
    {
    if (params.data_init_cond == DATA_INIT::GAUSSIAN) {geom.set_dual_2form_values(gaussian, progvars.fields_arr[QVAR], i);}
    if (params.data_init_cond == DATA_INIT::VORTICES) {geom.set_dual_2form_values(vortices, progvars.fields_arr[QVAR], i);}
    if (params.data_init_cond == DATA_INIT::SQUARE)   {geom.set_dual_2form_values(square,   progvars.fields_arr[QVAR], i);}
    }

    if (params.wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_primal_1form_values(uniform_x_wind,     constvars.fields_arr[UVAR], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_Y    ) {geom.set_primal_1form_values(uniform_y_wind,     constvars.fields_arr[UVAR], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (params.wind_init_cond == WIND_INIT::UNIFORM_XY   ) {geom.set_primal_1form_values(uniform_xy_wind,    constvars.fields_arr[UVAR], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (params.wind_init_cond == WIND_INIT::DEFORMATIONAL) {geom.set_primal_1form_values(deformational_wind, constvars.fields_arr[UVAR], 0, LINE_INTEGRAL_TYPE::NORMAL);}


}

#endif
