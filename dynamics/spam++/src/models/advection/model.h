#ifndef _MODEL_H_
#define _MODEL_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "divergence.h"
#include "finitevolume.h"
#include "geometry.h"
#include "params.h"

// Number of variables
uint constexpr nprognostic = 1;
uint constexpr nconstant = 1;
uint constexpr ndiagnostic = 1;

// Initial conditions related variables and functions
#define A_GAUSSIAN 1.
#define AX_GAUSSIAN 0.1
#define AY_GAUSSIAN 0.1
#define AZ_GAUSSIAN 0.1
#define C_UNIFORM_WIND 1.
#define XC_GAUSSIAN 25.0
#define YC_GAUSSIAN 25.0
#define ZC_GAUSSIAN 25.0

real YAKL_INLINE gaussian(real x)                 { return A_GAUSSIAN * exp(-AX_GAUSSIAN * pow(x-XC_GAUSSIAN,2.)); }
real YAKL_INLINE gaussian(real x, real y)         { return A_GAUSSIAN * exp(-AX_GAUSSIAN * pow(x-XC_GAUSSIAN,2.)) * exp(-AY_GAUSSIAN * pow(y-YC_GAUSSIAN,2.)); }
real YAKL_INLINE gaussian(real x, real y, real z) { return A_GAUSSIAN * exp(-AX_GAUSSIAN * pow(x-XC_GAUSSIAN,2.)) * exp(-AY_GAUSSIAN * pow(y-YC_GAUSSIAN,2.)) * exp(-AZ_GAUSSIAN * pow(z-ZC_GAUSSIAN,2.)); }

// FIX THIS
real YAKL_INLINE vortices(real x)                 {return x;}
real YAKL_INLINE vortices(real x, real y)         {return x;}
real YAKL_INLINE vortices(real x, real y, real z) {return x;}

// FIX THIS
real YAKL_INLINE square(real x)                 {return x;}
real YAKL_INLINE square(real x, real y)         {return x;}
real YAKL_INLINE square(real x, real y, real z) {return x;}

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

template<uint ndims> void set_model_specific_params(Parameters &params)
{
  if (data_init_cond == DATA_INIT::GAUSSIAN)
  {
  params.xlen = 50.0;
  params.xc = 25.0;
  if (ndims>=2)
  {
    params.ylen = 50.0;
    params.yc = 25.0;
  }
  if (ndims ==3)
  {
  params.zlen = 50.0;
  params.zc = 25.0;
  }
  params.etime = 0.0;
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
   if (reconstruction_type == RECONSTRUCTION_TYPE::UFV && reconstruction_order == 1)
   { ufv1_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology, *this->geom); }

   if (reconstruction_type == RECONSTRUCTION_TYPE::CFV && reconstruction_order == 2)
   { cfv2_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, *this->topology, *this->geom);}

   this->diag_exchange->exchange_variable_set(diagnostic_vars);

   //compute D (qrecon U)
   if (differential_order == 2)
   { divergence2<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology); }
   if (differential_order == 4)
   { divergence4<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology); }
   if (differential_order == 6)
   { divergence6<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology); }
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

template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (VariableSet<1, nprog> &progvars, VariableSet<1, nconst> &constvars, Geometry<1, nquadx, nquady, nquadz> &geom)
{

    for (int i=0; i<nqdofs; i++)
    {
    if (data_init_cond == DATA_INIT::GAUSSIAN) {geom.set_1form_values(gaussian, progvars.fields_arr[0], i);}
    if (data_init_cond == DATA_INIT::VORTICES) {geom.set_1form_values(vortices, progvars.fields_arr[0], i);}
    if (data_init_cond == DATA_INIT::SQUARE)   {geom.set_1form_values(square,   progvars.fields_arr[0], i);}
    }
    if (wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_0form_values(uniform_x_wind,     constvars.fields_arr[0], 0);}
}

template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (VariableSet<2, nprog> &progvars, VariableSet<2, nconst> &constvars, Geometry<2, nquadx, nquady, nquadz> &geom)
{
    for (int i=0; i<nqdofs; i++)
    {
    if (data_init_cond == DATA_INIT::GAUSSIAN) {geom.set_2form_values(gaussian, progvars.fields_arr[0], i);}
    if (data_init_cond == DATA_INIT::VORTICES) {geom.set_2form_values(vortices, progvars.fields_arr[0], i);}
    if (data_init_cond == DATA_INIT::SQUARE)   {geom.set_2form_values(square,   progvars.fields_arr[0], i);}
    }
    if (wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_1form_values(uniform_x_wind,     constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (wind_init_cond == WIND_INIT::UNIFORM_Y    ) {geom.set_1form_values(uniform_y_wind,     constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (wind_init_cond == WIND_INIT::UNIFORM_XY   ) {geom.set_1form_values(uniform_xy_wind,    constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}
    if (wind_init_cond == WIND_INIT::DEFORMATIONAL) {geom.set_1form_values(deformational_wind, constvars.fields_arr[0], 0, LINE_INTEGRAL_TYPE::NORMAL);}

}


template <int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (VariableSet<3, nprog> &progvars, VariableSet<3, nconst> &constvars, Geometry<3, nquadx, nquady, nquadz> &geom)
{
    for (int i=0; i<nqdofs; i++)
    {
    if (data_init_cond == DATA_INIT::GAUSSIAN) {geom.set_3form_values(gaussian, progvars.fields_arr[0], i);}
    if (data_init_cond == DATA_INIT::SQUARE)   {geom.set_3form_values(square,   progvars.fields_arr[0], i);}
    }
    if (wind_init_cond == WIND_INIT::UNIFORM_X    ) {geom.set_2form_values(uniform_x_wind,     constvars.fields_arr[0], 0);}
    if (wind_init_cond == WIND_INIT::UNIFORM_Y    ) {geom.set_2form_values(uniform_y_wind,     constvars.fields_arr[0], 0);}
    if (wind_init_cond == WIND_INIT::UNIFORM_Z    ) {geom.set_2form_values(uniform_z_wind,     constvars.fields_arr[0], 0);}
    if (wind_init_cond == WIND_INIT::DEFORMATIONAL) {geom.set_2form_values(deformational_wind, constvars.fields_arr[0], 0);}
}

#endif
