#ifndef _ADVECTION_H_
#define _ADVECTION_H_


#include "common.h"
#include "topology.h"
#include "variable_sets.h"
#include "divergence.h"
#include "finitevolume.h"
#include "geometry.h"
#include "util.h"

// Initial conditions related variables and functions

real xc, yc, zc; // UNIFORM RECT SPECIFIC...
real xlen, ylen, zlen; // UNIFORM RECT SPECIFIC...

real A = 1.;
real ax = 0.1;
real ay = 0.1;
real az = 0.1;
real c = 1.;

enum class DATA_INIT { GAUSSIAN, VORTICES, SQUARE };
enum class WIND_INIT { UNIFORM_X, UNIFORM_Y, UNIFORM_Z, DEFORMATIONAL };

// MAYBE OVERLOAD STUFF?
// + DO MORE CLEVER STUFF IN GEOMETRY?

// HOW WE DO PROPERLY TREAT THIS
// gaussian sees y=0,z=0 in 1D, z=0 in 2D
//real YAKL_INLINE gaussian(real x, real y, real z) { return A * exp(-ax * pow(x-xc,2.)); }
real YAKL_INLINE gaussian(real x, real y, real z) { return A * exp(-ax * pow(x-xc,2.)) * exp(-ay * pow(y-yc,2.)); }

// FIX THIS
real YAKL_INLINE vortices(real x, real y, real z) {
  return x;
}

// FIX THIS
real YAKL_INLINE square(real x, real y, real z) {
  return x;
}

real YAKL_INLINE uniform_x_wind(real x, real y, real z) {
  return c;
}

real YAKL_INLINE uniform_y_wind(real x, real y, real z) {
  return x;
}

real YAKL_INLINE uniform_z_wind(real x, real y, real z) {
  return x;
}

real YAKL_INLINE deformational_wind(real x, real y, real z) {
  return x;
}



// *********** COMPILE TIME CONSTANTS ************** //
// EVENTUALLY THESE SHOULD ALL BE COMPILE TIME/PRE-PROCESSOR FLAGS...
// WITH REASONABLE DEFAULTS!

// Initial conditions
DATA_INIT data_init_cond = DATA_INIT::GAUSSIAN;
WIND_INIT wind_init_cond = WIND_INIT::UNIFORM_X;

// Number of variables
uint constexpr nprognostic = 1;
uint constexpr nconstant = 1;
uint constexpr ndiagnostic = 1;
uint constexpr nqdofs = 1;


// *******   Tendencies   ***********//

// EVENTUALLY NEED TO FIGURE OUT HOW TO MAKE THIS SWITCHABLE/SUBCLASSABLE...

template <uint ndims, uint nprog, uint nconst, uint ndiag> class Tendencies {
public:

  const Topology<ndims> *topology;
  ExchangeSet<ndims, ndiag> *diag_exchange;

  bool is_initialized;

   Tendencies() {
     this->is_initialized = false;
     std::cout << "CREATED TENDENCIES\n";
   }

  void initialize(const Topology<ndims> &topo, ExchangeSet<ndims, ndiag> &diag_exchange)
  {
    this->topology = &topo;
    this->diag_exchange = &diag_exchange;
    this->is_initialized = true;
  }

  void compute_rhs(const VariableSet<ndims, nconst> &const_vars, VariableSet<ndims, nprog> &x, VariableSet<ndims, ndiag> &diagnostic_vars, VariableSet<ndims, nprog> &xtend)
  {
    std::cout << "adv tend\n";
// EVENTUALLY HERE WE SHOULD TAKE RECON TYPE, RECON ORDER AND DIFFERENTIAL ORDER INTO ACCOUNT!

   //compute reconstructions
   //ufv1_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology);
   cfv2_recon<ndims, nqdofs>(diagnostic_vars.fields_arr[0].data, x.fields_arr[0].data, *this->topology);
   this->diag_exchange->exchange_variable_set(diagnostic_vars);

   //compute D (qrecon U)
   divergence2<ndims, nqdofs>(xtend.fields_arr[0].data, diagnostic_vars.fields_arr[0].data, const_vars.fields_arr[0].data, *this->topology);
 }

};

template <uint ndims, uint nprog, uint nconst, uint ndiag> class AdvectionTendencies : public Tendencies<ndims, nprog, nconst, ndiag> {

public:

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

template <int ndims, int nprog, int nconst, int ndiag, int nquadx, int nquady, int nquadz> void set_initial_conditions (VariableSet<ndims, nprog> &progvars, VariableSet<ndims, nconst> &constvars, Geometry<ndims, nquadx, nquady, nquadz> &geom)
{

// set advected quantity
for (int i=0; i<nqdofs; i++)
{
if (data_init_cond == DATA_INIT::GAUSSIAN) { set_n_form<ndims, nquadx, nquady, nquadz>(gaussian, progvars.fields_arr[0], geom, i); }
if (data_init_cond == DATA_INIT::VORTICES) { set_n_form<ndims, nquadx, nquady, nquadz>(vortices, progvars.fields_arr[0], geom, i); }
if (data_init_cond == DATA_INIT::SQUARE) { set_n_form<ndims, nquadx, nquady, nquadz>(square, progvars.fields_arr[0], geom, i); }
}

// set wind
if (wind_init_cond == WIND_INIT::UNIFORM_X) { set_n_minus_1_form<ndims, nquadx, nquady, nquadz>(uniform_x_wind, constvars.fields_arr[0], geom, 0); }
if (wind_init_cond == WIND_INIT::UNIFORM_Y) { set_n_minus_1_form<ndims, nquadx, nquady, nquadz>(uniform_y_wind, constvars.fields_arr[0], geom, 0); }
if (wind_init_cond == WIND_INIT::UNIFORM_Z) { set_n_minus_1_form<ndims, nquadx, nquady, nquadz>(uniform_z_wind, constvars.fields_arr[0], geom, 0); }
if (wind_init_cond == WIND_INIT::DEFORMATIONAL) { set_n_minus_1_form<ndims, nquadx, nquady, nquadz>(deformational_wind, constvars.fields_arr[0], geom, 0); }

}




#endif
