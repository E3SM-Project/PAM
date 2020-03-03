
#include "advection.h"



// *******   Parameters   ***********//

// FIX THIS TO ACTUALLY READ A FILE, EVENTUALLY
void readParamsFile(std::string inFile) {
  // Topology
  nx = 50;
  ny = 1;
  nz = 1;

  // Time stepping
  dt = 1.0;
  Nsteps = 100;
  Nout = 5;

  // EVENTUALLY THIS SHOULD BE SET BY THE INITIAL CONDITION CHOICE!
  xlen = 50.0;
  ylen = 0.0;
  zlen = 0.0;
  xc = 25.0;
  yc = 0.0;
  zc = 0.0;
};


// *******   Tendencies   ***********//
void Tendencies::compute_rhs(VariableSet &const_vars, VariableSet &x, VariableSet &diagnostic_vars, VariableSet &xtend, Topology &topology) {

// EVENTUALLY HERE WE SHOULD TAKE RECON TYPE, RECON ORDER AND DIFFERENTIAL ORDER INTO ACCOUNT!

  //compute reconstructions
  fv1_recon<nqdofs>(diagnostic_vars.field_arr[0].data, x.field_arr[0].data, topology);
  diagnostic_vars.exchange();

  //compute D (qrecon U)
  divergence2<nqdofs>(xtend.field_arr[0], diagnostic_vars.field_arr[0].data, const_vars.field_arr[0], topology);
};


// *******   VariableSet Initialization   ***********//

void initialize_variables(Topology &topo, prognostic_vars, constant_vars, diagnostic_vars)
{
  prognostic_vars.initialize("x", true);
  constant_vars.initialize("cons", true);
  diagnostic_vars.initialize("diag", true);

  if (ndims == 1) {
    prognostic_vars.field_arr[0].intialize(topology, "q", 0, nqdofs, 0, 0);
    prognostic_vars.exchange_arr[0].initialize(topology, 0, nqdofs, 0, 0);
    constant_vars.field_arr[0].intialize(topology, "u", 1, 0, 0, 0);
    constant_vars.exchange_arr[0].initialize(topology, 1, 0, 0, 0);
    diagnostic_vars.field_arr[0].intialize(topology, "qrecon", nqdofs, 0, 0, 0);
    diagnostic_vars.exchange_arr[0].initialize(topology, nqdofs, 0, 0, 0);
  }

  if (ndims == 2) {
    prognostic_vars.field_arr[0].intialize(topology, "q", 0, 0, nqdofs, 0);
    prognostic_vars.exchange_arr[0].initialize(topology, 0, 0, nqdofs, 0);
    constant_vars.field_arr[0].intialize(topology, "u", 0, 1, 0, 0);
    constant_vars.exchange_arr[0].initialize(topology, 0, 1, 0, 0);
    diagnostic_vars.field_arr[0].intialize(topology, "qrecon", 0, nqdofs, 0, 0);
    diagnostic_vars.exchange_arr[0].initialize(topology, 0, nqdofs, 0, 0);
  }

  if (ndims == 3) {
    prognostic_vars.field_arr[0].intialize(topology, "q", 0, 0, 0, nqdofs);
    prognostic_vars.exchange_arr[0].initialize(topology, 0, 0, 0, nqdofs);
    constant_vars.field_arr[0].intialize(topology, "u", 0, 0, 1, 0);
    constant_vars.exchange_arr[0].initialize(topology, 0, 0, 1, 0);
    diagnostic_vars.field_arr[0].intialize(topology, "qrecon", 0, 0, nqdofs, 0);
    diagnostic_vars.exchange_arr[0].initialize(topology, 0, 0, nqdofs, 0);
  }

}



  // *******   Initial Conditions   ***********//

void set_initial_conditions(VariableSet &constant_vars, VariableSet &prognostic_vars)

{

// set advected quantity
for (int i=0; i<nqdofs; i++)
{
if (data_init_cond == DATA_INIT_GAUSSIAN) { set_n_form(gaussian, prognostic_vars.field_arr[0], i); }
if (data_init_cond == DATA_INIT_VORTICES) { set_n_form(vortices, prognostic_vars.field_arr[0], i); }
if (data_init_cond == DATA_INIT_SQUARE) { set_n_form(square, prognostic_vars.field_arr[0], i); }
}

// set wind
if (wind_init_cond == WIND_INIT_UNIFORM_X) { set_n_minus_1_form(uniform_x_wind, constant_vars.field_arr[0], 0); }
if (wind_init_cond == WIND_INIT_UNIFORM_Y) { set_n_minus_1_form(uniform_y_wind, constant_vars.field_arr[0], 0); }
if (wind_init_cond == WIND_INIT_UNIFORM_Z) { set_n_minus_1_form(uniform_z_wind, constant_vars.field_arr[0], 0); }
if (wind_init_cond == WIND_INIT_DEFORMATIONAL) { set_n_minus_1_form(deformational_wind, constant_vars.field_arr[0], 0); }

// Do a boundary exchange
prognostic_vars.exchange();
constant_vars.exchange();

}
