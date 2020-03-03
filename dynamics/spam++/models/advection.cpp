
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
void ModelTendencies::compute_rhs(VariableSet &const_vars, VariableSet &x, VariableSet &diagnostic_vars, VariableSet &xtend, Topology &topology) {

// EVENTUALLY HERE WE SHOULD TAKE RECON TYPE, RECON ORDER AND DIFFERENTIAL ORDER INTO ACCOUNT!

  //compute reconstructions
  fv1_recon( diagnostic_vars.field_arr[0].data, x.field_arr[0].data, nqdofs, topology);
  diagnostic_vars.exchange();

  //compute D (qrecon U)
  divergence2(xtend.q, diagnostic_vars.field_arr[0].data, const_vars.field_arr[0], nqdofs, topology);
};


// *******   VariableSets   ***********//

  void DiagnosticVars::create(Topology &topology, bool create_exchange = true) {
// MOVE MOST OF THIS LOGIC INTO A VARIABLESET CONSTRUCTOR?
    baseName = "";
    num_fields = 1;
    //CREATE FIELDS_ARR and EXCHANGE_ARR
    if (create_exchange) {}



    if (ndims == 1) {
      field_arr[0].intialize(topology, "qrecon", nqdofs, 0);
      if (create_exchange) { exchange_arr[0].initialize(topology, nqdofs, 0); }
    }
    if (ndims == 2) {
      field_arr[0].intialize(topology, "qrecon", 0, nqdofs, 0);
      if (create_exchange) { exchange_arr[0].initialize(topology, 0, nqdofs, 0); }
    }
    if (ndims == 3) {
      field_arr[0].intialize(topology, "qrecon", 0, 0, nqdofs, 0);
      if (create_exchange) { exchange_arr[0].initialize(topology, 0, 0, nqdofs, 0); }
    }
    }


  void ConstantVars::create(Topology &topology, bool create_exchange = true) {
    if (ndims == 1) {
      field_arr[0].intialize(topology, "u", 1, 0);
      if (create_exchange) { exchange_arr[0].initialize(topology, 1, 0); }
    }
    if (ndims == 2) {
      field_arr[0].intialize(topology, "u", 0, 1, 0);
      if (create_exchange) { exchange_arr[0].initialize(topology, 0, 1, 0); }
    }
    if (ndims == 3) {
      field_arr[0].intialize(topology, "u", 0, 0, 1, 0);
      if (create_exchange) { exchange_arr[0].initialize(topology, 0, 0, 1, 0); }
    }
  }



  void PrognosticVars::create(Topology &topology, Params &params, bool create_exchange = true) {
    if (ndims == 1) {
      field_arr[0].intialize(topology, "q", 0, nqdofs);
      if (create_exchange) { exchange_arr[0].initialize(topology, 0, nqdofs); }
    }
    if (ndims == 2) {
      field_arr[0].intialize(topology, "q", 0, 0, nqdofs);
      if (create_exchange) { exchange_arr[0].initialize(topology, 0, 0, nqdofs); }
    }
    if (ndims == 3) {
      field_arr[0].intialize(topology, "q", 0, 0, 0, nqdofs);
      if (create_exchange) { exchange_arr[0].initialize(topology, 0, 0, 0, nqdofs); }
    }
  }




  // *******   Initial Conditions   ***********//

void set_initial_conditions(VariableSet &constant_vars, VariableSet &prognostic_vars)

{

// set advected quantity
if (data_init_cond == DATA_INIT_GAUSSIAN) { set_n_form(gaussian, prognostic_vars.field_arr(0)); }
if (data_init_cond == DATA_INIT_VORTICES) { set_n_form(vortices, prognostic_vars.field_arr(0)); }
if (data_init_cond == DATA_INIT_SQUARE) { set_n_form(square, prognostic_vars.field_arr(0)); }

// set wind
if (wind_init_cond == WIND_INIT_UNIFORM_X) { set_n_minus_1_form(uniform_x_wind, constant_vars.field_arr(0)); }
if (wind_init_cond == WIND_INIT_UNIFORM_Y) { set_n_minus_1_form(uniform_y_wind, constant_vars.field_arr(0)); }
if (wind_init_cond == WIND_INIT_UNIFORM_Z) { set_n_minus_1_form(uniform_z_wind, constant_vars.field_arr(0)); }
if (wind_init_cond == WIND_INIT_DEFORMATIONAL) { set_n_minus_1_form(deformational_wind, constant_vars.field_arr(0); }

// Do a boundary exchange
prognostic_vars.exchange();
constant_vars.exchange();

}
