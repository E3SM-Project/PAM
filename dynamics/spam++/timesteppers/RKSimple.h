

#ifndef _RKSIMPLE_H_
#define _RKSIMPLE_H_


#include "common.h"
#include "topology.h"
#include "variableset.h"
#include "tendencies.h"
#include "STDLIB.h"

class RKSimpleTimeIntegrator {

public:

  int nstages;
  std::vector<real> stage_coeffs;
  VariableSet xtend;
  VariableSet xtemp;
  VariableSet x;
  Topology topology;
  Tendencies tendencies;
  VariableSet const_vars;
  VariableSet diagnostic_vars;

  real dt;

  void stepForward();
  void initialize(Tendencies &tend, VariableSet &xvars, Topology &topo, VariableSet &consts, VariableSet &diagnostics, real dtt);

};

class KG4 : RKSimpleTimeIntegrator {};

#endif
