

#ifndef _RKSIMPLE_H_
#define _RKSIMPLE_H_


#include "common.h"
#include "topology.h"
#include "variableset.h"
#include "tendencies.h"

template<int nstages> class RKSimpleTimeIntegrator {

public:

  SArray<real, nstages> stage_coeffs;
  VariableSet xtend;
  VariableSet xtemp;
  VariableSet x;
  Topology topology;
  Tendencies tendencies;
  VariableSet const_vars;
  VariableSet diagnostic_vars;

  void stepForward(real dt);
  void initialize(Tendencies &tend, VariableSet &xvars, Topology &topo, VariableSet &consts, VariableSet &diagnostics);
};

#endif
