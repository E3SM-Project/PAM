

#include "RKSimple.h"

  void RKSimpleTimeIntegrator::initialize(Tendencies &tend, VariableSet &xvars, Topology &topo, VariableSet &consts, VariableSet &diagnostics, real dt) {
    xvars.clone(xtemp, "xtemp", true);
    xvars.clone(xtend, "xtend", false);
    x = xvars;

    topology = topo;
    tendencies = tend;
    const_vars = consts;
    diagnostic_vars = diagnostics;

    dt = dtt;
  }

  void RKSimpleTimeIntegrator::stepForward() {

    tendencies.compute_rhs(const_vars, x, diagnostic_vars, xtend, topology, params);
    xtemp.waxpy(dt * stage_coeffs(0), xtend, x);

    for (int i=1; i<nstages; i++)
    {
      xtemp.exchange();
      tendencies.compute_rhs(constant_vars, xtemp, diagnostic_vars, xtend, topology, params);
      xtemp.waxpy(dt * stage_coeffs(i), xtend, x);
    }
    // THIS COPY CAN BE AVOIDED IF WE ARE CLEVER ie swap x and xtemp
    // Would require being careful with IO also?
    x.copy(xtemp);
    x.exchange();
  }


  void KG4::initialize(Tendencies &tend, VariableSet &xvars, Topology &topo, VariableSet &consts, VariableSet &diagnostics, Params &params) {
    super.initialize(tend, xvars, topo, consts, diagnostics, params);
    stage_coeffs = { 0.25, 1./3., 1./2., 1.};
    nstages = 4;
  }
