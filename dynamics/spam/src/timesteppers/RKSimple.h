
#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"

#define NSTAGESMAX 6

class RKSimpleTimeIntegrator : public TimeIntegrator {

public:
  using TimeIntegrator::TimeIntegrator;

  SArray<real, 1, NSTAGESMAX> stage_coeffs;
  FieldSet<nprognostic> xtend;
  FieldSet<nprognostic> xtemp;
  FieldSet<nprognostic> *x;
  Tendencies *tendencies;
  FieldSet<nconstant> *const_vars;
  FieldSet<nauxiliary> *auxiliary_vars;
  uint nstages;

  void set_stage_coefficients(ModelParameters &params);

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {
    this->xtemp.initialize(xvars, "xtemp");
    this->xtend.initialize(xvars, "xtend");
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;
    set_stage_coefficients(params);
    this->is_initialized = true;
  }

  void stepForward(real dt) override {

    this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                  *this->auxiliary_vars, this->xtend);
    this->xtemp.waxpy(-1. * dt * this->stage_coeffs(0), this->xtend, *this->x);

    for (int i = 1; i < nstages; i++) {
      this->xtemp.exchange();
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xtemp,
                                    *this->auxiliary_vars, this->xtend);
      this->xtemp.waxpy(-1. * dt * this->stage_coeffs(i), this->xtend,
                        *this->x);
    }
    // THIS COPY CAN BE AVOIDED IF WE ARE CLEVER ie swap x and xtemp
    // Would require being careful with IO also?
    this->x->copy(this->xtemp);
    this->x->exchange();
  }
};

void RKSimpleTimeIntegrator::set_stage_coefficients(ModelParameters &params) {
  if (params.tstype == "kgrk4") {
    this->nstages = 4;
    this->stage_coeffs(0) = 1. / 4.;
    this->stage_coeffs(1) = 1. / 3.;
    this->stage_coeffs(2) = 1. / 2.;
    this->stage_coeffs(3) = 1.;
  }
}
