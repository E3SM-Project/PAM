
#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"

namespace pamc {

int constexpr NSTAGESMAX = 10;

class KGRKTimeIntegrator : public TimeIntegrator {

public:
  using TimeIntegrator::TimeIntegrator;

  SArray<real, 1, NSTAGESMAX> stage_coeffs;
  FieldSet<nprognostic> xtend;
  FieldSet<nprognostic> xtemp;
  uint nstages;

  void set_stage_coefficients(ModelParameters &params);

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {
    TimeIntegrator::initialize(params, tend, linsys, xvars, consts, auxiliarys);

    this->xtemp.initialize(xvars, "xtemp");
    this->xtend.initialize(xvars, "xtend");
    set_stage_coefficients(params);
    this->is_initialized = true;
  }

  void step_forward(real dt) override {

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

void KGRKTimeIntegrator::set_stage_coefficients(ModelParameters &params) {
  if (params.tstype == "kgrk2") {
    this->nstages = 2;
    this->stage_coeffs(0) = 1. / 2.;
    this->stage_coeffs(1) = 1. / 1.;
  } else if (params.tstype == "kgrk3") {
    this->nstages = 3;
    this->stage_coeffs(0) = 1. / 3.;
    this->stage_coeffs(1) = 1. / 2.;
    this->stage_coeffs(2) = 1. / 1.;
  } else if (params.tstype == "kgrk4") {
    this->nstages = 4;
    this->stage_coeffs(0) = 1. / 4.;
    this->stage_coeffs(1) = 1. / 3.;
    this->stage_coeffs(2) = 1. / 2.;
    this->stage_coeffs(3) = 1. / 1.;
  } else if (params.tstype == "kgrk5") {
    this->nstages = 5;
    this->stage_coeffs(0) = 1. / 5.;
    this->stage_coeffs(1) = 1. / 5.;
    this->stage_coeffs(2) = 1. / 3.;
    this->stage_coeffs(3) = 1. / 2.;
    this->stage_coeffs(4) = 1. / 1.;
  } else if (params.tstype == "kgrk6") {
    this->nstages = 6;
    this->stage_coeffs(0) = 1. / 6.;
    this->stage_coeffs(1) = 2. / 15.;
    this->stage_coeffs(2) = 1. / 4.;
    this->stage_coeffs(3) = 1. / 3.;
    this->stage_coeffs(4) = 1. / 2.;
    this->stage_coeffs(5) = 1. / 1.;
  } else if (params.tstype == "kgrk7") {
    this->nstages = 7;
    this->stage_coeffs(0) = 1. / 7.;
    this->stage_coeffs(1) = 2. / 21.;
    this->stage_coeffs(2) = 1. / 5.;
    this->stage_coeffs(3) = 8. / 35.;
    this->stage_coeffs(4) = 1. / 3.;
    this->stage_coeffs(5) = 1. / 2.;
    this->stage_coeffs(6) = 1. / 1.;
  } else if (params.tstype == "kgrk8") {
    this->nstages = 8;
    this->stage_coeffs(0) = 1. / 8.;
    this->stage_coeffs(1) = 1. / 14.;
    this->stage_coeffs(2) = 1. / 6.;
    this->stage_coeffs(3) = 1. / 6.;
    this->stage_coeffs(4) = 1. / 4.;
    this->stage_coeffs(5) = 1. / 3.;
    this->stage_coeffs(6) = 1. / 2.;
    this->stage_coeffs(7) = 1. / 1.;
  } else if (params.tstype == "kgrk9") {
    this->nstages = 9;
    this->stage_coeffs(0) = 1. / 9.;
    this->stage_coeffs(1) = 1. / 18.;
    this->stage_coeffs(2) = 1. / 7.;
    this->stage_coeffs(3) = 8. / 63.;
    this->stage_coeffs(4) = 1. / 5.;
    this->stage_coeffs(5) = 5. / 21.;
    this->stage_coeffs(6) = 1. / 3.;
    this->stage_coeffs(7) = 1. / 2.;
    this->stage_coeffs(8) = 1. / 1.;
  } else if (params.tstype == "kgrk10") {
    this->nstages = 10;
    this->stage_coeffs(0) = 1. / 10.;
    this->stage_coeffs(1) = 2. / 45.;
    this->stage_coeffs(2) = 1. / 8.;
    this->stage_coeffs(3) = 1. / 10.;
    this->stage_coeffs(4) = 1. / 6.;
    this->stage_coeffs(5) = 9. / 50.;
    this->stage_coeffs(6) = 1. / 4.;
    this->stage_coeffs(7) = 1. / 3.;
    this->stage_coeffs(8) = 1. / 2.;
    this->stage_coeffs(9) = 1. / 1.;
  }
}
} // namespace pamc
