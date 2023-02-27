#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"

class SSPKKTimeIntegrator : public TimeIntegrator {

public:
  using TimeIntegrator::TimeIntegrator;

  FieldSet<nprognostic> F1;
  FieldSet<nprognostic> x1;
  FieldSet<nprognostic> F2;
  FieldSet<nprognostic> x2;
  FieldSet<nprognostic> x3;
  FieldSet<nprognostic> F3;
  FieldSet<nprognostic> *x;
  Tendencies *tendencies;
  FieldSet<nconstant> *const_vars;
  FieldSet<nauxiliary> *auxiliary_vars;

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {
    this->x1.initialize(xvars, "x1");
    this->F1.initialize(xvars, "F1");
    this->x2.initialize(xvars, "x2");
    this->F2.initialize(xvars, "F2");
    if (tstype == "ssprk3") {
      this->x3.initialize(xvars, "x3");
      this->F3.initialize(xvars, "F3");
    }
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;
    this->is_ssp = true;
    this->is_initialized = true;
  }

  void step_forward(real dt) override {

    this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                  *this->auxiliary_vars, this->F1);
    this->x1.waxpy(-1. * dt, this->F1, *this->x);
    this->x1.exchange();
    this->tendencies->compute_rhs(dt, *this->const_vars, this->x1,
                                  *this->auxiliary_vars, this->F2);

    if (tstype == "ssprk2") {
      this->x2.waxpbypcz(0.5, 0.5, -0.5 * dt, *this->x, this->x1, this->F2);
      this->x->copy(this->x2);
    }

    if (tstype == "ssprk3") {
      this->x2.waxpbypcz(0.75, 0.25, -0.25 * dt, *this->x, this->x1, this->F2);
      this->x2.exchange();
      this->tendencies->compute_rhs(dt, *this->const_vars, this->x2,
                                    *this->auxiliary_vars, this->F3);
      this->x3.waxpbypcz(1. / 3., 2. / 3., -2. / 3. * dt, *this->x, this->x2,
                         this->F3);
      this->x->copy(this->x3);
    }

    this->x->exchange();
  }
};
