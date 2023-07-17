#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"

namespace pamc {

class SSPRKTimeIntegrator : public TimeIntegrator {

public:
  using TimeIntegrator::TimeIntegrator;

  int nstages;

  FieldSet<nprognostic> F;
  std::vector<FieldSet<nprognostic>> xstage;

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {
    TimeIntegrator::initialize(params, tend, linsys, xvars, consts, auxiliarys);

    if (tstype == "ssprk2") {
      this->nstages = 2;
    } else if (tstype == "ssprk3") {
      this->nstages = 3;
    } else if (tstype == "ssprk34") {
      this->nstages = 4;
    }
    this->xstage.resize(nstages);

    for (int stage = 0; stage < nstages; ++stage) {
      this->xstage[stage].initialize(xvars, "x" + std::to_string(stage));
    }
    this->F.initialize(xvars, "F");

    this->is_ssp = true;
    this->is_initialized = true;
  }

  void step_forward(real dt) override {

    if (tstype == "ssprk2") {
      this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                    *this->auxiliary_vars, this->F);
      this->xstage[0].waxpy(-1. * dt, this->F, *this->x);
      this->xstage[0].exchange();

      this->tendencies->compute_rhs(dt, *this->const_vars, this->xstage[0],
                                    *this->auxiliary_vars, this->F);

      this->xstage[1].waxpbypcz(0.5, 0.5, -0.5 * dt, *this->x, this->xstage[0],
                                this->F);
      this->x->copy(this->xstage[1]);
    }

    if (tstype == "ssprk3") {
      this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                    *this->auxiliary_vars, this->F);
      this->xstage[0].waxpy(-1. * dt, this->F, *this->x);
      this->xstage[0].exchange();

      this->tendencies->compute_rhs(dt, *this->const_vars, this->xstage[0],
                                    *this->auxiliary_vars, this->F);

      this->xstage[1].waxpbypcz(0.75, 0.25, -0.25 * dt, *this->x,
                                this->xstage[0], this->F);
      this->xstage[1].exchange();
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xstage[1],
                                    *this->auxiliary_vars, this->F);
      this->xstage[2].waxpbypcz(1. / 3., 2. / 3., -2. / 3. * dt, *this->x,
                                this->xstage[1], this->F);
      this->x->copy(this->xstage[2]);
    }

    if (tstype == "ssprk34") {
      this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                    *this->auxiliary_vars, this->F);

      this->xstage[0].waxpy(1. / 2. * dt, this->F, *this->x);
      this->xstage[0].exchange();
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xstage[0],
                                    *this->auxiliary_vars, this->F);

      this->xstage[1].waxpy(-1. / 2. * dt, this->F, this->xstage[0]);
      this->xstage[1].exchange();
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xstage[1],
                                    *this->auxiliary_vars, this->F);

      this->xstage[2].waxpbypcz(2. / 3., 1. / 3., -1. / 6. * dt, *this->x,
                                this->xstage[1], this->F);
      this->xstage[2].exchange();
      this->tendencies->compute_rhs(dt, *this->const_vars, this->xstage[2],
                                    *this->auxiliary_vars, this->F);

      this->xstage[3].waxpy(-1. / 2. * dt, this->F, this->xstage[2]);
      this->xstage[3].exchange();

      this->x->copy(this->xstage[3]);
    }

    this->x->exchange();
  }
};
} // namespace pamc
