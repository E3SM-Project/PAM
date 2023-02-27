#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"

class LSRKTimeIntegrator : public TimeIntegrator {

public:
  using TimeIntegrator::TimeIntegrator;

  FieldSet<nprognostic> dx;
  FieldSet<nprognostic> *x;
  Tendencies *tendencies;
  FieldSet<nconstant> *const_vars;
  FieldSet<nauxiliary> *auxiliary_vars;
  std::vector<real> rka;
  std::vector<real> rkb;
  int nstages;

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {
    this->dx.initialize(xvars, "dx");
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;
    this->is_ssp = false;
    this->is_initialized = true;

    if (tstype == "lsrk12") {
      nstages = 12;
      this->rka.resize(nstages);
      this->rkb.resize(nstages);
      this->rka = {0,
                   -0.0923311242368072,
                   -0.9441056581158819,
                   -4.3271273247576394,
                   -2.1557771329026072,
                   -0.9770727190189062,
                   -0.7581835342571139,
                   -1.7977525470825499,
                   -2.6915667972700770,
                   -4.6466798960268143,
                   -0.1539613783825189,
                   -0.5943293901830616};

      this->rkb = {0.0650008435125904, 0.0161459902249842, 0.5758627178358159,
                   0.1649758848361671, 0.3934619494248182, 0.0443509641602719,
                   0.2074504268408778, 0.6914247433015102, 0.3766646883450449,
                   0.0757190350155483, 0.2027862031054088, 0.2167029365631842};
    }
  }

  void step_forward(real dt) override {
    for (int stage = 0; stage < nstages; ++stage) {
      // dx *= rka[stage];
      dx.wscal(this->rka[stage], this->dx);
      this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                    *this->auxiliary_vars, this->dx,
                                    ADD_MODE::ADD);
      // x += rkb[stage] * dt * dx;
      this->x->waxpy(this->rkb[stage] * (-dt), this->dx, *this->x);
      this->x->exchange();
    }
  }
};
