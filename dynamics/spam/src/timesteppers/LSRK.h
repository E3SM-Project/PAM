#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"

namespace pamc {

class LSRKTimeIntegrator : public TimeIntegrator {

public:
  using TimeIntegrator::TimeIntegrator;

  FieldSet<nprognostic> dx;
  std::vector<real> rka;
  std::vector<real> rkb;
  int nstages;

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {
    TimeIntegrator::initialize(params, tend, linsys, xvars, consts, auxiliarys);

    this->dx.initialize(xvars, "dx");
    this->is_ssp = false;
    this->is_initialized = true;

    if (tstype == "lsrk5") {
      this->nstages = 5;
      this->rka.resize(nstages);
      this->rkb.resize(nstages);
      this->rka = {
          0., -567301805773. / 1357537059087., -2404267990393. / 2016746695238.,
          -3550918686646. / 2091501179385., -1275806237668. / 842570457699.};
      this->rkb = {
          1432997174477. / 9575080441755., 5161836677717. / 13612068292357.,
          1720146321549. / 2090206949498., 3134564353537. / 4481467310338.,
          2277821191437. / 14882151754819.};

    } else if (tstype == "lsrk12") {
      this->nstages = 12;
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
    } else if (tstype == "lsrk13") {
      this->nstages = 13;
      this->rka.resize(nstages);
      this->rkb.resize(nstages);
      this->rka = {0,
                   -0.6160178650170565,
                   -0.4449487060774118,
                   -1.0952033345276178,
                   -1.2256030785959187,
                   -0.2740182222332805,
                   -0.0411952089052647,
                   -0.1797084899153560,
                   -1.1771530652064288,
                   -0.4078831463120878,
                   -0.8295636426191777,
                   -4.7895970584252288,
                   -0.6606671432964504};
      this->rkb = {0.0271990297818803, 0.1772488819905108, 0.0378528418949694,
                   0.6086431830142991, 0.2154313974316100, 0.2066152563885843,
                   0.0415864076069797, 0.0219891884310925, 0.9893081222650993,
                   0.0063199019859826, 0.3749640721105318, 1.6080235151003195,
                   0.0961209123818189};
    } else if (tstype == "lsrk14") {
      this->nstages = 14;
      this->rka.resize(nstages);
      this->rkb.resize(nstages);
      this->rka = {0,
                   -0.7188012108672410,
                   -0.7785331173421570,
                   -0.0053282796654044,
                   -0.8552979934029281,
                   -3.9564138245774565,
                   -1.5780575380587385,
                   -2.0837094552574054,
                   -0.7483334182761610,
                   -0.7032861106563359,
                   0.0013917096117681,
                   -0.0932075369637460,
                   -0.9514200470875948,
                   -7.1151571693922548};

      this->rkb = {0.0367762454319673, 0.3136296607553959, 0.1531848691869027,
                   0.0030097086818182, 0.3326293790646110, 0.2440251405350864,
                   0.3718879239592277, 0.6204126221582444, 0.1524043173028741,
                   0.0760894927419266, 0.0077604214040978, 0.0024647284755382,
                   0.0780348340049386, 5.5059777270269628};
    } else {
      throw std::runtime_error("Unknown lsrk method");
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
} // namespace pamc
