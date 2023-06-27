#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"

real norm(FieldSet<nprognostic> &x) {
  // note that this function assumes that x halos have been exchanged
  real accum = 0;
  for (auto f : x.fields_arr) {
    accum = std::max(yakl::intrinsics::maxval(yakl::intrinsics::abs(f.data)),
                     accum);
  }
  return accum;
}

class TimeIntegrator {
public:
  bool is_initialized = false;
  bool is_ssp = false;
  bool is_semi_implicit = false;
  std::string tstype;

  virtual void initialize(ModelParameters &params, Tendencies &tend,
                          LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                          FieldSet<nconstant> &consts,
                          FieldSet<nauxiliary> &auxiliarys) = 0;

  virtual void step_forward(real dt) = 0;

  TimeIntegrator(const std::string tstype) : tstype(tstype) {}
  virtual ~TimeIntegrator() = default;
};

class SemiImplicitTimeIntegrator : public TimeIntegrator {
public:
  using TimeIntegrator::TimeIntegrator;
  int monitor_convergence;
  // 0 = do not monitor (does si_max_iters iterations)
  // 1 = computes initial and final residual but still does si_max_iter
  // iterations 2 = iterates until convergence or si_max_iter is reached

  int verbosity_level;
  // 0 = do not print
  // 1 = print initial and final
  // 2 = print every iteration

  int max_iters;
  int nquad;

  bool two_point_discrete_gradient;

  real tol;
  std::vector<real> quad_pts;
  std::vector<real> quad_wts;

  virtual void initialize(ModelParameters &params, Tendencies &tend,
                          LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                          FieldSet<nconstant> &consts,
                          FieldSet<nauxiliary> &auxiliarys) override {

    this->tol = params.si_tolerance;
    this->two_point_discrete_gradient = params.si_two_point_discrete_gradient;
    this->monitor_convergence = params.si_monitor_convergence;
    this->verbosity_level = params.si_max_iters;
    this->max_iters = params.si_max_iters;
    this->nquad = params.si_nquad;

    this->quad_pts.resize(nquad);
    this->quad_wts.resize(nquad);
    set_ref_quad_pts_wts(this->quad_pts, this->quad_wts, nquad);

    this->is_semi_implicit = true;
  }
};
