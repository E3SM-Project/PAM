#pragma once

#include "common.h"
#include "field_sets.h"
#include "model.h"

namespace pamc {

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
  bool is_semi_implicit;
  std::string tstype;

  Tendencies *tendencies;

  FieldSet<nprognostic> *x;
  FieldSet<nconstant> *const_vars;
  FieldSet<nauxiliary> *auxiliary_vars;

  virtual void initialize(ModelParameters &params, Tendencies &tend,
                          FieldSet<nprognostic> &xvars,
                          FieldSet<nconstant> &consts,
                          FieldSet<nauxiliary> &auxiliarys) {
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;
  }

  virtual void step_forward(real dt) = 0;

  TimeIntegrator(const std::string tstype, bool is_semi_implicit = false)
      : tstype(tstype), is_semi_implicit(is_semi_implicit) {}
  virtual ~TimeIntegrator() = default;
};

class SemiImplicitTimeIntegrator : public TimeIntegrator {
public:
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

  SemiImplicitTimeIntegrator(const std::string tstype)
      : TimeIntegrator(tstype, true) {}

  void compute_discrete_gradient(real dt, FieldSet<nprognostic> &xn,
                                 FieldSet<nconstant> &const_vars,
                                 FieldSet<nauxiliary> &auxiliary_vars,
                                 FieldSet<nprognostic> &xm) {

    if (this->two_point_discrete_gradient) {
      this->tendencies->compute_two_point_discrete_gradient(
          dt, const_vars, *this->x, xn, auxiliary_vars);
    } else {
      xm.waxpby(1 - this->quad_pts[0], this->quad_pts[0], *this->x, xn);
      xm.exchange();
      this->tendencies->compute_functional_derivatives(
          dt, const_vars, xm, auxiliary_vars, this->quad_wts[0],
          ADD_MODE::REPLACE);
      for (int m = 1; m < nquad; ++m) {
        xm.waxpby(1 - this->quad_pts[m], this->quad_pts[m], *this->x, xn);
        xm.exchange();
        this->tendencies->compute_functional_derivatives(
            dt, const_vars, xm, auxiliary_vars, this->quad_wts[m],
            ADD_MODE::ADD);
      }
    }
  }

  virtual void initialize(ModelParameters &params, Tendencies &tend,
                          FieldSet<nprognostic> &xvars,
                          FieldSet<nconstant> &consts,
                          FieldSet<nauxiliary> &auxiliarys) override {
    TimeIntegrator::initialize(params, tend, xvars, consts, auxiliarys);

    this->tol = params.si_tolerance;
    this->two_point_discrete_gradient = params.si_two_point_discrete_gradient;
    this->monitor_convergence = params.si_monitor_convergence;
    this->verbosity_level = params.si_verbosity_level;
    this->max_iters = params.si_max_iters;
    this->nquad = params.si_nquad;

    this->quad_pts.resize(nquad);
    this->quad_wts.resize(nquad);
    set_ref_quad_pts_wts(this->quad_pts, this->quad_wts, nquad);
  }
};
} // namespace pamc
