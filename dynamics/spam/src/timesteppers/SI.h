#pragma once

#include "common.h"
#include "exchange.h"
#include "field_sets.h"
#include "model.h"
#include "topology.h"
#include <sstream>

constexpr int si_verbosity_level = 2;

template <uint nquad> class SITimeIntegrator {

public:
  int step;
  real tol;
  real avg_iters;
  SArray<real, 1, nquad> quad_pts;
  SArray<real, 1, nquad> quad_wts;
  FieldSet<nprognostic> *x;
  FieldSet<nprognostic> dx;
  FieldSet<nprognostic> xn;
  FieldSet<nprognostic> xm;
  Tendencies *tendencies;
  LinearSystem *linear_system;
  FieldSet<nconstant> *const_vars;
  FieldSet<nauxiliary> *auxiliary_vars;
  ExchangeSet<nprognostic> *x_exchange;

  bool is_initialized;
  SITimeIntegrator();
  SITimeIntegrator(const SITimeIntegrator &) = delete;
  SITimeIntegrator &operator=(const SITimeIntegrator &) = delete;
  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts, FieldSet<nauxiliary> &auxiliarys,
                  ExchangeSet<nprognostic> &prog_exch);
  void stepForward(real dt);
};

template <uint nquad> SITimeIntegrator<nquad>::SITimeIntegrator() {
  this->is_initialized = false;
}

template <uint nquad>
void SITimeIntegrator<nquad>::initialize(ModelParameters &params,
                                         Tendencies &tend, LinearSystem &linsys,
                                         FieldSet<nprognostic> &xvars,
                                         FieldSet<nconstant> &consts,
                                         FieldSet<nauxiliary> &auxiliarys,
                                         ExchangeSet<nprognostic> &prog_exch) {

  set_ref_quad_pts_wts(this->quad_pts, this->quad_wts);

  this->dx.initialize(xvars, "dx");
  this->xn.initialize(xvars, "xn");
  this->xm.initialize(xvars, "xm");
  this->x = &xvars;
  this->tendencies = &tend;
  this->linear_system = &linsys;
  this->const_vars = &consts;
  this->auxiliary_vars = &auxiliarys;
  this->x_exchange = &prog_exch;

  this->tol = params.si_tolerance;
  this->step = 0;
  this->avg_iters = 0;

  this->is_initialized = true;
}

real norm(FieldSet<nprognostic> &x) {
  // note that this function assumes that x halos have been exchanged
  real accum = 0;
  for (auto f : x.fields_arr) {
    accum = std::max(yakl::intrinsics::maxval(yakl::intrinsics::abs(f.data)),
                     accum);
  }
  return accum;
}

template <uint nquad> void SITimeIntegrator<nquad>::stepForward(real dt) {

  this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                *this->auxiliary_vars, this->dx);
  this->xn.copy(*this->x);
  // store residual in xm
  this->xm.waxpbypcz(-1, 1, -dt, this->xn, *this->x, this->dx);
  this->x_exchange->exchange_variable_set(this->xm);

  int iter = 0;
  int maxiters = 50;

  real res_norm = norm(xm);
  real initial_norm = res_norm;
  bool converged = false;

  if (si_verbosity_level > 0) {
    std::stringstream msg;
    msg << "Starting Newton iteration, step = " << step
        << ", initial residual = " << res_norm;
    std::cout << msg.str() << std::endl;
  }
  while (true) {
    if (res_norm / initial_norm < this->tol) {
      converged = true;
      break;
    } else if (iter > maxiters) {
      break;
    }

    this->linear_system->solve(dt, this->xm, *this->const_vars,
                               *this->auxiliary_vars, this->dx);

    this->xn.waxpy(1, this->dx, this->xn);

    this->xm.waxpby(1 - this->quad_pts(0), this->quad_pts(0), *this->x,
                    this->xn);
    this->x_exchange->exchange_variable_set(this->xm);
    this->tendencies->compute_functional_derivatives(
        ADD_MODE::REPLACE, this->quad_wts(0), dt, *this->const_vars, this->xm,
        *this->auxiliary_vars);

    for (int m = 1; m < nquad; ++m) {
      this->xm.waxpby(1 - this->quad_pts(m), this->quad_pts(m), *this->x,
                      this->xn);
      this->x_exchange->exchange_variable_set(this->xm);
      this->tendencies->compute_functional_derivatives(
          ADD_MODE::ADD, this->quad_wts(m), dt, *this->const_vars, this->xm,
          *this->auxiliary_vars);
    }

    this->xm.waxpby(0.5_fp, 0.5_fp, *this->x, this->xn);
    this->x_exchange->exchange_variable_set(this->xm);

    this->tendencies->apply_symplectic(dt, *this->const_vars, this->xm,
                                       *this->auxiliary_vars, this->dx);

    // store residual in xm
    this->xm.waxpbypcz(-1, 1, -dt, this->xn, *this->x, this->dx);
    this->x_exchange->exchange_variable_set(this->xm);
    res_norm = norm(xm);

    iter++;

    if (si_verbosity_level > 1) {
      std::stringstream msg;
      msg << "Iter: " << iter << " "
          << " " << res_norm;
      std::cout << msg.str() << std::endl;
    }
  }

  this->avg_iters *= step;
  this->step += 1;
  this->avg_iters += iter;
  this->avg_iters /= step;

  if (si_verbosity_level > 0) {
    std::stringstream msg;
    if (converged) {
      msg << "Newton solve converged in " << iter << " iters.\n";
    } else {
      msg << "!!! Newton solve failed to converge in " << iter << " iters.\n";
      msg << "!!! Solver tolerance: " << this->tol << "\n";
      msg << "!!! Achieved tolerance: " << res_norm / initial_norm << "\n";
    }
    msg << "Iters avg: " << this->avg_iters;
    std::cout << msg.str() << std::endl;
  }

  this->x->copy(this->xn);
  this->x_exchange->exchange_variable_set(*this->x);
}
