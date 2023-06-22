#pragma once

#include "common.h"
#include "exchange.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"
#include <sstream>

template <uint nquad> class SIFixedTimeIntegrator : public TimeIntegrator {

public:
  using TimeIntegrator::TimeIntegrator;
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
  FieldSet<nconstant> *const_vars;
  FieldSet<nauxiliary> *auxiliary_vars;

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {

    set_ref_quad_pts_wts(this->quad_pts, this->quad_wts);

    this->dx.initialize(xvars, "dx");
    this->xn.initialize(xvars, "xn");
    this->xm.initialize(xvars, "xm");
    this->x = &xvars;
    this->tendencies = &tend;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;

    this->tol = params.si_tolerance;
    this->step = 0;
    this->avg_iters = 0;

    this->is_semi_implicit = true;
    this->is_initialized = true;
  }

  // evaluates -dt * J((x^n + x^v) / 2) dHtilde/dx(x^n, x^v), stores in dx
  void evaluate_fixed_point_rhs(real dt) {
    this->xm.waxpby(1 - this->quad_pts(0), this->quad_pts(0), *this->x,
                    this->xn);
    this->xm.exchange();
    this->tendencies->compute_functional_derivatives(
        dt, *this->const_vars, this->xm, *this->auxiliary_vars,
        this->quad_wts(0), ADD_MODE::REPLACE);

    for (int m = 1; m < nquad; ++m) {
      this->xm.waxpby(1 - this->quad_pts(m), this->quad_pts(m), *this->x,
                      this->xn);
      this->xm.exchange();
      this->tendencies->compute_functional_derivatives(
          dt, *this->const_vars, this->xm, *this->auxiliary_vars,
          this->quad_wts(m), ADD_MODE::ADD);
    }

    this->xm.waxpby(0.5_fp, 0.5_fp, *this->x, this->xn);
    this->xm.exchange();

    this->tendencies->apply_symplectic(dt, *this->const_vars, this->xm,
                                       *this->auxiliary_vars, this->dx,
                                       ADD_MODE::REPLACE);
    this->tendencies->add_pressure_perturbation(
        dt, *this->const_vars, this->xm, *this->auxiliary_vars, this->dx);
  }

  void step_forward(real dt) override {

    this->xn.copy(*this->x);
    this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                  *this->auxiliary_vars, this->dx);

    real res_norm;
    real initial_res_norm;
    if (si_monitor_convergence > 0) {
      this->dx.exchange();
      initial_res_norm = dt * norm(this->dx);
    }

    if (si_verbosity_level > 0) {
      std::stringstream msg;
      msg << "Starting fixed-point iteration, step = " << step
          << ", initial residual = " << initial_res_norm;
      std::cout << msg.str() << std::endl;
    }

    bool converged = false;
    int iter = 0;
    while (true) {
      this->xn.waxpby(1, -dt, *this->x, this->dx);

      if (converged) {
        break;
      }

      iter++;

      if (iter >= si_max_iters) {
        break;
      }

      evaluate_fixed_point_rhs(dt);

      if (si_monitor_convergence > 1) {
        this->xm.waxpbypcz(1, -1, dt, this->xn, *this->x, this->dx);
        this->xm.exchange();
        res_norm = norm(xm);

        if (res_norm / initial_res_norm < this->tol) {
          converged = true;
        }
      }

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
      evaluate_fixed_point_rhs(dt);

      this->xm.waxpbypcz(1, -1, dt, this->xn, *this->x, this->dx);
      this->xm.exchange();
      res_norm = norm(xm);

      std::stringstream msg;
      if (converged) {
        msg << "Fixed-point iteration converged in " << iter << " iters.\n";
      } else {
        msg << "!!! Fixed-point iteration failed to converge in " << iter
            << " iters.\n";
        msg << "!!! Solver tolerance: " << this->tol << "\n";
        msg << "!!! Achieved tolerance: " << res_norm / initial_res_norm
            << "\n";
      }
      msg << "Iters avg: " << this->avg_iters;
      std::cout << msg.str() << std::endl;
    }

    this->x->copy(this->xn);
    this->x->exchange();
  }
};