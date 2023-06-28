#pragma once

#include "common.h"
#include "exchange.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"
#include <sstream>

class SIFixedTimeIntegrator : public SemiImplicitTimeIntegrator {

public:
  using SemiImplicitTimeIntegrator::SemiImplicitTimeIntegrator;
  int step;
  real avg_iters;
  FieldSet<nprognostic> dx;
  FieldSet<nprognostic> xn;
  FieldSet<nprognostic> xm;

  void initialize(ModelParameters &params, Tendencies &tend,
                  LinearSystem &linsys, FieldSet<nprognostic> &xvars,
                  FieldSet<nconstant> &consts,
                  FieldSet<nauxiliary> &auxiliarys) override {

    SemiImplicitTimeIntegrator::initialize(params, tend, linsys, xvars, consts,
                                           auxiliarys);

    this->dx.initialize(xvars, "dx");
    this->xn.initialize(xvars, "xn");
    this->xm.initialize(xvars, "xm");

    this->step = 0;
    this->avg_iters = 0;

    this->is_initialized = true;
  }

  // evaluates -dt * J((x^n + x^v) / 2) dHtilde/dx(x^n, x^v), stores in dx
  void evaluate_fixed_point_rhs(real dt) {
    compute_discrete_gradient(dt, this->xn, *this->const_vars,
                              *this->auxiliary_vars, this->xm);

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
    if (monitor_convergence > 0) {
      this->dx.exchange();
      initial_res_norm = dt * norm(this->dx);
    }

    if (verbosity_level > 0) {
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

      if (iter >= max_iters) {
        break;
      }

      evaluate_fixed_point_rhs(dt);

      if (monitor_convergence > 1) {
        this->xm.waxpbypcz(1, -1, dt, this->xn, *this->x, this->dx);
        this->xm.exchange();
        res_norm = norm(xm);

        if (res_norm / initial_res_norm < this->tol) {
          converged = true;
        }
      }

      if (verbosity_level > 1) {
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

    if (monitor_convergence > 0) {
      evaluate_fixed_point_rhs(dt);

      this->xm.waxpbypcz(1, -1, dt, this->xn, *this->x, this->dx);
      this->xm.exchange();
      res_norm = norm(xm);

      if (res_norm / initial_res_norm < this->tol) {
        converged = true;
      }
    }

    if (verbosity_level > 0) {
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
