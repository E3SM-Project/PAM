#pragma once

#include "common.h"
#include "exchange.h"
#include "field_sets.h"
#include "model.h"
#include "time_integrator.h"
#include "topology.h"
#include <sstream>

template <uint nquad> class SINewtonTimeIntegrator : public TimeIntegrator {

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
  LinearSystem *linear_system;
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
    this->linear_system = &linsys;
    this->const_vars = &consts;
    this->auxiliary_vars = &auxiliarys;

    this->tol = params.si_tolerance;
    this->step = 0;
    this->avg_iters = 0;

    this->is_semi_implicit = true;
    this->is_initialized = true;
  }

  void step_forward(real dt) override {

    this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                  *this->auxiliary_vars, this->dx);
    this->xn.copy(*this->x);
    // store residual in xm
    this->xm.waxpbypcz(-1, 1, -dt, this->xn, *this->x, this->dx);
    this->xm.exchange();

    int iter = 0;
    int maxiters = 50;

    real res_norm;
    real initial_res_norm;
    if (si_monitor_convergence > 0) {
      res_norm = norm(xm);
      initial_res_norm = res_norm;
    }

    bool converged = false;

    if (si_verbosity_level > 0) {
      std::stringstream msg;
      msg << "Starting Newton iteration, step = " << step
          << ", initial residual = " << initial_res_norm;
      std::cout << msg.str() << std::endl;
    }
    while (true) {
      if (si_monitor_convergence > 1 &&
          res_norm / initial_res_norm < this->tol) {
        converged = true;
        break;
      }
      if (iter >= si_max_iters) {
        break;
      }

      this->linear_system->solve(dt, this->xm, *this->const_vars,
                                 *this->auxiliary_vars, this->dx);

      this->xn.waxpy(1, this->dx, this->xn);

      if (si_compute_functional_derivatives_quadrature) {
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
      } else {
        this->xn.exchange();
        this->tendencies->compute_functional_derivatives_two_point(
            dt, *this->const_vars, *this->x, this->xn, *this->auxiliary_vars);
      }

      this->xm.waxpby(0.5_fp, 0.5_fp, *this->x, this->xn);
      this->xm.exchange();

      this->tendencies->apply_symplectic(dt, *this->const_vars, this->xm,
                                         *this->auxiliary_vars, this->dx,
                                         ADD_MODE::REPLACE);

      // store residual in xm
      this->xm.waxpbypcz(-1, 1, -dt, this->xn, *this->x, this->dx);
      this->xm.exchange();

      if (si_monitor_convergence > 1) {
        res_norm = norm(xm);
      }

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

    if (si_verbosity_level == 1) {
      res_norm = norm(xm);
    }

    if (si_verbosity_level > 0) {
      std::stringstream msg;
      if (converged) {
        msg << "Newton solve converged in " << iter << " iters.\n";
      } else {
        msg << "!!! Newton solve failed to converge in " << iter << " iters.\n";
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