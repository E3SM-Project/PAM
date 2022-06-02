#pragma once

#include "common.h"
#include "exchange.h"
#include "field_sets.h"
#include "model.h"
#include "topology.h"

class SITimeIntegrator {

public:
  FieldSet<nprognostic> F;
  FieldSet<nprognostic> F1;
  FieldSet<nprognostic> F2;
  FieldSet<nprognostic> *x;
  FieldSet<nprognostic> dx;
  FieldSet<nprognostic> xn;
  FieldSet<nprognostic> x1;
  FieldSet<nprognostic> x2;
  FieldSet<nprognostic> xm;
  FieldSet<nprognostic> res;
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

SITimeIntegrator::SITimeIntegrator() {
  this->is_initialized = false;
  std::cout << "CREATED SI\n";
}

void SITimeIntegrator::initialize(ModelParameters &params, Tendencies &tend,
                                  LinearSystem &linsys,
                                  FieldSet<nprognostic> &xvars,
                                  FieldSet<nconstant> &consts,
                                  FieldSet<nauxiliary> &auxiliarys,
                                  ExchangeSet<nprognostic> &prog_exch) {

  this->F.initialize(xvars, "F");
  this->F1.initialize(xvars, "F1");
  this->F2.initialize(xvars, "F2");
  this->dx.initialize(xvars, "dx");
  this->xn.initialize(xvars, "xn");
  this->xm.initialize(xvars, "xm");
  this->x1.initialize(xvars, "x1");
  this->x2.initialize(xvars, "x2");
  this->res.initialize(xvars, "res");
  this->x = &xvars;
  this->tendencies = &tend;
  this->linear_system = &linsys;
  this->const_vars = &consts;
  this->auxiliary_vars = &auxiliarys;
  this->x_exchange = &prog_exch;

  // this->tendencies->initialize_linsolver(params, tend, prog_exch, xvars,
  // consts, auxiliarys);

  this->is_initialized = true;
}

// real slownorm(FieldSet<nprognostic> &x, ExchangeSet<nprognostic> &x_exch)
//{
//   x_exch.exchange_variable_set(x);
//   real accum = 0;
//   for (auto f : x.fields_arr)
//   {
//     accum = std::max(yakl::intrinsics::maxval(yakl::intrinsics::abs(f.data)),
//     accum);
//   }
//
//   return accum;
// }
real slownorm(FieldSet<nprognostic> &x, ExchangeSet<nprognostic> &x_exch) {
  real accum = 0;

  for (auto f : x.fields_arr) {
    auto topo = f.topology;
    real4d trim("trim", f.total_dofs, f._nz, topo->n_cells_y, topo->n_cells_x);
    int is = topo->is;
    int js = topo->js;
    int ks = topo->ks;
    parallel_for(
        "Compute trimmed density",
        // SimpleBounds<4>(f.total_dofs, topo->nl, topo->n_cells_y,
        SimpleBounds<4>(1, topo->nl, topo->n_cells_y, topo->n_cells_x),
        YAKL_LAMBDA(int l, int k, int j, int i) {
          trim(l, k, j, i) = f.data(l, k + ks, j + js, i + is, 0);
        });
    accum =
        std::max(yakl::intrinsics::maxval(yakl::intrinsics::abs(trim)), accum);
  }

  return accum;
}

void SITimeIntegrator::stepForward(real dt) {

  real gamma1 = 0.5_fp * (1._fp - 1._fp / std::sqrt(3._fp));
  real gamma2 = 0.5_fp * (1._fp + 1._fp / std::sqrt(3._fp));

  this->tendencies->compute_rhs(dt, *this->const_vars, *this->x,
                                *this->auxiliary_vars, this->F);
  this->xn.copy(*this->x);
  this->res.waxpbypcz(-1, 1, -dt, this->xn, *this->x, this->F);

  int iter = 0;

  real res_norm = slownorm(res, *this->x_exchange);
  real initial_norm = res_norm;
  std::cout << "Initial res norm: " << res_norm << std::endl;
  while (true) {
    if (res_norm / initial_norm < 1e-10 || iter > 200) {
      break;
    }

    this->x_exchange->exchange_variable_set(this->xn);
    this->x_exchange->exchange_variable_set(this->res);
    this->linear_system->solve(dt, this->res, *this->const_vars,
                               *this->auxiliary_vars, this->dx);

    this->xn.waxpy(1, this->dx, this->xn);

    this->xm.waxpby(0.5_fp, 0.5_fp, *this->x, this->xn);
    this->x_exchange->exchange_variable_set(this->xm);
    this->x1.waxpby(1 - gamma1, gamma1, *this->x, this->xn);
    this->x2.waxpby(1 - gamma2, gamma2, *this->x, this->xn);

    this->x_exchange->exchange_variable_set(this->x1);
    this->tendencies->compute_functional_derivatives(
        dt, *this->const_vars, this->x1, *this->auxiliary_vars);
    this->tendencies->apply_symplectic(dt, *this->const_vars, this->xm,
                                       *this->auxiliary_vars, this->F1);

    this->x_exchange->exchange_variable_set(this->x2);
    this->tendencies->compute_functional_derivatives(
        dt, *this->const_vars, this->x2, *this->auxiliary_vars);
    this->tendencies->apply_symplectic(dt, *this->const_vars, this->xm,
                                       *this->auxiliary_vars, this->F2);

    this->F.waxpby(0.5_fp, 0.5_fp, this->F1, this->F2);
    this->res.waxpbypcz(-1, 1, -dt, this->xn, *this->x, this->F);

    res_norm = slownorm(res, *this->x_exchange);

    iter++;

    std::cout << "Iter: " << iter << " " << res_norm << std::endl;
  }

  this->x->copy(this->xn);
  this->x_exchange->exchange_variable_set(*this->x);
}