// clang-format off
#include "pam_const.h"
#include "common.h"
#include <array>
#include <iostream>
// clang-format on

namespace pamc {
uint constexpr ndims = 2;
uint constexpr nprognostic = 0;
uint constexpr nconstant = 0;
uint constexpr nauxiliary = 0;
uint constexpr ndiagnostic = 0;
uint constexpr ntracers_dycore = 0;
} // namespace pamc

#include "params.h"

namespace pamc {
struct ModelParameters : public Parameters {
  // std::string initdataStr;
  std::string tracerdataStr[ntracers_dycore + GPU_PAD];
  bool dycore_tracerpos[ntracers_dycore + GPU_PAD];
  // bool acoustic_balance;
};
} // namespace pamc

#include "exchange.h"
#include "fields.h"
#include "geometry.h"
#include "topology.h"

namespace pamc {
Parallel parallel_stub(int nx, int ny) {
  Parallel par;
  par.nx_glob = par.nx = nx;
  par.ny_glob = par.ny = ny;
  par.nz = 1;

  par.nprocx = par.nprocy = par.nranks = 1;
  par.masterproc = par.myrank = par.px = par.py = 0;

  par.i_beg = 0;
  par.i_end = nx - 1;
  par.j_beg = 0;
  par.j_end = ny - 1;

  par.halox = maxhalosize;
  par.haloy = maxhalosize;
  par.x_neigh(0) = par.x_neigh(1) = par.y_neigh(0) = par.y_neigh(1) = 0;
  par.ll_neigh = par.ur_neigh = par.ul_neigh = par.lr_neigh = 0;

  par.nens = 1;

  return par;
}

struct PeriodicUnitSquare {
  Topology primal_topology;
  Topology dual_topology;

  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

  bool is_initialized = false;

  void initialize(int nx, int ny) {
    ModelParameters params;

    params.nx_glob = nx;
    params.ny_glob = ny;
    params.xlen = 1;
    params.xc = 0;
    params.ylen = 1;
    params.yc = 0;

    Parallel par = parallel_stub(nx, ny);

    primal_topology.initialize(par, true);
    dual_topology.initialize(par, false);

    primal_geometry.initialize(primal_topology, params);
    dual_geometry.initialize(dual_topology, params);

    is_initialized = true;
  }

  PeriodicUnitSquare(int nx, int ny) { initialize(nx, ny); }

  template <int deg> Field create_straight_form() {
    static Exchange exchange;
    exchange.initialize(primal_topology, deg, 1, 1);
    Field field;
    field.initialize(primal_topology, &exchange, "test field", deg, 1, 1);
    return field;
  }

  template <int deg> Field create_twisted_form() {
    static Exchange exchange;
    exchange.initialize(dual_topology, deg, 1, 1);
    Field field;
    field.initialize(dual_topology, &exchange, "test field", deg, 1, 1);
    return field;
  }

  real compute_Linf_error(const Field &f1, const Field &f2) {
    Field error;
    error.initialize(f1, "error");

    int is = error.topology.is;
    int js = error.topology.js;
    int ks = error.topology.ks;

    SArray<real, 1, 2> scale;
    real dx = primal_geometry.dx;
    real dy = primal_geometry.dy;

    if (f1.basedof == 0) {
      scale(0) = 1;
    }
    if (f1.basedof == 1) {
      if (f1.topology.primal) {
        scale(0) = dx;
        scale(1) = dy;
      } else {
        scale(0) = dy;
        scale(1) = dx;
      }
    }
    if (f1.basedof == 2) {
      scale(0) = dx * dy;
    }

    error.set(0);
    parallel_for(
        SimpleBounds<4>(error.total_dofs, error.topology.nl,
                        error.topology.n_cells_y, error.topology.n_cells_x),
        YAKL_LAMBDA(int l, int k, int j, int i) {
          error.data(l, ks + k, js + j, is + i, 0) =
              (abs(f1.data(l, ks + k, js + j, is + i, 0) -
                   f2.data(l, ks + k, js + j, is + i, 0))) /
              scale(l);
        });

    return yakl::intrinsics::maxval(error.data);
  }
};

template <int nlevels> struct ConvergenceTest {
  std::string name;
  std::array<real, nlevels> errs;
  std::array<real, nlevels - 1> rates;

  template <typename F1, typename F2>
  ConvergenceTest(const std::string &name, F1 error_fun, F2 ic_fun)
      : name(name) {
    for (int l = 0; l < nlevels; ++l) {
      int np = 8 * std::pow(2, l - 1);
      errs[l] = error_fun(np, ic_fun);
    }

    for (int l = 0; l < nlevels; ++l) {
      if (l < nlevels - 1) {
        rates[l] = std::log2(errs[l] / errs[l + 1]);
      }
    }
  }

  template <typename F1>
  ConvergenceTest(const std::string &name, F1 error_fun) : name(name) {
    for (int l = 0; l < nlevels; ++l) {
      int np = 8 * std::pow(2, l - 1);
      errs[l] = error_fun(np);
    }

    for (int l = 0; l < nlevels; ++l) {
      if (l < nlevels - 1) {
        rates[l] = std::log2(errs[l] / errs[l + 1]);
      }
    }
  }

  void print_errors_and_rates() const {
    std::cout << "Errors" << std::endl;
    for (int l = 0; l < nlevels; ++l) {
      std::cout << l << " " << errs[l] << std::endl;
    }
    std::cout << "Rates" << std::endl;
    for (int l = 0; l < nlevels - 1; ++l) {
      std::cout << l << " " << rates[l] << std::endl;
    }
  }

  void check_rate(real expected_rate, real atol) const {
    auto rate = rates[nlevels - 2];
    if (std::abs(rate - expected_rate) > atol || std::isnan(rate)) {
      std::cout << "Failed convergence test: " << name << std::endl;
      std::cout << "Achieved rate: " << rate << " "
                << "Expected rate: " << expected_rate << std::endl;
      print_errors_and_rates();
      exit(-1);
    }
  }
};
} // namespace pamc
