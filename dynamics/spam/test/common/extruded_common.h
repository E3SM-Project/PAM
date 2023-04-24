// clang-format off
#include "pam_const.h"
#include "common.h"
#include <array>
#include <iostream>
// clang-format on

uint constexpr nprognostic = 0;
uint constexpr nconstant = 0;
uint constexpr nauxiliary = 0;
uint constexpr ndiagnostic = 0;
uint constexpr ntracers_dycore = 0;

#include "params.h"

struct ModelParameters : public Parameters {
  // std::string initdataStr;
  std::string tracerdataStr[ntracers_dycore];
  bool dycore_tracerpos[ntracers_dycore];
  // bool acoustic_balance;
  bool uniform_vertical;
  real2d zint;
};

#include "exchange.h"
#include "fields.h"
#include "geometry.h"
#include "params.h"
#include "topology.h"

Parallel parallel_stub(int nx, int ny, int nz) {
  Parallel par;
  par.nx_glob = par.nx = nx;
  par.ny_glob = par.ny = ny;

  par.nz = nz;

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

struct ExtrudedUnitSquare {
  Topology primal_topology;
  Topology dual_topology;

  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

  ExtrudedUnitSquare(int nx, int ny, int nz, bool uniform_vertical) {
    ModelParameters params;

    params.nx_glob = nx;
    params.ny_glob = ny;
    params.xlen = 1;
    params.ylen = 1;
    params.xc = 0;
    params.yc = 0;
    params.nz_dual = nz;
    params.uniform_vertical = uniform_vertical;

    params.zint = real2d("zint", nz + 1, 1);
    const real dz = 1.0 / (nz - 1);
    parallel_for(
        "compute zint", nz + 1, YAKL_LAMBDA(int k) {
          if (k == 0) {
            params.zint(k, 0) = 0;
          } else if (k == nz) {
            params.zint(k, 0) = 1;
          } else {
            params.zint(k, 0) = k * dz - dz / 2;
          }
          if (!uniform_vertical) {
            real z = params.zint(k, 0);
            params.zint(k, 0) = (-std::tanh(-0.5) + std::tanh((z - 0.5_fp))) /
                                (2 * std::tanh(0.5_fp));
          }
        });

    Parallel par = parallel_stub(nx, ny, nz);

    primal_topology.initialize(par, true);
    dual_topology.initialize(par, false);

    primal_geometry.initialize(primal_topology, params);
    dual_geometry.initialize(dual_topology, params);
  }

  ExtrudedUnitSquare(int nx, int nz, bool uniform_vertical)
      : ExtrudedUnitSquare(nx, 1, nz, uniform_vertical) {}

  template <int deg, int edeg> Field create_straight_form() {
    static Exchange exchange;
    exchange.initialize(primal_topology, deg, edeg, 1);
    Field field;
    field.initialize(primal_topology, &exchange, "test field", deg, edeg, 1);
    return field;
  }

  template <int deg, int edeg> Field create_twisted_form() {
    static Exchange exchange;
    exchange.initialize(dual_topology, deg, edeg, 1);
    Field field;
    field.initialize(dual_topology, &exchange, "test field", deg, edeg, 1);
    return field;
  }

  real compute_Linf_error(const Field &f1, const Field &f2,
                          bool compute_boundary_error = true) {
    Field error;
    error.initialize(f1, "error");

    int is = error.topology.is;
    int js = error.topology.js;
    int ks = error.topology.ks;

    error.set(0);

    int nz, koff;
    if (compute_boundary_error) {
      nz = error._nz;
      koff = 0;
    } else {
      nz = error._nz - 2;
      koff = 1;
    }

    parallel_for(
        SimpleBounds<4>(error.total_dofs, nz, error.topology.n_cells_y,
                        error.topology.n_cells_x),
        YAKL_LAMBDA(int l, int k, int j, int i) {
          real scale = 1;
          real dx = primal_geometry.dx;
          real dy = primal_geometry.dy;

          real dz = f1.topology.primal ? primal_geometry.dz(ks + ks + koff, 0)
                                       : dual_geometry.dz(ks + ks + koff, 0);

          if (f1.basedof >= 1) {
            scale *= dx;
          }
          if (f1.basedof >= 2) {
            scale *= dy;
          }
          if (f1.extdof == 1) {
            scale *= dz;
          }
          error.data(l, ks + k + koff, js + j, is + i, 0) =
              (abs(f1.data(l, ks + k + koff, js + j, is + i, 0) -
                   f2.data(l, ks + k + koff, js + j, is + i, 0)));

          error.data(l, ks + k + koff, js + j, is + i, 0) /= scale;
        });

    return yakl::intrinsics::maxval(error.data);
  }
};

template <int nlevels> struct ConvergenceTest {
  std::string name;
  std::array<real, nlevels> errs;
  std::array<real, nlevels - 1> rates;

  template <typename F1, typename F2>
  ConvergenceTest(const std::string &name, bool uniform_vertical, F1 error_fun,
                  F2 ic_fun)
      : name(name) {
    for (int l = 0; l < nlevels; ++l) {
      int np = 8 * std::pow(2, l - 1);
      errs[l] = error_fun(np, uniform_vertical, ic_fun);
    }

    for (int l = 0; l < nlevels; ++l) {
      if (l < nlevels - 1) {
        rates[l] = std::log2(errs[l] / errs[l + 1]);
      }
    }
  }

  template <typename F1>
  ConvergenceTest(const std::string &name, bool uniform_vertical, F1 error_fun)
      : name(name) {
    for (int l = 0; l < nlevels; ++l) {
      int np = 8 * std::pow(2, l - 1);
      errs[l] = error_fun(np, uniform_vertical);
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
