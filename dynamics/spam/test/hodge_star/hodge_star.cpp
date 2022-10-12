#include "layer_common.h"
#include "hodge_star.h"

real YAKL_INLINE fun_x(real x, real y) { return sin(2 * M_PI * x); }

real YAKL_INLINE fun_y(real x, real y) { return sin(2 * M_PI * y); }

real YAKL_INLINE fun_xy(real x, real y) {
  real sx = sin(2 * M_PI * x);
  real sy = sin(2 * M_PI * y);
  return sx * sy;
}

vec<2> YAKL_INLINE vecfun_x(real x, real y) {
  real sx = sin(2 * M_PI * x);

  vec<2> vvec;
  vvec.u = sx;
  vvec.v = 0;
  return vvec;
}

vec<2> YAKL_INLINE vecfun_y(real x, real y) {
  real sy = sin(2 * M_PI * y);

  vec<2> vvec;
  vvec.u = 0;
  vvec.v = sy;
  return vvec;
}

vec<2> YAKL_INLINE vecfun_xy(real x, real y) {
  real sx = sin(2 * M_PI * x);
  real sy = sin(2 * M_PI * y);

  vec<2> vvec;
  vvec.u = sx * sy;
  vvec.v = sx * sy * sx * sy;
  return vvec;
}

template <int diff_ord>
real compute_I_error(int np, real (*ic_fun)(real, real)) {
  PeriodicUnitSquare square(np, np);

  auto tw2 = square.create_twisted_form<2>();
  square.dual_geometry.set_2form_values(ic_fun, tw2, 0);

  auto st0 = square.create_straight_form<0>();
  auto st0_expected = square.create_straight_form<0>();
  square.primal_geometry.set_0form_values(ic_fun, st0_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    tw2.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_I<1, diff_ord>(st0.data, tw2.data, square.primal_geometry,
                                 square.dual_geometry, pis, pjs, pks, i, j, k,
                                 0);
        });
  }

  real errf = square.compute_Linf_error(st0_expected, st0);
  return errf;
}

void test_I_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x =
        ConvergenceTest<nlevels>("I 2 x", compute_I_error<diff_order>, fun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y =
        ConvergenceTest<nlevels>("I 2 y", compute_I_error<diff_order>, fun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy =
        ConvergenceTest<nlevels>("I 2 xy", compute_I_error<diff_order>, fun_xy);
    conv_xy.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x =
        ConvergenceTest<nlevels>("I 4 x", compute_I_error<diff_order>, fun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y =
        ConvergenceTest<nlevels>("I 4 y", compute_I_error<diff_order>, fun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy =
        ConvergenceTest<nlevels>("I 4 xy", compute_I_error<diff_order>, fun_xy);
    conv_xy.check_rate(4, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x =
        ConvergenceTest<nlevels>("I 6 x", compute_I_error<diff_order>, fun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y =
        ConvergenceTest<nlevels>("I 6 y", compute_I_error<diff_order>, fun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy =
        ConvergenceTest<nlevels>("I 6 xy", compute_I_error<diff_order>, fun_xy);
    conv_xy.check_rate(4, atol);
  }
}

template <int diff_ord>
real compute_J_error(int np, real (*ic_fun)(real, real)) {
  PeriodicUnitSquare square(np, np);

  auto st2 = square.create_straight_form<2>();
  square.primal_geometry.set_2form_values(ic_fun, st2, 0);

  auto tw0 = square.create_twisted_form<0>();
  auto tw0_expected = square.create_twisted_form<0>();
  square.dual_geometry.set_0form_values(ic_fun, tw0_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st2.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_J<1, diff_ord>(tw0.data, st2.data, square.primal_geometry,
                                 square.dual_geometry, dis, djs, dks, i, j, k,
                                 0);
        });
  }

  real errf = square.compute_Linf_error(tw0_expected, tw0);
  return errf;
}

void test_J_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x =
        ConvergenceTest<nlevels>("I 2 x", compute_J_error<diff_order>, fun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y =
        ConvergenceTest<nlevels>("I 2 y", compute_J_error<diff_order>, fun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy =
        ConvergenceTest<nlevels>("I 2 xy", compute_J_error<diff_order>, fun_xy);
    conv_xy.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x =
        ConvergenceTest<nlevels>("I 4 x", compute_J_error<diff_order>, fun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y =
        ConvergenceTest<nlevels>("I 4 y", compute_J_error<diff_order>, fun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy =
        ConvergenceTest<nlevels>("I 4 xy", compute_J_error<diff_order>, fun_xy);
    conv_xy.check_rate(4, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x =
        ConvergenceTest<nlevels>("I 6 x", compute_J_error<diff_order>, fun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y =
        ConvergenceTest<nlevels>("I 6 y", compute_J_error<diff_order>, fun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy =
        ConvergenceTest<nlevels>("I 6 xy", compute_J_error<diff_order>, fun_xy);
    conv_xy.check_rate(4, atol);
  }
}

template <int diff_ord>
real compute_H_error(int np, vec<2> (*ic_fun)(real, real)) {
  PeriodicUnitSquare square(np, np);

  auto st1 = square.create_straight_form<1>();
  square.primal_geometry.set_1form_values(ic_fun, st1, 0,
                                          LINE_INTEGRAL_TYPE::TANGENT);

  auto tw1 = square.create_twisted_form<1>();
  auto tw1_expected = square.create_twisted_form<1>();
  square.dual_geometry.set_1form_values(ic_fun, tw1_expected, 0,
                                        LINE_INTEGRAL_TYPE::NORMAL);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st1.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H<1, diff_ord>(tw1.data, st1.data, square.primal_geometry,
                                 square.dual_geometry, dis, djs, dks, i, j, k,
                                 0);
        });
  }

  real errf = square.compute_Linf_error(tw1_expected, tw1);
  return errf;
}

void test_H_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>("H 2 x", compute_H_error<diff_order>,
                                           vecfun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>("H 2 y", compute_H_error<diff_order>,
                                           vecfun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H 2 xy", compute_H_error<diff_order>, vecfun_xy);
    conv_xy.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>("H 4 x", compute_H_error<diff_order>,
                                           vecfun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>("H 4 y", compute_H_error<diff_order>,
                                           vecfun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H 4 xy", compute_H_error<diff_order>, vecfun_xy);
    conv_xy.check_rate(2, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>("H 6 x", compute_H_error<diff_order>,
                                           vecfun_x);
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>("H 6 y", compute_H_error<diff_order>,
                                           vecfun_y);
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H 6 xy", compute_H_error<diff_order>, vecfun_xy);
    conv_xy.check_rate(2, atol);
  }
}

int main() {
  yakl::init();

  test_I_convergence();
  test_J_convergence();
  test_H_convergence();

  yakl::finalize();
}
