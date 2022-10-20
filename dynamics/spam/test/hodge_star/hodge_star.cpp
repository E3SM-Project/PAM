// clang-format off
#include "layer_common.h"
#include "hodge_star.h"
// clang-format on

struct fun_x {
  real YAKL_INLINE operator()(real x, real y) const {
    return sin(2 * M_PI * x);
  }
};

struct fun_y {
  real YAKL_INLINE operator()(real x, real y) const {
    return sin(2 * M_PI * y);
  }
};

struct fun_xy {
  real YAKL_INLINE operator()(real x, real y) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    return sx * sy;
  }
};

struct vecfun_x {
  vec<2> YAKL_INLINE operator()(real x, real y) const {
    real sx = sin(2 * M_PI * x);

    vec<2> vvec;
    vvec.u = sx;
    vvec.v = 0;
    return vvec;
  }
};

struct vecfun_y {
  vec<2> YAKL_INLINE operator()(real x, real y) const {
    real sy = sin(2 * M_PI * y);

    vec<2> vvec;
    vvec.u = 0;
    vvec.v = sy;
    return vvec;
  }
};

struct vecfun_xy {
  vec<2> YAKL_INLINE operator()(real x, real y) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);

    vec<2> vvec;
    vvec.u = sx * sy;
    vvec.v = sx * sy * sx * sy;
    return vvec;
  }
};

template <int diff_ord, class F> real compute_H2bar_error(int np, F ic_fun) {
  PeriodicUnitSquare square(np, 2 * np);

  auto tw2 = square.create_twisted_form<2>();
  square.dual_geometry.set_2form_values(ic_fun, tw2, 0);

  auto st0 = square.create_straight_form<0>();
  auto st0_expected = square.create_straight_form<0>();
  square.primal_geometry.set_0form_values(ic_fun, st0_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw2.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H2bar<1, diff_ord>(st0.data, tw2.data, square.primal_geometry,
                                     square.dual_geometry, pis, pjs, pks, i, j,
                                     k, 0);
        });
  }

  real errf = square.compute_Linf_error(st0_expected, st0);
  return errf;
}

void test_H2bar_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2bar 2 x", compute_H2bar_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H2bar 2 y", compute_H2bar_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H2bar 2 xy", compute_H2bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2bar 4 x", compute_H2bar_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H2bar 4 y", compute_H2bar_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H2bar 4 xy", compute_H2bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(4, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2bar 6 x", compute_H2bar_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H2bar 6 y", compute_H2bar_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H2bar 6 xy", compute_H2bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(4, atol);
  }
}

template <int diff_ord, class F> real compute_H2_error(int np, F ic_fun) {
  PeriodicUnitSquare square(np, 2 * np);

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
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H2<1, diff_ord>(tw0.data, st2.data, square.primal_geometry,
                                  square.dual_geometry, dis, djs, dks, i, j, k,
                                  0);
        });
  }

  real errf = square.compute_Linf_error(tw0_expected, tw0);
  return errf;
}

void test_H2_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2 2 x", compute_H2_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H2 2 y", compute_H2_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H2 2 xy", compute_H2_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2 4 x", compute_H2_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H2 4 y", compute_H2_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H2 4 xy", compute_H2_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(4, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2 6 x", compute_H2_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H2 6 y", compute_H2_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H2 6 xy", compute_H2_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(4, atol);
  }
}

template <int diff_ord, class F> real compute_H0_error(int np, F ic_fun) {
  PeriodicUnitSquare square(np, 2 * np);

  auto st0 = square.create_straight_form<0>();
  square.primal_geometry.set_0form_values(ic_fun, st0, 0);

  auto tw2 = square.create_twisted_form<2>();
  auto tw2_expected = square.create_twisted_form<2>();
  square.dual_geometry.set_2form_values(ic_fun, tw2_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st0.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H0<1, diff_ord>(tw2.data, st0.data, square.primal_geometry,
                                  square.dual_geometry, dis, djs, dks, i, j, k,
                                  0);
        });
  }

  real errf = square.compute_Linf_error(tw2_expected, tw2);
  return errf;
}

void test_H0_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H0 2 x", compute_H0_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H0 2 y", compute_H0_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H0 2 xy", compute_H0_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(diff_order, atol);
  }
}

template <int diff_ord, class F> real compute_H0bar_error(int np, F ic_fun) {
  PeriodicUnitSquare square(np, 2 * np);

  auto tw0 = square.create_twisted_form<0>();
  square.dual_geometry.set_0form_values(ic_fun, tw0, 0);

  auto st2 = square.create_straight_form<2>();
  auto st2_expected = square.create_straight_form<2>();
  square.primal_geometry.set_2form_values(ic_fun, st2_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw0.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H0bar<1, diff_ord>(st2.data, tw0.data, square.primal_geometry,
                                     square.dual_geometry, pis, pjs, pks, i, j,
                                     k, 0);
        });
  }

  real errf = square.compute_Linf_error(st2_expected, st2);
  return errf;
}

void test_H0bar_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H0bar 2 x", compute_H0bar_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H0bar 2 y", compute_H0bar_error<diff_order, fun_y>, fun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H0bar 2 xy", compute_H0bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(diff_order, atol);
  }
}

template <int diff_ord, class F> real compute_H1_error(int np, F ic_fun) {
  PeriodicUnitSquare square(np, 2 * np);

  auto st1 = square.create_straight_form<1>();
  square.primal_geometry.set_1form_values(ic_fun, st1, 0,
                                          LINE_INTEGRAL_TYPE::TANGENT);

  auto tw1 = square.create_twisted_form<1>();
  auto tw1_expected = square.create_twisted_form<1>();
  square.dual_geometry.set_1form_values(ic_fun, tw1_expected, 0,
                                        LINE_INTEGRAL_TYPE::NORMAL);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st1.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H1<1, diff_ord>(tw1.data, st1.data, square.primal_geometry,
                                  square.dual_geometry, dis, djs, dks, i, j, k,
                                  0);
        });
  }

  real errf = square.compute_Linf_error(tw1_expected, tw1);
  return errf;
}

void test_H1_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H1 2 x", compute_H1_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H1 2 y", compute_H1_error<diff_order, vecfun_y>, vecfun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H1 2 xy", compute_H1_error<diff_order, vecfun_xy>, vecfun_xy{});
    conv_xy.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "H1 4 x", compute_H1_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H1 4 y", compute_H1_error<diff_order, vecfun_y>, vecfun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H1 4 xy", compute_H1_error<diff_order, vecfun_xy>, vecfun_xy{});
    conv_xy.check_rate(2, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "H1 6 x", compute_H1_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_y = ConvergenceTest<nlevels>(
        "H1 6 y", compute_H1_error<diff_order, vecfun_y>, vecfun_y{});
    conv_y.check_rate(diff_order, atol);
    auto conv_xy = ConvergenceTest<nlevels>(
        "H1 6 xy", compute_H1_error<diff_order, vecfun_xy>, vecfun_xy{});
    conv_xy.check_rate(2, atol);
  }
}

int main() {
  yakl::init();

  test_H0_convergence();
  test_H0bar_convergence();

  test_H1_convergence();

  test_H2_convergence();
  test_H2bar_convergence();

  yakl::finalize();
}
