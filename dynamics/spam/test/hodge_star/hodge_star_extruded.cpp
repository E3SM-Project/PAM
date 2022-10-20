// clang-format off
#include "extruded_common.h"
#include "hodge_star.h"
// clang-format on

struct fun_x {
  real YAKL_INLINE operator()(real x, real z) const {
    return sin(2 * M_PI * x);
  }
};

struct fun_z {
  real YAKL_INLINE operator()(real x, real z) const {
    return sin(2 * M_PI * z);
  }
};

struct fun_xz {
  real YAKL_INLINE operator()(real x, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = sin(2 * M_PI * z);
    return sx * sz;
  }
};

struct vecfun_x {
  vecext<2> YAKL_INLINE operator()(real x, real z) const {
    real sx = sin(2 * M_PI * x);

    vecext<2> vvec;
    vvec.u = sx;
    vvec.w = 0;
    return vvec;
  }
};

struct vecfun_z {
  vecext<2> YAKL_INLINE operator()(real x, real z) const {
    real sz = sin(2 * M_PI * z);

    vecext<2> vvec;
    vvec.u = 0;
    vvec.w = sz;
    return vvec;
  }
};

struct vecfun_xz {
  vecext<2> YAKL_INLINE operator()(real x, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = sin(2 * M_PI * z);

    vecext<2> vvec;
    vvec.u = sx * sz;
    vvec.w = sx * sz * sz * sz;
    return vvec;
  }
};

template <int diff_ord, class F> real compute_Iext_error(int np, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto tw11 = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_11form_values(ic_fun, tw11, 0);

  auto st00 = square.create_straight_form<0, 0>();
  auto st00_expected = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(ic_fun, st00_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw11.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Iext<1, diff_ord, vert_diff_ord>(
              st00.data, tw11.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st00_expected, st00);
  return errf;
}

void test_Iext_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "Iext 2 x", compute_Iext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Iext 2 z", compute_Iext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Iext 2 xz", compute_Iext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "Iext 4 x", compute_Iext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Iext 4 z", compute_Iext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Iext 4 xz", compute_Iext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "Iext 6 x", compute_Iext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Iext 6 z", compute_Iext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Iext 6 xz", compute_Iext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int diff_ord, class F> real compute_Jext_error(int np, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto st11 = square.create_straight_form<1, 1>();
  square.primal_geometry.set_11form_values(ic_fun, st11, 0);

  auto tw00 = square.create_twisted_form<0, 0>();
  auto tw00_expected = square.create_twisted_form<0, 0>();
  square.dual_geometry.set_00form_values(ic_fun, tw00_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st11.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Jext<1, diff_ord, vert_diff_ord>(
              tw00.data, st11.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw00_expected, tw00, false);
  return errf;
}

void test_Jext_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "Jext 2 x", compute_Jext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Jext 2 z", compute_Jext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Jext 2 xz", compute_Jext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(2, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "Jext 4 x", compute_Jext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Jext 4 z", compute_Jext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Jext 4 xz", compute_Jext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(2, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "Jext 6 x", compute_Jext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Jext 6 z", compute_Jext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Jext 6 xz", compute_Jext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(2, atol);
  }
}

template <int diff_ord, class F> real compute_H1ext_error(int np, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto st10 = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(ic_fun, st10, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  auto tw01 = square.create_twisted_form<0, 1>();
  auto tw01_expected = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_01form_values(ic_fun, tw01_expected, 0,
                                         LINE_INTEGRAL_TYPE::NORMAL);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st10.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H1ext<1, diff_ord>(
              tw01.data, st10.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw01_expected, tw01);
  return errf;
}

void test_H1ext_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H1ext 2 x", compute_H1ext_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H1ext 2 xz", compute_H1ext_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;

    auto conv_x = ConvergenceTest<nlevels>(
        "H1ext 4 x", compute_H1ext_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);

    auto conv_xz = ConvergenceTest<nlevels>(
        "H1ext 4 xz", compute_H1ext_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;

    auto conv_x = ConvergenceTest<nlevels>(
        "H1ext 6 x", compute_H1ext_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);

    auto conv_xz = ConvergenceTest<nlevels>(
        "H1ext 6 xz", compute_H1ext_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int vdiff_ord, class F> real compute_H1vert_error(int np, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto st01 = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(ic_fun, st01, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  auto tw10 = square.create_twisted_form<1, 0>();
  auto tw10_expected = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_10form_values(ic_fun, tw10_expected, 0,
                                         LINE_INTEGRAL_TYPE::NORMAL);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st01.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H1vert<1, vdiff_ord>(tw10.data, st01.data, square.primal_geometry,
                                   square.dual_geometry, dis, djs, dks, i, j,
                                   k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw10_expected, tw10, false);
  return errf;
}

void test_H1vert_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "H1vert 2 z", compute_H1vert_error<vert_diff_ord, vecfun_z>, vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H1vert 2 xz", compute_H1vert_error<vert_diff_ord, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(vert_diff_ord, atol);
  }
}

int main() {
  yakl::init();
  test_Iext_convergence();
  test_Jext_convergence();
  test_H1ext_convergence();
  test_H1vert_convergence();
  yakl::finalize();
}
