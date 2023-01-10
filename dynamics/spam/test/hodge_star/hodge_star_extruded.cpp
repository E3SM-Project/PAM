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
  real YAKL_INLINE operator()(real x, real z) const { return z * z * z; }
};

struct fun_xz {
  real YAKL_INLINE operator()(real x, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = z * z * z;
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
    real sz = z * z * z;

    vecext<2> vvec;
    vvec.u = 0;
    vvec.w = sz;
    return vvec;
  }
};

struct vecfun_xz {
  vecext<2> YAKL_INLINE operator()(real x, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = z * z * z;

    vecext<2> vvec;
    vvec.u = sz * sx;
    vvec.w = sx * sz * sz * sz;
    return vvec;
  }
};

template <int diff_ord, class F>
real compute_H0_ext_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(ic_fun, st00, 0);

  auto tw11 = square.create_twisted_form<1, 1>();
  auto tw11_expected = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_11form_values(ic_fun, tw11_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st00.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H0_ext<1, diff_ord, vert_diff_ord>(
              tw11.data, st00.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw11_expected, tw11);
  return errf;
}

void test_H0_ext_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.15;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H0_ext 2 x", uniform_vertical, compute_H0_ext_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H0_ext 2 z", uniform_vertical, compute_H0_ext_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H0_ext 2 xz", uniform_vertical,
        compute_H0_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_H0bar_ext_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw00 = square.create_twisted_form<0, 0>();
  square.dual_geometry.set_00form_values(ic_fun, tw00, 0);

  auto st11 = square.create_straight_form<1, 1>();
  auto st11_expected = square.create_straight_form<1, 1>();
  square.primal_geometry.set_11form_values(ic_fun, st11_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw00.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H0bar_ext<1, diff_ord, vert_diff_ord>(
              st11.data, tw00.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st11_expected, st11);
  return errf;
}

void test_H0bar_ext_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H0bar_ext 2 x", uniform_vertical,
        compute_H0bar_ext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H0bar_ext 2 z", uniform_vertical,
        compute_H0bar_ext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H0bar_ext 2 xz", uniform_vertical,
        compute_H0bar_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_H1_ext_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

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
          compute_H1_ext<1, diff_ord>(
              tw01.data, st10.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw01_expected, tw01);
  return errf;
}

void test_H1_ext_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H1_ext 2 x", uniform_vertical,
        compute_H1_ext_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H1_ext 2 xz", uniform_vertical,
        compute_H1_ext_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;

    auto conv_x = ConvergenceTest<nlevels>(
        "H1_ext 4 x", uniform_vertical,
        compute_H1_ext_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);

    auto conv_xz = ConvergenceTest<nlevels>(
        "H1_ext 4 xz", uniform_vertical,
        compute_H1_ext_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;

    auto conv_x = ConvergenceTest<nlevels>(
        "H1_ext 6 x", uniform_vertical,
        compute_H1_ext_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);

    auto conv_xz = ConvergenceTest<nlevels>(
        "H1_ext 6 xz", uniform_vertical,
        compute_H1_ext_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int vdiff_ord, class F>
real compute_H1_vert_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

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
          compute_H1_vert<1, vdiff_ord>(
              tw10.data, st01.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw10_expected, tw10, false);
  return errf;
}

void test_H1_vert_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "H1_vert 2 z", uniform_vertical,
        compute_H1_vert_error<vert_diff_ord, vecfun_z>, vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H1_vert 2 xz", uniform_vertical,
        compute_H1_vert_error<vert_diff_ord, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(vert_diff_ord, atol);
  }
}

template <int diff_ord, class F>
real compute_H1bar_ext_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw01 = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_01form_values(ic_fun, tw01, 0,
                                         LINE_INTEGRAL_TYPE::NORMAL);

  auto st10 = square.create_straight_form<1, 0>();
  auto st10_expected = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(ic_fun, st10_expected, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw01.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          st10_expected.data(0, k + pks, j + pjs, i + pis, 0) *= -1;
          compute_H1bar_ext<1, diff_ord>(
              st10.data, tw01.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st10_expected, st10);
  return errf;
}

void test_H1bar_ext_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H1bar_ext 2 x", uniform_vertical,
        compute_H1bar_ext_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H1bar_ext 2 xz", uniform_vertical,
        compute_H1bar_ext_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int vdiff_ord, class F>
real compute_H1bar_vert_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw10 = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_10form_values(ic_fun, tw10, 0,
                                         LINE_INTEGRAL_TYPE::NORMAL);

  auto st01 = square.create_straight_form<0, 1>();
  auto st01_expected = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(ic_fun, st01_expected, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw10.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          st01_expected.data(0, k + pks, j + pjs, i + pis, 0) *= -1;
          compute_H1bar_vert<1, vdiff_ord>(
              st01.data, tw10.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st01_expected, st01);
  return errf;
}

void test_H1bar_vert_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "H1bar_vert 2 z", uniform_vertical,
        compute_H1bar_vert_error<vert_diff_ord, vecfun_z>, vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H1bar_vert 2 xz", uniform_vertical,
        compute_H1bar_vert_error<vert_diff_ord, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(vert_diff_ord, atol);
  }
}

template <int diff_ord, class F>
real compute_H2bar_ext_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

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
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              st00.data, tw11.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st00_expected, st00);
  return errf;
}

void test_H2bar_ext_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2bar_ext 2 x", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H2bar_ext 2 z", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H2bar_ext 2 xz", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2bar_ext 4 x", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H2bar_ext 4 z", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H2bar_ext 4 xz", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2bar_ext 6 x", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H2bar_ext 6 z", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H2bar_ext 6 xz", uniform_vertical,
        compute_H2bar_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_H2_ext_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

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
          compute_H2_ext<1, diff_ord, vert_diff_ord>(
              tw00.data, st11.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw00_expected, tw00, false);
  return errf;
}

void test_H2_ext_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2_ext 2 x", uniform_vertical, compute_H2_ext_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H2_ext 2 z", uniform_vertical, compute_H2_ext_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H2_ext 2 xz", uniform_vertical,
        compute_H2_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(2, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2_ext 4 x", uniform_vertical, compute_H2_ext_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H2_ext 4 z", uniform_vertical, compute_H2_ext_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H2_ext 4 xz", uniform_vertical,
        compute_H2_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(2, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "H2_ext 6 x", uniform_vertical, compute_H2_ext_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H2_ext 6 z", uniform_vertical, compute_H2_ext_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H2_ext 6 xz", uniform_vertical,
        compute_H2_ext_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(2, atol);
  }
}

int main() {
  yakl::init();

  for (bool uniform_vertical : {true, false}) {
    test_H0_ext_convergence(uniform_vertical);
    test_H0bar_ext_convergence(uniform_vertical);

    test_H1_ext_convergence(uniform_vertical);
    test_H1_vert_convergence(uniform_vertical);
    test_H1bar_ext_convergence(uniform_vertical);
    test_H1bar_vert_convergence(uniform_vertical);

    test_H2_ext_convergence(uniform_vertical);
    test_H2bar_ext_convergence(uniform_vertical);
  }

  yakl::finalize();
}
