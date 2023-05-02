// clang-format off
unsigned constexpr ndims = 1;
#include "extruded_common.h"
#include "hodge_star.h"
// clang-format on

struct fun_x {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    return sin(2 * M_PI * x);
  }
};

struct fun_z {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    return z * z * z;
  }
};

struct fun_xz {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = z * z * z;
    return sx * sz;
  }
};

struct vecfun_x {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);

    VecXYZ vvec;
    vvec.u = sx;
    vvec.w = 0;
    return vvec;
  }
};

struct vecfun_z {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real sz = z * z * z;

    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = sz;
    return vvec;
  }
};

struct vecfun_xz {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = z * z * z;

    VecXYZ vvec;
    vvec.u = sz * sx;
    vvec.w = sx * sz * sz * sz;
    return vvec;
  }
};

template <int diff_ord, class F>
real compute_H00_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(ic_fun, st00, 0);

  auto tw11 = square.create_twisted_form<1, 1>();
  auto tw11_expected = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_n1form_values(ic_fun, tw11_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st00.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H00<1, diff_ord, vert_diff_ord>(
              tw11.data, st00.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw11_expected, tw11);
  return errf;
}

void test_H00_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.15;

  {
    const int diff_order = 2;
    auto conv_x =
        ConvergenceTest<nlevels>("H00 2 x", uniform_vertical,
                                 compute_H00_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("H00 2 z", uniform_vertical,
                                 compute_H00_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H00 2 xz", uniform_vertical, compute_H00_error<diff_order, fun_xz>,
        fun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_H00bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw00 = square.create_twisted_form<0, 0>();
  square.dual_geometry.set_00form_values(ic_fun, tw00, 0);

  auto st11 = square.create_straight_form<1, 1>();
  auto st11_expected = square.create_straight_form<1, 1>();
  square.primal_geometry.set_n1form_values(ic_fun, st11_expected, 0);

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
          compute_H00bar<1, diff_ord, vert_diff_ord>(
              st11.data, tw00.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st11_expected, st11);
  return errf;
}

void test_H00bar_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H00bar 2 x", uniform_vertical, compute_H00bar_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H00bar 2 z", uniform_vertical, compute_H00bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H00bar 2 xz", uniform_vertical,
        compute_H00bar_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_H10_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st10 = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(ic_fun, st10, 0);

  auto tw01 = square.create_twisted_form<0, 1>();
  auto tw01_expected = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_nm11form_values(ic_fun, tw01_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st10.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H10<1, diff_ord>(tw01.data, st10.data, square.primal_geometry,
                                   square.dual_geometry, dis, djs, dks, i, j, k,
                                   0);
        });
  }

  real errf = square.compute_Linf_error(tw01_expected, tw01);
  return errf;
}

void test_H10_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "H10 2 x", uniform_vertical, compute_H10_error<diff_order, vecfun_x>,
        vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H10 2 xz", uniform_vertical, compute_H10_error<diff_order, vecfun_xz>,
        vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;

    auto conv_x = ConvergenceTest<nlevels>(
        "H10 4 x", uniform_vertical, compute_H10_error<diff_order, vecfun_x>,
        vecfun_x{});
    conv_x.check_rate(diff_order, atol);

    auto conv_xz = ConvergenceTest<nlevels>(
        "H10 4 xz", uniform_vertical, compute_H10_error<diff_order, vecfun_xz>,
        vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;

    auto conv_x = ConvergenceTest<nlevels>(
        "H10 6 x", uniform_vertical, compute_H10_error<diff_order, vecfun_x>,
        vecfun_x{});
    conv_x.check_rate(diff_order, atol);

    auto conv_xz = ConvergenceTest<nlevels>(
        "H10 6 xz", uniform_vertical, compute_H10_error<diff_order, vecfun_xz>,
        vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int vdiff_ord, class F>
real compute_H01_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st01 = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(ic_fun, st01, 0);

  auto tw10 = square.create_twisted_form<1, 0>();
  auto tw10_expected = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_n0form_values(ic_fun, tw10_expected, 0);

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
          compute_H01<1, vdiff_ord>(
              tw10.data, st01.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw10_expected, tw10, false);
  return errf;
}

void test_H01_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "H01 2 z", uniform_vertical, compute_H01_error<vert_diff_ord, vecfun_z>,
        vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "H01 2 xz", uniform_vertical,
        compute_H01_error<vert_diff_ord, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(vert_diff_ord, atol);
  }
}

template <int diff_ord, class F>
real compute_Hnm11bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw01 = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_nm11form_values(ic_fun, tw01, 0);

  auto st10 = square.create_straight_form<1, 0>();
  auto st10_expected = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(ic_fun, st10_expected, 0);

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
          compute_Hnm11bar<1, diff_ord>(
              st10.data, tw01.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st10_expected, st10);
  return errf;
}

void test_Hnm11bar_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "Hnm11bar 2 x", uniform_vertical,
        compute_Hnm11bar_error<diff_order, vecfun_x>, vecfun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hnm11bar 2 xz", uniform_vertical,
        compute_Hnm11bar_error<diff_order, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int vdiff_ord, class F>
real compute_Hn0bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw10 = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_n0form_values(ic_fun, tw10, 0);

  auto st01 = square.create_straight_form<0, 1>();
  auto st01_expected = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(ic_fun, st01_expected, 0);

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
          compute_Hn0bar<1, vdiff_ord>(
              st01.data, tw10.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st01_expected, st01);
  return errf;
}

void test_Hn0bar_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn0bar 2 z", uniform_vertical,
        compute_Hn0bar_error<vert_diff_ord, vecfun_z>, vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hn0bar 2 xz", uniform_vertical,
        compute_Hn0bar_error<vert_diff_ord, vecfun_xz>, vecfun_xz{});
    conv_xz.check_rate(vert_diff_ord, atol);
  }
}

template <int diff_ord, class F>
real compute_Hn1bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw11 = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_n1form_values(ic_fun, tw11, 0);

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
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              st00.data, tw11.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st00_expected, st00);
  return errf;
}

void test_Hn1bar_convergence(bool uniform_vertical) {
  const int nlevels = 6;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>(
        "Hn1bar 2 x", uniform_vertical, compute_Hn1bar_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn1bar 2 z", uniform_vertical, compute_Hn1bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hn1bar 2 xz", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>(
        "Hn1bar 4 x", uniform_vertical, compute_Hn1bar_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn1bar 4 z", uniform_vertical, compute_Hn1bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hn1bar 4 xz", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>(
        "Hn1bar 6 x", uniform_vertical, compute_Hn1bar_error<diff_order, fun_x>,
        fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn1bar 6 z", uniform_vertical, compute_Hn1bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hn1bar 6 xz", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xz>, fun_xz{});
    conv_xz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_Hn1_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st11 = square.create_straight_form<1, 1>();
  square.primal_geometry.set_n1form_values(ic_fun, st11, 0);

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
          compute_Hn1<1, diff_ord, vert_diff_ord>(
              tw00.data, st11.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw00_expected, tw00, false);
  return errf;
}

void test_Hn1_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_x =
        ConvergenceTest<nlevels>("Hn1 2 x", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("Hn1 2 z", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hn1 2 xz", uniform_vertical, compute_Hn1_error<diff_order, fun_xz>,
        fun_xz{});
    conv_xz.check_rate(2, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x =
        ConvergenceTest<nlevels>("Hn1 4 x", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("Hn1 4 z", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hn1 4 xz", uniform_vertical, compute_Hn1_error<diff_order, fun_xz>,
        fun_xz{});
    conv_xz.check_rate(2, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x =
        ConvergenceTest<nlevels>("Hn1 6 x", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_x>, fun_x{});
    conv_x.check_rate(diff_order, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("Hn1 6 z", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xz = ConvergenceTest<nlevels>(
        "Hn1 6 xz", uniform_vertical, compute_Hn1_error<diff_order, fun_xz>,
        fun_xz{});
    conv_xz.check_rate(2, atol);
  }
}

int main() {
  yakl::init();

  for (bool uniform_vertical : {true, false}) {
    test_H00_convergence(uniform_vertical);
    test_H00bar_convergence(uniform_vertical);

    test_H10_convergence(uniform_vertical);
    test_H01_convergence(uniform_vertical);
    test_Hnm11bar_convergence(uniform_vertical);
    test_Hn0bar_convergence(uniform_vertical);

    test_Hn1_convergence(uniform_vertical);
    test_Hn1bar_convergence(uniform_vertical);
  }

  yakl::finalize();
}
