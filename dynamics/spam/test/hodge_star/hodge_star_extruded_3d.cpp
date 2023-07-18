// clang-format off
unsigned constexpr ndims = 2;
#include "extruded_common.h"
#include "hodge_star.h"
// clang-format on

using namespace pamc;

struct fun_xy {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    return sx * sy;
  }
};

struct fun_z {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    return z * z * z;
  }
};

struct fun_xyz {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real fz = z * z * z;
    return sx * sy * fz;
  }
};

struct vecfun_xy {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);

    VecXYZ vvec;
    vvec.u = sx;
    vvec.v = sy;
    vvec.w = 0;
    return vvec;
  }
};

struct vecfun_z {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real fz = z * z * z;

    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = fz;
    return vvec;
  }
};

struct vecfun_xyz {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real fz = z * z * z;

    VecXYZ vvec;
    vvec.u = sx * sy * fz;
    vvec.w = sx * sx * sx * sy * sy * fz;
    return vvec;
  }
};

template <int diff_ord, class F>
real compute_H00_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(ic_fun, st00, 0);

  auto tw21 = square.create_twisted_form<2, 1>();
  auto tw21_expected = square.create_twisted_form<2, 1>();
  square.dual_geometry.set_n1form_values(ic_fun, tw21_expected, 0);

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
              tw21.data, st00.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw21_expected, tw21);
  return errf;
}

void test_H00_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_xy = ConvergenceTest<nlevels>(
        "H00 2 xy", uniform_vertical, compute_H00_error<diff_order, fun_xy>,
        fun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("H00 2 z", uniform_vertical,
                                 compute_H00_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "H00 2 xyz", uniform_vertical, compute_H00_error<diff_order, fun_xyz>,
        fun_xyz{});
    conv_xyz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_H00bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto tw00 = square.create_twisted_form<0, 0>();
  square.dual_geometry.set_00form_values(ic_fun, tw00, 0);

  auto st21 = square.create_straight_form<2, 1>();
  auto st21_expected = square.create_straight_form<2, 1>();
  square.primal_geometry.set_n1form_values(ic_fun, st21_expected, 0);

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
              st21.data, tw00.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st21_expected, st21);
  return errf;
}

void test_H00bar_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_xy = ConvergenceTest<nlevels>(
        "H00bar 2 x", uniform_vertical,
        compute_H00bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "H00bar 2 z", uniform_vertical, compute_H00bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "H00bar 2 xyz", uniform_vertical,
        compute_H00bar_error<diff_order, fun_xyz>, fun_xyz{});
    conv_xyz.check_rate(1, atol);
  }
}

template <int diff_ord, class F>
real compute_H10_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st10 = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(ic_fun, st10, 0);

  auto tw11 = square.create_twisted_form<1, 1>();
  auto tw11_expected = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_nm11form_values(ic_fun, tw11_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st10.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H10<1, diff_ord>(tw11.data, st10.data, square.primal_geometry,
                                   square.dual_geometry, dis, djs, dks, i, j, k,
                                   0);
        });
  }

  real errf = square.compute_Linf_error(tw11_expected, tw11);
  return errf;
}

void test_H10_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_xy = ConvergenceTest<nlevels>(
        "H10 2 xy", uniform_vertical, compute_H10_error<diff_order, vecfun_xy>,
        vecfun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "H10 2 xyz", uniform_vertical,
        compute_H10_error<diff_order, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;

    auto conv_xy = ConvergenceTest<nlevels>(
        "H10 4 xy", uniform_vertical, compute_H10_error<diff_order, vecfun_xy>,
        vecfun_xy{});
    conv_xy.check_rate(diff_order, atol);

    auto conv_xyz = ConvergenceTest<nlevels>(
        "H10 4 xyz", uniform_vertical,
        compute_H10_error<diff_order, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;

    auto conv_xy = ConvergenceTest<nlevels>(
        "H10 6 xy", uniform_vertical, compute_H10_error<diff_order, vecfun_xy>,
        vecfun_xy{});
    conv_xy.check_rate(diff_order, atol);

    auto conv_xyz = ConvergenceTest<nlevels>(
        "H10 6 xyz", uniform_vertical,
        compute_H10_error<diff_order, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(1, atol);
  }
}

template <int vdiff_ord, class F>
real compute_H01_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st01 = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(ic_fun, st01, 0);

  auto tw20 = square.create_twisted_form<2, 0>();
  auto tw20_expected = square.create_twisted_form<2, 0>();
  square.dual_geometry.set_n0form_values(ic_fun, tw20_expected, 0);

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
              tw20.data, st01.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw20_expected, tw20, false);
  return errf;
}

void test_H01_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "H01 2 z", uniform_vertical, compute_H01_error<vert_diff_ord, vecfun_z>,
        vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "H01 2 xyz", uniform_vertical,
        compute_H01_error<vert_diff_ord, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(vert_diff_ord, atol);
  }
}

template <int diff_ord, class F>
real compute_Hnm11_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st11 = square.create_straight_form<1, 1>();
  square.primal_geometry.set_nm11form_values(ic_fun, st11, 0);

  auto tw10 = square.create_twisted_form<1, 0>();
  auto tw10_expected = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_10form_values(ic_fun, tw10_expected, 0);

  int dis = square.primal_topology.is;
  int djs = square.primal_topology.js;
  int dks = square.primal_topology.ks;

  {
    st11.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          for (int d = 0; d < ndims; ++d) {
            tw10_expected.data(d, k + 1 + dks, j + djs, i + dis, 0) *= -1;
          }
          compute_Hnm11<1, diff_ord>(
              tw10.data, st11.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
  }

  real errf = square.compute_Linf_error(tw10_expected, tw10, false);
  return errf;
}

void test_Hnm11_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hnm11 2 xy", uniform_vertical,
        compute_Hnm11_error<diff_order, vecfun_xy>, vecfun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hnm11 2 xyz", uniform_vertical,
        compute_Hnm11_error<diff_order, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(2, atol);
  }
}

template <int vdiff_ord, class F>
real compute_Hn0_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st20 = square.create_straight_form<2, 0>();
  square.primal_geometry.set_n0form_values(ic_fun, st20, 0);

  auto tw01 = square.create_twisted_form<0, 1>();
  auto tw01_expected = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_01form_values(ic_fun, tw01_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st20.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          tw01_expected.data(0, k + dks, j + djs, i + dis, 0) *= -1;
          compute_Hn0<1, vdiff_ord>(
              tw01.data, st20.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw01_expected, tw01, false);
  return errf;
}

void test_Hn0_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.14;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn0 2 z", uniform_vertical, compute_Hn0_error<vert_diff_ord, vecfun_z>,
        vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn0 2 xyz", uniform_vertical,
        compute_Hn0_error<vert_diff_ord, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(vert_diff_ord, atol);
  }
}

template <int diff_ord, class F>
real compute_Hnm11bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto tw11 = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_nm11form_values(ic_fun, tw11, 0);

  auto st10 = square.create_straight_form<1, 0>();
  auto st10_expected = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(ic_fun, st10_expected, 0);

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
          for (int d = 0; d < ndims; ++d) {
            st10_expected.data(d, k + pks, j + pjs, i + pis, 0) *= -1;
          }
          compute_Hnm11bar<1, diff_ord>(
              st10.data, tw11.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st10_expected, st10);
  return errf;
}

void test_Hnm11bar_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.11;

  {
    const int diff_order = 2;
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hnm11bar 2 xy", uniform_vertical,
        compute_Hnm11bar_error<diff_order, vecfun_xy>, vecfun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hnm11bar 2 xyz", uniform_vertical,
        compute_Hnm11bar_error<diff_order, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(1, atol);
  }
}

template <int vdiff_ord, class F>
real compute_Hn0bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto tw20 = square.create_twisted_form<2, 0>();
  square.dual_geometry.set_n0form_values(ic_fun, tw20, 0);

  auto st01 = square.create_straight_form<0, 1>();
  auto st01_expected = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(ic_fun, st01_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw20.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          st01_expected.data(0, k + pks, j + pjs, i + pis, 0) *= -1;
          compute_Hn0bar<1, vdiff_ord>(
              st01.data, tw20.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st01_expected, st01);
  return errf;
}

void test_Hn0bar_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.14;

  {
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn0bar 2 z", uniform_vertical,
        compute_Hn0bar_error<vert_diff_ord, vecfun_z>, vecfun_z{});
    conv_z.check_rate(vert_diff_ord, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn0bar 2 xyz", uniform_vertical,
        compute_Hn0bar_error<vert_diff_ord, vecfun_xyz>, vecfun_xyz{});
    conv_xyz.check_rate(vert_diff_ord, atol);
  }
}

template <int diff_ord, class F>
real compute_Hn1_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st21 = square.create_straight_form<2, 1>();
  square.primal_geometry.set_n1form_values(ic_fun, st21, 0);

  auto tw00 = square.create_twisted_form<0, 0>();
  auto tw00_expected = square.create_twisted_form<0, 0>();
  square.dual_geometry.set_00form_values(ic_fun, tw00_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st21.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hn1<1, diff_ord, vert_diff_ord>(
              tw00.data, st21.data, square.primal_geometry,
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
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hn1 2 xy", uniform_vertical, compute_Hn1_error<diff_order, fun_xy>,
        fun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("Hn1 2 z", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn1 2 xyz", uniform_vertical, compute_Hn1_error<diff_order, fun_xyz>,
        fun_xyz{});
    conv_xyz.check_rate(2, atol);
  }

  {
    const int diff_order = 4;
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hn1 4 xy", uniform_vertical, compute_Hn1_error<diff_order, fun_xy>,
        fun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("Hn1 4 z", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn1 4 xyz", uniform_vertical, compute_Hn1_error<diff_order, fun_xyz>,
        fun_xyz{});
    conv_xyz.check_rate(2, atol);
  }

  {
    const int diff_order = 6;
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hn1 6 xy", uniform_vertical, compute_Hn1_error<diff_order, fun_xy>,
        fun_xy{});
    conv_xy.check_rate(4, atol);
    auto conv_z =
        ConvergenceTest<nlevels>("Hn1 6 z", uniform_vertical,
                                 compute_Hn1_error<diff_order, fun_z>, fun_z{});
    conv_z.check_rate(2, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn1 6 xyz", uniform_vertical, compute_Hn1_error<diff_order, fun_xyz>,
        fun_xyz{});
    conv_xyz.check_rate(2, atol);
  }
}

template <int diff_ord, class F>
real compute_Hn1bar_error(int np, bool uniform_vertical, F ic_fun) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto tw21 = square.create_twisted_form<2, 1>();
  square.dual_geometry.set_n1form_values(ic_fun, tw21, 0);

  auto st00 = square.create_straight_form<0, 0>();
  auto st00_expected = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(ic_fun, st00_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  {
    tw21.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              st00.data, tw21.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st00_expected, st00);
  return errf;
}

void test_Hn1bar_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.13;

  {
    const int diff_order = 2;
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hn1bar 2 xy", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn1bar 2 z", uniform_vertical, compute_Hn1bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn1bar 2 xyz", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xyz>, fun_xyz{});
    conv_xyz.check_rate(1, atol);
  }

  {
    const int diff_order = 4;
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hn1bar 4 xy", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(diff_order, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn1bar 4 z", uniform_vertical, compute_Hn1bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn1bar 4 xyz", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xyz>, fun_xyz{});
    conv_xyz.check_rate(1, atol);
  }

  {
    const int diff_order = 6;
    auto conv_xy = ConvergenceTest<nlevels>(
        "Hn1bar 6 xy", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xy>, fun_xy{});
    conv_xy.check_rate(4, atol);
    auto conv_z = ConvergenceTest<nlevels>(
        "Hn1bar 6 z", uniform_vertical, compute_Hn1bar_error<diff_order, fun_z>,
        fun_z{});
    conv_z.check_rate(1, atol);
    auto conv_xyz = ConvergenceTest<nlevels>(
        "Hn1bar 6 xyz", uniform_vertical,
        compute_Hn1bar_error<diff_order, fun_xyz>, fun_xyz{});
    conv_xyz.check_rate(1, atol);
  }
}

int main() {
  yakl::init();

  for (bool uniform_vertical : {true, false}) {
    test_H00_convergence(uniform_vertical);
    test_H00bar_convergence(uniform_vertical);

    test_H10_convergence(uniform_vertical);
    test_H01_convergence(uniform_vertical);

    test_Hnm11_convergence(uniform_vertical);
    test_Hn0_convergence(uniform_vertical);

    test_Hnm11bar_convergence(uniform_vertical);
    test_Hn0bar_convergence(uniform_vertical);

    test_Hn1_convergence(uniform_vertical);
    test_Hn1bar_convergence(uniform_vertical);
  }

  yakl::finalize();
}
