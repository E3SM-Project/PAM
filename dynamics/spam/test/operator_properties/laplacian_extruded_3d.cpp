// clang-format off
unsigned constexpr ndims = 2;
#include "extruded_common.h"
#include "hodge_star.h"
#include "ext_deriv.h"
// clang-format on

using namespace pamc;

struct fun {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real fz = z * z * (z - 1) * (z - 1);
    return sx * sy * fz;
  }
};

struct lap_fun {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real fz = z * z * (z - 1) * (z - 1);
    real d2fz = 2 * (z - 1) * (z - 1) + 2 * z * z + 8 * z * (z - 1);
    return -8 * M_PI * M_PI * sx * sy * fz + sx * sy * d2fz;
  }
};

struct vecfun {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real cx = cos(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real cy = cos(2 * M_PI * y);
    real fz = z * z * (z - 1) * (z - 1);

    VecXYZ vvec;
    vvec.u = cx * sy * fz;
    vvec.v = sx * cy * fz;
    vvec.w = sx * sy * fz;
    // vvec.w = 0;

    // vvec.u = fz;
    // vvec.v = 0;
    // vvec.w = 0;
    return vvec;
  }
};

struct lap_vecfun {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real cx = cos(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real cy = cos(2 * M_PI * y);
    real fz = z * z * (z - 1) * (z - 1);
    real d2fz = 2 * (z - 1) * (z - 1) + 2 * z * z + 8 * z * (z - 1);

    VecXYZ vvec;
    vvec.u = (-8 * M_PI * M_PI * fz + d2fz) * cx * sy;
    vvec.v = (-8 * M_PI * M_PI * fz + d2fz) * sx * cy;
    vvec.w = (-8 * M_PI * M_PI * fz + d2fz) * sx * sy;

    // vvec.u = d2fz;
    // vvec.v = 0;
    // vvec.w = 0;
    return vvec;
  }
};

template <int diff_ord, int vert_diff_ord>
real compute_straight_00form_laplacian_error(int np, bool uniform_vertical) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(fun{}, st00, 0);

  auto lap_st00 = square.create_straight_form<0, 0>();
  auto lap_st00_expected = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(lap_fun{}, lap_st00_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st00.exchange();

    auto tmp_st10 = square.create_straight_form<1, 0>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0<1>(tmp_st10.data, st00.data, pis, pjs, pks, i, j, k, 0);
        });

    tmp_st10.exchange();
    auto tmp_st01 = square.create_straight_form<0, 1>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0_vert<1>(tmp_st01.data, st00.data, pis, pjs, pks, i, j, k,
                             0);
        });
    tmp_st01.exchange();

    auto tmp_tw11 = square.create_twisted_form<1, 1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H10<1, diff_ord>(tmp_tw11.data, tmp_st10.data,
                                   square.primal_geometry, square.dual_geometry,
                                   dis, djs, dks, i, j, k, 0);
        });
    tmp_tw11.exchange();

    auto tmp_tw20 = square.create_twisted_form<2, 0>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H01<1, diff_ord>(tmp_tw20.data, tmp_st01.data,
                                   square.primal_geometry, square.dual_geometry,
                                   dis, djs, dks, i, j, k + 1, 0);
        });
    tmp_tw20.exchange();
    tmp_tw20.set_bnd(0);

    auto tmp_tw21 = square.create_twisted_form<2, 1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Dnm1bar<1>(tmp_tw21.data, tmp_tw11.data, dis, djs, dks, i, j,
                             k, 0);
          compute_Dnm1bar_vert<1, ADD_MODE::ADD>(tmp_tw21.data, tmp_tw20.data,
                                                 dis, djs, dks, i, j, k, 0);
        });
    tmp_tw21.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              lap_st00.data, tmp_tw21.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(lap_st00_expected, lap_st00, false);
  return errf;
}

template <int diff_ord, int vert_diff_ord>
real compute_vel_laplacian_error(int np, bool uniform_vertical) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto st10 = square.create_straight_form<1, 0>();
  auto st01 = square.create_straight_form<0, 1>();
  square.primal_geometry.set_10form_values(vecfun{}, st10, 0);
  square.primal_geometry.set_01form_values(vecfun{}, st01, 0);

  auto lap_st10 = square.create_straight_form<1, 0>();
  auto lap_st01 = square.create_straight_form<0, 1>();
  auto lap_st10_expected = square.create_straight_form<1, 0>();
  auto lap_st01_expected = square.create_straight_form<0, 1>();
  square.primal_geometry.set_10form_values(lap_vecfun{}, lap_st10_expected, 0);
  square.primal_geometry.set_01form_values(lap_vecfun{}, lap_st01_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st10.exchange();
    st01.exchange();

    // *d*d
    auto tmp_st11 = square.create_straight_form<1, 1>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D1_ext<1>(tmp_st11.data, st10.data, st01.data, pis, pjs, pks,
                            i, j, k, 0);
        });
    tmp_st11.exchange();

    auto tmp_st20 = square.create_straight_form<2, 0>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D1<1>(tmp_st20.data, st10.data, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st20.exchange();

    auto tmp_tw10 = square.create_twisted_form<1, 0>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hnm11<1, diff_ord>(
              tmp_tw10.data, tmp_st11.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });
    tmp_tw10.exchange();
    tmp_tw10.set_bnd(0);

    auto tmp_tw01 = square.create_twisted_form<0, 1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hn0<1, diff_ord>(tmp_tw01.data, tmp_st20.data,
                                   square.primal_geometry, square.dual_geometry,
                                   dis, djs, dks, i, j, k, 0);
        });
    tmp_tw01.exchange();
    // std::cout << "TEST: "
    //           << yakl::intrinsics::minval(tmp_tw01.data)
    //           << " " << yakl::intrinsics::maxval(tmp_tw01.data) << std::endl;

    auto tmp_tw11 = square.create_twisted_form<1, 1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D1bar_ext<1>(tmp_tw11.data, tmp_tw10.data, tmp_tw01.data, dis,
                               djs, dks, i, j, k, 0);
        });
    tmp_tw11.exchange();

    auto tmp_tw20 = square.create_twisted_form<2, 0>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D1bar<1>(tmp_tw20.data, tmp_tw10.data, dis, djs, dks, i, j, k,
                           0);
        });
    tmp_tw20.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hnm11bar<1, diff_ord>(
              lap_st10.data, tmp_tw11.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
          for (int d = 0; d < ndims; ++d) {
            lap_st10.data(d, k + pks, j + pjs, i + pis, 0) *= -1;
          }
        });

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hn0bar<1, vert_diff_ord>(
              lap_st01.data, tmp_tw20.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
          lap_st01.data(0, k + pks, j + pjs, i + pis, 0) *= -1;
        });

    // d*d*
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H10<1, diff_ord>(tmp_tw11.data, st10.data,
                                   square.primal_geometry, square.dual_geometry,
                                   dis, djs, dks, i, j, k, 0);
        });
    tmp_tw11.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H01<1, vert_diff_ord>(
              tmp_tw20.data, st01.data, square.primal_geometry,
              square.dual_geometry, dis, djs, dks, i, j, k + 1, 0);
        });

    tmp_tw20.exchange();
    tmp_tw20.set_bnd(0);

    auto tmp_tw21 = square.create_twisted_form<2, 1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Dnm1bar<1>(tmp_tw21.data, tmp_tw11.data, dis, djs, dks, i, j,
                             k, 0);
          compute_Dnm1bar_vert<1, ADD_MODE::ADD>(tmp_tw21.data, tmp_tw20.data,
                                                 dis, djs, dks, i, j, k, 0);
        });
    tmp_tw21.exchange();

    auto tmp_st00 = square.create_straight_form<0, 0>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              tmp_st00.data, tmp_tw21.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st00.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0<1, ADD_MODE::ADD>(lap_st10.data, tmp_st00.data, pis, pjs,
                                       pks, i, j, k, 0);
        });

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0_vert<1, ADD_MODE::ADD>(lap_st01.data, tmp_st00.data, pis,
                                            pjs, pks, i, j, k, 0);
        });
  }

  real errf_st10 =
      square.compute_Linf_error(lap_st10_expected, lap_st10, false);
  real errf_st01 =
      square.compute_Linf_error(lap_st01_expected, lap_st01, false);

  return errf_st10; // + errf_st01;
}

template <int diff_ord, int vert_diff_ord>
real compute_twisted_21form_laplacian_error(int np, bool uniform_vertical) {
  ExtrudedUnitSquare square(np, 7 * np / 8, 9 * np / 8, uniform_vertical);

  auto tw21 = square.create_twisted_form<2, 1>();
  square.dual_geometry.set_n1form_values(fun{}, tw21, 0);

  auto lap_tw21 = square.create_twisted_form<2, 1>();
  auto lap_tw21_expected = square.create_twisted_form<2, 1>();
  square.dual_geometry.set_n1form_values(lap_fun{}, lap_tw21_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    tw21.exchange();

    auto tmp_st00 = square.create_straight_form<0, 0>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              tmp_st00.data, tw21.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st00.exchange();

    auto tmp_st10 = square.create_straight_form<1, 0>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0<1>(tmp_st10.data, tmp_st00.data, pis, pjs, pks, i, j, k,
                        0);
        });
    tmp_st10.exchange();

    auto tmp_st01 = square.create_straight_form<0, 1>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0_vert<1>(tmp_st01.data, tmp_st00.data, pis, pjs, pks, i, j,
                             k, 0);
        });
    tmp_st01.exchange();

    auto tmp_tw11 = square.create_twisted_form<1, 1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H10<1, diff_ord>(tmp_tw11.data, tmp_st10.data,
                                   square.primal_geometry, square.dual_geometry,
                                   dis, djs, dks, i, j, k, 0);
        });
    tmp_tw11.exchange();

    auto tmp_tw20 = square.create_twisted_form<2, 0>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni - 2,
                        square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H01<1, diff_ord>(tmp_tw20.data, tmp_st01.data,
                                   square.primal_geometry, square.dual_geometry,
                                   dis, djs, dks, i, j, k + 1, 0);
        });
    tmp_tw20.exchange();
    tmp_tw20.set_bnd(0);

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Dnm1bar<1>(lap_tw21.data, tmp_tw11.data, dis, djs, dks, i, j,
                             k, 0);
          compute_Dnm1bar_vert<1, ADD_MODE::ADD>(lap_tw21.data, tmp_tw20.data,
                                                 dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(lap_tw21_expected, lap_tw21, false);
  return errf;
}

void test_laplacian_convergence(bool uniform_vertical) {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_ord = 2;
    const int vert_diff_ord = 2;
    auto conv_st00 = ConvergenceTest<nlevels>(
        "straight 00 form laplacian", uniform_vertical,
        compute_straight_00form_laplacian_error<diff_ord, vert_diff_ord>);
    conv_st00.check_rate(vert_diff_ord, atol);
  }

  {
    const int diff_ord = 2;
    const int vert_diff_ord = 2;
    auto conv_vel = ConvergenceTest<nlevels>(
        "velocity laplacian", uniform_vertical,
        compute_vel_laplacian_error<diff_ord, vert_diff_ord>);
    conv_vel.check_rate(vert_diff_ord, atol);
  }

  {
    const int diff_ord = 2;
    const int vert_diff_ord = 2;
    auto conv_tw21 = ConvergenceTest<nlevels>(
        "twisted 21 form laplacian", uniform_vertical,
        compute_twisted_21form_laplacian_error<diff_ord, vert_diff_ord>);
    conv_tw21.check_rate(vert_diff_ord, atol);
  }
}

int main() {
  yakl::init();

  for (bool uniform_vertical : {true, false}) {
    test_laplacian_convergence(uniform_vertical);
  }
  yakl::finalize();
}
