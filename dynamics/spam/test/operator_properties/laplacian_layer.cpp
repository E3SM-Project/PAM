// clang-format off
#include "layer_common.h"
#include "hodge_star.h"
#include "ext_deriv.h"
// clang-format on

struct fun {
  real YAKL_INLINE operator()(real x, real y) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    return sx * sy;
  }
};

struct lap_fun {
  real YAKL_INLINE operator()(real x, real y) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    return -8 * M_PI * M_PI * sx * sy;
  }
};

struct vecfun {
  VecXY YAKL_INLINE operator()(real x, real y) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real cx = cos(2 * M_PI * x);
    real cy = cos(2 * M_PI * y);

    VecXY vvec;
    vvec.u = sx * sy;
    vvec.v = cx * cy;
    return vvec;
  }
};

struct lap_vecfun {
  VecXY YAKL_INLINE operator()(real x, real y) const {
    real sx = sin(2 * M_PI * x);
    real sy = sin(2 * M_PI * y);
    real cx = cos(2 * M_PI * x);
    real cy = cos(2 * M_PI * y);

    VecXY vvec;
    vvec.u = -8 * M_PI * M_PI * sx * sy;
    vvec.v = -8 * M_PI * M_PI * cx * cy;
    return vvec;
  }
};

template <int diff_ord> real compute_straight_0form_laplacian_error(int np) {
  PeriodicUnitSquare square(np, 2 * np);

  auto st0 = square.create_straight_form<0>();
  square.primal_geometry.set_0form_values(fun{}, st0, 0);

  auto lap_st0 = square.create_straight_form<0>();
  auto lap_st0_expected = square.create_straight_form<0>();
  square.primal_geometry.set_0form_values(lap_fun{}, lap_st0_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st0.exchange();

    auto tmp_st1 = square.create_straight_form<1>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0<1>(tmp_st1.data, st0.data, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st1.exchange();

    auto tmp_tw1 = square.create_twisted_form<1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H1<1, diff_ord>(tmp_tw1.data, tmp_st1.data,
                                  square.primal_geometry, square.dual_geometry,
                                  dis, djs, dks, i, j, k, 0);
        });
    tmp_tw1.exchange();

    auto tmp_tw2 = square.create_twisted_form<2>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Dnm1bar<1>(tmp_tw2.data, tmp_tw1.data, dis, djs, dks, i, j, k,
                             0);
        });
    tmp_tw2.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H2bar<1, diff_ord>(
              lap_st0.data, tmp_tw2.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(lap_st0_expected, lap_st0);
  return errf;
}

template <int diff_ord> real compute_straight_1form_laplacian_error(int np) {
  PeriodicUnitSquare square(np, 2 * np);

  auto st1 = square.create_straight_form<1>();
  square.primal_geometry.set_1form_values(vecfun{}, st1, 0,
                                          LINE_INTEGRAL_TYPE::TANGENT);

  auto lap_st1 = square.create_straight_form<1>();
  auto lap_st1_expected = square.create_straight_form<1>();
  square.primal_geometry.set_1form_values(lap_vecfun{}, lap_st1_expected, 0,
                                          LINE_INTEGRAL_TYPE::TANGENT);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    st1.exchange();

    // *d*d
    auto tmp_st2 = square.create_straight_form<2>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D1<1>(tmp_st2.data, st1.data, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st2.exchange();

    auto tmp_tw0 = square.create_twisted_form<0>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H2<1, diff_ord>(tmp_tw0.data, tmp_st2.data,
                                  square.primal_geometry, square.dual_geometry,
                                  dis, djs, dks, i, j, k, 0);
        });
    tmp_tw0.exchange();

    auto tmp_tw1 = square.create_twisted_form<1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0bar<1>(tmp_tw1.data, tmp_tw0.data, dis, djs, dks, i, j, k,
                           0);
        });
    tmp_tw1.exchange();

    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H1bar<1, diff_ord>(
              lap_st1.data, tmp_tw1.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });

    // d*d*
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H1<1, diff_ord>(tmp_tw1.data, st1.data,
                                  square.primal_geometry, square.dual_geometry,
                                  pis, pjs, pks, i, j, k, 0);
        });
    tmp_tw1.exchange();

    auto tmp_tw2 = square.create_twisted_form<2>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Dnm1bar<1>(tmp_tw2.data, tmp_tw1.data, dis, djs, dks, i, j, k,
                             0);
        });
    tmp_tw2.exchange();

    auto tmp_st0 = square.create_straight_form<0>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H2bar<1, diff_ord>(
              tmp_st0.data, tmp_tw2.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st0.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0<1, ADD_MODE::ADD>(lap_st1.data, tmp_st0.data, dis, djs,
                                       dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(lap_st1_expected, lap_st1);
  return errf;
}

template <int diff_ord> real compute_twisted_2form_laplacian_error(int np) {
  PeriodicUnitSquare square(np, 2 * np);

  auto tw2 = square.create_twisted_form<2>();
  square.dual_geometry.set_2form_values(fun{}, tw2, 0);

  auto lap_tw2 = square.create_twisted_form<2>();
  auto lap_tw2_expected = square.create_twisted_form<2>();
  square.dual_geometry.set_2form_values(lap_fun{}, lap_tw2_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;

  {
    tw2.exchange();

    auto tmp_st0 = square.create_straight_form<0>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H2bar<1, diff_ord>(
              tmp_st0.data, tw2.data, square.primal_geometry,
              square.dual_geometry, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st0.exchange();

    auto tmp_st1 = square.create_straight_form<1>();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0<1>(tmp_st1.data, tmp_st0.data, pis, pjs, pks, i, j, k, 0);
        });
    tmp_st1.exchange();

    auto tmp_tw1 = square.create_twisted_form<1>();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_H1<1, diff_ord>(tmp_tw1.data, tmp_st1.data,
                                  square.primal_geometry, square.dual_geometry,
                                  dis, djs, dks, i, j, k, 0);
        });
    tmp_tw1.exchange();

    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Dnm1bar<1>(lap_tw2.data, tmp_tw1.data, dis, djs, dks, i, j, k,
                             0);
        });
  }

  real errf = square.compute_Linf_error(lap_tw2_expected, lap_tw2);
  return errf;
}

void test_laplacian_convergence() {
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_st0 = ConvergenceTest<nlevels>(
        "straight 0 form laplacian",
        compute_straight_0form_laplacian_error<diff_order>);
    conv_st0.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 2;
    auto conv_st1 = ConvergenceTest<nlevels>(
        "straight 1 form laplacian",
        compute_straight_1form_laplacian_error<diff_order>);
    conv_st1.check_rate(diff_order, atol);
  }

  {
    const int diff_order = 2;
    auto conv_tw2 = ConvergenceTest<nlevels>(
        "twisted 2 form laplacian",
        compute_twisted_2form_laplacian_error<diff_order>);
    conv_tw2.check_rate(diff_order, atol);
  }
}

int main() {
  yakl::init();
  test_laplacian_convergence();
  yakl::finalize();
}
