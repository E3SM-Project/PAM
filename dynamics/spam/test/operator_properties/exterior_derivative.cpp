// clang-format off
#include "layer_common.h"
#include "ext_deriv.h"
// clang-format on

real YAKL_INLINE fun(real x, real y) {
  real sx = sin(2 * M_PI * x);
  real sy = sin(2 * M_PI * y);
  return sx * sy;
}
vec<2> YAKL_INLINE grad_fun(real x, real y) {
  vec<2> vvec;
  vvec.u = 2 * M_PI * cos(2 * M_PI * x) * sin(2 * M_PI * y);
  vvec.v = 2 * M_PI * sin(2 * M_PI * x) * cos(2 * M_PI * y);
  return vvec;
}

vec<2> YAKL_INLINE vecfun(real x, real y) {
  vec<2> vvec;
  vvec.u = sin(2 * M_PI * x) * sin(2 * M_PI * y);
  vvec.v = sin(2 * M_PI * x) * cos(2 * M_PI * y);
  return vvec;
}
real YAKL_INLINE div_vecfun(real x, real y) {
  return 2 * M_PI *
         (cos(2 * M_PI * x) * sin(2 * M_PI * y) -
          sin(2 * M_PI * x) * sin(2 * M_PI * y));
}

real YAKL_INLINE curl_vecfun(real x, real y) {
  return 2 * M_PI *
         (cos(2 * M_PI * x) * cos(2 * M_PI * y) -
          sin(2 * M_PI * x) * cos(2 * M_PI * y));
}

void test_D1(int np, real atol) {
  PeriodicUnitSquare square(np, 2 * np);

  auto st0 = square.create_straight_form<0>();
  square.primal_geometry.set_0form_values(fun, st0, 0);

  auto st1 = square.create_straight_form<1>();
  auto st1_expected = square.create_straight_form<1>();
  square.primal_geometry.set_1form_values(grad_fun, st1_expected, 0,
                                          LINE_INTEGRAL_TYPE::TANGENT);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;
  {
    st0.exchange();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          SArray<real, 1, 1> c;
          c(0) = 1;
          compute_cwD1<1>(st1.data, c, st0.data, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st1_expected, st1);

  if (errf > atol) {
    std::cout << "Exactness of D1 failed, error = " << errf << " tol = " << atol
              << std::endl;
    exit(-1);
  }
}

void test_D2(int np, real atol) {
  PeriodicUnitSquare square(np, 2 * np);

  auto st1 = square.create_straight_form<1>();
  square.primal_geometry.set_1form_values(vecfun, st1, 0,
                                          LINE_INTEGRAL_TYPE::TANGENT);

  auto st2 = square.create_straight_form<2>();
  auto st2_expected = square.create_straight_form<2>();
  square.primal_geometry.set_2form_values(curl_vecfun, st2_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;
  {
    st1.exchange();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D2<1>(st2.data, st1.data, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st2_expected, st2);

  if (errf > atol) {
    std::cout << "Exactness of D2 failed, error = " << errf << " tol = " << atol
              << std::endl;
    exit(-1);
  }
}

void test_Dbar2(int np, real atol) {
  PeriodicUnitSquare square(np, 2 * np);

  auto tw1 = square.create_twisted_form<1>();
  square.dual_geometry.set_1form_values(vecfun, tw1, 0,
                                        LINE_INTEGRAL_TYPE::NORMAL);

  auto tw2 = square.create_twisted_form<2>();
  auto tw2_expected = square.create_twisted_form<2>();
  square.dual_geometry.set_2form_values(div_vecfun, tw2_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;
  {
    tw1.exchange();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_cwDbar2<1>(tw2.data, 1, tw1.data, dis, djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw2_expected, tw2);

  if (errf > atol) {
    std::cout << "Exactness of Dbar2 failed, error = " << errf
              << " tol = " << atol << std::endl;
    exit(-1);
  }
}

int main() {
  yakl::init();
  real atol = 500 * std::numeric_limits<real>::epsilon();
  test_D1(33, atol);
  test_Dbar2(33, atol);
  test_D2(33, atol);
  yakl::finalize();
}
