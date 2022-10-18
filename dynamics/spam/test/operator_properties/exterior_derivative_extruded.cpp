// clang-format off
#include "extruded_common.h"
#include "ext_deriv.h"
// clang-format on

struct fun {
  real YAKL_INLINE operator()(real x, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = sin(2 * M_PI * z);
    return sx * sz;
  }
};

struct grad_fun {
  vecext<2> YAKL_INLINE operator()(real x, real z) const {
    vecext<2> vvec;
    vvec.u = 2 * M_PI * cos(2 * M_PI * x) * sin(2 * M_PI * z);
    vvec.w = 2 * M_PI * sin(2 * M_PI * x) * cos(2 * M_PI * z);
    return vvec;
  }
};

struct vecfun {
  vecext<2> YAKL_INLINE operator()(real x, real z) const {
    vecext<2> vvec;
    vvec.u = sin(2 * M_PI * x) * sin(2 * M_PI * z);
    vvec.w = sin(2 * M_PI * x) * cos(2 * M_PI * z);
    return vvec;
  }
};

struct div_vecfun {
  real YAKL_INLINE operator()(real x, real z) const {
    return 2 * M_PI *
           (cos(2 * M_PI * x) * sin(2 * M_PI * z) -
            sin(2 * M_PI * x) * sin(2 * M_PI * z));
  }
};

struct curl_vecfun {
  real YAKL_INLINE operator()(real x, real z) const {
    return 2 * M_PI *
           (cos(2 * M_PI * x) * cos(2 * M_PI * z) -
            sin(2 * M_PI * x) * cos(2 * M_PI * z));
  }
};

void test_D1(int np, real atol) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(fun{}, st00, 0);

  auto st10 = square.create_straight_form<1, 0>();
  auto st10_expected = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(grad_fun{}, st10_expected, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;
  {
    st00.exchange();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          SArray<real, 1, 1> c;
          c(0) = 1;
          compute_cwD1<1>(st10.data, c, st00.data, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st10_expected, st10);

  if (errf > atol) {
    std::cout << "Exactness of D1 failed, error = " << errf << " tol = " << atol
              << std::endl;
    exit(-1);
  }
}

void test_Dv(int np, real atol) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(fun{}, st00, 0);

  auto st01 = square.create_straight_form<0, 1>();
  auto st01_expected = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(grad_fun{}, st01_expected, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  auto ones = square.create_twisted_form<1, 0>();
  ones.set(1);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;
  {
    st00.exchange();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.ni,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_wDv<1>(st01.data, ones.data, st00.data, pis, pjs, pks, i, j,
                         k, 0);
        });
  }

  real errf = square.compute_Linf_error(st01_expected, st01);

  if (errf > atol) {
    std::cout << "Exactness of Dv failed, error = " << errf << " tol = " << atol
              << std::endl;
    exit(-1);
  }
}

void test_Dxz(int np, real atol) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto st10 = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(vecfun{}, st10, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  auto st01 = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(vecfun{}, st01, 0,
                                           LINE_INTEGRAL_TYPE::TANGENT);

  auto st11 = square.create_straight_form<1, 1>();
  auto st11_expected = square.create_straight_form<1, 1>();
  square.primal_geometry.set_11form_values(curl_vecfun{}, st11_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;
  {
    st01.exchange();
    st10.exchange();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_Dxz<1>(st11.data, st10.data, st01.data, pis, pjs, pks, i, j,
                         k, 0);
        });
  }

  real errf = square.compute_Linf_error(st11_expected, st11);

  if (errf > atol) {
    std::cout << "Exactness of Dxz failed, error = " << errf
              << " tol = " << atol << std::endl;
    exit(-1);
  }
}

void test_Dbar2_and_Dvbar(int np, real atol) {
  ExtrudedUnitSquare square(np, 2 * np);

  auto tw01 = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_01form_values(vecfun{}, tw01, 0,
                                         LINE_INTEGRAL_TYPE::NORMAL);

  auto tw10 = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_10form_values(vecfun{}, tw10, 0,
                                         LINE_INTEGRAL_TYPE::NORMAL);

  auto tw11 = square.create_twisted_form<1, 1>();
  auto tw11_expected = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_11form_values(div_vecfun{}, tw11_expected, 0);

  auto ones = square.create_twisted_form<1, 0>();
  ones.set(1);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;
  {
    tw01.exchange();
    tw10.exchange();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_cwDbar2<1>(tw11.data, 1, tw01.data, dis, djs, dks, i, j, k,
                             0);
          compute_wDvbar<1, ADD_MODE::ADD>(tw11.data, ones.data, tw10.data, dis,
                                           djs, dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw11_expected, tw11);

  if (errf > atol) {
    std::cout << "Exactness of Dbar2 and Dvbar failed, error = " << errf
              << " tol = " << atol << std::endl;
    exit(-1);
  }
}

int main() {
  yakl::init();
  real atol = 500 * std::numeric_limits<real>::epsilon();
  test_D1(33, atol);
  test_Dv(33, atol);
  test_Dbar2_and_Dvbar(33, atol);
  test_Dxz(33, atol);
  yakl::finalize();
}
