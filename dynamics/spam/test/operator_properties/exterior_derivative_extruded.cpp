// clang-format off
unsigned constexpr ndims = 1;
#include "extruded_common.h"
#include "ext_deriv.h"
// clang-format on

using namespace pamc;

struct fun {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    real sx = sin(2 * M_PI * x);
    real sz = sin(2 * M_PI * z);
    return sx * sz;
  }
};

struct grad_fun {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    VecXYZ vvec;
    vvec.u = 2 * M_PI * cos(2 * M_PI * x) * sin(2 * M_PI * z);
    vvec.w = 2 * M_PI * sin(2 * M_PI * x) * cos(2 * M_PI * z);
    return vvec;
  }
};

struct vecfun {
  VecXYZ YAKL_INLINE operator()(real x, real y, real z) const {
    VecXYZ vvec;
    vvec.u = sin(2 * M_PI * x) * sin(2 * M_PI * z);
    vvec.w = sin(2 * M_PI * x) * cos(2 * M_PI * z);
    return vvec;
  }
};

struct div_vecfun {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    return 2 * M_PI *
           (cos(2 * M_PI * x) * sin(2 * M_PI * z) -
            sin(2 * M_PI * x) * sin(2 * M_PI * z));
  }
};

struct curl_vecfun {
  real YAKL_INLINE operator()(real x, real y, real z) const {
    return 2 * M_PI *
           (sin(2 * M_PI * x) * cos(2 * M_PI * z) -
            cos(2 * M_PI * x) * cos(2 * M_PI * z));
  }
};

void test_D0(int np, bool uniform_vertical, real atol) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(fun{}, st00, 0);

  auto st10 = square.create_straight_form<1, 0>();
  auto st10_expected = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(grad_fun{}, st10_expected, 0);

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
          compute_D0<1>(st10.data, st00.data, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st10_expected, st10);

  if (errf > atol) {
    std::cout << "Exactness of D0 failed, error = " << errf << " tol = " << atol
              << std::endl;
    exit(-1);
  }
}

void test_D0_vert(int np, bool uniform_vertical, real atol) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st00 = square.create_straight_form<0, 0>();
  square.primal_geometry.set_00form_values(fun{}, st00, 0);

  auto st01 = square.create_straight_form<0, 1>();
  auto st01_expected = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(grad_fun{}, st01_expected, 0);

  int pis = square.primal_topology.is;
  int pjs = square.primal_topology.js;
  int pks = square.primal_topology.ks;
  {
    st00.exchange();
    parallel_for(
        SimpleBounds<3>(square.primal_topology.nl,
                        square.primal_topology.n_cells_y,
                        square.primal_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0_vert<1>(st01.data, st00.data, pis, pjs, pks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(st01_expected, st01);

  if (errf > atol) {
    std::cout << "Exactness of D0_vert failed, error = " << errf
              << " tol = " << atol << std::endl;
    exit(-1);
  }
}

void test_D0bar(int np, bool uniform_vertical, real atol) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw00 = square.create_twisted_form<0, 0>();
  square.dual_geometry.set_00form_values(fun{}, tw00, 0);

  auto tw10 = square.create_twisted_form<1, 0>();
  auto tw10_expected = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_10form_values(grad_fun{}, tw10_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;
  {
    tw00.exchange();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.ni, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0bar_ext<1>(tw10.data, tw00.data, dis, djs, dks, i, j, k, 0);
          tw10.data(0, k + dks, j + djs, i + dis, 0) *= -1;
        });
  }

  real errf = square.compute_Linf_error(tw10_expected, tw10);

  if (errf > atol) {
    std::cout << "Exactness of D0bar failed, error = " << errf
              << " tol = " << atol << std::endl;
    exit(-1);
  }
}

void test_D0bar_vert(int np, bool uniform_vertical, real atol) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw00 = square.create_twisted_form<0, 0>();
  square.dual_geometry.set_00form_values(fun{}, tw00, 0);

  auto tw01 = square.create_twisted_form<0, 1>();
  auto tw01_expected = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_01form_values(grad_fun{}, tw01_expected, 0);

  int dis = square.dual_topology.is;
  int djs = square.dual_topology.js;
  int dks = square.dual_topology.ks;
  {
    tw00.exchange();
    parallel_for(
        SimpleBounds<3>(square.dual_topology.nl, square.dual_topology.n_cells_y,
                        square.dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D0bar_vert<1>(tw01.data, tw00.data, dis, djs, dks, i, j, k,
                                0);
          tw01.data(0, k + dks, j + djs, i + dis, 0) *= -1;
        });
  }

  real errf = square.compute_Linf_error(tw01_expected, tw01);

  if (errf > atol) {
    std::cout << "Exactness of D0bar_vert failed, error = " << errf
              << " tol = " << atol << std::endl;
    exit(-1);
  }
}

void test_D1_ext(int np, bool uniform_vertical, real atol) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto st10 = square.create_straight_form<1, 0>();
  square.primal_geometry.set_10form_values(vecfun{}, st10, 0);

  auto st01 = square.create_straight_form<0, 1>();
  square.primal_geometry.set_01form_values(vecfun{}, st01, 0);

  auto st11 = square.create_straight_form<1, 1>();
  auto st11_expected = square.create_straight_form<1, 1>();
  square.primal_geometry.set_n1form_values(curl_vecfun{}, st11_expected, 0);

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
          compute_D1_ext<1>(st11.data, st10.data, st01.data, pis, pjs, pks, i,
                            j, k, 0);
          st11.data(0, k + pks, j + pjs, i + pis, 0) *= -1;
        });
  }

  real errf = square.compute_Linf_error(st11_expected, st11);

  if (errf > atol) {
    std::cout << "Exactness of D1_ext failed, error = " << errf
              << " tol = " << atol << std::endl;
    exit(-1);
  }
}

void test_Dnm1bar_and_Dnm1bar_vert(int np, bool uniform_vertical, real atol) {
  ExtrudedUnitSquare square(np, 2 * np, uniform_vertical);

  auto tw01 = square.create_twisted_form<0, 1>();
  square.dual_geometry.set_nm11form_values(vecfun{}, tw01, 0);

  auto tw10 = square.create_twisted_form<1, 0>();
  square.dual_geometry.set_n0form_values(vecfun{}, tw10, 0);

  auto tw11 = square.create_twisted_form<1, 1>();
  auto tw11_expected = square.create_twisted_form<1, 1>();
  square.dual_geometry.set_n1form_values(div_vecfun{}, tw11_expected, 0);

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
          compute_Dnm1bar<1>(tw11.data, tw01.data, dis, djs, dks, i, j, k, 0);
          compute_Dnm1bar_vert<1, ADD_MODE::ADD>(tw11.data, tw10.data, dis, djs,
                                                 dks, i, j, k, 0);
        });
  }

  real errf = square.compute_Linf_error(tw11_expected, tw11);

  if (errf > atol) {
    std::cout << "Exactness of Dnm1bar and Dnm1bar_vert failed, error = "
              << errf << " tol = " << atol << std::endl;
    exit(-1);
  }
}

int main() {
  yakl::init();
  real atol = 500 * std::numeric_limits<real>::epsilon();

  for (bool uniform_vertical : {true, false}) {
    test_D0(33, uniform_vertical, atol);
    test_D0_vert(33, uniform_vertical, atol);
    test_D0bar(33, uniform_vertical, atol);
    test_D0bar_vert(33, uniform_vertical, atol);
    test_Dnm1bar_and_Dnm1bar_vert(33, uniform_vertical, atol);
    test_D1_ext(33, uniform_vertical, atol);
  }
  yakl::finalize();
}
