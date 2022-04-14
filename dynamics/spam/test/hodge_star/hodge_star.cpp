
#include <iostream>
#include <array>
#include "pam_const.h"

uint constexpr ndims = 2;
uint constexpr nprognostic = 0;
uint constexpr nconstant = 0;
uint constexpr nauxiliary = 0;
uint constexpr ndiagnostic = 0;

#include "common.h"

struct ModelParameters : public Parameters {};

#include "params.h"
#include "fields.h"
#include "topology.h"
#include "geometry.h"
#include "exchange.h"
#include "model.h"
#include "hodge_star.h"



real YAKL_INLINE fun_x(real x, real y)
{
  return sin(2 * M_PI * x);
}

real YAKL_INLINE fun_y(real x, real y)
{
  return sin(2 * M_PI * y);
}

real YAKL_INLINE fun_xy(real x, real y)
{
  real sx = sin(2 * M_PI * x);
  real sy = sin(2 * M_PI * y);
  return sx * sy;
}


vec<2> YAKL_INLINE vecfun_x(real x, real y)
{
  real sx = sin(2 * M_PI * x);

  vec<2> vvec;
  vvec.u = sx;
  vvec.v = 0;
  return vvec;
}

vec<2> YAKL_INLINE vecfun_y(real x, real y)
{
  real sy = sin(2 * M_PI * y);

  vec<2> vvec;
  vvec.u = 0;
  vvec.v = sy;
  return vvec;
}

vec<2> YAKL_INLINE vecfun_xy(real x, real y)
{
  real sx = sin(2 * M_PI * x);
  real sy = sin(2 * M_PI * y);

  vec<2> vvec;
  vvec.u = sx * sy;
  vvec.v = sx * sy * sx * sy;
  return vvec;
}

template<int nlevels>
struct ConvergenceTest
{
  std::string name;
  std::array<real, nlevels> errs;
  std::array<real, nlevels-1> rates;

  template<typename F1, typename F2>
  ConvergenceTest(const std::string &name, F1 error_fun, F2 ic_fun) : name(name)
  {
    for (int l = 0; l < nlevels; ++l)
    {
      int np = 8 * std::pow(2, l - 1);
      errs[l] = error_fun(np, ic_fun);
    }

    for (int l = 0; l < nlevels; ++l)
    {
      if (l < nlevels - 1)
      {
        rates[l] = std::log2(errs[l] / errs[l+1]);
      }
    }

  }

  void print_errors_and_rates() const
  {
    std::cout << "Errors" << std::endl;
    for (int l = 0; l < nlevels; ++l)
    {
      std::cout << l << " " << errs[l] << std::endl;
    }
    std::cout << "Rates" << std::endl;
    for (int l = 0; l < nlevels-1; ++l)
    {
      std::cout << l << " " << rates[l] << std::endl;
    }
  }

  void check_rate(real expected_rate, real atol) const
  {
    auto rate = rates[nlevels - 2];
    if (std::abs(rate - expected_rate) > atol)
    {
      std::cout << "Failed convergence test: " << name << std::endl;
      std::cout << "Achieved rate: " << rate << " "
                << "Expected rate: " << expected_rate << std::endl;
      print_errors_and_rates();
      exit(-1);
    }
  }
};


Parallel parallel_stub(int nx, int ny)
{
  Parallel par;
  par.nx_glob = par.nx = nx;
  par.ny_glob = par.ny = ny;
  par.nz = 1;

  par.nprocx = par.nprocy = par.nranks = 1;
  par.masterproc = par.myrank = par.px = par.py = 0;

  par.i_beg = 0;
  par.i_end = nx - 1;
  par.j_beg = 0;
  par.j_end = ny - 1;

  par.halox = maxhalosize;
  par.haloy = maxhalosize;
  par.x_neigh(0) = par.x_neigh(1) = par.y_neigh(0) = par.y_neigh(1) = 0;
  par.ll_neigh = par.ur_neigh = par.ul_neigh = par.lr_neigh = 0;

  par.nens = 1;

  return par;
}

struct PeriodicUnitSquare
{
  Topology primal_topology;
  Topology dual_topology;

  UniformRectangularStraightGeometry primal_geometry;
  UniformRectangularTwistedGeometry dual_geometry;

  PeriodicUnitSquare(int nx, int ny)
  {
    ModelParameters params;

    params.nx_glob = nx;
    params.ny_glob = ny;
    params.xlen = 1;
    params.xc = 0;
    params.ylen = 1;
    params.yc = 0;

    Parallel par = parallel_stub(nx, ny);

    primal_topology.initialize(par, true);
    dual_topology.initialize(par, false);

    primal_geometry.initialize(primal_topology, params);

    dual_geometry.initialize(dual_topology, params);
  }
};

template<int diff_ord>
real compute_I_error(int np, real (*ic_fun)(real, real))
{
    PeriodicUnitSquare square(np, np);

    // set up one twisted 2-form
    Field tw2;
    tw2.initialize(square.dual_topology, "twisted 2-form", ndims, 1, 1);
    Exchange tw2_exchng;
    tw2_exchng.initialize(square.dual_topology, ndims, 1, 1);
    square.dual_geometry.set_2form_values(ic_fun, tw2, 0);

    // set up three straight 0-forms
    // result, expected, and error
    Field st0, st0_expected, st0_err;
    st0.initialize(square.primal_topology, "result straight 0-form", 0, 1, 1);
    st0_expected.initialize(square.primal_topology, "expected straight 0-form", 0, 1, 1);
    st0_err.initialize(square.primal_topology, "error straight 0-form", 0, 1, 1);
    square.primal_geometry.set_0form_values(ic_fun, st0_expected, 0);


    int pis = square.primal_topology.is;
    int pjs = square.primal_topology.js;
    int pks = square.primal_topology.ks;

    int dis = square.dual_topology.is;
    int djs = square.dual_topology.js;
    int dks = square.dual_topology.ks;

    {
      tw2_exchng.exchange_field(tw2);

      //this shouldnt be necessary ?
      Geometry *pgeom = &square.primal_geometry;
      Geometry *dgeom = &square.dual_geometry;

      parallel_for(Bounds<3>(square.primal_topology.nl,
                             square.primal_topology.n_cells_y,
                             square.primal_topology.n_cells_x),
      YAKL_LAMBDA(int k, int j, int i)
      {
          compute_I<1, diff_ord>(st0.data, tw2.data, *pgeom, *dgeom, pis, pjs, pks, i, j, k, 0);
      });
    }

    st0_err.set(0);
    parallel_for(Bounds<3>(square.primal_topology.nl,
                           square.primal_topology.n_cells_y,
                           square.primal_topology.n_cells_x),
    YAKL_LAMBDA(int k, int j, int i)
    {
      st0_err.data(0, pks + k, pjs + j, pis + i) =
        abs(st0.data(0, pks + k, pjs + j, pis + i) - st0_expected.data(0, pks + k, pjs + j, pis + i));
    });

    real errf = yakl::intrinsics::maxval(st0_err.data);

    return errf;
}

void test_I_convergence()
{
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>("I 2 x", compute_I_error<diff_order>, fun_x);
    conv_x.check_rate(2, atol);
    auto conv_y = ConvergenceTest<nlevels>("I 2 y", compute_I_error<diff_order>, fun_y);
    conv_y.check_rate(2, atol);
    auto conv_xy = ConvergenceTest<nlevels>("I 2 xy", compute_I_error<diff_order>, fun_xy);
    conv_xy.check_rate(2, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>("I 4 x", compute_I_error<diff_order>, fun_x);
    conv_x.check_rate(4, atol);
    auto conv_y = ConvergenceTest<nlevels>("I 4 y", compute_I_error<diff_order>, fun_y);
    conv_y.check_rate(4, atol);
    auto conv_xy = ConvergenceTest<nlevels>("I 4 xy", compute_I_error<diff_order>, fun_xy);
    conv_xy.check_rate(4, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>("I 6 x", compute_I_error<diff_order>, fun_x);
    conv_x.check_rate(6, atol);
    auto conv_y = ConvergenceTest<nlevels>("I 6 y", compute_I_error<diff_order>, fun_y);
    conv_y.check_rate(6, atol);
    auto conv_xy = ConvergenceTest<nlevels>("I 6 xy", compute_I_error<diff_order>, fun_xy);
    conv_xy.check_rate(4, atol);
  }
}

template<int diff_ord>
real compute_H_error(int np, vec<2> (*ic_fun)(real, real))
{
    PeriodicUnitSquare square(np, np);

    // set up one straight 1-form
    Field st1;
    st1.initialize(square.primal_topology, "straight 1-form", 1, 1, 1);
    Exchange st1_exchng;
    st1_exchng.initialize(square.primal_topology, 1, 1, 1);
    square.primal_geometry.set_1form_values(ic_fun, st1, 0, LINE_INTEGRAL_TYPE::TANGENT);

    // set up three twisted 1-forms
    // result, expected, and error
    Field tw1, tw1_expected, tw1_err;
    tw1.initialize(square.dual_topology, "result twisted 1-form", 1, 1, 1);
    tw1_expected.initialize(square.dual_topology, "expected twisted 1-form", 1, 1, 1);
    tw1_err.initialize(square.dual_topology, "error twisted 1-form", 1, 1, 1);
    square.dual_geometry.set_1form_values(ic_fun, tw1_expected, 0,  LINE_INTEGRAL_TYPE::NORMAL);


    int pis = square.primal_topology.is;
    int pjs = square.primal_topology.js;
    int pks = square.primal_topology.ks;

    int dis = square.dual_topology.is;
    int djs = square.dual_topology.js;
    int dks = square.dual_topology.ks;


    {
      st1_exchng.exchange_field(st1);

      //this shouldnt be necessary ?
      Geometry *pgeom = &square.primal_geometry;
      Geometry *dgeom = &square.dual_geometry;

      parallel_for(Bounds<3>(square.dual_topology.nl,
                             square.dual_topology.n_cells_y,
                             square.dual_topology.n_cells_x),
      YAKL_LAMBDA(int k, int j, int i)
      {
        compute_H<1, diff_ord>(tw1.data, st1.data, *pgeom, *dgeom, dis, djs, dks, i, j, k, 0);
      });
    }

    tw1_err.set(0);
    parallel_for(Bounds<3>(square.dual_topology.nl,
                           square.dual_topology.n_cells_y,
                           square.dual_topology.n_cells_x),
    YAKL_LAMBDA(int k, int j, int i)
    {
      tw1_err.data(0, pks + k, pjs + j, pis + i) =
        abs(tw1.data(0, pks + k, pjs + j, pis + i) - tw1_expected.data(0, pks + k, pjs + j, pis + i));
      tw1_err.data(1, pks + k, pjs + j, pis + i) =
        abs(tw1.data(1, pks + k, pjs + j, pis + i) - tw1_expected.data(1, pks + k, pjs + j, pis + i));
    });


    real errf = yakl::intrinsics::maxval(tw1_err.data);

    return errf;
}


void test_H_convergence()
{
  const int nlevels = 5;
  const real atol = 0.1;

  {
    const int diff_order = 2;
    auto conv_x = ConvergenceTest<nlevels>("H 2 x", compute_H_error<diff_order>, vecfun_x);
    conv_x.check_rate(3, atol);
    auto conv_y = ConvergenceTest<nlevels>("H 2 y", compute_H_error<diff_order>, vecfun_y);
    conv_y.check_rate(3, atol);
    auto conv_xy = ConvergenceTest<nlevels>("H 2 xy", compute_H_error<diff_order>, vecfun_xy);
    conv_xy.check_rate(3, atol);
  }

  {
    const int diff_order = 4;
    auto conv_x = ConvergenceTest<nlevels>("H 4 x", compute_H_error<diff_order>, vecfun_x);
    conv_x.check_rate(5, atol);
    auto conv_y = ConvergenceTest<nlevels>("H 4 y", compute_H_error<diff_order>, vecfun_y);
    conv_y.check_rate(5, atol);
    auto conv_xy = ConvergenceTest<nlevels>("H 4 xy", compute_H_error<diff_order>, vecfun_xy);
    conv_xy.check_rate(3, atol);
  }

  {
    const int diff_order = 6;
    auto conv_x = ConvergenceTest<nlevels>("H 6 x", compute_H_error<diff_order>, vecfun_x);
    conv_x.check_rate(7, atol);
    auto conv_y = ConvergenceTest<nlevels>("H 6 y", compute_H_error<diff_order>, vecfun_y);
    conv_y.check_rate(7, atol);
    auto conv_xy = ConvergenceTest<nlevels>("H 6 xy", compute_H_error<diff_order>, vecfun_xy);
    conv_xy.check_rate(3, atol);
  }
}


int main()
{
  yakl::init();

  test_I_convergence();
  test_H_convergence();

  yakl::finalize();
}
