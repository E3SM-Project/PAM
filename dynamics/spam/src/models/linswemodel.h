#pragma once

#include "common.h"
#include "model.h"
#include "stats.h"

#include "ext_deriv.h"
#include "fct.h"
#include "hamiltonian.h"
#include "hodge_star.h"
#include "recon.h"
#include "wedge.h"

real constexpr _H = 10;
real constexpr _G = 10;

Hamiltonian_Hk Hk;
Hamiltonian_SWE_Hs Hs;
ThermoPotential thermo;

real YAKL_INLINE sol2d_h(real x, real y, real t) {
  real xx = x + 50;
  real yy = y + 50;
  real kx = 2 * M_PI / 100;
  real ky = 2 * M_PI / 100;
  real const c = std::sqrt(_H * _G);
  real alpha = c * std::sqrt(kx * kx + ky * ky);
  return sin(kx * xx) * sin(ky * yy) * cos(alpha * t);
}
vec<2> YAKL_INLINE sol2d_v(real x, real y, real t) {
  vec<2> vvec;

  real xx = x + 50;
  real yy = y + 50;
  real kx = 2 * M_PI / 100;
  real ky = 2 * M_PI / 100;
  real const c = std::sqrt(_H * _G);
  real alpha = c * std::sqrt(kx * kx + ky * ky);
  vvec.u =
      -alpha / (2 * _H * kx) * cos(kx * xx) * sin(ky * yy) * sin(alpha * t);
  vvec.v =
      -alpha / (2 * _H * ky) * sin(kx * xx) * cos(ky * yy) * sin(alpha * t);
  return vvec;
}

real YAKL_INLINE gauss(real x) {
  real r = yakl::abs(x) / 10;
  return exp(-2 * r * r);
}
real YAKL_INLINE gaussian1d_h(real x, real y, real t) {
  real const c = std::sqrt(_H * _G);
  return (gauss(x - c * t) + gauss(x + c * t)) / 2;
}

vec<2> YAKL_INLINE gaussian1d_v(real x, real y, real t) {
  real const c = std::sqrt(_H * _G);
  vec<2> vvec;
  vvec.u = c * (gauss(x - c * t) - gauss(x + c * t)) / (2 * _H);
  vvec.v = 0;
  return vvec;
}

// *******   Diagnostics   ***********//

class ModelDiagnostics : public Diagnostics {
public:
  void compute_diag(ModelParameters &params, real time,
                    const VariableSet<nconstant> &const_vars,
                    VariableSet<nprognostic> &x,
                    VariableSet<ndiagnostic> &diagnostic_vars) {

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    if (params.initdataStr == "gaussian1d") {
      auto h_f = [=](real x, real y) { return gaussian1d_h(x, y, time); };
      auto v_f = [=](real x, real y) { return gaussian1d_v(x, y, time); };
      dual_geometry->set_2form_values(
          h_f, diagnostic_vars.fields_arr[HECDIAGVAR], 0);
      primal_geometry->set_0form_values(
          h_f, diagnostic_vars.fields_arr[HE0DIAGVAR], 0);
      primal_geometry->set_1form_values(v_f,
                                        diagnostic_vars.fields_arr[VEDIAGVAR],
                                        0, LINE_INTEGRAL_TYPE::TANGENT);
    } else if (params.initdataStr == "sol2d") {
      auto h_f = [=](real x, real y) { return sol2d_h(x, y, time); };
      auto v_f = [=](real x, real y) { return sol2d_v(x, y, time); };
      dual_geometry->set_2form_values(
          h_f, diagnostic_vars.fields_arr[HECDIAGVAR], 0);
      primal_geometry->set_0form_values(
          h_f, diagnostic_vars.fields_arr[HE0DIAGVAR], 0);
      primal_geometry->set_1form_values(v_f,
                                        diagnostic_vars.fields_arr[VEDIAGVAR],
                                        0, LINE_INTEGRAL_TYPE::TANGENT);
    }

    parallel_for(
        "Compute H0 Diag",
        Bounds<4>(primal_topology->nl, primal_topology->n_cells_y,
                  primal_topology->n_cells_x, primal_topology->nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_I<ndensity, diff_ord>(
              diagnostic_vars.fields_arr[H0DIAGVAR].data,
              x.fields_arr[HVAR].data, *this->primal_geometry,
              *this->dual_geometry, pis, pjs, pks, i, j, k, n);
        });
  }
};

// *******   Tendencies   ***********//

class ModelTendencies : public Tendencies {
public:
  void initialize(ModelParameters &params, const Topology &primal_topo,
                  const Topology &dual_topo, Geometry &primal_geom,
                  Geometry &dual_geom, ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nconstant> &const_exchange) {

    Tendencies::initialize(params, primal_topo, dual_topo, primal_geom,
                           dual_geom, aux_exchange, const_exchange);
    Hk.initialize(*this->primal_geometry, *this->dual_geometry);
    Hs.initialize(thermo, *this->primal_geometry, *this->dual_geometry);
  }

  void compute_constants(VariableSet<nconstant> &const_vars,
                         VariableSet<nprognostic> &x) {}

  void YAKL_INLINE compute_rhs(real dt, VariableSet<nconstant> &const_vars,
                               VariableSet<nprognostic> &x,
                               VariableSet<nauxiliary> &auxiliary_vars,
                               VariableSet<nprognostic> &xtend) {

    const auto G = const_vars.fields_arr[GRAVVAR].data;
    const auto H = const_vars.fields_arr[HREFVAR].data;

    const auto v = x.fields_arr[VVAR].data;
    const auto h = x.fields_arr[HVAR].data;

    auto vtend = xtend.fields_arr[VVAR].data;
    auto htend = xtend.fields_arr[HVAR].data;

    auto u = auxiliary_vars.fields_arr[UVAR].data;
    auto h0 = auxiliary_vars.fields_arr[H0VAR].data;

    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

    parallel_for(
        "Compute h0",
        Bounds<4>(primal_topology->nl, primal_topology->n_cells_y,
                  primal_topology->n_cells_x, primal_topology->nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_I<ndensity, diff_ord>(h0, h, *this->primal_geometry,
                                        *this->dual_geometry, pis, pjs, pks, i,
                                        j, k, n);
        });

    parallel_for(
        "Compute u",
        Bounds<4>(dual_topology->nl, dual_topology->n_cells_y,
                  dual_topology->n_cells_x, dual_topology->nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H<1, diff_ord>(u, v, *this->primal_geometry,
                                 *this->dual_geometry, dis, djs, dks, i, j, k,
                                 n);
        });

    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[H0VAR].exchange_field(
        auxiliary_vars.fields_arr[H0VAR]);

    parallel_for(
        "Compute v tend",
        Bounds<4>(primal_topology->nl, primal_topology->n_cells_y,
                  primal_topology->n_cells_x, primal_topology->nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD1<ndensity>(vtend, G, h0, pis, pjs, pks, i, j, k, n);
        });

    parallel_for(
        "Compute h tend",
        Bounds<4>(dual_topology->nl, dual_topology->n_cells_y,
                  dual_topology->n_cells_x, dual_topology->nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDbar2<ndensity>(htend, H, u, dis, djs, dks, i, j, k, n);
        });
  }
};

// *******   Statistics   ***********//

class ModelStats : public Stats {
public:
  real3d TEarr, KEarr, PEarr;
  void initialize(ModelParameters &params, Parallel &par,
                  const Topology &primal_topo, const Topology &dual_topo,
                  Geometry &primal_geom, Geometry &dual_geom) {
    Stats::initialize(params, par, primal_topo, dual_topo, primal_geom,
                      dual_geom);

    std::cout << "statsize: " << this->statsize << std::endl;

    this->stats_arr[ESTAT].initialize("energy", 3, this->statsize, this->nens,
                                      this->masterproc);
    this->TEarr =
        real3d("TE", this->dual_topology->nl, this->dual_topology->n_cells_y,
               this->dual_topology->n_cells_x);
    this->KEarr =
        real3d("KE", this->dual_topology->nl, this->dual_topology->n_cells_y,
               this->dual_topology->n_cells_x);
    this->PEarr =
        real3d("PE", this->dual_topology->nl, this->dual_topology->n_cells_y,
               this->dual_topology->n_cells_x);
  }

  void compute(VariableSet<nprognostic> &progvars,
               VariableSet<nconstant> &constvars, int tind) {

    for (int n = 0; n < nens; n++) {

      SArray<real, 1, 3> elocal, eglobal;

      elocal(0) = 0.;
      elocal(1) = 0.;
      elocal(2) = 0.;
      eglobal(0) = 0.;
      eglobal(1) = 0.;
      eglobal(2) = 0.;

      int dis = dual_topology->is;
      int djs = dual_topology->js;
      int dks = dual_topology->ks;

      parallel_for(
          "Compute energy stats",
          Bounds<3>(dual_topology->nl, dual_topology->n_cells_y,
                    dual_topology->n_cells_x),
          YAKL_LAMBDA(int k, int j, int i) {
            real KE, PE;
            KE = Hk.compute_KE(progvars.fields_arr[VVAR].data,
                               constvars.fields_arr[HREFCVAR].data, dis, djs, dks,
                               i, j, k, n);
            PE = Hs.compute_PE(progvars.fields_arr[H0VAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, k, n);
            TEarr(k, j, i) = KE + PE;
            KEarr(k, j, i) = KE;
            PEarr(k, j, i) = PE;
          });

      elocal(0) = yakl::intrinsics::sum(TEarr);
      elocal(1) = yakl::intrinsics::sum(KEarr);
      elocal(2) = yakl::intrinsics::sum(PEarr);


    this->ierr = MPI_Ireduce( &elocal, &eglobal, 3, REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[ESTAT]);

    this->ierr = MPI_Waitall(nstats, this->Req, this->Status);


    if (masterproc)
    {

       this->stats_arr[ESTAT].data(0,tind,n) = eglobal(0);
       this->stats_arr[ESTAT].data(1,tind,n) = eglobal(1);
       this->stats_arr[ESTAT].data(2,tind,n) = eglobal(2);

    }
    }
  }
};

// *******   VariableSet Initialization   ***********//
void initialize_variables(
    const Topology &ptopo, const Topology &dtopo,
    SArray<int, 2, nprognostic, 3> &prog_ndofs_arr,
    SArray<int, 2, nconstant, 3> &const_ndofs_arr,
    SArray<int, 2, nauxiliary, 3> &aux_ndofs_arr,
    SArray<int, 2, ndiagnostic, 3> &diag_ndofs_arr,
    std::array<std::string, nprognostic> &prog_names_arr,
    std::array<std::string, nconstant> &const_names_arr,
    std::array<std::string, nauxiliary> &aux_names_arr,
    std::array<std::string, ndiagnostic> &diag_names_arr,
    std::array<const Topology *, nprognostic> &prog_topo_arr,
    std::array<const Topology *, nconstant> &const_topo_arr,
    std::array<const Topology *, nauxiliary> &aux_topo_arr,
    std::array<const Topology *, ndiagnostic> &diag_topo_arr) {

  // primal grid represents straight quantities, dual grid twisted quantities

  // v, h
  prog_topo_arr[VVAR] = &ptopo;
  prog_topo_arr[HVAR] = &dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[HVAR] = "h";
  set_dofs_arr(prog_ndofs_arr, VVAR, 1, 1, 1); // v = straight 1-form
  set_dofs_arr(prog_ndofs_arr, HVAR, ndims, 1,
               ndensity); // h = twisted n-form

  aux_topo_arr[UVAR] = &dtopo;
  aux_names_arr[UVAR] = "u";
  set_dofs_arr(aux_ndofs_arr, UVAR, ndims - 1, 1, 1); // U = twisted (n-1)-form

  aux_topo_arr[H0VAR] = &ptopo;
  aux_names_arr[H0VAR] = "dens0";
  set_dofs_arr(aux_ndofs_arr, H0VAR, 0, 1, ndensity);

  const_topo_arr[GRAVVAR] = &dtopo;
  const_names_arr[GRAVVAR] = "grav";
  set_dofs_arr(const_ndofs_arr, GRAVVAR, ndims - 1, 1, 1);

  const_topo_arr[HREFVAR] = &dtopo;
  const_names_arr[HREFVAR] = "ref height";
  set_dofs_arr(const_ndofs_arr, HREFVAR, ndims - 1, 1, 1);
  
  const_topo_arr[HREFCVAR] = &dtopo;
  const_names_arr[HREFCVAR] = "ref height cell";
  set_dofs_arr(const_ndofs_arr, HREFCVAR, ndims, 1, 1);
  
  const_topo_arr[HSVAR] = &dtopo;
  const_names_arr[HSVAR] = "surf height";
  set_dofs_arr(const_ndofs_arr, HSVAR, ndims, 1, 1);

  diag_topo_arr[H0DIAGVAR] = &ptopo;
  diag_topo_arr[HE0DIAGVAR] = &ptopo;
  diag_topo_arr[HECDIAGVAR] = &dtopo;
  diag_topo_arr[VEDIAGVAR] = &ptopo;

  diag_names_arr[H0DIAGVAR] = "H0 diag";
  diag_names_arr[HE0DIAGVAR] = "HE0 diag";
  diag_names_arr[HECDIAGVAR] = "HEC diag";
  diag_names_arr[VEDIAGVAR] = "VE diag";

  set_dofs_arr(diag_ndofs_arr, H0DIAGVAR, 0, 1, 1);
  set_dofs_arr(diag_ndofs_arr, HE0DIAGVAR, 0, 1, 1);
  set_dofs_arr(diag_ndofs_arr, HECDIAGVAR, ndims, 1, 1);
  set_dofs_arr(diag_ndofs_arr, VEDIAGVAR, 1, 1, 1);
}

//***************** Set Initial Conditions ***************************//

void readModelParamsFile(std::string inFile, ModelParameters &params,
                         Parallel &par, int nz) {
  readParamsFile(inFile, params, par, nz);

  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  // Read the data initialization options
  params.initdataStr = config["initData"].as<std::string>();
}

void set_domain_sizes_ic(ModelParameters &params, std::string initData) {
  params.zlen = 1.0;
  params.zc = 0.5;
  if (initData == "sol2d" || initData == "gaussian1d") {
    params.xlen = 100.0;
    params.ylen = 100.0;
    params.xc = 0.0;
    params.yc = 0.0;
  }
}

void set_initial_conditions(ModelParameters &params,
                            VariableSet<nprognostic> &progvars,
                            VariableSet<nconstant> &constvars,
                            Geometry &primal_geom, Geometry &dual_geom) {

  Hs.set_parameters(_G);
  constvars.fields_arr[HREFVAR].set(_H);
  constvars.fields_arr[HREFCVAR].set(_H);
  constvars.fields_arr[HSVAR].set(0);
  
  constvars.fields_arr[GRAVVAR].set(_G);

  if (params.initdataStr == "gaussian1d") {
    std::cout << "IC: gaussian1d "
              << "\n";
    auto h_f = [=](real x, real y) { return gaussian1d_h(x, y, 0); };
    auto v_f = [=](real x, real y) { return gaussian1d_v(x, y, 0); };
    dual_geom.set_2form_values(h_f, progvars.fields_arr[HVAR], 0);
    primal_geom.set_1form_values(v_f, progvars.fields_arr[VVAR], 0,
                                 LINE_INTEGRAL_TYPE::TANGENT);
  } else if (params.initdataStr == "sol2d") {
    std::cout << "IC: sol2 "
              << "\n";
    auto h_f = [=](real x, real y) { return sol2d_h(x, y, 0); };
    auto v_f = [=](real x, real y) { return sol2d_v(x, y, 0); };
    dual_geom.set_2form_values(h_f, progvars.fields_arr[HVAR], 0);
    primal_geom.set_1form_values(v_f, progvars.fields_arr[VVAR], 0,
                                 LINE_INTEGRAL_TYPE::TANGENT);
  }
}
