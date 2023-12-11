#pragma once

#include "common.h"
#include "model.h"
#include "profiles.h"
#include "refstate.h"
#include "stats.h"

#include "ext_deriv.h"
#include "fct.h"
#include "hamiltonian.h"
#include "hodge_star.h"
#include "recon.h"
#include "thermo.h"
#include "variableset.h"
#include "wedge.h"

namespace pamc {
// *******   Diagnostics   ***********//

struct TotalDensityDiagnostic : public Diagnostic {
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    name = "total_dens";
    topology = dgeom.topology;
    dofs_arr = {1, 1, 1};
    Diagnostic::initialize(pgeom, dgeom, equations);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &dual_topology = dual_geometry.topology;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    const auto &densvar = x.fields_arr[DENSVAR].data;
    const auto &varset = this->equations->varset;

    parallel_for(
        "Compute total density diagnostic",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          const real total_dens =
              varset.get_total_density(densvar, k, j, i, dks, djs, dis, n);

          field.data(0, k + dks, j + djs, i + dis, n) = total_dens;
        });
  }
};

class Dens0Diagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    // concentration 0-forms for dens
    name = "densl";
    topology = pgeom.topology;
    dofs_arr = {
        0, 0,
        VariableSet::ndensity_prognostic}; // densldiag = straight (0,0)-form
    Diagnostic::initialize(pgeom, dgeom, equations);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "Compute DENS0 DIAG",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Hn1bar<VariableSet::ndensity_prognostic, diff_ord,
                         vert_diff_ord>(
              field.data, x.fields_arr[DENSVAR].data, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
        });
  }
};

class QHZDiagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    name = "QHZl";
    topology = dgeom.topology;
    dofs_arr = {ndims - 1, 0, 1}; // // Qldiag = twisted (n-1,0)-form
    Diagnostic::initialize(pgeom, dgeom, equations);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(PVPE, equations->PVPE);
    parallel_for(
        "Compute Q0 DIAG",
        SimpleBounds<4>(dual_topology.ni - 3, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_qhz(field.data, x.fields_arr[VVAR].data,
                           x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data,
                           const_vars.fields_arr[CORIOLISHZVAR].data, dis, djs,
                           dks, i, j, k + 2, n);
        });

    // Bottom is k=1 and top is k=dual_topology.ni-2
    parallel_for(
        "Compute Q0 DIAG TOP/BOTTOM",
        SimpleBounds<3>(dual_topology.n_cells_y, dual_topology.n_cells_x,
                        dual_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          PVPE.compute_qhz_bottom(field.data, x.fields_arr[VVAR].data,
                                  x.fields_arr[WVAR].data,
                                  x.fields_arr[DENSVAR].data,
                                  const_vars.fields_arr[CORIOLISHZVAR].data,
                                  dis, djs, dks, i, j, 1, n);
          PVPE.compute_qhz_top(field.data, x.fields_arr[VVAR].data,
                               x.fields_arr[WVAR].data,
                               x.fields_arr[DENSVAR].data,
                               const_vars.fields_arr[CORIOLISHZVAR].data, dis,
                               djs, dks, i, j, dual_topology.ni - 2, n);
        });

    field.set_bnd(0.0);
  }
};

class QXYDiagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    name = "QXYl";
    topology = dgeom.topology;
    dofs_arr = {0, 1, 1}; // // Qxydiag = twisted (0,1)-form
    Diagnostic::initialize(pgeom, dgeom, equations);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(PVPE, equations->PVPE);
    parallel_for(
        "Compute QXYl DIAG",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_qxy(field.data, x.fields_arr[VVAR].data,
                           x.fields_arr[DENSVAR].data,
                           const_vars.fields_arr[CORIOLISXYVAR].data, dis, djs,
                           dks, i, j, k, n);
          real dz = dual_geometry.get_dz(k + dks, n);
          field.data(0, k + dks, j + djs, i + dis, n) *= dz;
        });
  }
};

void add_model_diagnostics(
    std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
  diagnostics.emplace_back(std::make_unique<Dens0Diagnostic>());
  diagnostics.emplace_back(std::make_unique<QHZDiagnostic>());
  if (ndims > 1) {
    diagnostics.emplace_back(std::make_unique<QXYDiagnostic>());
  }
  diagnostics.emplace_back(std::make_unique<TotalDensityDiagnostic>());
}

struct TotalDensityFunctor {
  const VariableSet &varset;
  YAKL_INLINE real operator()(const real5d &densvar, int d, int k, int j, int i,
                              int n) const {
    return varset.get_total_density(densvar, k, j, i, 0, 0, 0, n);
  };
  YAKL_INLINE real operator()(const real3d &densvar, int d, int k,
                              int n) const {
    return varset.get_total_density(densvar, k, 0, n);
  };
};

// *******   Tendencies   ***********//

class ModelTendencies : public ExtrudedTendencies {
  real scalar_horiz_diffusion_coeff;
  real scalar_vert_diffusion_coeff;
  real velocity_vort_horiz_diffusion_coeff;
  real velocity_vort_vert_diffusion_coeff;
  real velocity_div_horiz_diffusion_coeff;
  real velocity_div_vert_diffusion_coeff;
  real scalar_diffusion_subtract_refstate;
  real velocity_diffusion_subtract_refstate;
  bool force_refstate_hydrostatic_balance;
  bool check_anelastic_constraint;
  using VS = VariableSet;

public:
  void initialize(ModelParameters &params, Equations &equations,
                  LinearSystem &linear_system,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {

    ExtrudedTendencies::initialize(params, equations, linear_system,
                                   primal_geom, dual_geom);
    scalar_horiz_diffusion_coeff = params.scalar_horiz_diffusion_coeff;
    scalar_vert_diffusion_coeff = params.scalar_vert_diffusion_coeff;
    velocity_vort_horiz_diffusion_coeff =
        params.velocity_vort_horiz_diffusion_coeff;
    velocity_vort_vert_diffusion_coeff =
        params.velocity_vort_vert_diffusion_coeff;
    velocity_div_horiz_diffusion_coeff =
        params.velocity_div_horiz_diffusion_coeff;
    velocity_div_vert_diffusion_coeff =
        params.velocity_div_vert_diffusion_coeff;

    scalar_diffusion_subtract_refstate =
        params.scalar_diffusion_subtract_refstate;
    velocity_diffusion_subtract_refstate =
        params.velocity_diffusion_subtract_refstate;
    force_refstate_hydrostatic_balance =
        params.force_refstate_hydrostatic_balance;
    check_anelastic_constraint = params.check_anelastic_constraint;
  }

  void convert_dynamics_to_coupler_state(PamCoupler &coupler,
                                         const FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars,
                                         bool couple_wind,
                                         bool couple_wind_exact_inverse) {
    convert_dynamics_to_coupler_densities(equations->varset, coupler, prog_vars,
                                          const_vars);

    if (couple_wind) {
      convert_dynamics_to_coupler_wind(equations->varset, coupler, prog_vars,
                                       const_vars, couple_wind_exact_inverse);
    }
  }

  void convert_dynamics_to_coupler_state_staggered(
      PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars,
      const FieldSet<nconstant> &const_vars) {
    convert_dynamics_to_coupler_densities(equations->varset, coupler, prog_vars,
                                          const_vars);
    convert_dynamics_to_coupler_staggered_wind(equations->varset, coupler,
                                               prog_vars, const_vars);
  }

  void convert_coupler_to_dynamics_state(PamCoupler &coupler,
                                         FieldSet<nprognostic> &prog_vars,
                                         FieldSet<nauxiliary> &auxiliary_vars,
                                         FieldSet<nconstant> &const_vars,
                                         bool couple_wind,
                                         bool couple_wind_exact_inverse) {
    convert_coupler_to_dynamics_densities(equations->varset, coupler, prog_vars,
                                          const_vars);
    if (couple_wind) {
      convert_coupler_to_dynamics_wind(equations->varset, coupler, prog_vars,
                                       const_vars, couple_wind_exact_inverse);
    }
    prog_vars.exchange();
    const_vars.exchange();
    auxiliary_vars.exchange();

#if defined PAMC_AN || defined PAMC_MAN
    if (couple_wind) {
      project_to_anelastic(const_vars, prog_vars, auxiliary_vars);
    }
#endif
  }

  void convert_coupler_to_dynamics_state_staggered(
      PamCoupler &coupler, FieldSet<nprognostic> &prog_vars,
      FieldSet<nauxiliary> &auxiliary_vars, FieldSet<nconstant> &const_vars) {
    convert_coupler_to_dynamics_densities(equations->varset, coupler, prog_vars,
                                          const_vars);
    convert_coupler_to_dynamics_staggered_wind(equations->varset, coupler,
                                               prog_vars, const_vars);
#if defined(PAMC_AN) || defined(PAMC_MAN)
    project_to_anelastic(const_vars, prog_vars, auxiliary_vars);
#endif
  }

  void compute_constants(FieldSet<nconstant> &const_vars,
                         FieldSet<nprognostic> &x) override {}

  void compute_dens0(real5d dens0var, const real5d densvar) {

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    const auto &refstate = this->equations->reference_state;
    const auto &varset = this->equations->varset;

    YAKL_SCOPE(refdens, refstate.dens.data);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);

    const auto subtract_refstate_f =
        YAKL_LAMBDA(const real5d &densvar, int d, int k, int j, int i, int n) {
      return densvar(d, k, j, i, n) - refdens(d, k, n);
    };

    parallel_for(
        "Compute Dens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
#if defined PAMC_AN || defined PAMC_MAN
          for (int d = 0; d < VS::ndensity_prognostic; ++d) {
            dens0var(d, k + pks, j + pjs, i + pis, n) =
                (densvar(d, k + pks, j + pjs, i + pis, n) -
                 refdens(d, k + pks, n)) /
                varset.get_total_density(densvar, k, j, i, pks, pjs, pis, n);
          }
#else
          compute_Hn1bar<VS::ndensity_prognostic, diff_ord, vert_diff_ord>(
              subtract_refstate_f, dens0var, densvar, primal_geometry,
              dual_geometry, pis, pjs, pks, i, j, k, n);
#endif
        });
  }

  void compute_U(real5d Uvar, const real5d Vvar) {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    parallel_for(
        "Compute Uvar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H10<1, diff_ord>(Uvar, Vvar, primal_geometry, dual_geometry,
                                   dis, djs, dks, i, j, k, n);
        });
  }

  void compute_UW(real5d UWvar, const real5d Wvar) {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    parallel_for(
        "Compute UWVar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H01<1, vert_diff_ord>(UWvar, Wvar, primal_geometry,
                                        dual_geometry, dis, djs, dks, i, j,
                                        k + 1, n);
        });
  }

  void compute_F_and_FW(real5d Fvar, real5d FWvar, real5d densvar, real5d Vvar,
                        real5d Wvar) {

    const auto &dual_topology = dual_geometry.topology;
    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(Hk, this->equations->Hk);
    parallel_for(
        "Compute F and Fw",
        SimpleBounds<4>(dual_topology.ni, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> he, u;
          real hew, uw, K2;
          Hk.compute_he_U_and_K(he, hew, u, uw, K2, densvar, Vvar, Wvar, pis,
                                pjs, pks, i, j, k, n);

          for (int d = 0; d < ndims; ++d) {
            Fvar(d, pks + k, pjs + j, pis + i, n) = he(d) * u(d);
          }
          FWvar(0, pks + k, pjs + j, pis + i, n) = hew * uw;
        });
  }

  void compute_FT_and_FTW(real5d FTvar, real5d FTWvar, real5d Fvar,
                          real5d FWvar, optional_real5d opt_FTxyvar) {

    const auto &primal_topology = primal_geometry.topology;
    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "Compute FTvar",
        SimpleBounds<4>(primal_topology.ni - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Wxz_u(FTvar, FWvar, pis, pjs, pks, i, j, k + 1, n);
        });
    parallel_for(
        "Compute FTvar bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          compute_Wxz_u_bottom(FTvar, FWvar, pis, pjs, pks, i, j, 0, n);
          compute_Wxz_u_top(FTvar, FWvar, pis, pjs, pks, i, j,
                            primal_topology.ni - 1, n);
        });
    parallel_for(
        "Compute FTWvar",
        SimpleBounds<4>(primal_topology.nl - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Wxz_w(FTWvar, Fvar, pis, pjs, pks, i, j, k + 1, n);
        });
    parallel_for(
        "Compute FTWvar bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          compute_Wxz_w_bottom(FTWvar, Fvar, pis, pjs, pks, i, j, 0, n);
          compute_Wxz_w_top(FTWvar, Fvar, pis, pjs, pks, i, j,
                            primal_topology.nl - 1, n);
        });

    if (ndims > 1) {
      auto FTxyvar = opt_FTxyvar.value();
      parallel_for(
          "Compute FTXY",
          SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_W(FTxyvar, Fvar, pis, pjs, pks, i, j, k, n);
          });
    }
  }

  void compute_q_and_f(real5d qhzvar, real5d fhzvar, const real5d Vvar,
                       const real5d Wvar, const real5d densvar,
                       const real5d coriolishzvar, optional_real5d opt_qxyvar,
                       optional_real5d opt_fxyvar,
                       optional_real5d opt_coriolisxyvar) {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(PVPE, this->equations->PVPE);
    parallel_for(
        "Compute Qhz, Fhz",
        SimpleBounds<4>(dual_topology.ni - 4, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_qhzfhz(qhzvar, fhzvar, Vvar, Wvar, densvar,
                              coriolishzvar, dis, djs, dks, i, j, k + 2, n);
        });
    parallel_for(
        "Compute Qhz, Fhz bnd",
        SimpleBounds<3>(dual_topology.n_cells_y, dual_topology.n_cells_x,
                        dual_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          PVPE.compute_qhzfhz_bottom(qhzvar, fhzvar, Vvar, Wvar, densvar,
                                     coriolishzvar, dis, djs, dks, i, j, 1, n);
          PVPE.compute_qhzfhz_top(qhzvar, fhzvar, Vvar, Wvar, densvar,
                                  coriolishzvar, dis, djs, dks, i, j,
                                  dual_topology.ni - 2, n);
        });

    if (ndims > 1) {
      auto qxyvar = opt_qxyvar.value();
      auto fxyvar = opt_fxyvar.value();
      auto coriolisxyvar = opt_coriolisxyvar.value();
      parallel_for(
          "Compute Qxy, Fxy",
          SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                          dual_topology.n_cells_x, dual_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            PVPE.compute_qxyfxy(qxyvar, fxyvar, Vvar, Wvar, densvar,
                                coriolisxyvar, dis, djs, dks, i, j, k, n);
          });
    }
  }

  void compute_edge_reconstructions_uniform(
      real5d densedgereconvar, real5d densvertedgereconvar,
      real5d qhzedgereconvar, real5d qhzvertedgereconvar,
      real5d coriolishzedgereconvar, real5d coriolishzvertedgereconvar,
      const real5d dens0var, const real5d qhzvar, const real5d fhzvar,
      optional_real5d opt_qxyedgereconvar, optional_real5d opt_qxyvar,
      optional_real5d opt_coriolisxyedgereconvar, optional_real5d opt_fxyvar) {
    yakl::timer_start("compute_edge_reconstructions_uniform");

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(dual_wenoRecon, this->dual_wenoRecon);
    YAKL_SCOPE(dual_to_gll, this->dual_to_gll);
    YAKL_SCOPE(dual_wenoIdl, this->dual_wenoIdl);
    YAKL_SCOPE(dual_wenoSigma, this->dual_wenoSigma);

    YAKL_SCOPE(dual_vert_wenoRecon, this->dual_vert_wenoRecon);
    YAKL_SCOPE(dual_vert_to_gll, this->dual_vert_to_gll);
    YAKL_SCOPE(dual_vert_wenoIdl, this->dual_vert_wenoIdl);
    YAKL_SCOPE(dual_vert_wenoSigma, this->dual_vert_wenoSigma);

    parallel_for(
        "ComputeDensEdgeRecon",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_edge_recon<VS::ndensity_prognostic,
                                     dual_reconstruction_type,
                                     dual_reconstruction_order>(
              densedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
          compute_twisted_vert_edge_recon_uniform<
              VS::ndensity_prognostic, dual_vert_reconstruction_type,
              dual_vert_reconstruction_order>(
              densvertedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl,
              dual_vert_wenoSigma);
        });

    YAKL_SCOPE(primal_wenoRecon, this->primal_wenoRecon);
    YAKL_SCOPE(primal_to_gll, this->primal_to_gll);
    YAKL_SCOPE(primal_wenoIdl, this->primal_wenoIdl);
    YAKL_SCOPE(primal_wenoSigma, this->primal_wenoSigma);

    YAKL_SCOPE(primal_vert_wenoRecon, this->primal_vert_wenoRecon);
    YAKL_SCOPE(primal_vert_to_gll, this->primal_vert_to_gll);
    YAKL_SCOPE(primal_vert_wenoIdl, this->primal_vert_wenoIdl);
    YAKL_SCOPE(primal_vert_wenoSigma, this->primal_vert_wenoSigma);

    YAKL_SCOPE(coriolis_wenoRecon, this->coriolis_wenoRecon);
    YAKL_SCOPE(coriolis_to_gll, this->coriolis_to_gll);
    YAKL_SCOPE(coriolis_wenoIdl, this->coriolis_wenoIdl);
    YAKL_SCOPE(coriolis_wenoSigma, this->coriolis_wenoSigma);

    YAKL_SCOPE(coriolis_vert_wenoRecon, this->coriolis_vert_wenoRecon);
    YAKL_SCOPE(coriolis_vert_to_gll, this->coriolis_vert_to_gll);
    YAKL_SCOPE(coriolis_vert_wenoIdl, this->coriolis_vert_wenoIdl);
    YAKL_SCOPE(coriolis_vert_wenoSigma, this->coriolis_vert_wenoSigma);

    parallel_for(
        "ComputeQhzEdgeRecon",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_straight_hz_edge_recon<1, reconstruction_type,
                                         reconstruction_order>(
              qhzedgereconvar, qhzvar, pis, pjs, pks, i, j, k, n,
              primal_wenoRecon, primal_to_gll, primal_wenoIdl,
              primal_wenoSigma);
          compute_straight_hz_vert_edge_recon_uniform<
              ndims, vert_reconstruction_type, vert_reconstruction_order>(
              qhzvertedgereconvar, qhzvar, pis, pjs, pks, i, j, k, n,
              primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl,
              primal_vert_wenoSigma);
          compute_straight_hz_edge_recon<1, coriolis_reconstruction_type,
                                         coriolis_reconstruction_order>(
              coriolishzedgereconvar, fhzvar, pis, pjs, pks, i, j, k, n,
              coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
              coriolis_wenoSigma);
          compute_straight_hz_vert_edge_recon_uniform<
              ndims, coriolis_vert_reconstruction_type,
              coriolis_vert_reconstruction_order>(
              coriolishzvertedgereconvar, fhzvar, pis, pjs, pks, i, j, k, n,
              coriolis_vert_wenoRecon, coriolis_vert_to_gll,
              coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
        });

    if (ndims > 1) {
      auto qxyvar = opt_qxyvar.value();
      auto qxyedgereconvar = opt_qxyedgereconvar.value();
      auto fxyvar = opt_fxyvar.value();
      auto coriolisxyedgereconvar = opt_coriolisxyedgereconvar.value();

      parallel_for(
          "ComputeQxyEdgeRecon",
          SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_straight_edge_recon<1, reconstruction_type,
                                        reconstruction_order>(
                qxyedgereconvar, qxyvar, pis, pjs, pks, i, j, k, n,
                primal_wenoRecon, primal_to_gll, primal_wenoIdl,
                primal_wenoSigma);
            compute_straight_edge_recon<1, coriolis_reconstruction_type,
                                        coriolis_reconstruction_order>(
                coriolisxyedgereconvar, fxyvar, pis, pjs, pks, i, j, k, n,
                coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
                coriolis_wenoSigma);
          });
    }
    yakl::timer_stop("compute_edge_reconstructions_uniform");
  }

  void compute_edge_reconstructions_variable(
      real5d densedgereconvar, real5d densvertedgereconvar,
      real5d qhzedgereconvar, real5d qhzvertedgereconvar,
      real5d coriolishzedgereconvar, real5d coriolishzvertedgereconvar,
      const real5d dens0var, const real5d qhzvar, const real5d fhzvar,
      optional_real5d opt_qxyedgereconvar, optional_real5d opt_qxyvar,
      optional_real5d opt_coriolisxyedgereconvar, optional_real5d opt_fxyvar) {
    yakl::timer_start("compute_edge_reconstructions_variable");

    if (vert_reconstruction_type == RECONSTRUCTION_TYPE::WENO ||
        dual_vert_reconstruction_type == RECONSTRUCTION_TYPE::WENO ||
        coriolis_vert_reconstruction_type == RECONSTRUCTION_TYPE::WENO) {
      throw std::runtime_error("classical WENO not supported with vertically "
                               "variable grid, you may want WENOFUNC");
    }

    if ((vert_reconstruction_type == RECONSTRUCTION_TYPE::CFV &&
         vert_reconstruction_order > 1) ||
        (dual_vert_reconstruction_type == RECONSTRUCTION_TYPE::CFV &&
         dual_vert_reconstruction_order > 1) ||
        (coriolis_vert_reconstruction_type == RECONSTRUCTION_TYPE::CFV &&
         coriolis_vert_reconstruction_order > 1)) {
      throw std::runtime_error("CFV reconstruction only supports order 1 with "
                               "vertically variable grid");
    }

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(dual_wenoRecon, this->dual_wenoRecon);
    YAKL_SCOPE(dual_to_gll, this->dual_to_gll);
    YAKL_SCOPE(dual_wenoIdl, this->dual_wenoIdl);
    YAKL_SCOPE(dual_wenoSigma, this->dual_wenoSigma);

    YAKL_SCOPE(dual_vert_coefs_to_gll_arr, this->dual_vert_coefs_to_gll_arr);
    YAKL_SCOPE(dual_vert_sten_to_gll_arr, this->dual_vert_sten_to_gll_arr);
    YAKL_SCOPE(dual_vert_sten_to_coefs_arr, this->dual_vert_sten_to_coefs_arr);
    YAKL_SCOPE(dual_vert_weno_recon_lower_arr,
               this->dual_vert_weno_recon_lower_arr);
    YAKL_SCOPE(dual_vert_wenoIdl, this->dual_vert_wenoIdl);
    YAKL_SCOPE(dual_vert_wenoSigma, this->dual_vert_wenoSigma);

    parallel_for(
        "ComputeDensEdgeRecon",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, dual_vert_reconstruction_order, 2>
              dual_vert_coefs_to_gll;
          SArray<real, 2, dual_vert_reconstruction_order, 2>
              dual_vert_sten_to_gll;
          SArray<real, 2, dual_vert_reconstruction_order,
                 dual_vert_reconstruction_order>
              dual_vert_sten_to_coefs;
          SArray<real, 3, (dual_vert_reconstruction_order - 1) / 2 + 1,
                 (dual_vert_reconstruction_order - 1) / 2 + 1,
                 (dual_vert_reconstruction_order - 1) / 2 + 1>
              dual_vert_weno_recon_lower;

          for (int h = 0; h < dual_vert_reconstruction_order; h++) {
            for (int g = 0; g < 2; g++) {
              dual_vert_coefs_to_gll(h, g) =
                  dual_vert_coefs_to_gll_arr(k, h, g, n);
              dual_vert_sten_to_gll(h, g) =
                  dual_vert_sten_to_gll_arr(k, h, g, n);
            }
          }

          for (int h1 = 0; h1 < dual_vert_reconstruction_order; h1++) {
            for (int h2 = 0; h2 < dual_vert_reconstruction_order; h2++) {
              dual_vert_sten_to_coefs(h1, h2) =
                  dual_vert_sten_to_coefs_arr(k, h1, h2, n);
            }
          }

          for (int h1 = 0; h1 < (dual_vert_reconstruction_order - 1) / 2 + 1;
               h1++) {
            for (int h2 = 0; h2 < (dual_vert_reconstruction_order - 1) / 2 + 1;
                 h2++) {
              for (int h3 = 0;
                   h3 < (dual_vert_reconstruction_order - 1) / 2 + 1; h3++) {
                dual_vert_weno_recon_lower(h1, h2, h3) =
                    dual_vert_weno_recon_lower_arr(k, h1, h2, h3, n);
              }
            }
          }

          compute_twisted_edge_recon<VS::ndensity_prognostic,
                                     dual_reconstruction_type,
                                     dual_reconstruction_order>(
              densedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);

          compute_twisted_vert_edge_recon_variable<
              VS::ndensity_prognostic, dual_vert_reconstruction_type,
              dual_vert_reconstruction_order>(
              densvertedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_vert_coefs_to_gll, dual_vert_sten_to_gll,
              dual_vert_sten_to_coefs, dual_vert_weno_recon_lower,
              dual_vert_wenoIdl, dual_vert_wenoSigma);
        });

    YAKL_SCOPE(primal_wenoRecon, this->primal_wenoRecon);
    YAKL_SCOPE(primal_to_gll, this->primal_to_gll);
    YAKL_SCOPE(primal_wenoIdl, this->primal_wenoIdl);
    YAKL_SCOPE(primal_wenoSigma, this->primal_wenoSigma);

    YAKL_SCOPE(primal_vert_coefs_to_gll_arr,
               this->primal_vert_coefs_to_gll_arr);
    YAKL_SCOPE(primal_vert_sten_to_gll_arr, this->primal_vert_sten_to_gll_arr);
    YAKL_SCOPE(primal_vert_sten_to_coefs_arr,
               this->primal_vert_sten_to_coefs_arr);
    YAKL_SCOPE(primal_vert_weno_recon_lower_arr,
               this->primal_vert_weno_recon_lower_arr);
    YAKL_SCOPE(primal_vert_wenoIdl, this->primal_vert_wenoIdl);
    YAKL_SCOPE(primal_vert_wenoSigma, this->primal_vert_wenoSigma);

    YAKL_SCOPE(coriolis_wenoRecon, this->coriolis_wenoRecon);
    YAKL_SCOPE(coriolis_to_gll, this->coriolis_to_gll);
    YAKL_SCOPE(coriolis_wenoIdl, this->coriolis_wenoIdl);
    YAKL_SCOPE(coriolis_wenoSigma, this->coriolis_wenoSigma);

    YAKL_SCOPE(coriolis_vert_coefs_to_gll_arr,
               this->coriolis_vert_coefs_to_gll_arr);
    YAKL_SCOPE(coriolis_vert_sten_to_gll_arr,
               this->coriolis_vert_sten_to_gll_arr);
    YAKL_SCOPE(coriolis_vert_sten_to_coefs_arr,
               this->coriolis_vert_sten_to_coefs_arr);
    YAKL_SCOPE(coriolis_vert_weno_recon_lower_arr,
               this->coriolis_vert_weno_recon_lower_arr);
    YAKL_SCOPE(coriolis_vert_wenoIdl, this->coriolis_vert_wenoIdl);
    YAKL_SCOPE(coriolis_vert_wenoSigma, this->coriolis_vert_wenoSigma);

    parallel_for(
        "ComputeQEdgeRecon",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, vert_reconstruction_order, 2> primal_vert_to_gll;
          SArray<real, 3, vert_reconstruction_order, vert_reconstruction_order,
                 vert_reconstruction_order>
              primal_vert_wenoRecon;

          SArray<real, 2, coriolis_vert_reconstruction_order, 2>
              coriolis_vert_to_gll;
          SArray<real, 3, coriolis_vert_reconstruction_order,
                 coriolis_vert_reconstruction_order,
                 coriolis_vert_reconstruction_order>
              coriolis_vert_wenoRecon;

          SArray<real, 2, coriolis_vert_reconstruction_order, 2>
              coriolis_vert_coefs_to_gll;
          SArray<real, 2, coriolis_vert_reconstruction_order, 2>
              coriolis_vert_sten_to_gll;
          SArray<real, 2, coriolis_vert_reconstruction_order,
                 coriolis_vert_reconstruction_order>
              coriolis_vert_sten_to_coefs;
          SArray<real, 3, (coriolis_vert_reconstruction_order - 1) / 2 + 1,
                 (coriolis_vert_reconstruction_order - 1) / 2 + 1,
                 (coriolis_vert_reconstruction_order - 1) / 2 + 1>
              coriolis_vert_weno_recon_lower;

          SArray<real, 2, vert_reconstruction_order, 2>
              primal_vert_coefs_to_gll;
          SArray<real, 2, vert_reconstruction_order, 2> primal_vert_sten_to_gll;
          SArray<real, 2, vert_reconstruction_order, vert_reconstruction_order>
              primal_vert_sten_to_coefs;
          SArray<real, 3, (vert_reconstruction_order - 1) / 2 + 1,
                 (vert_reconstruction_order - 1) / 2 + 1,
                 (vert_reconstruction_order - 1) / 2 + 1>
              primal_vert_weno_recon_lower;

          for (int h = 0; h < vert_reconstruction_order; h++) {
            for (int g = 0; g < 2; g++) {
              primal_vert_coefs_to_gll(h, g) =
                  primal_vert_coefs_to_gll_arr(k, h, g, n);
              primal_vert_sten_to_gll(h, g) =
                  primal_vert_sten_to_gll_arr(k, h, g, n);
            }
          }

          for (int h1 = 0; h1 < vert_reconstruction_order; h1++) {
            for (int h2 = 0; h2 < vert_reconstruction_order; h2++) {
              primal_vert_sten_to_coefs(h1, h2) =
                  primal_vert_sten_to_coefs_arr(k, h1, h2, n);
            }
          }

          for (int h1 = 0; h1 < (vert_reconstruction_order - 1) / 2 + 1; h1++) {
            for (int h2 = 0; h2 < (vert_reconstruction_order - 1) / 2 + 1;
                 h2++) {
              for (int h3 = 0; h3 < (vert_reconstruction_order - 1) / 2 + 1;
                   h3++) {
                primal_vert_weno_recon_lower(h1, h2, h3) =
                    primal_vert_weno_recon_lower_arr(k, h1, h2, h3, n);
              }
            }
          }

          for (int h = 0; h < coriolis_vert_reconstruction_order; h++) {
            for (int g = 0; g < 2; g++) {
              coriolis_vert_coefs_to_gll(h, g) =
                  coriolis_vert_coefs_to_gll_arr(k, h, g, n);
              coriolis_vert_sten_to_gll(h, g) =
                  coriolis_vert_sten_to_gll_arr(k, h, g, n);
            }
          }

          for (int h1 = 0; h1 < coriolis_vert_reconstruction_order; h1++) {
            for (int h2 = 0; h2 < coriolis_vert_reconstruction_order; h2++) {
              coriolis_vert_sten_to_coefs(h1, h2) =
                  coriolis_vert_sten_to_coefs_arr(k, h1, h2, n);
            }
          }

          for (int h1 = 0;
               h1 < (coriolis_vert_reconstruction_order - 1) / 2 + 1; h1++) {
            for (int h2 = 0;
                 h2 < (coriolis_vert_reconstruction_order - 1) / 2 + 1; h2++) {
              for (int h3 = 0;
                   h3 < (coriolis_vert_reconstruction_order - 1) / 2 + 1;
                   h3++) {
                coriolis_vert_weno_recon_lower(h1, h2, h3) =
                    coriolis_vert_weno_recon_lower_arr(k, h1, h2, h3, n);
              }
            }
          }

          compute_straight_hz_edge_recon<1, reconstruction_type,
                                         reconstruction_order>(
              qhzedgereconvar, qhzvar, pis, pjs, pks, i, j, k, n,
              primal_wenoRecon, primal_to_gll, primal_wenoIdl,
              primal_wenoSigma);
          compute_straight_hz_vert_edge_recon_variable<
              ndims, vert_reconstruction_type, vert_reconstruction_order>(
              qhzvertedgereconvar, qhzvar, pis, pjs, pks, i, j, k, n,
              primal_vert_coefs_to_gll, primal_vert_sten_to_gll,
              primal_vert_sten_to_coefs, primal_vert_weno_recon_lower,
              primal_vert_wenoIdl, primal_vert_wenoSigma);
          compute_straight_hz_edge_recon<1, coriolis_reconstruction_type,
                                         coriolis_reconstruction_order>(
              coriolishzedgereconvar, fhzvar, pis, pjs, pks, i, j, k, n,
              coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
              coriolis_wenoSigma);
          compute_straight_hz_vert_edge_recon_variable<
              ndims, coriolis_vert_reconstruction_type,
              coriolis_vert_reconstruction_order>(
              coriolishzvertedgereconvar, fhzvar, pis, pjs, pks, i, j, k, n,
              coriolis_vert_coefs_to_gll, coriolis_vert_sten_to_gll,
              coriolis_vert_sten_to_coefs, coriolis_vert_weno_recon_lower,
              coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
        });

    if (ndims > 1) {
      auto qxyvar = opt_qxyvar.value();
      auto qxyedgereconvar = opt_qxyedgereconvar.value();
      auto fxyvar = opt_fxyvar.value();
      auto coriolisxyedgereconvar = opt_coriolisxyedgereconvar.value();

      parallel_for(
          "ComputeQxyEdgeRecon",
          SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_straight_edge_recon<1, reconstruction_type,
                                        reconstruction_order>(
                qxyedgereconvar, qxyvar, pis, pjs, pks, i, j, k, n,
                primal_wenoRecon, primal_to_gll, primal_wenoIdl,
                primal_wenoSigma);
            compute_straight_edge_recon<1, coriolis_reconstruction_type,
                                        coriolis_reconstruction_order>(
                coriolisxyedgereconvar, fxyvar, pis, pjs, pks, i, j, k, n,
                coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
                coriolis_wenoSigma);
          });
    }
    yakl::timer_stop("compute_edge_reconstructions_variable");
  }

  void compute_recons(
      real5d densreconvar, real5d densvertreconvar, real5d qhzreconvar,
      real5d qhzvertreconvar, real5d coriolishzreconvar,
      real5d coriolishzvertreconvar, const real5d densedgereconvar,
      const real5d densvertedgereconvar, const real5d qhzedgereconvar,
      const real5d qhzvertedgereconvar, const real5d coriolishzedgereconvar,
      const real5d coriolishzvertedgereconvar, const real5d densvar,
      const real5d Fvar, const real5d FWvar, const real5d FTvar,
      const real5d FTWvar, real tanh_upwind_coeff,
      optional_real5d opt_qxyreconvar, optional_real5d opt_qxyedgereconvar,
      optional_real5d opt_coriolisxyreconvar,
      optional_real5d opt_coriolisxyedgereconvar, optional_real5d opt_FTxyvar) {
    yakl::timer_start("compute_recons");

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(varset, this->equations->varset);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);

    const auto &refstate = this->equations->reference_state;

    parallel_for(
        "ComputeDensRECON",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_recon<VS::ndensity_prognostic,
                                dual_reconstruction_type>(
              densreconvar, densedgereconvar, dual_geometry, Fvar,
              tanh_upwind_coeff, dis, djs, dks, i, j, k, n);

#if defined PAMC_AN || defined PAMC_MAN
          // add reference state
          for (int d = 0; d < ndims; d++) {
            for (int l = 0; l < VS::ndensity_prognostic; l++) {
              densreconvar(d + l * ndims, k + dks, j + djs, i + dis, n) +=
                  refstate.q_pi.data(l, k + pks, n);
            }
          }
          for (int d = 0; d < ndims; d++) {
            densreconvar(d + varset.dens_id_mass * ndims, k + dks, j + djs,
                         i + dis, n) = 1;
          }
#else
          const auto total_density_f = TotalDensityFunctor{varset};

          SArray<real, 1, ndims> he;
          SArray<real, 1, 1> dens0_ik;
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_ik, densvar, primal_geometry,
              dual_geometry, pis, pjs, pks, i, j, k, n);

          SArray<real, 1, 1> dens0_im1;
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_im1, densvar, primal_geometry,
              dual_geometry, pis, pjs, pks, i - 1, j, k, n);

          he(0) = 0.5_fp * (dens0_ik(0) + dens0_im1(0));
          if (ndims > 1) {
            SArray<real, 1, 1> dens0_jm1;
            compute_Hn1bar<1, diff_ord, vert_diff_ord>(
                total_density_f, dens0_jm1, densvar, primal_geometry,
                dual_geometry, pis, pjs, pks, i, j - 1, k, n);
            he(1) = 0.5_fp * (dens0_ik(0) + dens0_jm1(0));
          }
          // scale twisted recons and add reference state
          for (int d = 0; d < ndims; d++) {
            for (int l = 0; l < VS::ndensity_prognostic; l++) {
              densreconvar(d + l * ndims, k + dks, j + djs, i + dis, n) +=
                  refstate.rho_pi.data(0, k + pks, n) *
                  refstate.q_pi.data(l, k + pks, n);
              densreconvar(d + l * ndims, k + dks, j + djs, i + dis, n) /=
                  he(d);
            }
          }
#endif
        });

    parallel_for(
        "ComputeDensVertRECON",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_vert_recon<VS::ndensity_prognostic,
                                     dual_vert_reconstruction_type>(
              densvertreconvar, densvertedgereconvar, dual_geometry, FWvar,
              tanh_upwind_coeff, dis, djs, dks, i, j, k + 1, n);

#if defined PAMC_AN || defined PAMC_MAN
          // add reference state
          for (int l = 0; l < VS::ndensity_prognostic; l++) {
            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) +=
                refstate.q_di.data(l, k + dks + 1, n);
          }
          densvertreconvar(varset.dens_id_mass, k + dks + 1, j + djs, i + dis,
                           n) = 1;
#else
          const auto total_density_f = TotalDensityFunctor{varset};
          SArray<real, 1, 1> dens0_kp1;
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_kp1, densvar, primal_geometry,
              dual_geometry, pis, pjs, pks, i, j, k + 1, n);
          SArray<real, 1, 1> dens0_ik;
          compute_Hn1bar<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_ik, densvar, primal_geometry,
              dual_geometry, pis, pjs, pks, i, j, k, n);

          real hew = 0.5_fp * (dens0_kp1(0) + dens0_ik(0));
          // scale twisted recons and add reference state
          for (int l = 0; l < VS::ndensity_prognostic; l++) {
            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) +=
                refstate.rho_di.data(0, k + dks + 1, n) *
                refstate.q_di.data(l, k + dks + 1, n);
            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) /= hew;
          }
#endif
        });

    parallel_for(
        "ComputeQhzRECON",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_straight_hz_recon<1, reconstruction_type>(
              qhzreconvar, qhzedgereconvar, primal_geometry, FTWvar,
              tanh_upwind_coeff, pis, pjs, pks, i, j, k, n);
          compute_straight_hz_recon<1, coriolis_reconstruction_type>(
              coriolishzreconvar, coriolishzedgereconvar, primal_geometry,
              FTWvar, tanh_upwind_coeff, pis, pjs, pks, i, j, k, n);
        });
    parallel_for(
        "ComputeQhzVERTRECON",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_straight_hz_vert_recon<1, vert_reconstruction_type>(
              qhzvertreconvar, qhzvertedgereconvar, primal_geometry, FTvar,
              tanh_upwind_coeff, pis, pjs, pks, i, j, k, n);
          compute_straight_hz_vert_recon<1, coriolis_vert_reconstruction_type>(
              coriolishzvertreconvar, coriolishzvertedgereconvar,
              primal_geometry, FTvar, tanh_upwind_coeff, pis, pjs, pks, i, j, k,
              n);
        });

    if (ndims > 1) {
      auto qxyedgereconvar = opt_qxyedgereconvar.value();
      auto qxyreconvar = opt_qxyreconvar.value();
      auto coriolisxyedgereconvar = opt_coriolisxyedgereconvar.value();
      auto coriolisxyreconvar = opt_coriolisxyreconvar.value();
      auto FTxyvar = opt_FTxyvar.value();
      parallel_for(
          "ComputeQxyRECON",
          SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_straight_recon<1, reconstruction_type>(
                qxyreconvar, qxyedgereconvar, primal_geometry, FTxyvar,
                tanh_upwind_coeff, pis, pjs, pks, i, j, k, n);
            compute_straight_recon<1, coriolis_reconstruction_type>(
                coriolisxyreconvar, coriolisxyedgereconvar, primal_geometry,
                FTxyvar, tanh_upwind_coeff, pis, pjs, pks, i, j, k, n);
          });
    }
    yakl::timer_stop("compute_recons");
  }

  void add_scalar_diffusion(real scalar_horiz_diffusion_coeff,
                            real scalar_vert_diffusion_coeff,
                            real5d denstendvar, const real5d densvar,
                            const real5d dens0var, const real5d Fdiffvar,
                            const real5d FWdiffvar,
                            FieldSet<nauxiliary> &auxiliary_vars) {
    yakl::timer_start("add_scalar_diffusion");

    YAKL_SCOPE(subtract_refstate, this->scalar_diffusion_subtract_refstate);

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(varset, this->equations->varset);
    YAKL_SCOPE(q_di, this->equations->reference_state.q_di.data);
    YAKL_SCOPE(q_pi, this->equations->reference_state.q_pi.data);
    YAKL_SCOPE(refdens, this->equations->reference_state.dens.data);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);

    parallel_for(
        "Scalar diffusion init",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          const real total_dens =
              varset.get_total_density(densvar, k, j, i, dks, djs, dis, n);

          for (int d = 0; d < VS::ndensity_diffused; ++d) {
            int dens_id = varset.diffused_dens_ids(d);
            real dens0 = densvar(dens_id, k + dks, j + djs, i + dis, n);
            dens0 /= total_dens;
            if (subtract_refstate) {
              dens0 -= q_pi(dens_id, k + pks, n);
            }
            dens0var(d, k + pks, j + pjs, i + pis, n) = dens0;
          }
        });
    // Further optimization idea: no need to exchange everything here
    auxiliary_vars.exchange({DENS0VAR});

    parallel_for(
        "Scalar diffusion horz flux",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, VS::ndensity_diffused, ndims> dens_grad;
          compute_D0<VS::ndensity_diffused>(dens_grad, dens0var, pis, pjs, pks,
                                            i, j, k, n);

          SArray<real, 1, ndims> H10_diag;
          H10_diagonal(H10_diag, primal_geometry, dual_geometry, pis, pjs, pks,
                       i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            for (int l = 0; l < VS::ndensity_diffused; ++l) {
              Fdiffvar(ndims * l + d, k + dks, j + djs, i + dis, n) =
                  dens_grad(l, d) * H10_diag(d);
            }
          }
        });
    auxiliary_vars.exchange({FDIFFVAR});

    parallel_for(
        "Scalar diffusion vert flux",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, VS::ndensity_diffused> dens_vert_grad;
          compute_D0_vert<VS::ndensity_diffused>(dens_vert_grad, dens0var, pis,
                                                 pjs, pks, i, j, k, n);
          const real H01_diag = H01_diagonal(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k + 1, n);
          for (int d = 0; d < VS::ndensity_diffused; ++d) {
            FWdiffvar(d, k + 1 + dks, j + djs, i + dis, n) =
                dens_vert_grad(d) * H01_diag;
          }
        });
    auxiliary_vars.exchange({FWDIFFVAR});
    auxiliary_vars.fields_arr[FWDIFFVAR].set_bnd(0.0);

    parallel_for(
        "Scalar diffusion tendency",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, VS::ndensity_diffused> hdiv;
          compute_Dnm1bar<VS::ndensity_diffused>(hdiv, Fdiffvar, dis, djs, dks,
                                                 i, j, k, n);
          SArray<real, 1, VS::ndensity_diffused> vdiv;
          compute_Dnm1bar_vert<VS::ndensity_diffused>(vdiv, FWdiffvar, dis, djs,
                                                      dks, i, j, k, n);

          const real rho =
              varset.get_total_density(densvar, k, j, i, pks, pjs, pis, n);
          const real Hn1bar_diag = Hn1bar_diagonal(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, n);

          for (int d = 0; d < VS::ndensity_diffused; ++d) {
            const real diff_tend =
                -scalar_horiz_diffusion_coeff * rho * Hn1bar_diag * hdiv(d) +
                -scalar_vert_diffusion_coeff * rho * Hn1bar_diag * vdiv(d);
            int dens_id = varset.diffused_dens_ids(d);
            denstendvar(dens_id, k + pks, j + pjs, i + pis, n) += diff_tend;
          }
        });

    yakl::timer_stop("add_scalar_diffusion");
  }

  void add_velocity_diffusion_2d(real velocity_vort_horiz_diffusion_coeff,
                                 real velocity_vort_vert_diffusion_coeff,
                                 real velocity_div_horiz_diffusion_coeff,
                                 real velocity_div_vert_diffusion_coeff,
                                 real5d Vtendvar, real5d Wtendvar,
                                 const real5d Vvar, const real5d Wvar,
                                 const real5d qhzedgereconvar,
                                 const real5d qhzvar, const real5d dens0var,
                                 const real5d Kvar, const real5d Fvar,
                                 const real5d FWvar, const real5d FTWvar,
                                 FieldSet<nauxiliary> &auxiliary_vars) {
    yakl::timer_start("add_velocity_diffusion");

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(refv, this->equations->reference_state.v.data);
    YAKL_SCOPE(subtract_refstate, this->velocity_diffusion_subtract_refstate);

    parallel_for(
        "Velocity diffusion 0",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          for (int d = 0; d < ndims; ++d) {
            FTWvar(d, k + pks, j + pjs, i + pis, n) =
                Vvar(d, k + pks, j + pjs, i + pis, n);
            if (subtract_refstate) {
              FTWvar(d, k + pks, j + pjs, i + pis, n) -= refv(d, k + pks, n);
            }
          }
        });
    auxiliary_vars.exchange({FTWVAR});

    // *d
    parallel_for(
        "Velocity diffusion 1",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> Dv;
          compute_D1_ext<1>(Dv, FTWvar, Wvar, pis, pjs, pks, i, j, k, n);
          const real Hn1_diag = Hn1_diagonal(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k + 1, n);
          qhzvar(0, k + 1 + dks, j + djs, i + dis, n) = Hn1_diag * Dv(0, 0);
        });
    auxiliary_vars.exchange({QHZVAR});
    auxiliary_vars.fields_arr[QHZVAR].set_bnd(0.0);

    // *d*
    parallel_for(
        "Velocity diffusion 2",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H10<1, diffusion_diff_ord>(Fvar, FTWvar, primal_geometry,
                                             dual_geometry, dis, djs, dks, i, j,
                                             k, n);
        });
    auxiliary_vars.exchange({FVAR});

    parallel_for(
        "Velocity diffusion 3",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H01<1, vert_diffusion_diff_ord>(FWvar, Wvar, primal_geometry,
                                                  dual_geometry, dis, djs, dks,
                                                  i, j, k + 1, n);
        });
    auxiliary_vars.exchange({FWVAR});
    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);

    parallel_for(
        "Velocity diffusion 4",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, 1> Dhorz;
          SArray<real, 1, 1> Dvert;
          compute_Dnm1bar<1>(Dhorz, Fvar, dis, djs, dks, i, j, k, n);
          compute_Dnm1bar_vert<1>(Dvert, FWvar, dis, djs, dks, i, j, k, n);
          const real Hn1bar_diag = Hn1bar_diagonal(
              primal_geometry, dual_geometry, dis, djs, dks, i, j, k, n);
          dens0var(0, k + pks, j + pjs, i + pis, n) =
              (Dhorz(0) + Dvert(0)) * Hn1bar_diag;
        });
    auxiliary_vars.exchange({DENS0VAR});

    parallel_for(
        "Velocity diffusion 5",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          // *d(*d) = vort
          const real Dvert =
              compute_D0bar_vert<1>(qhzvar, dis, djs, dks, i, j, k, n);
          SArray<real, 1, ndims> Hnm11bar_diag;
          Hnm11bar_diagonal(Hnm11bar_diag, primal_geometry, dual_geometry, dis,
                            djs, dks, i, j, k, n);
          Vtendvar(0, k + pks, j + pjs, i + pis, n) -=
              velocity_vort_horiz_diffusion_coeff * Dvert * Hnm11bar_diag(0);

          // d(*d*) = div
          SArray<real, 2, 1, ndims> vdiff;
          compute_D0<1>(vdiff, dens0var, pis, pjs, pks, i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            Vtendvar(d, pks + k, pjs + j, pis + i, n) -=
                velocity_div_horiz_diffusion_coeff * vdiff(0, d);
          }
        });

    parallel_for(
        "Velocity diffusion 6",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          // *d(*d) = vort
          SArray<real, 1, ndims> Dhorz;
          compute_D0bar_ext<1>(Dhorz, qhzvar, dis, djs, dks, i, j, k + 1, n);
          const real Hn0bar_diag = Hn0bar_diagonal(
              primal_geometry, dual_geometry, dis, djs, dks, i, j, k, n);
          Wtendvar(0, k + pks, j + pjs, i + pis, n) -=
              velocity_vort_vert_diffusion_coeff * Dhorz(0) * Hn0bar_diag;

          // d(*d*) = div
          SArray<real, 1, 1> wdiff;
          compute_D0_vert<1>(wdiff, dens0var, pis, pjs, pks, i, j, k, n);
          Wtendvar(0, pks + k, pjs + j, pis + i, n) -=
              velocity_div_vert_diffusion_coeff * wdiff(0);
        });

    yakl::timer_stop("add_velocity_diffusion");
  }

  void add_velocity_diffusion_3d(
      real velocity_vort_horiz_diffusion_coeff,
      real velocity_vort_vert_diffusion_coeff,
      real velocity_div_horiz_diffusion_coeff,
      real velocity_div_vert_diffusion_coeff, real5d Vtendvar, real5d Wtendvar,
      const real5d Vvar, const real5d Wvar, const real5d qhzedgereconvar,
      const real5d qhzvar, const real5d dens0var, const real5d Kvar,
      const real5d Fvar, const real5d FWvar, const real5d FTWvar,
      const real5d qxyedgereconvar, const real5d qxyvar,
      FieldSet<nauxiliary> &auxiliary_vars) {
    yakl::timer_start("add_velocity_diffusion");

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(refv, this->equations->reference_state.v.data);
    YAKL_SCOPE(subtract_refstate, this->velocity_diffusion_subtract_refstate);

    parallel_for(
        "Velocity diffusion 0",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          for (int d = 0; d < ndims; ++d) {
            FTWvar(d, k + pks, j + pjs, i + pis, n) =
                Vvar(d, k + pks, j + pjs, i + pis, n);
            if (subtract_refstate) {
              FTWvar(d, k + pks, j + pjs, i + pis, n) -= refv(d, k + pks, n);
            }
          }
        });
    auxiliary_vars.exchange({FTWVAR});

    // *d*d = vort
    parallel_for(
        "Velocity diffusion 1",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D1_ext<1>(qhzedgereconvar, FTWvar, Wvar, pis, pjs, pks, i, j,
                            k, n);
        });
    auxiliary_vars.exchange({QHZEDGERECONVAR});

    parallel_for(
        "Velocity diffusion 1",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D1<1>(qxyedgereconvar, FTWvar, pis, pjs, pks, i, j, k, n);
        });
    auxiliary_vars.exchange({QXYEDGERECONVAR});

    parallel_for(
        "Velocity diffusion 2",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hnm11<1, diffusion_diff_ord>(qhzvar, qhzedgereconvar,
                                               primal_geometry, dual_geometry,
                                               dis, djs, dks, i, j, k + 1, n);
        });
    auxiliary_vars.exchange({QHZVAR});
    auxiliary_vars.fields_arr[QHZVAR].set_bnd(0.0);

    parallel_for(
        "Velocity diffusion 2",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hn0<1, diffusion_diff_ord>(qxyvar, qxyedgereconvar,
                                             primal_geometry, dual_geometry,
                                             dis, djs, dks, i, j, k, n);
        });
    auxiliary_vars.exchange({QXYVAR});

    parallel_for(
        SimpleBounds<3>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D1bar_ext<1>(Fvar, qhzvar, qxyvar, dis, djs, dks, i, j, k, 0);
        });
    auxiliary_vars.exchange({FVAR});

    parallel_for(
        SimpleBounds<3>(dual_topology.ni, dual_topology.n_cells_y,
                        dual_topology.n_cells_x),
        YAKL_LAMBDA(int k, int j, int i) {
          compute_D1bar<1>(FWvar, qhzvar, dis, djs, dks, i, j, k, 0);
        });
    auxiliary_vars.exchange({FWVAR});

    parallel_for(
        "Velocity diffusion 5",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> vdiff;
          compute_Hnm11bar<1, diffusion_diff_ord>(vdiff, Fvar, primal_geometry,
                                                  dual_geometry, pis, pjs, pks,
                                                  i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            Vtendvar(d, k + pks, j + pjs, i + pis, n) +=
                velocity_vort_horiz_diffusion_coeff * vdiff(d);
          }
        });

    parallel_for(
        "Velocity diffusion 6",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, 1> wdiff;
          compute_Hn0bar<1, vert_diffusion_diff_ord>(
              wdiff, FWvar, primal_geometry, dual_geometry, pis, pjs, pks, i, j,
              k, n);
          Wtendvar(0, k + pks, j + pjs, i + pis, n) +=
              velocity_vort_vert_diffusion_coeff * wdiff(0);
        });

    // d*d* = div
    parallel_for(
        "Velocity diffusion 7",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H10<1, diffusion_diff_ord>(Fvar, FTWvar, primal_geometry,
                                             dual_geometry, dis, djs, dks, i, j,
                                             k, n);
        });
    auxiliary_vars.exchange({FVAR});

    parallel_for(
        "Velocity diffusion 8",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H01<1, vert_diffusion_diff_ord>(FWvar, Wvar, primal_geometry,
                                                  dual_geometry, dis, djs, dks,
                                                  i, j, k + 1, n);
        });
    auxiliary_vars.exchange({FWVAR});
    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);

    parallel_for(
        "Velocity diffusion 9",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Dnm1bar<1>(Kvar, Fvar, dis, djs, dks, i, j, k, n);
          compute_Dnm1bar_vert<1, ADD_MODE::ADD>(Kvar, FWvar, dis, djs, dks, i,
                                                 j, k, n);
        });
    auxiliary_vars.exchange({KVAR});

    parallel_for(
        "Velocity diffusion 10",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hn1bar<1, diffusion_diff_ord, vert_diffusion_diff_ord>(
              dens0var, Kvar, primal_geometry, dual_geometry, pis, pjs, pks, i,
              j, k, n);
        });
    auxiliary_vars.exchange({DENS0VAR});

    parallel_for(
        "Velocity diffusion 11",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> vdiff;
          compute_D0<1>(vdiff, dens0var, pis, pjs, pks, i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            Vtendvar(d, pks + k, pjs + j, pis + i, n) -=
                velocity_div_horiz_diffusion_coeff * vdiff(0, d);
          }
        });

    parallel_for(
        "Velocity diffusion 12",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, 1, 1> wdiff;
          compute_D0_vert<1>(wdiff, dens0var, pis, pjs, pks, i, j, k, n);
          Wtendvar(0, pks + k, pjs + j, pis + i, n) -=
              velocity_div_vert_diffusion_coeff * wdiff(0);
        });

    yakl::timer_stop("add_velocity_diffusion");
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void compute_tendencies(
      real5d denstendvar, real5d Vtendvar, real5d Wtendvar,
      const real5d densreconvar, const real5d densvertreconvar,
      const real5d qhzreconvar, const real5d qhzvertreconvar,
      const real5d coriolishzreconvar, const real5d coriolishzvertreconvar,
      const real5d Bvar, const real5d Fvar, const real5d FWvar,
      optional_real5d opt_qxyreconvar, optional_real5d opt_coriolisxyreconvar) {
    yakl::timer_start("compute_tendencies");

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    const auto &refstate = this->equations->reference_state;
    YAKL_SCOPE(active_dens_ids, this->equations->varset.active_dens_ids);
    YAKL_SCOPE(force_refstate_hydrostatic_balance,
               this->force_refstate_hydrostatic_balance);

    real5d qxyreconvar, coriolisxyreconvar;
    if (ndims > 1) {
      qxyreconvar = opt_qxyreconvar.value();
      coriolisxyreconvar = opt_coriolisxyreconvar.value();
    }

    parallel_for(
        "Compute Wtend",
        SimpleBounds<4>(primal_topology.nl - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD0_vert<VS::ndensity_active, addmode>(
              Wtendvar, densvertreconvar, active_dens_ids, Bvar, pis, pjs, pks,
              i, j, k + 1, n);
          if (force_refstate_hydrostatic_balance) {
            compute_wD0_vert<VS::ndensity_active, ADD_MODE::ADD>(
                Wtendvar, refstate.q_di.data, active_dens_ids, refstate.B.data,
                pis, pjs, pks, i, j, k + 1, n);
          }
          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, qhzreconvar,
                                               qhzvertreconvar, Fvar, pis, pjs,
                                               pks, i, j, k + 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_w_nonEC<1, ADD_MODE::ADD>(
                Wtendvar, qhzreconvar, Fvar, pis, pjs, pks, i, j, k + 1, n);
          }

          compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, coriolishzreconvar,
                                             coriolishzvertreconvar, Fvar, pis,
                                             pjs, pks, i, j, k + 1, n);
          if (ndims > 1) {
            if (qf_choice == QF_MODE::EC) {
              compute_Qyz_w_EC<1, ADD_MODE::ADD>(Wtendvar, qhzreconvar,
                                                 qhzvertreconvar, Fvar, pis,
                                                 pjs, pks, i, j, k + 1, n);
            }
            if (qf_choice == QF_MODE::NOEC) {
              compute_Qyz_w_nonEC<1, ADD_MODE::ADD>(
                  Wtendvar, qhzreconvar, Fvar, pis, pjs, pks, i, j, k + 1, n);
            }

            compute_Qyz_w_EC<1, ADD_MODE::ADD>(Wtendvar, coriolishzreconvar,
                                               coriolishzvertreconvar, Fvar,
                                               pis, pjs, pks, i, j, k + 1, n);
          }
        });

    parallel_for(
        "Compute Wtend Bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          compute_wD0_vert<VS::ndensity_active, addmode>(
              Wtendvar, densvertreconvar, active_dens_ids, Bvar, pis, pjs, pks,
              i, j, 0, n);
          compute_wD0_vert<VS::ndensity_active, addmode>(
              Wtendvar, densvertreconvar, active_dens_ids, Bvar, pis, pjs, pks,
              i, j, primal_topology.nl - 1, n);
          if (force_refstate_hydrostatic_balance) {
            compute_wD0_vert<VS::ndensity_active, ADD_MODE::ADD>(
                Wtendvar, refstate.q_di.data, active_dens_ids, refstate.B.data,
                pis, pjs, pks, i, j, 0, n);
            compute_wD0_vert<VS::ndensity_active, ADD_MODE::ADD>(
                Wtendvar, refstate.q_di.data, active_dens_ids, refstate.B.data,
                pis, pjs, pks, i, j, primal_topology.nl - 1, n);
          }
          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(
                Wtendvar, qhzreconvar, qhzvertreconvar, Fvar, pis, pjs, pks, i,
                j, 0, n);
            compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(
                Wtendvar, qhzreconvar, qhzvertreconvar, Fvar, pis, pjs, pks, i,
                j, primal_topology.nl - 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_w_nonEC_bottom<1, ADD_MODE::ADD>(
                Wtendvar, qhzreconvar, Fvar, pis, pjs, pks, i, j, 0, n);
            compute_Qxz_w_nonEC_top<1, ADD_MODE::ADD>(
                Wtendvar, qhzreconvar, Fvar, pis, pjs, pks, i, j,
                primal_topology.nl - 1, n);
          }

          compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(
              Wtendvar, coriolishzreconvar, coriolishzvertreconvar, Fvar, pis,
              pjs, pks, i, j, 0, n);
          compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(
              Wtendvar, coriolishzreconvar, coriolishzvertreconvar, Fvar, pis,
              pjs, pks, i, j, primal_topology.nl - 1, n);

          if (ndims > 1) {
            if (qf_choice == QF_MODE::EC) {
              compute_Qyz_w_EC_bottom<1, ADD_MODE::ADD>(
                  Wtendvar, qhzreconvar, qhzvertreconvar, Fvar, pis, pjs, pks,
                  i, j, 0, n);
              compute_Qyz_w_EC_top<1, ADD_MODE::ADD>(
                  Wtendvar, qhzreconvar, qhzvertreconvar, Fvar, pis, pjs, pks,
                  i, j, primal_topology.nl - 1, n);
            }
            if (qf_choice == QF_MODE::NOEC) {
              compute_Qyz_w_nonEC_bottom<1, ADD_MODE::ADD>(
                  Wtendvar, qhzreconvar, Fvar, pis, pjs, pks, i, j, 0, n);
              compute_Qyz_w_nonEC_top<1, ADD_MODE::ADD>(
                  Wtendvar, qhzreconvar, Fvar, pis, pjs, pks, i, j,
                  primal_topology.nl - 1, n);
            }

            compute_Qyz_w_EC_bottom<1, ADD_MODE::ADD>(
                Wtendvar, coriolishzreconvar, coriolishzvertreconvar, Fvar, pis,
                pjs, pks, i, j, 0, n);
            compute_Qyz_w_EC_top<1, ADD_MODE::ADD>(
                Wtendvar, coriolishzreconvar, coriolishzvertreconvar, Fvar, pis,
                pjs, pks, i, j, primal_topology.nl - 1, n);
          }
        });

    parallel_for(
        "Compute Vtend",
        SimpleBounds<4>(primal_topology.ni - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD0<VS::ndensity_active, addmode>(Vtendvar, densreconvar,
                                                    active_dens_ids, Bvar, pis,
                                                    pjs, pks, i, j, k + 1, n);
          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar, qhzreconvar,
                                               qhzvertreconvar, FWvar, pis, pjs,
                                               pks, i, j, k + 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_u_nonEC<1, ADD_MODE::ADD>(Vtendvar, qhzvertreconvar,
                                                  FWvar, pis, pjs, pks, i, j,
                                                  k + 1, n);
          }

          compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar, coriolishzreconvar,
                                             coriolishzvertreconvar, FWvar, pis,
                                             pjs, pks, i, j, k + 1, n);

          if (ndims > 1) {
            if (qf_choice == QF_MODE::EC) {
              compute_Qyz_v_EC<1, ADD_MODE::ADD>(Vtendvar, qhzreconvar,
                                                 qhzvertreconvar, FWvar, pis,
                                                 pjs, pks, i, j, k + 1, n);
              compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, qxyreconvar, Fvar, pis,
                                             pjs, pks, i, j, k + 1, n);
            }
            if (qf_choice == QF_MODE::NOEC) {
              compute_Qyz_v_nonEC<1, ADD_MODE::ADD>(Vtendvar, qhzvertreconvar,
                                                    FWvar, pis, pjs, pks, i, j,
                                                    k + 1, n);
              compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, qxyreconvar, Fvar,
                                                pis, pjs, pks, i, j, k + 1, n);
            }

            compute_Qyz_v_EC<1, ADD_MODE::ADD>(Vtendvar, coriolishzreconvar,
                                               coriolishzvertreconvar, FWvar,
                                               pis, pjs, pks, i, j, k + 1, n);
            compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, coriolisxyreconvar, Fvar,
                                           pis, pjs, pks, i, j, k + 1, n);
          }
        });
    parallel_for(
        "Compute Vtend Bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          compute_wD0<VS::ndensity_active, addmode>(Vtendvar, densreconvar,
                                                    active_dens_ids, Bvar, pis,
                                                    pjs, pks, i, j, 0, n);
          compute_wD0<VS::ndensity_active, addmode>(
              Vtendvar, densreconvar, active_dens_ids, Bvar, pis, pjs, pks, i,
              j, primal_topology.ni - 1, n);

          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(
                Vtendvar, qhzreconvar, qhzvertreconvar, FWvar, pis, pjs, pks, i,
                j, 0, n);
            compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(
                Vtendvar, qhzreconvar, qhzvertreconvar, FWvar, pis, pjs, pks, i,
                j, primal_topology.ni - 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_u_nonEC_bottom<1, ADD_MODE::ADD>(
                Vtendvar, qhzvertreconvar, FWvar, pis, pjs, pks, i, j, 0, n);
            compute_Qxz_u_nonEC_top<1, ADD_MODE::ADD>(
                Vtendvar, qhzvertreconvar, FWvar, pis, pjs, pks, i, j,
                primal_topology.ni - 1, n);
          }

          compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(
              Vtendvar, coriolishzreconvar, coriolishzvertreconvar, FWvar, pis,
              pjs, pks, i, j, 0, n);
          compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(
              Vtendvar, coriolishzreconvar, coriolishzvertreconvar, FWvar, pis,
              pjs, pks, i, j, primal_topology.ni - 1, n);

          if (ndims > 1) {
            if (qf_choice == QF_MODE::EC) {
              compute_Qyz_v_EC_bottom<1, ADD_MODE::ADD>(
                  Vtendvar, qhzreconvar, qhzvertreconvar, FWvar, pis, pjs, pks,
                  i, j, 0, n);
              compute_Qyz_v_EC_top<1, ADD_MODE::ADD>(
                  Vtendvar, qhzreconvar, qhzvertreconvar, FWvar, pis, pjs, pks,
                  i, j, primal_topology.ni - 1, n);
              compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, qxyreconvar, Fvar, pis,
                                             pjs, pks, i, j, 0, n);
              compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, qxyreconvar, Fvar, pis,
                                             pjs, pks, i, j,
                                             primal_topology.ni - 1, n);
            }
            if (qf_choice == QF_MODE::NOEC) {
              compute_Qyz_v_nonEC_bottom<1, ADD_MODE::ADD>(
                  Vtendvar, qhzvertreconvar, FWvar, pis, pjs, pks, i, j, 0, n);
              compute_Qyz_v_nonEC_top<1, ADD_MODE::ADD>(
                  Vtendvar, qhzvertreconvar, FWvar, pis, pjs, pks, i, j,
                  primal_topology.ni - 1, n);
              compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, qxyreconvar, Fvar,
                                                pis, pjs, pks, i, j, 0, n);
              compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, qxyreconvar, Fvar,
                                                pis, pjs, pks, i, j,
                                                primal_topology.ni - 1, n);
            }

            compute_Qyz_v_EC_bottom<1, ADD_MODE::ADD>(
                Vtendvar, coriolishzreconvar, coriolishzvertreconvar, FWvar,
                pis, pjs, pks, i, j, 0, n);
            compute_Qyz_v_EC_top<1, ADD_MODE::ADD>(
                Vtendvar, coriolishzreconvar, coriolishzvertreconvar, FWvar,
                pis, pjs, pks, i, j, primal_topology.ni - 1, n);
            compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, coriolisxyreconvar, Fvar,
                                           pis, pjs, pks, i, j, 0, n);
            compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, coriolisxyreconvar, Fvar,
                                           pis, pjs, pks, i, j,
                                           primal_topology.ni - 1, n);
          }
        });

    parallel_for(
        "Compute Dens Tend",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDnm1bar<VS::ndensity_prognostic, addmode>(
              denstendvar, densreconvar, Fvar, dis, djs, dks, i, j, k, n);
          compute_wDnm1bar_vert<VS::ndensity_prognostic, ADD_MODE::ADD>(
              denstendvar, densvertreconvar, FWvar, dis, djs, dks, i, j, k, n);
        });

    yakl::timer_stop("compute_tendencies");
  }

  real compute_max_anelastic_constraint(FieldSet<nprognostic> &x,
                                        FieldSet<nauxiliary> &auxiliary_vars,
                                        bool has_f_and_fw = false) override {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    const auto mfvar = auxiliary_vars.fields_arr[MFVAR].data;

    const auto &densvar = x.fields_arr[DENSVAR].data;
    const auto &Vvar = x.fields_arr[VVAR].data;
    const auto &Wvar = x.fields_arr[WVAR].data;

    const auto &Fvar = auxiliary_vars.fields_arr[FVAR].data;
    const auto &FWvar = auxiliary_vars.fields_arr[FWVAR].data;

    YAKL_SCOPE(Hk, this->equations->Hk);
    YAKL_SCOPE(varset, this->equations->varset);
    YAKL_SCOPE(rho_pi, this->equations->reference_state.rho_pi.data);
    YAKL_SCOPE(rho_di, this->equations->reference_state.rho_di.data);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);

    if (!has_f_and_fw) {
      parallel_for(
          "Compute F and Fw",
          SimpleBounds<4>(dual_topology.ni, dual_topology.n_cells_y,
                          dual_topology.n_cells_x, dual_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            SArray<real, 1, ndims> he, u;

            // no need to compute K2 here but this is only for diagnostic and
            // maybe compiler can optimize it away ?
            real hew, uw, K2;
            Hk.compute_he_U_and_K(he, hew, u, uw, K2, densvar, Vvar, Wvar, pis,
                                  pjs, pks, i, j, k, n);

            for (int d = 0; d < ndims; ++d) {
              Fvar(d, pks + k, pjs + j, pis + i, n) = he(d) * u(d);
            }
            FWvar(0, pks + k, pjs + j, pis + i, n) = hew * uw;
          });
      auxiliary_vars.exchange({FVAR, FWVAR});
    } else {
    }
    parallel_for(
        "Compute anelastic constraint",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Dnm1bar<1>(mfvar, Fvar, dis, djs, dks, i, j, k, n);
          compute_Dnm1bar_vert<1, ADD_MODE::ADD>(mfvar, FWvar, dis, djs, dks, i,
                                                 j, k, n);
          const real total_dens =
              varset.get_total_density(densvar, k, j, i, dks, djs, dis, n);
          mfvar(0, k + dks, j + djs, i + dis, n) /= total_dens;
        });

    auxiliary_vars.exchange({MFVAR});
    const auto divvar = auxiliary_vars.fields_arr[MFVAR].data.slice<4>(
        0, yakl::COLON, yakl::COLON, yakl::COLON, yakl::COLON);

    return yakl::intrinsics::maxval(yakl::intrinsics::abs(divvar));
  }

  void compute_functional_derivatives(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x,
      FieldSet<nauxiliary> &auxiliary_vars, real fac = 1,
      ADD_MODE addmode = ADD_MODE::REPLACE) override {
    yakl::timer_start("compute_functional_derivatives");

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    const auto &densvar = x.fields_arr[DENSVAR].data;
    const auto &Vvar = x.fields_arr[VVAR].data;
    const auto &Wvar = x.fields_arr[WVAR].data;

    const auto &Fvar = auxiliary_vars.fields_arr[FVAR].data;
    const auto &FWvar = auxiliary_vars.fields_arr[FWVAR].data;
    const auto &Kvar = auxiliary_vars.fields_arr[KVAR].data;
    const auto &Bvar = auxiliary_vars.fields_arr[BVAR].data;
    const auto &HSvar = const_vars.fields_arr[HSVAR].data;

    YAKL_SCOPE(Hk, this->equations->Hk);
    YAKL_SCOPE(Hs, this->equations->Hs);
    YAKL_SCOPE(rho_pi, this->equations->reference_state.rho_pi.data);
    YAKL_SCOPE(rho_di, this->equations->reference_state.rho_di.data);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);

    parallel_for(
        "Functional derivatives",
        SimpleBounds<4>(dual_topology.ni, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> he, u;
          real hew, uw, K2;
          Hk.compute_he_U_and_K(he, hew, u, uw, K2, densvar, Vvar, Wvar, pis,
                                pjs, pks, i, j, k, n);

          Kvar(0, k + pks, j + pjs, i + pis, n) = K2;

          if (addmode == ADD_MODE::ADD) {
            for (int d = 0; d < ndims; ++d) {
              Fvar(d, pks + k, pjs + j, pis + i, n) += fac * he(d) * u(d);
            }
            FWvar(0, pks + k, pjs + j, pis + i, n) += fac * hew * uw;
          } else if (addmode == ADD_MODE::REPLACE) {
            for (int d = 0; d < ndims; ++d) {
              Fvar(d, pks + k, pjs + j, pis + i, n) = fac * he(d) * u(d);
            }
            FWvar(0, pks + k, pjs + j, pis + i, n) = fac * hew * uw;
          }

          if (k < primal_topology.ni) {
            if (addmode == ADD_MODE::ADD) {
              Hs.compute_dHsdx<ADD_MODE::ADD>(Bvar, densvar, HSvar, pis, pjs,
                                              pks, i, j, k, n, fac);
            } else if (addmode == ADD_MODE::REPLACE) {
              Hs.compute_dHsdx<ADD_MODE::REPLACE>(Bvar, densvar, HSvar, pis,
                                                  pjs, pks, i, j, k, n, fac);
            }
          }
        });
    auxiliary_vars.exchange({KVAR});

    parallel_for(
        "Add K to B",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j, k,
                                            n, fac);
        });

    auxiliary_vars.exchange({BVAR, FVAR, FWVAR});

#if defined PAMC_AN || defined PAMC_MAN
    if (this->check_anelastic_constraint) {
      real max_div = compute_max_anelastic_constraint(x, auxiliary_vars, true);
      std::cout << "Anelastic constraint: " << max_div << std::endl;
    }
#endif
    yakl::timer_stop("compute_functional_derivatives");
  }

  void compute_two_point_discrete_gradient(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x1,
      FieldSet<nprognostic> &x2,
      FieldSet<nauxiliary> &auxiliary_vars) override {
    compute_two_point_discrete_gradient_impl(
        this->equations->Hs, this->equations->thermo, dt, const_vars, x1, x2,
        auxiliary_vars);
  }

  template <class HamilT = Hamiltonian, class ThermoT = ThermoPotential>
  void compute_two_point_discrete_gradient_impl(
      HamilT Hs, ThermoT thermo, real dt, FieldSet<nconstant> &const_vars,
      FieldSet<nprognostic> &x1, FieldSet<nprognostic> &x2,
      FieldSet<nauxiliary> &auxiliary_vars) {
    if constexpr (!two_point_discrete_gradient_implemented_v<HamilT, ThermoT>) {
      throw std::runtime_error(
          "two point discrete gradient not implemented for this combination of "
          "hamiltonian and thermodynamics");
    } else {
      const auto &primal_topology = primal_geometry.topology;
      const auto &dual_topology = dual_geometry.topology;

      int pis = primal_topology.is;
      int pjs = primal_topology.js;
      int pks = primal_topology.ks;

      const auto &densvar1 = x1.fields_arr[DENSVAR].data;
      const auto &Vvar1 = x1.fields_arr[VVAR].data;
      const auto &Wvar1 = x1.fields_arr[WVAR].data;

      const auto &densvar2 = x2.fields_arr[DENSVAR].data;
      const auto &Vvar2 = x2.fields_arr[VVAR].data;
      const auto &Wvar2 = x2.fields_arr[WVAR].data;

      const auto &Fvar = auxiliary_vars.fields_arr[FVAR].data;
      const auto &FWvar = auxiliary_vars.fields_arr[FWVAR].data;
      const auto &Kvar = auxiliary_vars.fields_arr[KVAR].data;
      const auto &Bvar = auxiliary_vars.fields_arr[BVAR].data;
      const auto &HSvar = const_vars.fields_arr[HSVAR].data;

      YAKL_SCOPE(Hk, this->equations->Hk);
      YAKL_SCOPE(primal_geometry, this->primal_geometry);
      YAKL_SCOPE(dual_geometry, this->dual_geometry);
      parallel_for(
          "Functional derivatives",
          SimpleBounds<4>(dual_topology.ni, dual_topology.n_cells_y,
                          dual_topology.n_cells_x, dual_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            SArray<real, 1, ndims> he_1, u_1;
            real hew_1, uw_1, K2_1;
            Hk.compute_he_U_and_K(he_1, hew_1, u_1, uw_1, K2_1, densvar1, Vvar1,
                                  Wvar1, pis, pjs, pks, i, j, k, n);

            SArray<real, 1, ndims> he_2, u_2;
            real hew_2, uw_2, K2_2;
            Hk.compute_he_U_and_K(he_2, hew_2, u_2, uw_2, K2_2, densvar2, Vvar2,
                                  Wvar2, pis, pjs, pks, i, j, k, n);

            Kvar(0, k + pks, j + pjs, i + pis, n) = 0.5_fp * (K2_1 + K2_2);

            for (int d = 0; d < ndims; ++d) {
              Fvar(d, pks + k, pjs + j, pis + i, n) =
                  0.25_fp * (he_1(d) + he_2(d)) * (u_1(d) + u_2(d));
            }
            FWvar(0, pks + k, pjs + j, pis + i, n) =
                0.25_fp * (hew_1 + hew_2) * (uw_1 + uw_2);

            if (k < primal_topology.ni) {
              Hs.compute_dHsdx_two_point(thermo, Bvar, densvar1, densvar2,
                                         HSvar, pis, pjs, pks, i, j, k, n);
            }
          });
      auxiliary_vars.exchange({KVAR});

      parallel_for(
          "Add K to B",
          SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j,
                                              k, n);
          });

      auxiliary_vars.exchange({BVAR, FVAR, FWVAR});
    }
  }

  void apply_symplectic(real dt, FieldSet<nconstant> &const_vars,
                        FieldSet<nprognostic> &x,
                        FieldSet<nauxiliary> &auxiliary_vars,
                        FieldSet<nprognostic> &xtend,
                        ADD_MODE addmode = ADD_MODE::REPLACE,
                        bool needs_to_recompute_F = true) override {
    yakl::timer_start("apply_symplectic");

    const auto &dual_topology = dual_geometry.topology;

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);

    auxiliary_vars.exchange({DENS0VAR});

    if (needs_to_recompute_F) {
      compute_F_and_FW(auxiliary_vars.fields_arr[F2VAR].data,
                       auxiliary_vars.fields_arr[FW2VAR].data,
                       x.fields_arr[DENSVAR].data, x.fields_arr[VVAR].data,
                       x.fields_arr[WVAR].data);

      auxiliary_vars.exchange({F2VAR, FW2VAR});
      auxiliary_vars.fields_arr[FW2VAR].set_bnd(0.0);
    }

    compute_FT_and_FTW(
        auxiliary_vars.fields_arr[FTVAR].data,
        auxiliary_vars.fields_arr[FTWVAR].data,
        needs_to_recompute_F ? auxiliary_vars.fields_arr[F2VAR].data
                             : auxiliary_vars.fields_arr[FVAR].data,
        needs_to_recompute_F ? auxiliary_vars.fields_arr[FW2VAR].data
                             : auxiliary_vars.fields_arr[FWVAR].data,
        ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[FTXYVAR].data}
                  : std::nullopt);

    auxiliary_vars.exchange({FTVAR, FTWVAR});
    if (ndims > 1) {
      auxiliary_vars.exchange({FTXYVAR});
    }

    compute_q_and_f(
        auxiliary_vars.fields_arr[QHZVAR].data,
        auxiliary_vars.fields_arr[FHZVAR].data, x.fields_arr[VVAR].data,
        x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data,
        const_vars.fields_arr[CORIOLISHZVAR].data,
        ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[QXYVAR].data}
                  : std::nullopt,
        ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[FXYVAR].data}
                  : std::nullopt,
        ndims > 1 ? optional_real5d{const_vars.fields_arr[CORIOLISXYVAR].data}
                  : std::nullopt);

    auxiliary_vars.exchange({QHZVAR, FHZVAR});
    auxiliary_vars.fields_arr[QHZVAR].set_bnd(0.0);
    auxiliary_vars.fields_arr[FHZVAR].set_bnd(0.0);
    if (ndims > 1) {
      auxiliary_vars.exchange({QXYVAR, FXYVAR});
    }

    // Compute densrecon, densvertrecon, qrecon and frecon
    if (dual_geometry.uniform_vertical) {
      compute_edge_reconstructions_uniform(
          auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QHZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QHZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[QHZVAR].data,
          auxiliary_vars.fields_arr[FHZVAR].data,
          ndims > 1
              ? optional_real5d{auxiliary_vars.fields_arr[QXYEDGERECONVAR].data}
              : std::nullopt,
          ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[QXYVAR].data}
                    : std::nullopt,
          ndims > 1 ? optional_real5d{auxiliary_vars
                                          .fields_arr[CORIOLISXYEDGERECONVAR]
                                          .data}
                    : std::nullopt,
          ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[FXYVAR].data}
                    : std::nullopt);
    } else {
      compute_edge_reconstructions_variable(
          auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QHZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QHZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[QHZVAR].data,
          auxiliary_vars.fields_arr[FHZVAR].data,
          ndims > 1
              ? optional_real5d{auxiliary_vars.fields_arr[QXYEDGERECONVAR].data}
              : std::nullopt,
          ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[QXYVAR].data}
                    : std::nullopt,
          ndims > 1 ? optional_real5d{auxiliary_vars
                                          .fields_arr[CORIOLISXYEDGERECONVAR]
                                          .data}
                    : std::nullopt,
          ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[FXYVAR].data}
                    : std::nullopt);
    }

    auxiliary_vars.exchange({DENSEDGERECONVAR, DENSVERTEDGERECONVAR,
                             QHZEDGERECONVAR, QHZVERTEDGERECONVAR,
                             CORIOLISHZEDGERECONVAR,
                             CORIOLISHZVERTEDGERECONVAR});
    if (ndims > 1) {
      auxiliary_vars.exchange({QXYEDGERECONVAR, CORIOLISXYEDGERECONVAR});
    }

    compute_recons(
        auxiliary_vars.fields_arr[DENSRECONVAR].data,
        auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
        auxiliary_vars.fields_arr[QHZRECONVAR].data,
        auxiliary_vars.fields_arr[QHZVERTRECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISHZRECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISHZVERTRECONVAR].data,
        auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[QHZEDGERECONVAR].data,
        auxiliary_vars.fields_arr[QHZVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISHZEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISHZVERTEDGERECONVAR].data,
        x.fields_arr[DENSVAR].data,
        needs_to_recompute_F ? auxiliary_vars.fields_arr[F2VAR].data
                             : auxiliary_vars.fields_arr[FVAR].data,
        needs_to_recompute_F ? auxiliary_vars.fields_arr[FW2VAR].data
                             : auxiliary_vars.fields_arr[F2VAR].data,
        auxiliary_vars.fields_arr[FTVAR].data,
        auxiliary_vars.fields_arr[FTWVAR].data, this->tanh_upwind_coeff,
        ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[QXYRECONVAR].data}
                  : std::nullopt,
        ndims > 1
            ? optional_real5d{auxiliary_vars.fields_arr[QXYEDGERECONVAR].data}
            : std::nullopt,
        ndims > 1
            ? optional_real5d{auxiliary_vars.fields_arr[CORIOLISXYRECONVAR]
                                  .data}
            : std::nullopt,
        ndims > 1
            ? optional_real5d{auxiliary_vars.fields_arr[CORIOLISXYEDGERECONVAR]
                                  .data}
            : std::nullopt,
        ndims > 1 ? optional_real5d{auxiliary_vars.fields_arr[FTXYVAR].data}
                  : std::nullopt);

    auxiliary_vars.exchange({DENSRECONVAR, DENSVERTRECONVAR, QHZRECONVAR,
                             QHZVERTRECONVAR, CORIOLISHZRECONVAR,
                             CORIOLISHZVERTRECONVAR});

    if (ndims > 1) {
      auxiliary_vars.exchange({QXYRECONVAR, CORIOLISXYRECONVAR});
    }

    // Compute fct
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(dens_pos, equations->varset.dens_pos);
    parallel_for(
        "ComputeEdgeFlux",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_edgefluxes<VS::ndensity_prognostic>(
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[DENSRECONVAR].data,
              auxiliary_vars.fields_arr[FVAR].data, dens_pos, dis, djs, dks, i,
              j, k, n);

          if (k < dual_topology.ni - 2) {
            compute_vertedgefluxes<VS::ndensity_prognostic>(
                auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data,
                auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                auxiliary_vars.fields_arr[FWVAR].data, dens_pos, dis, djs, dks,
                i, j, k + 1, n);
          }
        });
    auxiliary_vars.exchange({EDGEFLUXVAR, VERTEDGEFLUXVAR});

    parallel_for(
        "ComputeMf",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Mfext<VS::ndensity_prognostic>(
              auxiliary_vars.fields_arr[MFVAR].data,
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, dt, dens_pos,
              dis, djs, dks, i, j, k, n);
        });

    auxiliary_vars.exchange({MFVAR});

    parallel_for(
        "Apply Phi",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          apply_Phi<VS::ndensity_prognostic>(
              auxiliary_vars.fields_arr[DENSRECONVAR].data,
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data,
              dens_pos, dis, djs, dks, i, j, k, n);
          if (k < dual_topology.ni - 2) {
            apply_Phivert<VS::ndensity_prognostic>(
                auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data,
                auxiliary_vars.fields_arr[MFVAR].data,
                x.fields_arr[DENSVAR].data, dens_pos, dis, djs, dks, i, j,
                k + 1, n);
          }
        });

    auxiliary_vars.exchange({DENSRECONVAR, DENSVERTRECONVAR});

    // Compute tendencies
    if (addmode == ADD_MODE::REPLACE) {
      compute_tendencies<ADD_MODE::REPLACE>(
          xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
          xtend.fields_arr[WVAR].data,
          auxiliary_vars.fields_arr[DENSRECONVAR].data,
          auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
          auxiliary_vars.fields_arr[QHZRECONVAR].data,
          auxiliary_vars.fields_arr[QHZVERTRECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZRECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZVERTRECONVAR].data,
          auxiliary_vars.fields_arr[BVAR].data,
          auxiliary_vars.fields_arr[FVAR].data,
          auxiliary_vars.fields_arr[FWVAR].data,
          ndims > 1
              ? optional_real5d{auxiliary_vars.fields_arr[QXYRECONVAR].data}
              : std::nullopt,
          ndims > 1
              ? optional_real5d{auxiliary_vars.fields_arr[CORIOLISXYRECONVAR]
                                    .data}
              : std::nullopt);
    }
    if (addmode == ADD_MODE::ADD) {
      compute_tendencies<ADD_MODE::ADD>(
          xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
          xtend.fields_arr[WVAR].data,
          auxiliary_vars.fields_arr[DENSRECONVAR].data,
          auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
          auxiliary_vars.fields_arr[QHZRECONVAR].data,
          auxiliary_vars.fields_arr[QHZVERTRECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZRECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISHZVERTRECONVAR].data,
          auxiliary_vars.fields_arr[BVAR].data,
          auxiliary_vars.fields_arr[FVAR].data,
          auxiliary_vars.fields_arr[FWVAR].data,
          ndims > 1
              ? optional_real5d{auxiliary_vars.fields_arr[QXYRECONVAR].data}
              : std::nullopt,
          ndims > 1
              ? optional_real5d{auxiliary_vars.fields_arr[CORIOLISXYRECONVAR]
                                    .data}
              : std::nullopt);
    }

    if (scalar_horiz_diffusion_coeff > 0 or scalar_vert_diffusion_coeff > 0) {
      add_scalar_diffusion(
          scalar_horiz_diffusion_coeff, scalar_vert_diffusion_coeff,
          xtend.fields_arr[DENSVAR].data, x.fields_arr[DENSVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[FDIFFVAR].data,
          auxiliary_vars.fields_arr[FWDIFFVAR].data, auxiliary_vars);
    }
    if (velocity_vort_horiz_diffusion_coeff > 0 or
        velocity_vort_vert_diffusion_coeff > 0 or
        velocity_div_horiz_diffusion_coeff > 0 or
        velocity_div_vert_diffusion_coeff > 0) {
      if (ndims == 1) {
        add_velocity_diffusion_2d(
            velocity_vort_horiz_diffusion_coeff,
            velocity_vort_vert_diffusion_coeff,
            velocity_div_horiz_diffusion_coeff,
            velocity_div_vert_diffusion_coeff, xtend.fields_arr[VVAR].data,
            xtend.fields_arr[WVAR].data, x.fields_arr[VVAR].data,
            x.fields_arr[WVAR].data,
            auxiliary_vars.fields_arr[QHZEDGERECONVAR].data,
            auxiliary_vars.fields_arr[QHZVAR].data,
            auxiliary_vars.fields_arr[DENS0VAR].data,
            auxiliary_vars.fields_arr[KVAR].data,
            auxiliary_vars.fields_arr[FVAR].data,
            auxiliary_vars.fields_arr[FWVAR].data,
            auxiliary_vars.fields_arr[FTWVAR].data, auxiliary_vars);
      } else {
        add_velocity_diffusion_3d(
            velocity_vort_horiz_diffusion_coeff,
            velocity_vort_vert_diffusion_coeff,
            velocity_div_horiz_diffusion_coeff,
            velocity_div_vert_diffusion_coeff, xtend.fields_arr[VVAR].data,
            xtend.fields_arr[WVAR].data, x.fields_arr[VVAR].data,
            x.fields_arr[WVAR].data,
            auxiliary_vars.fields_arr[QHZEDGERECONVAR].data,
            auxiliary_vars.fields_arr[QHZVAR].data,
            auxiliary_vars.fields_arr[DENS0VAR].data,
            auxiliary_vars.fields_arr[KVAR].data,
            auxiliary_vars.fields_arr[FVAR].data,
            auxiliary_vars.fields_arr[FWVAR].data,
            auxiliary_vars.fields_arr[FTWVAR].data,
            auxiliary_vars.fields_arr[QXYEDGERECONVAR].data,
            auxiliary_vars.fields_arr[QXYVAR].data, auxiliary_vars);
      }
    }

    yakl::timer_stop("apply_symplectic");
  }

#if defined PAMC_AN || defined PAMC_MAN
  void project_to_anelastic(FieldSet<nconstant> &const_vars,
                            FieldSet<nprognostic> &x,
                            FieldSet<nauxiliary> &auxiliary_vars) override {
    // this happens to do what we want, note x used for xtend
    add_pressure_perturbation(1, const_vars, x, auxiliary_vars, x);
  }

  void add_pressure_perturbation(real dt, FieldSet<nconstant> &const_vars,
                                 FieldSet<nprognostic> &x,
                                 FieldSet<nauxiliary> &auxiliary_vars,
                                 FieldSet<nprognostic> &xtend) override {
    yakl::timer_start("add_pressure_perturbation");
    linear_system->solve(dt, xtend, const_vars, auxiliary_vars, xtend);
    yakl::timer_stop("add_pressure_perturbation");
  }
#endif

  void clip_negative_densities(FieldSet<nprognostic> &x) override {
    const auto &dual_topology = dual_geometry.topology;

    const auto densvar = x.fields_arr[DENSVAR].data;
    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(dens_pos, equations->varset.dens_pos);
    parallel_for(
        "Remove negatives",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          for (int d = 0; d < VS::ndensity_prognostic; ++d) {
            if (dens_pos(d)) {
              densvar(d, k + dks, j + djs, i + dis, n) =
                  std::max(0._fp, densvar(d, k + dks, j + djs, i + dis, n));
            }
          }
        });
  }
};

// *******   Linear system   ***********//
class CompressibleVelocityLinearSystem : public LinearSystem {

  yakl::RealFFT1D<real> fftv_x;
  yakl::RealFFT1D<real> fftw_x;
  // yakl::RealFFT1D<real> fftv_y;
  // yakl::RealFFT1D<real> fftw_y;

  int nxf, nyf;

  real4d Blin_coeff;
  real4d v_transform;
  real4d w_transform;
  complex4d complex_vrhs;
  complex4d complex_wrhs;
  complex5d complex_vcoeff;

  complex4d tri_l;
  complex4d tri_d;
  complex4d tri_u;

  complex4d tri_c;

  using VS = VariableSet;

public:
  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  Equations &equations) override {

    if (ndims > 1) {
      throw std::runtime_error(
          "velocity-based linear solver not implemented in 3d");
    }

    LinearSystem::initialize(params, primal_geom, dual_geom, equations);

    const auto &primal_topology = primal_geom.topology;

    auto pni = primal_topology.ni;
    auto pnl = primal_topology.nl;
    auto nx = primal_topology.n_cells_x;
    auto ny = primal_topology.n_cells_y;
    auto nens = primal_topology.nens;

    this->nxf = nx + 2 - nx % 2;
    this->nyf = ndims > 1 ? ny + 2 - ny % 2 : ny;

    this->Blin_coeff = real4d("Blin coeff", VS::ndensity_dycore,
                              VS::ndensity_dycore, pni, nens);

    v_transform = real4d("v transform", pni, nyf, nxf, nens);
    w_transform = real4d("w transform", pnl, nyf, nxf, nens);
    yakl::memset(v_transform, 0);
    yakl::memset(w_transform, 0);

    complex_vrhs = complex4d("complex vrhs", pni, ny, nx, nens);
    complex_wrhs = complex4d("complex wrhs", pnl, ny, nx, nens);

    fftv_x.init(v_transform, 2, nx);
    fftw_x.init(w_transform, 2, nx);
    // fftv_y.init(v_transform, 1, ny);
    // fftw_y.init(w_transform, 1, ny);

    complex_vcoeff =
        complex5d("complex vcoeff", 1 + VS::ndensity_dycore, pni, ny, nx, nens);

    tri_d = complex4d("tri d", pnl, ny, nx, nens);
    tri_l = complex4d("tri l", pnl, ny, nx, nens);
    tri_u = complex4d("tri u", pnl, ny, nx, nens);
    tri_c = complex4d("tri c", pnl, ny, nx, nens);
  }

  virtual void compute_coefficients(real dt) override {
    const auto &refstate = this->equations->reference_state;

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    auto n_cells_x = dual_topology.n_cells_x;
    auto n_cells_y = dual_topology.n_cells_y;
    auto nens = dual_topology.nens;
    auto nl = primal_topology.nl;
    auto ni = primal_topology.ni;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    real dtf = dt / 2;
    real dtf2 = dt * dt / 4;

    const auto &rho_pi = refstate.rho_pi.data;
    const auto &q_pi = refstate.q_pi.data;
    const auto &rho_di = refstate.rho_di.data;
    const auto &q_di = refstate.q_di.data;
    const auto &Nsq_pi = refstate.Nsq_pi.data;

    YAKL_SCOPE(thermo, this->equations->thermo);
    YAKL_SCOPE(grav, this->equations->Hs.g);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(Blin_coeff, this->Blin_coeff);
    YAKL_SCOPE(tri_l, this->tri_l);
    YAKL_SCOPE(tri_d, this->tri_d);
    YAKL_SCOPE(tri_u, this->tri_u);
    YAKL_SCOPE(complex_vcoeff, this->complex_vcoeff);

    parallel_for(
        "Compute Blin_coeff",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          real z = primal_geometry.zint(k + primal_topology.ks, n);

          real Rd = thermo.cst.Rd;
          real pr = thermo.cst.pr;
          real gamma_d = thermo.cst.gamma_d;
          real Cpd = thermo.cst.Cpd;
          real Cvd = thermo.cst.Cvd;
          real grav2 = grav * grav;

          real rho_ref = rho_pi(0, k + pks, n);
          real alpha_ref = 1 / rho_ref;
          real s_ref = q_pi(1, k + pks, n);

          real rho_ref2 = rho_ref * rho_ref;
          real p_ref = thermo.solve_p(rho_ref, s_ref, 0, 0, 0, 0);

          real dpds_ref =
              thermo.compute_dpdentropic_var(alpha_ref, s_ref, 0, 0, 0, 0);
          real dpds_ref2 = dpds_ref * dpds_ref;

          real Nref2 = Nsq_pi(0, k + pks, n);
          real cref = thermo.compute_soundspeed(alpha_ref, s_ref, 0, 0, 0, 0);
          real cref2 = cref * cref;

          real b0_rho = (cref2 * rho_ref - dpds_ref * s_ref) / rho_ref2;
          real b0_s =
              dpds_ref / rho_ref - dpds_ref2 * s_ref / (cref2 * rho_ref2) -
              dpds_ref2 * grav2 * s_ref / (Nref2 * cref2 * cref2 * rho_ref2);
          real b0_S = b0_s / rho_ref;
          b0_rho -= s_ref / rho_ref * b0_s;

          real b1_rho = dpds_ref / rho_ref2;
          real b1_s = dpds_ref2 * (Nref2 * cref2 + grav2) /
                      (Nref2 * cref2 * cref2 * rho_ref2);
          real b1_S = b1_s / rho_ref;
          b1_rho -= s_ref / rho_ref * b1_s;

          Blin_coeff(0, 0, k, n) = b0_rho;
          Blin_coeff(0, 1, k, n) = b0_S;
          Blin_coeff(1, 0, k, n) = b1_rho;
          Blin_coeff(1, 1, k, n) = b1_S;
        });

    parallel_for(
        "compute vcoeff",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> fD0Dbar;
          fourier_cwD0Dnm1bar(fD0Dbar, 1, i, j, k, dual_topology.n_cells_x,
                              dual_topology.n_cells_y, dual_topology.ni);

          real fH2bar = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          SArray<real, 1, ndims> fH1;
          fourier_H10<diff_ord>(fH1, primal_geometry, dual_geometry, pis, pjs,
                                pks, i, j, k, 0, n_cells_x, n_cells_y,
                                dual_topology.ni);
          SArray<complex, 1, ndims> fD0;
          fourier_cwD0(fD0, 1, i, j, k, dual_topology.n_cells_x,
                       dual_topology.n_cells_y, dual_topology.ni);

          real he = rho_pi(0, k + pks, n);

          real c1 = 1;
          for (int d1 = 0; d1 < VS::ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < VS::ndensity_dycore; ++d2) {
              c1 -= dtf2 * fH2bar * fH1(0) * fD0Dbar(0) * he *
                    q_pi(d1, k + pks, n) * q_pi(d2, k + pks, n) *
                    Blin_coeff(d1, d2, k, n);
            }
          }

          complex_vcoeff(0, k, j, i, n) = 1 / c1;
          for (int d1 = 0; d1 < VS::ndensity_dycore; ++d1) {
            complex cd1 = 0;
            for (int d2 = 0; d2 < VS::ndensity_dycore; ++d2) {
              cd1 += fD0(0) * dtf2 * fH2bar * q_pi(d2, k + pks, n) *
                     Blin_coeff(d2, d1, k, n);
            }
            complex_vcoeff(1 + d1, k, j, i, n) = cd1 / c1;
          }
        });

    parallel_for(
        "Compute vertical tridiag",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real fH2bar_k = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fH2bar_kp1 = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          real gamma_fac_kp2 = rho_di(0, k + dks + 2, n) *
                               H01_diagonal(primal_geometry, dual_geometry, pis,
                                            pjs, pks, i, j, k + 2, n);
          real gamma_fac_kp1 = rho_di(0, k + dks + 1, n) *
                               H01_diagonal(primal_geometry, dual_geometry, pis,
                                            pjs, pks, i, j, k + 1, n);
          real gamma_fac_k = rho_di(0, k + dks, n) *
                             H01_diagonal(primal_geometry, dual_geometry, pis,
                                          pjs, pks, i, j, k, n);

          tri_u(k, j, i, n) = 0;
          tri_d(k, j, i, n) = 1;
          tri_l(k, j, i, n) = 0;

          for (int d1 = 0; d1 < VS::ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < VS::ndensity_dycore; ++d2) {
              real alpha_kp1 = q_di(d1, k + dks + 1, n);

              real beta_kp1 = fH2bar_kp1 * Blin_coeff(d1, d2, k + 1, n);
              real beta_k = fH2bar_k * Blin_coeff(d1, d2, k, n);

              real gamma_kp2 = gamma_fac_kp2 * q_di(d2, k + dks + 2, n);
              real gamma_kp1 = gamma_fac_kp1 * q_di(d2, k + dks + 1, n);
              real gamma_k = gamma_fac_k * q_di(d2, k + dks, n);

              tri_u(k, j, i, n) += -dtf2 * alpha_kp1 * beta_kp1 * gamma_kp2;
              tri_d(k, j, i, n) +=
                  dtf2 * alpha_kp1 * (beta_kp1 + beta_k) * gamma_kp1;
              tri_l(k, j, i, n) += -dtf2 * alpha_kp1 * beta_k * gamma_k;
            }
          }
        });

    parallel_for(
        "Compute horizontal tridiag",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real fH2bar_k = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fH2bar_kp1 = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          real gamma_fac_kp2 = rho_di(0, k + dks + 2, n) *
                               H01_diagonal(primal_geometry, dual_geometry, pis,
                                            pjs, pks, i, j, k + 2, n);
          real gamma_fac_kp1 = rho_di(0, k + dks + 1, n) *
                               H01_diagonal(primal_geometry, dual_geometry, pis,
                                            pjs, pks, i, j, k + 1, n);
          real gamma_fac_k = rho_di(0, k + dks, n) *
                             H01_diagonal(primal_geometry, dual_geometry, pis,
                                          pjs, pks, i, j, k, n);

          SArray<real, 1, ndims> fH1_kp1_a;
          SArray<real, 1, ndims> fH1_k_a;
          fourier_H10<diff_ord>(fH1_kp1_a, primal_geometry, dual_geometry, pis,
                                pjs, pks, i, j, k + 1, 0, n_cells_x, n_cells_y,
                                dual_topology.ni);
          fourier_H10<diff_ord>(fH1_k_a, primal_geometry, dual_geometry, pis,
                                pjs, pks, i, j, k, 0, n_cells_x, n_cells_y,
                                dual_topology.ni);
          real fH1h_kp1 = fH1_kp1_a(0);
          real fH1h_k = fH1_k_a(0);

          complex fDnm1bar_kp1 = fourier_Dnm1bar(1, i, j, k + 1, n_cells_x,
                                                 n_cells_y, dual_topology.ni);
          complex fDnm1bar_k = fourier_Dnm1bar(1, i, j, k, n_cells_x, n_cells_y,
                                               dual_topology.ni);

          real he_kp1 = rho_pi(0, k + pks + 1, n);
          real he_k = rho_pi(0, k + pks, n);

          for (int d1 = 0; d1 < VS::ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < VS::ndensity_dycore; ++d2) {
              for (int d3 = 0; d3 < VS::ndensity_dycore; ++d3) {

                real alpha_kp1 = dtf2 * q_di(d1, k + dks + 1, n);
                complex beta_kp1 = fH2bar_kp1 * Blin_coeff(d1, d2, k + 1, n) *
                                   q_pi(d2, k + dks + 1, n) * fDnm1bar_kp1 *
                                   he_kp1 * fH1h_kp1;
                complex beta_k = fH2bar_k * Blin_coeff(d1, d2, k, n) *
                                 q_pi(d2, k + pks, n) * fDnm1bar_k * he_k *
                                 fH1h_k;

                real gamma_kp2 = gamma_fac_kp2 * q_di(d3, k + dks + 2, n);
                real gamma_kp1 = gamma_fac_kp1 * q_di(d3, k + dks + 1, n);
                real gamma_k = gamma_fac_k * q_di(d3, k + dks, n);

                complex vc_kp1 = complex_vcoeff(1 + d3, k + 1, j, i, n);
                complex vc_k = complex_vcoeff(1 + d3, k, j, i, n);

                tri_u(k, j, i, n) += -alpha_kp1 * beta_kp1 * vc_kp1 * gamma_kp2;
                tri_d(k, j, i, n) +=
                    alpha_kp1 * (beta_kp1 * vc_kp1 + beta_k * vc_k) * gamma_kp1;
                tri_l(k, j, i, n) += -alpha_kp1 * beta_k * vc_k * gamma_k;
              }
            }
          }
        });
  }

  virtual void solve(real dt, FieldSet<nprognostic> &rhs,
                     FieldSet<nconstant> &const_vars,
                     FieldSet<nauxiliary> &auxiliary_vars,
                     FieldSet<nprognostic> &solution) override {
    yakl::timer_start("linear_solve");

    const auto &refstate = this->equations->reference_state;

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    // fourier

    auto sol_dens = solution.fields_arr[DENSVAR].data;
    auto sol_v = solution.fields_arr[VVAR].data;
    auto sol_w = solution.fields_arr[WVAR].data;

    auto rhs_v = rhs.fields_arr[VVAR].data;
    auto rhs_w = rhs.fields_arr[WVAR].data;
    auto rhs_dens = rhs.fields_arr[DENSVAR].data;

    auto rhs_dens0 = auxiliary_vars.fields_arr[DENS0VAR].data;
    auto bvar = auxiliary_vars.fields_arr[BVAR].data;
    auto wtendvar = auxiliary_vars.fields_arr[BVAR].data;
    auto uvar = auxiliary_vars.fields_arr[UVAR].data;
    auto uwvar = auxiliary_vars.fields_arr[UWVAR].data;
    auto fvar = auxiliary_vars.fields_arr[FVAR].data;
    auto fwvar = auxiliary_vars.fields_arr[FWVAR].data;

    auto n_cells_x = dual_topology.n_cells_x;
    auto n_cells_y = dual_topology.n_cells_y;
    auto nens = dual_topology.nens;
    auto nl = primal_topology.nl;
    auto ni = primal_topology.ni;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    real dtf = dt / 2;
    real dtf2 = dt * dt / 4;

    const auto &rho_pi = refstate.rho_pi.data;
    const auto &q_pi = refstate.q_pi.data;
    const auto &rho_di = refstate.rho_di.data;
    const auto &q_di = refstate.q_di.data;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(tri_l, this->tri_l);
    YAKL_SCOPE(tri_d, this->tri_d);
    YAKL_SCOPE(tri_u, this->tri_u);
    YAKL_SCOPE(tri_c, this->tri_c);
    YAKL_SCOPE(complex_vcoeff, this->complex_vcoeff);
    YAKL_SCOPE(v_transform, this->v_transform);
    YAKL_SCOPE(w_transform, this->w_transform);
    YAKL_SCOPE(complex_vrhs, this->complex_vrhs);
    YAKL_SCOPE(complex_wrhs, this->complex_wrhs);
    YAKL_SCOPE(Blin_coeff, this->Blin_coeff);

    parallel_for(
        "Prepare rhs 1 - B",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, VS::ndensity_dycore> rhs0;
          compute_Hn1bar<VS::ndensity_dycore, diff_ord, vert_diff_ord>(
              rhs0, rhs_dens, primal_geometry, dual_geometry, pis, pjs, pks, i,
              j, k, n);

          for (int d1 = 0; d1 < VS::ndensity_dycore; ++d1) {
            real b_d1 = 0;
            for (int d2 = 0; d2 < VS::ndensity_dycore; ++d2) {
              b_d1 -= dtf * Blin_coeff(d1, d2, k, n) * rhs0(d2);
            }
            bvar(d1, pks + k, pjs + j, pis + i, n) = b_d1;
          }
        });

    auxiliary_vars.exchange({BVAR});

    parallel_for(
        "Prepare rhs 2 - compute v/w transform",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> mod_v;
          compute_wD0<VS::ndensity_dycore>(mod_v, q_pi, bvar, pis, pjs, pks, i,
                                           j, k, n);
          v_transform(k, j, i, n) =
              rhs_v(0, k + pks, j + pjs, i + pis, n) + mod_v(0);
          if (k < primal_topology.nl) {
            real mod_w = compute_wD0_vert<VS::ndensity_dycore>(
                q_di, bvar, pis, pjs, pks, i, j, k, n);
            w_transform(k, j, i, n) =
                rhs_w(0, k + pks, j + pjs, i + pis, n) + mod_w;
          }
        });

    yakl::timer_start("ffts");
    fftv_x.forward_real(v_transform);
    fftw_x.forward_real(w_transform);
    // fftv_y.forward_real(v_transform);
    // fftw_y.forward_real(w_transform);
    yakl::timer_stop("ffts");

    parallel_for(
        "Transform result to complex",
        yakl::c::Bounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                           {0, nxf - 1, 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real v_real = v_transform(k, j, i, n);
          real v_imag = v_transform(k, j, i + 1, n);
          complex_vrhs(k, j, i / 2, n) = complex(v_real, v_imag);
          if (k < primal_topology.nl) {
            real w_real = w_transform(k, j, i, n);
            real w_imag = w_transform(k, j, i + 1, n);
            complex_wrhs(k, j, i / 2, n) = complex(w_real, w_imag);
          }
        });

    parallel_for(
        "Modify wrhs",
        Bounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                  {0, (nxf - 1) / 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          int ik = i / 2;
          int jk = j / 2;

          complex vc0_kp1 =
              complex_vcoeff(0, k + 1, j, i, n) * complex_vrhs(k + 1, j, i, n);
          complex vc0_k =
              complex_vcoeff(0, k, j, i, n) * complex_vrhs(k, j, i, n);

          real fH2bar_k = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fH2bar_kp1 = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          SArray<real, 1, ndims> fH1_kp1_a;
          SArray<real, 1, ndims> fH1_k_a;
          fourier_H10<diff_ord>(fH1_kp1_a, primal_geometry, dual_geometry, pis,
                                pjs, pks, i, j, k + 1, 0, n_cells_x, n_cells_y,
                                dual_topology.ni);
          fourier_H10<diff_ord>(fH1_k_a, primal_geometry, dual_geometry, pis,
                                pjs, pks, i, j, k, 0, n_cells_x, n_cells_y,
                                dual_topology.ni);
          real fH1h_kp1 = fH1_kp1_a(0);
          real fH1h_k = fH1_k_a(0);

          complex fDnm1bar_kp1 = fourier_Dnm1bar(1, i, j, k + 1, n_cells_x,
                                                 n_cells_y, dual_topology.ni);
          complex fDnm1bar_k = fourier_Dnm1bar(1, i, j, k, n_cells_x, n_cells_y,
                                               dual_topology.ni);

          real he_kp1 = rho_pi(0, k + pks + 1, n);
          real he_k = rho_pi(0, k + pks, n);

          for (int d1 = 0; d1 < VS::ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < VS::ndensity_dycore; ++d2) {
              real alpha_kp1 = dtf2 * q_di(d1, k + dks + 1, n);
              complex beta_kp1 = fH2bar_kp1 * Blin_coeff(d1, d2, k + 1, n) *
                                 q_pi(d2, k + pks + 1, n) * fDnm1bar_kp1 *
                                 he_kp1 * fH1h_kp1;
              complex beta_k = fH2bar_k * Blin_coeff(d1, d2, k, n) *
                               q_pi(d2, k + pks, n) * fDnm1bar_k * he_k *
                               fH1h_k;
              complex_wrhs(k, j, i, n) +=
                  alpha_kp1 * (beta_kp1 * vc0_kp1 - beta_k * vc0_k);
            }
          }
        });

    parallel_for(
        "Tridiagonal solve",
        Bounds<3>(primal_topology.n_cells_y, {0, (nxf - 1) / 2},
                  primal_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          int nz = primal_topology.nl;

          tri_c(0, j, i, n) = tri_u(0, j, i, n) / tri_d(0, j, i, n);
          for (int k = 1; k < nz - 1; ++k) {
            tri_c(k, j, i, n) =
                tri_u(k, j, i, n) /
                (tri_d(k, j, i, n) - tri_l(k, j, i, n) * tri_c(k - 1, j, i, n));
          }
          complex_wrhs(0, j, i, n) /= tri_d(0, j, i, n);
          for (int k = 1; k < nz; ++k) {
            complex_wrhs(k, j, i, n) =
                (complex_wrhs(k, j, i, n) -
                 tri_l(k, j, i, n) * complex_wrhs(k - 1, j, i, n)) /
                (tri_d(k, j, i, n) - tri_l(k, j, i, n) * tri_c(k - 1, j, i, n));
          }
          for (int k = nz - 2; k >= 0; --k) {
            complex_wrhs(k, j, i, n) -=
                tri_c(k, j, i, n) * complex_wrhs(k + 1, j, i, n);
          }
        });

    parallel_for(
        "Compute vhat",
        Bounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                  {0, (nxf - 1) / 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          complex w_kp1;
          if (k < primal_topology.ni - 1) {
            w_kp1 = complex_wrhs(k, j, i, n);
          } else {
            w_kp1 = 0;
          }
          complex w_k;
          if (k > 0) {
            w_k = complex_wrhs(k - 1, j, i, n);
          } else {
            w_k = 0;
          }

          real gamma_fac_kp1 = rho_di(0, k + dks + 1, n) *
                               H01_diagonal(primal_geometry, dual_geometry, pis,
                                            pjs, pks, i, j, k + 1, n);
          real gamma_fac_k = rho_di(0, k + dks, n) *
                             H01_diagonal(primal_geometry, dual_geometry, pis,
                                          pjs, pks, i, j, k, n);

          complex_vrhs(k, j, i, n) *= complex_vcoeff(0, k, j, i, n);
          for (int d1 = 0; d1 < VS::ndensity_dycore; ++d1) {
            complex_vrhs(k, j, i, n) +=
                complex_vcoeff(1 + d1, k, j, i, n) *
                (gamma_fac_kp1 * q_di(d1, k + dks + 1, n) * w_kp1 -
                 gamma_fac_k * q_di(d1, k + dks, n) * w_k);
          }
        });

    parallel_for(
        "Complex solution to real",
        yakl::c::Bounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                           {0, nxf - 1, 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          v_transform(k, j, i, n) = complex_vrhs(k, j, i / 2, n).real();
          v_transform(k, j, i + 1, n) = complex_vrhs(k, j, i / 2, n).imag();
          if (k < primal_topology.nl) {
            w_transform(k, j, i, n) = complex_wrhs(k, j, i / 2, n).real();
            w_transform(k, j, i + 1, n) = complex_wrhs(k, j, i / 2, n).imag();
          }
        });

    yakl::timer_start("ffts");
    fftv_x.inverse_real(v_transform);
    fftw_x.inverse_real(w_transform);
    // fftv_y.inverse_real(v_transform);
    // fftw_y.inverse_real(w_transform);
    yakl::timer_stop("ffts");

    parallel_for(
        "Store v/w solution into array with halos",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          sol_v(0, k + pks, j + pjs, i + pis, n) = v_transform(k, j, i, n);
          if (k < primal_topology.nl) {
            sol_w(0, k + pks, j + pjs, i + pis, n) = w_transform(k, j, i, n);
          }
        });

    solution.exchange({VVAR, WVAR});

    // recover densities

    YAKL_SCOPE(Hk, this->equations->Hk);
    parallel_for(
        "Recover densities 1 - F/Fw",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> u;
          compute_H10<1, diff_ord>(u, sol_v, primal_geometry, dual_geometry,
                                   dis, djs, dks, i, j, k, n);

          fvar(0, k + dks, j + djs, i + dis, n) =
              u(0, 0) * rho_pi(0, k + pks, n);

          if (k < dual_topology.ni - 2) {
            SArray<real, 1, 1> uw;
            compute_H01(uw, sol_w, primal_geometry, dual_geometry, dis, djs,
                        dks, i, j, k + 1, n);
            fwvar(0, k + 1 + dks, j + djs, i + dis, n) =
                uw(0) * rho_di(0, k + dks + 1, n);
          }
        });

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({FVAR, FWVAR});

    parallel_for(
        "Recover densities 2 - Dnm1bar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDnm1bar<VS::ndensity_prognostic>(sol_dens, q_pi, fvar, dis,
                                                    djs, dks, i, j, k, n);
          compute_wDnm1bar_vert<VS::ndensity_prognostic, ADD_MODE::ADD>(
              sol_dens, q_di, fwvar, dis, djs, dks, i, j, k, n);
          for (int d = 0; d < VS::ndensity_prognostic; ++d) {
            sol_dens(d, k + dks, j + djs, i + dis, n) *= -dt / 2;
            sol_dens(d, k + dks, j + djs, i + dis, n) +=
                rhs_dens(d, pks + k, pjs + j, pis + i, n);
          }
        });

    yakl::timer_stop("linear_solve");
  }
};

struct PressureLinearSystem : public LinearSystem {
  yakl::RealFFT1D<real> fftp_x;
  yakl::RealFFT1D<real> fftp_y;

  bool is_initialized = false;

  int nxf, nyf;
  real4d p_transform;

  real4d tri_l;
  real4d tri_d;
  real4d tri_u;
  real4d tri_c;

  using VS = VariableSet;

public:
  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  Equations &equations) override {

    LinearSystem::initialize(params, primal_geom, dual_geom, equations);

    const auto &primal_topology = primal_geom.topology;

    auto pni = primal_topology.ni;
    auto pnl = primal_topology.nl;
    auto nx = primal_topology.n_cells_x;
    auto ny = primal_topology.n_cells_y;
    auto nens = primal_topology.nens;

    this->nxf = nx + 2 - nx % 2;
    this->nyf = ndims > 1 ? ny + 2 - ny % 2 : ny;

    p_transform = real4d("p transform", pni, nyf, nxf, nens);
    yakl::memset(p_transform, 0);

    fftp_x.init(p_transform, 2, nx);
    if (ndims > 1) {
      fftp_y.init(p_transform, 1, ny);
    }

    tri_d = real4d("tri d", pni, nyf, nxf, nens);
    tri_l = real4d("tri l", pni, nyf, nxf, nens);
    tri_u = real4d("tri u", pni, nyf, nxf, nens);
    tri_c = real4d("tri c", pni, nyf, nxf, nens);

    this->is_initialized = true;
  }

  virtual void prepare_pressure_rhs(real dt, FieldSet<nprognostic> &rhs,
                                    FieldSet<nconstant> &const_vars,
                                    FieldSet<nauxiliary> &auxiliary_vars) = 0;

  virtual void solve_for_pressure(real dt, FieldSet<nprognostic> &rhs,
                                  FieldSet<nconstant> &const_vars,
                                  FieldSet<nauxiliary> &auxiliary_vars) = 0;

  virtual void update_velocity(real dt, FieldSet<nprognostic> &rhs,
                               FieldSet<nconstant> &const_vars,
                               FieldSet<nauxiliary> &auxiliary_vars,
                               FieldSet<nprognostic> &solution) = 0;

  virtual void update_densities(real dt, FieldSet<nprognostic> &rhs,
                                FieldSet<nconstant> &const_vars,
                                FieldSet<nauxiliary> &auxiliary_vars,
                                FieldSet<nprognostic> &solution) = 0;

  void solve(real dt, FieldSet<nprognostic> &rhs,
             FieldSet<nconstant> &const_vars,
             FieldSet<nauxiliary> &auxiliary_vars,
             FieldSet<nprognostic> &solution) override {

    prepare_pressure_rhs(dt, rhs, const_vars, auxiliary_vars);
    solve_for_pressure(dt, rhs, const_vars, auxiliary_vars);
    update_velocity(dt, rhs, const_vars, auxiliary_vars, solution);
    update_densities(dt, rhs, const_vars, auxiliary_vars, solution);
  }
};

struct AnelasticLinearSystem : PressureLinearSystem {
  int kfix;

  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  Equations &equations) override {
    PressureLinearSystem::initialize(params, primal_geom, dual_geom, equations);
    this->kfix = primal_geom.topology.ni / 2;
  }

  void compute_coefficients(real dt) override {
    const auto &refstate = this->equations->reference_state;

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    auto n_cells_x = dual_topology.n_cells_x;
    auto n_cells_y = dual_topology.n_cells_y;
    auto nens = dual_topology.nens;
    auto nl = primal_topology.nl;
    auto ni = primal_topology.ni;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    const auto &rho_pi = refstate.rho_pi.data;
    const auto &rho_di = refstate.rho_di.data;

    YAKL_SCOPE(thermo, this->equations->thermo);
    YAKL_SCOPE(grav, this->equations->Hs.g);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(tri_l, this->tri_l);
    YAKL_SCOPE(tri_d, this->tri_d);
    YAKL_SCOPE(tri_u, this->tri_u);
    YAKL_SCOPE(kfix, this->kfix);

    parallel_for(
        "Anelastic set coeffs",
        SimpleBounds<4>(primal_topology.ni, nyf, nxf, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          int ik = i / 2;
          int jk = j / 2;

          SArray<real, 1, ndims> fH1;
          fourier_H10<diff_ord>(fH1, primal_geometry, dual_geometry, pis, pjs,
                                pks, ik, jk, k, 0, dual_topology.n_cells_x,
                                dual_topology.n_cells_y, dual_topology.ni);

          SArray<real, 1, ndims> fD0Dbar;
          fourier_cwD0Dnm1bar(fD0Dbar, 1, ik, jk, k, dual_topology.n_cells_x,
                              dual_topology.n_cells_y, dual_topology.ni);

          tri_l(k, j, i, n) = 0;
          tri_d(k, j, i, n) = 0;
          for (int d = 0; d < ndims; ++d) {
            tri_d(k, j, i, n) += fH1(d) * fD0Dbar(d) * rho_pi(0, k + pks, n);
          }
          tri_u(k, j, i, n) = 0;

          const real h_k = rho_di(0, k + dks, n) *
                           H01_diagonal(primal_geometry, dual_geometry, pis,
                                        pjs, pks, i, j, k, n);
          const real h_kp1 = rho_di(0, k + dks + 1, n) *
                             H01_diagonal(primal_geometry, dual_geometry, pis,
                                          pjs, pks, i, j, k + 1, n);

          tri_u(k, j, i, n) += h_kp1;
          tri_l(k, j, i, n) += h_k;

          if (k == 0) {
            tri_d(k, j, i, n) += -h_kp1;
          } else if (k == (primal_topology.ni - 1)) {
            tri_d(k, j, i, n) += -h_k;
          } else {
            tri_d(k, j, i, n) += -(h_kp1 + h_k);
          }

          // the tridiagonal system that we need to solve is formally singular
          // because of Neumann conditons on both boundaries. To avoid issues
          // with direct solve in the vertical we fix the horizontal mean of
          // pressure at one vertical level
          if (ik == 0 && jk == 0 && k == kfix) {
            tri_d(k, j, i, n) = 1;
            tri_u(k, j, i, n) = 0;
            tri_l(k, j, i, n) = 0;
          }
        });
  }

  void prepare_pressure_rhs(real dt, FieldSet<nprognostic> &rhs,
                            FieldSet<nconstant> &const_vars,
                            FieldSet<nauxiliary> &auxiliary_vars) override {

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const auto &refstate = this->equations->reference_state;
    const auto &rho_pi = refstate.rho_pi.data;
    const auto &rho_di = refstate.rho_di.data;

    const auto Fvar = auxiliary_vars.fields_arr[FVAR].data;
    const auto FWvar = auxiliary_vars.fields_arr[FWVAR].data;
    const auto mfvar = auxiliary_vars.fields_arr[MFVAR].data;
    const auto rhs_v = rhs.fields_arr[VVAR].data;
    const auto rhs_w = rhs.fields_arr[WVAR].data;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    parallel_for(
        "Anelastic rhs 1",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> u;
          compute_H10<1, diff_ord>(u, rhs_v, primal_geometry, dual_geometry,
                                   dis, djs, dks, i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            Fvar(d, k + dks, j + djs, i + dis, n) =
                u(0, d) * rho_pi(0, k + pks, n);
          }

          if (k < dual_topology.ni - 2) {
            SArray<real, 1, 1> uw;
            compute_H01(uw, rhs_w, primal_geometry, dual_geometry, dis, djs,
                        dks, i, j, k + 1, n);
            FWvar(0, k + 1 + dks, j + djs, i + dis, n) =
                uw(0) * rho_di(0, k + dks + 1, n);
          }
        });
    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({FVAR, FWVAR});

    parallel_for(
        "Anelastic rhs 2",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Dnm1bar<1>(mfvar, Fvar, dis, djs, dks, i, j, k, n);
          compute_Dnm1bar_vert<1, ADD_MODE::ADD>(mfvar, FWvar, dis, djs, dks, i,
                                                 j, k, n);
        });
  }

  void solve_for_pressure(real dt, FieldSet<nprognostic> &rhs,
                          FieldSet<nconstant> &const_vars,
                          FieldSet<nauxiliary> &auxiliary_vars) override {

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(p_transform, this->p_transform);
    YAKL_SCOPE(tri_l, this->tri_l);
    YAKL_SCOPE(tri_d, this->tri_d);
    YAKL_SCOPE(tri_u, this->tri_u);
    YAKL_SCOPE(tri_c, this->tri_c);
    YAKL_SCOPE(kfix, this->kfix);

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const auto mfvar = auxiliary_vars.fields_arr[MFVAR].data;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    parallel_for(
        "Anelastic pressure solve 1",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          p_transform(k, j, i, n) = -mfvar(0, k + pks, j + pjs, i + pis, n);
        });

    yakl::timer_start("ffts");
    fftp_x.forward_real(p_transform);
    if (ndims > 1) {
      fftp_y.forward_real(p_transform);
    }
    yakl::timer_stop("ffts");

    parallel_for(
        "Anelastic tridiagonal solve",
        Bounds<3>(nyf, nxf, primal_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          // set the horizontal mean of pressure to zero at k = kfix
          int ik = i / 2;
          int jk = j / 2;
          if (ik == 0 && jk == 0) {
            p_transform(kfix, j, i, n) = 0;
          }
          int nz = primal_topology.ni;
          tri_c(0, j, i, n) = tri_u(0, j, i, n) / tri_d(0, j, i, n);
          for (int k = 1; k < nz - 1; ++k) {
            tri_c(k, j, i, n) =
                tri_u(k, j, i, n) /
                (tri_d(k, j, i, n) - tri_l(k, j, i, n) * tri_c(k - 1, j, i, n));
          }
          p_transform(0, j, i, n) /= tri_d(0, j, i, n);
          for (int k = 1; k < nz; ++k) {
            p_transform(k, j, i, n) =
                (p_transform(k, j, i, n) -
                 tri_l(k, j, i, n) * p_transform(k - 1, j, i, n)) /
                (tri_d(k, j, i, n) - tri_l(k, j, i, n) * tri_c(k - 1, j, i, n));
          }
          for (int k = nz - 2; k >= 0; --k) {
            p_transform(k, j, i, n) -=
                tri_c(k, j, i, n) * p_transform(k + 1, j, i, n);
          }
        });

    yakl::timer_start("ffts");
    fftp_x.inverse_real(p_transform);
    if (ndims > 1) {
      fftp_y.inverse_real(p_transform);
    }
    yakl::timer_stop("ffts");
  }

  void update_velocity(real dt, FieldSet<nprognostic> &rhs,
                       FieldSet<nconstant> &const_vars,
                       FieldSet<nauxiliary> &auxiliary_vars,
                       FieldSet<nprognostic> &solution) override {

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(p_transform, this->p_transform);
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const auto Bvar = auxiliary_vars.fields_arr[BVAR].data;
    const auto sol_v = solution.fields_arr[VVAR].data;
    const auto sol_w = solution.fields_arr[WVAR].data;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    // need to store p into a field with halos
    parallel_for(
        "Anelastic update v 1",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Bvar(0, k + pks, j + pjs, i + pis, n) = p_transform(k, j, i, n);
        });
    auxiliary_vars.exchange({BVAR});

    parallel_for(
        "Anelastic update v 2",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D0_vert<1, ADD_MODE::ADD>(sol_w, Bvar, pis, pjs, pks, i, j, k,
                                            n);
        });

    parallel_for(
        "Anelastic update v 3",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D0<1, ADD_MODE::ADD>(sol_v, Bvar, pis, pjs, pks, i, j, k, n);
        });
  }

  void update_densities(real dt, FieldSet<nprognostic> &rhs,
                        FieldSet<nconstant> &const_vars,
                        FieldSet<nauxiliary> &auxiliary_vars,
                        FieldSet<nprognostic> &solution) override {}
};

struct CompressiblePressureLinearSystem : PressureLinearSystem {
  real3d linp_coeff;

  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  Equations &equations) override {
    PressureLinearSystem::initialize(params, primal_geom, dual_geom, equations);
    linp_coeff =
        real3d("linp coeff", VS::ndensity_active, primal_geometry.topology.ni,
               primal_geometry.topology.nens);
  }

  void compute_coefficients(real dt) override {
    const auto &refstate = this->equations->reference_state;
    const auto &varset = this->equations->varset;

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    auto n_cells_x = dual_topology.n_cells_x;
    auto n_cells_y = dual_topology.n_cells_y;
    auto nens = dual_topology.nens;
    auto nl = primal_topology.nl;
    auto ni = primal_topology.ni;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    real alpha = dt / 2;

    const auto &rho_pi = refstate.rho_pi.data;
    const auto &q_pi = refstate.q_pi.data;
    const auto &rho_di = refstate.rho_di.data;
    const auto &q_di = refstate.q_di.data;
    const auto &pres_pi = refstate.pres_pi.data;
    const auto &pres_di = refstate.pres_di.data;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(tri_l, this->tri_l);
    YAKL_SCOPE(tri_d, this->tri_d);
    YAKL_SCOPE(tri_u, this->tri_u);
    YAKL_SCOPE(linp_coeff, this->linp_coeff);

    parallel_for(
        "pres coeffs",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          SArray<real, 1, VS::ndensity_active> p_coeff;
          varset.linear_pressure_coeffs(p_coeff, k, pks, n);
          for (int d = 0; d < VS::ndensity_active; ++d) {
            linp_coeff(d, k, n) = p_coeff(d);
          }
        });

    parallel_for(
        "Compressible linear system coefficients",
        SimpleBounds<4>(primal_topology.ni, nyf, nxf, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          int ik = i / 2;
          int jk = j / 2;

          SArray<real, 1, ndims> fH1;
          fourier_H10<diff_ord>(fH1, primal_geometry, dual_geometry, pis, pjs,
                                pks, ik, jk, k, 0, dual_topology.n_cells_x,
                                dual_topology.n_cells_y, dual_topology.ni);

          SArray<real, 1, ndims> fD0Dbar;
          fourier_cwD0Dnm1bar(fD0Dbar, 1, ik, jk, k, dual_topology.n_cells_x,
                              dual_topology.n_cells_y, dual_topology.ni);

          real fHn1bar = fourier_Hn1bar<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, ik, jk, k, 0,
              dual_topology.n_cells_x, dual_topology.n_cells_y,
              dual_topology.nl);

          tri_d(k, j, i, n) = 1;
          tri_l(k, j, i, n) = 0;
          tri_u(k, j, i, n) = 0;

          for (int dd = 0; dd < VS::ndensity_active; ++dd) {
            for (int d = 0; d < ndims; ++d) {
              tri_d(k, j, i, n) -= alpha * alpha * linp_coeff(dd, k, n) *
                                   fHn1bar * fH1(d) * fD0Dbar(d) *
                                   q_pi(dd, k + pks, n);
            }
          }

          for (int d = 0; d < VS::ndensity_active; ++d) {
            const real gamma_kp1 = 1;
            const real gamma_k = 1;
            const real gamma_km1 = 1;

            const real beta_kp1 = q_di(d, k + 1 + dks, n) *
                                  H01_diagonal(primal_geometry, dual_geometry,
                                               pis, pjs, pks, i, j, k + 1, n);

            const real beta_k = q_di(d, k + dks, n) *
                                H01_diagonal(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k, n);

            // TODO: This is more tricky with higher order hodge stars !!!
            const real alpha_k = -alpha * alpha * fHn1bar * linp_coeff(d, k, n);

            tri_u(k, j, i, n) += alpha_k * beta_kp1 * gamma_kp1;
            tri_l(k, j, i, n) += alpha_k * beta_k * gamma_km1;

            if (k == 0) {
              tri_d(k, j, i, n) += -alpha_k * beta_kp1 * gamma_k;
            } else if (k == (primal_topology.ni - 1)) {
              tri_d(k, j, i, n) += -alpha_k * beta_k * gamma_k;
            } else {
              tri_d(k, j, i, n) += -alpha_k * (beta_kp1 + beta_k) * gamma_k;
            }
          }
        });
  }

  void prepare_pressure_rhs(real dt, FieldSet<nprognostic> &rhs,
                            FieldSet<nconstant> &const_vars,
                            FieldSet<nauxiliary> &auxiliary_vars) override {

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(p_transform, this->p_transform);
    YAKL_SCOPE(linp_coeff, this->linp_coeff);
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const auto &refstate = this->equations->reference_state;
    const auto &rho_pi = refstate.rho_pi.data;
    const auto &rho_di = refstate.rho_di.data;
    const auto &q_pi = refstate.q_pi.data;
    const auto &q_di = refstate.q_di.data;

    const auto Fvar = auxiliary_vars.fields_arr[FVAR].data;
    const auto FWvar = auxiliary_vars.fields_arr[FWVAR].data;
    const auto Bvar = auxiliary_vars.fields_arr[BVAR].data;
    const auto mfvar = auxiliary_vars.fields_arr[MFVAR].data;
    const auto rhs_v = rhs.fields_arr[VVAR].data;
    const auto rhs_w = rhs.fields_arr[WVAR].data;
    const auto rhs_dens = rhs.fields_arr[DENSVAR].data;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    parallel_for(
        "Linear solve rhs 1",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> u;
          compute_H10<1, diff_ord>(u, rhs_v, primal_geometry, dual_geometry,
                                   dis, djs, dks, i, j, k, n);

          for (int d = 0; d < ndims; ++d) {
            Fvar(d, k + dks, j + djs, i + dis, n) =
                u(0, d) * rho_pi(0, k + pks, n);
          }

          if (k < dual_topology.ni - 2) {
            SArray<real, 1, 1> uw;
            compute_H01(uw, rhs_w, primal_geometry, dual_geometry, dis, djs,
                        dks, i, j, k + 1, n);
            FWvar(0, k + 1 + dks, j + djs, i + dis, n) =
                uw(0) * rho_di(0, k + dks + 1, n);
          }
        });
    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({FVAR, FWVAR});

    parallel_for(
        "Linear solve rhs 2",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, VS::ndensity_active> tend;
          SArray<real, 1, VS::ndensity_active> vtend;

          compute_wDnm1bar<VS::ndensity_active>(tend, q_pi, Fvar, dis, djs, dks,
                                                i, j, k, n);
          compute_wDnm1bar_vert<VS::ndensity_active>(vtend, q_di, FWvar, dis,
                                                     djs, dks, i, j, k, n);

          for (int d = 0; d < VS::ndensity_active; ++d) {
            mfvar(d, k + dks, j + djs, i + dis, n) =
                rhs_dens(d, k + dks, j + djs, i + dis, n);
            mfvar(d, k + dks, j + djs, i + dis, n) -= 0.5_fp * dt * tend(d);
            mfvar(d, k + dks, j + djs, i + dis, n) -= 0.5_fp * dt * vtend(d);
          }
        });
    auxiliary_vars.exchange({MFVAR});

    parallel_for(
        "Linear solve rhs 3",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hn1bar<VS::ndensity_active, diff_ord, vert_diff_ord>(
              Bvar, mfvar, primal_geometry, dual_geometry, pis, pjs, pks, i, j,
              k, n);
        });

    parallel_for(
        "Linear solve rhs 4",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          p_transform(k, j, i, n) = 0;
          for (int d = 0; d < VS::ndensity_active; ++d) {
            p_transform(k, j, i, n) +=
                linp_coeff(d, k, n) * Bvar(d, k + pks, j + pjs, i + pis, n);
          }
        });
  }

  void solve_for_pressure(real dt, FieldSet<nprognostic> &rhs,
                          FieldSet<nconstant> &const_vars,
                          FieldSet<nauxiliary> &auxiliary_vars) override {

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(p_transform, this->p_transform);
    YAKL_SCOPE(tri_c, this->tri_c);
    YAKL_SCOPE(tri_l, this->tri_l);
    YAKL_SCOPE(tri_u, this->tri_u);
    YAKL_SCOPE(tri_d, this->tri_d);
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const auto mfvar = auxiliary_vars.fields_arr[MFVAR].data;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    yakl::timer_start("ffts");
    fftp_x.forward_real(p_transform);
    if (ndims > 1) {
      fftp_y.forward_real(p_transform);
    }
    yakl::timer_stop("ffts");

    parallel_for(
        "Anelastic tridiagonal solve",
        Bounds<3>(nyf, nxf, primal_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          int ik = i / 2;
          int jk = j / 2;

          int nz = primal_topology.ni;
          tri_c(0, j, i, n) = tri_u(0, j, i, n) / tri_d(0, j, i, n);
          for (int k = 1; k < nz - 1; ++k) {
            tri_c(k, j, i, n) =
                tri_u(k, j, i, n) /
                (tri_d(k, j, i, n) - tri_l(k, j, i, n) * tri_c(k - 1, j, i, n));
          }
          p_transform(0, j, i, n) /= tri_d(0, j, i, n);
          for (int k = 1; k < nz; ++k) {
            p_transform(k, j, i, n) =
                (p_transform(k, j, i, n) -
                 tri_l(k, j, i, n) * p_transform(k - 1, j, i, n)) /
                (tri_d(k, j, i, n) - tri_l(k, j, i, n) * tri_c(k - 1, j, i, n));
          }
          for (int k = nz - 2; k >= 0; --k) {
            p_transform(k, j, i, n) -=
                tri_c(k, j, i, n) * p_transform(k + 1, j, i, n);
          }
        });

    yakl::timer_start("ffts");
    fftp_x.inverse_real(p_transform);
    if (ndims > 1) {
      fftp_y.inverse_real(p_transform);
    }
    yakl::timer_stop("ffts");
  }

  void update_velocity(real dt, FieldSet<nprognostic> &rhs,
                       FieldSet<nconstant> &const_vars,
                       FieldSet<nauxiliary> &auxiliary_vars,
                       FieldSet<nprognostic> &solution) override {

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(p_transform, this->p_transform);
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const auto &refstate = this->equations->reference_state;
    const auto &rho_pi = refstate.rho_pi.data;
    const auto &rho_di = refstate.rho_di.data;

    const auto Bvar = auxiliary_vars.fields_arr[BVAR].data;
    const auto rhs_v = rhs.fields_arr[VVAR].data;
    const auto rhs_w = rhs.fields_arr[WVAR].data;
    const auto sol_v = solution.fields_arr[VVAR].data;
    const auto sol_w = solution.fields_arr[WVAR].data;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    parallel_for(
        "Anelastic - store p",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Bvar(0, k + pks, j + pjs, i + pis, n) = p_transform(k, j, i, n);
        });
    auxiliary_vars.exchange({BVAR});

    parallel_for(
        "Anelastic - add pressure gradient W",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, 1> dpdz;
          compute_D0_vert<1>(dpdz, Bvar, pis, pjs, pks, i, j, k, n);
          real rho = rho_di(0, k + 1 + dks, n);
          sol_w(0, k + pks, j + pjs, i + pis, n) =
              rhs_w(0, k + pks, j + pjs, i + pis, n) -
              0.5_fp * dt * dpdz(0) / rho;
        });

    parallel_for(
        "Anelastic - add pressure gradient V",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> dpdh;
          compute_D0<1>(dpdh, Bvar, pis, pjs, pks, i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            real rho = rho_pi(0, k + pks, n);
            sol_v(d, k + pks, j + pjs, i + pis, n) =
                rhs_v(d, k + pks, j + pjs, i + pis, n) -
                0.5_fp * dt * dpdh(0, d) / rho;
          }
        });
  }

  void update_densities(real dt, FieldSet<nprognostic> &rhs,
                        FieldSet<nconstant> &const_vars,
                        FieldSet<nauxiliary> &auxiliary_vars,
                        FieldSet<nprognostic> &solution) override {

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const auto &refstate = this->equations->reference_state;
    const auto &rho_pi = refstate.rho_pi.data;
    const auto &rho_di = refstate.rho_di.data;
    const auto &q_pi = refstate.q_pi.data;
    const auto &q_di = refstate.q_di.data;

    const auto rhs_dens = rhs.fields_arr[DENSVAR].data;
    const auto sol_v = solution.fields_arr[VVAR].data;
    const auto sol_w = solution.fields_arr[WVAR].data;
    const auto sol_dens = solution.fields_arr[DENSVAR].data;
    const auto Fvar = auxiliary_vars.fields_arr[FVAR].data;
    const auto FWvar = auxiliary_vars.fields_arr[FWVAR].data;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    solution.exchange({VVAR});
    solution.exchange({WVAR});

    parallel_for(
        "Recover densities 1 - F/Fw",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> u;
          compute_H10<1, diff_ord>(u, sol_v, primal_geometry, dual_geometry,
                                   dis, djs, dks, i, j, k, n);

          for (int d = 0; d < ndims; ++d) {
            Fvar(d, k + dks, j + djs, i + dis, n) =
                u(0, d) * rho_pi(0, k + pks, n);
          }

          if (k < dual_topology.ni - 2) {
            SArray<real, 1, 1> uw;
            compute_H01(uw, sol_w, primal_geometry, dual_geometry, dis, djs,
                        dks, i, j, k + 1, n);
            FWvar(0, k + 1 + dks, j + djs, i + dis, n) =
                uw(0) * rho_di(0, k + dks + 1, n);
          }
        });

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({FVAR, FWVAR});
    yakl::memset(sol_dens, 0);

    parallel_for(
        "Recover densities 2 - Dnm1bar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDnm1bar<VS::ndensity_prognostic>(sol_dens, q_pi, Fvar, dis,
                                                    djs, dks, i, j, k, n);
          compute_wDnm1bar_vert<VS::ndensity_prognostic, ADD_MODE::ADD>(
              sol_dens, q_di, FWvar, dis, djs, dks, i, j, k, n);
          for (int d = 0; d < VS::ndensity_prognostic; ++d) {
            sol_dens(d, k + dks, j + djs, i + dis, n) *= -dt / 2;
            sol_dens(d, k + dks, j + djs, i + dis, n) +=
                rhs_dens(d, pks + k, pjs + j, pis + i, n);
          }
        });
  }
};

std::unique_ptr<LinearSystem>
model_linear_system(const ModelParameters &params) {
  if (VariableSet::compressible) {
    if (params.linear_system == "velocity") {
      return std::make_unique<CompressibleVelocityLinearSystem>();
    } else if (params.linear_system == "pressure") {
      return std::make_unique<CompressiblePressureLinearSystem>();
      //} else if (params.linear_system == "pressure_gravity") {
      //  return std::make_unique<CompressiblePressureGravityLinearSystem>();
    } else {
      throw std::runtime_error("Unknown linear system choice");
    }
  } else {
    return std::make_unique<AnelasticLinearSystem>();
  }
}

// *******   Statistics   ***********//

class ModelStats : public Stats {
  using VS = VariableSet;

public:
  real3d TEarr, KEarr, PEarr, IEarr, PENSarr, trimmed_density;
  real4d PVarr;

  void initialize(ModelParameters &params, Parallel &par,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom, Equations &equations) {
    Stats::initialize(params, par, primal_geom, dual_geom, equations);
    this->stats_arr[DENSSTAT].initialize("mass", VS::ndensity_prognostic,
                                         this->statsize, this->nens,
                                         this->masterproc);
    this->stats_arr[DENSMAXSTAT].initialize("densmax", VS::ndensity_prognostic,
                                            this->statsize, this->nens,
                                            this->masterproc);
    this->stats_arr[DENSMINSTAT].initialize("densmin", VS::ndensity_prognostic,
                                            this->statsize, this->nens,
                                            this->masterproc);
    this->stats_arr[ESTAT].initialize("energy", 4, this->statsize, this->nens,
                                      this->masterproc);
    this->stats_arr[PVSTAT].initialize("pv", ndims > 1 ? 3 : 1, this->statsize,
                                       this->nens, this->masterproc);
    this->stats_arr[PESTAT].initialize("pens", 1, this->statsize, this->nens,
                                       this->masterproc);

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    TEarr = real3d("TE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    KEarr = real3d("KE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    IEarr = real3d("IE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    PEarr = real3d("PE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    PVarr = real4d("PV", ndims > 1 ? 3 : 1, primal_topology.nl,
                   primal_topology.n_cells_y, primal_topology.n_cells_x);
    PENSarr = real3d("PENS", primal_topology.nl, primal_topology.n_cells_y,
                     primal_topology.n_cells_x);
    trimmed_density = real3d("trimmed_density", dual_topology.nl,
                             dual_topology.n_cells_y, dual_topology.n_cells_x);
  }

  void compute(FieldSet<nprognostic> &progvars, FieldSet<nconstant> &constvars,
               int tind) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    for (int n = 0; n < nens; n++) {

      SArray<real, 1, VS::ndensity_prognostic> masslocal, massglobal;
      SArray<real, 1, VS::ndensity_prognostic> densmaxlocal, densmaxglobal;
      SArray<real, 1, VS::ndensity_prognostic> densminlocal, densminglobal;
      SArray<real, 1, (ndims > 1 ? 3 : 1)> pvlocal, pvglobal;
      SArray<real, 1, 4> elocal, eglobal;
      SArray<real, 1, 1> pelocal, peglobal;

      pvlocal(0) = 0.;
      pvglobal(0) = 0.;
      pelocal(0) = 0.;
      peglobal(0) = 0.;
      elocal(0) = 0.;
      elocal(1) = 0.;
      elocal(2) = 0.;
      elocal(3) = 0.;
      eglobal(0) = 0.;
      eglobal(1) = 0.;
      eglobal(2) = 0.;
      eglobal(3) = 0.;
      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        masslocal(l) = 0.;
        massglobal(l) = 0.;
      }
      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        densmaxlocal(l) = 0.;
        densmaxglobal(l) = 0.;
      }
      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        densminlocal(l) = 0.;
        densminglobal(l) = 0.;
      }

      int dis = dual_topology.is;
      int djs = dual_topology.js;
      int dks = dual_topology.ks;

      YAKL_SCOPE(Hk, equations->Hk);
      YAKL_SCOPE(Hs, equations->Hs);
      parallel_for(
          "Compute energetics stats",
          SimpleBounds<3>(dual_topology.nl - 2, dual_topology.n_cells_y,
                          dual_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int k, int j, int i) {
            real KE, PE, IE;
            KE = Hk.compute_KE(progvars.fields_arr[VVAR].data,
                               progvars.fields_arr[WVAR].data,
                               progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, k + 1, n);
            PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, k + 1, n);
            IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, k + 1, n);
            TEarr(k + 1, j, i) = KE + PE + IE;
            KEarr(k + 1, j, i) = KE;
            PEarr(k + 1, j, i) = PE;
            IEarr(k + 1, j, i) = IE;
          });
      parallel_for(
          "Compute energetics stats bnd",
          SimpleBounds<2>(dual_topology.n_cells_y, dual_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int j, int i) {
            real KE, PE, IE;
            KE = Hk.compute_KE_bottom(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, 0, n);
            PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, 0, n);
            IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, 0, n);
            TEarr(0, j, i) = KE + PE + IE;
            KEarr(0, j, i) = KE;
            PEarr(0, j, i) = PE;
            IEarr(0, j, i) = IE;
            KE = Hk.compute_KE_top(progvars.fields_arr[VVAR].data,
                                   progvars.fields_arr[WVAR].data,
                                   progvars.fields_arr[DENSVAR].data, dis, djs,
                                   dks, i, j, dual_topology.nl - 1, n);
            PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, dual_topology.nl - 1, n);
            IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, dual_topology.nl - 1, n);
            TEarr(dual_topology.nl - 1, j, i) = KE + PE + IE;
            KEarr(dual_topology.nl - 1, j, i) = KE;
            PEarr(dual_topology.nl - 1, j, i) = PE;
            IEarr(dual_topology.nl - 1, j, i) = IE;
          });

      elocal(0) = yakl::intrinsics::sum(TEarr);
      elocal(1) = yakl::intrinsics::sum(KEarr);
      elocal(2) = yakl::intrinsics::sum(PEarr);
      elocal(3) = yakl::intrinsics::sum(IEarr);

      int pis = primal_topology.is;
      int pjs = primal_topology.js;
      int pks = primal_topology.ks;

      YAKL_SCOPE(PVPE, equations->PVPE);
      parallel_for(
          "Compute PV/PE stats",
          SimpleBounds<3>(primal_topology.nl - 2, primal_topology.n_cells_y,
                          primal_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int k, int j, int i) {
            pvpe_extruded vals_pvpe;
            vals_pvpe = PVPE.compute_PVPE(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISHZVAR].data,
                constvars.fields_arr[CORIOLISXYVAR].data, pis, pjs, pks, i, j,
                k + 1, n);

            for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
              PVarr(d, k + 1, j, i) = vals_pvpe.pv(d);
            }
            PENSarr(k + 1, j, i) = vals_pvpe.pe;
          });
      parallel_for(
          "Compute PV/PE stats bnd",
          SimpleBounds<2>(primal_topology.n_cells_y, primal_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int j, int i) {
            pvpe_extruded vals_pvpe;
            vals_pvpe = PVPE.compute_PVPE_bottom(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISHZVAR].data,
                constvars.fields_arr[CORIOLISXYVAR].data, pis, pjs, pks, i, j,
                0, n);
            for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
              PVarr(d, 0, j, i) = vals_pvpe.pv(d);
            }
            PENSarr(0, j, i) = vals_pvpe.pe;
            vals_pvpe = PVPE.compute_PVPE_top(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISHZVAR].data,
                constvars.fields_arr[CORIOLISXYVAR].data, pis, pjs, pks, i, j,
                primal_topology.nl - 1, n);
            for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
              PVarr(d, primal_topology.nl - 1, j, i) = vals_pvpe.pv(d);
            }
            PENSarr(primal_topology.nl - 1, j, i) = vals_pvpe.pe;
          });

      pvlocal(0) = yakl::intrinsics::sum(
          PVarr.slice<3>(0, yakl::COLON, yakl::COLON, yakl::COLON));
      if (ndims > 1) {
        pvlocal(1) = yakl::intrinsics::sum(
            PVarr.slice<3>(1, yakl::COLON, yakl::COLON, yakl::COLON));
        pvlocal(2) = yakl::intrinsics::sum(
            PVarr.slice<3>(2, yakl::COLON, yakl::COLON, yakl::COLON));
      }
      pelocal(0) = yakl::intrinsics::sum(PENSarr);

      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        parallel_for(
            "Compute trimmed density",
            SimpleBounds<3>(dual_topology.nl, dual_topology.n_cells_y,
                            dual_topology.n_cells_x),
            YAKL_CLASS_LAMBDA(int k, int j, int i) {
              trimmed_density(k, j, i) = progvars.fields_arr[DENSVAR].data(
                  l, k + dks, j + djs, i + dis, n);
            });

        masslocal(l) = yakl::intrinsics::sum(trimmed_density);
        densmaxlocal(l) = yakl::intrinsics::maxval(trimmed_density);
        densminlocal(l) = yakl::intrinsics::minval(trimmed_density);
      }

      // MPI sum/min/max
      this->ierr = MPI_Ireduce(&masslocal, &massglobal, VS::ndensity_prognostic,
                               PAMC_MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,
                               &this->Req[DENSSTAT]);
      this->ierr = MPI_Ireduce(&densmaxlocal, &densmaxglobal,
                               VS::ndensity_prognostic, PAMC_MPI_REAL, MPI_MAX,
                               0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
      this->ierr = MPI_Ireduce(&densminlocal, &densminglobal,
                               VS::ndensity_prognostic, PAMC_MPI_REAL, MPI_MIN,
                               0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
      this->ierr =
          MPI_Ireduce(&pvlocal, &pvglobal, ndims > 1 ? 3 : 1, PAMC_MPI_REAL,
                      MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
      this->ierr = MPI_Ireduce(&pelocal, &peglobal, 1, PAMC_MPI_REAL, MPI_SUM,
                               0, MPI_COMM_WORLD, &this->Req[PESTAT]);
      this->ierr = MPI_Ireduce(&elocal, &eglobal, 4, PAMC_MPI_REAL, MPI_SUM, 0,
                               MPI_COMM_WORLD, &this->Req[ESTAT]);

      this->ierr = MPI_Waitall(nstats, this->Req, this->Status);

      if (masterproc) {
        for (int l = 0; l < VS::ndensity_prognostic; l++) {
          this->stats_arr[DENSSTAT].data(l, tind, n) = massglobal(l);
          this->stats_arr[DENSMAXSTAT].data(l, tind, n) = densmaxglobal(l);
          this->stats_arr[DENSMINSTAT].data(l, tind, n) = densminglobal(l);
        }

        this->stats_arr[ESTAT].data(0, tind, n) = eglobal(0);
        this->stats_arr[ESTAT].data(1, tind, n) = eglobal(1);
        this->stats_arr[ESTAT].data(2, tind, n) = eglobal(2);
        this->stats_arr[ESTAT].data(3, tind, n) = eglobal(3);
        for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
          this->stats_arr[PVSTAT].data(d, tind, n) = pvglobal(d);
        }
        this->stats_arr[PESTAT].data(0, tind, n) = peglobal(0);
      }
    }
  }
};

// *******   FieldSet Initialization   ***********//
void initialize_variables(
    const Topology &ptopo, const Topology &dtopo,
    std::array<FieldDescription, nprognostic> &prog_desc_arr,
    std::array<FieldDescription, nconstant> &const_desc_arr,
    std::array<FieldDescription, nauxiliary> &aux_desc_arr) {

  using VS = VariableSet;

  // primal grid represents straight quantities, dual grid twisted quantities
  // ndims is the BASEDIM size!

  // v, w, dens
  prog_desc_arr[VVAR] = {"v", ptopo, 1, 0, 1}; // v = straight (1,0)-form
  prog_desc_arr[WVAR] = {"w", ptopo, 0, 1, 1}; // w = straight (0,1)-form
  prog_desc_arr[DENSVAR] = {
      "dens", dtopo, ndims, 1,
      VS::ndensity_prognostic}; // dens = twisted (n,1)-form

  // hs
  const_desc_arr[HSVAR] = {"hs", dtopo, ndims, 1, 1}; // hs = twisted (n,1)-form
  const_desc_arr[CORIOLISHZVAR] = {"coriolishz", ptopo, 1, 1,
                                   1}; // f = straight (1,1)-form

  // functional derivatives = F, B, K, he, U
  aux_desc_arr[FVAR] = {"F", dtopo, ndims - 1, 1, 1}; // F = twisted
                                                      // (n-1,1)-form

  aux_desc_arr[BVAR] = {"B", ptopo, 0, 0,
                        VS::ndensity_active}; // B = straight (0,0)-form

  aux_desc_arr[KVAR] = {"K", dtopo, ndims, 1, 1}; // K = twisted (n,1)-form

  aux_desc_arr[UVAR] = {"U", dtopo, ndims - 1, 1, 1}; // U = twisted
                                                      // (n-1,1)-form

  aux_desc_arr[FWVAR] = {"Fw", dtopo, ndims, 0, 1}; // Fw = twisted (n,0)-form

  aux_desc_arr[UWVAR] = {"Uw", dtopo, ndims, 0, 1}; // Uw = twisted (n,0)-form

  // second copy on F/Fw, needed for the semi-implicit solver because
  // the symplecitic operator gets evaluated at a different state
  // and needs its own F and Fw
  aux_desc_arr[F2VAR] = {"F2", dtopo, ndims - 1, 1, 1}; // F2 = twisted
                                                        // (n-1,1)-form
  aux_desc_arr[FW2VAR] = {"Fw2", dtopo, ndims, 0,
                          1}; // Fw2 = twisted (n,0)-form

  aux_desc_arr[FTVAR] = {"FT", ptopo, 0, 1, ndims}; // FT
  aux_desc_arr[FTWVAR] = {"FTW", ptopo, 1, 0, 1};   // FTW

  // dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_desc_arr[DENS0VAR] = {
      "dens0", ptopo, 0, 0,
      VS::ndensity_prognostic}; // dens0 = straight (0,0)-form
  aux_desc_arr[DENSEDGERECONVAR] = {
      "densedgerecon", dtopo, ndims, 1,
      2 * ndims * VS::ndensity_prognostic}; // densedgerecon lives on dual
                                            // cells, associated with F
  aux_desc_arr[DENSRECONVAR] = {"densrecon", dtopo, ndims - 1, 1,
                                VS::ndensity}; // densrecon lives on horiz dual
                                               // edges, associated with F
  aux_desc_arr[DENSVERTEDGERECONVAR] = {
      "densvertedgerecon", dtopo, ndims, 1,
      2 * VS::ndensity_prognostic}; // densedgerecon lives on dual cells,
                                    // associated with Fw
  aux_desc_arr[DENSVERTRECONVAR] = {
      "densvertrecon", dtopo, ndims, 0,
      VS::ndensity}; // densvertrecon lives on vert dual
                     // edges, associated with Fw

  // fct stuff- Phi, Mf, edgeflux
  aux_desc_arr[MFVAR] = {"Mf", dtopo, ndims, 1, VS::ndensity_prognostic};
  aux_desc_arr[EDGEFLUXVAR] = {"edgeflux", dtopo, ndims - 1, 1,
                               VS::ndensity_prognostic};
  aux_desc_arr[VERTEDGEFLUXVAR] = {"vertedgeflux", dtopo, ndims, 0,
                                   VS::ndensity_prognostic};

  // Q stuff
  aux_desc_arr[QHZVAR] = {"QHZ", dtopo, ndims - 1, 0,
                          1}; // in 2d QHZ=QXZ00 = twisted (0,0)-form,
                              // in 3d QHZ=(QXZ10,QYZ10) lives on twisted edges
  aux_desc_arr[QHZRECONVAR] = {
      "qhzrecon", ptopo, 0, 1,
      ndims}; // qhzrecon lives on vert primal edges, associated with w
  aux_desc_arr[QHZEDGERECONVAR] = {
      "qhzedgerecon", ptopo, 1, 1,
      2 * 1}; // qhzedgerecon lives on primal cells in 2d or primal faces in 3d,
              // associated with Fw/w
  aux_desc_arr[QHZVERTRECONVAR] = {
      "qhzvertrecon", ptopo, 1, 0,
      1}; // qhzsvertrecon lives on horiz primal edges, associated with v
  aux_desc_arr[QHZVERTEDGERECONVAR] = {
      "qhzvertedgerecon", ptopo, 1, 1,
      2 * 1}; // qhzvertedgerecon lives on primal cells 2d or primal faces in
              // 3d, associated with F/v
  aux_desc_arr[QHZFLUXVAR] = {
      "qhzflux", ptopo, 0, 1,
      1}; // qhzflux lives on vert primal edges, associated with w
  aux_desc_arr[QHZVERTFLUXVAR] = {
      "qhzvertflux", ptopo, 1, 0,
      1}; // qhzvertflux lives on horiz primal edges, associated with v

  aux_desc_arr[FHZVAR] = {"fhz", dtopo, ndims - 1, 0,
                          1}; // in 2d fhz=fxz00 is a twisted 0,0 form
                              // in 3d fhz=(fxz10,fyz10) lives on twisted edges
  aux_desc_arr[CORIOLISHZRECONVAR] = {
      "coriolishzrecon", ptopo, 0, 1,
      ndims}; // coriolishzrecon lives on vert primal edges, associated with w
  aux_desc_arr[CORIOLISHZEDGERECONVAR] = {
      "coriolishzedgerecon", ptopo, 1, 1,
      2 * 1}; // coriolishzedgerecon lives on primal cells in 2d or primal faces
              // in 3d, associated with Fw/w
  aux_desc_arr[CORIOLISHZVERTRECONVAR] = {
      "coriolishzvertrecon", ptopo, 1, 0,
      1}; // coriolishzsvertrecon lives on horiz primal edges, associated with v
  aux_desc_arr[CORIOLISHZVERTEDGERECONVAR] = {
      "coriolishzvertedgerecon", ptopo, 1, 1,
      2 * 1}; // coriolishzvertedgerecon lives on primal cells in 2d or primal
              // faces in 3d, associated with F/v

  // diffusion variables
  aux_desc_arr[FDIFFVAR] = {
      "F", dtopo, ndims - 1, 1,
      VS::ndensity_diffused}; // diffusive F = twisted (n-1,1)-form
  aux_desc_arr[FWDIFFVAR] = {
      "Fw", dtopo, ndims, 0,
      VS::ndensity_diffused}; // diffusive Fw = twisted (n,0)-form

  // additonal Q and coriolis variables for 3d
  if (ndims > 1) {
    aux_desc_arr[QXYVAR] = {"QXY", dtopo, 0, ndims - 1,
                            1}; // QXY lives on twisted edges
    aux_desc_arr[QXYRECONVAR] = {"qxyrecon", ptopo, 1, 0,
                                 2}; // qxyrecon lives on primal edges
    aux_desc_arr[QXYEDGERECONVAR] = {
        "qxyedgerecon", ptopo, 2, 0,
        4}; // qxyedgerecon lives on primal vertical faces

    aux_desc_arr[FXYVAR] = {"FXY", dtopo, 0, ndims - 1,
                            1}; // FXY lives on twisted edges
    aux_desc_arr[CORIOLISXYRECONVAR] = {
        "coriolisxyrecon", ptopo, 1, 0,
        2}; // coriolisxyrecon lives on primal edges
    aux_desc_arr[CORIOLISXYEDGERECONVAR] = {
        "corliolisxyedgerecon", ptopo, 2, 0,
        4}; // coriolisxyedgerecon lives on primal vertical faces

    aux_desc_arr[FTXYVAR] = {"FTXY", ptopo, 1, 0, 1}; // FTXY
    const_desc_arr[CORIOLISXYVAR] = {"coriolisxy", ptopo, 2, 0, 1};
  }
}

std::unique_ptr<TestCase> make_coupled_test_case(PamCoupler &coupler);
void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance);

void read_model_params_file(std::string inFile, ModelParameters &params,
                            Parallel &par, PamCoupler &coupler,
                            std::unique_ptr<TestCase> &testcase) {
#ifdef PAM_STANDALONE
  // Read common parameters
  int nz = coupler.get_nz();
  readParamsFile(inFile, params, par, nz);

  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  params.acoustic_balance = config["balance_initial_density"].as<bool>(false);
  params.uniform_vertical = (config["vcoords"].as<std::string>() == "uniform");
  // Read diffusion coefficients
  params.scalar_horiz_diffusion_coeff =
      config["scalar_horiz_diffusion_coeff"].as<real>(0);
  params.scalar_vert_diffusion_coeff =
      config["scalar_vert_diffusion_coeff"].as<real>(0);
  params.velocity_vort_horiz_diffusion_coeff =
      config["velocity_vort_horiz_diffusion_coeff"].as<real>(0);
  params.velocity_vort_vert_diffusion_coeff =
      config["velocity_vort_vert_diffusion_coeff"].as<real>(0);
  params.velocity_div_horiz_diffusion_coeff =
      config["velocity_div_horiz_diffusion_coeff"].as<real>(0);
  params.velocity_div_vert_diffusion_coeff =
      config["velocity_div_vert_diffusion_coeff"].as<real>(0);

  params.scalar_diffusion_subtract_refstate =
      config["scalar_diffusion_subtract_refstate"].as<bool>(true);
  params.velocity_diffusion_subtract_refstate =
      config["velocity_diffusion_subtract_refstate"].as<bool>(false);
  // Read the data initialization options
  params.init_data = config["init_data"].as<std::string>();
  params.force_refstate_hydrostatic_balance =
      config["force_refstate_hydrostatic_balance"].as<bool>(false);
  params.check_anelastic_constraint =
      config["check_anelastic_constraint"].as<bool>(false);

  params.linear_system = config["linear_system"].as<std::string>("pressure");

  for (int i = 0; i < ntracers_dycore; i++) {
    params.init_dycore_tracer[i] =
        config["init_dycore_tracer" + std::to_string(i)].as<std::string>(
            "constant");
    params.dycore_tracer_pos[i] =
        config["dycore_tracer" + std::to_string(i) + "_pos"].as<bool>(false);
  }

  // Store vertical cell interface heights in the data manager
  auto &dm = coupler.get_data_manager_device_readonly();
  params.zint = dm.get<real const, 2>("vertical_interface_height");

  params.ylen = 1.0;
  params.yc = 0.5;

  testcase_from_string(testcase, params.init_data, params.acoustic_balance);
#endif
}

void read_model_params_coupler(ModelParameters &params, Parallel &par,
                               PamCoupler &coupler,
                               std::unique_ptr<TestCase> &testcase) {

  // Read common parameters
  read_params_coupler(params, par, coupler);

  params.acoustic_balance = false;
  params.uniform_vertical = false;
  params.scalar_horiz_diffusion_coeff = 0;
  params.scalar_vert_diffusion_coeff = 0;
  params.velocity_vort_horiz_diffusion_coeff = 0;
  params.velocity_vort_vert_diffusion_coeff = 0;
  params.velocity_div_horiz_diffusion_coeff = 0;
  params.velocity_div_vert_diffusion_coeff = 0;
  params.scalar_diffusion_subtract_refstate = true;
  params.velocity_diffusion_subtract_refstate = false;
  params.init_data = "coupler";
  params.force_refstate_hydrostatic_balance = true;
  params.check_anelastic_constraint = false;

  // Store vertical cell interface heights in the data manager
  auto &dm = coupler.get_data_manager_device_readonly();
  params.zint = dm.get<real const, 2>("vertical_interface_height");

  params.ylen = 1.0;
  params.yc = 0.5;

  testcase = make_coupled_test_case(coupler);
}

void check_and_print_model_parameters(const ModelParameters &params,
                                      const Parallel &par,
                                      bool verbose = false) {

  check_and_print_parameters(params, par, verbose);

  if (verbose) {
    serial_print("IC: " + params.init_data, par.masterproc);
    serial_print("acoustically balanced: " +
                     std::to_string(params.acoustic_balance),
                 par.masterproc);
    serial_print("scalar_horiz_diffusion_coeff: " +
                     std::to_string(params.scalar_horiz_diffusion_coeff),
                 par.masterproc);
    serial_print("scalar_vert_diffusion_coeff: " +
                     std::to_string(params.scalar_vert_diffusion_coeff),
                 par.masterproc);
    serial_print("velocity_vort_horiz_diffusion_coeff: " +
                     std::to_string(params.velocity_vort_horiz_diffusion_coeff),
                 par.masterproc);
    serial_print("velocity_vort_vert_diffusion_coeff: " +
                     std::to_string(params.velocity_vort_vert_diffusion_coeff),
                 par.masterproc);
    serial_print("velocity_div_horiz_diffusion_coeff: " +
                     std::to_string(params.velocity_div_horiz_diffusion_coeff),
                 par.masterproc);
    serial_print("velocity_div_vert_diffusion_coeff: " +
                     std::to_string(params.velocity_div_vert_diffusion_coeff),
                 par.masterproc);
    serial_print("scalar_diffusion_subtract_refstate: " +
                     std::to_string(params.scalar_diffusion_subtract_refstate),
                 par.masterproc);
    serial_print(
        "velocity_diffusion_subtract_refstate: " +
            std::to_string(params.velocity_diffusion_subtract_refstate),
        par.masterproc);

    for (int i = 0; i < ntracers_dycore; i++) {
      serial_print("Dycore Tracer" + std::to_string(i) +
                       " IC: " + params.init_dycore_tracer[i],
                   par.masterproc);
    }
  }
}

//***************** Test Cases ***************************//

// Universal

real YAKL_INLINE isentropic_T(real z, real theta0, real g,
                              const ThermoPotential &thermo) {
  return theta0 - z * g / thermo.cst.Cpd;
}

real YAKL_INLINE isentropic_p(real z, real theta0, real g,
                              const ThermoPotential &thermo) {
  return thermo.cst.pr * pow(isentropic_T(z, theta0, g, thermo) / theta0,
                             1. / thermo.cst.kappa_d);
}

real YAKL_INLINE isentropic_rho(real z, real theta0, real g,
                                const ThermoPotential &thermo) {
  real p = isentropic_p(z, theta0, g, thermo);
  real T = isentropic_T(z, theta0, g, thermo);
  real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
  return 1._fp / alpha;
}

real YAKL_INLINE isothermal_zdep(real z, real var_s, real T_ref, real g,
                                 const ThermoPotential &thermo) {
  real Rd = thermo.cst.Rd;
  real delta = g / (Rd * T_ref);
  return var_s * exp(-delta * z);
}

real YAKL_INLINE const_stability_p(real z, real N, real g, real ps, real Ts,
                                   const ThermoPotential &thermo) {
  real S = N * N / g;
  real G = g / (thermo.cst.Cpd * Ts * S);
  return ps * pow(1 - G * (1 - exp(-S * z)), 1 / thermo.cst.kappa_d);
}
real YAKL_INLINE const_stability_T(real z, real N, real g, real Ts,
                                   const ThermoPotential &thermo) {
  real S = N * N / g;
  real G = g / (thermo.cst.Cpd * Ts * S);
  return Ts * exp(S * z) * (1 - G * (1 - exp(-S * z)));
}

real YAKL_INLINE linear_ellipsoid(real x, real z, real x0, real z0, real xrad,
                                  real zrad, real amp) {
  real xn = (x - x0) / xrad;
  real zn = (z - z0) / zrad;
  real dist = sqrt(xn * xn + zn * zn);
  return amp * std::max(1._fp - dist, 0._fp);
}

real YAKL_INLINE flat_geop(real z, real g) { return g * z; }

// Returns saturation vapor pressure
real YAKL_INLINE saturation_vapor_pressure(real temp) {
  real tc = temp - 273.15_fp;
  return 610.94_fp * exp(17.625_fp * tc / (243.04_fp + tc));
}

real YAKL_INLINE saturation_mixing_ratio(real T, real p) {
  return 380 / p * exp(17.27 * (T - 273) / (T - 36));
}

template <class T> class SWETestCase : public TestCase, public T {
public:
  using T::g;

  using T::Lx;
  using T::Ly;
  using T::Lz;
  using T::xc;
  using T::yc;

  using T::h_f;
#ifdef PAMC_TSWE
  using T::S_f;
#endif
  using T::coriolis_f;
  using T::v_f;

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
    params.ylen = T::Ly;
    params.yc = T::yc;
  }

  std::array<real, 3> get_domain() const override { return {Lx, Ly, Lz}; }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return h_f(x, y, z); },
        progvars.fields_arr[DENSVAR], 0);
#ifdef PAMC_TSWE
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return S_f(x, y, z); },
        progvars.fields_arr[DENSVAR], 1);
#endif
    primal_geom.set_10form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[VVAR], 0);
    primal_geom.set_01form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[WVAR], 0);

    if (ndims == 1) {
      primal_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) { return coriolis_f(x, y, z).v; },
          constvars.fields_arr[CORIOLISHZVAR], 0);
    } else {
      primal_geom.set_nm11form_values(
          YAKL_LAMBDA(real x, real y, real z) { return coriolis_f(x, y, z); },
          constvars.fields_arr[CORIOLISHZVAR], 0);

      primal_geom.set_n0form_values(
          YAKL_LAMBDA(real x, real y, real z) { return coriolis_f(x, y, z); },
          constvars.fields_arr[CORIOLISXYVAR], 0);
    }

    // YAKL_SCOPE(tracer_f, this->tracer_f);
    // for (int i = 0; i < ntracers_dycore; i++) {
    //   dual_geom.set_11form_values(
    //       YAKL_LAMBDA(real x, real y) {
    //         return h_f(x, y) * tracer_f(i)->compute(x, y, Lx, Ly, xc, yc);
    //       },
    //       progvars.fields_arr[DENSVAR], i + ndensity_dycore);
    // }
    this->equations->Hs.set_parameters(g);
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = this->equations->reference_state;

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    const int pks = primal_topology.ks;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(varset, equations->varset);
    YAKL_SCOPE(Hs, equations->Hs);

    // real href = T::H0;
    real href = 0;

    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return 0; }, refstate.geop, 0);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return href; }, refstate.dens,
        varset.dens_id_mass);

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return href; }, refstate.rho_pi,
        varset.dens_id_mass);
    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return 1; }, refstate.q_pi, varset.dens_id_mass);

    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return href; }, refstate.rho_di,
        varset.dens_id_mass);
    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return 1; }, refstate.q_di, varset.dens_id_mass);
  }
};

template <class T> class EulerTestCase : public TestCase, public T {
  using VS = VariableSet;

public:
  using T::g;
  using T::Lx;
  using T::Lz;
  using T::xc;
  using T::zc;

  using T::entropicvar_f;

  using T::refentropicdensity_f;
  using T::refnsq_f;
  using T::refrho_f;
  using T::v_f;

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
#ifdef PAMC_AN
    return refrho_f(z, thermo);
#else
    return T::rho_f(x, y, z, thermo);
#endif
  }

  std::array<real, 3> get_domain() const override {
    real Ly = 1;
    if constexpr (T::max_ndims > 1) {
      Ly = T::Ly;
    }
    return {Lx, Ly, Lz};
  }

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
    if constexpr (T::max_ndims > 1) {
      params.ylen = T::Ly;
      params.yc = T::yc;
    }
  }

  void add_diagnostics(
      std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {
    T::add_diagnostics(diagnostics);
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(varset, equations->varset);
#ifndef PAMC_AN
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return rho_f(x, y, z, thermo); },
        progvars.fields_arr[DENSVAR], varset.dens_id_mass);
#endif
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) {
          return rho_f(x, y, z, thermo) * entropicvar_f(x, y, z, thermo);
        },
        progvars.fields_arr[DENSVAR], varset.dens_id_entr);

    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return flat_geop(z, g); },
        constvars.fields_arr[HSVAR], 0);

    primal_geom.set_10form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[VVAR], 0);
    primal_geom.set_01form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[WVAR], 0);

    YAKL_SCOPE(tracers, this->tracers);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return rho_f(x, y, z, thermo) *
                   TracerFunctor{}(tracers(i), x, z, Lx, Lz, xc, zc);
          },
          progvars.fields_arr[DENSVAR], i + VS::ndensity_dycore_prognostic);
    }
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = this->equations->reference_state;

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    const int pks = primal_topology.ks;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(varset, equations->varset);
    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(Hs, equations->Hs);

    Hs.set_parameters(g);

    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return flat_geop(z, g); }, refstate.geop, 0);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.dens,
        varset.dens_id_mass);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.dens, varset.dens_id_entr);

    parallel_for(
        "compute rho_pi and unscaled q_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          SArray<real, 1, VS::ndensity> dens0;
          compute_Hn1bar<VS::ndensity, vert_diff_ord>(
              dens0, refstate.dens.data, primal_geom, dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k + pks, n) = dens0(varset.dens_id_mass);
          for (int d = 0; d < VS::ndensity; ++d) {
            refstate.q_pi.data(d, k + pks, n) = dens0(d);
          }
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refnsq_f(z, thermo); }, refstate.Nsq_pi,
        0);

    parallel_for(
        "compute rho_di and unscaled q_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, dks, n) = refstate.rho_pi.data(0, pks, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1 + dks, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1 + pks, n);
          } else {
            refstate.rho_di.data(0, k + dks, n) =
                0.5_fp * (refstate.rho_pi.data(0, k + pks, n) +
                          refstate.rho_pi.data(0, k - 1 + pks, n));
          }

          for (int d = 0; d < VS::ndensity; ++d) {
            if (k == 0) {
              refstate.q_di.data(d, dks, n) = refstate.q_pi.data(d, pks, n);
            } else if (k == dual_topology.ni - 1) {
              refstate.q_di.data(d, dual_topology.ni - 1 + dks, n) =
                  refstate.q_pi.data(d, primal_topology.ni - 1 + pks, n);
            } else {
              refstate.q_di.data(d, k + dks, n) =
                  0.5_fp * (refstate.q_pi.data(d, k + pks, n) +
                            refstate.q_pi.data(d, k - 1 + pks, n));
            }
          }
        });

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_pi.data(l, k + pks, n) /=
                refstate.rho_pi.data(0, k + pks, n);
          }
        });

    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_di.data(l, k + dks, n) /=
                refstate.rho_di.data(0, k + dks, n);
          }
        });

    parallel_for(
        "Compute refstate B",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          Hs.compute_dHsdx(refstate.B.data, refstate.dens.data,
                           refstate.geop.data, pks, k, n, -1);
        });

    parallel_for(
        "Compute refstate pres_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_pi.data(0, k + pks, n);
          const real entropicvar =
              refstate.q_pi.data(varset.dens_id_entr, k + pks, n);
          refstate.pres_pi.data(0, k + pks, n) =
              thermo.solve_p(rho, entropicvar, 1, 0, 0, 0);
        });

    parallel_for(
        "Compute refstate pres_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_di.data(0, k + dks, n);
          const real entropicvar =
              refstate.q_di.data(varset.dens_id_entr, k + dks, n);
          refstate.pres_di.data(0, k + dks, n) =
              thermo.solve_p(rho, entropicvar, 1, 0, 0, 0);
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return T::refnsq_f(z, thermo); }, refstate.Nsq_pi,
        0);
  }
};

template <class T> class MoistEulerTestCase : public TestCase, public T {
  using VS = VariableSet;

public:
  using T::g;
  using T::Lx;
  using T::Lz;
  using T::xc;
  using T::zc;

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
#ifdef PAMC_MAN
    return T::refrho_f(z, thermo);
#elif defined PAMC_MCErhod || defined PAMC_MCErhodp
    return T::rhod_f(x, y, z, thermo);
#else
    return T::rho_f(x, y, z, thermo);
#endif
  }

  std::array<real, 3> get_domain() const override {
    real Ly = 1;
    if constexpr (T::max_ndims > 1) {
      Ly = T::Ly;
    }
    return {Lx, Ly, Lz};
  }

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
    if constexpr (T::max_ndims > 1) {
      params.ylen = T::Ly;
      params.yc = T::yc;
    }
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return flat_geop(z, g); },
        constvars.fields_arr[HSVAR], 0);

    if constexpr (T::needs_special_init) {
      T::initialize(equations, progvars, constvars, primal_geom, dual_geom);
    } else {
      YAKL_SCOPE(thermo, equations->thermo);
      YAKL_SCOPE(varset, equations->varset);

#ifndef PAMC_MAN
      dual_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return T::rho_f(x, y, z, thermo);
          },
          progvars.fields_arr[DENSVAR], varset.dens_id_mass);
#endif
      dual_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return T::rho_f(x, y, z, thermo) *
                   T::entropicvar_f(x, y, z, thermo);
          },
          progvars.fields_arr[DENSVAR], varset.dens_id_entr);

      dual_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return T::rhov_f(x, y, z, thermo);
          },
          progvars.fields_arr[DENSVAR], varset.dens_id_vap);

      YAKL_SCOPE(tracers, this->tracers);
      for (int i = 0; i < ntracers_dycore; i++) {
        dual_geom.set_n1form_values(
            YAKL_LAMBDA(real x, real y, real z) {
              return T::rho_f(x, y, z, thermo) *
                     TracerFunctor{}(tracers(i), x, z, Lx, Lz, xc, zc);
            },
            progvars.fields_arr[DENSVAR], i + VS::ndensity_dycore_prognostic);
      }
    }
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = this->equations->reference_state;

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    const int pks = primal_topology.ks;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(Hs, equations->Hs);
    YAKL_SCOPE(varset, equations->varset);

    Hs.set_parameters(g);

    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return flat_geop(z, g); }, refstate.geop, 0);

    if constexpr (T::needs_special_init) {
      T::initialize_refstate(equations, primal_geom, dual_geom);
    } else {
      dual_geom.set_profile_n1form_values(
          YAKL_LAMBDA(real z) { return T::refrho_f(z, thermo); }, refstate.dens,
          varset.dens_id_mass);
      dual_geom.set_profile_n1form_values(
          YAKL_LAMBDA(real z) { return T::refentropicdensity_f(z, thermo); },
          refstate.dens, varset.dens_id_entr);
      dual_geom.set_profile_n1form_values(
          YAKL_LAMBDA(real z) { return T::refrhov_f(z, thermo); },
          refstate.dens, varset.dens_id_vap);
    }

    parallel_for(
        "compute rho_pi and unscaled q_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const auto total_density_f = TotalDensityFunctor{varset};
          SArray<real, 1, 1> rho0;
          compute_Hn1bar<1, vert_diff_ord>(total_density_f, rho0,
                                           refstate.dens.data, primal_geom,
                                           dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k + pks, n) = rho0(0);

          SArray<real, 1, VS::ndensity> dens0;
          compute_Hn1bar<VS::ndensity, vert_diff_ord>(
              dens0, refstate.dens.data, primal_geom, dual_geom, pks, k, n);
          for (int d = 0; d < VS::ndensity; ++d) {
            refstate.q_pi.data(d, k + pks, n) = dens0(d);
          }
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return T::refnsq_f(z, thermo); }, refstate.Nsq_pi,
        0);

    parallel_for(
        "compute rho_di and unscaled q_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, dks, n) = refstate.rho_pi.data(0, pks, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1 + dks, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1 + pks, n);
          } else {
            refstate.rho_di.data(0, k + dks, n) =
                0.5_fp * (refstate.rho_pi.data(0, k + pks, n) +
                          refstate.rho_pi.data(0, k - 1 + pks, n));
          }

          for (int d = 0; d < VS::ndensity; ++d) {
            if (k == 0) {
              refstate.q_di.data(d, dks, n) = refstate.q_pi.data(d, pks, n);
            } else if (k == dual_topology.ni - 1) {
              refstate.q_di.data(d, dual_topology.ni - 1 + dks, n) =
                  refstate.q_pi.data(d, primal_topology.ni - 1 + pks, n);
            } else {
              refstate.q_di.data(d, k + dks, n) =
                  0.5_fp * (refstate.q_pi.data(d, k + pks, n) +
                            refstate.q_pi.data(d, k - 1 + pks, n));
            }
          }
        });

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_pi.data(l, k + pks, n) /=
                refstate.rho_pi.data(0, k + pks, n);
          }
        });
    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_di.data(l, k + dks, n) /=
                refstate.rho_di.data(0, k + dks, n);
          }
        });

    parallel_for(
        "Compute refstate B",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          Hs.compute_dHsdx(refstate.B.data, refstate.dens.data,
                           refstate.geop.data, pks, k, n, -1);
        });

    parallel_for(
        "Compute refstate pres_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_pi.data(0, k + pks, n);
          const real entropicvar =
              refstate.q_pi.data(varset.dens_id_entr, k + pks, n);
          const real qv = refstate.q_pi.data(varset.dens_id_vap, k + pks, n);
          refstate.pres_pi.data(0, k + pks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });

    parallel_for(
        "Compute refstate pres_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_di.data(0, k + dks, n);
          const real entropicvar =
              refstate.q_di.data(varset.dens_id_entr, k + dks, n);
          const real qv = refstate.q_di.data(varset.dens_id_vap, k + dks, n);
          refstate.pres_di.data(0, k + dks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return T::refnsq_f(z, thermo); }, refstate.Nsq_pi,
        0);
  }

  void add_diagnostics(
      std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {
    T::add_diagnostics(diagnostics);
  }
};

class CoupledTestCase : public TestCase {
  using VS = VariableSet;

public:
  PamCoupler &coupler;
  bool set_coupler_state = false;

  CoupledTestCase(PamCoupler &coupler) : coupler(coupler) {}

  std::array<real, 3> get_domain() const override { return {0, 1, 0}; }

  void set_domain(ModelParameters &params) override {
    params.xlen = coupler.get_xlen();
    params.xc = 0;
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    auto &varset = this->equations->varset;
    const real g = this->equations->Hs.g;
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return flat_geop(z, g); },
        constvars.fields_arr[HSVAR], 0);

    convert_coupler_to_dynamics_densities(varset, coupler, progvars, constvars);
    convert_coupler_to_dynamics_wind(varset, coupler, progvars, constvars,
                                     false);
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {

    auto &refstate = this->equations->reference_state;
    const auto &varset = this->equations->varset;
    auto &thermo = this->equations->thermo;
    const auto &primal_topology = primal_geom.topology;
    const auto &dual_topology = dual_geom.topology;

    // Set thermo constants based on coupler values
    // Need a better way to separate fundamental vs derived constants
    // Also I think only P3 defines all of these

    thermo.cst.Rd = coupler.get_option<real>("R_d");
    thermo.cst.Rv = coupler.get_option<real>("R_v");
    thermo.cst.pr = coupler.get_option<real>("p0");
    thermo.cst.Cpd = coupler.get_option<real>("cp_d");
    thermo.cst.Cvd = thermo.cst.Cpd - thermo.cst.Rd;
    thermo.cst.Cpv = coupler.get_option<real>("cp_v");

    thermo.cst.Lvr = coupler.get_option<real>("latvap");
    thermo.cst.Lfr = coupler.get_option<real>("latice");

    thermo.cst.gamma_d = thermo.cst.Cpd / thermo.cst.Cvd;
    thermo.cst.kappa_d = thermo.cst.Rd / thermo.cst.Cpd;
    thermo.cst.delta_d = thermo.cst.Rd / thermo.cst.Cvd;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readonly();

    auto dm_ref_dens_dry = dm.get<real const, 2>("ref_density_dry");
    auto dm_ref_dens_vap = dm.get<real const, 2>("ref_density_vapor");
    auto dm_ref_dens_liq = dm.get<real const, 2>("ref_density_liq");
    auto dm_ref_dens_ice = dm.get<real const, 2>("ref_density_ice");
    auto dm_ref_temp = dm.get<real const, 2>("ref_temp");

    const real grav = coupler.get_option<real>("grav");
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return flat_geop(z, grav); }, refstate.geop, 0);

    // sets dens and unscaled q_pi
    parallel_for(
        "Coupled reference state 1",
        SimpleBounds<2>(dual_topology.nl, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int n) {
          const real temp = dm_ref_temp(k, n);
          const real dens_dry = dm_ref_dens_dry(k, n);
          const real dens_vap = dm_ref_dens_vap(k, n);
          const real dens_liq = dm_ref_dens_liq(k, n);
          const real dens_ice = dm_ref_dens_ice(k, n);
          const real dens = dens_dry + dens_vap;

          const real qd = dens_dry / dens;
          const real qv = dens_vap / dens;
          const real ql = dens_liq / dens;
          const real qi = dens_ice / dens;

          const real alpha = 1.0_fp / dens;
          const real entropic_var = thermo.compute_entropic_var_from_alpha_T(
              alpha, temp, qd, qv, ql, qi);

          const real dual_volume =
              dual_geom.get_area_n1entity(k + dks, djs, dis, n);
          refstate.dens.data(varset.dens_id_mass, k + dks, n) =
              dens * dual_volume;
          refstate.dens.data(varset.dens_id_entr, k + dks, n) =
              entropic_var * dens * dual_volume;
          refstate.dens.data(varset.dens_id_vap, k + dks, n) =
              dens_vap * dual_volume;

          refstate.q_pi.data(varset.dens_id_mass, k + dks, n) = dens;
          refstate.q_pi.data(varset.dens_id_entr, k + dks, n) =
              dens * entropic_var;
          refstate.q_pi.data(varset.dens_id_vap, k + dks, n) = dens_vap;
        });

    parallel_for(
        "compute unscaled q_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int d = 0; d < VS::ndensity; ++d) {
            if (k == 0) {
              refstate.q_di.data(d, dks, n) = refstate.q_pi.data(d, pks, n);
            } else if (k == dual_topology.ni - 1) {
              refstate.q_di.data(d, dual_topology.ni - 1 + dks, n) =
                  refstate.q_pi.data(d, primal_topology.ni - 1 + pks, n);
            } else {
              real q_km1 = refstate.q_pi.data(d, k - 1 + pks, n);
              real q_k = refstate.q_pi.data(d, k + pks, n);

              refstate.q_di.data(d, k + dks, n) =
                  q_km1 + (q_k - q_km1) *
                              (dual_geom.zint(k + pks, n) -
                               primal_geom.zint(k - 1 + pks, n)) /
                              primal_geom.dz(k - 1 + pks, n);
            }
          }
        });

    parallel_for(
        "compute rho_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const auto total_density_f = TotalDensityFunctor{varset};
          SArray<real, 1, 1> rho0;
          compute_Hn1bar<1, vert_diff_ord>(total_density_f, rho0,
                                           refstate.dens.data, primal_geom,
                                           dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k + pks, n) = rho0(0);
        });

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_pi.data(l, k + pks, n) /=
                refstate.rho_pi.data(0, k + pks, n);
          }
#if defined PAMC_AN || defined PAMC_MAN
          refstate.q_pi.data(varset.dens_id_mass, k + pks, n) = 1;
#endif
        });

    parallel_for(
        "compute rho_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, dks, n) = refstate.rho_pi.data(0, pks, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1 + dks, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1 + pks, n);
          } else {
            real rho_k = refstate.rho_pi.data(0, k + pks, n);
            real rho_km1 = refstate.rho_pi.data(0, k - 1 + pks, n);

            refstate.rho_di.data(0, k + dks, n) =
                rho_km1 + (rho_k - rho_km1) *
                              (dual_geom.zint(k + pks, n) -
                               primal_geom.zint(k - 1 + pks, n)) /
                              primal_geom.dz(k - 1 + pks, n);
          }
        });
    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_di.data(l, k + dks, n) /=
                refstate.rho_di.data(0, k + dks, n);
          }
#if defined PAMC_AN || defined PAMC_MAN
          refstate.q_di.data(varset.dens_id_mass, k + pks, n) = 1;
#endif
        });

    YAKL_SCOPE(Hs, equations->Hs);
    Hs.set_parameters(grav);

    parallel_for(
        "Compute refstate B",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          Hs.compute_dHsdx(refstate.B.data, refstate.dens.data,
                           refstate.geop.data, pks, k, n, -1);
        });

    parallel_for(
        "compute Nsq",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real Rd = thermo.cst.Rd;
          const real Rv = thermo.cst.Rv;
          const real Cpd = thermo.cst.Cpd;
          const real Lvr = thermo.cst.Lvr;
          const real eta = Rv / Rd;

          real T_kp, T_km, T;
          real rv_kp, rv_km, rv;
          real dz;

          if (k == 0) {
            T_km = dm_ref_temp(k, n);
            T_kp = dm_ref_temp(k + 1, n);
            T = T_km;

            rv_kp = dm_ref_dens_vap(k + 1, n) / dm_ref_dens_dry(k + 1, n);
            rv_km = dm_ref_dens_vap(k, n) / dm_ref_dens_dry(k, n);
            rv = rv_km;

            dz = primal_geom.dz(k + pks, n);
          } else if (k == primal_topology.ni - 1) {
            T_km = dm_ref_temp(k - 1, n);
            T_kp = dm_ref_temp(k, n);
            T = dm_ref_temp(k, n);

            rv_kp = dm_ref_dens_vap(k, n) / dm_ref_dens_dry(k, n);
            rv_km = dm_ref_dens_vap(k - 1, n) / dm_ref_dens_dry(k - 1, n);
            rv = rv_kp;

            dz = primal_geom.dz(k - 1 + pks, n);
          } else {
            T = dm_ref_temp(k, n);
            T_km = dm_ref_temp(k - 1, n);
            T_kp = dm_ref_temp(k + 1, n);

            rv_km = dm_ref_dens_vap(k - 1, n) / dm_ref_dens_dry(k - 1, n);
            rv_kp = dm_ref_dens_vap(k + 1, n) / dm_ref_dens_dry(k + 1, n);
            rv = dm_ref_dens_vap(k, n) / dm_ref_dens_dry(k, n);

            dz = primal_geom.dz(k + pks, n) + primal_geom.dz(k - 1 + pks, n);
          }
          real dTdz = (T_kp - T_km) / dz;
          real drvdz = (rv_kp - rv_km) / dz;

          real Tv = T * (1 + eta * rv) / (1 + rv);
          real es = saturation_vapor_pressure(T);
          real rsw = (es / (Rd * T) - 1) * Rd / Rv;
          real qsw = rsw / (1 + rsw);

          real D1w = 1 + (1 + eta * rsw) * Lvr * qsw / (Rd * Tv);
          real D2w = 1 + (1 + eta * rsw) * Lvr * Lvr * qsw / (Cpd * Rv * T * T);
          real gamma_m = grav / Cpd * D1w / D2w;

          refstate.Nsq_pi.data(0, k + pks, n) =
              grav / T * D1w * (dTdz + gamma_m) - grav / (1 + rv) * drvdz;
        });

    parallel_for(
        "Compute refstate pres_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_pi.data(0, k + pks, n);
          const real entropicvar =
              refstate.q_pi.data(varset.dens_id_entr, k + pks, n);
          const real qv = refstate.q_pi.data(varset.dens_id_vap, k + pks, n);
          refstate.pres_pi.data(0, k + pks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });

    parallel_for(
        "Compute refstate pres_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_di.data(0, k + dks, n);
          const real entropicvar =
              refstate.q_di.data(varset.dens_id_entr, k + dks, n);
          const real qv = refstate.q_di.data(varset.dens_id_vap, k + dks, n);
          refstate.pres_di.data(0, k + dks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });
  }

  void add_diagnostics(
      std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {}
};

struct DoubleVortex {
  enum class PLANE { XZ, YZ, XY };

  static PLANE constexpr plane = PLANE::XZ;
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 5000000._fp;
  static real constexpr Ly = 5000000._fp;
  static real constexpr Lz = 5000000._fp;
  static real constexpr coriolis = 0.00006147_fp;
  static real constexpr H0 = 750.0_fp;
  static real constexpr ox = 0.1_fp;
  static real constexpr oy = 0.1_fp;
  static real constexpr sigmax = 3._fp / 40._fp * Lx;
  static real constexpr sigmay = 3._fp / 40._fp * Ly;
  static real constexpr dh = 75.0_fp;
  static real constexpr xc1 = (0.5_fp - ox) * Lx;
  static real constexpr yc1 = (0.5_fp - oy) * Ly;
  static real constexpr xc2 = (0.5_fp + ox) * Lx;
  static real constexpr yc2 = (0.5_fp + oy) * Ly;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr yc = 0.5_fp * Ly;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr c = 0.05_fp;
  static real constexpr a = 1.0_fp / 3.0_fp;
  static real constexpr D = 0.5_fp * Lx;

  static std::pair<real, real> YAKL_INLINE get_plane_coords(real x, real y,
                                                            real z) {
    if (plane == PLANE::XY) {
      return {x, y};
    }
    if (plane == PLANE::XZ) {
      return {x, z};
    }
    if (plane == PLANE::YZ) {
      return {y, z};
    }
  }

  static VecXYZ YAKL_INLINE coriolis_f(real x, real y, real z) {
    VecXYZ v = {0, 0, 0};
    if (plane == PLANE::XY) {
      v.w = coriolis;
    }
    if (plane == PLANE::XZ) {
      v.v = ndims > 1 ? -coriolis : coriolis;
    }
    if (plane == PLANE::YZ) {
      v.u = coriolis;
    }

    return v;
  }

  static real YAKL_INLINE h_f(real xx, real yy, real zz) {
    auto [x, y] = get_plane_coords(xx, yy, zz);

    real xprime1 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc1));
    real yprime1 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc1));
    real xprime2 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc2));
    real yprime2 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc2));
    real xprimeprime1 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc1));
    real yprimeprime1 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc1));
    real xprimeprime2 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc2));
    real yprimeprime2 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc2));

    return H0 - dh * (exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
                      exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)) -
                      4._fp * pi * sigmax * sigmay / Lx / Ly);
  }

  static VecXYZ YAKL_INLINE v_f(real xx, real yy, real zz) {
    auto [x, y] = get_plane_coords(xx, yy, zz);

    real xprime1 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc1));
    real yprime1 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc1));
    real xprime2 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc2));
    real yprime2 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc2));
    real xprimeprime1 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc1));
    real yprimeprime1 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc1));
    real xprimeprime2 =
        Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc2));
    real yprimeprime2 =
        Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc2));

    real u1 =
        -g * dh / coriolis / sigmay *
        (yprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         yprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));
    real u2 =
        g * dh / coriolis / sigmax *
        (xprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         xprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));

    VecXYZ vvec = {0, 0, 0};

    if (plane == PLANE::XY) {
      vvec.u = u1;
      vvec.v = u2;
    }
    if (plane == PLANE::XZ) {
      vvec.u = u1;
      vvec.w = u2;
    }
    if (plane == PLANE::YZ) {
      vvec.v = u1;
      vvec.w = u2;
    }
    return vvec;
  }

  static real YAKL_INLINE S_f(real xx, real yy, real zz) {
    auto [x, y] = get_plane_coords(xx, yy, zz);
    real sval =
        g * (1._fp + c * exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) /
                             (a * a * D * D)));
    return sval * h_f(xx, yy, zz);
  }
};

// contains defualts
struct TestCaseInit {
  static int constexpr max_ndims = 1;
  static bool constexpr needs_special_init = false;
};

template <bool acoustic_balance> struct RisingBubble : TestCaseInit {
  static int constexpr max_ndims = 2;
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 1000._fp;
  static real constexpr Ly = 1000._fp;
  static real constexpr Lz = 1500._fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr yc = 0.5_fp * Ly;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr theta0 = 300.0_fp;
  static real constexpr bzc = 350._fp;
  static real constexpr dss = 0.5_fp;
  static real constexpr rc = 250._fp;
  static real constexpr rh0 = 0.8_fp;
  static real constexpr T_ref = 300.0_fp;
  static real constexpr N_ref = 0.0001;

  static real YAKL_INLINE refnsq_f(real z, const ThermoPotential &thermo) {
    return N_ref * N_ref;
  }

  static real YAKL_INLINE refp_f(real z, const ThermoPotential &thermo) {
    return const_stability_p(z, N_ref, g, thermo.cst.pr, theta0, thermo);
  }

  static real YAKL_INLINE refT_f(real z, const ThermoPotential &thermo) {
    return const_stability_T(z, N_ref, g, theta0, thermo);
  }

  static real YAKL_INLINE refrho_f(real z, const ThermoPotential &thermo) {
    real p = refp_f(z, thermo);
    real T = refT_f(z, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE refentropicdensity_f(real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(z, thermo);
    real T_ref = refT_f(z, thermo);
    real p_ref = refp_f(z, thermo);

    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rho_b = isentropic_rho(z, theta0, g, thermo);
    if (acoustic_balance) {
      real theta = entropicvar_f(x, y, z, thermo);
      return rho_b * theta0 / theta;
    } else {
      return rho_b;
    }
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real r;
    if (ndims == 1) {
      r = sqrt((x - xc) * (x - xc) + (z - bzc) * (z - bzc));
    } else {
      r = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc) +
               (z - bzc) * (z - bzc));
    }
    real dtheta = (r < rc) ? dss * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.v = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct TwoBubbles : TestCaseInit {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 1000._fp;
  static real constexpr Lz = 1000._fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr theta0 = 303.15_fp;

  static real constexpr A1 = 0.5_fp;
  static real constexpr a1 = 150;
  static real constexpr s1 = 50;
  static real constexpr x1 = 500;
  static real constexpr z1 = 300;

  static real constexpr A2 = -0.15_fp;
  static real constexpr a2 = 0;
  static real constexpr s2 = 50;
  static real constexpr x2 = 560;
  static real constexpr z2 = 640;

  static real constexpr T_ref = 303.15_fp;
  static real constexpr N_ref = 0.0001;

  static real YAKL_INLINE refnsq_f(real z, const ThermoPotential &thermo) {
    return N_ref * N_ref;
  }

  static real YAKL_INLINE refp_f(real z, const ThermoPotential &thermo) {
    return const_stability_p(z, N_ref, g, thermo.cst.pr, theta0, thermo);
  }

  static real YAKL_INLINE refT_f(real z, const ThermoPotential &thermo) {
    return const_stability_T(z, N_ref, g, theta0, thermo);
  }

  static real YAKL_INLINE refrho_f(real z, const ThermoPotential &thermo) {
    real p = refp_f(z, thermo);
    real T = refT_f(z, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE refentropicdensity_f(real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(z, thermo);
    real T_ref = refT_f(z, thermo);
    real p_ref = refp_f(z, thermo);
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    return isentropic_rho(z, theta0, g, thermo);
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);

    real dtheta = 0;

    real r1 = sqrt((x - x1) * (x - x1) + (z - z1) * (z - z1));
    if (r1 <= a1) {
      dtheta += A1;
    } else {
      dtheta += A1 * exp(-(r1 - a1) * (r1 - a1) / (s1 * s1));
    }

    real r2 = sqrt((x - x2) * (x - x2) + (z - z2) * (z - z2));
    if (r2 <= a2) {
      dtheta += A2;
    } else {
      dtheta += A2 * exp(-(r2 - a2) * (r2 - a2) / (s2 * s2));
    }

    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct DensityCurrent : TestCaseInit {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 51.2e3;
  static real constexpr Lz = 6400;
  static real constexpr xc = 0;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr theta0 = 300.0_fp;
  static real constexpr bxc = 0;
  static real constexpr bzc = 3000;
  static real constexpr bxr = 4000;
  static real constexpr bzr = 2000;
  static real constexpr dss = -15;
  static real constexpr T_ref = 300.0_fp;
  static real constexpr N_ref = 0.0001;

  static real YAKL_INLINE refnsq_f(real z, const ThermoPotential &thermo) {
    return N_ref * N_ref;
  }

  static real YAKL_INLINE refp_f(real z, const ThermoPotential &thermo) {
    return const_stability_p(z, N_ref, g, thermo.cst.pr, theta0, thermo);
  }

  static real YAKL_INLINE refT_f(real z, const ThermoPotential &thermo) {
    return const_stability_T(z, N_ref, g, theta0, thermo);
  }

  static real YAKL_INLINE refrho_f(real z, const ThermoPotential &thermo) {
    real p = refp_f(z, thermo);
    real T = refT_f(z, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE refentropicdensity_f(real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(z, thermo);
    real T_ref = refT_f(z, thermo);
    real p_ref = refp_f(z, thermo);
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rho_b = isentropic_rho(z, theta0, g, thermo);
    return rho_b;
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real r = sqrt((x - bxc) * (x - bxc) / (bxr * bxr) +
                  (z - bzc) * (z - bzc) / (bzr * bzr));
    real dtheta = (r < 1) ? dss * 0.5_fp * (1._fp + cos(pi * r)) : 0._fp;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct MoistRisingBubble : public RisingBubble<false> {

  static real YAKL_INLINE rhov_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {
    real r = (x - xc) * (x - xc) + (z - bzc) * (z - bzc);
    if (ndims > 1) {
      r += (y - yc) * (y - yc);
    }
    r = std::sqrt(r);

    real rh = (r < rc) ? rh0 * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
    real Th = isentropic_T(z, theta0, g, thermo);
    real svp = saturation_vapor_pressure(Th);
    real pv = svp * rh;
    return pv / (thermo.cst.Rv * Th);
  }

  static real YAKL_INLINE rhod_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rhod = rhod_f(x, y, z, thermo);
    real rhov = rhov_f(x, y, z, thermo);
    return rhod + rhov;
  }

  static real YAKL_INLINE refrhov_f(real z, const ThermoPotential &thermo) {
    return 0;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct LargeRisingBubble : TestCaseInit {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 20000._fp;
  static real constexpr Lz = 20000._fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr theta0 = 300.0_fp;
  static real constexpr bzc = 2000._fp;
  static real constexpr xrad = 2000._fp;
  static real constexpr zrad = 2000._fp;
  static real constexpr amp_theta = 2.0_fp;
  static real constexpr amp_vapor = 0.8_fp;
  static real constexpr Cpv = 1859._fp;
  static real constexpr Cpd = 1003._fp;
  static real constexpr N_ref = 0.0001;

  static real YAKL_INLINE refnsq_f(real z, const ThermoPotential &thermo) {
    return N_ref * N_ref;
  }

  static real YAKL_INLINE refp_f(real z, const ThermoPotential &thermo) {
    return const_stability_p(z, N_ref, g, thermo.cst.pr, theta0, thermo);
  }

  static real YAKL_INLINE refT_f(real z, const ThermoPotential &thermo) {
    return const_stability_T(z, N_ref, g, theta0, thermo);
  }

  static real YAKL_INLINE refrho_f(real z, const ThermoPotential &thermo) {
    real p = refp_f(z, thermo);
    real T = refT_f(z, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE refentropicdensity_f(real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(z, thermo);
    real T_ref = refT_f(z, thermo);
    real p_ref = refp_f(z, thermo);
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    return isentropic_rho(z, theta0, g, thermo);
  }
  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {

    real p = isentropic_p(z, theta0, g, thermo);
    real T0 = isentropic_T(z, theta0, g, thermo);
    real dtheta = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_theta);
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T0 + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct MoistLargeRisingBubble : LargeRisingBubble {

  static real YAKL_INLINE rhod_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE rhov_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {

    real pert = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_vapor);
    real Th = isentropic_T(z, theta0, g, thermo);
    real svp = saturation_vapor_pressure(Th);
    real pv = svp * pert;
    return pv / (thermo.cst.Rv * Th);
  }

  static real YAKL_INLINE refrhov_f(real z, const ThermoPotential &thermo) {
    return 0;
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rhod = rhod_f(x, y, z, thermo);
    real rhov = rhov_f(x, y, z, thermo);
    return rhod + rhov;
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {

    real p = isentropic_p(z, theta0, g, thermo);
    real T0 = isentropic_T(z, theta0, g, thermo);
    real dtheta = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_theta);
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T0 + dT, 1, 0, 0, 0);
  }
};

template <bool add_perturbation> struct GravityWave : TestCaseInit {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 300e3_fp;
  static real constexpr Lz = 10e3_fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr d = 5e3_fp;
  static real constexpr T_ref = 250._fp;
  static real constexpr u_0 = 20._fp;
  static real constexpr x_c = 100e3_fp;
  static real constexpr p_s = 1e5_fp;
  static real constexpr dT_max = 0.01_fp;

  static real YAKL_INLINE refnsq_f(real z, const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real gamma_d = thermo.cst.gamma_d;
    real N2 = (gamma_d - 1) / gamma_d * g * g / (Rd * T_ref);
    return N2;
  }

  static real YAKL_INLINE refrho_f(real z, const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real rho_s = p_s / (Rd * T_ref);
    real rho_ref = isothermal_zdep(z, rho_s, T_ref, g, thermo);
    return rho_ref;
  }

  static real YAKL_INLINE refp_f(real z, const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real rho_ref = refrho_f(z, thermo);
    real p_ref = Rd * rho_ref * T_ref;
    return p_ref;
  }

  static real YAKL_INLINE refentropicdensity_f(real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(z, thermo);
    real p_ref = refp_f(z, thermo);
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;

    real delta = g / (Rd * T_ref);
    real rho_s = p_s / (Rd * T_ref);

    real rho_ref = refrho_f(z, thermo);

    real dT_b = dT_max * exp(-pow((x - x_c) / d, 2)) * sin(pi * z / Lz);
    real dT = exp(delta * z / 2) * dT_b;

    real drho_b = -rho_s * dT_b / T_ref;
    real drho = exp(-delta * z / 2) * drho_b;

    real T = T_ref;
    real rho = rho_ref;

    if (add_perturbation) {
      T += dT;
      rho += drho;
    }

    return rho;
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;

    real delta = g / (Rd * T_ref);
    real rho_s = p_s / (Rd * T_ref);

    real rho_ref = isothermal_zdep(z, rho_s, T_ref, g, thermo);

    real dT_b = dT_max * exp(-pow((x - x_c) / d, 2)) * std::sin(pi * z / Lz);
    real dT = exp(delta * z / 2) * dT_b;

    real drho_b = -rho_s * dT_b / T_ref;
    real drho = exp(-delta * z / 2) * drho_b;

    real T = T_ref;
    if (add_perturbation) {
      T += dT;
    }
    real rho = rho_ref;

    // real p = rho * Rd * T;
    real p_ref = isothermal_zdep(z, p_s, T_ref, g, thermo);
    real dp = Rd * T_ref * drho + Rd * rho_ref * dT;
    real p = p_ref;

    if (add_perturbation) {
      rho += drho;
      p += dp;
    }

    return thermo.compute_entropic_var_from_p_T(p, T, 1, 0, 0, 0);
  }

  static real YAKL_INLINE entropicdensity_f(real x, real y, real z,
                                            const ThermoPotential &thermo) {
    return rho_f(x, y, z, thermo) * entropicvar_f(x, y, z, thermo);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = u_0;
    vvec.w = 0;
    return vvec;
  }

  static real YAKL_INLINE rhoexact_f(real x, real y, real z, real t,
                                     const ThermoPotential &thermo) {
    real rho = refrho_f(z, thermo);
    if (add_perturbation) {
      rho += sum_series(x, z, t, thermo).drho;
    }
    return rho;
  }

  static real YAKL_INLINE entropicdensityexact_f(
      real x, real y, real z, real t, const ThermoPotential &thermo) {

    const auto sol = sum_series(x, z, t, thermo);

    real p_ref = isothermal_zdep(z, p_s, T_ref, g, thermo);
    real rho_ref = refrho_f(z, thermo);

    real rho = rho_ref;
    real p = p_ref;
    real T = T_ref;

    if (add_perturbation) {
      rho += sol.drho;
      p += sol.dp;
      T += sol.dT;
    }

    return rho * thermo.compute_entropic_var_from_p_T(p, T, 1, 0, 0, 0);
  }

  static real YAKL_INLINE Texact_f(real x, real y, real z, real t,
                                   const ThermoPotential &thermo) {
    real T = T_ref;
    if (add_perturbation) {
      T += sum_series(x, z, t, thermo).dT;
    }
    return T;
  }

  static VecXYZ YAKL_INLINE vexact_f(real x, real y, real z, real t,
                                     const ThermoPotential &thermo) {
    VecXYZ v;

    const auto sol = sum_series(x, z, t, thermo);
    v.u = u_0;
    v.w = 0;
    if (add_perturbation) {
      v.u += sol.du;
      v.w += sol.dw;
    }
    return v;
  }

  struct SeriesSolution {
    real drho;
    real dp;
    real dT;
    real du;
    real dv;
    real dw;
  };

  static SeriesSolution YAKL_INLINE sum_series(real x, real z, real t,
                                               const ThermoPotential &thermo) {
    SeriesSolution sol;

    real x_c = GravityWave::x_c;
    real Lz = GravityWave::Lz;
    real g = GravityWave::g;
    real p_s = GravityWave::p_s;
    real T_ref = GravityWave::T_ref;

    complex im(0, 1);

    real Rd = thermo.cst.Rd;
    real cv_d = thermo.cst.Cvd;
    real cp_d = thermo.cst.Cpd;

    real xp = x - u_0 * t;
    real delta = g / (Rd * T_ref);
    real delta2 = delta * delta;
    real delta4 = delta2 * delta2;
    real c_s = sqrt(cp_d / cv_d * Rd * T_ref);
    real c_s2 = c_s * c_s;
    real rho_s = p_s / (Rd * T_ref);
    real f = 0;
    real f2 = f * f;

    complex drho_b = 0;
    complex dp_b = 0;
    complex du_b = 0;
    complex dv_b = 0;
    complex dw_b = 0;
    for (int m = -1; m < 2; m += 2) {
      for (int n = -100; n < 101; n++) {
        real k_x = 2 * pi * n / Lx;
        real k_z = pi * m / Lz;

        real k_x2 = k_x * k_x;
        real k_z2 = k_z * k_z;

        real p_1 = c_s2 * (k_x2 + k_z2 + delta2 / 4) + f2;
        real q_1 =
            g * k_x2 * (c_s2 * delta - g) + c_s2 * f2 * (k_z2 + delta2 / 4);

        real alpha = sqrt(p_1 / 2 - sqrt(p_1 * p_1 / 4 - q_1));
        real alpha2 = alpha * alpha;
        real beta = sqrt(p_1 / 2 + sqrt(p_1 * p_1 / 4 - q_1));
        real beta2 = beta * beta;

        real fac1 = 1 / (beta * beta - alpha * alpha);
        real L_m1 =
            (-std::cos(alpha * t) / alpha2 + std::cos(beta * t) / beta2) *
                fac1 +
            1 / (alpha2 * beta2);
        real L_0 =
            (std::sin(alpha * t) / alpha - std::sin(beta * t) / beta) * fac1;
        real L_1 = (std::cos(alpha * t) - std::cos(beta * t)) * fac1;
        real L_2 =
            (-alpha * std::sin(alpha * t) + beta * std::sin(beta * t)) * fac1;
        real L_3 =
            (-alpha2 * std::cos(alpha * t) + beta2 * std::cos(beta * t)) * fac1;

        if (alpha == 0) {
          L_m1 = (beta2 * t * t - 1 + std::cos(beta * t)) / (beta2 * beta2);
          L_0 = (beta * t - std::sin(beta * t)) / (beta2 * beta);
        }

        complex drhot_b0 = -rho_s / T_ref * dT_max / sqrt(pi) * d / Lx *
                           exp(-d * d * k_x2 / 4) * exp(-im * k_x * x_c) * k_z *
                           Lz / (2._fp * im);

        complex drhot_b =
            (L_3 + (p_1 + g * (im * k_z - delta / 2)) * L_1 +
             (c_s2 * (k_z2 + delta2 / 4) + g * (im * k_z - delta / 2)) * f2 *
                 L_m1) *
            drhot_b0;

        complex dpt_b = -(g - c_s2 * (im * k_z + delta / 2)) *
                        (L_1 + f2 * L_m1) * g * drhot_b0;

        complex dut_b = im * k_x * (g - c_s2 * (im * k_z + delta / 2)) * L_0 *
                        g * drhot_b0 / rho_s;

        complex dvt_b = -f * im * k_x * (g - c_s2 * (im * k_z + delta / 2)) *
                        L_m1 * g * drhot_b0 / rho_s;

        complex dwt_b =
            -(L_2 + (f2 + c_s2 * k_x2) * L_0) * g * drhot_b0 / rho_s;

        complex expfac = exp(im * (k_x * xp + k_z * z));

        drho_b += drhot_b * expfac;
        dp_b += dpt_b * expfac;
        du_b += dut_b * expfac;
        dv_b += dvt_b * expfac;
        dw_b += dwt_b * expfac;
      }
    }
    complex dT_b = T_ref * (dp_b / p_s - drho_b / rho_s);

    sol.drho = exp(-delta * z / 2) * drho_b.real();
    sol.dp = exp(-delta * z / 2) * dp_b.real();
    sol.dT = exp(delta * z / 2) * dT_b.real();
    sol.du = exp(delta * z / 2) * du_b.real();
    sol.dv = exp(delta * z / 2) * dv_b.real();
    sol.dw = exp(delta * z / 2) * dw_b.real();

    return sol;
  }

  struct ExactDensityDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "dense";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return rhoexact_f(x, y, z, time, thermo);
          },
          field, 0);

      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return entropicdensityexact_f(x, y, z, time, thermo);
          },
          field, 1);
    }
  };

  struct ExactTemperatureDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "Te";
      topology = pgeom.topology;
      dofs_arr = {0, 0, 1};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      dual_geometry.set_00form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return Texact_f(x, y, z, time, thermo);
          },
          field, 0);
    }
  };

  struct ExactWDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "we";
      topology = pgeom.topology;
      dofs_arr = {0, 1, 1};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      primal_geometry.set_01form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return vexact_f(x, y, z, time, thermo);
          },
          field, 0);
    }
  };

  struct BackgroundDensityDiagnostic : public Diagnostic {

    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "densb";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) { return refrho_f(z, thermo); },
          field, 0);

      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return refentropicdensity_f(z, thermo);
          },
          field, 1);
    }
  };

  struct TemperatureDiagnostic : public Diagnostic {
    static constexpr bool linear_T = true;
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "T";
      topology = pgeom.topology;
      dofs_arr = {0, 0, 1};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      const auto &primal_topology = primal_geometry.topology;

      int pis = primal_topology.is;
      int pjs = primal_topology.js;
      int pks = primal_topology.ks;

      const auto &densvar = x.fields_arr[DENSVAR].data;
      const auto &refstate = this->equations->reference_state;
      const auto &refdens = refstate.dens.data;

      YAKL_SCOPE(thermo, equations->thermo);
      YAKL_SCOPE(varset, equations->varset);
      parallel_for(
          "Compute T diagnostic",
          SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
            real alpha = varset.get_alpha(densvar, k, j, i, pks, pjs, pis, n);
            real entropic_var =
                varset.get_entropic_var(densvar, k, j, i, pks, pjs, pis, n);

            real T;
            if (!linear_T) {
              T = thermo.compute_T_from_alpha(alpha, entropic_var, 1, 0, 0, 0);
            } else {
              real Rd = thermo.cst.Rd;
              real refalpha = varset.get_alpha(refdens, k, pks, n);
              real refentropic_var =
                  varset.get_entropic_var(refdens, k, pks, n);

              real rho = 1 / alpha;
              real refrho = 1 / refalpha;

              real p = thermo.solve_p(rho, entropic_var, 1, 0, 0, 0);
              real refp = thermo.solve_p(refrho, refentropic_var, 1, 0, 0, 0);

              real drho = rho - refrho;
              real dp = p - refp;

              T = T_ref + dp / (refrho * Rd) -
                  refp / (refrho * refrho * Rd) * drho;
            }

            field.data(0, k + pks, j + pjs, i + pis, n) = T;
          });
    }
  };

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
    diagnostics.emplace_back(std::make_unique<ExactDensityDiagnostic>());
    diagnostics.emplace_back(std::make_unique<ExactWDiagnostic>());
    diagnostics.emplace_back(std::make_unique<ExactTemperatureDiagnostic>());
    diagnostics.emplace_back(std::make_unique<BackgroundDensityDiagnostic>());
    diagnostics.emplace_back(std::make_unique<TemperatureDiagnostic>());
  }
};

struct Supercell : TestCaseInit {
  static int constexpr max_ndims = 2;
  static bool constexpr needs_special_init = true;
  static real constexpr g = 9.81;
  static real constexpr Lx = 168e3;
  static real constexpr Ly = 168e3;
  static real constexpr Lz = 20e3;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr yc = 0.5_fp * Ly;
  static real constexpr zc = 0.5_fp * Lz;

  static real constexpr xbc = 0.5_fp * Lx;
  static real constexpr ybc = 0.5_fp * Ly;
  static real constexpr zbc = 1.5e3;

  static real constexpr rx = 10e3;
  static real constexpr ry = 10e3;
  static real constexpr rz = 1.5e3;

  static real constexpr dtht = 3;
  static real constexpr tht_0 = 300;
  static real constexpr z_tr = 12e3;
  static real constexpr tht_tr = 343;
  static real constexpr T_tr = 213;

  static real constexpr z_s = 5e3;
  static real constexpr U_s = 30;
  static real constexpr U_c = 15;
  static real constexpr dz_u = 1e3;

  static real constexpr N_ref = 0.011;

  static real YAKL_INLINE refnsq_f(real z, const ThermoPotential &thermo) {
    return N_ref * N_ref;
  }

  static real YAKL_INLINE refp_f(real z, const ThermoPotential &thermo) {
    return const_stability_p(z, N_ref, g, thermo.cst.pr, tht_0, thermo);
  }

  static real YAKL_INLINE refT_f(real z, const ThermoPotential &thermo) {
    return const_stability_T(z, N_ref, g, tht_0, thermo);
  }

  static real YAKL_INLINE refrho_f(real z, const ThermoPotential &thermo) {
    real p = refp_f(z, thermo);
    real T = refT_f(z, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE refentropicdensity_f(real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(z, thermo);
    real T_ref = refT_f(z, thermo);
    real p_ref = refp_f(z, thermo);
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE refrhov_f(real z, const ThermoPotential &thermo) {
    return 0;
  }

  static real YAKL_INLINE tht_f(real z, const ThermoPotential &thermo) {
    const real Cpd = thermo.cst.Cpd;
    if (z <= z_tr) {
      return tht_0 + (tht_tr - tht_0) * pow(z / z_tr, 5. / 4.);
    } else {
      return tht_tr * exp(g / (Cpd * T_tr) * (z - z_tr));
    }
  }

  static real YAKL_INLINE hum_f(real z, const ThermoPotential &thermo) {
    if (z <= z_tr) {
      return 1 - 0.75 * pow(z / z_tr, 5. / 4.);
    } else {
      return 0.25;
    }
  }

  static real YAKL_INLINE tht_perturb_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real dx = (x - xbc) / rx;
    real dy = (y - ybc) / ry;
    real dz = (z - zbc) / rz;
    real r = sqrt(dx * dx + dy * dy + dz * dz);
    return r < 1 ? dtht * pow(cos(pi * r / 2), 2) : 0;
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;

    real u_e;
    if (z < z_s - dz_u) {
      u_e = U_s * z / z_s - U_c;
    } else if (abs(z - z_s) <= dz_u) {
      u_e = (-4. / 5 + 3 * z / z_s - 5. / 4 * pow(z / z_s, 2)) * U_s - U_c;
    } else if (z > z_s + dz_u) {
      u_e = U_s - U_c;
    }
    vvec.u = u_e;
    return vvec;
  }

  static void initialize_refstate(Equations *equations,
                                  const Geometry<Straight> &primal_geometry,
                                  const Geometry<Twisted> &dual_geometry) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    real2d thtv("thtv", primal_topology.ni + 2 * primal_topology.mirror_halo,
                primal_topology.nens);
    real2d qv("qv", primal_topology.ni + 2 * primal_topology.mirror_halo,
              primal_topology.nens);
    real2d exner("exner", primal_topology.ni + 2 * primal_topology.mirror_halo,
                 primal_topology.nens);

    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(Hs, equations->Hs);
    YAKL_SCOPE(varset, equations->varset);
    YAKL_SCOPE(refstate, equations->reference_state);

    thermo.cst.Rd = 287.;
    thermo.cst.Rv = 461;
    thermo.cst.pr = 1e5;
    thermo.cst.Cpd = 1003;
    thermo.cst.Cvd = thermo.cst.Cpd - thermo.cst.Rd;
    thermo.cst.Cpv = 1859;

    thermo.cst.gamma_d = thermo.cst.Cpd / thermo.cst.Cvd;
    thermo.cst.kappa_d = thermo.cst.Rd / thermo.cst.Cpd;
    thermo.cst.delta_d = thermo.cst.Rd / thermo.cst.Cvd;

    const real pr = thermo.cst.pr;
    const real Cpd = thermo.cst.Cpd;
    const real kappa_d = thermo.cst.kappa_d;

    const int nonlinear_iters = 10;
    const real max_qv = 0.014;
    const real veps = thermo.cst.Rv / thermo.cst.Rd - 1;

    parallel_for(
        "Supercell set densities", primal_topology.nens, YAKL_LAMBDA(int n) {
          // set intial thtv to tht
          for (int k = 0; k < primal_topology.ni; ++k) {
            real z = primal_geometry.zint(k + pks, n);
            thtv(k + pks, n) = tht_f(z, thermo);
          }

          // iterate to solve nonlinear problem
          for (int m = 0; m < nonlinear_iters; ++m) {
            // compute exner
            exner(0 + pks, n) = 1;
            for (int k = 1; k < primal_topology.ni; ++k) {
              // TODO: change to phi
              real dz = primal_geometry.dz(k + pks, n);
              exner(k + pks, n) =
                  exner(k - 1 + pks, n) -
                  g / (Cpd * (thtv(k - 1 + pks, n) + thtv(k + pks, n)) / 2) *
                      dz;
            }

            // update thtv
            for (int k = 0; k < primal_topology.ni; ++k) {
              real z = primal_geometry.zint(k + pks, n);
              real p = pr * std::pow(exner(k + pks, n), 1. / kappa_d);
              real T = tht_f(z, thermo) * exner(k + pks, n);
              real qvs = saturation_mixing_ratio(T, p);
              qv(k + pks, n) = std::min(qvs * hum_f(z, thermo), max_qv);
              real tht = tht_f(z, thermo);
              thtv(k + pks, n) = tht * (1 + veps * qv(k + pks, n));
            }
          }
        });

    parallel_for(
        "Supercell set densitiy profiles",
        SimpleBounds<2>(dual_topology.nl, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          real dual_volume =
              dual_geometry.get_area_n1entity(k + dks, djs, dis, n);
          real z = primal_geometry.zint(k + pks, n);

          real p = pr * std::pow(exner(k + pks, n), 1. / kappa_d);
          real rho = p / (thermo.cst.Rd * exner(k + pks, n) * thtv(k + pks, n));

          refstate.dens.data(varset.dens_id_mass, k + dks, n) =
              rho * dual_volume;
          refstate.dens.data(varset.dens_id_entr, k + dks, n) =
              rho * thtv(k + pks, n) * dual_volume;
          refstate.dens.data(varset.dens_id_vap, k + dks, n) =
              rho * qv(k + pks, n) * dual_volume;
        });

    primal_geometry.set_profile_10form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        refstate.v, 0);
  }

  static void initialize(Equations *equations, FieldSet<nprognostic> &progvars,
                         FieldSet<nconstant> &constvars,
                         const Geometry<Straight> &primal_geometry,
                         const Geometry<Twisted> &dual_geometry) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(varset, equations->varset);
    YAKL_SCOPE(refdens, equations->reference_state.dens.data);

    auto densvar = progvars.fields_arr[DENSVAR].data;
    parallel_for(
        "Supercell initial condition",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
#ifndef PAMC_MAN
          densvar(varset.dens_id_mass, k + dks, j + djs, i + dis, n) =
              refdens(varset.dens_id_mass, k + dks, n);
#endif
          densvar(varset.dens_id_entr, k + dks, j + djs, i + dis, n) =
              refdens(varset.dens_id_entr, k + dks, n);

          // add perturbation
          CoordsXYZ ll_corner;
          dual_geometry.get_ll_corner(ll_corner, k, j, i, n);
          const real x = ll_corner.x + dual_geometry.dx / 2;
          const real y = ll_corner.y + dual_geometry.dy / 2;
          const real z = ll_corner.z + dual_geometry.dz(k + dks, n) / 2;
          densvar(varset.dens_id_entr, k + dks, j + djs, i + dis, n) +=
              tht_perturb_f(x, y, z, thermo) *
              refdens(varset.dens_id_mass, k + dks, n);

          densvar(varset.dens_id_vap, k + dks, j + djs, i + dis, n) =
              refdens(varset.dens_id_vap, k + dks, n);
        });

    primal_geometry.set_10form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[VVAR], 0);
    primal_geometry.set_01form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[WVAR], 0);
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance) {
  if (name == "gravitywave") {
    testcase = std::make_unique<EulerTestCase<GravityWave<true>>>();
  } else if (name == "twobubbles") {
    testcase = std::make_unique<EulerTestCase<TwoBubbles>>();
  } else if (name == "risingbubble") {
    if (acoustic_balance) {
      testcase = std::make_unique<EulerTestCase<RisingBubble<true>>>();
    } else {
      testcase = std::make_unique<EulerTestCase<RisingBubble<false>>>();
    }
  } else if (name == "densitycurrent") {
    testcase = std::make_unique<EulerTestCase<DensityCurrent>>();
  } else if (name == "moistrisingbubble") {
    testcase = std::make_unique<MoistEulerTestCase<MoistRisingBubble>>();
  } else if (name == "supercell") {
    testcase = std::make_unique<MoistEulerTestCase<Supercell>>();
  } else if (name == "largerisingbubble") {
    testcase = std::make_unique<EulerTestCase<LargeRisingBubble>>();
  } else if (name == "moistlargerisingbubble") {
    testcase = std::make_unique<MoistEulerTestCase<MoistLargeRisingBubble>>();
  } else if (name == "doublevortex") {
    testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  } else {
    throw std::runtime_error("unknown test case");
  }
}

#ifdef PAM_STANDALONE
void testcase_from_config(std::unique_ptr<TestCase> &testcase,
                          const YAML::Node &config) {
  const std::string name = config["init_data"].as<std::string>();
  const bool acoustic_balance =
      config["balance_initial_density"].as<bool>(false);
  testcase_from_string(testcase, name, acoustic_balance);
}
#endif

std::unique_ptr<TestCase> make_coupled_test_case(PamCoupler &coupler) {
  return std::make_unique<CoupledTestCase>(coupler);
}
} // namespace pamc
