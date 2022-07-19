#pragma once

#include "common.h"
#include "model.h"
#include "stats.h"

#include "ext_deriv.h"
#include "fct.h"
#include "hamiltonian.h"
#include "hodge_star.h"
#include "recon.h"
#include "thermo.h"
#include "variableset.h"
#include "wedge.h"

// *******   Functionals/Hamiltonians   ***********//

Functional_PVPE_extruded PVPE;
Hamiltonian_Hk_extruded Hk;

VariableSet varset;
#ifdef _SWE
Hamiltonian_SWE_Hs Hs;
#elif _TSWE
Hamiltonian_TSWE_Hs Hs;
#elif _CE
Hamiltonian_CE_Hs Hs;
#elif _MCErho
Hamiltonian_MCE_Hs Hs;
#elif _MCErhod
Hamiltonian_MCE_Hs Hs;
#elif _CEp
Hamiltonian_CE_p_Hs Hs;
#elif _MCErhop
Hamiltonian_MCE_p_Hs Hs;
#elif _MCErhodp
Hamiltonian_MCE_p_Hs Hs;
#endif
// ADD ANELASTIC + MOIST ANELASTIC

ThermoPotential thermo;

// *******   Diagnostics   ***********//

class Dens0Diagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom) override {
    // concentration 0-forms for dens
    name = "densl";
    topology = pgeom.topology;
    dofs_arr = {0, 0, ndensity}; // densldiag = straight (0,0)-form
    Diagnostic::initialize(pgeom, dgeom);
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
          compute_Iext<ndensity, diff_ord, vert_diff_ord>(
              field.data, x.fields_arr[DENSVAR].data, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
        });
  }
};

class QXZ0Diagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom) override {
    name = "QXZl";
    topology = dgeom.topology;
    dofs_arr = {0, 0, 1}; // // Qldiag = twisted (0,0)-form
    Diagnostic::initialize(pgeom, dgeom);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(PVPE, ::PVPE);
    parallel_for(
        "Compute Q0 DIAG",
        SimpleBounds<4>(dual_topology.ni - 3, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_qxz0(field.data, x.fields_arr[VVAR].data,
                            x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data,
                            const_vars.fields_arr[CORIOLISXZVAR].data, dis, djs,
                            dks, i, j, k + 2, n);
        });

    // Bottom is k=1 and top is k=dual_topology.ni-2
    parallel_for(
        "Compute Q0 DIAG TOP/BOTTOM",
        SimpleBounds<3>(dual_topology.n_cells_y, dual_topology.n_cells_x,
                        dual_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          PVPE.compute_qxz0_bottom(field.data, x.fields_arr[VVAR].data,
                                   x.fields_arr[WVAR].data,
                                   x.fields_arr[DENSVAR].data,
                                   const_vars.fields_arr[CORIOLISXZVAR].data,
                                   dis, djs, dks, i, j, 1, n);
          PVPE.compute_qxz0_top(field.data, x.fields_arr[VVAR].data,
                                x.fields_arr[WVAR].data,
                                x.fields_arr[DENSVAR].data,
                                const_vars.fields_arr[CORIOLISXZVAR].data, dis,
                                djs, dks, i, j, dual_topology.ni - 2, n);
        });

    field.set_bnd(0.0);
  }
};

void add_model_diagnostics(
    std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
  diagnostics.emplace_back(std::make_unique<Dens0Diagnostic>());
  diagnostics.emplace_back(std::make_unique<QXZ0Diagnostic>());
}

// *******   Tendencies   ***********//

class ModelTendencies : public ExtrudedTendencies {
public:
  void initialize(PamCoupler &coupler, ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nconstant> &const_exchange) {

    ExtrudedTendencies::initialize(params, primal_geom, dual_geom, aux_exchange,
                                   const_exchange);
    varset.initialize(coupler, params, thermo, this->primal_geometry,
                      this->dual_geometry);
    PVPE.initialize(varset);
    Hk.initialize(varset, this->primal_geometry, this->dual_geometry);
    Hs.initialize(thermo, varset, this->primal_geometry, this->dual_geometry);
  }

  void
  convert_dynamics_to_coupler_state(PamCoupler &coupler,
                                    const FieldSet<nprognostic> &prog_vars,
                                    const FieldSet<nconstant> &const_vars) {
    varset.convert_dynamics_to_coupler_state(coupler, prog_vars, const_vars);
  }
  void
  convert_coupler_to_dynamics_state(PamCoupler &coupler,
                                    FieldSet<nprognostic> &prog_vars,
                                    const FieldSet<nconstant> &const_vars) {
    varset.convert_coupler_to_dynamics_state(coupler, prog_vars, const_vars);
  }

  void compute_constants(FieldSet<nconstant> &const_vars,
                         FieldSet<nprognostic> &x) {}

  void compute_functional_derivatives_and_diagnostic_quantities_I(
      real5d Uvar, real5d UWvar, real5d qxz0var, real5d fxz0var,
      real5d dens0var, const real5d Vvar, const real5d Wvar,
      const real5d densvar, const real5d coriolisxzvar) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute Dens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Iext<ndensity, diff_ord, vert_diff_ord>(
              dens0var, densvar, this->primal_geometry, this->dual_geometry,
              pis, pjs, pks, i, j, k, n);
        });

    parallel_for(
        "Compute Uvar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Hext<1, diff_ord>(Uvar, Vvar, this->primal_geometry,
                                    this->dual_geometry, dis, djs, dks, i, j, k,
                                    n);
        });

    parallel_for(
        "Compute UWVar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Hv<1, vert_diff_ord>(UWvar, Wvar, this->primal_geometry,
                                       this->dual_geometry, dis, djs, dks, i, j,
                                       k + 1, n);
        });

    YAKL_SCOPE(PVPE, ::PVPE);
    parallel_for(
        "Compute Q0, F0",
        SimpleBounds<4>(dual_topology.ni - 3, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_qxz0fxz0(qxz0var, fxz0var, Vvar, Wvar, densvar,
                                coriolisxzvar, dis, djs, dks, i, j, k + 2, n);
        });
    parallel_for(
        "Compute Q0, F0 bnd",
        SimpleBounds<3>(dual_topology.n_cells_y, dual_topology.n_cells_x,
                        dual_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          PVPE.compute_qxz0fxz0_bottom(qxz0var, fxz0var, Vvar, Wvar, densvar,
                                       coriolisxzvar, dis, djs, dks, i, j, 1,
                                       n);
          PVPE.compute_qxz0fxz0_top(qxz0var, fxz0var, Vvar, Wvar, densvar,
                                    coriolisxzvar, dis, djs, dks, i, j,
                                    dual_topology.ni - 2, n);
        });
  }

  void compute_functional_derivatives_and_diagnostic_quantities_II(
      real5d Fvar, real5d FWvar, real5d Kvar, real5d HEvar, real5d HEWvar,
      const real5d Vvar, const real5d Uvar, const real5d Wvar,
      const real5d UWvar, const real5d dens0var) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    // THIS WILL NEED SOME SLIGHT MODIFICATIONS FOR CASE OF NON-ZERO UWVAR_B IE
    // BOUNDARY FLUXES BUT FOR NOW IT IS FINE SINCE UWVAR=0 on BND AND THEREFORE
    // K COMPUTATIONS IGNORE IT
    YAKL_SCOPE(Hk, ::Hk);
    parallel_for(
        "Compute Fvar, Kvar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F(Fvar, HEvar, Uvar, dens0var, dis, djs, dks, i, j, k, n);
          Hk.compute_K(Kvar, Vvar, Uvar, Wvar, UWvar, dis, djs, dks, i, j, k,
                       n);
        });
    parallel_for(
        "Compute FWvar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_Fw(FWvar, HEWvar, UWvar, dens0var, dis, djs, dks, i, j,
                        k + 1, n);
        });
  }

  void compute_functional_derivatives_and_diagnostic_quantities_III(
      real5d Bvar, real5d FTvar, real5d FTWvar, const real5d Fvar,
      const real5d Uvar, const real5d FWvar, const real5d UWvar,
      const real5d Kvar, const real5d densvar, const real5d HSvar) {

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
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
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
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          compute_Wxz_w_bottom(FTWvar, Fvar, pis, pjs, pks, i, j, 0, n);
          compute_Wxz_w_top(FTWvar, Fvar, pis, pjs, pks, i, j,
                            primal_topology.nl - 1, n);
        });

    YAKL_SCOPE(Hk, ::Hk);
    YAKL_SCOPE(Hs, ::Hs);
    parallel_for(
        "Compute Bvar",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hs.compute_dHsdx(Bvar, densvar, HSvar, pis, pjs, pks, i, j, k, n);
          Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j, k,
                                            n);
        });
  }

  void compute_edge_reconstructions(
      real5d densedgereconvar, real5d densvertedgereconvar,
      real5d qxzedgereconvar, real5d qxzvertedgereconvar,
      real5d coriolisxzedgereconvar, real5d coriolisxzvertedgereconvar,
      const real5d dens0var, const real5d qxz0var, const real5d fxz0var) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "ComputeDensEdgeRecon",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_edge_recon<ndensity, dual_reconstruction_type,
                                     dual_reconstruction_order>(
              densedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
          compute_twisted_vert_edge_recon<ndensity,
                                          dual_vert_reconstruction_type,
                                          dual_vert_reconstruction_order>(
              densvertedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl,
              dual_vert_wenoSigma);
        });

    parallel_for(
        "ComputeQEdgeRecon",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_straight_xz_edge_recon<1, reconstruction_type,
                                         reconstruction_order>(
              qxzedgereconvar, qxz0var, pis, pjs, pks, i, j, k, n,
              primal_wenoRecon, primal_to_gll, primal_wenoIdl,
              primal_wenoSigma);
          compute_straight_xz_vert_edge_recon<1, vert_reconstruction_type,
                                              vert_reconstruction_order>(
              qxzvertedgereconvar, qxz0var, pis, pjs, pks, i, j, k, n,
              primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl,
              primal_vert_wenoSigma);
          compute_straight_xz_edge_recon<1, coriolis_reconstruction_type,
                                         coriolis_reconstruction_order>(
              coriolisxzedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
              coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
              coriolis_wenoSigma);
          compute_straight_xz_vert_edge_recon<
              1, coriolis_vert_reconstruction_type,
              coriolis_vert_reconstruction_order>(
              coriolisxzvertedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
              coriolis_vert_wenoRecon, coriolis_vert_to_gll,
              coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
        });
  }

  void compute_recons(
      real5d densreconvar, real5d densvertreconvar, real5d qxzreconvar,
      real5d qxzvertreconvar, real5d coriolisxzreconvar,
      real5d coriolisxzvertreconvar, const real5d densedgereconvar,
      const real5d densvertedgereconvar, const real5d qxzedgereconvar,
      const real5d qxzvertedgereconvar, const real5d coriolisxzedgereconvar,
      const real5d coriolisxzvertedgereconvar, const real5d HEvar,
      const real5d HEWvar, const real5d Uvar, const real5d UWvar,
      const real5d FTvar, const real5d FTWvar) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "ComputeDensRECON",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_recon<ndensity, dual_reconstruction_type>(
              densreconvar, densedgereconvar, Uvar, dis, djs, dks, i, j, k, n);
          // scale twisted recons
          for (int d = 0; d < ndims; d++) {
            for (int l = 0; l < ndensity; l++) {
              densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) =
                  densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) /
                  HEvar(d, k + dks, j + djs, i + dis, n);
            }
          }
        });

    parallel_for(
        "ComputeDensVertRECON",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_vert_recon<ndensity, dual_vert_reconstruction_type>(
              densvertreconvar, densvertedgereconvar, UWvar, dis, djs, dks, i,
              j, k + 1, n);
          // scale twisted recons
          for (int l = 0; l < ndensity; l++) {
            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) =
                densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) /
                HEWvar(0, k + dks + 1, j + djs, i + dis, n);
          }
        });

    parallel_for(
        "ComputeQRECON",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_straight_xz_recon<1, reconstruction_type>(
              qxzreconvar, qxzedgereconvar, FTWvar, pis, pjs, pks, i, j, k, n);
          compute_straight_xz_recon<1, coriolis_reconstruction_type>(
              coriolisxzreconvar, coriolisxzedgereconvar, FTWvar, pis, pjs, pks,
              i, j, k, n);
        });
    parallel_for(
        "ComputeQVERTRECON",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_straight_xz_vert_recon<1, vert_reconstruction_type>(
              qxzvertreconvar, qxzvertedgereconvar, FTvar, pis, pjs, pks, i, j,
              k, n);
          compute_straight_xz_vert_recon<1, coriolis_vert_reconstruction_type>(
              coriolisxzvertreconvar, coriolisxzvertedgereconvar, FTvar, pis,
              pjs, pks, i, j, k, n);
        });
  }

  void
  compute_tendencies(real5d denstendvar, real5d Vtendvar, real5d Wtendvar,
                     const real5d densreconvar, const real5d densvertreconvar,
                     const real5d qxzreconvar, const real5d qxzvertreconvar,
                     const real5d coriolisxzreconvar,
                     const real5d coriolisxzvertreconvar, const real5d Bvar,
                     const real5d Fvar, const real5d FWvar, const real5d Phivar,
                     const real5d Phivertvar) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute Wtend",
        SimpleBounds<4>(primal_topology.nl - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDv_fct<ndensity>(Wtendvar, densvertreconvar, Phivertvar,
                                    Bvar, pis, pjs, pks, i, j, k + 1, n);
          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, qxzreconvar,
                                               qxzvertreconvar, Fvar, pis, pjs,
                                               pks, i, j, k + 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_w_nonEC<1, ADD_MODE::ADD>(
                Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, k + 1, n);
          }
          compute_Qxz_w_EC<1, ADD_MODE::ADD>(Wtendvar, coriolisxzreconvar,
                                             coriolisxzvertreconvar, Fvar, pis,
                                             pjs, pks, i, j, k + 1, n);
        });
    parallel_for(
        "Compute Wtend Bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          compute_wDv_fct<ndensity>(Wtendvar, densvertreconvar, Phivertvar,
                                    Bvar, pis, pjs, pks, i, j, 0, n);
          compute_wDv_fct<ndensity>(Wtendvar, densvertreconvar, Phivertvar,
                                    Bvar, pis, pjs, pks, i, j,
                                    primal_topology.nl - 1, n);
          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(
                Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i,
                j, 0, n);
            compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(
                Wtendvar, qxzreconvar, qxzvertreconvar, Fvar, pis, pjs, pks, i,
                j, primal_topology.nl - 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_w_nonEC_bottom<1, ADD_MODE::ADD>(
                Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j, 0, n);
            compute_Qxz_w_nonEC_top<1, ADD_MODE::ADD>(
                Wtendvar, qxzreconvar, Fvar, pis, pjs, pks, i, j,
                primal_topology.nl - 1, n);
          }
          compute_Qxz_w_EC_bottom<1, ADD_MODE::ADD>(
              Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis,
              pjs, pks, i, j, 0, n);
          compute_Qxz_w_EC_top<1, ADD_MODE::ADD>(
              Wtendvar, coriolisxzreconvar, coriolisxzvertreconvar, Fvar, pis,
              pjs, pks, i, j, primal_topology.nl - 1, n);
        });

    parallel_for(
        "Compute Vtend",
        SimpleBounds<4>(primal_topology.ni - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD1_fct<ndensity>(Vtendvar, densreconvar, Phivar, Bvar, pis,
                                    pjs, pks, i, j, k + 1, n);
          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar, qxzreconvar,
                                               qxzvertreconvar, FWvar, pis, pjs,
                                               pks, i, j, k + 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_u_nonEC<1, ADD_MODE::ADD>(Vtendvar, qxzvertreconvar,
                                                  FWvar, pis, pjs, pks, i, j,
                                                  k + 1, n);
          }
          compute_Qxz_u_EC<1, ADD_MODE::ADD>(Vtendvar, coriolisxzreconvar,
                                             coriolisxzvertreconvar, FWvar, pis,
                                             pjs, pks, i, j, k + 1, n);
        });
    parallel_for(
        "Compute Vtend Bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          compute_wD1_fct<ndensity>(Vtendvar, densreconvar, Phivar, Bvar, pis,
                                    pjs, pks, i, j, 0, n);
          compute_wD1_fct<ndensity>(Vtendvar, densreconvar, Phivar, Bvar, pis,
                                    pjs, pks, i, j, primal_topology.ni - 1, n);
          if (qf_choice == QF_MODE::EC) {
            compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(
                Vtendvar, qxzreconvar, qxzvertreconvar, FWvar, pis, pjs, pks, i,
                j, 0, n);
            compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(
                Vtendvar, qxzreconvar, qxzvertreconvar, FWvar, pis, pjs, pks, i,
                j, primal_topology.ni - 1, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Qxz_u_nonEC_bottom<1, ADD_MODE::ADD>(
                Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j, 0, n);
            compute_Qxz_u_nonEC_top<1, ADD_MODE::ADD>(
                Vtendvar, qxzvertreconvar, FWvar, pis, pjs, pks, i, j,
                primal_topology.ni - 1, n);
          }
          compute_Qxz_u_EC_bottom<1, ADD_MODE::ADD>(
              Vtendvar, coriolisxzreconvar, coriolisxzvertreconvar, FWvar, pis,
              pjs, pks, i, j, 0, n);
          compute_Qxz_u_EC_top<1, ADD_MODE::ADD>(
              Vtendvar, coriolisxzreconvar, coriolisxzvertreconvar, FWvar, pis,
              pjs, pks, i, j, primal_topology.ni - 1, n);
        });

    parallel_for(
        "Compute Dens Tend",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDbar2_fct<ndensity>(denstendvar, densreconvar, Phivar, Fvar,
                                       dis, djs, dks, i, j, k, n);
          compute_wDvbar_fct<ndensity, ADD_MODE::ADD>(
              denstendvar, densvertreconvar, Phivertvar, FWvar, dis, djs, dks,
              i, j, k, n);
        });
  }

  void compute_rhs(real dt, FieldSet<nconstant> &const_vars,
                   FieldSet<nprognostic> &x,
                   FieldSet<nauxiliary> &auxiliary_vars,
                   FieldSet<nprognostic> &xtend) {

    const auto &dual_topology = dual_geometry.topology;

    // Compute U, W, q0, dens0
    compute_functional_derivatives_and_diagnostic_quantities_I(
        auxiliary_vars.fields_arr[UVAR].data,
        auxiliary_vars.fields_arr[UWVAR].data,
        auxiliary_vars.fields_arr[QXZ0VAR].data,
        auxiliary_vars.fields_arr[FXZ0VAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data, x.fields_arr[VVAR].data,
        x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data,
        const_vars.fields_arr[CORIOLISXZVAR].data);

    auxiliary_vars.fields_arr[QXZ0VAR].set_bnd(0.0);
    auxiliary_vars.fields_arr[FXZ0VAR].set_bnd(0.0);
    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
        auxiliary_vars.fields_arr[UWVAR]);
    this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
        auxiliary_vars.fields_arr[DENS0VAR]);
    this->aux_exchange->exchanges_arr[QXZ0VAR].exchange_field(
        auxiliary_vars.fields_arr[QXZ0VAR]);
    this->aux_exchange->exchanges_arr[FXZ0VAR].exchange_field(
        auxiliary_vars.fields_arr[FXZ0VAR]);

    // Compute K, F, FW, he, hew
    compute_functional_derivatives_and_diagnostic_quantities_II(
        auxiliary_vars.fields_arr[FVAR].data,
        auxiliary_vars.fields_arr[FWVAR].data,
        auxiliary_vars.fields_arr[KVAR].data,
        auxiliary_vars.fields_arr[HEVAR].data,
        auxiliary_vars.fields_arr[HEWVAR].data, x.fields_arr[VVAR].data,
        auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[WVAR].data,
        auxiliary_vars.fields_arr[UWVAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data);

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[FVAR].exchange_field(
        auxiliary_vars.fields_arr[FVAR]);
    this->aux_exchange->exchanges_arr[FWVAR].exchange_field(
        auxiliary_vars.fields_arr[FWVAR]);
    this->aux_exchange->exchanges_arr[KVAR].exchange_field(
        auxiliary_vars.fields_arr[KVAR]);
    this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
        auxiliary_vars.fields_arr[HEVAR]);
    this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(
        auxiliary_vars.fields_arr[HEWVAR]);

    // Compute FT, B
    compute_functional_derivatives_and_diagnostic_quantities_III(
        auxiliary_vars.fields_arr[BVAR].data,
        auxiliary_vars.fields_arr[FTVAR].data,
        auxiliary_vars.fields_arr[FTWVAR].data,
        auxiliary_vars.fields_arr[FVAR].data,
        auxiliary_vars.fields_arr[UVAR].data,
        auxiliary_vars.fields_arr[FWVAR].data,
        auxiliary_vars.fields_arr[UWVAR].data,
        auxiliary_vars.fields_arr[KVAR].data, x.fields_arr[DENSVAR].data,
        const_vars.fields_arr[HSVAR].data);
    // auxiliary_vars.fields_arr[KVAR].data,
    // auxiliary_vars.fields_arr[DENS0VAR].data,
    // const_vars.fields_arr[HSVAR].data);

    this->aux_exchange->exchanges_arr[BVAR].exchange_field(
        auxiliary_vars.fields_arr[BVAR]);
    this->aux_exchange->exchanges_arr[FTVAR].exchange_field(
        auxiliary_vars.fields_arr[FTVAR]);
    this->aux_exchange->exchanges_arr[FTWVAR].exchange_field(
        auxiliary_vars.fields_arr[FTWVAR]);

    // Compute densrecon, densvertrecon, qrecon and frecon
    compute_edge_reconstructions(
        auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[QXZEDGERECONVAR].data,
        auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data,
        auxiliary_vars.fields_arr[QXZ0VAR].data,
        auxiliary_vars.fields_arr[FXZ0VAR].data);

    this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(
        auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
    this->aux_exchange->exchanges_arr[DENSVERTEDGERECONVAR].exchange_field(
        auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR]);
    this->aux_exchange->exchanges_arr[QXZEDGERECONVAR].exchange_field(
        auxiliary_vars.fields_arr[QXZEDGERECONVAR]);
    this->aux_exchange->exchanges_arr[QXZEDGERECONVAR].exchange_field(
        auxiliary_vars.fields_arr[QXZEDGERECONVAR]);
    this->aux_exchange->exchanges_arr[CORIOLISXZVERTEDGERECONVAR]
        .exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR]);
    this->aux_exchange->exchanges_arr[CORIOLISXZVERTEDGERECONVAR]
        .exchange_field(auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR]);

    compute_recons(auxiliary_vars.fields_arr[DENSRECONVAR].data,
                   auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                   auxiliary_vars.fields_arr[QXZRECONVAR].data,
                   auxiliary_vars.fields_arr[QXZVERTRECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISXZRECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR].data,
                   auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[QXZEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[HEVAR].data,
                   auxiliary_vars.fields_arr[HEWVAR].data,
                   auxiliary_vars.fields_arr[UVAR].data,
                   auxiliary_vars.fields_arr[UWVAR].data,
                   auxiliary_vars.fields_arr[FTVAR].data,
                   auxiliary_vars.fields_arr[FTWVAR].data);

    this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[DENSRECONVAR]);
    this->aux_exchange->exchanges_arr[DENSVERTRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[DENSVERTRECONVAR]);
    this->aux_exchange->exchanges_arr[QXZRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[QXZRECONVAR]);
    this->aux_exchange->exchanges_arr[QXZVERTRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[QXZVERTRECONVAR]);
    this->aux_exchange->exchanges_arr[CORIOLISXZRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[CORIOLISXZRECONVAR]);
    this->aux_exchange->exchanges_arr[CORIOLISXZVERTRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR]);

    // Compute fct
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "ComputeEdgeFlux",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_edgefluxes<ndensity>(
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[DENSRECONVAR].data,
              auxiliary_vars.fields_arr[FVAR].data, dis, djs, dks, i, j, k, n);
        });
    parallel_for(
        "ComputeVertEdgeFlux",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_vertedgefluxes<ndensity>(
              auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
              auxiliary_vars.fields_arr[FWVAR].data, dis, djs, dks, i, j, k, n);
        });
    this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(
        auxiliary_vars.fields_arr[EDGEFLUXVAR]);
    this->aux_exchange->exchanges_arr[VERTEDGEFLUXVAR].exchange_field(
        auxiliary_vars.fields_arr[VERTEDGEFLUXVAR]);

    parallel_for(
        "ComputeMf",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Mfext<ndensity>(
              auxiliary_vars.fields_arr[MFVAR].data,
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, dt, dis, djs,
              dks, i, j, k, n);
        });

    this->aux_exchange->exchanges_arr[MFVAR].exchange_field(
        auxiliary_vars.fields_arr[MFVAR]);

    parallel_for(
        "ComputePhiVert",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Phivert<ndensity>(
              auxiliary_vars.fields_arr[PHIVERTVAR].data,
              auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[MFVAR].data, x.fields_arr[DENSVAR].data,
              dis, djs, dks, i, j, k + 1, n);
        });
    parallel_for(
        "ComputePhi",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Phi<ndensity>(auxiliary_vars.fields_arr[PHIVAR].data,
                                auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
                                auxiliary_vars.fields_arr[MFVAR].data,
                                x.fields_arr[DENSVAR].data, dis, djs, dks, i, j,
                                k, n);
        });

    // Don't do FCT for non-FCT vars
    for (int l = 0; l < ndensity; l++) {
      // if (not varset.dens_pos(l))
      if (!varset.dens_pos[l]) {
        auxiliary_vars.fields_arr[PHIVAR].set(l, 1.0);
        auxiliary_vars.fields_arr[PHIVERTVAR].set(l, 1.0);
      }
    }

    this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(
        auxiliary_vars.fields_arr[PHIVAR]);
    this->aux_exchange->exchanges_arr[PHIVERTVAR].exchange_field(
        auxiliary_vars.fields_arr[PHIVERTVAR]);

    // Compute tendencies
    compute_tendencies(xtend.fields_arr[DENSVAR].data,
                       xtend.fields_arr[VVAR].data, xtend.fields_arr[WVAR].data,
                       auxiliary_vars.fields_arr[DENSRECONVAR].data,
                       auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                       auxiliary_vars.fields_arr[QXZRECONVAR].data,
                       auxiliary_vars.fields_arr[QXZVERTRECONVAR].data,
                       auxiliary_vars.fields_arr[CORIOLISXZRECONVAR].data,
                       auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR].data,
                       auxiliary_vars.fields_arr[BVAR].data,
                       auxiliary_vars.fields_arr[FVAR].data,
                       auxiliary_vars.fields_arr[FWVAR].data,
                       auxiliary_vars.fields_arr[PHIVAR].data,
                       auxiliary_vars.fields_arr[PHIVERTVAR].data);
  }
};

// *******   Statistics   ***********//

class ModelStats : public Stats {
public:
  real3d TEarr, KEarr, PEarr, IEarr, PVarr, PENSarr, trimmed_density;

  void initialize(ModelParameters &params, Parallel &par,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    Stats::initialize(params, par, primal_geom, dual_geom);
    this->stats_arr[DENSSTAT].initialize("mass", ndensity, this->statsize,
                                         this->nens, this->masterproc);
    this->stats_arr[DENSMAXSTAT].initialize("densmax", ndensity, this->statsize,
                                            this->nens, this->masterproc);
    this->stats_arr[DENSMINSTAT].initialize("densmin", ndensity, this->statsize,
                                            this->nens, this->masterproc);
    this->stats_arr[ESTAT].initialize("energy", 4, this->statsize, this->nens,
                                      this->masterproc);
    this->stats_arr[PVSTAT].initialize("pv", 1, this->statsize, this->nens,
                                       this->masterproc);
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
    PVarr = real3d("PV", primal_topology.nl, primal_topology.n_cells_y,
                   primal_topology.n_cells_x);
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

      SArray<real, 1, ndensity> masslocal, massglobal;
      SArray<real, 1, ndensity> densmaxlocal, densmaxglobal;
      SArray<real, 1, ndensity> densminlocal, densminglobal;
      SArray<real, 1, 1> pvlocal, pvglobal;
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
      for (int l = 0; l < ndensity; l++) {
        masslocal(l) = 0.;
        massglobal(l) = 0.;
      }
      for (int l = 0; l < ndensity; l++) {
        densmaxlocal(l) = 0.;
        densmaxglobal(l) = 0.;
      }
      for (int l = 0; l < ndensity; l++) {
        densminlocal(l) = 0.;
        densminglobal(l) = 0.;
      }

      int dis = dual_topology.is;
      int djs = dual_topology.js;
      int dks = dual_topology.ks;

      YAKL_SCOPE(Hk, ::Hk);
      YAKL_SCOPE(Hs, ::Hs);
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

      YAKL_SCOPE(PVPE, ::PVPE);
      parallel_for(
          "Compute PV/PE stats",
          SimpleBounds<3>(primal_topology.nl - 2, primal_topology.n_cells_y,
                          primal_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int k, int j, int i) {
            pvpe vals_pvpe;
            vals_pvpe = PVPE.compute_PVPE(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j,
                k + 1, n);
            PVarr(k + 1, j, i) = vals_pvpe.pv;
            PENSarr(k + 1, j, i) = vals_pvpe.pe;
          });
      parallel_for(
          "Compute PV/PE stats bnd",
          SimpleBounds<2>(primal_topology.n_cells_y, primal_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int j, int i) {
            pvpe vals_pvpe;
            vals_pvpe = PVPE.compute_PVPE_bottom(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j,
                0, n);
            PVarr(0, j, i) = vals_pvpe.pv;
            PENSarr(0, j, i) = vals_pvpe.pe;
            vals_pvpe = PVPE.compute_PVPE_top(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISXZVAR].data, pis, pjs, pks, i, j,
                primal_topology.nl - 1, n);
            PVarr(primal_topology.nl - 1, j, i) = vals_pvpe.pv;
            PENSarr(primal_topology.nl - 1, j, i) = vals_pvpe.pe;
          });

      pvlocal(0) = yakl::intrinsics::sum(PVarr);
      pelocal(0) = yakl::intrinsics::sum(PENSarr);

      for (int l = 0; l < ndensity; l++) {
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
      this->ierr =
          MPI_Ireduce(&masslocal, &massglobal, ndensity, REAL_MPI, MPI_SUM, 0,
                      MPI_COMM_WORLD, &this->Req[DENSSTAT]);
      this->ierr =
          MPI_Ireduce(&densmaxlocal, &densmaxglobal, ndensity, REAL_MPI,
                      MPI_MAX, 0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
      this->ierr =
          MPI_Ireduce(&densminlocal, &densminglobal, ndensity, REAL_MPI,
                      MPI_MIN, 0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
      this->ierr = MPI_Ireduce(&pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0,
                               MPI_COMM_WORLD, &this->Req[PVSTAT]);
      this->ierr = MPI_Ireduce(&pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0,
                               MPI_COMM_WORLD, &this->Req[PESTAT]);
      this->ierr = MPI_Ireduce(&elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0,
                               MPI_COMM_WORLD, &this->Req[ESTAT]);

      this->ierr = MPI_Waitall(nstats, this->Req, this->Status);

      if (masterproc) {
        for (int l = 0; l < ndensity; l++) {
          this->stats_arr[DENSSTAT].data(l, tind, n) = massglobal(l);
          this->stats_arr[DENSMAXSTAT].data(l, tind, n) = densmaxglobal(l);
          this->stats_arr[DENSMINSTAT].data(l, tind, n) = densminglobal(l);
        }

        this->stats_arr[ESTAT].data(0, tind, n) = eglobal(0);
        this->stats_arr[ESTAT].data(1, tind, n) = eglobal(1);
        this->stats_arr[ESTAT].data(2, tind, n) = eglobal(2);
        this->stats_arr[ESTAT].data(3, tind, n) = eglobal(3);
        this->stats_arr[PVSTAT].data(0, tind, n) = pvglobal(0);
        this->stats_arr[PESTAT].data(0, tind, n) = peglobal(0);
      }
    }
  }
};

// *******   FieldSet Initialization   ***********//
void initialize_variables(const Topology &ptopo, const Topology &dtopo,
                          SArray<int, 2, nprognostic, 3> &prog_ndofs_arr,
                          SArray<int, 2, nconstant, 3> &const_ndofs_arr,
                          SArray<int, 2, nauxiliary, 3> &aux_ndofs_arr,
                          std::array<std::string, nprognostic> &prog_names_arr,
                          std::array<std::string, nconstant> &const_names_arr,
                          std::array<std::string, nauxiliary> &aux_names_arr,
                          std::array<Topology, nprognostic> &prog_topo_arr,
                          std::array<Topology, nconstant> &const_topo_arr,
                          std::array<Topology, nauxiliary> &aux_topo_arr) {

  // primal grid represents straight quantities, dual grid twisted quantities
  // ndims is the BASEDIM size!

  // v, w, dens
  prog_topo_arr[VVAR] = ptopo;
  prog_topo_arr[WVAR] = ptopo;
  prog_topo_arr[DENSVAR] = dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[WVAR] = "w";
  prog_names_arr[DENSVAR] = "dens";
  set_dofs_arr(prog_ndofs_arr, VVAR, 1, 0, 1); // v = straight (1,0)-form
  set_dofs_arr(prog_ndofs_arr, WVAR, 0, 1, 1); // w = straight (0,1)-form
  set_dofs_arr(prog_ndofs_arr, DENSVAR, ndims, 1,
               ndensity); // dens = twisted (n,1)-form

  // hs
  const_topo_arr[HSVAR] = dtopo;
  const_names_arr[HSVAR] = "hs";
  set_dofs_arr(const_ndofs_arr, HSVAR, ndims, 1, 1); // hs = twisted (n,1)-form
  const_topo_arr[CORIOLISXZVAR] = ptopo;
  const_names_arr[CORIOLISXZVAR] = "coriolisxz";
  set_dofs_arr(const_ndofs_arr, CORIOLISXZVAR, 1, 1,
               1); // f = straight (1,1)-form

  // functional derivatives = F, B, K, he, U
  aux_topo_arr[BVAR] = ptopo;
  aux_topo_arr[FVAR] = dtopo;
  aux_topo_arr[UVAR] = dtopo;
  aux_topo_arr[HEVAR] = dtopo;
  aux_topo_arr[FWVAR] = dtopo;
  aux_topo_arr[UWVAR] = dtopo;
  aux_topo_arr[HEWVAR] = dtopo;
  aux_topo_arr[KVAR] = dtopo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  aux_names_arr[FWVAR] = "Fw";
  aux_names_arr[UWVAR] = "Uw";
  aux_names_arr[HEWVAR] = "hew";
  set_dofs_arr(aux_ndofs_arr, BVAR, 0, 0, ndensity);  // B = straight (0,0)-form
  set_dofs_arr(aux_ndofs_arr, KVAR, ndims, 1, 1);     // K = twisted (n,1)-form
  set_dofs_arr(aux_ndofs_arr, FVAR, ndims - 1, 1, 1); // F = twisted
                                                      // (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, UVAR, ndims - 1, 1, 1); // U = twisted
                                                      // (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, HEVAR, ndims - 1, 1,
               1); // he lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, FWVAR, ndims, 0, 1); // Fw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, UWVAR, ndims, 0, 1); // Uw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, HEWVAR, ndims, 0,
               1); // hew lives on vert dual edges, associated with Fw

  // dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = dtopo;
  aux_topo_arr[DENSEDGERECONVAR] = dtopo;
  aux_topo_arr[DENSVERTRECONVAR] = dtopo;
  aux_topo_arr[DENSVERTEDGERECONVAR] = dtopo;
  aux_topo_arr[DENS0VAR] = ptopo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSVERTRECONVAR] = "densvertrecon";
  aux_names_arr[DENSVERTEDGERECONVAR] = "densvertedgerecon";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  set_dofs_arr(
      aux_ndofs_arr, DENSRECONVAR, ndims - 1, 1,
      ndensity); // densrecon lives on horiz dual edges, associated with F
  set_dofs_arr(
      aux_ndofs_arr, DENSEDGERECONVAR, ndims, 1,
      2 * ndims *
          ndensity); // densedgerecon lives on dual cells, associated with F
  set_dofs_arr(
      aux_ndofs_arr, DENSVERTRECONVAR, ndims, 0,
      ndensity); // densvertrecon lives on vert dual edges, associated with Fw
  set_dofs_arr(
      aux_ndofs_arr, DENSVERTEDGERECONVAR, ndims, 1,
      2 * ndensity); // densedgerecon lives on dual cells, associated with Fw
  set_dofs_arr(aux_ndofs_arr, DENS0VAR, 0, 0,
               ndensity); // dens0 = straight (0,0)-form

  // fct stuff- Phi, Mf, edgeflux
  aux_topo_arr[PHIVAR] = dtopo;
  aux_topo_arr[PHIVERTVAR] = dtopo;
  aux_topo_arr[MFVAR] = dtopo;
  aux_topo_arr[EDGEFLUXVAR] = dtopo;
  aux_topo_arr[VERTEDGEFLUXVAR] = dtopo;
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[PHIVERTVAR] = "PhiVert";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  aux_names_arr[VERTEDGEFLUXVAR] = "vertedgeflux";
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims - 1, 1, ndensity);
  set_dofs_arr(aux_ndofs_arr, PHIVERTVAR, ndims, 0, ndensity);
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ndensity);
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims - 1, 1, ndensity);
  set_dofs_arr(aux_ndofs_arr, VERTEDGEFLUXVAR, ndims, 0, ndensity);

  // Q stuff
  aux_topo_arr[QXZ0VAR] = dtopo;
  aux_names_arr[QXZ0VAR] = "QXZ0";
  set_dofs_arr(aux_ndofs_arr, QXZ0VAR, 0, 0, 1); // Q0 = twisted (0,0)-form
  aux_topo_arr[QXZRECONVAR] = ptopo;
  aux_topo_arr[QXZEDGERECONVAR] = ptopo;
  aux_topo_arr[QXZVERTRECONVAR] = ptopo;
  aux_topo_arr[QXZVERTEDGERECONVAR] = ptopo;
  aux_topo_arr[QXZFLUXVAR] = ptopo;
  aux_topo_arr[QXZVERTFLUXVAR] = ptopo;
  aux_names_arr[QXZRECONVAR] = "qxzrecon";
  aux_names_arr[QXZEDGERECONVAR] = "qxzedgerecon";
  aux_names_arr[QXZVERTRECONVAR] = "qxzvertrecon";
  aux_names_arr[QXZVERTEDGERECONVAR] = "qxzvertedgerecon";
  aux_names_arr[QXZFLUXVAR] = "qxzflux";
  aux_names_arr[QXZVERTFLUXVAR] = "qxzvertflux";
  set_dofs_arr(aux_ndofs_arr, QXZRECONVAR, 0, 1,
               1); // qxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(
      aux_ndofs_arr, QXZEDGERECONVAR, ndims, 1,
      2 * 1); // qxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(
      aux_ndofs_arr, QXZVERTRECONVAR, 1, 0,
      1); // qxzsvertrecon lives on horiz primal edges, associated with v
  set_dofs_arr(
      aux_ndofs_arr, QXZVERTEDGERECONVAR, ndims, 1,
      2 * 1); // qxzvertedgerecon lives on primal cells, associated with F/v
  set_dofs_arr(aux_ndofs_arr, QXZFLUXVAR, 0, 1,
               1); // qxzflux lives on vert primal edges, associated with w
  set_dofs_arr(aux_ndofs_arr, QXZVERTFLUXVAR, 1, 0,
               1); // qxzvertflux lives on horiz primal edges, associated with v

  aux_topo_arr[FTVAR] = ptopo;
  aux_topo_arr[FTWVAR] = ptopo;
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[FTWVAR] = "FTW";
  set_dofs_arr(aux_ndofs_arr, FTVAR, 1, 0,
               1); // FT = straight (1,0)-form ie Fw at v pts
  set_dofs_arr(aux_ndofs_arr, FTWVAR, 0, 1,
               1); // FTW = straight (0,1)-form ie F at w pts

  aux_topo_arr[FXZ0VAR] = dtopo;
  aux_topo_arr[CORIOLISXZRECONVAR] = ptopo;
  aux_topo_arr[CORIOLISXZEDGERECONVAR] = ptopo;
  aux_topo_arr[CORIOLISXZVERTRECONVAR] = ptopo;
  aux_topo_arr[CORIOLISXZVERTEDGERECONVAR] = ptopo;
  aux_names_arr[FXZ0VAR] = "fxz0";
  aux_names_arr[CORIOLISXZRECONVAR] = "coriolisxzrecon";
  aux_names_arr[CORIOLISXZEDGERECONVAR] = "coriolisxzedgerecon";
  aux_names_arr[CORIOLISXZVERTRECONVAR] = "coriolisxzvertrecon";
  aux_names_arr[CORIOLISXZVERTEDGERECONVAR] = "coriolisxzvertedgerecon";
  set_dofs_arr(aux_ndofs_arr, FXZ0VAR, 0, 0, 1); // fxz0 is a twisted 0,0 form
  set_dofs_arr(
      aux_ndofs_arr, CORIOLISXZRECONVAR, 0, 1,
      1); // coriolisxzrecon lives on vert primal edges, associated with w
  set_dofs_arr(
      aux_ndofs_arr, CORIOLISXZEDGERECONVAR, ndims, 1,
      2 * 1); // coriolisxzedgerecon lives on primal cells, associated with Fw/w
  set_dofs_arr(
      aux_ndofs_arr, CORIOLISXZVERTRECONVAR, 1, 0,
      1); // coriolisxzsvertrecon lives on horiz primal edges, associated with v
  set_dofs_arr(aux_ndofs_arr, CORIOLISXZVERTEDGERECONVAR, ndims, 1,
               2 * 1); // coriolisxzvertedgerecon lives on primal cells,
                       // associated with F/v

  // #if defined _AN || defined _MAN
  // aux_topo_arr[PVAR] = ptopo; //p = straight 0-form
  // aux_names_arr[PVAR] = "p";
  // set_dofs_arr(aux_ndofs_arr, PVAR, 0, 1, 1);  //p = straight 0-form
  // #endif
}

void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance);

void readModelParamsFile(std::string inFile, ModelParameters &params,
                         Parallel &par, int nz,
                         std::unique_ptr<TestCase> &testcase) {

  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  params.acoustic_balance = config["balance_initial_density"].as<bool>(false);

  // Read the data initialization options
  params.initdataStr = config["initData"].as<std::string>();
  testcase_from_string(testcase, params.initdataStr, params.acoustic_balance);

  serial_print("IC: " + params.initdataStr, par.masterproc);
  serial_print("acoustically balanced: " +
                   std::to_string(params.acoustic_balance),
               par.masterproc);

  for (int i = 0; i < ntracers_dycore; i++) {
    params.tracerdataStr[i] =
        config["initTracer" + std::to_string(i)].as<std::string>();
    params.dycore_tracerpos[i] =
        config["initTracerPos" + std::to_string(i)].as<bool>();
    serial_print("Dycore Tracer" + std::to_string(i) +
                     " IC: " + params.tracerdataStr[i],
                 par.masterproc);
  }

  testcase->set_tracers(params);

  params.ylen = 1.0;
  params.yc = 0.5;
  testcase->set_domain(params);

  readParamsFile(inFile, params, par, nz);
}

//***************** Test Cases ***************************//

// Universal

real YAKL_INLINE isentropic_T(real x, real z, real theta0, real g,
                              const ThermoPotential &thermo) {
  return theta0 - z * g / thermo.cst.Cpd;
}

real YAKL_INLINE isentropic_p(real x, real z, real theta0, real g,
                              const ThermoPotential &thermo) {
  return thermo.cst.pr * pow(isentropic_T(x, z, theta0, g, thermo) / theta0,
                             1. / thermo.cst.kappa_d);
}

real YAKL_INLINE isentropic_rho(real x, real z, real theta0, real g,
                                const ThermoPotential &thermo) {
  real p = isentropic_p(x, z, theta0, g, thermo);
  real T = isentropic_T(x, z, theta0, g, thermo);
  real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
  return 1._fp / alpha;
}

real YAKL_INLINE linear_ellipsoid(real x, real z, real x0, real z0, real xrad,
                                  real zrad, real amp) {
  real xn = (x - x0) / xrad;
  real zn = (z - z0) / zrad;
  real dist = sqrt(xn * xn + zn * zn);
  return amp * std::max(1._fp - dist, 0._fp);
}

real YAKL_INLINE flat_geop(real x, real z, real g) { return g * z; }

// Returns saturation vapor pressure
real YAKL_INLINE saturation_vapor_pressure(real temp) {
  real tc = temp - 273.15_fp;
  return 610.94_fp * exp(17.625_fp * tc / (243.04_fp + tc));
}

template <class T> class SWETestCase : public TestCase, public T {
public:
  using T::g;

  using T::Lx;
  using T::Ly;
  using T::xc;
  using T::yc;

  using T::h_f;
#ifdef _TSWE
  using T::S_f;
#endif
  using T::coriolis_f;
  using T::v_f;

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.zlen = Ly;
    params.xc = xc;
    params.zc = yc;
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real y) { return h_f(x, y); },
        progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real y) { return S_f(x, y); },
        progvars.fields_arr[DENSVAR], 1);
#endif
    primal_geom.set_10form_values(
        YAKL_LAMBDA(real x, real y) { return v_f(x, y); },
        progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    primal_geom.set_01form_values(
        YAKL_LAMBDA(real x, real y) { return v_f(x, y); },
        progvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);

    // RESTORE ONCE PRIMAL GRID CFV RECON IS FIXED IN PARALLEL
    // primal_geom.set_11form_values(
    //      YAKL_LAMBDA(real x, real y) { return coriolis_f(x, y); },
    //      constvars.fields_arr[CORIOLISVAR], 0);

    YAKL_SCOPE(tracer_f, this->tracer_f);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_11form_values(
          YAKL_LAMBDA(real x, real y) {
            return h_f(x, y) * tracer_f(i)->compute(x, y, Lx, Ly, xc, yc);
          },
          progvars.fields_arr[DENSVAR], i + ndensity_dycore);
    }
    Hs.set_parameters(g);
  }
};

template <class T> class EulerTestCase : public TestCase, public T {
public:
  using T::g;
  using T::Lx;
  using T::Lz;
  using T::xc;
  using T::zc;

  using T::entropicvar_f;
  using T::rho_f;

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.zlen = Lz;
    params.xc = xc;
    params.zc = zc;
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    YAKL_SCOPE(thermo, ::thermo);
    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) { return rho_f(x, z, thermo); },
        progvars.fields_arr[DENSVAR], 0);
    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) {
          return rho_f(x, z, thermo) * entropicvar_f(x, z, thermo);
        },
        progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) { return flat_geop(x, z, g); },
        constvars.fields_arr[HSVAR], 0);

    YAKL_SCOPE(tracer_f, this->tracer_f);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return rho_f(x, z, thermo) *
                   tracer_f(i)->compute(x, z, Lx, Lz, xc, zc);
          },
          progvars.fields_arr[DENSVAR], i + ndensity_dycore);
    }
  }
};

template <class T> class MoistEulerTestCase : public TestCase, public T {
public:
  using T::g;
  using T::Lx;
  using T::Lz;
  using T::xc;
  using T::zc;

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
#if defined _MCErhod || defined _MCErhodp
    return T::rhod_f(x, z, thermo);
#else
    return T::rho_f(x, z, thermo);
#endif
  }
  using T::entropicdensity_f;

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.zlen = Lz;
    params.xc = xc;
    params.zc = zc;
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    YAKL_SCOPE(thermo, ::thermo);

    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) { return rho_f(x, z, thermo); },
        progvars.fields_arr[DENSVAR], 0);
    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) { return entropicdensity_f(x, z, thermo); },
        progvars.fields_arr[DENSVAR], 1);
    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) { return flat_geop(x, z, g); },
        constvars.fields_arr[HSVAR], 0);

    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) { return T::rhov_f(x, z, thermo); },
        progvars.fields_arr[DENSVAR], varset.dm_id_vap + ndensity_nophysics);

    YAKL_SCOPE(tracer_f, this->tracer_f);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return rho_f(x, z, thermo) *
                   tracer_f(i)->compute(x, z, Lx, Lz, xc, zc);
          },
          progvars.fields_arr[DENSVAR], i + ndensity_dycore);
    }
  }
};

struct DoubleVortex {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 5000000._fp;
  static real constexpr Ly = 5000000._fp;
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
  static real constexpr c = 0.05_fp;
  static real constexpr a = 1.0_fp / 3.0_fp;
  static real constexpr D = 0.5_fp * Lx;

  static real YAKL_INLINE coriolis_f(real x, real y) { return coriolis; }

  static real YAKL_INLINE h_f(real x, real y) {
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

  static vecext<2> YAKL_INLINE v_f(real x, real y) {
    vecext<2> vvec;

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

    vvec.u =
        -g * dh / coriolis / sigmay *
        (yprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         yprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));
    vvec.w =
        g * dh / coriolis / sigmax *
        (xprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         xprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));
    return vvec;
  }

  static real YAKL_INLINE S_f(real x, real y) {
    // real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. * M_PI
    // / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
    real sval =
        g * (1._fp + c * exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) /
                             (a * a * D * D)));
    // real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
    // real sval = g;
    // real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35 *
    // Ly
    // && y < 0.65 * Ly ) ? 1. : 0.));
    return sval * h_f(x, y);
  }
};

template <bool acoustic_balance> struct RisingBubble {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 1000._fp;
  static real constexpr Lz = 1500._fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr theta0 = 300.0_fp;
  static real constexpr bzc = 350._fp;
  static real constexpr dss = 0.5_fp;
  static real constexpr rc = 250._fp;
  static real constexpr rh0 = 0.8_fp;

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    real rho_b = isentropic_rho(x, z, theta0, g, thermo);
    if (acoustic_balance) {
      real theta = entropicvar_f(x, z, thermo);
      return rho_b * theta0 / theta;
    } else {
      return rho_b;
    }
  }

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(x, z, theta0, g, thermo);
    real T = isentropic_T(x, z, theta0, g, thermo);
    real r = sqrt((x - xc) * (x - xc) + (z - bzc) * (z - bzc));
    real dtheta = (r < rc) ? dss * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var(p, T + dT, 0, 0, 0, 0);
  }
};

struct MoistRisingBubble : public RisingBubble<false> {

  static real YAKL_INLINE rhov_f(real x, real z,
                                 const ThermoPotential &thermo) {
    real r = sqrt((x - xc) * (x - xc) + (z - bzc) * (z - bzc));
    real rh = (r < rc) ? rh0 * (1._fp + cos(pi * r / rc)) : 0._fp;
    real Th = isentropic_T(x, z, theta0, g, thermo);
    real svp = saturation_vapor_pressure(Th);
    real pv = svp * rh;
    return pv / (thermo.cst.Rv * Th);
  }

  static real YAKL_INLINE rhod_f(real x, real z,
                                 const ThermoPotential &thermo) {
    real p = isentropic_p(x, z, theta0, g, thermo);
    real T = isentropic_T(x, z, theta0, g, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    real rhod = rhod_f(x, z, thermo);
    real rhov = rhov_f(x, z, thermo);
    return rhod + rhov;
  }

  static real YAKL_INLINE entropicdensity_f(real x, real z,
                                            const ThermoPotential &thermo) {
    real p = isentropic_p(x, z, theta0, g, thermo);
    real T = isentropic_T(x, z, theta0, g, thermo);
    real r = sqrt((x - xc) * (x - xc) + (z - bzc) * (z - bzc));
    real dtheta = (r < rc) ? dss * 0.5_fp * (1. + cos(pi * r / rc)) : 0._fp;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    real theta = thermo.compute_entropic_var(p, T + dT, 1, 0, 0, 0);
    return theta * rho_f(x, z, thermo);
  }
};

struct LargeRisingBubble {
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

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    return isentropic_rho(x, z, theta0, g, thermo);
  }
  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {

    real p = isentropic_p(x, z, theta0, g, thermo);
    real T0 = isentropic_T(x, z, theta0, g, thermo);
    real dtheta = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_theta);
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var(p, T0 + dT, 0, 0, 0, 0);
  }
};

struct MoistLargeRisingBubble : LargeRisingBubble {

  static real YAKL_INLINE rhod_f(real x, real z,
                                 const ThermoPotential &thermo) {
    real p = isentropic_p(x, z, theta0, g, thermo);
    real T = isentropic_T(x, z, theta0, g, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE rhov_f(real x, real z,
                                 const ThermoPotential &thermo) {

    real pert = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_vapor);
    real Th = isentropic_T(x, z, theta0, g, thermo);
    real svp = saturation_vapor_pressure(Th);
    real pv = svp * pert;
    return pv / (thermo.cst.Rv * Th);
  }

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    real rhod = rhod_f(x, z, thermo);
    real rhov = rhov_f(x, z, thermo);
    return rhod + rhov;
  }

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {

    real p = isentropic_p(x, z, theta0, g, thermo);
    real T0 = isentropic_T(x, z, theta0, g, thermo);
    real dtheta = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_theta);
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var(p, T0 + dT, 1, 0, 0, 0);
  }

  static real YAKL_INLINE entropicdensity_f(real x, real z,
                                            const ThermoPotential &thermo) {
    return entropicvar_f(x, z, thermo) * rho_f(x, z, thermo);
  }
};

void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance) {
  if (name == "doublevortex") {
    testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  } else if (name == "risingbubble") {
    if (acoustic_balance) {
      testcase = std::make_unique<EulerTestCase<RisingBubble<true>>>();
    } else {
      testcase = std::make_unique<EulerTestCase<RisingBubble<false>>>();
    }
  } else if (name == "moistrisingbubble") {
    testcase = std::make_unique<MoistEulerTestCase<MoistRisingBubble>>();
  } else if (name == "largerisingbubble") {
    testcase = std::make_unique<EulerTestCase<LargeRisingBubble>>();
  } else if (name == "moistlargerisingbubble") {
    testcase = std::make_unique<MoistEulerTestCase<MoistLargeRisingBubble>>();
  } else {
    throw std::runtime_error("unknown test case");
  }
}
