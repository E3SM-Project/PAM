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

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "pocketfft_hdronly.h"

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
  void initialize(const Topology &ptopo, const Topology &dtopo,
                  const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom) override {
    // concentration 0-forms for dens
    name = "densl";
    topology = ptopo;
    dofs_arr = {0, 0, ndensity}; // densldiag = straight (0,0)-form
    Diagnostic::initialize(ptopo, dtopo, pgeom, dgeom);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

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
  void initialize(const Topology &ptopo, const Topology &dtopo,
                  const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom) override {
    name = "QXZl";
    topology = dtopo;
    dofs_arr = {0, 0, 1}; // // Qldiag = twisted (0,0)-form
    Diagnostic::initialize(ptopo, dtopo, pgeom, dgeom);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

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
                  Topology &primal_topo, Topology &dual_topo,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nconstant> &const_exchange) {

    ExtrudedTendencies::initialize(params, primal_topo, dual_topo, primal_geom,
                                   dual_geom, aux_exchange, const_exchange);
    varset.initialize(coupler, params, thermo, this->primal_topology,
                      this->dual_topology, this->primal_geometry,
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
                         FieldSet<nprognostic> &x) override {}

  void compute_dens0(real5d dens0var, const real5d densvar) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "Compute Dens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Iext<ndensity, diff_ord, vert_diff_ord>(
              dens0var, densvar, this->primal_geometry, this->dual_geometry,
              pis, pjs, pks, i, j, k, n);
        });
  }

  void compute_U(real5d Uvar, const real5d Vvar) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute Uvar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Hext<1, diff_ord>(Uvar, Vvar, this->primal_geometry,
                                    this->dual_geometry, dis, djs, dks, i, j, k,
                                    n);
        });
  }

  void compute_UW(real5d UWvar, const real5d Wvar) {
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute UWVar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Hv<1, vert_diff_ord>(UWvar, Wvar, this->primal_geometry,
                                       this->dual_geometry, dis, djs, dks, i, j,
                                       k + 1, n);
        });
  }

  void compute_q0f0(real5d qxz0var, real5d fxz0var, const real5d Vvar,
                    const real5d Wvar, const real5d densvar,
                    const real5d coriolisxzvar) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

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

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void compute_F_FW_and_K(real fac, real5d Fvar, real5d FWvar, real5d Kvar,
                          const real5d Vvar, const real5d Uvar,
                          const real5d Wvar, const real5d UWvar,
                          const real5d dens0var) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    // THIS WILL NEED SOME SLIGHT MODIFICATIONS FOR CASE OF NON-ZERO UWVAR_B IE
    // BOUNDARY FLUXES BUT FOR NOW IT IS FINE SINCE UWVAR=0 on BND AND THEREFORE
    // K COMPUTATIONS IGNORE IT
    YAKL_SCOPE(Hk, ::Hk);
    parallel_for(
        "Compute Fvar, Kvar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F<addmode>(Fvar, Uvar, dens0var, dis, djs, dks, i, j, k, n,
                                fac);
          Hk.compute_K(Kvar, Vvar, Uvar, Wvar, UWvar, dis, djs, dks, i, j, k,
                       n);
        });
    parallel_for(
        "Compute FWvar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_Fw<addmode>(FWvar, UWvar, dens0var, dis, djs, dks, i, j,
                                 k + 1, n, fac);
        });
  }

  void compute_F_FW_and_he(real5d Fvar, real5d FWvar, real5d HEvar,
                           real5d HEWvar, const real5d Vvar, const real5d Uvar,
                           const real5d Wvar, const real5d UWvar,
                           const real5d dens0var) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    // THIS WILL NEED SOME SLIGHT MODIFICATIONS FOR CASE OF NON-ZERO UWVAR_B IE
    // BOUNDARY FLUXES BUT FOR NOW IT IS FINE SINCE UWVAR=0 on BND AND THEREFORE
    // K COMPUTATIONS IGNORE IT
    YAKL_SCOPE(Hk, ::Hk);
    parallel_for(
        "Compute Fvar he",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F_and_he(Fvar, HEvar, Uvar, dens0var, dis, djs, dks, i, j,
                              k, n);
        });
    parallel_for(
        "Compute FWvar he",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_Fw_and_he(FWvar, HEWvar, UWvar, dens0var, dis, djs, dks, i,
                               j, k + 1, n);
        });
  }

  void compute_FT_and_FTW(real5d FTvar, real5d FTWvar, const real5d Fvar,
                          const real5d FWvar) {

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
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void compute_B(real fac, real5d Bvar, const real5d Kvar, const real5d densvar,
                 const real5d HSvar) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(Hk, ::Hk);
    YAKL_SCOPE(Hs, ::Hs);
    parallel_for(
        "Compute Bvar",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hs.compute_dHsdx<addmode>(Bvar, densvar, HSvar, pis, pjs, pks, i, j,
                                    k, n, fac);
          Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j, k,
                                            n, fac);
        });
  }

  void compute_edge_reconstructions(
      real5d densedgereconvar, real5d densvertedgereconvar,
      real5d qxzedgereconvar, real5d qxzvertedgereconvar,
      real5d coriolisxzedgereconvar, real5d coriolisxzvertedgereconvar,
      const real5d dens0var, const real5d qxz0var, const real5d fxz0var) {

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

  // LINEAR
  // void YAKL_INLINE compute_functional_derivatives(
  //    ADD_MODE addmode, real fac, real dt, FieldSet<nconstant> &const_vars,
  //    FieldSet<nprognostic> &x, FieldSet<nauxiliary> &auxiliary_vars) override
  //    {

  //  int pis = primal_topology.is;
  //  int pjs = primal_topology.js;
  //  int pks = primal_topology.ks;

  //  int dis = dual_topology.is;
  //  int djs = dual_topology.js;
  //  int dks = dual_topology.ks;

  //  auto densvar = x.fields_arr[DENSVAR].data;
  //  auto vvar = x.fields_arr[VVAR].data;
  //  auto wvar = x.fields_arr[WVAR].data;
  //
  //
  //  auto refdensvar = const_vars.fields_arr[REFDENSVAR].data;
  //  auto refnsq0var = const_vars.fields_arr[REFNSQ0VAR].data;
  //  auto refdens0var = const_vars.fields_arr[REFDENS0VAR].data;

  //  auto dens0var = auxiliary_vars.fields_arr[DENS0VAR].data;
  //  auto uvar = auxiliary_vars.fields_arr[UVAR].data;
  //  auto uwvar = auxiliary_vars.fields_arr[UWVAR].data;
  //  auto bvar = auxiliary_vars.fields_arr[BVAR].data;
  //
  //  auto fvar = auxiliary_vars.fields_arr[FVAR].data;
  //  auto fwvar = auxiliary_vars.fields_arr[FWVAR].data;
  //  auto hevar = auxiliary_vars.fields_arr[HEVAR].data;
  //  auto hewvar = auxiliary_vars.fields_arr[HEWVAR].data;

  //  auto densedgereconvar = auxiliary_vars.fields_arr[DENSEDGERECONVAR].data;
  //  auto densvertedgereconvar =
  //  auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data; auto densreconvar =
  //  auxiliary_vars.fields_arr[DENSRECONVAR].data; auto densvertreconvar =
  //  auxiliary_vars.fields_arr[DENSVERTRECONVAR].data;

  //  parallel_for(
  //      "Compute Dens0var",
  //      SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
  //                      primal_topology.n_cells_x, primal_topology.nens),
  //      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
  //        compute_Iext<ndensity, diff_ord, vert_diff_ord>(
  //            dens0var, densvar, this->primal_geometry, this->dual_geometry,
  //            pis, pjs, pks, i, j, k, n);
  //      });

  //  parallel_for(
  //      "Compute Uvar",
  //      SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
  //        compute_Hext<1, diff_ord>(uvar, vvar, this->primal_geometry,
  //                                  this->dual_geometry, dis, djs, dks, i, j,
  //                                  k, n);
  //      });

  //  parallel_for(
  //      "Compute UWVar",
  //      SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
  //        compute_Hv<1, vert_diff_ord>(uwvar, wvar, this->primal_geometry,
  //                                     this->dual_geometry, dis, djs, dks, i,
  //                                     j, k + 1, n);
  //      });

  //  auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
  //  this->aux_exchange->exchanges_arr[UVAR].exchange_field(
  //      auxiliary_vars.fields_arr[UVAR]);
  //  this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
  //      auxiliary_vars.fields_arr[UWVAR]);
  //  this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
  //      auxiliary_vars.fields_arr[DENS0VAR]);

  //  parallel_for(
  //      "Compute F",
  //      SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        Hk.compute_F_and_he(fvar, hevar, uvar, refdens0var, dis, djs, dks,
  //        i, j, k, n);
  //  });

  //  parallel_for(
  //      "Compute FWvar",
  //      SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        Hk.compute_Fw_and_he(fwvar, hewvar, uwvar, refdens0var, dis, djs,
  //        dks, i, j,
  //                      k + 1, n);
  //      });
  //
  //  auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
  //  this->aux_exchange->exchanges_arr[FVAR].exchange_field(
  //      auxiliary_vars.fields_arr[FVAR]);
  //  this->aux_exchange->exchanges_arr[FWVAR].exchange_field(
  //      auxiliary_vars.fields_arr[FWVAR]);
  //  this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
  //      auxiliary_vars.fields_arr[HEVAR]);
  //  this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(
  //      auxiliary_vars.fields_arr[HEWVAR]);

  //  parallel_for(
  //      "Compute Bvar",
  //      SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
  //                      primal_topology.n_cells_x, primal_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {

  //        real Rd = thermo.cst.Rd;
  //        real pr = thermo.cst.pr;
  //        real gamma_d = thermo.cst.gamma_d;
  //        real Cpd = thermo.cst.Cpd;
  //        real Cvd = thermo.cst.Cvd;
  //        real grav = 9.80616_fp;
  //        real grav2 = grav * grav;
  //
  //        real rho_ref = refdens0var(0, pks + k, pjs + j, pis + i, n);
  //        real Tht_ref = refdens0var(1, pks + k, pjs + j, pis + i, n);
  //        real tht_ref = Tht_ref / rho_ref;

  //        //real rho = rho_ref + dens0var(0, pks + k, pjs + j, pis + i, n);
  //        //real Tht = Tht_ref + dens0var(1, pks + k, pjs + j, pis + i, n);
  //        real rho = dens0var(0, pks + k, pjs + j, pis + i, n);
  //        real Tht = dens0var(1, pks + k, pjs + j, pis + i, n);

  //        real tht = Tht / rho_ref - (rho - rho_ref) / rho_ref * tht_ref;
  //        real rho_ref2 = rho_ref * rho_ref;
  //        real p_ref = pr * pow(Rd * Tht_ref / pr, gamma_d);
  //        real dpdtht_ref = gamma_d * p_ref / tht_ref;
  //        real dpdtht_ref2 = dpdtht_ref * dpdtht_ref;
  //        real Nref2 = refnsq0var(0, pks + k, pjs + j, pis + i, n);
  //        real cref2 = gamma_d * p_ref / rho_ref;

  //        real drho = rho - rho_ref;
  //        real dtht = tht - tht_ref;

  //        real b0_drho = (cref2 * rho_ref - dpdtht_ref * tht_ref) / rho_ref2;
  //        real b0_dtht = dpdtht_ref / rho_ref - dpdtht_ref2 * tht_ref / (cref2
  //        * rho_ref2) -
  //                      dpdtht_ref2 * grav2 * tht_ref / (Nref2 * cref2 * cref2
  //                      * rho_ref2);

  //        real b0 = drho * b0_drho + dtht * b0_dtht;

  //        real b1_drho = dpdtht_ref / rho_ref2;
  //        real b1_dtht = dpdtht_ref2 * (Nref2 * cref2 + grav2) / (Nref2 *
  //        cref2 * cref2 * rho_ref2);

  //        real b1 = drho * b1_drho + dtht * b1_dtht;

  //        bvar(0, pks + k, pjs + j, pis + i, n) = b0;
  //        bvar(1, pks + k, pjs + j, pis + i, n) = b1;
  //      });
  //
  //  this->aux_exchange->exchanges_arr[BVAR].exchange_field(
  //      auxiliary_vars.fields_arr[BVAR]);
  //}

  //// LINEAR
  // void YAKL_INLINE apply_symplectic(real dt, FieldSet<nconstant> &const_vars,
  //                  FieldSet<nprognostic> &x,
  //                  FieldSet<nauxiliary> &auxiliary_vars,
  //                  FieldSet<nprognostic> &xtend) override {

  //  int pis = primal_topology.is;
  //  int pjs = primal_topology.js;
  //  int pks = primal_topology.ks;

  //  int dis = dual_topology.is;
  //  int djs = dual_topology.js;
  //  int dks = dual_topology.ks;

  //  auto densvar = x.fields_arr[DENSVAR].data;
  //  auto vvar = x.fields_arr[VVAR].data;
  //  auto wvar = x.fields_arr[WVAR].data;
  //
  //  auto denstendvar = xtend.fields_arr[DENSVAR].data;
  //  auto vtendvar = xtend.fields_arr[VVAR].data;
  //  auto wtendvar = xtend.fields_arr[WVAR].data;
  //
  //  auto refdensvar = const_vars.fields_arr[REFDENSVAR].data;
  //  auto refnsq0var = const_vars.fields_arr[REFNSQ0VAR].data;
  //  auto refdens0var = const_vars.fields_arr[REFDENS0VAR].data;

  //  auto dens0var = auxiliary_vars.fields_arr[DENS0VAR].data;
  //  auto uvar = auxiliary_vars.fields_arr[UVAR].data;
  //  auto uwvar = auxiliary_vars.fields_arr[UWVAR].data;
  //  auto bvar = auxiliary_vars.fields_arr[BVAR].data;
  //
  //  auto fvar = auxiliary_vars.fields_arr[FVAR].data;
  //  auto fwvar = auxiliary_vars.fields_arr[FWVAR].data;
  //  auto fvar2 = auxiliary_vars.fields_arr[FVAR2].data;
  //  auto fwvar2 = auxiliary_vars.fields_arr[FWVAR2].data;
  //  auto hevar = auxiliary_vars.fields_arr[HEVAR].data;
  //  auto hewvar = auxiliary_vars.fields_arr[HEWVAR].data;

  //  auto densedgereconvar = auxiliary_vars.fields_arr[DENSEDGERECONVAR].data;
  //  auto densvertedgereconvar =
  //  auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data; auto densreconvar =
  //  auxiliary_vars.fields_arr[DENSRECONVAR].data; auto densvertreconvar =
  //  auxiliary_vars.fields_arr[DENSVERTRECONVAR].data;

  //  parallel_for(
  //      "Compute Dens0var",
  //      SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
  //                      primal_topology.n_cells_x, primal_topology.nens),
  //      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
  //        compute_Iext<ndensity, diff_ord, vert_diff_ord>(
  //            dens0var, densvar, this->primal_geometry, this->dual_geometry,
  //            pis, pjs, pks, i, j, k, n);
  //      });

  //  parallel_for(
  //      "Compute Uvar",
  //      SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
  //        compute_Hext<1, diff_ord>(uvar, vvar, this->primal_geometry,
  //                                  this->dual_geometry, dis, djs, dks, i, j,
  //                                  k, n);
  //      });

  //  parallel_for(
  //      "Compute UWVar",
  //      SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
  //        compute_Hv<1, vert_diff_ord>(uwvar, wvar, this->primal_geometry,
  //                                     this->dual_geometry, dis, djs, dks, i,
  //                                     j, k + 1, n);
  //      });

  //  auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
  //  this->aux_exchange->exchanges_arr[UVAR].exchange_field(
  //      auxiliary_vars.fields_arr[UVAR]);
  //  this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
  //      auxiliary_vars.fields_arr[UWVAR]);
  //  this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
  //      auxiliary_vars.fields_arr[DENS0VAR]);

  //  parallel_for(
  //      "Compute F",
  //      SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        Hk.compute_F_and_he(fvar2, hevar, uvar, refdens0var, dis, djs, dks,
  //        i, j, k, n);
  //  });

  //  parallel_for(
  //      "Compute FWvar",
  //      SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        Hk.compute_Fw_and_he(fwvar2, hewvar, uwvar, refdens0var, dis, djs,
  //        dks, i, j,
  //                      k + 1, n);
  //      });
  //
  //  auxiliary_vars.fields_arr[FWVAR2].set_bnd(0.0);
  //  this->aux_exchange->exchanges_arr[FVAR2].exchange_field(
  //      auxiliary_vars.fields_arr[FVAR2]);
  //  this->aux_exchange->exchanges_arr[FWVAR2].exchange_field(
  //      auxiliary_vars.fields_arr[FWVAR2]);
  //  this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
  //      auxiliary_vars.fields_arr[HEVAR]);
  //  this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(
  //      auxiliary_vars.fields_arr[HEWVAR]);

  //  parallel_for(
  //      "ComputeDensEdgeRecon",
  //      SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
  //        compute_twisted_edge_recon<ndensity, dual_reconstruction_type,
  //                                   dual_reconstruction_order>(
  //            densedgereconvar, refdens0var, dis, djs, dks, i, j, k, n,
  //            dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
  //        compute_twisted_vert_edge_recon<ndensity,
  //                                        dual_vert_reconstruction_type,
  //                                        dual_vert_reconstruction_order>(
  //            densvertedgereconvar, refdens0var, dis, djs, dks, i, j, k, n,
  //            dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl,
  //            dual_vert_wenoSigma);
  //      });

  //  this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(
  //      auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
  //  this->aux_exchange->exchanges_arr[DENSVERTEDGERECONVAR].exchange_field(
  //      auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR]);

  //  parallel_for(
  //      "ComputeDens recon",
  //      SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        compute_twisted_recon<ndensity, dual_reconstruction_type>(
  //            densreconvar, densedgereconvar, uvar, dis, djs, dks, i, j, k,
  //            n);

  //        // scale twisted recons
  //        for (int d = 0; d < ndims; d++) {
  //          for (int l = 0; l < ndensity; l++) {
  //            densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) =
  //                densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n)
  //                / hevar(d, k + dks, j + djs, i + dis, n);
  //          }
  //        }
  //      });

  //  parallel_for(
  //      "ComputeDensVertRECON",
  //      SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        compute_twisted_vert_recon<ndensity, dual_vert_reconstruction_type>(
  //            densvertreconvar, densvertedgereconvar, uwvar, dis, djs, dks, i,
  //            j, k + 1, n);

  //        // scale twisted recons
  //        for (int l = 0; l < ndensity; l++) {
  //          densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) =
  //              densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) /
  //              hewvar(0, k + dks + 1, j + djs, i + dis, n);
  //        }
  //      });

  //  this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(
  //      auxiliary_vars.fields_arr[DENSRECONVAR]);
  //  this->aux_exchange->exchanges_arr[DENSVERTRECONVAR].exchange_field(
  //      auxiliary_vars.fields_arr[DENSVERTRECONVAR]);

  //  parallel_for(
  //      "Compute Wtend",
  //      SimpleBounds<4>(primal_topology.nl - 2, primal_topology.n_cells_y,
  //                      primal_topology.n_cells_x, primal_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        compute_wDv<ndensity>(wtendvar, densvertreconvar,
  //                                  bvar, pis, pjs, pks, i, j, k + 1, n);
  //      });

  //  parallel_for(
  //      "Compute Wtend Bnd",
  //      SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
  //                      primal_topology.nens),
  //      YAKL_CLASS_LAMBDA(int j, int i, int n) {
  //        compute_wDv<ndensity>(wtendvar, densvertreconvar,
  //                                  bvar, pis, pjs, pks, i, j, 0, n);
  //        compute_wDv<ndensity>(wtendvar, densvertreconvar,
  //                                  bvar, pis, pjs, pks, i, j,
  //                                  primal_topology.nl - 1, n);
  //      });
  //
  //  parallel_for(
  //      "Compute Vtend",
  //      SimpleBounds<4>(primal_topology.ni - 2, primal_topology.n_cells_y,
  //                      primal_topology.n_cells_x, primal_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        compute_wD1<ndensity>(vtendvar, densreconvar, bvar, pis,
  //                                  pjs, pks, i, j, k + 1, n);
  //  });

  //  parallel_for(
  //      "Compute Vtend Bnd",
  //      SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
  //                      primal_topology.nens),
  //      YAKL_CLASS_LAMBDA(int j, int i, int n) {
  //        compute_wD1<ndensity>(vtendvar, densreconvar, bvar, pis,
  //                                  pjs, pks, i, j, 0, n);
  //        compute_wD1<ndensity>(vtendvar, densreconvar, bvar, pis,
  //                                  pjs, pks, i, j, primal_topology.ni - 1,
  //                                  n);
  //      });

  //  parallel_for(
  //      "Compute Dens Tend",
  //      SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
  //                      dual_topology.n_cells_x, dual_topology.nens),
  //      YAKL_LAMBDA(int k, int j, int i, int n) {
  //        compute_wDbar2<ndensity>(denstendvar, densreconvar, fvar,
  //                                     dis, djs, dks, i, j, k, n);
  //        compute_wDvbar<ndensity, ADD_MODE::ADD>(
  //            denstendvar, densvertreconvar, fwvar, dis, djs, dks,
  //            i, j, k, n);
  //      });
  //}

  // NONLINEAR
  void YAKL_INLINE compute_functional_derivatives(
      ADD_MODE addmode, real fac, real dt, FieldSet<nconstant> &const_vars,
      FieldSet<nprognostic> &x, FieldSet<nauxiliary> &auxiliary_vars) override {

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U(auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[VVAR].data);
    compute_UW(auxiliary_vars.fields_arr[UWVAR].data, x.fields_arr[WVAR].data);
    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
        auxiliary_vars.fields_arr[UWVAR]);
    this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
        auxiliary_vars.fields_arr[DENS0VAR]);

    if (addmode == ADD_MODE::ADD) {
      compute_F_FW_and_K<ADD_MODE::ADD>(
          fac, auxiliary_vars.fields_arr[FVAR].data,
          auxiliary_vars.fields_arr[FWVAR].data,
          auxiliary_vars.fields_arr[KVAR].data, x.fields_arr[VVAR].data,
          auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[WVAR].data,
          auxiliary_vars.fields_arr[UWVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data);
    } else if (addmode == ADD_MODE::REPLACE) {
      compute_F_FW_and_K<ADD_MODE::REPLACE>(
          fac, auxiliary_vars.fields_arr[FVAR].data,
          auxiliary_vars.fields_arr[FWVAR].data,
          auxiliary_vars.fields_arr[KVAR].data, x.fields_arr[VVAR].data,
          auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[WVAR].data,
          auxiliary_vars.fields_arr[UWVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data);
    }

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[FVAR].exchange_field(
        auxiliary_vars.fields_arr[FVAR]);
    this->aux_exchange->exchanges_arr[FWVAR].exchange_field(
        auxiliary_vars.fields_arr[FWVAR]);
    this->aux_exchange->exchanges_arr[KVAR].exchange_field(
        auxiliary_vars.fields_arr[KVAR]);

    if (addmode == ADD_MODE::ADD) {
      compute_B<ADD_MODE::ADD>(fac, auxiliary_vars.fields_arr[BVAR].data,
                               auxiliary_vars.fields_arr[KVAR].data,
                               x.fields_arr[DENSVAR].data,
                               const_vars.fields_arr[HSVAR].data);
    } else if (addmode == ADD_MODE::REPLACE) {
      compute_B<ADD_MODE::REPLACE>(fac, auxiliary_vars.fields_arr[BVAR].data,
                                   auxiliary_vars.fields_arr[KVAR].data,
                                   x.fields_arr[DENSVAR].data,
                                   const_vars.fields_arr[HSVAR].data);
    }
    this->aux_exchange->exchanges_arr[BVAR].exchange_field(
        auxiliary_vars.fields_arr[BVAR]);
  }

  // NONLINEAR
  void YAKL_INLINE apply_symplectic(real dt, FieldSet<nconstant> &const_vars,
                                    FieldSet<nprognostic> &x,
                                    FieldSet<nauxiliary> &auxiliary_vars,
                                    FieldSet<nprognostic> &xtend) override {

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U(auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[VVAR].data);
    compute_UW(auxiliary_vars.fields_arr[UWVAR].data, x.fields_arr[WVAR].data);
    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);

    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
        auxiliary_vars.fields_arr[UWVAR]);
    this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
        auxiliary_vars.fields_arr[DENS0VAR]);

    compute_F_FW_and_he(
        auxiliary_vars.fields_arr[FVAR2].data,
        auxiliary_vars.fields_arr[FWVAR2].data,
        auxiliary_vars.fields_arr[HEVAR].data,
        auxiliary_vars.fields_arr[HEWVAR].data, x.fields_arr[VVAR].data,
        auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[WVAR].data,
        auxiliary_vars.fields_arr[UWVAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data);

    auxiliary_vars.fields_arr[FWVAR2].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[FVAR2].exchange_field(
        auxiliary_vars.fields_arr[FVAR2]);
    this->aux_exchange->exchanges_arr[FWVAR2].exchange_field(
        auxiliary_vars.fields_arr[FWVAR2]);

    this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
        auxiliary_vars.fields_arr[HEVAR]);
    this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(
        auxiliary_vars.fields_arr[HEWVAR]);

    compute_q0f0(auxiliary_vars.fields_arr[QXZ0VAR].data,
                 auxiliary_vars.fields_arr[FXZ0VAR].data,
                 x.fields_arr[VVAR].data, x.fields_arr[WVAR].data,
                 x.fields_arr[DENSVAR].data,
                 const_vars.fields_arr[CORIOLISXZVAR].data);

    auxiliary_vars.fields_arr[QXZ0VAR].set_bnd(0.0);
    auxiliary_vars.fields_arr[FXZ0VAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[QXZ0VAR].exchange_field(
        auxiliary_vars.fields_arr[QXZ0VAR]);
    this->aux_exchange->exchanges_arr[FXZ0VAR].exchange_field(
        auxiliary_vars.fields_arr[FXZ0VAR]);

    compute_FT_and_FTW(auxiliary_vars.fields_arr[FTVAR].data,
                       auxiliary_vars.fields_arr[FTWVAR].data,
                       auxiliary_vars.fields_arr[FVAR2].data,
                       auxiliary_vars.fields_arr[FWVAR2].data);

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

  void compute_linrhs(real dt, FieldSet<nconstant> &const_vars,
                      FieldSet<nprognostic> &x,
                      FieldSet<nauxiliary> &auxiliary_vars,
                      FieldSet<nprognostic> &xtend) override {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    auto densvar = x.fields_arr[DENSVAR].data;
    auto vvar = x.fields_arr[VVAR].data;
    auto wvar = x.fields_arr[WVAR].data;

    auto denstendvar = xtend.fields_arr[DENSVAR].data;
    auto vtendvar = xtend.fields_arr[VVAR].data;
    auto wtendvar = xtend.fields_arr[WVAR].data;

    auto refdensvar = const_vars.fields_arr[REFDENSVAR].data;
    auto refnsq0var = const_vars.fields_arr[REFNSQ0VAR].data;
    auto refdens0var = const_vars.fields_arr[REFDENS0VAR].data;

    auto dens0var = auxiliary_vars.fields_arr[DENS0VAR].data;
    auto uvar = auxiliary_vars.fields_arr[UVAR].data;
    auto uwvar = auxiliary_vars.fields_arr[UWVAR].data;
    auto bvar = auxiliary_vars.fields_arr[BVAR].data;

    auto fvar = auxiliary_vars.fields_arr[FVAR].data;
    auto fwvar = auxiliary_vars.fields_arr[FWVAR].data;
    auto hevar = auxiliary_vars.fields_arr[HEVAR].data;
    auto hewvar = auxiliary_vars.fields_arr[HEWVAR].data;

    auto densedgereconvar = auxiliary_vars.fields_arr[DENSEDGERECONVAR].data;
    auto densvertedgereconvar =
        auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data;
    auto densreconvar = auxiliary_vars.fields_arr[DENSRECONVAR].data;
    auto densvertreconvar = auxiliary_vars.fields_arr[DENSVERTRECONVAR].data;
    
    auto bdens0var = const_vars.fields_arr[BDENS0VAR].data;
    
    auto refhevar = const_vars.fields_arr[REFHEVAR].data;
    auto refhewvar = const_vars.fields_arr[REFHEWVAR].data;
    auto refdensedgereconvar = const_vars.fields_arr[REFDENSEDGERECONVAR].data;
    auto refdensvertedgereconvar =
        const_vars.fields_arr[REFDENSVERTEDGERECONVAR].data;
    auto refdensreconvar = const_vars.fields_arr[REFDENSRECONVAR].data;
    auto refdensvertreconvar = const_vars.fields_arr[REFDENSVERTRECONVAR].data;

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
          compute_Hext<1, diff_ord>(uvar, vvar, this->primal_geometry,
                                    this->dual_geometry, dis, djs, dks, i, j, k,
                                    n);
        });

    parallel_for(
        "Compute UWVar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Hv<1, vert_diff_ord>(uwvar, wvar, this->primal_geometry,
                                       this->dual_geometry, dis, djs, dks, i, j,
                                       k + 1, n);
        });

    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
        auxiliary_vars.fields_arr[UWVAR]);
    this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
        auxiliary_vars.fields_arr[DENS0VAR]);

    parallel_for(
        "Compute F",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F_and_he(fvar, hevar, uvar, refdens0var, dis, djs, dks, i,
                              j, k, n);
        });

    parallel_for(
        "Compute FWvar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_Fw_and_he(fwvar, hewvar, uwvar, refdens0var, dis, djs, dks,
                               i, j, k + 1, n);
        });

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[FVAR].exchange_field(
        auxiliary_vars.fields_arr[FVAR]);
    this->aux_exchange->exchanges_arr[FWVAR].exchange_field(
        auxiliary_vars.fields_arr[FWVAR]);
    this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
        auxiliary_vars.fields_arr[HEVAR]);
    this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(
        auxiliary_vars.fields_arr[HEWVAR]);

    parallel_for(
        "Compute Bvar",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {

          real rho = dens0var(0, pks + k, pjs + j, pis + i, n);
          real Tht = dens0var(1, pks + k, pjs + j, pis + i, n);
          
          real b0_rho = bdens0var(0, pks + k, pjs + j, pis + i, n);
          real b0_Tht = bdens0var(1, pks + k, pjs + j, pis + i, n);
          real b1_rho = bdens0var(2, pks + k, pjs + j, pis + i, n);
          real b1_Tht = bdens0var(3, pks + k, pjs + j, pis + i, n);

          real b0 = rho * b0_rho + Tht * b0_Tht;
          real b1 = rho * b1_rho + Tht * b1_Tht;

          bvar(0, pks + k, pjs + j, pis + i, n) = b0;
          bvar(1, pks + k, pjs + j, pis + i, n) = b1;
        });

    this->aux_exchange->exchanges_arr[BVAR].exchange_field(
        auxiliary_vars.fields_arr[BVAR]);

    //parallel_for(
    //    "ComputeDensEdgeRecon",
    //    SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
    //      compute_twisted_edge_recon<ndensity, dual_reconstruction_type,
    //                                 dual_reconstruction_order>(
    //          densedgereconvar, refdens0var, dis, djs, dks, i, j, k, n,
    //          dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
    //      compute_twisted_vert_edge_recon<ndensity,
    //                                      dual_vert_reconstruction_type,
    //                                      dual_vert_reconstruction_order>(
    //          densvertedgereconvar, refdens0var, dis, djs, dks, i, j, k, n,
    //          dual_vert_wenoRecon, dual_vert_to_gll, dual_vert_wenoIdl,
    //          dual_vert_wenoSigma);
    //    });

    //this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(
    //    auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
    //this->aux_exchange->exchanges_arr[DENSVERTEDGERECONVAR].exchange_field(
    //    auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR]);

    //parallel_for(
    //    "ComputeDens recon",
    //    SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      compute_twisted_recon<ndensity, dual_reconstruction_type>(
    //          densreconvar, densedgereconvar, uvar, dis, djs, dks, i, j, k, n);

    //      // scale twisted recons
    //      for (int d = 0; d < ndims; d++) {
    //        for (int l = 0; l < ndensity; l++) {
    //          densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) =
    //              densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) /
    //              hevar(d, k + dks, j + djs, i + dis, n);
    //        }
    //      }
    //    });

    //parallel_for(
    //    "ComputeDensVertRECON",
    //    SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      compute_twisted_vert_recon<ndensity, dual_vert_reconstruction_type>(
    //          densvertreconvar, densvertedgereconvar, uwvar, dis, djs, dks, i,
    //          j, k + 1, n);

    //      // scale twisted recons
    //      for (int l = 0; l < ndensity; l++) {
    //        densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) =
    //            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) /
    //            hewvar(0, k + dks + 1, j + djs, i + dis, n);
    //      }
    //    });

    //this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(
    //    auxiliary_vars.fields_arr[DENSRECONVAR]);
    //this->aux_exchange->exchanges_arr[DENSVERTRECONVAR].exchange_field(
    //    auxiliary_vars.fields_arr[DENSVERTRECONVAR]);

    parallel_for(
        "Compute Wtend",
        SimpleBounds<4>(primal_topology.nl - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDv<ndensity>(wtendvar, refdensvertreconvar, bvar, pis, pjs, pks,
                                i, j, k + 1, n);
        });

    parallel_for(
        "Compute Wtend Bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          compute_wDv<ndensity>(wtendvar, refdensvertreconvar, bvar, pis, pjs, pks,
                                i, j, 0, n);
          compute_wDv<ndensity>(wtendvar, refdensvertreconvar, bvar, pis, pjs, pks,
                                i, j, primal_topology.nl - 1, n);
        });

    parallel_for(
        "Compute Vtend",
        SimpleBounds<4>(primal_topology.ni - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD1<ndensity>(vtendvar, refdensreconvar, bvar, pis, pjs, pks, i,
                                j, k + 1, n);
        });

    parallel_for(
        "Compute Vtend Bnd",
        SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                        primal_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          compute_wD1<ndensity>(vtendvar, refdensreconvar, bvar, pis, pjs, pks, i,
                                j, 0, n);
          compute_wD1<ndensity>(vtendvar, refdensreconvar, bvar, pis, pjs, pks, i,
                                j, primal_topology.ni - 1, n);
        });

    parallel_for(
        "Compute Dens Tend",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDbar2<ndensity>(denstendvar, refdensreconvar, fvar, dis, djs,
                                   dks, i, j, k, n);
          compute_wDvbar<ndensity, ADD_MODE::ADD>(
              denstendvar, refdensvertreconvar, fwvar, dis, djs, dks, i, j, k, n);
        });
  }
};

int offset(const Field &field) {
  auto topo = field.topology;

  auto ndofs = field.total_dofs;
  auto nz = field._nz;

  auto halosize_x = topo.halosize_x;
  auto halosize_y = topo.halosize_y;

  auto n_cells_x = topo.n_cells_x;
  auto n_cells_y = topo.n_cells_y;

  auto nx = n_cells_x + 2 * halosize_x;
  auto ny = n_cells_y + 2 * halosize_y;

  auto nens = topo.nens;

  int ntot = ndofs * nz * n_cells_y * n_cells_x * nens;
  return nens * (halosize_x + nx * (halosize_y + ny * (topo.mirror_halo)));
}

int realsize(const Field &field) {
  auto topo = field.topology;

  auto ndofs = field.total_dofs;
  auto nz = field._nz;

  auto n_cells_x = topo.n_cells_x;
  auto n_cells_y = topo.n_cells_y;

  auto nens = topo.nens;

  return ndofs * nz * n_cells_y * n_cells_x * nens;
}

int realsize(const FieldSet<nprognostic> &vars) {
  int accum = 0;
  for (auto f : vars.fields_arr) {
    accum += realsize(f);
  }

  return accum;
}

real getindex(const Field &field, int ix) {
  auto topo = field.topology;

  auto ndofs = field.total_dofs;
  auto nz = field._nz;
  auto n_cells_x = topo.n_cells_x;
  auto n_cells_y = topo.n_cells_y;
  auto nens = topo.nens;

  int n = ix % nens;
  int i = (ix / nens) % n_cells_x;
  int j = (ix / nens / n_cells_x) % n_cells_y;
  int k = (ix / nens / n_cells_x / n_cells_y) % nz;
  int d = ix / nens / n_cells_x / n_cells_y / nz;

  return field.data(d, topo.ks + k, topo.js + j, topo.is + i, n);
}

real getindex(const FieldSet<nprognostic> &vars, int ix) {
  int s = 0;
  while (ix / realsize(vars.fields_arr[s]) > 0) {
    ix -= realsize(vars.fields_arr[s]);
    s++;
  }
  auto field = vars.fields_arr[s];
  return getindex(field, ix);
}

void setindex(Field &field, int ix, real val) {
  auto topo = field.topology;

  auto ndofs = field.total_dofs;
  auto nz = field._nz;
  auto n_cells_x = topo.n_cells_x;
  auto n_cells_y = topo.n_cells_y;
  auto nens = topo.nens;

  int n = ix % nens;
  int i = (ix / nens) % n_cells_x;
  int j = (ix / nens / n_cells_x) % n_cells_y;
  int k = (ix / nens / n_cells_x / n_cells_y) % nz;
  int d = ix / nens / n_cells_x / n_cells_y / nz;

  field.data(d, topo.ks + k, topo.js + j, topo.is + i, n) = val;
}

void setindex(FieldSet<nprognostic> &vars, int ix, real val) {
  int s = 0;
  while (ix / realsize(vars.fields_arr[s]) > 0) {
    ix -= realsize(vars.fields_arr[s]);
    s++;
  }
  auto field = vars.fields_arr[s];
  setindex(field, ix, val);
}

void set_one(FieldSet<nprognostic> &vars, int ix) {
  for (auto f : vars.fields_arr) {
    f.set(0);
  }
  setindex(vars, ix, 1);
}

// *******   Linear system   ***********//
class ModelLinearSystem : public LinearSystem {
  Eigen::SparseLU<Eigen::SparseMatrix<real>> lu;
  Eigen::VectorXd eigen_rhs;
  Eigen::VectorXd eigen_sol;
  
  complex5d complex_v;
  complex5d complex_vrhs;
  complex5d complex_w;
  complex5d complex_wrhs;
  complex5d complex_wrhs2;
  complex5d complex_dens;
  complex5d complex_densrhs;
  
  complex5d complex_v_halo;
  complex5d complex_vcoeff_halo;
  complex5d complex_vrhs_halo;
  complex5d complex_w_halo;
  complex5d complex_uw_halo;
  complex5d complex_wrhs_halo;
  complex5d complex_dens_halo;
  complex5d complex_dens0_halo;
  complex5d complex_densrhs_halo;
  
  complex1d tri_l;
  complex1d tri_d;
  complex1d tri_u;
  complex1d tri_rhs;
  
  real5d vtend;
  real5d wtend;
  real5d wtend2;
  real5d denstend;

public:
  void initialize(ModelParameters &params, Tendencies *tend,
                  FieldSet<nprognostic> &x, FieldSet<nconstant> &const_vars,
                  FieldSet<nauxiliary> &auxiliary_vars,
                  ExchangeSet<nprognostic> &prog_exchange) override {
    LinearSystem::initialize(params, tend, x, const_vars, auxiliary_vars,
                             prog_exchange);

    auto N = realsize(x);

    this->eigen_rhs = Eigen::VectorXd(N);
    this->eigen_sol = Eigen::VectorXd(N);

    FieldSet<nprognostic> F;
    FieldSet<nprognostic> Q;

    F.initialize(x, "F");
    Q.initialize(x, "Q");

    //FieldSet<nprognostic> Z;
    //FieldSet<nprognostic> FZ;
    //Z.initialize(x, "Z");
    //FZ.initialize(x, "FZ");
    //tend->compute_linrhs(params.dtcrm, const_vars, Z, auxiliary_vars, FZ);

    yakl::timer_start("compute coefficients");
    using TT = Eigen::Triplet<real>;
    std::vector<Eigen::Triplet<real>> coefficients;
    for (int i = 0; i < N; ++i) {
      set_one(Q, i);
      prog_exchange.exchange_variable_set(Q);
      tend->compute_linrhs(params.dtcrm, const_vars, Q, auxiliary_vars, F);
      //F.waxpy(-1, FZ, F);
      for (int j = 0; j < N; ++j) {
        auto val = 0.5 * params.dtcrm * getindex(F, j);
        if (i == j) {
          val += 1;
        }
        // auto val = getindex(F, j);
        if (std::abs(val) > 1e-10) {
          coefficients.push_back(TT(j, i, val));
        }
      }
    }
    yakl::timer_stop("compute coefficients");

    Eigen::SparseMatrix<real> L(N, N);
    L.setFromTriplets(coefficients.begin(), coefficients.end());

    // checks
    // Eigen::VectorXd ex(N);
    // for (int i = 0; i < N; ++i) {
    //  ex(i) = getindex(x, i);
    //}
    // auto eLx = L * ex;
    //
    // prog_exchange.exchange_variable_set(x);
    // tend->compute_linrhs(params.dtcrm, const_vars, x, auxiliary_vars, F);
    // F.waxpy(-1, FZ, F);

    // for (int i = 0; i < N; ++i) {
    //   real Lx = getindex(F, i);

    //  std::cout << i << " " << eLx(i) - Lx << std::endl;
    //}
    // exit(1);

    yakl::timer_start("LU");
    // this->lu.analyzePattern(L);
    // this->lu.factorize(L);
    this->lu.compute(L);
    yakl::timer_stop("LU");


    // fourier initialize

    auto pni = primal_topology.ni;
    auto pnl = primal_topology.nl;
    auto dni = dual_topology.ni;
    auto dnl = dual_topology.nl;
    auto nx = primal_topology.n_cells_x;
    auto ny = primal_topology.n_cells_y;
    auto nens = primal_topology.nens;

    auto nxh = nx + 2 * primal_topology.halosize_x;
    auto nyh = nx + 2 * primal_topology.halosize_y;
    auto pnih = pni + 2 * primal_topology.mirror_halo;
    auto pnlh = pnl + 2 * primal_topology.mirror_halo;
    auto dnih = dni + 2 * primal_topology.mirror_halo;
    auto dnlh = dnl + 2 * primal_topology.mirror_halo;

    complex_v = complex5d("complex v", 1, pni, ny,
                          nx, nens);
    complex_vrhs = complex5d("complex vrhs", 1, pni, ny,
                          nx, nens);
    complex_w = complex5d("complex w", 1, pnl, ny,
                          nx, nens);
    complex_wrhs = complex5d("complex wrhs", 1, pnl, ny,
                          nx, nens);
    complex_wrhs2 = complex5d("complex wrhs", 1, pnl, ny,
                          nx, nens);
    complex_dens = complex5d("complex dens", ndensity, dnl, ny,
                          nx, nens);
    complex_densrhs = complex5d("complex densrhs", ndensity, dnl, ny,
                          nx, nens);

    complex_v_halo = complex5d("complex v halo", 1, pnih, nyh,
                          nxh, nens);
    complex_vcoeff_halo = complex5d("complex vcoeff halo", 3, pnih, nyh,
                          nxh, nens);
    complex_vrhs_halo = complex5d("complex vrhs halo", 1, pnih, nyh,
                          nxh, nens);
    complex_w_halo = complex5d("complex w halo", 1, pnlh, nyh,
                          nxh, nens);
    complex_uw_halo = complex5d("complex uw halo", 1, pnlh + 2, nyh,
                          nxh, nens);
    complex_wrhs_halo = complex5d("complex wrhs halo", 1, pnlh, nyh,
                          nxh, nens);
    complex_dens_halo = complex5d("complex dens halo", ndensity, dnlh, nyh,
                          nxh, nens);
    complex_dens0_halo = complex5d("complex dens0 halo", ndensity, dnlh, nyh,
                          nxh, nens);
    complex_densrhs_halo = complex5d("complex densrhs halo", ndensity, dnlh, nyh,
                          nxh, nens);
    
    vtend = real5d("complex v halo", 1, pnih, nyh,
                   nxh, nens);
    wtend = real5d("complex w halo", 1, pnlh, nyh,
                    nxh, nens);
    wtend2 = real5d("complex w halo", 1, pnlh, nyh,
                    nxh, nens);
    denstend = real5d("complex v halo", ndensity, dnlh, nyh,
                   nxh, nens);
    
    tri_d = complex1d("tri d", pnl);
    tri_rhs = complex1d("tri rhs", pnl);
    tri_l = complex1d("tri l", pnl-1);
    tri_u = complex1d("tri u", pnl-1);
  }

  virtual void YAKL_INLINE solve(real dt, FieldSet<nprognostic> &rhs,
                                 FieldSet<nconstant> &const_vars,
                                 FieldSet<nauxiliary> &auxiliary_vars,
                                 FieldSet<nprognostic> &solution) override {

    yakl::timer_start("linsolve");
    auto N = realsize(rhs);
    for (int i = 0; i < N; ++i) {
      this->eigen_rhs(i) = getindex(rhs, i);
    }

    this->eigen_sol = this->lu.solve(this->eigen_rhs);

    for (int i = 0; i < N; ++i) {
      setindex(solution, i, this->eigen_sol(i));
    }
    
    //int pis = primal_topology.is;
    //int pjs = primal_topology.js;
    //int pks = primal_topology.ks;

    //int dis = dual_topology.is;
    //int djs = dual_topology.js;
    //int dks = dual_topology.ks;

    //auto vvar = solution.fields_arr[VVAR].data;
    //auto wvar = solution.fields_arr[WVAR].data;

    //auto densvar = solution.fields_arr[DENSVAR].data;
    //auto denstendvar = rhs.fields_arr[DENSVAR].data;
    //auto refdens0var = const_vars.fields_arr[REFDENS0VAR].data;

    //auto uvar = auxiliary_vars.fields_arr[UVAR].data;
    //auto uwvar = auxiliary_vars.fields_arr[UWVAR].data;

    //auto fvar = auxiliary_vars.fields_arr[FVAR].data;
    //auto fwvar = auxiliary_vars.fields_arr[FWVAR].data;
    //
    //auto refhevar = const_vars.fields_arr[REFHEVAR].data;
    //auto refhewvar = const_vars.fields_arr[REFHEWVAR].data;
    //auto refdensreconvar = const_vars.fields_arr[REFDENSRECONVAR].data;
    //auto refdensvertreconvar = const_vars.fields_arr[REFDENSVERTRECONVAR].data;
    //
    //parallel_for(
    //    "Compute Uvar",
    //    SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      compute_Hext<1, diff_ord>(uvar, vvar, this->primal_geometry,
    //                                this->dual_geometry, dis, djs, dks, i, j, k,
    //                                n);
    //    });

    //parallel_for(
    //    "Compute UWVar",
    //    SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      compute_Hv<1, vert_diff_ord>(uwvar, wvar, this->primal_geometry,
    //                                   this->dual_geometry, dis, djs, dks, i, j,
    //                                   k + 1, n);
    //    });

    //auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
    //this->aux_exchange->exchanges_arr[UVAR].exchange_field(
    //    auxiliary_vars.fields_arr[UVAR]);
    //this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
    //    auxiliary_vars.fields_arr[UWVAR]);
    //this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
    //    auxiliary_vars.fields_arr[DENS0VAR]);

    //parallel_for(
    //    "Compute F",
    //    SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      Hk.compute_F_and_he(fvar, refhevar, uvar, refdens0var, dis, djs, dks, i,
    //                          j, k, n);
    //    });

    //parallel_for(
    //    "Compute FWvar",
    //    SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      Hk.compute_Fw_and_he(fwvar, refhewvar, uwvar, refdens0var, dis, djs, dks,
    //                           i, j, k + 1, n);
    //    });

    //auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    //this->aux_exchange->exchanges_arr[FVAR].exchange_field(
    //    auxiliary_vars.fields_arr[FVAR]);
    //this->aux_exchange->exchanges_arr[FWVAR].exchange_field(
    //    auxiliary_vars.fields_arr[FWVAR]);
    //this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
    //    auxiliary_vars.fields_arr[HEVAR]);
    //this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(
    //    auxiliary_vars.fields_arr[HEWVAR]);
    //
    //parallel_for(
    //    "Compute Dens Tend",
    //    SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      compute_wDbar2<ndensity>(denstendvar, refdensreconvar, fvar, dis, djs,
    //                               dks, i, j, k, n);
    //      real term2 = denstendvar(0, k + dks, j + djs, i + dis, n);
    //      compute_wDvbar<ndensity>(
    //          denstendvar, refdensvertreconvar, fwvar, dis, djs, dks, i, j, k, n);
    //      real term3 = denstendvar(0, k + dks, j + djs, i + dis, n);
    //      term2 *= dt / 2;
    //      term3 *= dt / 2;
    //      real term1 = densvar(0, k + dks, j + djs, i + dis, n); 
    //      //denstendvar(0, k + dks, j + djs, i + dis, n) *= dt / 2;
    //      //denstendvar(0, k + dks, j + djs, i + dis, n) += densvar(0, k + dks, j + djs, i + dis, n);
    //      if (i == 0) {
    //      std::cout << k << " " << term1 << " " << term2 << " " << term3 << " " << term1 + term2 + term3 << std::endl;
    //      }
    //    });

    // fourier
    auto refhevar = const_vars.fields_arr[REFHEVAR].data;
    auto refhewvar = const_vars.fields_arr[REFHEWVAR].data;
    auto refdensreconvar = const_vars.fields_arr[REFDENSRECONVAR].data;
    auto bdens0var = const_vars.fields_arr[BDENS0VAR].data;
    auto refdensvertreconvar = const_vars.fields_arr[REFDENSVERTRECONVAR].data;
    auto refdens0var = const_vars.fields_arr[REFDENS0VAR].data;

    //auto denstendvar = rhs.fields_arr[DENSVAR].data;


    //
    //auto refhevar = const_vars.fields_arr[REFHEVAR].data;
    //auto refhewvar = const_vars.fields_arr[REFHEWVAR].data;
    //auto refdensreconvar = const_vars.fields_arr[REFDENSRECONVAR].data;
    //auto refdensvertreconvar = const_vars.fields_arr[REFDENSVERTRECONVAR].data;

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

    
    parallel_for(
        "Compute Uvar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hext<1, diff_ord>(uvar, sol_v, this->primal_geometry,
                                    this->dual_geometry, dis, djs, dks, i, j, k,
                                    n);
        });

    parallel_for(
        "Compute UWVar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hv<1, vert_diff_ord>(uwvar, sol_w, this->primal_geometry,
                                       this->dual_geometry, dis, djs, dks, i, j,
                                       k + 1, n);
        });

    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[UWVAR].exchange_field(
        auxiliary_vars.fields_arr[UWVAR]);
    this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
        auxiliary_vars.fields_arr[DENS0VAR]);

    parallel_for(
        "Compute F",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F_and_he(fvar, refhevar, uvar, refdens0var, dis, djs, dks, i,
                              j, k, n);
        });

    parallel_for(
        "Compute FWvar",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_Fw_and_he(fwvar, refhewvar, uwvar, refdens0var, dis, djs, dks,
                               i, j, k + 1, n);
        });

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    this->aux_exchange->exchanges_arr[FVAR].exchange_field(
        auxiliary_vars.fields_arr[FVAR]);
    this->aux_exchange->exchanges_arr[FWVAR].exchange_field(
        auxiliary_vars.fields_arr[FWVAR]);
    this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
        auxiliary_vars.fields_arr[HEVAR]);
    this->aux_exchange->exchanges_arr[HEWVAR].exchange_field(
        auxiliary_vars.fields_arr[HEWVAR]);
    
    parallel_for(
        "Compute Dens Tend",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDvbar<ndensity>(
              denstend, refdensvertreconvar, fwvar, dis, djs, dks, i, j, k, n);
        });
    
    parallel_for(
        "Compute Dens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Iext<ndensity, diff_ord, vert_diff_ord>(
              rhs_dens0, denstend, this->primal_geometry, this->dual_geometry,
              pis, pjs, pks, i, j, k, n);
        });
    parallel_for(
        "Compute Bvar",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {

          real rho = rhs_dens0(0, pks + k, pjs + j, pis + i, n);
          real Tht = rhs_dens0(1, pks + k, pjs + j, pis + i, n);
          
          real b0_rho = bdens0var(0, pks + k, pjs + j, pis + i, n);
          real b0_Tht = bdens0var(1, pks + k, pjs + j, pis + i, n);
          real b1_rho = bdens0var(2, pks + k, pjs + j, pis + i, n);
          real b1_Tht = bdens0var(3, pks + k, pjs + j, pis + i, n);

          real b0 = rho * b0_rho + Tht * b0_Tht;
          real b1 = rho * b1_rho + Tht * b1_Tht;

          bvar(0, pks + k, pjs + j, pis + i, n) = -dtf2 * b0;
          bvar(1, pks + k, pjs + j, pis + i, n) = -dtf2 * b1;
        });
    this->aux_exchange->exchanges_arr[BVAR].exchange_field(
        auxiliary_vars.fields_arr[BVAR]);
    parallel_for(
        "Compute Wtend",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDv<ndensity>(wtend2, refdensvertreconvar, bvar, pis, pjs, pks,
                                i, j, k, n);
        });

    parallel_for(
        "Compute Dens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Iext<ndensity, diff_ord, vert_diff_ord>(
              rhs_dens0, rhs_dens, this->primal_geometry, this->dual_geometry,
              pis, pjs, pks, i, j, k, n);
        });
    parallel_for(
        "Compute Bvar",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {

          real rho = rhs_dens0(0, pks + k, pjs + j, pis + i, n);
          real Tht = rhs_dens0(1, pks + k, pjs + j, pis + i, n);
          
          real b0_rho = bdens0var(0, pks + k, pjs + j, pis + i, n);
          real b0_Tht = bdens0var(1, pks + k, pjs + j, pis + i, n);
          real b1_rho = bdens0var(2, pks + k, pjs + j, pis + i, n);
          real b1_Tht = bdens0var(3, pks + k, pjs + j, pis + i, n);

          real b0 = rho * b0_rho + Tht * b0_Tht;
          real b1 = rho * b1_rho + Tht * b1_Tht;

          bvar(0, pks + k, pjs + j, pis + i, n) = -dtf * b0;
          bvar(1, pks + k, pjs + j, pis + i, n) = -dtf * b1;
        });
    
    this->aux_exchange->exchanges_arr[BVAR].exchange_field(
        auxiliary_vars.fields_arr[BVAR]);


    parallel_for(
        "Compute Wtend",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDv<ndensity>(wtend, refdensvertreconvar, bvar, pis, pjs, pks,
                                i, j, k, n);
        });
    parallel_for(
        "Compute Vtend",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD1<ndensity>(vtend, refdensreconvar, bvar, pis, pjs, pks, i,
                                j, k, n);
        });


    real scale = 1.0 / (n_cells_x * n_cells_y);
    pocketfft::shape_t shape_v(5);
    pocketfft::shape_t shape_w(5);
    pocketfft::shape_t shape_dens(5);
    pocketfft::stride_t stride_v(5);
    pocketfft::stride_t stride_w(5);
    pocketfft::stride_t stride_dens(5);

    shape_v[0] = 1;
    shape_v[1] = ni;
    shape_v[2] = n_cells_y;
    shape_v[3] = n_cells_x;
    shape_v[4] = nens;
    
    shape_w[0] = 1;
    shape_w[1] = nl;
    shape_w[2] = n_cells_y; 
    shape_w[3] = n_cells_x; 
    shape_w[4] = nens;      
    
    shape_dens[0] = 2;
    shape_dens[1] = dual_topology.nl;
    shape_dens[2] = n_cells_y; 
    shape_dens[3] = n_cells_x; 
    shape_dens[4] = nens;      

    stride_v[4] = sizeof(complex);
    stride_w[4] = sizeof(complex);
    stride_dens[4] = sizeof(complex);
    for (int i = 3; i >= 0; i--) {
      stride_v[i] = stride_v[i + 1] * shape_v[i + 1];
      stride_w[i] = stride_w[i + 1] * shape_w[i + 1];
      stride_dens[i] = stride_dens[i + 1] * shape_dens[i + 1];
    }
    pocketfft::shape_t axes = {2, 3};
    
    parallel_for(
        "fft copy dens",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          for (int d = 0; d < ndensity; ++d) {
            complex_dens(d, k, j, i, n) = sol_dens(d, k + dks, j + djs, i + dis, n);
            complex_densrhs(d, k, j, i, n) = rhs_dens(d, k + dks, j + djs, i + dis, n);
          }
        });


    parallel_for(
        "fft copy v",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          complex_v(0, k, j, i, n) = sol_v(0, k + pks, j + pjs, i + pis, n);
          complex_vrhs(0, k, j, i, n) = rhs_v(0, k + pks, j + pjs, i + pis, n) +
                                        vtend(0, k + pks, j + pjs, i + pis, n);     
        });
    
    parallel_for(
        "fft copy w",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          complex_w(0, k, j, i, n) = sol_w(0, k + pks, j + pjs, i + pis, n);
          complex_wrhs(0, k, j, i, n) = rhs_w(0, k + pks, j + pjs, i + pis, n) +
                                        wtend(0, k + pks, j + pjs, i + pis, n); 
          complex_wrhs2(0, k, j, i, n) = wtend2(0, k + pks, j + pjs, i + pis, n); 
        });

    pocketfft::c2c(shape_v, stride_v, stride_v, axes, pocketfft::FORWARD,
                   complex_v.data(), complex_v.data(), 1._fp);
    pocketfft::c2c(shape_v, stride_v, stride_v, axes, pocketfft::FORWARD,
                   complex_vrhs.data(), complex_vrhs.data(), 1._fp);

    pocketfft::c2c(shape_w, stride_w, stride_w, axes, pocketfft::FORWARD,
                   complex_w.data(), complex_w.data(), 1._fp);
    pocketfft::c2c(shape_w, stride_w, stride_w, axes, pocketfft::FORWARD,
                   complex_wrhs.data(), complex_wrhs.data(), 1._fp);
    pocketfft::c2c(shape_w, stride_w, stride_w, axes, pocketfft::FORWARD,
                   complex_wrhs2.data(), complex_wrhs2.data(), 1._fp);
    
    pocketfft::c2c(shape_dens, stride_dens, stride_dens, axes, pocketfft::FORWARD,
                   complex_dens.data(), complex_dens.data(), 1._fp);
    pocketfft::c2c(shape_dens, stride_dens, stride_dens, axes, pocketfft::FORWARD,
                   complex_densrhs.data(), complex_densrhs.data(), 1._fp);
    
    parallel_for(
        "eh 1",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          complex_v_halo(0, k + pks, j + pjs, i + pis, n) = complex_v(0, k, j, i, n);
          complex_vrhs_halo(0, k + pks, j + pjs, i + pis, n) = complex_vrhs(0, k, j, i, n);
        });
    
    parallel_for(
        "eh 2",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          complex_w_halo(0, k + pks, j + pjs, i + pis, n) = complex_w(0, k, j, i, n);
          complex_wrhs_halo(0, k + pks, j + pjs, i + pis, n) = complex_wrhs(0, k, j, i, n);
        });


    parallel_for(
        "eh 3",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hv<1, vert_diff_ord>(complex_uw_halo, complex_w_halo, primal_geometry,
                                       dual_geometry, dis, djs, dks, i, j,
                                       k + 1, n);
          real hew = refhewvar(0, k + 1 + dks, j + djs, i + dis, n);
          complex_uw_halo(0, k + 1 + dks, j + djs, i + dis, n) *= hew;
        });
    
    parallel_for(
        "eh 5",
        SimpleBounds<3>(dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int j, int i, int n) {
          complex_uw_halo(0, 0 + dks, j + djs, i + dis, n) = 0;
          complex_uw_halo(0, dual_topology.ni + dks, j + djs, i + dis, n) = 0;
        });

    parallel_for(
        "eh 4",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDvbar<ndensity>(
              complex_dens_halo, refdensvertreconvar, complex_uw_halo, dis, djs, dks,
              i, j, k, n);
        });

    
    parallel_for(
        "compute vcoeff",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
            real b0_rho = bdens0var(0, pks + k, pjs + j, pis + i, n);
            real b0_Tht = bdens0var(1, pks + k, pjs + j, pis + i, n);
            real b1_rho = bdens0var(2, pks + k, pjs + j, pis + i, n);
            real b1_Tht = bdens0var(3, pks + k, pjs + j, pis + i, n);

            real refdens0recon = refdensreconvar(0, pks + k, pjs + j, pis + i, n);
            real refdens1recon = refdensreconvar(1, pks + k, pjs + j, pis + i, n);


            SArray<real, 1, ndims> fD1Dbar;
            fourier_cwD1Dbar2(fD1Dbar, 1.0, i, j, k, dual_topology.n_cells_x, dual_topology.n_cells_y,
                dual_topology.ni); 

            real fI = fourier_Iext<diff_ord>(primal_geometry, dual_geometry, pis, pjs, pks,
                                  i, j, k, 0, n_cells_x, n_cells_y, dual_topology.ni);

            SArray<real, 1, ndims> fH;
            fourier_Hext<diff_ord>(fH, primal_geometry, dual_geometry, pis, pjs, pks,
                              i, j, k, 0, n_cells_x, n_cells_y, dual_topology.ni);

            complex im(0, 1);
            complex fac = (2 * pi * i) / dual_topology.n_cells_x;
            //complex fD1 = exp(im * fac) - 1._fp;
            complex fD1 = 1._fp - exp(-im * fac);
            real he = refhevar(0, dks + k, djs + j, dis + i, n);

            real c1 = 1 - dtf2 * fI * fH(0) * fD1Dbar(0) * (
                           he * refdens0recon * refdens0recon * b0_rho
                         + he * refdens0recon * refdens1recon * b0_Tht
                         + he * refdens1recon * refdens0recon * b1_rho
                         + he * refdens1recon * refdens1recon * b1_Tht);

            complex c2 = -fD1 * dtf2 * fI * (refdens0recon * b0_rho + refdens1recon * b1_rho);
            complex c3 = -fD1 * dtf2 * fI * (refdens0recon * b0_Tht + refdens1recon * b1_Tht);
            
            complex c4 = -dtf * fD1 * fI * (refdens0recon * b0_rho + refdens1recon * b1_rho); 
            complex c5 = -dtf * fD1 * fI * (refdens0recon * b0_Tht + refdens1recon * b1_Tht); 

            complex lhs = c1 * complex_v(0, k, j, i, n) + 
                            c2 * complex_dens_halo(0, k + dks, j + djs, i + dis, n) +
                            c3 * complex_dens_halo(1, k + dks, j + djs, i + dis, n);
            //complex rhs = complex_vrhs(0, k, j, i, n) +
            //                c4 * complex_densrhs(0, k, j, i, n) +
            //                c5 * complex_densrhs(1, k, j, i, n);
            complex rhs = complex_vrhs(0, k, j, i, n);


            complex_vcoeff_halo(0, k + pks, j + pjs, i + pis, n) = rhs / c1;
            complex_vcoeff_halo(1, k + pks, j + pjs, i + pis, n) = -c2 / c1;
            complex_vcoeff_halo(2, k + pks, j + pjs, i + pis, n) = -c3 / c1;

      });
    
    parallel_for(
        "TEST",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {

          real fI_k = fourier_Iext<diff_ord>(primal_geometry, dual_geometry, pis, pjs, pks,
                                i, j, k, 0, n_cells_x, n_cells_y, dual_topology.ni);
          real fI_kp1 = fourier_Iext<diff_ord>(primal_geometry, dual_geometry, pis, pjs, pks,
                                i, j, k + 1, 0, n_cells_x, n_cells_y, dual_topology.ni);

          real alpha_kp1 = refdensvertreconvar(0, pks + k + 1, pjs + j, pis + i, n); 
          real beta_kp1 = fI_kp1 * bdens0var(0, pks + k + 1, pjs + j, pis + i, n);  
          real beta_k = fI_k * bdens0var(0, pks + k, pjs + j, pis + i, n);

          real gamma_kp2 = refdensvertreconvar(0, pks + k + 2, pjs + j, pis + i, n);
          real gamma_kp1 = refdensvertreconvar(0, pks + k + 1, pjs + j, pis + i, n);
          real gamma_k = refdensvertreconvar(0, pks + k, pjs + j, pis + i, n);

          complex uw_kp2 = complex_uw_halo(0, k + 2 + dks, j + djs, i + dis, n);
          complex uw_kp1 = complex_uw_halo(0, k + 1 + dks, j + djs, i + dis, n);
          complex uw_k = complex_uw_halo(0, k + dks, j + djs, i + dis, n);

          complex term1 = alpha_kp1 * (beta_kp1 * gamma_kp2 * uw_kp2 -
                          gamma_kp1 * (beta_kp1 + beta_k) * uw_kp1 + beta_k * gamma_k * uw_k);
          
          alpha_kp1 = refdensvertreconvar(0, pks + k + 1, pjs + j, pis + i, n); 
          beta_kp1 = fI_kp1 * bdens0var(1, pks + k + 1, pjs + j, pis + i, n);  
          beta_k = fI_k * bdens0var(1, pks + k, pjs + j, pis + i, n);

          gamma_kp2 = refdensvertreconvar(1, pks + k + 2, pjs + j, pis + i, n);
          gamma_kp1 = refdensvertreconvar(1, pks + k + 1, pjs + j, pis + i, n);
          gamma_k = refdensvertreconvar(1, pks + k, pjs + j, pis + i, n);
          
          complex term2 = alpha_kp1 * (beta_kp1 * gamma_kp2 * uw_kp2 -
                          gamma_kp1 * (beta_kp1 + beta_k) * uw_kp1 + beta_k * gamma_k * uw_k);
          
          alpha_kp1 = refdensvertreconvar(1, pks + k + 1, pjs + j, pis + i, n); 
          beta_kp1 = fI_kp1 * bdens0var(2, pks + k + 1, pjs + j, pis + i, n);  
          beta_k = fI_k *  bdens0var(2, pks + k, pjs + j, pis + i, n);

          gamma_kp2 = refdensvertreconvar(0, pks + k + 2, pjs + j, pis + i, n);
          gamma_kp1 = refdensvertreconvar(0, pks + k + 1, pjs + j, pis + i, n);
          gamma_k = refdensvertreconvar(0, pks + k, pjs + j, pis + i, n);
          
          complex term3 = alpha_kp1 * (beta_kp1 * gamma_kp2 * uw_kp2 -
                          gamma_kp1 * (beta_kp1 + beta_k) * uw_kp1 + beta_k * gamma_k * uw_k);
          
          alpha_kp1 = refdensvertreconvar(1, pks + k + 1, pjs + j, pis + i, n); 
          beta_kp1 = fI_kp1 * bdens0var(3, pks + k + 1, pjs + j, pis + i, n);  
          beta_k = fI_k * bdens0var(3, pks + k, pjs + j, pis + i, n);

          gamma_kp2 = refdensvertreconvar(1, pks + k + 2, pjs + j, pis + i, n);
          gamma_kp1 = refdensvertreconvar(1, pks + k + 1, pjs + j, pis + i, n);
          gamma_k = refdensvertreconvar(1, pks + k, pjs + j, pis + i, n);
          
          complex term4 = alpha_kp1 * (beta_kp1 * gamma_kp2 * uw_kp2 -
                          gamma_kp1 * (beta_kp1 + beta_k) * uw_kp1 + beta_k * gamma_k * uw_k);

          complex_wrhs2(0, k, j, i, n) = -dtf2 * (term1 + term2 + term3 + term4);
    });

    parallel_for(
        "TEST",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
            //complex c0 = complex_vcoeff_halo(0, pks + k, pjs + j, pis + i, n);
            //complex c1 = complex_vcoeff_halo(1, pks + k, pjs + j, pis + i, n);
            //complex c2 = complex_vcoeff_halo(2, pks + k, pjs + j, pis + i, n);


            //complex lhs = complex_v(0, k, j, i, n);
            //complex rhs = c0 +
            //                c1 * complex_dens_halo(0, k + pks, j + pjs, i + pis, n) +
            //                c2 * complex_dens_halo(1, k + pks, j + pjs, i + pis, n);
            //if (i == 4) {
            //  std::cout << k << " " << lhs - rhs << std::endl;
            //}
            
            complex lhs = complex_w(0, k, j, i, n) + complex_wrhs2(0, k, j, i, n);
            complex rhs = complex_wrhs(0, k, j, i, n);
            if (i == 0) {
              std::cout << k << " " << lhs << " " << rhs << " " << lhs - rhs << std::endl;
            }
      });

    exit(1);

    // prepare rhs

    //int dis = dual_topology.is;
    //int djs = dual_topology.js;
    //int dks = dual_topology.ks;

    //int pis = primal_topology.is;
    //int pjs = primal_topology.js;
    //int pks = primal_topology.ks;
    //
    //auto bdens0var = const_vars.fields_arr[BDENS0VAR].data;
    //auto refhevar = const_vars.fields_arr[REFHEVAR].data;
    //auto refhewvar = const_vars.fields_arr[REFHEWVAR].data;
    //auto refdensedgereconvar = const_vars.fields_arr[REFDENSEDGERECONVAR].data;
    //auto refdensvertedgereconvar =
    //    const_vars.fields_arr[REFDENSVERTEDGERECONVAR].data;
    //auto refdensreconvar = const_vars.fields_arr[REFDENSRECONVAR].data;
    //auto refdensvertreconvar = const_vars.fields_arr[REFDENSVERTRECONVAR].data;

    //auto rhs_densvar = rhs.fields_arr[DENSVAR].data;
    //auto rhs_vvar = rhs.fields_arr[VVAR].data;
    //auto rhs_wvar = rhs.fields_arr[WVAR].data;

    //auto bvar = auxiliary_vars.fields_arr[BVAR].data;
    //auto dens0var = auxiliary_vars.fields_arr[DENS0VAR].data;

    //parallel_for(
    //    "linsolve prepare rhs 1",
    //    SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
    //                    primal_topology.n_cells_x, primal_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      compute_Iext<ndensity, diff_ord, vert_diff_ord>(
    //          dens0var, rhs_densvar, primal_geometry, dual_geometry,
    //          pis, pjs, pks, i, j, k, n);
    //    });

    //this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
    //    auxiliary_vars.fields_arr[DENS0VAR]);

    //parallel_for(
    //    "linsolve prepare rhs 2",
    //    SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
    //                    primal_topology.n_cells_x, primal_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {

    //      real rho = dens0var(0, pks + k, pjs + j, pis + i, n);
    //      real Tht = dens0var(1, pks + k, pjs + j, pis + i, n);
    //      
    //      real b0_rho = bdens0var(0, pks + k, pjs + j, pis + i, n);
    //      real b0_Tht = bdens0var(1, pks + k, pjs + j, pis + i, n);
    //      real b1_rho = bdens0var(2, pks + k, pjs + j, pis + i, n);
    //      real b1_Tht = bdens0var(3, pks + k, pjs + j, pis + i, n);

    //      real b0 = rho * b0_rho + Tht * b0_Tht;
    //      real b1 = rho * b1_rho + Tht * b1_Tht;

    //      bvar(0, pks + k, pjs + j, pis + i, n) = b0;
    //      bvar(1, pks + k, pjs + j, pis + i, n) = b1;
    //    });

    //this->aux_exchange->exchanges_arr[BVAR].exchange_field(
    //    auxiliary_vars.fields_arr[BVAR]);

    //parallel_for(
    //    "linsolve prepare rhs 3",
    //    SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
    //                    primal_topology.n_cells_x, primal_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      real tmp = rhs_wvar(0, pks + k, pjs + j, pis + i, n);
    //      compute_wDv<ndensity>(rhs_wvar, refdensvertreconvar, bvar, pis, pjs, pks,
    //                            i, j, k, n);
    //      rhs_wvar(0, pks + k, pjs + j, pis + i, n) *= -dt / 2;
    //      rhs_wvar(0, pks + k, pjs + j, pis + i, n) += tmp;
    //    });

    //parallel_for(
    //    "linsolve prepare rhs 4",
    //    SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
    //                    primal_topology.n_cells_x, primal_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //      real tmp = rhs_vvar(0, pks + k, pjs + j, pis + i, n);
    //      compute_wD1<ndensity>(rhs_vvar, refdensreconvar, bvar, pis, pjs, pks, i,
    //                            j, k, n);
    //      rhs_vvar(0, pks + k, pjs + j, pis + i, n) *= -dt / 2;
    //      rhs_vvar(0, pks + k, pjs + j, pis + i, n) += tmp;
    //    });
    //

    //real dtf = dt * dt / 4;

    //// compute c1, c2, c3 in utilde = c1 * F_u + c2 * Dvbar(rho_0 * w) * c3 * Dvbar * (s_0 * w)
    //parallel_for(
    //    "compute c1, c2, c3",
    //    SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
    //                    dual_topology.n_cells_x, dual_topology.nens),
    //    YAKL_LAMBDA(int k, int j, int i, int n) {
    //        real b0_rho = bdens0var(0, pks + k, pjs + j, pis + i, n);
    //        real b0_Tht = bdens0var(1, pks + k, pjs + j, pis + i, n);
    //        real b1_rho = bdens0var(2, pks + k, pjs + j, pis + i, n);
    //        real b1_Tht = bdens0var(3, pks + k, pjs + j, pis + i, n);

    //        real he = hevar(d, dks + k, djs + j, dis + i)
    //        real tmp = 1 - dtf * he * refdens0recon * refdens0recon * b0_rho
    //                     - dtf * he * refdens0recon * refdens1recon * b0_Tht
    //                     - dtf * he * refdens1recon * refdens0recon * b1_rho
    //                     - dtf * he * refdens1recon * refdens1recon * b1_Tht;
    //        
    //        c1var(0, k + dks + 1, j + djs, i + dis, n) =
    //      }
    //    });
    

    yakl::timer_stop("linsolve");
  }
};

// *******   Statistics   ***********//

class ModelStats : public Stats {
public:
  real3d TEarr, KEarr, PEarr, IEarr, PVarr, PENSarr, trimmed_density;

  void initialize(ModelParameters &params, Parallel &par,
                  const Topology &primal_topo, const Topology &dual_topo,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    Stats::initialize(params, par, primal_topo, dual_topo, primal_geom,
                      dual_geom);
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

    this->TEarr =
        real3d("TE", this->dual_topology.nl, this->dual_topology.n_cells_y,
               this->dual_topology.n_cells_x);
    this->KEarr =
        real3d("KE", this->dual_topology.nl, this->dual_topology.n_cells_y,
               this->dual_topology.n_cells_x);
    this->IEarr =
        real3d("IE", this->dual_topology.nl, this->dual_topology.n_cells_y,
               this->dual_topology.n_cells_x);
    this->PEarr =
        real3d("PE", this->dual_topology.nl, this->dual_topology.n_cells_y,
               this->dual_topology.n_cells_x);
    this->PVarr =
        real3d("PV", this->primal_topology.nl, this->primal_topology.n_cells_y,
               this->primal_topology.n_cells_x);
    this->PENSarr = real3d("PENS", this->primal_topology.nl,
                           this->primal_topology.n_cells_y,
                           this->primal_topology.n_cells_x);
    this->trimmed_density =
        real3d("trimmed_density", this->dual_topology.nl,
               this->dual_topology.n_cells_y, this->dual_topology.n_cells_x);
  }

  void compute(FieldSet<nprognostic> &progvars, FieldSet<nconstant> &constvars,
               int tind) {

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
  aux_topo_arr[FVAR2] = dtopo;
  aux_topo_arr[UVAR] = dtopo;
  aux_topo_arr[HEVAR] = dtopo;
  aux_topo_arr[FWVAR] = dtopo;
  aux_topo_arr[FWVAR2] = dtopo;
  aux_topo_arr[UWVAR] = dtopo;
  aux_topo_arr[HEWVAR] = dtopo;
  aux_topo_arr[KVAR] = dtopo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[FVAR2] = "F2";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  aux_names_arr[FWVAR] = "Fw";
  aux_names_arr[FWVAR2] = "Fw2";
  aux_names_arr[UWVAR] = "Uw";
  aux_names_arr[HEWVAR] = "hew";
  set_dofs_arr(aux_ndofs_arr, BVAR, 0, 0, ndensity);  // B = straight (0,0)-form
  set_dofs_arr(aux_ndofs_arr, KVAR, ndims, 1, 1);     // K = twisted (n,1)-form
  set_dofs_arr(aux_ndofs_arr, FVAR, ndims - 1, 1, 1); // F = twisted
                                                      // (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, FVAR2, ndims - 1, 1, 1); // F = twisted
                                                       // (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, UVAR, ndims - 1, 1, 1);  // U = twisted
                                                       // (n-1,1)-form
  set_dofs_arr(aux_ndofs_arr, HEVAR, ndims - 1, 1,
               1); // he lives on horiz dual edges, associated with F
  set_dofs_arr(aux_ndofs_arr, FWVAR, ndims, 0, 1);  // Fw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, FWVAR2, ndims, 0, 1); // Fw = twisted (n,0)-form
  set_dofs_arr(aux_ndofs_arr, UWVAR, ndims, 0, 1);  // Uw = twisted (n,0)-form
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

  const_topo_arr[REFDENSVAR] = dtopo;
  const_names_arr[REFDENSVAR] = "refdens";
  set_dofs_arr(const_ndofs_arr, REFDENSVAR, ndims, 1,
               ndensity); // refdens = twisted (n,1)-form
  const_topo_arr[REFDENS0VAR] = ptopo;
  const_names_arr[REFDENS0VAR] = "refdens0";
  set_dofs_arr(const_ndofs_arr, REFDENS0VAR, 0, 0,
               ndensity); // refdens0 = straight (0,0)-form
  const_topo_arr[REFNSQ0VAR] = ptopo;
  const_names_arr[REFNSQ0VAR] = "refnsq0";
  set_dofs_arr(const_ndofs_arr, REFNSQ0VAR, 0, 0,
               ndensity); // refnsq0 = straight (0,0)-form
  
  const_topo_arr[BDENS0VAR] = ptopo;
  const_names_arr[BDENS0VAR] = "bdens0";
  set_dofs_arr(const_ndofs_arr, BDENS0VAR, 0, 0,
               ndensity * ndensity); // refdens0 = straight (0,0)-form

  const_topo_arr[REFHEVAR] = dtopo;
  const_names_arr[REFHEVAR] = "refhe";
  set_dofs_arr(const_ndofs_arr, REFHEVAR, ndims - 1, 1,
               1); // he lives on horiz dual edges, associated with F
  const_topo_arr[REFHEWVAR] = dtopo;
  const_names_arr[REFHEWVAR] = "refhew";
  set_dofs_arr(const_ndofs_arr, REFHEWVAR, ndims, 0,
               1); // hew lives on vert dual edges, associated with Fw

  const_topo_arr[REFDENSRECONVAR] = dtopo;
  const_topo_arr[REFDENSEDGERECONVAR] = dtopo;
  const_topo_arr[REFDENSVERTRECONVAR] = dtopo;
  const_topo_arr[REFDENSVERTEDGERECONVAR] = dtopo;
  const_names_arr[REFDENSVERTRECONVAR] = "refdensvertrecon";
  const_names_arr[REFDENSVERTEDGERECONVAR] = "refdensvertedgerecon";
  const_names_arr[REFDENSRECONVAR] = "refdensrecon";
  const_names_arr[REFDENSEDGERECONVAR] = "refdensedgerecon";
  set_dofs_arr(
      const_ndofs_arr, REFDENSRECONVAR, ndims - 1, 1,
      ndensity); // densrecon lives on horiz dual edges, associated with F
  set_dofs_arr(
      const_ndofs_arr, REFDENSEDGERECONVAR, ndims, 1,
      2 * ndims *
          ndensity); // densedgerecon lives on dual cells, associated with F
  set_dofs_arr(
      const_ndofs_arr, REFDENSVERTRECONVAR, ndims, 0,
      ndensity); // densvertrecon lives on vert dual edges, associated with Fw
  set_dofs_arr(
      const_ndofs_arr, REFDENSVERTEDGERECONVAR, ndims, 1,
      2 * ndensity); // densedgerecon lives on dual cells, associated with Fw
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
                              ExchangeSet<nconstant> &const_exchange,
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

  using T::refentropicdensity_f;
  using T::refnsq_f;
  using T::refrho_f;

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.zlen = Lz;
    params.xc = xc;
    params.zc = zc;
  }

  void add_diagnostics(
      std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {
    T::add_diagnostics(diagnostics);
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              ExchangeSet<nconstant> &const_exchange,
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
   
    // reference state stuff

    auto bdens0var = constvars.fields_arr[BDENS0VAR].data;
    auto refdensvar = constvars.fields_arr[REFDENSVAR].data;
    auto refnsq0var = constvars.fields_arr[REFNSQ0VAR].data;
    auto refdens0var = constvars.fields_arr[REFDENS0VAR].data;

    auto refhevar = constvars.fields_arr[REFHEVAR].data;
    auto refhewvar = constvars.fields_arr[REFHEWVAR].data;
    auto refdensedgereconvar = constvars.fields_arr[REFDENSEDGERECONVAR].data;
    auto refdensvertedgereconvar =
        constvars.fields_arr[REFDENSVERTEDGERECONVAR].data;
    auto refdensreconvar = constvars.fields_arr[REFDENSRECONVAR].data;
    auto refdensvertreconvar = constvars.fields_arr[REFDENSVERTRECONVAR].data;

    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) { return refrho_f(x, z, thermo); },
        constvars.fields_arr[REFDENSVAR], 0);
    dual_geom.set_11form_values(
        YAKL_LAMBDA(real x, real z) {
          return refentropicdensity_f(x, z, thermo);
        },
        constvars.fields_arr[REFDENSVAR], 1);

    primal_geom.set_00form_values(
        YAKL_LAMBDA(real x, real z) { return refnsq_f(x, z, thermo); },
        constvars.fields_arr[REFNSQ0VAR], 0);

    auto primal_topology = primal_geom.topology;
    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;
    
    auto dual_topology = dual_geom.topology;
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;
    
    const_exchange.exchanges_arr[REFDENSVAR].exchange_field(
        constvars.fields_arr[REFDENSVAR]);
    

    parallel_for(
        "Compute refdens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Iext<ndensity, diff_ord, vert_diff_ord>(
              refdens0var,
              refdensvar, primal_geom, dual_geom,
              pis, pjs, pks, i, j, k, n);
        });
    
    const_exchange.exchanges_arr[REFDENS0VAR].exchange_field(
        constvars.fields_arr[REFDENS0VAR]);
    
    parallel_for(
        "Compute refhe",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_he(refhevar, refdens0var, dis, djs, dks, i,
                              j, k, n);
        });

    parallel_for(
        "Compute refhew",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_hew(refhewvar, refdens0var, dis, djs, dks,
                               i, j, k + 1, n);
        });
    
    const_exchange.exchanges_arr[REFHEVAR].exchange_field(
        constvars.fields_arr[REFHEVAR]);
    const_exchange.exchanges_arr[REFHEWVAR].exchange_field(
        constvars.fields_arr[REFHEWVAR]);
    
    parallel_for(
        "Compute Bdens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          real Rd = thermo.cst.Rd;
          real pr = thermo.cst.pr;
          real gamma_d = thermo.cst.gamma_d;
          real Cpd = thermo.cst.Cpd;
          real Cvd = thermo.cst.Cvd;
          real grav = 9.80616_fp;
          real grav2 = grav * grav;

          real rho_ref = refdens0var(0, pks + k, pjs + j, pis + i, n);
          real Tht_ref = refdens0var(1, pks + k, pjs + j, pis + i, n);
          real tht_ref = Tht_ref / rho_ref;

          real rho_ref2 = rho_ref * rho_ref;
          real p_ref = pr * pow(Rd * Tht_ref / pr, gamma_d);
          real dpdtht_ref = gamma_d * p_ref / tht_ref;
          real dpdtht_ref2 = dpdtht_ref * dpdtht_ref;
          real Nref2 = refnsq0var(0, pks + k, pjs + j, pis + i, n);
          real cref2 = gamma_d * p_ref / rho_ref;

          real b0_rho = (cref2 * rho_ref - dpdtht_ref * tht_ref) / rho_ref2;
          real b0_tht = dpdtht_ref / rho_ref -
                         dpdtht_ref2 * tht_ref / (cref2 * rho_ref2) -
                         dpdtht_ref2 * grav2 * tht_ref /
                             (Nref2 * cref2 * cref2 * rho_ref2);
          real b0_Tht = b0_tht / rho_ref;
          b0_rho -= tht_ref / rho_ref * b0_tht;

          real b1_rho = dpdtht_ref / rho_ref2;
          real b1_tht = dpdtht_ref2 * (Nref2 * cref2 + grav2) /
                         (Nref2 * cref2 * cref2 * rho_ref2);
          real b1_Tht = b1_tht / rho_ref;
          b1_rho -= tht_ref / rho_ref * b1_tht;

          bdens0var(0, pks + k, pjs + j, pis + i, n) = b0_rho;
          bdens0var(1, pks + k, pjs + j, pis + i, n) = b0_Tht;
          bdens0var(2, pks + k, pjs + j, pis + i, n) = b1_rho;
          bdens0var(3, pks + k, pjs + j, pis + i, n) = b1_Tht;
        });

    constexpr int ord = 1;
    constexpr auto recon_type = RECONSTRUCTION_TYPE::CFV;

    SArray<real, 3, ord, ord, ord> wenoRecon;
    SArray<real, 2, ord, 2> to_gll;
    SArray<real, 1, (ord-1)/2 + 2> wenoIdl;
    real wenoSigma;
    
    parallel_for(
        "ComputeDensEdgeRecon",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_edge_recon<ndensity, recon_type,
                                     ord>(
              refdensedgereconvar, refdens0var, dis, djs, dks, i, j, k, n,
              wenoRecon, to_gll, wenoIdl, wenoSigma);
          compute_twisted_vert_edge_recon<ndensity,
                                          recon_type,
                                          ord>(
              refdensvertedgereconvar, refdens0var, dis, djs, dks, i, j, k, n,
              wenoRecon, to_gll, wenoIdl, wenoSigma);
        });

    const_exchange.exchanges_arr[REFDENSEDGERECONVAR].exchange_field(
      constvars.fields_arr[REFDENSEDGERECONVAR]);
    const_exchange.exchanges_arr[REFDENSVERTEDGERECONVAR].exchange_field(
       constvars.fields_arr[REFDENSVERTEDGERECONVAR]);

    // not used
    auto uvar = refdensvar;
    auto uwvar = refdensvar;

    parallel_for(
        "ComputeDens recon",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_recon<ndensity, recon_type>(
              refdensreconvar, refdensedgereconvar, uvar, dis, djs, dks, i, j, k, n);

          // scale twisted recons
          for (int d = 0; d < ndims; d++) {
            for (int l = 0; l < ndensity; l++) {
              refdensreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) =
                  refdensreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) /
                  refhevar(d, k + dks, j + djs, i + dis, n);
            }
          }
        });

    parallel_for(
        "ComputeDensVertRECON",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_vert_recon<ndensity, recon_type>(
              refdensvertreconvar, refdensvertedgereconvar, uwvar, dis, djs, dks, i,
              j, k + 1, n);

          // scale twisted recons
          for (int l = 0; l < ndensity; l++) {
            refdensvertreconvar(l, k + dks + 1, j + djs, i + dis, n) =
                refdensvertreconvar(l, k + dks + 1, j + djs, i + dis, n) /
                refhewvar(0, k + dks + 1, j + djs, i + dis, n);
          }
        });

    const_exchange.exchanges_arr[REFDENSRECONVAR].exchange_field(
        constvars.fields_arr[REFDENSRECONVAR]);
    const_exchange.exchanges_arr[REFDENSVERTRECONVAR].exchange_field(
        constvars.fields_arr[REFDENSVERTRECONVAR]);

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
                              ExchangeSet<nconstant> &const_exchange,
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

real YAKL_INLINE isothermal_zdep(real x, real z, real var_s, real T_ref, real g,
                                 const ThermoPotential &thermo) {
  real Rd = thermo.cst.Rd;
  real delta = g / (Rd * T_ref);
  return var_s * exp(-delta * z);
}

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
  static real constexpr T_ref = 300.0_fp;
  
  static real YAKL_INLINE refnsq_f(real x, real z,
                                   const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real gamma_d = thermo.cst.gamma_d;
    real N2 = (gamma_d - 1) / gamma_d * g * g / (Rd * T_ref);
    return N2;
  }

  static real YAKL_INLINE refrho_f(real x, real z,
                                   const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real rho_s = thermo.cst.pr / (Rd * T_ref);
    real rho_ref = isothermal_zdep(x, z, rho_s, T_ref, g, thermo);
    return rho_ref;
  }

  static real YAKL_INLINE refp_f(real x, real z,
                                 const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real rho_ref = refrho_f(x, z, thermo);
    real p_ref = Rd * rho_ref * T_ref;
    return p_ref;
  }

  static real YAKL_INLINE refentropicdensity_f(real x, real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(x, z, thermo);
    real p_ref = refp_f(x, z, thermo);
    return rho_ref * thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  }


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
    //real r = sqrt((z - bzc) * (z - bzc));
    //real dtheta = (r < rc) ? dss * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
    real dtheta = 0;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var(p, T + dT, 0, 0, 0, 0);
  }
  
  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
    //diagnostics.emplace_back(std::make_unique<ExactDensityDiagnostic>());
    //diagnostics.emplace_back(std::make_unique<BackgroundDensityDiagnostic>());
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


struct GravityWave {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 300e3_fp;
  static real constexpr Lz = 10e3_fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr d = 5e3_fp;
  static real constexpr T_ref = 250._fp;
  static real constexpr u_0 = 0._fp;
  static real constexpr x_c = 150e3_fp;
  static real constexpr p_s = 1e5_fp;
  // static real constexpr dT_max = 0.01_fp;
  static real constexpr dT_max = 1.0_fp;

  static real YAKL_INLINE refnsq_f(real x, real z,
                                   const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real gamma_d = thermo.cst.gamma_d;
    real N2 = (gamma_d - 1) / gamma_d * g * g / (Rd * T_ref);
    return N2;
  }

  static real YAKL_INLINE refrho_f(real x, real z,
                                   const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real rho_s = p_s / (Rd * T_ref);
    real rho_ref = isothermal_zdep(x, z, rho_s, T_ref, g, thermo);
    return rho_ref;
  }

  static real YAKL_INLINE refp_f(real x, real z,
                                 const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;
    real rho_ref = refrho_f(x, z, thermo);
    real p_ref = Rd * rho_ref * T_ref;
    return p_ref;
  }

  static real YAKL_INLINE refentropicdensity_f(real x, real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(x, z, thermo);
    real p_ref = refp_f(x, z, thermo);
    return rho_ref * thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;

    real delta = g / (Rd * T_ref);
    real rho_s = p_s / (Rd * T_ref);

    real rho_ref = refrho_f(x, z, thermo);

    real dT_b = dT_max * exp(-pow((x - x_c) / d, 2)) * sin(pi * z / Lz);
    real dT = exp(delta * z / 2) * dT_b;

    real drho_b = -rho_s * dT_b / T_ref;
    real drho = exp(-delta * z / 2) * drho_b;

    real T = (T_ref + dT);
    real rho = (rho_ref + drho);

    return rho;
    // return rho_ref;
  }

  static real YAKL_INLINE rhoexact_f(real x, real z, real t,
                                     const ThermoPotential &thermo) {
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
        real L_m1 = (-cos(alpha * t) / alpha2 + cos(beta * t) / beta2) * fac1 +
                    1 / (alpha2 * beta2);
        real L_0 = (sin(alpha * t) / alpha - sin(beta * t) / beta) * fac1;
        real L_1 = (cos(alpha * t) - cos(beta * t)) * fac1;
        real L_2 = (-alpha * sin(alpha * t) + alpha * sin(beta * t)) * fac1;
        real L_3 = (-alpha2 * cos(alpha * t) + beta2 * cos(beta * t)) * fac1;

        if (alpha == 0) {
          L_m1 = (beta2 * t * t - 1 + cos(beta * t)) / (beta2 * beta2);
          L_0 = (beta * t - sin(beta * t)) / (beta2 * beta);
        }

        complex drhot_b0 = -rho_s / T_ref * dT_max / sqrt(pi) * d / Lx *
                           exp(-d * d * k_x2 / 4) * exp(-im * k_x * x_c) * k_z *
                           Lz / (2._fp * im);

        complex drhot_b =
            (L_3 + (p_1 + g * (im * k_z - delta / 2)) * L_1 +
             (c_s2 * (k_z2 + delta2 / 4) + g * (im * k_z - delta / 2)) * f2 *
                 L_m1) *
            drhot_b0;

        complex expfac = exp(im * (k_x * xp + k_z * z));

        drho_b += drhot_b * expfac;
      }
    }

    real drho = exp(-delta * z / 2) * drho_b.real();
    real rho = refrho_f(x, z, thermo) + drho;

    return rho;
  }

  static real YAKL_INLINE entropicdensityexact_f(
      real x, real z, real t, const ThermoPotential &thermo) {
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
        real L_m1 = (-cos(alpha * t) / alpha2 + cos(beta * t) / beta2) * fac1 +
                    1 / (alpha2 * beta2);
        real L_0 = (sin(alpha * t) / alpha - sin(beta * t) / beta) * fac1;
        real L_1 = (cos(alpha * t) - cos(beta * t)) * fac1;
        real L_2 = (-alpha * sin(alpha * t) + alpha * sin(beta * t)) * fac1;
        real L_3 = (-alpha2 * cos(alpha * t) + beta2 * cos(beta * t)) * fac1;

        if (alpha == 0) {
          L_m1 = (beta2 * t * t - 1 + cos(beta * t)) / (beta2 * beta2);
          L_0 = (beta * t - sin(beta * t)) / (beta2 * beta);
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

        complex expfac = exp(im * (k_x * xp + k_z * z));

        drho_b += drhot_b * expfac;
        dp_b += dpt_b * expfac;
      }
    }

    complex dT_b = T_ref * (dp_b / p_s - drho_b / rho_s);

    real p_ref = isothermal_zdep(x, z, p_s, T_ref, g, thermo);
    real dp = exp(-delta * z / 2) * dp_b.real();
    real dT = exp(delta * z / 2) * dT_b.real();

    real drho = exp(-delta * z / 2) * drho_b.real();
    real rho = refrho_f(x, z, thermo) + drho;

    return rho *
           thermo.compute_entropic_var(p_ref + dp, T_ref + dT, 0, 0, 0, 0);
  }

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;

    real delta = g / (Rd * T_ref);
    real rho_s = p_s / (Rd * T_ref);

    real rho_ref = isothermal_zdep(x, z, rho_s, T_ref, g, thermo);

    real dT_b = dT_max * exp(-pow((x - x_c) / d, 2)) * sin(pi * z / Lz);
    real dT = exp(delta * z / 2) * dT_b;

    real drho_b = -rho_s * dT_b / T_ref;
    real drho = exp(-delta * z / 2) * drho_b;

    real T = (T_ref + dT);
    real rho = (rho_ref + drho);

    // real p = rho * Rd * T;
    real p_ref = isothermal_zdep(x, z, p_s, T_ref, g, thermo);
    real dp = Rd * T_ref * drho + Rd * rho_ref * dT;
    real p = p_ref + dp;

    return thermo.compute_entropic_var(p, T, 0, 0, 0, 0);
    // return thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  }

  static real YAKL_INLINE entropicdensity_f(real x, real z,
                                            const ThermoPotential &thermo) {
    return rho_f(x, z, thermo) * entropicvar_f(x, z, thermo);
  }

  static vecext<2> YAKL_INLINE v_f(real x, real y) {
    vecext<2> vvec;
    vvec.u = u_0;
    vvec.w = 0;
    return vvec;
  }

  struct ExactDensityDiagnostic : public Diagnostic {

    void initialize(const Topology &ptopo, const Topology &dtopo,
                    const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "dense";
      topology = dtopo;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(ptopo, dtopo, pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return rhoexact_f(x, z, time, thermo) - refrho_f(x, z, thermo);
          },
          field, 0);

      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return entropicdensityexact_f(x, z, time, thermo) -
                   refentropicdensity_f(x, z, thermo);
          },
          field, 1);
    }
  };

  struct BackgroundDensityDiagnostic : public Diagnostic {

    void initialize(const Topology &ptopo, const Topology &dtopo,
                    const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "densb";
      topology = dtopo;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(ptopo, dtopo, pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) { return refrho_f(x, z, thermo); }, field,
          0);

      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return refentropicdensity_f(x, z, thermo);
          },
          field, 1);
    }
  };

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
    diagnostics.emplace_back(std::make_unique<ExactDensityDiagnostic>());
    diagnostics.emplace_back(std::make_unique<BackgroundDensityDiagnostic>());
  }
};

void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance) {
  if (name == "doublevortex") {
    // testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  } else if (name == "gravitywave") {
    testcase = std::make_unique<EulerTestCase<GravityWave>>();
  } else if (name == "risingbubble") {
      if (acoustic_balance) {
        testcase = std::make_unique<EulerTestCase<RisingBubble<true>>>();
      } else {
        testcase = std::make_unique<EulerTestCase<RisingBubble<false>>>();
      }
    //} else if (name == "moistrisingbubble") {
    //  testcase = std::make_unique<MoistEulerTestCase<MoistRisingBubble>>();
    //} else if (name == "largerisingbubble") {
    //  testcase = std::make_unique<EulerTestCase<LargeRisingBubble>>();
    //} else if (name == "moistlargerisingbubble") {
    //  testcase =
    //  std::make_unique<MoistEulerTestCase<MoistLargeRisingBubble>>();
  } else {
    throw std::runtime_error("unknown test case");
  }
}
