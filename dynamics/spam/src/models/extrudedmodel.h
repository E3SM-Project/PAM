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
                  ReferenceState &refstate) {

    ExtrudedTendencies::initialize(params, primal_geom, dual_geom, refstate);
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
                         FieldSet<nprognostic> &x) override {}

  void compute_dens0(real5d dens0var, const real5d densvar) {

    const auto &primal_topology = primal_geometry.topology;

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

    const auto &dual_topology = dual_geometry.topology;

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

    const auto &dual_topology = dual_geometry.topology;

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

    const auto &dual_topology = dual_geometry.topology;

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

    const auto &dual_topology = dual_geometry.topology;

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

    const auto &dual_topology = dual_geometry.topology;

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
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void compute_B(real fac, real5d Bvar, const real5d Kvar, const real5d densvar,
                 const real5d HSvar) {

    const auto &primal_topology = primal_geometry.topology;

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
      real5d densedgereconvar, real5d densvertedgereconvar, real5d densvertedgereconvar2, 
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

          compute_twisted_vert_edge_recon<ndensity,
                                          RECONSTRUCTION_TYPE::CFV,
                                          dual_vert_reconstruction_order>(
              densvertedgereconvar2, dens0var, dis, djs, dks, i, j, k, n,
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
      real5d densreconvar, real5d densvertreconvar, real5d densvertreconvar2,real5d qxzreconvar,
      real5d qxzvertreconvar, real5d coriolisxzreconvar,
      real5d coriolisxzvertreconvar, const real5d densedgereconvar,
      const real5d densvertedgereconvar,
      const real5d densvertedgereconvar2,
      const real5d qxzedgereconvar,
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
          compute_twisted_vert_recon<ndensity, RECONSTRUCTION_TYPE::CFV>(
              densvertreconvar2, densvertedgereconvar2, UWvar, dis, djs, dks, i,
              j, k + 1, n);
          // scale twisted recons
          for (int l = 0; l < ndensity; l++) {
            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) =
                densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) /
                HEWvar(0, k + dks + 1, j + djs, i + dis, n);
          }
          for (int l = 0; l < ndensity; l++) {
            densvertreconvar2(l, k + dks + 1, j + djs, i + dis, n) =
                densvertreconvar2(l, k + dks + 1, j + djs, i + dis, n) /
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
                     const real5d densreconvar,
                     const real5d densvertreconvar,
                     const real5d densvertreconvar2,
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
          compute_wDv_fct<ndensity>(Wtendvar, densvertreconvar2, Phivertvar,
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
          compute_wDv_fct<ndensity>(Wtendvar, densvertreconvar2, Phivertvar,
                                    Bvar, pis, pjs, pks, i, j, 0, n);
          compute_wDv_fct<ndensity>(Wtendvar, densvertreconvar2, Phivertvar,
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

    std::cout << "wtend: " << yakl::intrinsics::maxval(yakl::intrinsics::abs(Wtendvar)) << std::endl;

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

  void compute_functional_derivatives(
      ADD_MODE addmode, real fac, real dt, FieldSet<nconstant> &const_vars,
      FieldSet<nprognostic> &x, FieldSet<nauxiliary> &auxiliary_vars) override {

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U(auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[VVAR].data);
    compute_UW(auxiliary_vars.fields_arr[UWVAR].data, x.fields_arr[WVAR].data);

    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({UVAR, UWVAR, DENS0VAR});

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
    auxiliary_vars.exchange({FVAR, FWVAR, KVAR});

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
    auxiliary_vars.exchange({BVAR});
  }

  void apply_symplectic(real dt, FieldSet<nconstant> &const_vars,
                        FieldSet<nprognostic> &x,
                        FieldSet<nauxiliary> &auxiliary_vars,
                        FieldSet<nprognostic> &xtend) override {

    const auto &dual_topology = dual_geometry.topology;

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U(auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[VVAR].data);
    compute_UW(auxiliary_vars.fields_arr[UWVAR].data, x.fields_arr[WVAR].data);
    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);

    auxiliary_vars.exchange({UVAR, UWVAR, DENS0VAR});

    compute_F_FW_and_he(
        auxiliary_vars.fields_arr[FVAR2].data,
        auxiliary_vars.fields_arr[FWVAR2].data,
        auxiliary_vars.fields_arr[HEVAR].data,
        auxiliary_vars.fields_arr[HEWVAR].data, x.fields_arr[VVAR].data,
        auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[WVAR].data,
        auxiliary_vars.fields_arr[UWVAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data);

    auxiliary_vars.fields_arr[FWVAR2].set_bnd(0.0);
    auxiliary_vars.exchange({FVAR2, FWVAR2, HEVAR, HEWVAR});

    compute_q0f0(auxiliary_vars.fields_arr[QXZ0VAR].data,
                 auxiliary_vars.fields_arr[FXZ0VAR].data,
                 x.fields_arr[VVAR].data, x.fields_arr[WVAR].data,
                 x.fields_arr[DENSVAR].data,
                 const_vars.fields_arr[CORIOLISXZVAR].data);

    auxiliary_vars.fields_arr[QXZ0VAR].set_bnd(0.0);
    auxiliary_vars.fields_arr[FXZ0VAR].set_bnd(0.0);
    auxiliary_vars.exchange({QXZ0VAR, FXZ0VAR});

    compute_FT_and_FTW(auxiliary_vars.fields_arr[FTVAR].data,
                       auxiliary_vars.fields_arr[FTWVAR].data,
                       auxiliary_vars.fields_arr[FVAR2].data,
                       auxiliary_vars.fields_arr[FWVAR2].data);

    auxiliary_vars.exchange({FTVAR, FTWVAR});

    // Compute densrecon, densvertrecon, qrecon and frecon
    compute_edge_reconstructions(
        auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR2].data,
        auxiliary_vars.fields_arr[QXZEDGERECONVAR].data,
        auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data,
        auxiliary_vars.fields_arr[QXZ0VAR].data,
        auxiliary_vars.fields_arr[FXZ0VAR].data);

    auxiliary_vars.exchange({DENSEDGERECONVAR, DENSVERTEDGERECONVAR, DENSVERTEDGERECONVAR2,
                             QXZEDGERECONVAR, QXZVERTEDGERECONVAR,
                             CORIOLISXZEDGERECONVAR,
                             CORIOLISXZVERTEDGERECONVAR});

    compute_recons(auxiliary_vars.fields_arr[DENSRECONVAR].data,
                   auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                   auxiliary_vars.fields_arr[DENSVERTRECONVAR2].data,
                   auxiliary_vars.fields_arr[QXZRECONVAR].data,
                   auxiliary_vars.fields_arr[QXZVERTRECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISXZRECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISXZVERTRECONVAR].data,
                   auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR2].data,
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

    auxiliary_vars.exchange({DENSRECONVAR, DENSVERTRECONVAR, QXZRECONVAR,
                             QXZVERTRECONVAR, CORIOLISXZRECONVAR,
                             CORIOLISXZVERTRECONVAR});

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
    auxiliary_vars.exchange({EDGEFLUXVAR, VERTEDGEFLUXVAR});

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

    auxiliary_vars.exchange({MFVAR});

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

    auxiliary_vars.exchange({PHIVAR, PHIVERTVAR});

    // Compute tendencies
    compute_tendencies(xtend.fields_arr[DENSVAR].data,
                       xtend.fields_arr[VVAR].data, xtend.fields_arr[WVAR].data,
                       auxiliary_vars.fields_arr[DENSRECONVAR].data,
                       auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                       auxiliary_vars.fields_arr[DENSVERTRECONVAR2].data,
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

struct ModelReferenceState : ReferenceState {
  real3d q_di;
  real3d q_pi;
  real3d rho_di;
  real3d rho_pi;
  real4d Blin_coeff;

  void initialize(const Topology &primal_topology,
                  const Topology &dual_topology) override {
    this->q_pi =
        real3d("refq_pi", ndensity, primal_topology.ni, primal_topology.nens);
    this->q_di =
        real3d("refq_di", ndensity, dual_topology.ni, dual_topology.nens);
    this->rho_pi =
        real3d("refrho_pi", 1, primal_topology.ni, primal_topology.nens);
    this->rho_di = real3d("refrho_di", 1, dual_topology.ni, dual_topology.nens);
    this->Blin_coeff = real4d("Blin coeff", ndensity, ndensity,
                              primal_topology.ni, primal_topology.nens);

    this->is_initialized = true;
  }
};

// *******   Linear system   ***********//
class ModelLinearSystem : public LinearSystem {
  ModelReferenceState *reference_state;

  yakl::RealFFT1D<real> fftv_x;
  yakl::RealFFT1D<real> fftw_x;
  // yakl::RealFFT1D<real> fftv_y;
  // yakl::RealFFT1D<real> fftw_y;

  int nxf, nyf;

  real4d v_transform;
  real4d w_transform;
  complex4d complex_vrhs;
  complex4d complex_wrhs;
  complex5d complex_vcoeff;

  complex4d tri_l;
  complex4d tri_d;
  complex4d tri_u;

  complex4d tri_c;

public:
  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  ReferenceState &refstate) override {

    LinearSystem::initialize(params, primal_geom, dual_geom, refstate);

    this->reference_state = static_cast<ModelReferenceState *>(&refstate);

    const auto &primal_topology = primal_geom.topology;

    auto pni = primal_topology.ni;
    auto pnl = primal_topology.nl;
    auto nx = primal_topology.n_cells_x;
    auto ny = primal_topology.n_cells_y;
    auto nens = primal_topology.nens;

    this->nxf = nx + 2 - nx % 2;
    this->nyf = ny + 2 - ny % 2;

    v_transform = real4d("v transform", pni, nyf, nxf, nens);
    w_transform = real4d("w transform", pnl, nyf, nxf, nens);
    complex_vrhs = complex4d("complex vrhs", pni, ny, nx, nens);
    complex_wrhs = complex4d("complex wrhs", pnl, ny, nx, nens);

    fftv_x.init(v_transform, 2, nx);
    fftw_x.init(w_transform, 2, nx);
    // fftv_y.init(v_transform, 1, ny);
    // fftw_y.init(w_transform, 1, ny);

    complex_vcoeff =
        complex5d("complex vcoeff", 1 + ndensity, pni, ny, nx, nens);

    tri_d = complex4d("tri d", pnl, ny, nx, nens);
    tri_l = complex4d("tri l", pnl, ny, nx, nens);
    tri_u = complex4d("tri u", pnl, ny, nx, nens);
    tri_c = complex4d("tri c", pnl, ny, nx, nens);
  }

  virtual void compute_coefficients(real dt) override {
    auto &refstate = *this->reference_state;

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

    parallel_for(
        "compute vcoeff",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> fD1Dbar;
          fourier_cwD1Dbar2(fD1Dbar, 1, i, j, k, dual_topology.n_cells_x,
                            dual_topology.n_cells_y, dual_topology.ni);

          real fI = fourier_Iext<diff_ord>(primal_geometry, dual_geometry, pis,
                                           pjs, pks, i, j, k, 0, n_cells_x,
                                           n_cells_y, dual_topology.ni);
          SArray<real, 1, ndims> fH;
          fourier_Hext<diff_ord>(fH, primal_geometry, dual_geometry, pis, pjs,
                                 pks, i, j, k, 0, n_cells_x, n_cells_y,
                                 dual_topology.ni);
          SArray<complex, 1, ndims> fD1;
          fourier_cwD1(fD1, 1, i, j, k, dual_topology.n_cells_x,
                       dual_topology.n_cells_y, dual_topology.ni);

          real he = refstate.rho_pi(0, k, n);

          real c1 = 1;
          for (int d1 = 0; d1 < ndensity; ++d1) {
            for (int d2 = 0; d2 < ndensity; ++d2) {
              c1 -= dtf2 * fI * fH(0) * fD1Dbar(0) * he *
                    refstate.q_pi(d1, k, n) * refstate.q_pi(d2, k, n) *
                    refstate.Blin_coeff(d1, d2, k, n);
            }
          }

          complex_vcoeff(0, k, j, i, n) = 1 / c1;
          for (int d1 = 0; d1 < ndensity; ++d1) {
            complex cd1 = 0;
            for (int d2 = 0; d2 < ndensity; ++d2) {
              cd1 += fD1(0) * dtf2 * fI * refstate.q_pi(d2, k, n) *
                     refstate.Blin_coeff(d2, d1, k, n);
            }
            complex_vcoeff(1 + d1, k, j, i, n) = cd1 / c1;
          }
        });

    parallel_for(
        "Compute vertical tridiag",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real fI_k = fourier_Iext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fI_kp1 = fourier_Iext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          real gamma_fac_kp2 = refstate.rho_di(0, k + 2, n) *
                               Hv_coeff(primal_geometry, dual_geometry, pis,
                                        pjs, pks, i, j, k + 2);
          real gamma_fac_kp1 = refstate.rho_di(0, k + 1, n) *
                               Hv_coeff(primal_geometry, dual_geometry, pis,
                                        pjs, pks, i, j, k + 1);
          real gamma_fac_k =
              refstate.rho_di(0, k, n) *
              Hv_coeff(primal_geometry, dual_geometry, pis, pjs, pks, i, j, k);

          tri_u(k, j, i, n) = 0;
          tri_d(k, j, i, n) = 1;
          tri_l(k, j, i, n) = 0;

          for (int d1 = 0; d1 < ndensity; ++d1) {
            for (int d2 = 0; d2 < ndensity; ++d2) {
              real alpha_kp1 = refstate.q_di(d1, k + 1, n);

              real beta_kp1 = fI_kp1 * refstate.Blin_coeff(d1, d2, k + 1, n);
              real beta_k = fI_k * refstate.Blin_coeff(d1, d2, k, n);

              real gamma_kp2 = gamma_fac_kp2 * refstate.q_di(d2, k + 2, n);
              real gamma_kp1 = gamma_fac_kp1 * refstate.q_di(d2, k + 1, n);
              real gamma_k = gamma_fac_k * refstate.q_di(d2, k, n);

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
          real fI_k = fourier_Iext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fI_kp1 = fourier_Iext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          real gamma_fac_kp2 = refstate.rho_di(0, k + 2, n) *
                               Hv_coeff(primal_geometry, dual_geometry, pis,
                                        pjs, pks, i, j, k + 2);
          real gamma_fac_kp1 = refstate.rho_di(0, k + 1, n) *
                               Hv_coeff(primal_geometry, dual_geometry, pis,
                                        pjs, pks, i, j, k + 1);
          real gamma_fac_k =
              refstate.rho_di(0, k, n) *
              Hv_coeff(primal_geometry, dual_geometry, pis, pjs, pks, i, j, k);

          SArray<real, 1, ndims> fH_kp1_a;
          SArray<real, 1, ndims> fH_k_a;
          fourier_Hext<diff_ord>(fH_kp1_a, primal_geometry, dual_geometry, pis,
                                 pjs, pks, i, j, k + 1, 0, n_cells_x, n_cells_y,
                                 dual_topology.ni);
          fourier_Hext<diff_ord>(fH_k_a, primal_geometry, dual_geometry, pis,
                                 pjs, pks, i, j, k, 0, n_cells_x, n_cells_y,
                                 dual_topology.ni);
          real fHh_kp1 = fH_kp1_a(0);
          real fHh_k = fH_k_a(0);

          complex fDbar2_kp1 = fourier_Dbar2(1, i, j, k + 1, n_cells_x,
                                             n_cells_y, dual_topology.ni);
          complex fDbar2_k =
              fourier_Dbar2(1, i, j, k, n_cells_x, n_cells_y, dual_topology.ni);

          real he_kp1 = refstate.rho_pi(0, k + 1, n);
          real he_k = refstate.rho_pi(0, k, n);

          for (int d1 = 0; d1 < ndensity; ++d1) {
            for (int d2 = 0; d2 < ndensity; ++d2) {
              for (int d3 = 0; d3 < ndensity; ++d3) {

                real alpha_kp1 = dtf2 * refstate.q_di(d1, k + 1, n);
                complex beta_kp1 =
                    fI_kp1 * refstate.Blin_coeff(d1, d2, k + 1, n) *
                    refstate.q_pi(d2, k + 1, n) * fDbar2_kp1 * he_kp1 * fHh_kp1;
                complex beta_k = fI_k * refstate.Blin_coeff(d1, d2, k, n) *
                                 refstate.q_pi(d2, k, n) * fDbar2_k * he_k *
                                 fHh_k;

                real gamma_kp2 = gamma_fac_kp2 * refstate.q_di(d3, k + 2, n);
                real gamma_kp1 = gamma_fac_kp1 * refstate.q_di(d3, k + 1, n);
                real gamma_k = gamma_fac_k * refstate.q_di(d3, k, n);

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

    auto &refstate = *this->reference_state;

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    yakl::timer_start("linsolve");

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

    parallel_for(
        "Prepare rhs 1 - I",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Iext<ndensity, diff_ord, vert_diff_ord>(
              sol_dens, rhs_dens, this->primal_geometry, this->dual_geometry,
              pis, pjs, pks, i, j, k, n);
        });

    parallel_for(
        "Prepare rhs 2 - B",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          for (int d1 = 0; d1 < ndensity; ++d1) {
            real b_d1 = 0;
            for (int d2 = 0; d2 < ndensity; ++d2) {
              b_d1 -= dtf * refstate.Blin_coeff(d1, d2, k, n) *
                      sol_dens(d2, pks + k, pjs + j, pis + i, n);
            }
            bvar(d1, pks + k, pjs + j, pis + i, n) = b_d1;
          }
        });

    auxiliary_vars.exchange({BVAR});

    parallel_for(
        "Prepare rhs 3 - wDv",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDv<ndensity>(sol_w, refstate.q_di, bvar, pis, pjs, pks, i, j,
                                k, n);
        });
    parallel_for(
        "Prepare rhs 4 - wD1",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD1<ndensity>(sol_v, refstate.q_pi, bvar, pis, pjs, pks, i, j,
                                k, n);
        });

    parallel_for(
        "Load v_transform",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          v_transform(k, j, i, n) = rhs_v(0, k + pks, j + pjs, i + pis, n) +
                                    sol_v(0, k + pks, j + pjs, i + pis, n);
        });

    parallel_for(
        "Load w_transform",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          w_transform(k, j, i, n) = rhs_w(0, k + pks, j + pjs, i + pis, n) +
                                    sol_w(0, k + pks, j + pjs, i + pis, n);
        });

    yakl::timer_start("fft fwd");
    fftv_x.forward_real(v_transform);
    fftw_x.forward_real(w_transform);
    // fftv_y.forward_real(v_transform);
    // fftw_y.forward_real(w_transform);
    yakl::timer_stop("fft fwd");

    parallel_for(
        "Compute complex vrhs",
        yakl::c::Bounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                           {0, nxf - 1, 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real v_real = v_transform(k, j, i, n);
          real v_imag = v_transform(k, j, i + 1, n);
          complex_vrhs(k, j, i / 2, n) = complex(v_real, v_imag);
        });

    parallel_for(
        "Compute complex wrhs",
        Bounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                  {0, nxf - 1, 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real w_real = w_transform(k, j, i, n);
          real w_imag = w_transform(k, j, i + 1, n);
          complex_wrhs(k, j, i / 2, n) = complex(w_real, w_imag);
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

          real fI_k = fourier_Iext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fI_kp1 = fourier_Iext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          SArray<real, 1, ndims> fH_kp1_a;
          SArray<real, 1, ndims> fH_k_a;
          fourier_Hext<diff_ord>(fH_kp1_a, primal_geometry, dual_geometry, pis,
                                 pjs, pks, i, j, k + 1, 0, n_cells_x, n_cells_y,
                                 dual_topology.ni);
          fourier_Hext<diff_ord>(fH_k_a, primal_geometry, dual_geometry, pis,
                                 pjs, pks, i, j, k, 0, n_cells_x, n_cells_y,
                                 dual_topology.ni);
          real fHh_kp1 = fH_kp1_a(0);
          real fHh_k = fH_k_a(0);

          complex fDbar2_kp1 = fourier_Dbar2(1, i, j, k + 1, n_cells_x,
                                             n_cells_y, dual_topology.ni);
          complex fDbar2_k =
              fourier_Dbar2(1, i, j, k, n_cells_x, n_cells_y, dual_topology.ni);

          real he_kp1 = refstate.rho_pi(0, k + 1, n);
          real he_k = refstate.rho_pi(0, k, n);

          for (int d1 = 0; d1 < ndensity; ++d1) {
            for (int d2 = 0; d2 < ndensity; ++d2) {
              real alpha_kp1 = dtf2 * refstate.q_di(d1, k + 1, n);
              complex beta_kp1 =
                  fI_kp1 * refstate.Blin_coeff(d1, d2, k + 1, n) *
                  refstate.q_pi(d2, k + 1, n) * fDbar2_kp1 * he_kp1 * fHh_kp1;
              complex beta_k = fI_k * refstate.Blin_coeff(d1, d2, k, n) *
                               refstate.q_pi(d2, k, n) * fDbar2_k * he_k *
                               fHh_k;
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

          real gamma_fac_kp1 = refstate.rho_di(0, k + 1, n) *
                               Hv_coeff(primal_geometry, dual_geometry, pis,
                                        pjs, pks, i, j, k + 1);
          real gamma_fac_k =
              refstate.rho_di(0, k, n) *
              Hv_coeff(primal_geometry, dual_geometry, pis, pjs, pks, i, j, k);

          complex_vrhs(k, j, i, n) *= complex_vcoeff(0, k, j, i, n);
          for (int d1 = 0; d1 < ndensity; ++d1) {
            complex_vrhs(k, j, i, n) +=
                complex_vcoeff(1 + d1, k, j, i, n) *
                (gamma_fac_kp1 * refstate.q_di(d1, k + 1, n) * w_kp1 -
                 gamma_fac_k * refstate.q_di(d1, k, n) * w_k);
          }
        });

    parallel_for(
        "Store complex vrhs",
        yakl::c::Bounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                           {0, nxf - 1, 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          v_transform(k, j, i, n) = complex_vrhs(k, j, i / 2, n).real();
          v_transform(k, j, i + 1, n) = complex_vrhs(k, j, i / 2, n).imag();
        });

    parallel_for(
        "Store complex wrhs",
        yakl::c::Bounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                           {0, nxf - 1, 2}, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          w_transform(k, j, i, n) = complex_wrhs(k, j, i / 2, n).real();
          w_transform(k, j, i + 1, n) = complex_wrhs(k, j, i / 2, n).imag();
        });

    yakl::timer_start("fft bwd");
    fftv_x.inverse_real(v_transform);
    fftw_x.inverse_real(w_transform);
    // fftv_y.inverse_real(v_transform);
    // fftw_y.inverse_real(w_transform);
    yakl::timer_stop("fft bwd");

    parallel_for(
        "Store w",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          sol_w(0, k + pks, j + pjs, i + pis, n) = w_transform(k, j, i, n);
        });
    parallel_for(
        "Store v",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          sol_v(0, k + pks, j + pjs, i + pis, n) = v_transform(k, j, i, n);
        });

    solution.exchange({VVAR, WVAR});

    // back out densities

    parallel_for(
        "Revover densities 1 - Hext",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hext<1, diff_ord>(uvar, sol_v, this->primal_geometry,
                                    this->dual_geometry, dis, djs, dks, i, j, k,
                                    n);
        });

    parallel_for(
        "Recover densities 2 - Hv",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Hv<1, vert_diff_ord>(uwvar, sol_w, this->primal_geometry,
                                       this->dual_geometry, dis, djs, dks, i, j,
                                       k + 1, n);
        });

    auxiliary_vars.fields_arr[UWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({UVAR, UWVAR});

    YAKL_SCOPE(Hk, ::Hk);
    parallel_for(
        "Recover densities 3 - F",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F(fvar, uvar, refstate.rho_pi, dis, djs, dks, i, j, k, n);
        });

    parallel_for(
        "Recover densities 4 - Fw",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_Fw(fwvar, uwvar, refstate.rho_di, dis, djs, dks, i, j,
                        k + 1, n);
        });

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({FVAR, FWVAR});

    parallel_for(
        "Recover densities 5 - finish",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDbar2<ndensity>(sol_dens, refstate.q_pi, fvar, dis, djs, dks,
                                   i, j, k, n);
          compute_wDvbar<ndensity, ADD_MODE::ADD>(
              sol_dens, refstate.q_di, fwvar, dis, djs, dks, i, j, k, n);
          for (int d = 0; d < ndensity; ++d) {
            sol_dens(d, k + dks, j + djs, i + dis, n) *= -dt / 2;
            sol_dens(d, k + dks, j + djs, i + dis, n) +=
                rhs_dens(d, pks + k, pjs + j, pis + i, n);
          }
        });

    yakl::timer_stop("linsolve");
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
void initialize_variables(
    const Topology &ptopo, const Topology &dtopo,
    std::array<FieldDescription, nprognostic> &prog_desc_arr,
    std::array<FieldDescription, nconstant> &const_desc_arr,
    std::array<FieldDescription, nauxiliary> &aux_desc_arr) {

  // primal grid represents straight quantities, dual grid twisted quantities
  // ndims is the BASEDIM size!

  // v, w, dens
  prog_desc_arr[VVAR] = {"v", ptopo, 1, 0, 1}; // v = straight (1,0)-form
  prog_desc_arr[WVAR] = {"w", ptopo, 0, 1, 1}; // w = straight (0,1)-form
  prog_desc_arr[DENSVAR] = {"dens", dtopo, ndims, 1,
                            ndensity}; // dens = twisted (n,1)-form

  // hs
  const_desc_arr[HSVAR] = {"hs", dtopo, ndims, 1, 1}; // hs = twisted (n,1)-form
  const_desc_arr[CORIOLISXZVAR] = {"coriolisxz", ptopo, 1, 1,
                                   1}; // f = straight (1,1)-form

  // functional derivatives = F, B, K, he, U
  aux_desc_arr[FVAR] = {"F", dtopo, ndims - 1, 1, 1};   // F = twisted
                                                        // (n-1,1)-form
  aux_desc_arr[FVAR2] = {"F2", dtopo, ndims - 1, 1, 1}; // F2 = twisted
                                                        // (n-1,1)-form
  aux_desc_arr[BVAR] = {"B", ptopo, 0, 0, ndensity}; // B = straight (0,0)-form

  aux_desc_arr[KVAR] = {"K", dtopo, ndims, 1, 1}; // K = twisted (n,1)-form

  aux_desc_arr[HEVAR] = {"he", dtopo, ndims - 1, 1,
                         1}; // he lives on horiz dual edges, associated with F

  aux_desc_arr[UVAR] = {"U", dtopo, ndims - 1, 1, 1}; // U = twisted
                                                      // (n-1,1)-form

  aux_desc_arr[FWVAR] = {"Fw", dtopo, ndims, 0, 1}; // Fw = twisted (n,0)-form
  aux_desc_arr[FWVAR2] = {"Fw2", dtopo, ndims, 0,
                          1}; // Fw2 = twisted (n,0)-form

  aux_desc_arr[HEWVAR] = {
      "hew", dtopo, ndims, 0,
      1}; // hew lives on vert dual edges, associated with Fw

  aux_desc_arr[UWVAR] = {"Uw", dtopo, ndims, 0, 1}; // Uw = twisted (n,0)-form

  // dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_desc_arr[DENS0VAR] = {"dens0", ptopo, 0, 0,
                            ndensity}; // dens0 = straight (0,0)-form
  aux_desc_arr[DENSEDGERECONVAR] = {
      "densedgerecon", dtopo, ndims, 1,
      2 * ndims *
          ndensity}; // densedgerecon lives on dual cells, associated with F
  aux_desc_arr[DENSRECONVAR] = {
      "densrecon", dtopo, ndims - 1, 1,
      ndensity}; // densrecon lives on horiz dual edges, associated with F
  aux_desc_arr[DENSVERTEDGERECONVAR] = {
      "densvertedgerecon", dtopo, ndims, 1,
      2 * ndensity}; // densedgerecon lives on dual cells, associated with Fw
  aux_desc_arr[DENSVERTEDGERECONVAR2] = {
      "densvertedgerecon", dtopo, ndims, 1,
      2 * ndensity}; // densedgerecon lives on dual cells, associated with Fw
  aux_desc_arr[DENSVERTRECONVAR] = {
      "densvertrecon", dtopo, ndims, 0,
      ndensity}; // densvertrecon lives on vert dual edges, associated with Fw
  aux_desc_arr[DENSVERTRECONVAR2] = {
      "densvertrecon", dtopo, ndims, 0,
      ndensity}; // densvertrecon lives on vert dual edges, associated with Fw

  // fct stuff- Phi, Mf, edgeflux
  aux_desc_arr[PHIVAR] = {"Phi", dtopo, ndims - 1, 1, ndensity};
  aux_desc_arr[MFVAR] = {"Mf", dtopo, ndims, 1, ndensity};
  aux_desc_arr[EDGEFLUXVAR] = {"edgeflux", dtopo, ndims - 1, 1, ndensity};
  aux_desc_arr[PHIVERTVAR] = {"PhiVert", dtopo, ndims, 0, ndensity};
  aux_desc_arr[VERTEDGEFLUXVAR] = {"vertedgeflux", dtopo, ndims, 0, ndensity};

  // Q stuff
  aux_desc_arr[QXZ0VAR] = {"QXZ0", dtopo, 0, 0, 1}; // Q0 = twisted (0,0)-form
  aux_desc_arr[QXZRECONVAR] = {
      "qxzrecon", ptopo, 0, 1,
      1}; // qxzrecon lives on vert primal edges, associated with w
  aux_desc_arr[QXZEDGERECONVAR] = {
      "qxzedgerecon", ptopo, ndims, 1,
      2 * 1}; // qxzedgerecon lives on primal cells, associated with Fw/w
  aux_desc_arr[QXZVERTRECONVAR] = {
      "qxzvertrecon", ptopo, 1, 0,
      1}; // qxzsvertrecon lives on horiz primal edges, associated with v
  aux_desc_arr[QXZVERTEDGERECONVAR] = {
      "qxzvertedgerecon", ptopo, ndims, 1,
      2 * 1}; // qxzvertedgerecon lives on primal cells, associated with F/v
  aux_desc_arr[QXZFLUXVAR] = {
      "qxzflux", ptopo, 0, 1,
      1}; // qxzflux lives on vert primal edges, associated with w
  aux_desc_arr[QXZVERTFLUXVAR] = {
      "qxzvertflux", ptopo, 1, 0,
      1}; // qxzvertflux lives on horiz primal edges, associated with v

  aux_desc_arr[FTVAR] = {"FT", ptopo, 1, 0,
                         1}; // FT = straight (1,0)-form ie Fw at v pts
  aux_desc_arr[FTWVAR] = {"FTW", ptopo, 0, 1,
                          1}; // FTW = straight (0,1)-form ie F at w pts

  aux_desc_arr[FXZ0VAR] = {"fxz0", dtopo, 0, 0,
                           1}; // fxz0 is a twisted 0,0 form
  aux_desc_arr[CORIOLISXZRECONVAR] = {
      "coriolisxzrecon", ptopo, 0, 1,
      1}; // coriolisxzrecon lives on vert primal edges, associated with w
  aux_desc_arr[CORIOLISXZEDGERECONVAR] = {
      "coriolisxzedgerecon", ptopo, ndims, 1,
      2 * 1}; // coriolisxzedgerecon lives on primal cells, associated with Fw/w
  aux_desc_arr[CORIOLISXZVERTRECONVAR] = {
      "coriolisxzvertrecon", ptopo, 1, 0,
      1}; // coriolisxzsvertrecon lives on horiz primal edges, associated with v
  aux_desc_arr[CORIOLISXZVERTEDGERECONVAR] = {
      "coriolisxzvertedgerecon", ptopo, ndims, 1,
      2 * 1}; // coriolisxzvertedgerecon lives on primal cells,
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

  using T::refentropicdensity_f;
  using T::refnsq_f;
  using T::refrho_f;
  using T::v_f;

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

    primal_geom.set_10form_values(
        YAKL_LAMBDA(real x, real y) { return v_f(x, y); },
        progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    primal_geom.set_01form_values(
        YAKL_LAMBDA(real x, real y) { return v_f(x, y); },
        progvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);

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

  void set_reference_state(ReferenceState &reference_state,
                           const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = static_cast<ModelReferenceState &>(reference_state);

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    YAKL_SCOPE(thermo, ::thermo);
    parallel_for(
        "Compute refstate rho_pi/q_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          real z = k * primal_geom.dz;
          refstate.rho_pi(0, k, n) = refrho_f(z, thermo);
          refstate.q_pi(0, k, n) =
              refrho_f(z, thermo) / refstate.rho_pi(0, k, n);
          refstate.q_pi(1, k, n) =
              refentropicdensity_f(z, thermo) / refstate.rho_pi(0, k, n);
        });

    parallel_for(
        "Compute refstate rho_di/q_di",
        SimpleBounds<2>(dual_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          real dz = dual_geom.dz;
          real z;
          if (k == 0) {
            z = 0;
          } else if (k == dual_topology.ni - 1) {
            z = dual_geom.Lz;
          } else {
            z = (k - 0.5_fp) * dual_geom.dz;
          }

          refstate.rho_di(0, k, n) = refrho_f(z, thermo);
          refstate.q_di(0, k, n) =
              refrho_f(z, thermo) / refstate.rho_di(0, k, n);
          refstate.q_di(1, k, n) =
              refentropicdensity_f(z, thermo) / refstate.rho_di(0, k, n);
        });

    parallel_for(
        "Compute Blin_coeff",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          real z = k * primal_geom.dz;

          real Rd = thermo.cst.Rd;
          real pr = thermo.cst.pr;
          real gamma_d = thermo.cst.gamma_d;
          real Cpd = thermo.cst.Cpd;
          real Cvd = thermo.cst.Cvd;
          real grav = g;
          real grav2 = grav * grav;

          real rho_ref = refrho_f(z, thermo);
          real alpha_ref = 1 / rho_ref;
          real S_ref = refentropicdensity_f(z, thermo);
          real s_ref = S_ref / rho_ref;

          real rho_ref2 = rho_ref * rho_ref;
          real p_ref = thermo.solve_p(rho_ref, s_ref, 0, 0, 0, 0);

          real dpds_ref =
              thermo.compute_dpdentropic_var(alpha_ref, s_ref, 0, 0, 0, 0);
          real dpds_ref2 = dpds_ref * dpds_ref;

          real Nref2 = refnsq_f(z, thermo);
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

          refstate.Blin_coeff(0, 0, k, n) = b0_rho;
          refstate.Blin_coeff(0, 1, k, n) = b0_S;
          refstate.Blin_coeff(1, 0, k, n) = b1_rho;
          refstate.Blin_coeff(1, 1, k, n) = b1_S;
        });
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
    //real dtheta = (r < rc) ? dss * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
    real dtheta = 0._fp;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var(p, T + dT, 0, 0, 0, 0);
  }

  static vecext<2> YAKL_INLINE v_f(real x, real y) {
    vecext<2> vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
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
    return rho_ref * thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  }

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

  static vecext<2> YAKL_INLINE v_f(real x, real y) {
    vecext<2> vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
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
    return rho_ref * thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;

    real delta = g / (Rd * T_ref);
    real rho_s = p_s / (Rd * T_ref);

    real rho_ref = refrho_f(z, thermo);

    real dT_b = dT_max * exp(-pow((x - x_c) / d, 2)) * sin(pi * z / Lz);
    real dT = exp(delta * z / 2) * dT_b;

    real drho_b = -rho_s * dT_b / T_ref;
    real drho = exp(-delta * z / 2) * drho_b;

    //real T = (T_ref + dT);
    //real rho = (rho_ref + drho);
    real T = (T_ref);
    real rho = (rho_ref);

    return rho;
  }

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;

    real delta = g / (Rd * T_ref);
    real rho_s = p_s / (Rd * T_ref);

    real rho_ref = isothermal_zdep(z, rho_s, T_ref, g, thermo);

    real dT_b = dT_max * exp(-pow((x - x_c) / d, 2)) * std::sin(pi * z / Lz);
    real dT = exp(delta * z / 2) * dT_b;

    real drho_b = -rho_s * dT_b / T_ref;
    real drho = exp(-delta * z / 2) * drho_b;

    real T = (T_ref + dT);
    real rho = (rho_ref + drho);

    // real p = rho * Rd * T;
    real p_ref = isothermal_zdep(z, p_s, T_ref, g, thermo);
    real dp = Rd * T_ref * drho + Rd * rho_ref * dT;
    real p = p_ref + dp;

    return thermo.compute_entropic_var(p, T, 0, 0, 0, 0);
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

  static real YAKL_INLINE rhoexact_f(real x, real z, real t,
                                     const ThermoPotential &thermo) {
    return refrho_f(z, thermo) + sum_series(x, z, t, thermo).drho;
  }

  static real YAKL_INLINE entropicdensityexact_f(
      real x, real z, real t, const ThermoPotential &thermo) {

    const auto sol = sum_series(x, z, t, thermo);

    real p_ref = isothermal_zdep(z, p_s, T_ref, g, thermo);
    real rho_ref = refrho_f(z, thermo);

    real rho = rho_ref + sol.drho;
    real p = p_ref + sol.dp;
    real T = T_ref + sol.dT;

    return rho * thermo.compute_entropic_var(p, T, 0, 0, 0, 0);
  }

  static real YAKL_INLINE Texact_f(real x, real z, real t,
                                   const ThermoPotential &thermo) {
    return T_ref + sum_series(x, z, t, thermo).dT;
  }

  static vecext<2> YAKL_INLINE vexact_f(real x, real z, real t,
                                        const ThermoPotential &thermo) {
    vecext<2> v;

    const auto sol = sum_series(x, z, t, thermo);
    v.u = u_0 + sol.du;
    v.w = sol.dw;
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

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {
    real Rd = thermo.cst.Rd;

    real delta = g / (Rd * T_ref);
    real rho_s = p_s / (Rd * T_ref);

    real rho_ref = isothermal_zdep(z, rho_s, T_ref, g, thermo);

    real dT_b = dT_max * exp(-pow((x - x_c) / d, 2)) * std::sin(pi * z / Lz);
    real dT = exp(delta * z / 2) * dT_b;

    real drho_b = -rho_s * dT_b / T_ref;
    real drho = exp(-delta * z / 2) * drho_b;

    //real T = (T_ref + dT);
    real T = (T_ref);
    real rho = (rho_ref + drho);

    // real p = rho * Rd * T;
    real p_ref = isothermal_zdep(z, p_s, T_ref, g, thermo);
    real dp = Rd * T_ref * drho + Rd * rho_ref * dT;
    //real p = p_ref + dp;
    real p = p_ref;

    return thermo.compute_entropic_var(p, T, 0, 0, 0, 0);
  }
>>>>>>> 0da22d4 (debug)

    return sol;
  }

  struct ExactDensityDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "dense";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, ::thermo);
      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return rhoexact_f(x, z, time, thermo);
          },
          field, 0);

      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return entropicdensityexact_f(x, z, time, thermo);
          },
          field, 1);
    }
  };

  struct ExactTemperatureDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "Te";
      topology = pgeom.topology;
      dofs_arr = {0, 0, 1};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, ::thermo);
      dual_geometry.set_00form_values(
          YAKL_LAMBDA(real x, real z) { return Texact_f(x, z, time, thermo); },
          field, 0);
    }
  };

  struct ExactWDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "we";
      topology = pgeom.topology;
      dofs_arr = {0, 1, 1};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, ::thermo);
      primal_geometry.set_01form_values(
          YAKL_LAMBDA(real x, real z) { return vexact_f(x, z, time, thermo); },
          field, 0, LINE_INTEGRAL_TYPE::TANGENT);
    }
  };

  struct BackgroundDensityDiagnostic : public Diagnostic {

    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "densb";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, ::thermo);
      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) { return refrho_f(z, thermo); }, field,
          0);

      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return refentropicdensity_f(z, thermo);
          },
          field, 1);
    }
  };

  struct TemperatureDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "T";
      topology = pgeom.topology;
      dofs_arr = {0, 0, 1};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      const auto &primal_topology = primal_geometry.topology;

      int pis = primal_topology.is;
      int pjs = primal_topology.js;
      int pks = primal_topology.ks;

      const auto &densvar = x.fields_arr[DENSVAR].data;

      YAKL_SCOPE(thermo, ::thermo);
      YAKL_SCOPE(varset, ::varset);
      parallel_for(
          "Compute T diagnostic",
          SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
            real alpha = varset.get_alpha(densvar, k, j, i, pks, pjs, pis, n);
            real entropic_var =
                varset.get_entropic_var(densvar, k, j, i, pks, pjs, pis, n);
            real T = thermo.compute_T(alpha, entropic_var, 1, 0, 0, 0);

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

struct BalancedState {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 10e3;
  static real constexpr Lz = 10e3;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr p_s = 1e5;
  static real constexpr T_ref = 300;
  static real constexpr theta0 = 300;
  static real constexpr N_ref = 0.0001;

  //static real YAKL_INLINE refnsq_f(real z, const ThermoPotential &thermo) {
  //  return N_ref * N_ref;
  //}

  //static real YAKL_INLINE refp_f(real z, const ThermoPotential &thermo) {
  //  return const_stability_p(z, N_ref, g, thermo.cst.pr, theta0, thermo);
  //}

  //static real YAKL_INLINE refT_f(real z, const ThermoPotential &thermo) {
  //  return const_stability_T(z, N_ref, g, theta0, thermo);
  //}

  //static real YAKL_INLINE refrho_f(real z, const ThermoPotential &thermo) {
  //  real p = refp_f(z, thermo);
  //  real T = refT_f(z, thermo);
  //  real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
  //  return 1._fp / alpha;
  //}

  //static real YAKL_INLINE refentropicvar_f(real z,
  //                                             const ThermoPotential &thermo) {
  //  real T_ref = refT_f(z, thermo);
  //  real p_ref = refp_f(z, thermo);
  //  return thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  //}

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
  
  static real YAKL_INLINE refentropicvar_f(real z, const ThermoPotential &thermo) {
    real p = isothermal_zdep(z, p_s, T_ref, g, thermo);
    return thermo.compute_entropic_var(p, T_ref, 0, 0, 0, 0);
  }

  static real YAKL_INLINE refentropicdensity_f(real z,
                                               const ThermoPotential &thermo) {
    real rho_ref = refrho_f(z, thermo);
    return rho_ref * refentropicvar_f(z, thermo);
  }

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    return refrho_f(z, thermo);
    //return isentropic_rho(x, z, theta0, g, thermo);
  }

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {
    return refentropicvar_f(z, thermo);
    //real p = isentropic_p(x, z, theta0, g, thermo);
    //real T = isentropic_T(x, z, theta0, g, thermo);
    //return thermo.compute_entropic_var(p, T, 0, 0, 0, 0);
  }

  static vecext<2> YAKL_INLINE v_f(real x, real y) {
    vecext<2> vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }
  
  struct BackgroundDensityDiagnostic : public Diagnostic {

    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "densb";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, ::thermo);
      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) { return refrho_f(z, thermo); }, field,
          0);

      dual_geometry.set_11form_values(
          YAKL_LAMBDA(real x, real z) {
            return refentropicdensity_f(z, thermo);
          },
          field, 1);
    }
  };

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
    diagnostics.emplace_back(std::make_unique<BackgroundDensityDiagnostic>());
  }
};

void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance) {
  if (name == "doublevortex") {
    testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  } else if (name == "gravitywave") {
    testcase = std::make_unique<EulerTestCase<GravityWave>>();
  } else if (name == "balancedstate") {
    testcase = std::make_unique<EulerTestCase<BalancedState>>();
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
