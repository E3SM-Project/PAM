#pragma once

#include "common.h"
#include "model.h"
#include "profiles.h"
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

struct ModelReferenceState : ReferenceState {
  Profile dens;
  Profile geop;
  Profile q_di;
  Profile q_pi;
  Profile rho_di;
  Profile rho_pi;
  Profile B;
  // should it also be a profile ?
  real4d Blin_coeff;

  void initialize(const Topology &primal_topology,
                  const Topology &dual_topology) override {

    this->dens.initialize(dual_topology, "ref dens", 1, 1, ndensity);
    this->geop.initialize(dual_topology, "ref geop", 1, 1, 1);
    this->rho_pi.initialize(primal_topology, "refrho_pi", 0, 0, 1);
    this->q_pi.initialize(primal_topology, "refq_pi", 0, 0, ndensity);
    this->rho_di.initialize(dual_topology, "refrho_di", 0, 0, 1);
    this->q_di.initialize(dual_topology, "refq_di", 0, 0, ndensity);
    this->B.initialize(dual_topology, "ref B", 1, 1, ndensity_B);

    this->Blin_coeff = real4d("Blin coeff", ndensity_dycore, ndensity_dycore,
                              primal_topology.ni, primal_topology.nens);

    this->is_initialized = true;
  }
};

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

  void compute(real time, const ReferenceState &reference_state,
               const FieldSet<nconstant> &const_vars,
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
          compute_H2bar_ext<ndensity, diff_ord, vert_diff_ord>(
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

  void compute(real time, const ReferenceState &reference_state,
               const FieldSet<nconstant> &const_vars,
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
  real entropicvar_diffusion_coeff;
  real velocity_diffusion_coeff;

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
    entropicvar_diffusion_coeff = params.entropicvar_diffusion_coeff;
    velocity_diffusion_coeff = params.velocity_diffusion_coeff;
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

    const auto &refstate =
        static_cast<ModelReferenceState &>(*this->reference_state);

    YAKL_SCOPE(refdens, refstate.dens.data);
    const auto subtract_refstate_f =
        YAKL_LAMBDA(const real5d &densvar, int d, int k, int j, int i, int n) {
      // have to subtract pks here because refstate doesn't have halos
      return densvar(d, k, j, i, n) - refdens(d, k - pks, n);
    };

    parallel_for(
        "Compute Dens0var",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_H2bar_ext<ndensity, diff_ord, vert_diff_ord>(
              subtract_refstate_f, dens0var, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
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
          compute_H1_ext<1, diff_ord>(Uvar, Vvar, this->primal_geometry,
                                      this->dual_geometry, dis, djs, dks, i, j,
                                      k, n);
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
          compute_H1_vert<1, vert_diff_ord>(UWvar, Wvar, this->primal_geometry,
                                            this->dual_geometry, dis, djs, dks,
                                            i, j, k + 1, n);
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

  void compute_edge_reconstructions_uniform(
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
          compute_twisted_vert_edge_recon_uniform<
              ndensity, dual_vert_reconstruction_type,
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
          compute_straight_xz_vert_edge_recon_uniform<
              1, vert_reconstruction_type, vert_reconstruction_order>(
              qxzvertedgereconvar, qxz0var, pis, pjs, pks, i, j, k, n,
              primal_vert_wenoRecon, primal_vert_to_gll, primal_vert_wenoIdl,
              primal_vert_wenoSigma);
          compute_straight_xz_edge_recon<1, coriolis_reconstruction_type,
                                         coriolis_reconstruction_order>(
              coriolisxzedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
              coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
              coriolis_wenoSigma);
          compute_straight_xz_vert_edge_recon_uniform<
              1, coriolis_vert_reconstruction_type,
              coriolis_vert_reconstruction_order>(
              coriolisxzvertedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
              coriolis_vert_wenoRecon, coriolis_vert_to_gll,
              coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
        });
  }

  void compute_edge_reconstructions_variable(
      real5d densedgereconvar, real5d densvertedgereconvar,
      real5d qxzedgereconvar, real5d qxzvertedgereconvar,
      real5d coriolisxzedgereconvar, real5d coriolisxzvertedgereconvar,
      const real5d dens0var, const real5d qxz0var, const real5d fxz0var) {

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

    parallel_for(
        "ComputeDensEdgeRecon",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
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

          compute_twisted_edge_recon<ndensity, dual_reconstruction_type,
                                     dual_reconstruction_order>(
              densedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);

          compute_twisted_vert_edge_recon_variable<
              ndensity, dual_vert_reconstruction_type,
              dual_vert_reconstruction_order>(
              densvertedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_vert_coefs_to_gll, dual_vert_sten_to_gll,
              dual_vert_sten_to_coefs, dual_vert_weno_recon_lower,
              dual_vert_wenoIdl, dual_vert_wenoSigma);
        });

    parallel_for(
        "ComputeQEdgeRecon",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
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

          compute_straight_xz_edge_recon<1, reconstruction_type,
                                         reconstruction_order>(
              qxzedgereconvar, qxz0var, pis, pjs, pks, i, j, k, n,
              primal_wenoRecon, primal_to_gll, primal_wenoIdl,
              primal_wenoSigma);
          compute_straight_xz_vert_edge_recon_variable<
              1, vert_reconstruction_type, vert_reconstruction_order>(
              qxzvertedgereconvar, qxz0var, pis, pjs, pks, i, j, k, n,
              primal_vert_coefs_to_gll, primal_vert_sten_to_gll,
              primal_vert_sten_to_coefs, primal_vert_weno_recon_lower,
              primal_vert_wenoIdl, primal_vert_wenoSigma);
          compute_straight_xz_edge_recon<1, coriolis_reconstruction_type,
                                         coriolis_reconstruction_order>(
              coriolisxzedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
              coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
              coriolis_wenoSigma);
          compute_straight_xz_vert_edge_recon_variable<
              1, coriolis_vert_reconstruction_type,
              coriolis_vert_reconstruction_order>(
              coriolisxzvertedgereconvar, fxz0var, pis, pjs, pks, i, j, k, n,
              coriolis_vert_coefs_to_gll, coriolis_vert_sten_to_gll,
              coriolis_vert_sten_to_coefs, coriolis_vert_weno_recon_lower,
              coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
        });
  }

  void compute_recons(
      real5d densreconvar, real5d densvertreconvar, real5d qxzreconvar,
      real5d qxzvertreconvar, real5d coriolisxzreconvar,
      real5d coriolisxzvertreconvar, const real5d densedgereconvar,
      const real5d densvertedgereconvar, const real5d qxzedgereconvar,
      const real5d qxzvertedgereconvar, const real5d coriolisxzedgereconvar,
      const real5d coriolisxzvertedgereconvar, const real5d densvar,
      const real5d Vvar, const real5d Wvar) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(varset, ::varset);

    const auto total_density_f =
        YAKL_LAMBDA(const real5d &densvar, int d, int k, int j, int i, int n) {
      return varset.get_total_density(densvar, k, j, i, 0, 0, 0, n);
    };

    const auto &refstate =
        static_cast<ModelReferenceState &>(*this->reference_state);

    parallel_for(
        "ComputeDensRECON",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_recon<ndensity, dual_reconstruction_type>(
              densreconvar, densedgereconvar, this->primal_geometry,
              this->dual_geometry, Vvar, dis, djs, dks, i, j, k, n);

          SArray<real, 1, 1> dens0_ik;
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_ik, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
          SArray<real, 1, 1> dens0_im1;
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_im1, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i - 1, j, k, n);

          real he = 0.5_fp * (dens0_ik(0) + dens0_im1(0));
          // scale twisted recons and add reference state
          for (int d = 0; d < ndims; d++) {
            for (int l = 0; l < ndensity; l++) {
              densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) +=
                  refstate.rho_pi.data(0, k, n) * refstate.q_pi.data(l, k, n);
              densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) /=
                  he;
            }
          }
        });

    parallel_for(
        "ComputeDensVertRECON",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_vert_recon<ndensity, dual_vert_reconstruction_type>(
              densvertreconvar, densvertedgereconvar, this->primal_geometry,
              this->dual_geometry, Wvar, dis, djs, dks, i, j, k + 1, n);

          SArray<real, 1, 1> dens0_kp1;
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_kp1, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k + 1, n);
          SArray<real, 1, 1> dens0_ik;
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_ik, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);

          real hew = 0.5_fp * (dens0_kp1(0) + dens0_ik(0));
          // scale twisted recons and add reference state
          for (int l = 0; l < ndensity; l++) {
            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) +=
                refstate.rho_di.data(0, k + 1, n) *
                refstate.q_di.data(l, k + 1, n);
            densvertreconvar(l, k + dks + 1, j + djs, i + dis, n) /= hew;
          }
        });

    parallel_for(
        "ComputeQRECON",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_straight_xz_recon<1, reconstruction_type>(
              qxzreconvar, qxzedgereconvar, this->primal_geometry,
              this->dual_geometry, Vvar, pis, pjs, pks, i, j, k, n);
          compute_straight_xz_recon<1, coriolis_reconstruction_type>(
              coriolisxzreconvar, coriolisxzedgereconvar, this->primal_geometry,
              this->dual_geometry, Vvar, pis, pjs, pks, i, j, k, n);
        });
    parallel_for(
        "ComputeQVERTRECON",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_straight_xz_vert_recon<1, vert_reconstruction_type>(
              qxzvertreconvar, qxzvertedgereconvar, this->primal_geometry,
              this->dual_geometry, Wvar, pis, pjs, pks, i, j, k, n);
          compute_straight_xz_vert_recon<1, coriolis_vert_reconstruction_type>(
              coriolisxzvertreconvar, coriolisxzvertedgereconvar,
              this->primal_geometry, this->dual_geometry, Wvar, pis, pjs, pks,
              i, j, k, n);
        });
  }

  void add_entropicvar_diffusion(real entropicvar_coeff, real5d denstendvar,
                                 const real5d densvar, const real5d dens0var,
                                 const real5d Kvar, const real5d qxzfluxvar,
                                 const real5d qxzvertfluxvar, const real5d Fvar,
                                 const real5d FWvar,
                                 FieldSet<nauxiliary> &auxiliary_vars) {

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(varset, ::varset);
    parallel_for(
        "Entropicvar diffusion 1",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          dens0var(0, k + pks, j + pjs, i + pis, n) =
              varset.get_entropic_var(densvar, k, j, i, pks, pjs, pis, n);
        });
    auxiliary_vars.exchange({DENS0VAR});

    parallel_for(
        "Entropicvar diffusion 2",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D0<1>(qxzvertfluxvar, dens0var, pis, pjs, pks, i, j, k, n);
        });
    auxiliary_vars.exchange({QXZVERTFLUXVAR});

    parallel_for(
        "Entropicvar diffusion 3",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D0_vert<1>(qxzfluxvar, dens0var, pis, pjs, pks, i, j, k, n);
        });
    auxiliary_vars.exchange({QXZFLUXVAR});

    parallel_for(
        "Entropicvar diffusion 4",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H1_ext<1, diffusion_diff_ord>(Fvar, qxzvertfluxvar,
                                                primal_geometry, dual_geometry,
                                                dis, djs, dks, i, j, k, n);
        });
    auxiliary_vars.exchange({FVAR});

    parallel_for(
        "Entropicvar diffusion 5",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H1_vert<1, diffusion_diff_ord>(FWvar, qxzfluxvar,
                                                 primal_geometry, dual_geometry,
                                                 dis, djs, dks, i, j, k + 1, 0);
        });
    auxiliary_vars.exchange({FWVAR});
    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);

    parallel_for(
        "Entropicvar diffusion 6",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D1bar<1>(Kvar, Fvar, dis, djs, dks, i, j, k, n);
          compute_D1bar_vert<1, ADD_MODE::ADD>(Kvar, FWvar, dis, djs, dks, i, j,
                                               k, n);
        });
    auxiliary_vars.exchange({KVAR});

    parallel_for(
        "Entropicvar diffusion 7",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, 1> sdiff;
          compute_H2bar_ext<1, diffusion_diff_ord, vert_diffusion_diff_ord>(
              sdiff, Kvar, primal_geometry, dual_geometry, pis, pjs, pks, i, j,
              k, n);
          sdiff(0) *= entropicvar_coeff;

          real rho = densvar(0, k + pks, j + pjs, i + pis, n);
          denstendvar(1, k + pks, j + pjs, i + pis, n) += -rho * sdiff(0);
        });
  }

  void add_velocity_diffusion(real velocity_coeff, real5d Vtendvar,
                              real5d Wtendvar, const real5d Vvar,
                              const real5d Wvar, const real5d qxzedgereconvar,
                              const real5d qxz0var, const real5d dens0var,
                              const real5d Kvar, const real5d Fvar,
                              const real5d FWvar,
                              FieldSet<nauxiliary> &auxiliary_vars) {

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    // *d*d

    parallel_for(
        "Velocity diffusion 1",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D1_ext<1>(qxzedgereconvar, Vvar, Wvar, pis, pjs, pks, i, j, k,
                            n);
        });
    auxiliary_vars.exchange({QXZEDGERECONVAR});

    parallel_for(
        "Velocity diffusion 2",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H2_ext<1, diffusion_diff_ord, vert_diffusion_diff_ord>(
              qxz0var, qxzedgereconvar, primal_geometry, dual_geometry, dis,
              djs, dks, i, j, k + 1, n);
        });
    auxiliary_vars.exchange({QXZ0VAR});
    auxiliary_vars.fields_arr[QXZ0VAR].set_bnd(0.0);

    parallel_for(
        "Velocity diffusion 3",
        SimpleBounds<4>(dual_topology.ni, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D0bar_ext<1>(FWvar, qxz0var, dis, djs, dks, i, j, k, n);
        });
    auxiliary_vars.exchange({FWVAR});

    parallel_for(
        "Velocity diffusion 4",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D0bar_vert<1>(Fvar, qxz0var, dis, djs, dks, i, j, k, n);
        });
    auxiliary_vars.exchange({FVAR});

    parallel_for(
        "Velocity diffusion 5",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> vdiff;
          compute_H1bar_ext<1, diffusion_diff_ord>(vdiff, Fvar, primal_geometry,
                                                   dual_geometry, pis, pjs, pks,
                                                   i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            Vtendvar(d, k + pks, j + pjs, i + pis, n) -=
                velocity_coeff * vdiff(d);
          }
        });

    parallel_for(
        "Velocity diffusion 6",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> wdiff;
          compute_H1bar_vert<1, vert_diffusion_diff_ord>(
              wdiff, FWvar, primal_geometry, dual_geometry, pis, pjs, pks, i, j,
              k, n);
          Wtendvar(0, k + pks, j + pjs, i + pis, n) -=
              velocity_coeff * wdiff(0);
        });

    // d*d*

    parallel_for(
        "Velocity diffusion 7",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H1_ext<1, diffusion_diff_ord>(Fvar, Vvar, primal_geometry,
                                                dual_geometry, dis, djs, dks, i,
                                                j, k, n);
        });
    auxiliary_vars.exchange({FVAR});

    parallel_for(
        "Velocity diffusion 8",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H1_vert<1, vert_diffusion_diff_ord>(
              FWvar, Wvar, primal_geometry, dual_geometry, dis, djs, dks, i, j,
              k + 1, n);
        });
    auxiliary_vars.exchange({FWVAR});
    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);

    parallel_for(
        "Velocity diffusion 9",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D1bar<1>(Kvar, Fvar, dis, djs, dks, i, j, k, n);
          compute_D1bar_vert<1, ADD_MODE::ADD>(Kvar, FWvar, dis, djs, dks, i, j,
                                               k, n);
        });
    auxiliary_vars.exchange({KVAR});

    parallel_for(
        "Velocity diffusion 10",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H2bar_ext<1, diffusion_diff_ord, vert_diffusion_diff_ord>(
              dens0var, Kvar, primal_geometry, dual_geometry, pis, pjs, pks, i,
              j, k, n);
        });
    auxiliary_vars.exchange({DENS0VAR});

    parallel_for(
        "Velocity diffusion 11",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> vdiff;
          compute_D0<1>(vdiff, dens0var, pis, pjs, pks, i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            Vtendvar(d, pks + k, pjs + j, pis + i, n) -=
                velocity_coeff * vdiff(d);
          }
        });

    parallel_for(
        "Velocity diffusion 12",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real wdiff = compute_D0_vert<1>(dens0var, pis, pjs, pks, i, j, k, n);
          Wtendvar(0, pks + k, pjs + j, pis + i, n) -= velocity_coeff * wdiff;
        });
  }

  void
  compute_tendencies(real5d denstendvar, real5d Vtendvar, real5d Wtendvar,
                     const real5d densreconvar, const real5d densvertreconvar,
                     const real5d qxzreconvar, const real5d qxzvertreconvar,
                     const real5d coriolisxzreconvar,
                     const real5d coriolisxzvertreconvar, const real5d Bvar,
                     const real5d Fvar, const real5d FWvar) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    const auto &refstate =
        static_cast<ModelReferenceState &>(*this->reference_state);

    parallel_for(
        "Compute Wtend",
        SimpleBounds<4>(primal_topology.nl - 2, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD0_vert<ndensity_B>(Wtendvar, densvertreconvar, Bvar, pis,
                                       pjs, pks, i, j, k + 1, n);
          if (force_refstate_hydrostatic_balance) {
            compute_wD0_vert<ndensity_B, ADD_MODE::ADD>(
                Wtendvar, refstate.q_di.data, refstate.B.data, pis, pjs, pks, i,
                j, k + 1, n);
          }
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
          compute_wD0_vert<ndensity_B>(Wtendvar, densvertreconvar, Bvar, pis,
                                       pjs, pks, i, j, 0, n);
          compute_wD0_vert<ndensity_B>(Wtendvar, densvertreconvar, Bvar, pis,
                                       pjs, pks, i, j, primal_topology.nl - 1,
                                       n);
          if (force_refstate_hydrostatic_balance) {
            compute_wD0_vert<ndensity_B, ADD_MODE::ADD>(
                Wtendvar, refstate.q_di.data, refstate.B.data, pis, pjs, pks, i,
                j, 0, n);
            compute_wD0_vert<ndensity_B, ADD_MODE::ADD>(
                Wtendvar, refstate.q_di.data, refstate.B.data, pis, pjs, pks, i,
                j, primal_topology.nl - 1, n);
          }
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
          compute_wD0<ndensity_B>(Vtendvar, densreconvar, Bvar, pis, pjs, pks,
                                  i, j, k + 1, n);
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
          compute_wD0<ndensity_B>(Vtendvar, densreconvar, Bvar, pis, pjs, pks,
                                  i, j, 0, n);
          compute_wD0<ndensity_B>(Vtendvar, densreconvar, Bvar, pis, pjs, pks,
                                  i, j, primal_topology.ni - 1, n);
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
          compute_wD1bar<ndensity>(denstendvar, densreconvar, Fvar, dis, djs,
                                   dks, i, j, k, n);
          compute_wD1bar_vert<ndensity, ADD_MODE::ADD>(
              denstendvar, densvertreconvar, FWvar, dis, djs, dks, i, j, k, n);
        });
  }

  void compute_functional_derivatives(
      ADD_MODE addmode, real fac, real dt, FieldSet<nconstant> &const_vars,
      FieldSet<nprognostic> &x, FieldSet<nauxiliary> &auxiliary_vars) override {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    const auto &densvar = x.fields_arr[DENSVAR].data;
    const auto &Vvar = x.fields_arr[VVAR].data;
    const auto &Wvar = x.fields_arr[WVAR].data;

    const auto &Fvar = auxiliary_vars.fields_arr[FVAR].data;
    const auto &FWvar = auxiliary_vars.fields_arr[FWVAR].data;
    const auto &Kvar = auxiliary_vars.fields_arr[KVAR].data;
    const auto &Bvar = auxiliary_vars.fields_arr[BVAR].data;
    const auto &HSvar = const_vars.fields_arr[HSVAR].data;

    YAKL_SCOPE(Hk, ::Hk);
    YAKL_SCOPE(Hs, ::Hs);
    YAKL_SCOPE(varset, ::varset);

    const auto total_density_f =
        YAKL_LAMBDA(const real5d &densvar, int d, int k, int j, int i, int n) {
      return varset.get_total_density(densvar, k, j, i, 0, 0, 0, n);
    };

    parallel_for(
        "Functional derivatives",
        SimpleBounds<4>(dual_topology.ni, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, 1> dens0_ik, dens0_im1, dens0_km1;
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_ik, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_im1, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i - 1, j, k, n);
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_km1, densvar, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k - 1, n);

          SArray<real, 1, ndims> u_ik;
          SArray<real, 1, 1> uw_ik;
          SArray<real, 1, ndims> u_ip1;
          SArray<real, 1, 1> uw_kp1;

          compute_H1_ext<1, diff_ord>(u_ik, Vvar, this->primal_geometry,
                                      this->dual_geometry, pis, pjs, pks, i, j,
                                      k, n);
          compute_H1_ext<1, diff_ord>(u_ip1, Vvar, this->primal_geometry,
                                      this->dual_geometry, pis, pjs, pks, i + 1,
                                      j, k, n);

          if (k == 0 || k == (dual_topology.ni - 1)) {
            uw_ik(0) = 0;
          } else {
            compute_H1_vert<1, vert_diff_ord>(
                uw_ik, Wvar, this->primal_geometry, this->dual_geometry, pis,
                pjs, pks, i, j, k, n);
          }

          if (k >= (dual_topology.ni - 2)) {
            uw_kp1(0) = 0;
          } else {
            compute_H1_vert<1, vert_diff_ord>(
                uw_kp1, Wvar, this->primal_geometry, this->dual_geometry, pis,
                pjs, pks, i, j, k + 1, n);
          }

          real K2 =
              0.5_fp * (Vvar(0, k + pks, j + pjs, i + pis, n) * u_ik(0) +
                        Vvar(0, k + pks, j + pjs, i + 1 + pis, n) * u_ip1(0));

          if (k < dual_topology.nl) {
            K2 +=
                0.5_fp * (Wvar(0, k - 1 + pks, j + pjs, i + pis, n) * uw_ik(0) +
                          Wvar(0, k + pks, j + pjs, i + pis, n) * uw_kp1(0));
          }

          Kvar(0, k + pks, j + pjs, i + pis, n) = 0.5_fp * K2;

          if (addmode == ADD_MODE::ADD) {
            Fvar(0, pks + k, pjs + j, pis + i, n) +=
                fac * 0.5_fp * (dens0_ik(0) + dens0_im1(0)) * u_ik(0);
            FWvar(0, pks + k, pjs + j, pis + i, n) +=
                fac * 0.5_fp * (dens0_ik(0) + dens0_km1(0)) * uw_ik(0);
          } else if (addmode == ADD_MODE::REPLACE) {
            Fvar(0, pks + k, pjs + j, pis + i, n) =
                fac * 0.5_fp * (dens0_ik(0) + dens0_im1(0)) * u_ik(0);
            FWvar(0, pks + k, pjs + j, pis + i, n) =
                fac * 0.5_fp * (dens0_ik(0) + dens0_km1(0)) * uw_ik(0);
          }

          if (addmode == ADD_MODE::ADD) {
            Hs.compute_dHsdx<ADD_MODE::ADD>(Bvar, densvar, HSvar, pis, pjs, pks,
                                            i, j, k, n, fac);
          } else if (addmode == ADD_MODE::REPLACE) {
            Hs.compute_dHsdx<ADD_MODE::REPLACE>(Bvar, densvar, HSvar, pis, pjs,
                                                pks, i, j, k, n, fac);
          }
        });
    auxiliary_vars.exchange({KVAR});

    parallel_for(
        "Add K to B",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j, k,
                                            n, fac);
        });

    auxiliary_vars.exchange({BVAR, FVAR, FWVAR});
  }

  void compute_functional_derivatives_two_point(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x1,
      FieldSet<nprognostic> &x2,
      FieldSet<nauxiliary> &auxiliary_vars) override {

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

    YAKL_SCOPE(Hk, ::Hk);
    YAKL_SCOPE(Hs, ::Hs);
    YAKL_SCOPE(varset, ::varset);

    const auto total_density_f =
        YAKL_LAMBDA(const real5d &densvar, int d, int k, int j, int i, int n) {
      return varset.get_total_density(densvar, k, j, i, 0, 0, 0, n);
    };

    parallel_for(
        "Functional derivatives",
        SimpleBounds<4>(dual_topology.ni, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, 1> dens0_ik_1, dens0_im1_1, dens0_km1_1;
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_ik_1, densvar1, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_im1_1, densvar1, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i - 1, j, k, n);
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_km1_1, densvar1, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k - 1, n);

          SArray<real, 1, 1> dens0_ik_2, dens0_im1_2, dens0_km1_2;
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_ik_2, densvar2, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_im1_2, densvar2, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i - 1, j, k, n);
          compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
              total_density_f, dens0_km1_2, densvar2, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k - 1, n);

          SArray<real, 1, ndims> u_ik_1;
          SArray<real, 1, 1> uw_ik_1;
          SArray<real, 1, ndims> u_ip1_1;
          SArray<real, 1, 1> uw_kp1_1;

          compute_H1_ext<1, diff_ord>(u_ik_1, Vvar1, this->primal_geometry,
                                      this->dual_geometry, pis, pjs, pks, i, j,
                                      k, n);
          compute_H1_ext<1, diff_ord>(u_ip1_1, Vvar1, this->primal_geometry,
                                      this->dual_geometry, pis, pjs, pks, i + 1,
                                      j, k, n);

          if (k == 0 || k == (dual_topology.ni - 1)) {
            uw_ik_1(0) = 0;
          } else {
            compute_H1_vert<1, vert_diff_ord>(
                uw_ik_1, Wvar1, this->primal_geometry, this->dual_geometry, pis,
                pjs, pks, i, j, k, n);
          }

          if (k >= (dual_topology.ni - 2)) {
            uw_kp1_1(0) = 0;
          } else {
            compute_H1_vert<1, vert_diff_ord>(
                uw_kp1_1, Wvar1, this->primal_geometry, this->dual_geometry,
                pis, pjs, pks, i, j, k + 1, n);
          }

          SArray<real, 1, ndims> u_ik_2;
          SArray<real, 1, 1> uw_ik_2;
          SArray<real, 1, ndims> u_ip1_2;
          SArray<real, 1, 1> uw_kp1_2;

          compute_H1_ext<1, diff_ord>(u_ik_2, Vvar2, this->primal_geometry,
                                      this->dual_geometry, pis, pjs, pks, i, j,
                                      k, n);
          compute_H1_ext<1, diff_ord>(u_ip1_2, Vvar2, this->primal_geometry,
                                      this->dual_geometry, pis, pjs, pks, i + 1,
                                      j, k, n);

          if (k == 0 || k == (dual_topology.ni - 1)) {
            uw_ik_2(0) = 0;
          } else {
            compute_H1_vert<1, vert_diff_ord>(
                uw_ik_2, Wvar2, this->primal_geometry, this->dual_geometry, pis,
                pjs, pks, i, j, k, n);
          }

          if (k >= (dual_topology.ni - 2)) {
            uw_kp1_2(0) = 0;
          } else {
            compute_H1_vert<1, vert_diff_ord>(
                uw_kp1_2, Wvar2, this->primal_geometry, this->dual_geometry,
                pis, pjs, pks, i, j, k + 1, n);
          }

          real K2_1 = 0.5_fp *
                      (Vvar1(0, k + pks, j + pjs, i + pis, n) * u_ik_1(0) +
                       Vvar1(0, k + pks, j + pjs, i + 1 + pis, n) * u_ip1_1(0));

          real K2_2 = 0.5_fp *
                      (Vvar2(0, k + pks, j + pjs, i + pis, n) * u_ik_2(0) +
                       Vvar2(0, k + pks, j + pjs, i + 1 + pis, n) * u_ip1_2(0));

          if (k < dual_topology.nl) {
            K2_1 += 0.5_fp *
                    (Wvar1(0, k - 1 + pks, j + pjs, i + pis, n) * uw_ik_1(0) +
                     Wvar1(0, k + pks, j + pjs, i + pis, n) * uw_kp1_1(0));

            K2_2 += 0.5_fp *
                    (Wvar2(0, k - 1 + pks, j + pjs, i + pis, n) * uw_ik_2(0) +
                     Wvar2(0, k + pks, j + pjs, i + pis, n) * uw_kp1_2(0));
          }

          Kvar(0, k + pks, j + pjs, i + pis, n) = 0.25_fp * (K2_1 + K2_2);

          Fvar(0, pks + k, pjs + j, pis + i, n) =
              0.125_fp *
              (dens0_ik_1(0) + dens0_im1_1(0) + dens0_ik_2(0) +
               dens0_im1_2(0)) *
              (u_ik_1(0) + u_ik_2(0));
          FWvar(0, pks + k, pjs + j, pis + i, n) =
              0.125_fp *
              (dens0_ik_1(0) + dens0_km1_1(0) + dens0_ik_2(0) +
               dens0_km1_2(0)) *
              (uw_ik_1(0) + uw_ik_2(0));

          Hs.compute_dHsdx_two_point(Bvar, densvar1, densvar2, HSvar, pis, pjs,
                                     pks, i, j, k, n);
        });
    auxiliary_vars.exchange({KVAR});

    parallel_for(
        "Add K to B",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_dKddens<ADD_MODE::ADD>(Bvar, Kvar, pis, pjs, pks, i, j, k,
                                            n);
        });

    auxiliary_vars.exchange({BVAR, FVAR, FWVAR});
  }

  void apply_symplectic(real dt, FieldSet<nconstant> &const_vars,
                        FieldSet<nprognostic> &x,
                        FieldSet<nauxiliary> &auxiliary_vars,
                        FieldSet<nprognostic> &xtend) override {

    const auto &dual_topology = dual_geometry.topology;

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);

    auxiliary_vars.exchange({DENS0VAR});

    compute_q0f0(auxiliary_vars.fields_arr[QXZ0VAR].data,
                 auxiliary_vars.fields_arr[FXZ0VAR].data,
                 x.fields_arr[VVAR].data, x.fields_arr[WVAR].data,
                 x.fields_arr[DENSVAR].data,
                 const_vars.fields_arr[CORIOLISXZVAR].data);

    auxiliary_vars.fields_arr[QXZ0VAR].set_bnd(0.0);
    auxiliary_vars.fields_arr[FXZ0VAR].set_bnd(0.0);
    auxiliary_vars.exchange({QXZ0VAR, FXZ0VAR});

    // Compute densrecon, densvertrecon, qrecon and frecon
    if (dual_geometry.uniform_vertical) {
      compute_edge_reconstructions_uniform(
          auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QXZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[QXZ0VAR].data,
          auxiliary_vars.fields_arr[FXZ0VAR].data);
    } else {
      compute_edge_reconstructions_variable(
          auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENSVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QXZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISXZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISXZVERTEDGERECONVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[QXZ0VAR].data,
          auxiliary_vars.fields_arr[FXZ0VAR].data);
    }

    auxiliary_vars.exchange({DENSEDGERECONVAR, DENSVERTEDGERECONVAR,
                             QXZEDGERECONVAR, QXZVERTEDGERECONVAR,
                             CORIOLISXZEDGERECONVAR,
                             CORIOLISXZVERTEDGERECONVAR});

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
                   x.fields_arr[DENSVAR].data, x.fields_arr[VVAR].data,
                   x.fields_arr[WVAR].data);

    auxiliary_vars.exchange({DENSRECONVAR, DENSVERTRECONVAR, QXZRECONVAR,
                             QXZVERTRECONVAR, CORIOLISXZRECONVAR,
                             CORIOLISXZVERTRECONVAR});

    // Compute fct
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(dens_pos, varset.dens_pos);
    parallel_for(
        "ComputeEdgeFlux",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_edgefluxes<ndensity>(
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[DENSRECONVAR].data,
              auxiliary_vars.fields_arr[FVAR].data, dens_pos, dis, djs, dks, i,
              j, k, n);

          if (k < dual_topology.ni - 2) {
            compute_vertedgefluxes<ndensity>(
                auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data,
                auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                auxiliary_vars.fields_arr[FWVAR].data, dens_pos, dis, djs, dks,
                i, j, k, n);
          }
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
              auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data, dt, dens_pos,
              dis, djs, dks, i, j, k, n);
        });

    auxiliary_vars.exchange({MFVAR});

    parallel_for(
        "Apply Phi",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          apply_Phi<ndensity>(auxiliary_vars.fields_arr[DENSRECONVAR].data,
                              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
                              auxiliary_vars.fields_arr[MFVAR].data,
                              x.fields_arr[DENSVAR].data, dens_pos, dis, djs,
                              dks, i, j, k, n);
          if (k < dual_topology.ni - 2) {
            apply_Phivert<ndensity>(
                auxiliary_vars.fields_arr[DENSVERTRECONVAR].data,
                auxiliary_vars.fields_arr[VERTEDGEFLUXVAR].data,
                auxiliary_vars.fields_arr[MFVAR].data,
                x.fields_arr[DENSVAR].data, dens_pos, dis, djs, dks, i, j,
                k + 1, n);
          }
        });

    auxiliary_vars.exchange({DENSRECONVAR, DENSVERTRECONVAR});

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
                       auxiliary_vars.fields_arr[FWVAR].data);

    if (entropicvar_diffusion_coeff > 0) {
      add_entropicvar_diffusion(
          entropicvar_diffusion_coeff, xtend.fields_arr[DENSVAR].data,
          x.fields_arr[DENSVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[KVAR].data,
          auxiliary_vars.fields_arr[QXZFLUXVAR].data,
          auxiliary_vars.fields_arr[QXZVERTFLUXVAR].data,
          auxiliary_vars.fields_arr[FVAR].data,
          auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars);
    }
    if (velocity_diffusion_coeff > 0) {
      add_velocity_diffusion(
          velocity_diffusion_coeff, xtend.fields_arr[VVAR].data,
          xtend.fields_arr[WVAR].data, x.fields_arr[VVAR].data,
          x.fields_arr[WVAR].data,
          auxiliary_vars.fields_arr[QXZEDGERECONVAR].data,
          auxiliary_vars.fields_arr[QXZ0VAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[KVAR].data,
          auxiliary_vars.fields_arr[FVAR].data,
          auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars);
    }
  }

  void remove_negative_densities(FieldSet<nprognostic> &x) override {
    const auto &dual_topology = dual_geometry.topology;

    const auto densvar = x.fields_arr[DENSVAR].data;
    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(dens_pos, varset.dens_pos);
    parallel_for(
        "Remove negatives",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          for (int d = 0; d < ndensity; ++d) {
            if (dens_pos(d)) {
              densvar(d, k + dks, j + djs, i + dis, n) =
                  std::max(0._fp, densvar(d, k + dks, j + djs, i + dis, n));
            }
          }
        });
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
    yakl::memset(v_transform, 0);
    yakl::memset(w_transform, 0);

    complex_vrhs = complex4d("complex vrhs", pni, ny, nx, nens);
    complex_wrhs = complex4d("complex wrhs", pnl, ny, nx, nens);

    fftv_x.init(v_transform, 2, nx);
    fftw_x.init(w_transform, 2, nx);
    // fftv_y.init(v_transform, 1, ny);
    // fftw_y.init(w_transform, 1, ny);

    complex_vcoeff =
        complex5d("complex vcoeff", 1 + ndensity_dycore, pni, ny, nx, nens);

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

    const auto &rho_pi = refstate.rho_pi.data;
    const auto &q_pi = refstate.q_pi.data;
    const auto &rho_di = refstate.rho_di.data;
    const auto &q_di = refstate.q_di.data;

    parallel_for(
        "compute vcoeff",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> fD0Dbar;
          fourier_cwD0D1bar(fD0Dbar, 1, i, j, k, dual_topology.n_cells_x,
                            dual_topology.n_cells_y, dual_topology.ni);

          real fH2bar = fourier_H2bar_ext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          SArray<real, 1, ndims> fH1;
          fourier_H1_ext<diff_ord>(fH1, primal_geometry, dual_geometry, pis,
                                   pjs, pks, i, j, k, 0, n_cells_x, n_cells_y,
                                   dual_topology.ni);
          SArray<complex, 1, ndims> fD0;
          fourier_cwD0(fD0, 1, i, j, k, dual_topology.n_cells_x,
                       dual_topology.n_cells_y, dual_topology.ni);

          real he = rho_pi(0, k, n);

          real c1 = 1;
          for (int d1 = 0; d1 < ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < ndensity_dycore; ++d2) {
              c1 -= dtf2 * fH2bar * fH1(0) * fD0Dbar(0) * he * q_pi(d1, k, n) *
                    q_pi(d2, k, n) * refstate.Blin_coeff(d1, d2, k, n);
            }
          }

          complex_vcoeff(0, k, j, i, n) = 1 / c1;
          for (int d1 = 0; d1 < ndensity_dycore; ++d1) {
            complex cd1 = 0;
            for (int d2 = 0; d2 < ndensity_dycore; ++d2) {
              cd1 += fD0(0) * dtf2 * fH2bar * q_pi(d2, k, n) *
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
          real fH2bar_k = fourier_H2bar_ext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fH2bar_kp1 = fourier_H2bar_ext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          real gamma_fac_kp2 = rho_di(0, k + 2, n) *
                               H1_vert_coeff(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k + 2, n);
          real gamma_fac_kp1 = rho_di(0, k + 1, n) *
                               H1_vert_coeff(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k + 1, n);
          real gamma_fac_k =
              rho_di(0, k, n) * H1_vert_coeff(primal_geometry, dual_geometry,
                                              pis, pjs, pks, i, j, k, n);

          tri_u(k, j, i, n) = 0;
          tri_d(k, j, i, n) = 1;
          tri_l(k, j, i, n) = 0;

          for (int d1 = 0; d1 < ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < ndensity_dycore; ++d2) {
              real alpha_kp1 = q_di(d1, k + 1, n);

              real beta_kp1 =
                  fH2bar_kp1 * refstate.Blin_coeff(d1, d2, k + 1, n);
              real beta_k = fH2bar_k * refstate.Blin_coeff(d1, d2, k, n);

              real gamma_kp2 = gamma_fac_kp2 * q_di(d2, k + 2, n);
              real gamma_kp1 = gamma_fac_kp1 * q_di(d2, k + 1, n);
              real gamma_k = gamma_fac_k * q_di(d2, k, n);

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
          real fH2bar_k = fourier_H2bar_ext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fH2bar_kp1 = fourier_H2bar_ext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          real gamma_fac_kp2 = rho_di(0, k + 2, n) *
                               H1_vert_coeff(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k + 2, n);
          real gamma_fac_kp1 = rho_di(0, k + 1, n) *
                               H1_vert_coeff(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k + 1, n);
          real gamma_fac_k =
              rho_di(0, k, n) * H1_vert_coeff(primal_geometry, dual_geometry,
                                              pis, pjs, pks, i, j, k, n);

          SArray<real, 1, ndims> fH1_kp1_a;
          SArray<real, 1, ndims> fH1_k_a;
          fourier_H1_ext<diff_ord>(fH1_kp1_a, primal_geometry, dual_geometry,
                                   pis, pjs, pks, i, j, k + 1, 0, n_cells_x,
                                   n_cells_y, dual_topology.ni);
          fourier_H1_ext<diff_ord>(fH1_k_a, primal_geometry, dual_geometry, pis,
                                   pjs, pks, i, j, k, 0, n_cells_x, n_cells_y,
                                   dual_topology.ni);
          real fH1h_kp1 = fH1_kp1_a(0);
          real fH1h_k = fH1_k_a(0);

          complex fD1bar_kp1 = fourier_D1bar(1, i, j, k + 1, n_cells_x,
                                             n_cells_y, dual_topology.ni);
          complex fD1bar_k =
              fourier_D1bar(1, i, j, k, n_cells_x, n_cells_y, dual_topology.ni);

          real he_kp1 = rho_pi(0, k + 1, n);
          real he_k = rho_pi(0, k, n);

          for (int d1 = 0; d1 < ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < ndensity_dycore; ++d2) {
              for (int d3 = 0; d3 < ndensity_dycore; ++d3) {

                real alpha_kp1 = dtf2 * q_di(d1, k + 1, n);
                complex beta_kp1 =
                    fH2bar_kp1 * refstate.Blin_coeff(d1, d2, k + 1, n) *
                    q_pi(d2, k + 1, n) * fD1bar_kp1 * he_kp1 * fH1h_kp1;
                complex beta_k = fH2bar_k * refstate.Blin_coeff(d1, d2, k, n) *
                                 q_pi(d2, k, n) * fD1bar_k * he_k * fH1h_k;

                real gamma_kp2 = gamma_fac_kp2 * q_di(d3, k + 2, n);
                real gamma_kp1 = gamma_fac_kp1 * q_di(d3, k + 1, n);
                real gamma_k = gamma_fac_k * q_di(d3, k, n);

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

    const auto &rho_pi = refstate.rho_pi.data;
    const auto &q_pi = refstate.q_pi.data;
    const auto &rho_di = refstate.rho_di.data;
    const auto &q_di = refstate.q_di.data;

    parallel_for(
        "Prepare rhs 1 - B",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndensity_dycore> rhs0;
          compute_H2bar_ext<ndensity_dycore, diff_ord, vert_diff_ord>(
              rhs0, rhs_dens, this->primal_geometry, this->dual_geometry, pis,
              pjs, pks, i, j, k, n);

          for (int d1 = 0; d1 < ndensity_dycore; ++d1) {
            real b_d1 = 0;
            for (int d2 = 0; d2 < ndensity_dycore; ++d2) {
              b_d1 -= dtf * refstate.Blin_coeff(d1, d2, k, n) * rhs0(d2);
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
          compute_wD0<ndensity_dycore>(mod_v, q_pi, bvar, pis, pjs, pks, i, j,
                                       k, n);
          v_transform(k, j, i, n) =
              rhs_v(0, k + pks, j + pjs, i + pis, n) + mod_v(0);
          if (k < primal_topology.nl) {
            real mod_w = compute_wD0_vert<ndensity_dycore>(q_di, bvar, pis, pjs,
                                                           pks, i, j, k, n);
            w_transform(k, j, i, n) =
                rhs_w(0, k + pks, j + pjs, i + pis, n) + mod_w;
          }
        });

    yakl::timer_start("fft fwd");
    fftv_x.forward_real(v_transform);
    fftw_x.forward_real(w_transform);
    // fftv_y.forward_real(v_transform);
    // fftw_y.forward_real(w_transform);
    yakl::timer_stop("fft fwd");

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

          real fH2bar_k = fourier_H2bar_ext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, 0,
              n_cells_x, n_cells_y, dual_topology.ni);
          real fH2bar_kp1 = fourier_H2bar_ext<diff_ord>(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k + 1, 0,
              n_cells_x, n_cells_y, dual_topology.ni);

          SArray<real, 1, ndims> fH1_kp1_a;
          SArray<real, 1, ndims> fH1_k_a;
          fourier_H1_ext<diff_ord>(fH1_kp1_a, primal_geometry, dual_geometry,
                                   pis, pjs, pks, i, j, k + 1, 0, n_cells_x,
                                   n_cells_y, dual_topology.ni);
          fourier_H1_ext<diff_ord>(fH1_k_a, primal_geometry, dual_geometry, pis,
                                   pjs, pks, i, j, k, 0, n_cells_x, n_cells_y,
                                   dual_topology.ni);
          real fH1h_kp1 = fH1_kp1_a(0);
          real fH1h_k = fH1_k_a(0);

          complex fD1bar_kp1 = fourier_D1bar(1, i, j, k + 1, n_cells_x,
                                             n_cells_y, dual_topology.ni);
          complex fD1bar_k =
              fourier_D1bar(1, i, j, k, n_cells_x, n_cells_y, dual_topology.ni);

          real he_kp1 = rho_pi(0, k + 1, n);
          real he_k = rho_pi(0, k, n);

          for (int d1 = 0; d1 < ndensity_dycore; ++d1) {
            for (int d2 = 0; d2 < ndensity_dycore; ++d2) {
              real alpha_kp1 = dtf2 * q_di(d1, k + 1, n);
              complex beta_kp1 =
                  fH2bar_kp1 * refstate.Blin_coeff(d1, d2, k + 1, n) *
                  q_pi(d2, k + 1, n) * fD1bar_kp1 * he_kp1 * fH1h_kp1;
              complex beta_k = fH2bar_k * refstate.Blin_coeff(d1, d2, k, n) *
                               q_pi(d2, k, n) * fD1bar_k * he_k * fH1h_k;
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

          real gamma_fac_kp1 = rho_di(0, k + 1, n) *
                               H1_vert_coeff(primal_geometry, dual_geometry,
                                             pis, pjs, pks, i, j, k + 1, n);
          real gamma_fac_k =
              rho_di(0, k, n) * H1_vert_coeff(primal_geometry, dual_geometry,
                                              pis, pjs, pks, i, j, k, n);

          complex_vrhs(k, j, i, n) *= complex_vcoeff(0, k, j, i, n);
          for (int d1 = 0; d1 < ndensity_dycore; ++d1) {
            complex_vrhs(k, j, i, n) +=
                complex_vcoeff(1 + d1, k, j, i, n) *
                (gamma_fac_kp1 * q_di(d1, k + 1, n) * w_kp1 -
                 gamma_fac_k * q_di(d1, k, n) * w_k);
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

    yakl::timer_start("fft bwd");
    fftv_x.inverse_real(v_transform);
    fftw_x.inverse_real(w_transform);
    // fftv_y.inverse_real(v_transform);
    // fftw_y.inverse_real(w_transform);
    yakl::timer_stop("fft bwd");

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

    YAKL_SCOPE(Hk, ::Hk);
    parallel_for(
        "Recover densities 1 - F/Fw",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 1, ndims> u;
          compute_H1_ext<1, diff_ord>(u, sol_v, this->primal_geometry,
                                      this->dual_geometry, dis, djs, dks, i, j,
                                      k, n);

          fvar(0, k + dks, j + djs, i + dis, n) = u(0) * rho_pi(0, k, n);

          if (k < dual_topology.ni - 2) {
            const real uw = compute_H1_vert(sol_w, this->primal_geometry,
                                            this->dual_geometry, dis, djs, dks,
                                            i, j, k + 1, n);
            fwvar(0, k + 1 + dks, j + djs, i + dis, n) =
                uw * rho_di(0, k + 1, n);
          }
        });

    auxiliary_vars.fields_arr[FWVAR].set_bnd(0.0);
    auxiliary_vars.exchange({FVAR, FWVAR});

    parallel_for(
        "Recover densities 2 - D1bar",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD1bar<ndensity>(sol_dens, q_pi, fvar, dis, djs, dks, i, j, k,
                                   n);
          compute_wD1bar_vert<ndensity, ADD_MODE::ADD>(
              sol_dens, q_di, fwvar, dis, djs, dks, i, j, k, n);
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
  aux_desc_arr[FVAR] = {"F", dtopo, ndims - 1, 1, 1}; // F = twisted
                                                      // (n-1,1)-form

  aux_desc_arr[BVAR] = {"B", ptopo, 0, 0,
                        ndensity_B}; // B = straight (0,0)-form

  aux_desc_arr[KVAR] = {"K", dtopo, ndims, 1, 1}; // K = twisted (n,1)-form

  aux_desc_arr[UVAR] = {"U", dtopo, ndims - 1, 1, 1}; // U = twisted
                                                      // (n-1,1)-form

  aux_desc_arr[FWVAR] = {"Fw", dtopo, ndims, 0, 1}; // Fw = twisted (n,0)-form

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
  aux_desc_arr[DENSVERTRECONVAR] = {
      "densvertrecon", dtopo, ndims, 0,
      ndensity}; // densvertrecon lives on vert dual edges, associated with Fw

  // fct stuff- Phi, Mf, edgeflux
  aux_desc_arr[MFVAR] = {"Mf", dtopo, ndims, 1, ndensity};
  aux_desc_arr[EDGEFLUXVAR] = {"edgeflux", dtopo, ndims - 1, 1, ndensity};
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
                         Parallel &par, PamCoupler &coupler,
                         std::unique_ptr<TestCase> &testcase) {

  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  int nz = coupler.get_nz();

  params.acoustic_balance = config["balance_initial_density"].as<bool>(false);
  params.uniform_vertical = (config["vcoords"].as<std::string>() == "uniform");

  // Read diffusion coefficients
  params.entropicvar_diffusion_coeff =
      config["entropicvar_diffusion_coeff"].as<real>(0);
  params.velocity_diffusion_coeff =
      config["velocity_diffusion_coeff"].as<real>(0);

  // Read the data initialization options
  params.initdataStr = config["initData"].as<std::string>();
  testcase_from_string(testcase, params.initdataStr, params.acoustic_balance, coupler);

  serial_print("IC: " + params.initdataStr, par.masterproc);
  serial_print("acoustically balanced: " +
                   std::to_string(params.acoustic_balance),
               par.masterproc);
  serial_print("entropicvar_diffusion_coeff: " +
                   std::to_string(params.entropicvar_diffusion_coeff),
               par.masterproc);
  serial_print("velocity_diffusion_coeff: " +
                   std::to_string(params.velocity_diffusion_coeff),
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

//THIS IS VERY 2D SPECIFIC, NEEDS TO CHANGE FOR 3D
  params.ylen = 1.0;
  params.yc = 0.5;
  testcase->set_domain(params, coupler);

  // Store vertical cell interface heights in the data manager
  auto &dm = coupler.get_data_manager_readonly();
  params.zint = dm.get<real const, 2>("vertical_interface_height");

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

// TODO: Restore this once there is an option to autogenerate a "uniform" grid
// in the vertical

// template <class T> class SWETestCase : public TestCase, public T {
// public:
//   using T::g;
//
//   using T::Lx;
//   using T::Ly;
//   using T::xc;
//   using T::yc;
//
//   using T::h_f;
// #ifdef _TSWE
//   using T::S_f;
// #endif
//   using T::coriolis_f;
//   using T::v_f;
//
//   void set_domain(ModelParameters &params) override {
//     params.xlen = Lx;
//     params.zlen = Ly;
//     params.xc = xc;
//   }
//
//   void set_initial_conditions(FieldSet<nprognostic> &progvars,
//                               FieldSet<nconstant> &constvars,
//                               const Geometry<Straight> &primal_geom,
//                               const Geometry<Twisted> &dual_geom) override {
//
//     dual_geom.set_11form_values(
//         YAKL_LAMBDA(real x, real y) { return h_f(x, y); },
//         progvars.fields_arr[DENSVAR], 0);
// #ifdef _TSWE
//     dual_geom.set_11form_values(
//         YAKL_LAMBDA(real x, real y) { return S_f(x, y); },
//         progvars.fields_arr[DENSVAR], 1);
// #endif
//     primal_geom.set_10form_values(
//         YAKL_LAMBDA(real x, real y) { return v_f(x, y); },
//         progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
//     primal_geom.set_01form_values(
//         YAKL_LAMBDA(real x, real y) { return v_f(x, y); },
//         progvars.fields_arr[WVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
//
//     // RESTORE ONCE PRIMAL GRID CFV RECON IS FIXED IN PARALLEL
//     // primal_geom.set_11form_values(
//     //      YAKL_LAMBDA(real x, real y) { return coriolis_f(x, y); },
//     //      constvars.fields_arr[CORIOLISVAR], 0);
//
//     YAKL_SCOPE(tracer_f, this->tracer_f);
//     for (int i = 0; i < ntracers_dycore; i++) {
//       dual_geom.set_11form_values(
//           YAKL_LAMBDA(real x, real y) {
//             return h_f(x, y) * tracer_f(i)->compute(x, y, Lx, Ly, xc, yc);
//           },
//           progvars.fields_arr[DENSVAR], i + ndensity_dycore);
//     }
//     Hs.set_parameters(g);
//   }
// };

class CouplerData : public TestCase
{
real g:
real Lx;
real xc;

std::array<real, 3> get_domain() const { return {Lx, 1, -1.}; }

void set_domain(ModelParameters &params) {

  g = coupler.get_grav();
  Lx = coupler.get_xlen();
  xc = Lx/2.0_fp;

  params.xlen = Lx;
  params.xc = xc;
}

void add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {};

void set_initial_conditions(FieldSet<nprognostic> &progvars,
                            FieldSet<nconstant> &constvars,
                            const Geometry<Straight> &primal_geom,
                            const Geometry<Twisted> &dual_geom) override {

dual_geom.set_11form_values(YAKL_LAMBDA(real x, real z) { return flat_geop(x, z, g); }, constvars.fields_arr[HSVAR], 0);

//Just copies data from coupler state
convert_coupler_to_dynamics_state(coupler, progvars,constvars);

}

//ADD THESE
void set_reference_state(ReferenceState &reference_state,
                         const Geometry<Straight> &primal_geom,
                         const Geometry<Twisted> &dual_geom) override {
  auto &refstate = static_cast<ModelReferenceState &>(reference_state);

//Needs the non-perturbed average coupler state

// Get GCM state to base reference state on
    auto &dm = coupler.get_data_manager_readwrite();

    auto gcm_rho_d = dm.get<real,2>("gcm_density_dry");
    //auto gcm_uvel  = dm.get<real,2>("gcm_uvel"       );
    //auto gcm_vvel  = dm.get<real,2>("gcm_vvel"       );
    //auto gcm_wvel  = dm.get<real,2>("gcm_wvel"       );
    auto gcm_temp  = dm.get<real,2>("gcm_temp"       );
    auto gcm_rho_v = dm.get<real,2>("gcm_water_vapor");

//LOTS TO ADD HERE
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

  std::array<real, 3> get_domain() const override { return {Lx, 1, Lz}; }

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
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

    const int pks = primal_topology.ks;

    YAKL_SCOPE(thermo, ::thermo);
    YAKL_SCOPE(Hs, ::Hs);

    dual_geom.set_profile_11form_values(
        YAKL_LAMBDA(real z) { return flat_geop(0, z, g); }, refstate.geop, 0);
    dual_geom.set_profile_11form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.dens, 0);
    dual_geom.set_profile_11form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.dens, 1);

    parallel_for(
        "compute rho_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          SArray<real, 1, 1> rho0;
          compute_H2bar_ext<1, vert_diff_ord>(
              rho0, refstate.dens.data, primal_geom, dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k, n) = rho0(0);
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.q_pi, 0);

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.q_pi, 1);

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < ndensity; ++l) {
            refstate.q_pi.data(l, k, n) /= refstate.rho_pi.data(0, k, n);
          }
        });

    parallel_for(
        "compute rho_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, 0, n) = refstate.rho_pi.data(0, 0, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1, n);
          } else {
            refstate.rho_di.data(0, k, n) =
                0.5_fp * (refstate.rho_pi.data(0, k, n) +
                          refstate.rho_pi.data(0, k - 1, n));
          }
        });

    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.q_di, 0);

    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.q_di, 1);

    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < ndensity; ++l) {
            refstate.q_di.data(l, k, n) /= refstate.rho_di.data(0, k, n);
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
        "Compute Blin_coeff",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          real z = primal_geom.zint(k + primal_topology.ks, n);

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
  using T::entropicvar_f;
  using T::refentropicdensity_f;
  using T::refnsq_f;
  using T::refrho_f;
  using T::refrhov_f;

  std::array<real, 3> get_domain() const override { return {Lx, 1, Lz}; }

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    YAKL_SCOPE(thermo, ::thermo);
    YAKL_SCOPE(varset, ::varset);

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

  void set_reference_state(ReferenceState &reference_state,
                           const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = static_cast<ModelReferenceState &>(reference_state);

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    const int pks = primal_topology.ks;

    YAKL_SCOPE(thermo, ::thermo);
    YAKL_SCOPE(Hs, ::Hs);
    YAKL_SCOPE(varset, ::varset);

    dual_geom.set_profile_11form_values(
        YAKL_LAMBDA(real z) { return flat_geop(0, z, g); }, refstate.geop, 0);
    dual_geom.set_profile_11form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.dens, 0);
    dual_geom.set_profile_11form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.dens, 1);
    dual_geom.set_profile_11form_values(
        YAKL_LAMBDA(real z) { return refrhov_f(z, thermo); }, refstate.dens, 2);

    const auto total_density_f =
        YAKL_LAMBDA(const real3d &densvar, int d, int k, int n) {
      return varset.get_total_density(densvar, k, 0, n);
    };

    parallel_for(
        "compute rho_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          SArray<real, 1, 1> rho0;
          compute_H2bar_ext<1, vert_diff_ord>(total_density_f, rho0,
                                              refstate.dens.data, primal_geom,
                                              dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k, n) = rho0(0);
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.q_pi, 0);

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.q_pi, 1);

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refrhov_f(z, thermo); }, refstate.q_pi, 2);

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < ndensity; ++l) {
            refstate.q_pi.data(l, k, n) /= refstate.rho_pi.data(0, k, n);
          }
        });

    parallel_for(
        "compute rho_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, 0, n) = refstate.rho_pi.data(0, 0, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1, n);
          } else {
            refstate.rho_di.data(0, k, n) =
                0.5_fp * (refstate.rho_pi.data(0, k, n) +
                          refstate.rho_pi.data(0, k - 1, n));
          }
        });

    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.q_di, 0);

    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.q_di, 1);

    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refrhov_f(z, thermo); }, refstate.q_di, 2);

    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < ndensity; ++l) {
            refstate.q_di.data(l, k, n) /= refstate.rho_di.data(0, k, n);
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
        "Compute Blin_coeff",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          real z = primal_geom.zint(k + primal_topology.ks, n);

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

// struct DoubleVortex {
//   static real constexpr g = 9.80616_fp;
//   static real constexpr Lx = 5000000._fp;
//   static real constexpr Ly = 5000000._fp;
//   static real constexpr coriolis = 0.00006147_fp;
//   static real constexpr H0 = 750.0_fp;
//   static real constexpr ox = 0.1_fp;
//   static real constexpr oy = 0.1_fp;
//   static real constexpr sigmax = 3._fp / 40._fp * Lx;
//   static real constexpr sigmay = 3._fp / 40._fp * Ly;
//   static real constexpr dh = 75.0_fp;
//   static real constexpr xc1 = (0.5_fp - ox) * Lx;
//   static real constexpr yc1 = (0.5_fp - oy) * Ly;
//   static real constexpr xc2 = (0.5_fp + ox) * Lx;
//   static real constexpr yc2 = (0.5_fp + oy) * Ly;
//   static real constexpr xc = 0.5_fp * Lx;
//   static real constexpr yc = 0.5_fp * Ly;
//   static real constexpr c = 0.05_fp;
//   static real constexpr a = 1.0_fp / 3.0_fp;
//   static real constexpr D = 0.5_fp * Lx;
//
//   static real YAKL_INLINE coriolis_f(real x, real y) { return coriolis; }
//
//   static real YAKL_INLINE h_f(real x, real y) {
//     real xprime1 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc1));
//     real yprime1 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc1));
//     real xprime2 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc2));
//     real yprime2 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc2));
//     real xprimeprime1 =
//         Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc1));
//     real yprimeprime1 =
//         Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc1));
//     real xprimeprime2 =
//         Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc2));
//     real yprimeprime2 =
//         Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc2));
//
//     return H0 - dh * (exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1))
//     +
//                       exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2))
//                       - 4._fp * pi * sigmax * sigmay / Lx / Ly);
//   }
//
//   static vecext<2> YAKL_INLINE v_f(real x, real y) {
//     vecext<2> vvec;
//
//     real xprime1 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc1));
//     real yprime1 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc1));
//     real xprime2 = Lx / (pi * sigmax) * sin(pi / Lx * (x - xc2));
//     real yprime2 = Ly / (pi * sigmay) * sin(pi / Ly * (y - yc2));
//     real xprimeprime1 =
//         Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc1));
//     real yprimeprime1 =
//         Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc1));
//     real xprimeprime2 =
//         Lx / (2.0_fp * pi * sigmax) * sin(2.0_fp * pi / Lx * (x - xc2));
//     real yprimeprime2 =
//         Ly / (2.0_fp * pi * sigmay) * sin(2.0_fp * pi / Ly * (y - yc2));
//
//     vvec.u =
//         -g * dh / coriolis / sigmay *
//         (yprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 *
//         yprime1)) +
//          yprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 *
//          yprime2)));
//     vvec.w =
//         g * dh / coriolis / sigmax *
//         (xprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 *
//         yprime1)) +
//          xprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 *
//          yprime2)));
//     return vvec;
//   }
//
//   static real YAKL_INLINE S_f(real x, real y) {
//     // real sval = g * (1. + c * sin(2. * M_PI / Lx * (x - xc)) * sin(2. *
//     M_PI
//     // / Ly * (y - yc)) * exp(-((x-xc)*(x-xc) + (y-yc)*(y-yc))/(a*a*D*D)));
//     real sval =
//         g * (1._fp + c * exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) /
//                              (a * a * D * D)));
//     // real sval = g * (1. + c * sin(2. * M_PI / Lx * (x- xc)));
//     // real sval = g;
//     // real sval = g * (1. + c * ((x > 0.35 * Lx && x < 0.65 * Lx && y > 0.35
//     *
//     // Ly
//     // && y < 0.65 * Ly ) ? 1. : 0.));
//     return sval * h_f(x, y);
//   }
// };

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

    return rho_ref * thermo.compute_entropic_var(p_ref, T_ref, 1, 0, 0, 0);
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
    real dtheta = (r < rc) ? dss * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
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

struct TwoBubbles {
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
    return rho_ref * thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    return isentropic_rho(x, z, theta0, g, thermo);
  }

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(x, z, theta0, g, thermo);
    real T = isentropic_T(x, z, theta0, g, thermo);

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

struct DensityCurrent {
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
    return rho_ref * thermo.compute_entropic_var(p_ref, T_ref, 0, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real z, const ThermoPotential &thermo) {
    real rho_b = isentropic_rho(x, z, theta0, g, thermo);
    return rho_b;
  }

  static real YAKL_INLINE entropicvar_f(real x, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(x, z, theta0, g, thermo);
    real T = isentropic_T(x, z, theta0, g, thermo);
    real r = sqrt((x - bxc) * (x - bxc) / (bxr * bxr) +
                  (z - bzc) * (z - bzc) / (bzr * bzr));
    real dtheta = (r < 1) ? dss * 0.5_fp * (1._fp + cos(pi * r)) : 0._fp;
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
    real rh = (r < rc) ? rh0 * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
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

  static real YAKL_INLINE refrhov_f(real z, const ThermoPotential &thermo) {
    return 0;
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

  static real YAKL_INLINE refrhov_f(real z, const ThermoPotential &thermo) {
    return 0;
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

template <bool add_perturbation> struct GravityWave {
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

    real T = T_ref;
    real rho = rho_ref;

    if (add_perturbation) {
      T += dT;
      rho += drho;
    }

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
    real rho = refrho_f(z, thermo);
    if (add_perturbation) {
      rho += sum_series(x, z, t, thermo).drho;
    }
    return rho;
  }

  static real YAKL_INLINE entropicdensityexact_f(
      real x, real z, real t, const ThermoPotential &thermo) {

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

    return rho * thermo.compute_entropic_var(p, T, 0, 0, 0, 0);
  }

  static real YAKL_INLINE Texact_f(real x, real z, real t,
                                   const ThermoPotential &thermo) {
    real T = T_ref;
    if (add_perturbation) {
      T += sum_series(x, z, t, thermo).dT;
    }
    return T;
  }

  static vecext<2> YAKL_INLINE vexact_f(real x, real z, real t,
                                        const ThermoPotential &thermo) {
    vecext<2> v;

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

// Diagnostics

  struct ExactDensityDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "dense";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const ReferenceState &reference_state,
                 const FieldSet<nconstant> &const_vars,
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

    void compute(real time, const ReferenceState &reference_state,
                 const FieldSet<nconstant> &const_vars,
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

    void compute(real time, const ReferenceState &reference_state,
                 const FieldSet<nconstant> &const_vars,
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

    void compute(real time, const ReferenceState &reference_state,
                 const FieldSet<nconstant> &const_vars,
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
    static constexpr bool linear_T = true;
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom) override {
      name = "T";
      topology = pgeom.topology;
      dofs_arr = {0, 0, 1};
      Diagnostic::initialize(pgeom, dgeom);
    }

    void compute(real time, const ReferenceState &reference_state,
                 const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      const auto &primal_topology = primal_geometry.topology;

      int pis = primal_topology.is;
      int pjs = primal_topology.js;
      int pks = primal_topology.ks;

      const auto &densvar = x.fields_arr[DENSVAR].data;
      const auto &refstate =
          static_cast<const ModelReferenceState &>(reference_state);
      const auto &refdens = refstate.dens.data;

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

            real T;
            if (!linear_T) {
              T = thermo.compute_T(alpha, entropic_var, 1, 0, 0, 0);
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

void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance, PamCoupler &coupler) {
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
  } else if (name == "largerisingbubble") {
    testcase = std::make_unique<EulerTestCase<LargeRisingBubble>>();
  } else if (name == "moistlargerisingbubble") {
    testcase = std::make_unique<MoistEulerTestCase<MoistLargeRisingBubble>>();
  }
  else if (name == "coupler") {
    testcase = std::make_unique<CouplerData>();
  }
  // else if (name == "doublevortex") {
  //  testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  //}
  else {
    throw std::runtime_error("unknown test case");
  }
}

//I don't think gets used anywhere?
//void testcase_from_config(std::unique_ptr<TestCase> &testcase,
//                          const YAML::Node &config) {
//  const std::string name = config["initData"].as<std::string>();
//  const bool acoustic_balance =
//      config["balance_initial_density"].as<bool>(false);
//  testcase_from_string(testcase, name, acoustic_balance);
//}
