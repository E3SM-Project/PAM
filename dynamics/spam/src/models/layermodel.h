#pragma once

#include "common.h"
#include "model.h"
#include "stats.h"

#include "ext_deriv.h"
#include "fct.h"
#include "hamiltonian.h"
#include "hodge_star.h"
#include "pocketfft_hdronly.h"
#include "recon.h"
#include "thermo.h"
#include "wedge.h"

// *******   Functionals/Hamiltonians   ***********//

Functional_PVPE PVPE;
Hamiltonian_Hk Hk;

VariableSet varset;
#ifdef _SWE
Hamiltonian_SWE_Hs Hs;
#elif _TSWE
Hamiltonian_TSWE_Hs Hs;
#endif

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
    dofs_arr = {0, 1, ndensity}; // densldiag = straight 0-form
    Diagnostic::initialize(ptopo, dtopo, pgeom, dgeom);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "Compute Dens0 Diag",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_I<ndensity, diff_ord>(
              field.data, x.fields_arr[DENSVAR].data, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
        });
  }
};

class Q0Diagnostic : public Diagnostic {
public:
  void initialize(const Topology &ptopo, const Topology &dtopo,
                  const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom) override {
    name = "q";
    topology = dtopo;
    dofs_arr = {0, 1, 1}; // qdiag = twisted 0-form
    Diagnostic::initialize(ptopo, dtopo, pgeom, dgeom);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(PVPE, ::PVPE);
    parallel_for(
        "Compute Q0 Diag",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_q0(field.data, x.fields_arr[VVAR].data,
                          x.fields_arr[DENSVAR].data,
                          const_vars.fields_arr[CORIOLISVAR].data, dis, djs,
                          dks, i, j, k, n);
        });
  }
};

void add_model_diagnostics(
    std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
  diagnostics.emplace_back(std::make_unique<Dens0Diagnostic>());
  diagnostics.emplace_back(std::make_unique<Q0Diagnostic>());
}

// *******   Tendencies   ***********//

class ModelTendencies : public Tendencies {
public:
  void initialize(PamCoupler &coupler, ModelParameters &params,
                  Topology &primal_topo, Topology &dual_topo,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nconstant> &const_exchange) {
    Tendencies::initialize(params, primal_topo, dual_topo, primal_geom,
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

  void YAKL_INLINE compute_dens0(real5d dens0var, const real5d densvar) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    // compute dens0var = I densvar
    parallel_for(
        "Compute Dens0",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_I<ndensity, diff_ord>(
              dens0var, densvar, this->primal_geometry, this->dual_geometry,
              pis, pjs, pks, i, j, k, n);
        });
  }

  void YAKL_INLINE compute_U(real5d Uvar, const real5d Vvar) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    // compute U = H v
    parallel_for(
        "Compute U",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H<1, diff_ord>(Uvar, Vvar, *this->primal_geometry,
                                 *this->dual_geometry, dis, djs, dks, i, j, k,
                                 n);
        });
  }

  void YAKL_INLINE compute_U_and_q0f0(real5d Uvar, real5d Q0var, real5d f0var,
                                      const real5d Vvar, const real5d densvar,
                                      const real5d coriolisvar) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(PVPE, ::PVPE);
    // compute U Hv, q0, f0
    parallel_for(
        "Compute U, Q0, F0",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_H<1, diff_ord>(Uvar, Vvar, this->primal_geometry,
                                 this->dual_geometry, dis, djs, dks, i, j, k,
                                 n);
          PVPE.compute_q0f0(Q0var, f0var, Vvar, densvar, coriolisvar, dis, djs,
                            dks, i, j, k, n);
        });
  }

  void YAKL_INLINE compute_F_and_K(real5d Fvar, real5d Kvar, const real5d Vvar,
                                   const real5d Uvar, const real5d dens0var) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(Hk, ::Hk);
    parallel_for(
        "Compute F/K",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F_and_K(Fvar, Kvar, Vvar, Uvar, dens0var, dis, djs, dks, i,
                             j, k, n);
        });
  }

  void YAKL_INLINE compute_F_and_he(real5d Fvar, real5d HEvar,
                                    const real5d Uvar, const real5d dens0var) {

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute F/he",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F_and_he(Fvar, HEvar, Uvar, dens0var, dis, djs, dks, i, j,
                              k, n);
        });
  }

  void YAKL_INLINE compute_B(real5d Bvar, const real5d Kvar,
                             const real5d densvar, const real5d HSvar) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(Hk, ::Hk);
    YAKL_SCOPE(Hs, ::Hs);
    parallel_for(
        "Compute B",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hs.compute_dHsdx(Bvar, densvar, HSvar, pis, pjs, pks, i, j, k, n);
          Hk.compute_dKddens(Bvar, Kvar, pis, pjs, pks, i, j, k, n);
        });
  }

  void YAKL_INLINE compute_FT(real5d FTvar, const real5d Fvar) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "Compute FT",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_W(FTvar, Fvar, pis, pjs, pks, i, j, k, n);
        });
  }

  void YAKL_INLINE compute_edge_reconstructions(
      real5d densedgereconvar, real5d Qedgereconvar, real5d fedgereconvar,
      const real5d dens0var, const real5d Q0var, const real5d f0var) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute straight edge recons",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_straight_edge_recon<1, reconstruction_type,
                                      reconstruction_order>(
              Qedgereconvar, Q0var, pis, pjs, pks, i, j, k, n, primal_wenoRecon,
              primal_to_gll, primal_wenoIdl, primal_wenoSigma);
          compute_straight_edge_recon<1, coriolis_reconstruction_type,
                                      coriolis_reconstruction_order>(
              fedgereconvar, f0var, pis, pjs, pks, i, j, k, n,
              coriolis_wenoRecon, coriolis_to_gll, coriolis_wenoIdl,
              coriolis_wenoSigma);
        });

    parallel_for(
        "Compute twisted edge recons",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_edge_recon<ndensity, dual_reconstruction_type,
                                     dual_reconstruction_order>(
              densedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
        });
  }

  void compute_recons(real5d densreconvar, real5d Qreconvar,
                      real5d Coriolisreconvar, const real5d densedgereconvar,
                      const real5d Qedgereconvar, const real5d fedgereconvar,
                      const real5d HEvar, const real5d FTvar,
                      const real5d Uvar) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute straight recons",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_straight_recon<1, reconstruction_type>(
              Qreconvar, Qedgereconvar, FTvar, pis, pjs, pks, i, j, k, n);
          compute_straight_recon<1, coriolis_reconstruction_type>(
              Coriolisreconvar, fedgereconvar, FTvar, pis, pjs, pks, i, j, k,
              n);
        });

    parallel_for(
        "Compute twisted recons",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_recon<ndensity, dual_reconstruction_type>(
              densreconvar, densedgereconvar, Uvar, dis, djs, dks, i, j, k, n);
          // scale primal recons
          for (int d = 0; d < ndims; d++) {
            for (int l = 0; l < ndensity; l++) {
              densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) =
                  densreconvar(l + d * ndensity, k + dks, j + djs, i + dis, n) /
                  HEvar(d, k + dks, j + djs, i + dis, n);
            }
          }
        });
  }

  void compute_tendencies(real5d denstendvar, real5d Vtendvar,
                          const real5d densreconvar, const real5d Qreconvar,
                          const real5d Coriolisreconvar, const real5d Bvar,
                          const real5d Fvar, const real5d Phivar) {

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute v tend",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD1_fct<ndensity>(Vtendvar, densreconvar, Phivar, Bvar, pis,
                                    pjs, pks, i, j, k, n);
          if (qf_choice == QF_MODE::EC) {
            compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, pis, pjs,
                                           pks, i, j, k, n);
          }
          if (qf_choice == QF_MODE::NOEC) {
            compute_Q_nonEC<1, ADD_MODE::ADD>(Vtendvar, Qreconvar, Fvar, pis,
                                              pjs, pks, i, j, k, n);
          }
          compute_Q_EC<1, ADD_MODE::ADD>(Vtendvar, Coriolisreconvar, Fvar, pis,
                                         pjs, pks, i, j, k, n);
        });

    parallel_for(
        "Compute dens tend",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wDbar2_fct<ndensity>(denstendvar, densreconvar, Phivar, Fvar,
                                       dis, djs, dks, i, j, k, n);
        });
  }

  void YAKL_INLINE compute_functional_derivatives(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x,
      FieldSet<nauxiliary> &auxiliary_vars) override {

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U(auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[VVAR].data);

    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
        auxiliary_vars.fields_arr[DENS0VAR]);

    compute_F_and_K(auxiliary_vars.fields_arr[FVAR].data,
                    auxiliary_vars.fields_arr[KVAR].data,
                    x.fields_arr[VVAR].data,
                    auxiliary_vars.fields_arr[UVAR].data,
                    auxiliary_vars.fields_arr[DENS0VAR].data);

    this->aux_exchange->exchanges_arr[FVAR].exchange_field(
        auxiliary_vars.fields_arr[FVAR]);
    this->aux_exchange->exchanges_arr[KVAR].exchange_field(
        auxiliary_vars.fields_arr[KVAR]);

    compute_B(auxiliary_vars.fields_arr[BVAR].data,
              auxiliary_vars.fields_arr[KVAR].data, x.fields_arr[DENSVAR].data,
              const_vars.fields_arr[HSVAR].data);

    this->aux_exchange->exchanges_arr[BVAR].exchange_field(
        auxiliary_vars.fields_arr[BVAR]);
  }

  void YAKL_INLINE apply_symplectic(real dt, FieldSet<nconstant> &const_vars,
                                    FieldSet<nprognostic> &x,
                                    FieldSet<nauxiliary> &auxiliary_vars,
                                    FieldSet<nprognostic> &xtend) override {
    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U_and_q0f0(auxiliary_vars.fields_arr[UVAR].data,
                       auxiliary_vars.fields_arr[Q0VAR].data,
                       auxiliary_vars.fields_arr[F0VAR].data,
                       x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data,
                       const_vars.fields_arr[CORIOLISVAR].data);

    this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
        auxiliary_vars.fields_arr[DENS0VAR]);
    this->aux_exchange->exchanges_arr[UVAR].exchange_field(
        auxiliary_vars.fields_arr[UVAR]);
    this->aux_exchange->exchanges_arr[Q0VAR].exchange_field(
        auxiliary_vars.fields_arr[Q0VAR]);
    this->aux_exchange->exchanges_arr[F0VAR].exchange_field(
        auxiliary_vars.fields_arr[F0VAR]);

    compute_F_and_he(auxiliary_vars.fields_arr[FVAR2].data,
                     auxiliary_vars.fields_arr[HEVAR].data,
                     auxiliary_vars.fields_arr[UVAR].data,
                     auxiliary_vars.fields_arr[DENS0VAR].data);

    this->aux_exchange->exchanges_arr[FVAR2].exchange_field(
        auxiliary_vars.fields_arr[FVAR2]);
    this->aux_exchange->exchanges_arr[HEVAR].exchange_field(
        auxiliary_vars.fields_arr[HEVAR]);

    compute_FT(auxiliary_vars.fields_arr[FTVAR].data,
               auxiliary_vars.fields_arr[FVAR2].data);

    this->aux_exchange->exchanges_arr[FTVAR].exchange_field(
        auxiliary_vars.fields_arr[FTVAR]);

    // Compute densrecon, qrecon and frecon
    compute_edge_reconstructions(
        auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
        auxiliary_vars.fields_arr[QEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data,
        auxiliary_vars.fields_arr[Q0VAR].data,
        auxiliary_vars.fields_arr[F0VAR].data);

    this->aux_exchange->exchanges_arr[DENSEDGERECONVAR].exchange_field(
        auxiliary_vars.fields_arr[DENSEDGERECONVAR]);
    this->aux_exchange->exchanges_arr[QEDGERECONVAR].exchange_field(
        auxiliary_vars.fields_arr[QEDGERECONVAR]);
    this->aux_exchange->exchanges_arr[CORIOLISEDGERECONVAR].exchange_field(
        auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR]);

    compute_recons(auxiliary_vars.fields_arr[DENSRECONVAR].data,
                   auxiliary_vars.fields_arr[QRECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
                   auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[QEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[HEVAR].data,
                   auxiliary_vars.fields_arr[FTVAR].data,
                   auxiliary_vars.fields_arr[UVAR].data);

    this->aux_exchange->exchanges_arr[DENSRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[DENSRECONVAR]);
    this->aux_exchange->exchanges_arr[QRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[QRECONVAR]);
    this->aux_exchange->exchanges_arr[CORIOLISRECONVAR].exchange_field(
        auxiliary_vars.fields_arr[CORIOLISRECONVAR]);

    // Compute fct

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    parallel_for(
        "Compute edgefluxes",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_edgefluxes<ndensity>(
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[DENSRECONVAR].data,
              auxiliary_vars.fields_arr[FVAR].data, dis, djs, dks, i, j, k, n);
        });
    this->aux_exchange->exchanges_arr[EDGEFLUXVAR].exchange_field(
        auxiliary_vars.fields_arr[EDGEFLUXVAR]);

    parallel_for(
        "Compute Mf",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Mf<ndensity>(auxiliary_vars.fields_arr[MFVAR].data,
                               auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt,
                               dis, djs, dks, i, j, k, n);
        });

    this->aux_exchange->exchanges_arr[MFVAR].exchange_field(
        auxiliary_vars.fields_arr[MFVAR]);

    parallel_for(
        "Compute Phi",
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
        for (int d = 0; d < ndims; ++d) {
          auxiliary_vars.fields_arr[PHIVAR].set(l + d * ndensity, 1.0);
        }
      }
    }

    this->aux_exchange->exchanges_arr[PHIVAR].exchange_field(
        auxiliary_vars.fields_arr[PHIVAR]);

    // Compute tendencies
    compute_tendencies(xtend.fields_arr[DENSVAR].data,
                       xtend.fields_arr[VVAR].data,
                       auxiliary_vars.fields_arr[DENSRECONVAR].data,
                       auxiliary_vars.fields_arr[QRECONVAR].data,
                       auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
                       auxiliary_vars.fields_arr[BVAR].data,
                       auxiliary_vars.fields_arr[FVAR].data,
                       auxiliary_vars.fields_arr[PHIVAR].data);
  }
};

// *******   Linear system   ***********//
typedef yakl::Array<complex, 2, yakl::memDevice, yakl::styleC> complex2d;
class ModelLinearSystem : public LinearSystem {
public:
  complex2d ch, cv0, cv1;

  void initialize(ModelParameters &params, Topology &primal_topo,
                  Topology &dual_topo, Geometry &primal_geom,
                  Geometry &dual_geom, ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nprognostic> &prog_exchange) override {
    LinearSystem::initialize(params, primal_topo, dual_topo, primal_geom,
                             dual_geom, aux_exchange, prog_exchange);

    ch = complex2d("ch", dual_topo.n_cells_y, dual_topo.n_cells_x);
    cv0 = complex2d("cv0", primal_topo.n_cells_y, primal_topo.n_cells_x);
    cv1 = complex2d("cv1", primal_topo.n_cells_y, primal_topo.n_cells_x);
  }

  virtual void YAKL_INLINE solve(real dt, FieldSet<nprognostic> &rhs,
                                 FieldSet<nconstant> &const_vars,
                                 FieldSet<nauxiliary> &auxiliary_vars,
                                 FieldSet<nprognostic> &solution) override {

    auto grav = Hs.g;
    // TODO: get this from somwhere
    real ref_height = 750;
    auto n_cells_x = dual_topology->n_cells_x;
    auto n_cells_y = dual_topology->n_cells_y;

    auto nl = dual_topology->nl;
    int pis = primal_topology->is;
    int pjs = primal_topology->js;
    int pks = primal_topology->ks;
    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;

    auto v_rhs = rhs.fields_arr[VVAR].data;
    auto h_rhs = rhs.fields_arr[DENSVAR].data;

    constexpr bool single_eq = true;

    real scale = 1.0 / (n_cells_x * n_cells_y);
    pocketfft::shape_t shape(2);
    pocketfft::stride_t stride(2);
    shape[0] = n_cells_y;
    shape[1] = n_cells_x;
    stride[0] = n_cells_x * sizeof(complex);
    stride[1] = sizeof(complex);
    pocketfft::shape_t axes = {0, 1};

    solution.copy(rhs);
    prog_exchange->exchange_variable_set(solution);

    if (single_eq) {
      auto u = auxiliary_vars.fields_arr[UVAR].data;
      auto h0 = auxiliary_vars.fields_arr[DENS0VAR].data;
      auto h_sol = solution.fields_arr[DENSVAR].data;

      parallel_for(
          "linsolve compute h0",
          SimpleBounds<4>(primal_topology->nl, primal_topology->n_cells_y,
                          primal_topology->n_cells_x, primal_topology->nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_I<ndensity_active, diff_ord>(
                h0, h_sol, *this->primal_geometry, *this->dual_geometry, pis,
                pjs, pks, i, j, k, n);
          });

      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
          auxiliary_vars.fields_arr[DENS0VAR]);

      parallel_for(
          "linsolve ???",
          SimpleBounds<4>(primal_topology->nl, primal_topology->n_cells_y,
                          primal_topology->n_cells_x, primal_topology->nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            SArray<real, 1, ndensity_active> c;
            c(0) = 0;
            for (int dof = 1; dof < ndensity_active; ++dof) {
              c(dof) = -0.25_fp * dt;
            }
            compute_cwD1<ndensity_active, ADD_MODE::ADD>(v_rhs, c, h0, pis, pjs,
                                                         pks, i, j, k, n);
          });

      prog_exchange->exchanges_arr[VVAR].exchange_field(rhs.fields_arr[VVAR]);

      parallel_for(
          "linsolve compute u",
          SimpleBounds<4>(dual_topology->nl, dual_topology->n_cells_y,
                          dual_topology->n_cells_x, dual_topology->nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_H<1, diff_ord>(u, v_rhs, *this->primal_geometry,
                                   *this->dual_geometry, dis, djs, dks, i, j, k,
                                   n);
          });
      this->aux_exchange->exchanges_arr[UVAR].exchange_field(
          auxiliary_vars.fields_arr[UVAR]);

      parallel_for(
          "linsolve prepare h rhs",
          SimpleBounds<4>(dual_topology->nl, dual_topology->n_cells_y,
                          dual_topology->n_cells_x, dual_topology->nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            real h = h_rhs(0, dks + k, djs + j, dis + i, n);
            compute_cwDbar2<1>(h_rhs, ref_height, u, dis, djs, dks, i, j, k, n);
            h_rhs(0, dks + k, djs + j, dis + i, n) *= -0.5_fp * dt;
            h_rhs(0, dks + k, djs + j, dis + i, n) += h;
          });
    }

    if (single_eq) {
      parallel_for(
          "fft copy in", SimpleBounds<2>(n_cells_y, n_cells_x),
          YAKL_LAMBDA(int j, int i) {
            ch(j, i) = h_rhs(0, dks, j + djs, i + dis, 0);
          });
    } else {
      parallel_for(
          "fft copy in", SimpleBounds<2>(n_cells_y, n_cells_x),
          YAKL_LAMBDA(int j, int i) {
            ch(j, i) = h_rhs(0, dks, j + djs, i + dis, 0);
            cv0(j, i) = v_rhs(0, pks, j + pjs, i + pis, 0);
            cv1(j, i) = v_rhs(1, pks, j + pjs, i + pis, 0);
          });
    }

    yakl::timer_start("fft fwd");
    // fftw_execute(fwd_h);
    pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD, ch.data(),
                   ch.data(), 1._fp);
    if (!single_eq) {
      pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD,
                     cv0.data(), cv0.data(), 1._fp);
      pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD,
                     cv1.data(), cv1.data(), 1._fp);
    }
    yakl::timer_stop("fft fwd");

    parallel_for(
        "fft invert", SimpleBounds<2>(n_cells_y, n_cells_x),
        YAKL_LAMBDA(int j, int i) {
          complex cI =
              fourier_I<diff_ord>(*primal_geometry, *dual_geometry, pis, pjs,
                                  pks, i, j, 0, 0, n_cells_x, n_cells_y, nl);
          SArray<complex, 1, ndims> cD1;
          fourier_cwD1(cD1, grav, i, j, 0, n_cells_x, n_cells_y, nl);
          SArray<complex, 1, ndims> cH;
          fourier_H<diff_ord>(cH, *primal_geometry, *dual_geometry, pis, pjs,
                              pks, i, j, 0, 0, n_cells_x, n_cells_y, nl);
          SArray<complex, 1, ndims> cD2bar;
          fourier_cwDbar2(cD2bar, ref_height, i, j, 0, n_cells_x, n_cells_y,
                          nl);

          complex hn = (ch(j, i) - 0.5_fp * dt * cD2bar(0) * cH(0) * cv0(j, i) -
                        0.5_fp * dt * cD2bar(1) * cH(1) * cv1(j, i));
          complex hd =
              (1._fp - 0.25_fp * dt * dt * cD1(0) * cD2bar(0) * cI * cH(0) -
               0.25_fp * dt * dt * cD1(1) * cD2bar(1) * cI * cH(1));

          if (single_eq) {
            ch(j, i) /= hd;
          } else {
            ch(j, i) = hn / hd;
            cv0(j, i) -= 0.5_fp * dt * cD1(0) * cI * ch(j, i);
            cv1(j, i) -= 0.5_fp * dt * cD1(1) * cI * ch(j, i);
          }
        });

    yakl::timer_start("fft bwd");
    pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD, ch.data(),
                   ch.data(), scale);
    if (!single_eq) {
      pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD,
                     cv0.data(), cv0.data(), scale);
      pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD,
                     cv1.data(), cv1.data(), scale);
    }
    yakl::timer_stop("fft bwd");

    auto h_sol = solution.fields_arr[DENSVAR].data;
    auto v_sol = solution.fields_arr[VVAR].data;

    parallel_for(
        "fft copy out", SimpleBounds<2>(n_cells_y, n_cells_x),
        YAKL_LAMBDA(int j, int i) {
          h_sol(0, dks, j + djs, i + dis, 0) = ch(j, i).real();

          if (!single_eq) {
            v_sol(0, dks, j + djs, i + dis, 0) = cv0(j, i).real();
            v_sol(1, dks, j + djs, i + dis, 0) = cv1(j, i).real();
          }
        });

    if (single_eq) {

      prog_exchange->exchanges_arr[DENSVAR].exchange_field(
          solution.fields_arr[DENSVAR]);

      auto h0 = auxiliary_vars.fields_arr[DENS0VAR].data;
      parallel_for(
          "linsolve compute h0",
          SimpleBounds<4>(primal_topology->nl, primal_topology->n_cells_y,
                          primal_topology->n_cells_x, primal_topology->nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_I<1, diff_ord>(h0, h_sol, *this->primal_geometry,
                                   *this->dual_geometry, pis, pjs, pks, i, j, k,
                                   n);
          });

      this->aux_exchange->exchanges_arr[DENS0VAR].exchange_field(
          auxiliary_vars.fields_arr[DENS0VAR]);

      parallel_for(
          "linsolve extract v from h",
          SimpleBounds<4>(primal_topology->nl, primal_topology->n_cells_y,
                          primal_topology->n_cells_x, primal_topology->nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            real v0 = v_rhs(0, dks, j + djs, i + dis, 0);
            real v1 = v_rhs(1, dks, j + djs, i + dis, 0);
            SArray<real, 1, ndensity_active> c;
            c(0) = grav;
            for (int dof = 1; dof < ndensity_active; ++dof) {
              c(dof) = 0.5_fp;
            }
            compute_cwD1<ndensity_active>(v_sol, c, h0, pis, pjs, pks, i, j, k,
                                          n);

            v_sol(0, dks, j + djs, i + dis, 0) *= -0.5_fp * dt;
            v_sol(1, dks, j + djs, i + dis, 0) *= -0.5_fp * dt;
            v_sol(0, dks, j + djs, i + dis, 0) += v0;
            v_sol(1, dks, j + djs, i + dis, 0) += v1;
          });
    }
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

      YAKL_SCOPE(Hs, ::Hs);
      YAKL_SCOPE(Hk, ::Hk);
      parallel_for(
          "Compute energy stats",
          SimpleBounds<3>(dual_topology.nl, dual_topology.n_cells_y,
                          dual_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int k, int j, int i) {
            real KE, PE, IE;
            KE = Hk.compute_KE(progvars.fields_arr[VVAR].data,
                               progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, k, n);
            PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, k, n);
            IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, k, n);
            TEarr(k, j, i) = KE + PE + IE;
            KEarr(k, j, i) = KE;
            PEarr(k, j, i) = PE;
            IEarr(k, j, i) = IE;
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
          "Compute pv/pens stats",
          SimpleBounds<3>(primal_topology.nl, primal_topology.n_cells_y,
                          primal_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int k, int j, int i) {
            pvpe vals_pvpe;
            vals_pvpe =
                PVPE.compute_PVPE(progvars.fields_arr[VVAR].data,
                                  progvars.fields_arr[DENSVAR].data,
                                  constvars.fields_arr[CORIOLISVAR].data, pis,
                                  pjs, pks, i, j, k, n);
            PVarr(k, j, i) = vals_pvpe.pv;
            PENSarr(k, j, i) = vals_pvpe.pe;
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

  // v, dens
  prog_topo_arr[VVAR] = ptopo;
  prog_topo_arr[DENSVAR] = dtopo;
  prog_names_arr[VVAR] = "v";
  prog_names_arr[DENSVAR] = "dens";
  set_dofs_arr(prog_ndofs_arr, VVAR, 1, 1, 1); // v = straight 1-form
  set_dofs_arr(prog_ndofs_arr, DENSVAR, ndims, 1,
               ndensity); // dens = twisted n-form

  // hs, coriolis
  const_topo_arr[HSVAR] = dtopo;
  const_topo_arr[CORIOLISVAR] = ptopo;
  const_names_arr[HSVAR] = "hs";
  const_names_arr[CORIOLISVAR] = "coriolis";
  set_dofs_arr(const_ndofs_arr, HSVAR, ndims, 1, 1);   // hs = twisted n-form
  set_dofs_arr(const_ndofs_arr, CORIOLISVAR, 2, 1, 1); // f = straight 2-form

  // functional derivatives = F, B, K, he, U
  aux_topo_arr[BVAR] = ptopo;
  aux_topo_arr[FVAR] = dtopo;
  aux_topo_arr[UVAR] = dtopo;
  aux_topo_arr[HEVAR] = dtopo;
  aux_topo_arr[KVAR] = dtopo;
  aux_names_arr[KVAR] = "K";
  aux_names_arr[BVAR] = "B";
  aux_names_arr[FVAR] = "F";
  aux_names_arr[FVAR2] = "F2";
  aux_names_arr[UVAR] = "U";
  aux_names_arr[HEVAR] = "he";
  set_dofs_arr(aux_ndofs_arr, BVAR, 0, 1, ndensity);   // B = straight 0-form
  set_dofs_arr(aux_ndofs_arr, KVAR, ndims, 1, 1);      // K = twisted n-form
  set_dofs_arr(aux_ndofs_arr, FVAR, ndims - 1, 1, 1);  // F = twisted (n-1)-form
  set_dofs_arr(aux_ndofs_arr, FVAR2, ndims - 1, 1, 1); // F = twisted (n-1)-form
  set_dofs_arr(aux_ndofs_arr, UVAR, ndims - 1, 1, 1);  // U = twisted (n-1)-form
  set_dofs_arr(aux_ndofs_arr, HEVAR, ndims - 1, 1,
               1); // he lives on dual edges, associated with F

  // dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_topo_arr[DENSRECONVAR] = dtopo;
  aux_topo_arr[DENSEDGERECONVAR] = dtopo;
  aux_topo_arr[DENS0VAR] = ptopo;
  aux_names_arr[DENS0VAR] = "dens0";
  aux_names_arr[DENSRECONVAR] = "densrecon";
  aux_names_arr[DENSEDGERECONVAR] = "densedgerecon";
  set_dofs_arr(aux_ndofs_arr, DENSRECONVAR, ndims - 1, 1,
               ndensity); // densrecon lives on dual edges, associated with F
  set_dofs_arr(
      aux_ndofs_arr, DENSEDGERECONVAR, ndims, 1,
      2 * ndims *
          ndensity); // densedgerecon lives on dual cells, associated with F
  set_dofs_arr(aux_ndofs_arr, DENS0VAR, 0, 1,
               ndensity); // dens0 = straight 0-form

  // dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon,
  // coriolisedgercon, coriolisrecon
  aux_topo_arr[FTVAR] = ptopo;
  aux_topo_arr[CORIOLISRECONVAR] = ptopo;
  aux_topo_arr[CORIOLISEDGERECONVAR] = ptopo;
  aux_topo_arr[Q0VAR] = dtopo;
  aux_topo_arr[F0VAR] = dtopo;
  aux_topo_arr[QRECONVAR] = ptopo;
  aux_topo_arr[QEDGERECONVAR] = ptopo;
  aux_names_arr[FTVAR] = "FT";
  aux_names_arr[CORIOLISRECONVAR] = "coriolisrecon";
  aux_names_arr[CORIOLISEDGERECONVAR] = "coriolisedgerecon";
  aux_names_arr[Q0VAR] = "q";
  aux_names_arr[F0VAR] = "f";
  aux_names_arr[QRECONVAR] = "qrecon";
  aux_names_arr[QEDGERECONVAR] = "qedgerecon";
  set_dofs_arr(aux_ndofs_arr, FTVAR, 1, 1, 1); // FT = straight 1-form
  set_dofs_arr(aux_ndofs_arr, Q0VAR, 0, 1, 1); // q0 = twisted 0-form
  set_dofs_arr(aux_ndofs_arr, F0VAR, 0, 1, 1); // f0 = twisted 0-form
  set_dofs_arr(aux_ndofs_arr, QRECONVAR, 1, 1,
               1); // qrecon lives on primal edges, associated with FT
  set_dofs_arr(aux_ndofs_arr, QEDGERECONVAR, 2, 1,
               4); // qedgerecon lives on primal cells
  set_dofs_arr(aux_ndofs_arr, CORIOLISRECONVAR, 1, 1,
               1); // coriolisrecon lives on primal edges, associated with FT
  set_dofs_arr(aux_ndofs_arr, CORIOLISEDGERECONVAR, 2, 1,
               4); // coriolisedgerecon lives on primal cells

  // fct stuff- Phi, Mf, edgeflux
  aux_topo_arr[PHIVAR] = dtopo;
  aux_topo_arr[MFVAR] = dtopo;
  aux_topo_arr[EDGEFLUXVAR] = dtopo;
  aux_names_arr[PHIVAR] = "Phi";
  aux_names_arr[MFVAR] = "Mf";
  aux_names_arr[EDGEFLUXVAR] = "edgeflux";
  set_dofs_arr(aux_ndofs_arr, PHIVAR, ndims - 1, 1, ndensity);
  set_dofs_arr(aux_ndofs_arr, MFVAR, ndims, 1, ndensity);
  set_dofs_arr(aux_ndofs_arr, EDGEFLUXVAR, ndims - 1, 1, ndensity);
}

void testcase_from_string(std::unique_ptr<TestCase> &testcase,
                          std::string name);

void readModelParamsFile(std::string inFile, ModelParameters &params,
                         Parallel &par, int nz,
                         std::unique_ptr<TestCase> &testcase) {

  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  // Read the data initialization options
  params.initdataStr = config["initData"].as<std::string>();
  testcase_from_string(testcase, params.initdataStr);

  serial_print("IC: " + params.initdataStr, par.masterproc);

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

  params.zlen = 1.0;
  params.zc = 0.5;
  testcase->set_domain(params);

  readParamsFile(inFile, params, par, nz);
}

//***************** Test Cases ***************************//

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
    params.ylen = Ly;
    params.xc = xc;
    params.yc = yc;
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    dual_geom.set_2form_values(
        YAKL_LAMBDA(real x, real y) { return h_f(x, y); },
        progvars.fields_arr[DENSVAR], 0);
#ifdef _TSWE
    dual_geom.set_2form_values(
        YAKL_LAMBDA(real x, real y) { return S_f(x, y); },
        progvars.fields_arr[DENSVAR], 1);
#endif
    primal_geom.set_1form_values(
        YAKL_LAMBDA(real x, real y) { return v_f(x, y); },
        progvars.fields_arr[VVAR], 0, LINE_INTEGRAL_TYPE::TANGENT);
    primal_geom.set_2form_values(
        YAKL_LAMBDA(real x, real y) { return coriolis_f(x, y); },
        constvars.fields_arr[CORIOLISVAR], 0);

    YAKL_SCOPE(tracer_f, this->tracer_f);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_2form_values(
          YAKL_LAMBDA(real x, real y) {
            return h_f(x, y) * tracer_f(i)->compute(x, y, Lx, Ly, xc, yc);
          },
          progvars.fields_arr[DENSVAR], i + ndensity_dycore);
    }
    Hs.set_parameters(g);
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

  static vec<2> YAKL_INLINE v_f(real x, real y) {
    vec<2> vvec;

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
    vvec.v =
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

void testcase_from_string(std::unique_ptr<TestCase> &testcase,
                          std::string name) {
  if (name == "doublevortex") {
    testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  } else {
    throw std::runtime_error("unknown test case");
  }
}
