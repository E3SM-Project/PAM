#pragma once

#include "common.h"
#include "model.h"
#include "refstate.h"
#include "stats.h"

#include "ext_deriv.h"
#include "fct.h"
#include "hamiltonian.h"
#include "hodge_star.h"
#include "recon.h"
#include "thermo.h"
#include "wedge.h"

// *******   Diagnostics   ***********//

class Dens0Diagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    // concentration 0-forms for dens
    name = "densl";
    topology = pgeom.topology;
    dofs_arr = {
        0, 1, VariableSet::ndensity_prognostic}; // densldiag = straight 0-form
    Diagnostic::initialize(pgeom, dgeom, equations);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &primal_topology = primal_geometry.topology;
    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "Compute Dens0 Diag",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_H2bar<VariableSet::ndensity_prognostic, diff_ord>(
              field.data, x.fields_arr[DENSVAR].data, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
        });
  }
};

class Q0Diagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    name = "q";
    topology = dgeom.topology;
    dofs_arr = {0, 1, 1}; // qdiag = twisted 0-form
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
  using VS = VariableSet;

  void initialize(ModelParameters &params, Equations &equations,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {

    Tendencies::initialize(params, equations, primal_geom, dual_geom);
  }

  void
  convert_dynamics_to_coupler_state(PamCoupler &coupler,
                                    const FieldSet<nprognostic> &prog_vars,
                                    const FieldSet<nconstant> &const_vars) {
    this->equations->varset.convert_dynamics_to_coupler_state(
        coupler, prog_vars, const_vars);
  }
  void convert_coupler_to_dynamics_state(PamCoupler &coupler,
                                         FieldSet<nprognostic> &prog_vars,
                                         FieldSet<nauxiliary> &auxiliary_vars,
                                         FieldSet<nconstant> &const_vars) {
    this->equations->varset.convert_coupler_to_dynamics_state(
        coupler, prog_vars, const_vars);
  }

  void compute_constants(FieldSet<nconstant> &const_vars,
                         FieldSet<nprognostic> &x) override {}

  void compute_dens0(real5d dens0var, const real5d densvar) {
    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    // compute dens0var = I densvar
    parallel_for(
        "Compute Dens0",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H2bar<VS::ndensity_prognostic, diff_ord>(
              dens0var, densvar, primal_geometry, dual_geometry, pis, pjs, pks,
              i, j, k, n);
        });
  }

  void compute_U(real5d Uvar, const real5d Vvar) {
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);

    // compute U = H1 v
    parallel_for(
        "Compute U",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H1<1, diff_ord>(Uvar, Vvar, primal_geometry, dual_geometry,
                                  dis, djs, dks, i, j, k, n);
        });
  }

  void compute_U_and_q0f0(real5d Uvar, real5d Q0var, real5d f0var,
                          const real5d Vvar, const real5d densvar,
                          const real5d coriolisvar) {
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(PVPE, this->equations->PVPE);
    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    // compute U H1v, q0, f0
    parallel_for(
        "Compute U, Q0, F0",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H1<1, diff_ord>(Uvar, Vvar, primal_geometry, dual_geometry,
                                  dis, djs, dks, i, j, k, n);
          PVPE.compute_q0f0(Q0var, f0var, Vvar, densvar, coriolisvar, dis, djs,
                            dks, i, j, k, n);
        });
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void compute_F_and_K(real fac, real5d Fvar, real5d Kvar, const real5d Vvar,
                       const real5d Uvar, const real5d dens0var) {
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(Hk, this->equations->Hk);
    parallel_for(
        "Compute F/K",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F_and_K<addmode>(Fvar, Kvar, Vvar, Uvar, dens0var, dis,
                                      djs, dks, i, j, k, n, fac);
        });
  }

  void compute_F_and_he(real5d Fvar, real5d HEvar, const real5d Uvar,
                        const real5d dens0var) {
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(Hk, this->equations->Hk);
    parallel_for(
        "Compute F/he",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hk.compute_F_and_he(Fvar, HEvar, Uvar, dens0var, dis, djs, dks, i, j,
                              k, n);
        });
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void compute_B(real fac, real5d Bvar, const real5d Kvar, const real5d densvar,
                 const real5d HSvar) {
    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    YAKL_SCOPE(Hk, this->equations->Hk);
    YAKL_SCOPE(Hs, this->equations->Hs);
    parallel_for(
        "Compute B",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          Hs.compute_dHsdx<addmode>(Bvar, densvar, HSvar, pis, pjs, pks, i, j,
                                    k, n, fac);
          Hk.compute_dKddens(Bvar, Kvar, pis, pjs, pks, i, j, k, n, fac);
        });
  }

  void compute_FT(real5d FTvar, const real5d Fvar) {
    const auto &primal_topology = primal_geometry.topology;

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

  void compute_edge_reconstructions(real5d densedgereconvar,
                                    real5d Qedgereconvar, real5d fedgereconvar,
                                    const real5d dens0var, const real5d Q0var,
                                    const real5d f0var) {
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(primal_wenoRecon, this->primal_wenoRecon);
    YAKL_SCOPE(primal_to_gll, this->primal_to_gll);
    YAKL_SCOPE(primal_wenoIdl, this->primal_wenoIdl);
    YAKL_SCOPE(primal_wenoSigma, this->primal_wenoSigma);

    YAKL_SCOPE(coriolis_wenoRecon, this->coriolis_wenoRecon);
    YAKL_SCOPE(coriolis_to_gll, this->coriolis_to_gll);
    YAKL_SCOPE(coriolis_wenoIdl, this->coriolis_wenoIdl);
    YAKL_SCOPE(coriolis_wenoSigma, this->coriolis_wenoSigma);

    parallel_for(
        "Compute straight edge recons",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
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

    YAKL_SCOPE(dual_wenoRecon, this->dual_wenoRecon);
    YAKL_SCOPE(dual_to_gll, this->dual_to_gll);
    YAKL_SCOPE(dual_wenoIdl, this->dual_wenoIdl);
    YAKL_SCOPE(dual_wenoSigma, this->dual_wenoSigma);

    parallel_for(
        "Compute twisted edge recons",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_edge_recon<VS::ndensity_prognostic,
                                     dual_reconstruction_type,
                                     dual_reconstruction_order>(
              densedgereconvar, dens0var, dis, djs, dks, i, j, k, n,
              dual_wenoRecon, dual_to_gll, dual_wenoIdl, dual_wenoSigma);
        });
  }

  void compute_recons(real5d densreconvar, real5d Qreconvar,
                      real5d Coriolisreconvar, const real5d densedgereconvar,
                      const real5d Qedgereconvar, const real5d fedgereconvar,
                      const real5d HEvar, const real5d FTvar,
                      const real5d Vvar) {
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);

    parallel_for(
        "Compute straight recons",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_straight_recon<1, reconstruction_type>(
              Qreconvar, Qedgereconvar, primal_geometry, FTvar, pis, pjs, pks,
              i, j, k, n);
          compute_straight_recon<1, coriolis_reconstruction_type>(
              Coriolisreconvar, fedgereconvar, primal_geometry, FTvar, pis, pjs,
              pks, i, j, k, n);
        });

    parallel_for(
        "Compute twisted recons",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_twisted_recon<VS::ndensity_prognostic,
                                dual_reconstruction_type>(
              densreconvar, densedgereconvar, dual_geometry, Vvar, dis, djs,
              dks, i, j, k, n);
          // scale primal recons
          for (int d = 0; d < ndims; d++) {
            for (int l = 0; l < VS::ndensity_prognostic; l++) {
              densreconvar(d + l * ndims, k + dks, j + djs, i + dis, n) =
                  densreconvar(d + l * ndims, k + dks, j + djs, i + dis, n) /
                  HEvar(d, k + dks, j + djs, i + dis, n);
            }
          }
        });
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void compute_tendencies(real5d denstendvar, real5d Vtendvar,
                          const real5d densreconvar, const real5d Qreconvar,
                          const real5d Coriolisreconvar, const real5d Bvar,
                          const real5d Fvar) {
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(active_dens_ids, this->equations->varset.active_dens_ids);

    parallel_for(
        "Compute v tend",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_wD0<VS::ndensity_active, addmode>(Vtendvar, densreconvar,
                                                    active_dens_ids, Bvar, pis,
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
          compute_wDnm1bar<VS::ndensity_prognostic, addmode>(
              denstendvar, densreconvar, Fvar, dis, djs, dks, i, j, k, n);
        });
  }

  void compute_functional_derivatives(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x,
      FieldSet<nauxiliary> &auxiliary_vars, real fac = 1,
      ADD_MODE addmode = ADD_MODE::REPLACE) override {
    const auto &dual_topology = dual_geometry.topology;

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U(auxiliary_vars.fields_arr[UVAR].data, x.fields_arr[VVAR].data);

    auxiliary_vars.exchange({UVAR, DENS0VAR});

    // doing it this way because virtual functions cannot be templates :(
    if (addmode == ADD_MODE::ADD) {
      compute_F_and_K<ADD_MODE::ADD>(fac, auxiliary_vars.fields_arr[FVAR].data,
                                     auxiliary_vars.fields_arr[KVAR].data,
                                     x.fields_arr[VVAR].data,
                                     auxiliary_vars.fields_arr[UVAR].data,
                                     auxiliary_vars.fields_arr[DENS0VAR].data);
    } else if (addmode == ADD_MODE::REPLACE) {
      compute_F_and_K<ADD_MODE::REPLACE>(
          fac, auxiliary_vars.fields_arr[FVAR].data,
          auxiliary_vars.fields_arr[KVAR].data, x.fields_arr[VVAR].data,
          auxiliary_vars.fields_arr[UVAR].data,
          auxiliary_vars.fields_arr[DENS0VAR].data);
    }

    auxiliary_vars.exchange({FVAR, KVAR});

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
                        FieldSet<nprognostic> &xtend,
                        ADD_MODE addmode = ADD_MODE::REPLACE) override {
    const auto &dual_topology = dual_geometry.topology;

    compute_dens0(auxiliary_vars.fields_arr[DENS0VAR].data,
                  x.fields_arr[DENSVAR].data);
    compute_U_and_q0f0(auxiliary_vars.fields_arr[UVAR].data,
                       auxiliary_vars.fields_arr[Q0VAR].data,
                       auxiliary_vars.fields_arr[F0VAR].data,
                       x.fields_arr[VVAR].data, x.fields_arr[DENSVAR].data,
                       const_vars.fields_arr[CORIOLISVAR].data);

    auxiliary_vars.exchange({DENS0VAR, UVAR, Q0VAR, F0VAR});

    compute_F_and_he(auxiliary_vars.fields_arr[FVAR2].data,
                     auxiliary_vars.fields_arr[HEVAR].data,
                     auxiliary_vars.fields_arr[UVAR].data,
                     auxiliary_vars.fields_arr[DENS0VAR].data);

    auxiliary_vars.exchange({FVAR2, HEVAR});

    compute_FT(auxiliary_vars.fields_arr[FTVAR].data,
               auxiliary_vars.fields_arr[FVAR2].data);

    auxiliary_vars.exchange({FTVAR});

    // Compute densrecon, qrecon and frecon
    compute_edge_reconstructions(
        auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
        auxiliary_vars.fields_arr[QEDGERECONVAR].data,
        auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
        auxiliary_vars.fields_arr[DENS0VAR].data,
        auxiliary_vars.fields_arr[Q0VAR].data,
        auxiliary_vars.fields_arr[F0VAR].data);

    auxiliary_vars.exchange(
        {DENSEDGERECONVAR, QEDGERECONVAR, CORIOLISEDGERECONVAR});

    compute_recons(auxiliary_vars.fields_arr[DENSRECONVAR].data,
                   auxiliary_vars.fields_arr[QRECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
                   auxiliary_vars.fields_arr[DENSEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[QEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[CORIOLISEDGERECONVAR].data,
                   auxiliary_vars.fields_arr[HEVAR].data,
                   auxiliary_vars.fields_arr[FTVAR].data,
                   x.fields_arr[VVAR].data);

    auxiliary_vars.exchange({DENSRECONVAR, QRECONVAR, CORIOLISRECONVAR});

    // Compute fct

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    YAKL_SCOPE(dens_pos, this->equations->varset.dens_pos);
    parallel_for(
        "Compute edgefluxes",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_edgefluxes<VS::ndensity_prognostic>(
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data,
              auxiliary_vars.fields_arr[DENSRECONVAR].data,
              auxiliary_vars.fields_arr[FVAR].data, dens_pos, dis, djs, dks, i,
              j, k, n);
        });

    auxiliary_vars.exchange({EDGEFLUXVAR});

    parallel_for(
        "Compute Mf",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_Mf<VS::ndensity_prognostic>(
              auxiliary_vars.fields_arr[MFVAR].data,
              auxiliary_vars.fields_arr[EDGEFLUXVAR].data, dt, dens_pos, dis,
              djs, dks, i, j, k, n);
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
        });

    auxiliary_vars.exchange({DENSRECONVAR});

    // Compute tendencies
    if (addmode == ADD_MODE::REPLACE) {
      compute_tendencies<ADD_MODE::REPLACE>(
          xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
          auxiliary_vars.fields_arr[DENSRECONVAR].data,
          auxiliary_vars.fields_arr[QRECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
          auxiliary_vars.fields_arr[BVAR].data,
          auxiliary_vars.fields_arr[FVAR].data);
    }
    if (addmode == ADD_MODE::ADD) {
      compute_tendencies<ADD_MODE::ADD>(
          xtend.fields_arr[DENSVAR].data, xtend.fields_arr[VVAR].data,
          auxiliary_vars.fields_arr[DENSRECONVAR].data,
          auxiliary_vars.fields_arr[QRECONVAR].data,
          auxiliary_vars.fields_arr[CORIOLISRECONVAR].data,
          auxiliary_vars.fields_arr[BVAR].data,
          auxiliary_vars.fields_arr[FVAR].data);
    }
  }
};

// *******   Linear system   ***********//
class ModelLinearSystem : public LinearSystem {
public:
  using VS = VariableSet;

  real4d dens_transform;
  yakl::RealFFT1D<real> fft_y;
  yakl::RealFFT1D<real> fft_x;
  int nxf, nyf;

  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  Equations &equations) override {
    LinearSystem::initialize(params, primal_geom, dual_geom, equations);

    const auto &dual_topo = dual_geom.topology;

    int nx = dual_topo.n_cells_x;
    this->nxf = nx + 2 - nx % 2;
    int ny = dual_topo.n_cells_y;
    this->nyf = ny + 2 - ny % 2;

    dens_transform =
        real4d("dens transform", dual_topo.nl, nyf, nxf, dual_topo.nens);
    yakl::memset(dens_transform, 0);

    fft_x.init(dens_transform, 2, nx);
    fft_y.init(dens_transform, 1, ny);
  }

  virtual void solve(real dt, FieldSet<nprognostic> &rhs,
                     FieldSet<nconstant> &const_vars,
                     FieldSet<nauxiliary> &auxiliary_vars,
                     FieldSet<nprognostic> &solution) override {

    const auto &dual_topology = dual_geometry.topology;
    const auto &primal_topology = primal_geometry.topology;
    const auto &refstate = this->equations->reference_state;

    yakl::timer_start("Linear solve");
    auto grav = this->equations->Hs.g;

    auto n_cells_x = dual_topology.n_cells_x;
    auto n_cells_y = dual_topology.n_cells_y;
    auto nl = dual_topology.nl;
    auto nens = dual_topology.nens;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;
    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    solution.copy(rhs);
    solution.exchange();

    auto v_rhs = rhs.fields_arr[VVAR].data;
    auto v_sol = solution.fields_arr[VVAR].data;
    auto dens_rhs = rhs.fields_arr[DENSVAR].data;
    auto dens_sol = solution.fields_arr[DENSVAR].data;
    auto U = auxiliary_vars.fields_arr[UVAR].data;
    auto dens0 = auxiliary_vars.fields_arr[DENS0VAR].data;

    YAKL_SCOPE(primal_geometry, this->primal_geometry);
    YAKL_SCOPE(dual_geometry, this->dual_geometry);
    YAKL_SCOPE(dens_transform, this->dens_transform);

    parallel_for(
        "compute dens0",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H2bar<VS::ndensity_active, diff_ord>(
              dens0, dens_sol, primal_geometry, dual_geometry, pis, pjs, pks, i,
              j, k, n);
        });

    auxiliary_vars.exchange({DENS0VAR});

    parallel_for(
        "prepare rhs 1",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, VS::ndensity_active, ndims> c;
          for (int d = 0; d < ndims; ++d) {
#ifdef _SWE
            c(0, d) = 0;
#elif _TSWE
            c(0, d) = 0.25_fp * grav * dt;
#endif
            for (int dof = 1; dof < VS::ndensity_active; ++dof) {
              c(dof, d) = -0.25_fp * dt;
            }
          }
          compute_wD0<VS::ndensity_active, ADD_MODE::ADD>(v_rhs, c, dens0, pis,
                                                          pjs, pks, i, j, k, n);
        });

    rhs.exchange({VVAR});

    parallel_for(
        "prepare rhs 2",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H1<1, diff_ord>(U, v_rhs, primal_geometry, dual_geometry, dis,
                                  djs, dks, i, j, k, n);
        });
    auxiliary_vars.exchange({UVAR});

    parallel_for(
        "prepare rhs 3",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real h = dens_rhs(0, dks + k, djs + j, dis + i, n);
          SArray<real, 3, 1, ndims, 2> c;
          for (int d = 0; d < ndims; ++d) {
            for (int n = 0; n < 2; ++n) {
              c(0, d, n) = refstate.ref_height;
            }
          }
          compute_wDnm1bar<1>(dens_rhs, c, U, dis, djs, dks, i, j, k, n);
          dens_rhs(0, dks + k, djs + j, dis + i, n) *= -0.5_fp * dt;
          dens_rhs(0, dks + k, djs + j, dis + i, n) += h;
        });

    parallel_for(
        "fft copy in",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          dens_transform(k, j, i, n) =
              dens_rhs(0, k + dks, j + djs, i + dis, n);
        });

    yakl::timer_start("fft fwd");
    fft_x.forward_real(dens_transform);
    fft_y.forward_real(dens_transform);
    yakl::timer_stop("fft fwd");

    parallel_for(
        "fft invert",
        SimpleBounds<4>(dual_topology.nl, nyf, nxf, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          int ik = i / 2;
          int jk = j / 2;

          real cH2bar = fourier_H2bar<diff_ord>(primal_geometry, dual_geometry,
                                                pis, pjs, pks, ik, jk, 0, 0,
                                                n_cells_x, n_cells_y, nl);

          SArray<real, 1, ndims> cH1;
          fourier_H1<diff_ord>(cH1, primal_geometry, dual_geometry, pis, pjs,
                               pks, ik, jk, 0, 0, n_cells_x, n_cells_y, nl);

          SArray<real, 1, ndims> cD0Dnm1bar;
          fourier_cwD0Dnm1bar(cD0Dnm1bar, grav * refstate.ref_height, ik, jk, 0,
                              n_cells_x, n_cells_y, nl);

          real hd =
              (1._fp - 0.25_fp * dt * dt * cD0Dnm1bar(0) * cH2bar * cH1(0) -
               0.25_fp * dt * dt * cD0Dnm1bar(1) * cH2bar * cH1(1));

          dens_transform(k, j, i, n) /= hd;
        });

    yakl::timer_start("fft bwd");
    fft_x.inverse_real(dens_transform);
    fft_y.inverse_real(dens_transform);
    yakl::timer_stop("fft bwd");

    parallel_for(
        "fft copy out",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real dens_old = dens_sol(0, k + dks, j + djs, i + dis, n);
          real dens_new = dens_transform(k, j, i, n);
          dens_sol(0, k + dks, j + djs, i + dis, n) = dens_new;
#ifdef _TSWE
          dens_sol(1, k + dks, j + djs, i + dis, n) -=
              grav * (dens_old - dens_new);
#endif
        });

    solution.exchange({DENSVAR});

    parallel_for(
        "compute dens0",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H2bar<1, diff_ord>(dens0, dens_sol, primal_geometry,
                                     dual_geometry, pis, pjs, pks, i, j, k, n);
        });

    auxiliary_vars.exchange({DENS0VAR});

    parallel_for(
        "extract v from dens",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          real v0 = v_rhs(0, k + dks, j + djs, i + dis, n);
          real v1 = v_rhs(1, k + dks, j + djs, i + dis, n);
          SArray<real, 2, VS::ndensity_active, ndims> c;
          for (int d = 0; d < ndims; ++d) {
#ifdef _SWE
            c(0, d) = grav;
#elif _TSWE
            c(0, d) = 0.5_fp * grav;
#endif
            for (int dof = 1; dof < VS::ndensity_active; ++dof) {
              c(dof, d) = 0.5_fp;
            }
            compute_wD0<VS::ndensity_active>(v_sol, c, dens0, pis, pjs, pks, i,
                                             j, k, n);
          }

          v_sol(0, k + dks, j + djs, i + dis, n) *= -0.5_fp * dt;
          v_sol(1, k + dks, j + djs, i + dis, n) *= -0.5_fp * dt;
          v_sol(0, k + dks, j + djs, i + dis, n) += v0;
          v_sol(1, k + dks, j + djs, i + dis, n) += v1;
        });

    yakl::timer_stop("Linear solve");
  }
};

// *******   Statistics   ***********//

class ModelStats : public Stats {
public:
  using VS = VariableSet;
  real3d TEarr, KEarr, PEarr, IEarr, PVarr, PENSarr, trimmed_density;

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
    this->stats_arr[PVSTAT].initialize("pv", 1, this->statsize, this->nens,
                                       this->masterproc);
    this->stats_arr[PESTAT].initialize("pens", 1, this->statsize, this->nens,
                                       this->masterproc);

    const auto &primal_topology = primal_geom.topology;
    const auto &dual_topology = dual_geom.topology;

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

      SArray<real, 1, VS::ndensity_prognostic> masslocal, massglobal;
      SArray<real, 1, VS::ndensity_prognostic> densmaxlocal, densmaxglobal;
      SArray<real, 1, VS::ndensity_prognostic> densminlocal, densminglobal;
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

      YAKL_SCOPE(Hs, equations->Hs);
      YAKL_SCOPE(Hk, equations->Hk);
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

      YAKL_SCOPE(PVPE, equations->PVPE);
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
                               REAL_MPI, MPI_SUM, 0, MPI_COMM_WORLD,
                               &this->Req[DENSSTAT]);
      this->ierr = MPI_Ireduce(&densmaxlocal, &densmaxglobal,
                               VS::ndensity_prognostic, REAL_MPI, MPI_MAX, 0,
                               MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
      this->ierr = MPI_Ireduce(&densminlocal, &densminglobal,
                               VS::ndensity_prognostic, REAL_MPI, MPI_MIN, 0,
                               MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
      this->ierr = MPI_Ireduce(&pvlocal, &pvglobal, 1, REAL_MPI, MPI_SUM, 0,
                               MPI_COMM_WORLD, &this->Req[PVSTAT]);
      this->ierr = MPI_Ireduce(&pelocal, &peglobal, 1, REAL_MPI, MPI_SUM, 0,
                               MPI_COMM_WORLD, &this->Req[PESTAT]);
      this->ierr = MPI_Ireduce(&elocal, &eglobal, 4, REAL_MPI, MPI_SUM, 0,
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

  using VS = VariableSet;

  // primal grid represents straight quantities, dual grid twisted quantities

  // v, dens
  prog_desc_arr[VVAR] = {"v", ptopo, 1, 1, 1}; // v = straight 1-form
  prog_desc_arr[DENSVAR] = {"dens", dtopo, ndims, 1,
                            VS::ndensity_prognostic}; // dens = twisted n-form

  // hs, coriolis
  const_desc_arr[HSVAR] = {"hs", dtopo, ndims, 1, 1}; // hs = twisted n-form
  const_desc_arr[CORIOLISVAR] = {"coriolis", ptopo, 2, 1,
                                 1}; // f = straight 2-form

  // functional derivatives = F, B, K, he, U
  aux_desc_arr[BVAR] = {"B", ptopo, 0, 1,
                        VS::ndensity_active};         // B = straight 0-form
  aux_desc_arr[FVAR] = {"F", dtopo, ndims - 1, 1, 1}; // F = twisted (n-1)-form
  aux_desc_arr[FVAR2] = {"F2", dtopo, ndims - 1, 1,
                         1};                      // F2 = twisted (n-1)-form
  aux_desc_arr[KVAR] = {"K", dtopo, ndims, 1, 1}; // K = twisted n-form
  aux_desc_arr[HEVAR] = {"he", dtopo, ndims - 1, 1,
                         1}; // he lives on dual edges, associated with F
  aux_desc_arr[UVAR] = {"U", dtopo, ndims - 1, 1, 1}; // U = twisted (n-1)-form

  // dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_desc_arr[DENS0VAR] = {"dens0", ptopo, 0, 1,
                            VS::ndensity_prognostic}; // dens0 = straight 0-form
  aux_desc_arr[DENSEDGERECONVAR] = {
      "densedgerecon", dtopo, ndims, 1,
      2 * ndims * VS::ndensity_prognostic}; // densedgerecon lives on dual
                                            // cells, associated with F
  aux_desc_arr[DENSRECONVAR] = {
      "densrecon", dtopo, ndims - 1, 1,
      VS::ndensity_prognostic}; // densrecon lives on dual edges, associated
                                // with F

  // dual grid reconstruction stuff- q0, f0, FT, qedgerecon, qrecon,
  // coriolisedgercon, coriolisrecon
  aux_desc_arr[Q0VAR] = {"q", dtopo, 0, 1, 1};  // q0 = twisted 0-form
  aux_desc_arr[F0VAR] = {"f", dtopo, 0, 1, 1};  // f0 = twisted 0-form
  aux_desc_arr[FTVAR] = {"FT", ptopo, 1, 1, 1}; // FT = straight 1-form
  aux_desc_arr[QEDGERECONVAR] = {"qedgerecon", ptopo, 2, 1,
                                 4}; // qedgerecon lives on primal cells
  aux_desc_arr[QRECONVAR] = {
      "qrecon", ptopo, 1, 1,
      1}; // qrecon lives on primal edges, associated with FT
  aux_desc_arr[CORIOLISEDGERECONVAR] = {
      "coriolisedgerecon", ptopo, 2, 1,
      4}; // coriolisedgerecon lives on primal cells
  aux_desc_arr[CORIOLISRECONVAR] = {
      "coriolisrecon", ptopo, 1, 1,
      1}; // coriolisrecon lives on primal edges, associated with FT

  // fct stuff- Mf, edgeflux
  aux_desc_arr[MFVAR] = {"Mf", dtopo, ndims, 1, VS::ndensity_prognostic};
  aux_desc_arr[EDGEFLUXVAR] = {"edgeflux", dtopo, ndims - 1, 1,
                               VS::ndensity_prognostic};
}

void testcase_from_string(std::unique_ptr<TestCase> &testcase,
                          std::string name);

void read_model_params_file(std::string inFile, ModelParameters &params,
                            Parallel &par, PamCoupler &coupler,
                            std::unique_ptr<TestCase> &testcase) {

  // Read common parameters
  int nz = coupler.get_nz();
  readParamsFile(inFile, params, par, nz);

  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  // Read the data initialization options
  params.initdataStr = config["initData"].as<std::string>();
  testcase_from_string(testcase, params.initdataStr);

  for (int i = 0; i < ntracers_dycore; i++) {
    params.tracerdataStr[i] =
        config["initTracer" + std::to_string(i)].as<std::string>("constant");
    params.dycore_tracerpos[i] =
        config["initTracerPos" + std::to_string(i)].as<bool>(false);
  }
}

void read_model_params_coupler(ModelParameters &params, Parallel &par,
                               PamCoupler &coupler,
                               std::unique_ptr<TestCase> &testcase) {}

void check_and_print_model_parameters(const ModelParameters &params,
                                      const Parallel &par) {

  check_and_print_parameters(params, par);

  serial_print("IC: " + params.initdataStr, par.masterproc);
  for (int i = 0; i < ntracers_dycore; i++) {
    serial_print("Dycore Tracer" + std::to_string(i) +
                     " IC: " + params.tracerdataStr[i],
                 par.masterproc);
  }
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

  std::array<real, 3> get_domain() const override { return {Lx, Ly, 1}; }

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

    YAKL_SCOPE(tracers, this->tracers);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_2form_values(
          YAKL_LAMBDA(real x, real y) {
            return h_f(x, y) *
                   TracerFunctor{}(tracers(i), x, y, Lx, Ly, xc, yc);
          },
          progvars.fields_arr[DENSVAR],
          i + VariableSet::ndensity_dycore_prognostic);
    }
    equations->Hs.set_parameters(g);
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = this->equations->reference_state;
    refstate.ref_height = T::ref_height;
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
  static real constexpr ref_height = H0;

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

  static VecXY YAKL_INLINE v_f(real x, real y) {
    VecXY vvec;

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

// real YAKL_INLINE bickley_tracer(real x, real y) {
//
//   return std::sin(bickley_constants.k * y);
// }

struct BickleyJet {
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 4 * pi;
  static real constexpr Ly = 4 * pi;
  static real constexpr eps = 0.1_fp;
  static real constexpr l = 0.5_fp;
  static real constexpr k = 0.5_fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr yc = 0.5_fp * Ly;
  static real constexpr ref_height = 1._fp;

  static real YAKL_INLINE coriolis_f(real x, real y) { return 0; }

  static real YAKL_INLINE h_f(real x, real y) { return 1; }

  static VecXY YAKL_INLINE v_f(real x, real y) {
    VecXY vvec;
    real U = std::pow(std::cosh(y), -2);
    real psi = std::exp(-std::pow(y + l / 10, 2) / (2 * l * l)) *
               std::cos(k * x) * std::cos(k * y);

    real u = psi * (k * std::tan(k * y) + y / (l * l));
    real v = -psi * k * std::tan(k * x);

    vvec.u = U + eps * u;
    vvec.v = eps * v;

    return vvec;
  }

  static real YAKL_INLINE S_f(real x, real y) { return g * h_f(x, y); }
};

void testcase_from_string(std::unique_ptr<TestCase> &testcase,
                          std::string name) {
  if (name == "doublevortex") {
    testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  } else if (name == "bickleyjet") {
    testcase = std::make_unique<SWETestCase<BickleyJet>>();
  } else {
    throw std::runtime_error("unknown test case");
  }
}

void testcase_from_config(std::unique_ptr<TestCase> &testcase,
                          const YAML::Node &config) {
  const std::string name = config["initData"].as<std::string>();
  testcase_from_string(testcase, name);
}
