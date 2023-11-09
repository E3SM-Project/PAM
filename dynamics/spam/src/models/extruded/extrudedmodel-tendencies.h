#pragma once

namespace pamc {

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
  real scalar_diffusion_coeff;
  real scalar_diffusion_subtract_refstate;
  real velocity_diffusion_coeff;
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
    scalar_diffusion_coeff = params.scalar_diffusion_coeff;
    scalar_diffusion_subtract_refstate =
        params.scalar_diffusion_subtract_refstate;
    velocity_diffusion_coeff = params.velocity_diffusion_coeff;
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
          "Compute FTvar",
          SimpleBounds<4>(primal_topology.ni - 2, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_Wyz_u(FTvar, FWvar, pis, pjs, pks, i, j, k + 1, n);
          });
      parallel_for(
          "Compute FTvar bnd",
          SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                          primal_topology.nens),
          YAKL_LAMBDA(int j, int i, int n) {
            compute_Wyz_u_bottom(FTvar, FWvar, pis, pjs, pks, i, j, 0, n);
            compute_Wyz_u_top(FTvar, FWvar, pis, pjs, pks, i, j,
                              primal_topology.ni - 1, n);
          });

      parallel_for(
          "Compute FTWvar",
          SimpleBounds<4>(primal_topology.nl - 2, primal_topology.n_cells_y,
                          primal_topology.n_cells_x, primal_topology.nens),
          YAKL_LAMBDA(int k, int j, int i, int n) {
            compute_Wyz_w(FTWvar, Fvar, pis, pjs, pks, i, j, k + 1, n);
          });
      parallel_for(
          "Compute FTWvar bnd",
          SimpleBounds<3>(primal_topology.n_cells_y, primal_topology.n_cells_x,
                          primal_topology.nens),
          YAKL_LAMBDA(int j, int i, int n) {
            compute_Wyz_w_bottom(FTWvar, Fvar, pis, pjs, pks, i, j, 0, n);
            compute_Wyz_w_top(FTWvar, Fvar, pis, pjs, pks, i, j,
                              primal_topology.nl - 1, n);
          });

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

  void add_scalar_diffusion(real scalar_diffusion_coeff, real5d denstendvar,
                            const real5d densvar, const real5d dens0var,
                            const real5d Fdiffvar, const real5d FWdiffvar,
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
              dens0 -= q_di(dens_id, k + dks, n);
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
          compute_Dnm1bar<1>(hdiv, Fdiffvar, dis, djs, dks, i, j, k, n);
          SArray<real, 1, VS::ndensity_diffused> vdiv;
          compute_Dnm1bar_vert<1>(vdiv, FWdiffvar, dis, djs, dks, i, j, k, n);

          const real rho =
              varset.get_total_density(densvar, k, j, i, pks, pjs, pis, n);
          const real Hn1bar_diag = Hn1bar_diagonal(
              primal_geometry, dual_geometry, pis, pjs, pks, i, j, k, n);

          for (int d = 0; d < VS::ndensity_diffused; ++d) {
            const real diff_tend = -scalar_diffusion_coeff * rho * Hn1bar_diag *
                                   (hdiv(d) + vdiv(d));
            int dens_id = varset.diffused_dens_ids(d);
            denstendvar(dens_id, k + pks, j + pjs, i + pis, n) += diff_tend;
          }
        });

    yakl::timer_stop("add_scalar_diffusion");
  }

  void add_velocity_diffusion_2d(
      real velocity_coeff, real5d Vtendvar, real5d Wtendvar, const real5d Vvar,
      const real5d Wvar, const real5d qhzedgereconvar, const real5d qhzvar,
      const real5d dens0var, const real5d Kvar, const real5d Fvar,
      const real5d FWvar, FieldSet<nauxiliary> &auxiliary_vars) {
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

    // *d
    parallel_for(
        "Velocity diffusion 1",
        SimpleBounds<4>(dual_topology.ni - 2, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          SArray<real, 2, 1, ndims> Dv;
          compute_D1_ext<1>(Dv, Vvar, Wvar, pis, pjs, pks, i, j, k, n);
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
          compute_H10<1, diffusion_diff_ord>(Fvar, Vvar, primal_geometry,
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
          // *d(*d)
          const real Dvert =
              compute_D0bar_vert<1>(qhzvar, dis, djs, dks, i, j, k, n);
          SArray<real, 1, ndims> Hnm11bar_diag;
          Hnm11bar_diagonal(Hnm11bar_diag, primal_geometry, dual_geometry, dis,
                            djs, dks, i, j, k, n);
          Vtendvar(0, k + pks, j + pjs, i + pis, n) -=
              velocity_coeff * Dvert * Hnm11bar_diag(0);

          // d(*d*)
          SArray<real, 2, 1, ndims> vdiff;
          compute_D0<1>(vdiff, dens0var, pis, pjs, pks, i, j, k, n);
          for (int d = 0; d < ndims; ++d) {
            Vtendvar(d, pks + k, pjs + j, pis + i, n) -=
                velocity_coeff * vdiff(0, d);
          }
        });

    parallel_for(
        "Velocity diffusion 6",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          // *d(*d)
          SArray<real, 1, ndims> Dhorz;
          compute_D0bar_ext<1>(Dhorz, qhzvar, dis, djs, dks, i, j, k + 1, n);
          const real Hn0bar_diag = Hn0bar_diagonal(
              primal_geometry, dual_geometry, dis, djs, dks, i, j, k, n);
          Wtendvar(0, k + pks, j + pjs, i + pis, n) -=
              velocity_coeff * Dhorz(0) * Hn0bar_diag;

          // d(*d*)
          SArray<real, 1, 1> wdiff;
          compute_D0_vert<1>(wdiff, dens0var, pis, pjs, pks, i, j, k, n);
          Wtendvar(0, pks + k, pjs + j, pis + i, n) -=
              velocity_coeff * wdiff(0);
        });

    yakl::timer_stop("add_velocity_diffusion");
  }

  void add_velocity_diffusion_3d(
      real velocity_coeff, real5d Vtendvar, real5d Wtendvar, const real5d Vvar,
      const real5d Wvar, const real5d qhzedgereconvar, const real5d qhzvar,
      const real5d dens0var, const real5d Kvar, const real5d Fvar,
      const real5d FWvar, const real5d qxyedgereconvar, const real5d qxyvar,
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

    // *d*d
    parallel_for(
        "Velocity diffusion 1",
        SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D1_ext<1>(qhzedgereconvar, Vvar, Wvar, pis, pjs, pks, i, j, k,
                            n);
        });
    auxiliary_vars.exchange({QHZEDGERECONVAR});

    parallel_for(
        "Velocity diffusion 1",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_D1<1>(qxyedgereconvar, Vvar, pis, pjs, pks, i, j, k, n);
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
                velocity_coeff * vdiff(d);
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
              velocity_coeff * wdiff(0);
        });

    // d*d*
    parallel_for(
        "Velocity diffusion 7",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          compute_H10<1, diffusion_diff_ord>(Fvar, Vvar, primal_geometry,
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
                velocity_coeff * vdiff(0, d);
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
              velocity_coeff * wdiff(0);
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
                             : auxiliary_vars.fields_arr[F2VAR].data,
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

    if (scalar_diffusion_coeff > 0) {
      add_scalar_diffusion(
          scalar_diffusion_coeff, xtend.fields_arr[DENSVAR].data,
          x.fields_arr[DENSVAR].data, auxiliary_vars.fields_arr[DENS0VAR].data,
          auxiliary_vars.fields_arr[FDIFFVAR].data,
          auxiliary_vars.fields_arr[FWDIFFVAR].data, auxiliary_vars);
    }
    if (velocity_diffusion_coeff > 0) {
      if (ndims == 1) {
        add_velocity_diffusion_2d(
            velocity_diffusion_coeff, xtend.fields_arr[VVAR].data,
            xtend.fields_arr[WVAR].data, x.fields_arr[VVAR].data,
            x.fields_arr[WVAR].data,
            auxiliary_vars.fields_arr[QHZEDGERECONVAR].data,
            auxiliary_vars.fields_arr[QHZVAR].data,
            auxiliary_vars.fields_arr[DENS0VAR].data,
            auxiliary_vars.fields_arr[KVAR].data,
            auxiliary_vars.fields_arr[FVAR].data,
            auxiliary_vars.fields_arr[FWVAR].data, auxiliary_vars);
      } else {
        add_velocity_diffusion_3d(
            velocity_diffusion_coeff, xtend.fields_arr[VVAR].data,
            xtend.fields_arr[WVAR].data, x.fields_arr[VVAR].data,
            x.fields_arr[WVAR].data,
            auxiliary_vars.fields_arr[QHZEDGERECONVAR].data,
            auxiliary_vars.fields_arr[QHZVAR].data,
            auxiliary_vars.fields_arr[DENS0VAR].data,
            auxiliary_vars.fields_arr[KVAR].data,
            auxiliary_vars.fields_arr[FVAR].data,
            auxiliary_vars.fields_arr[FWVAR].data,
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

  void remove_negative_densities(FieldSet<nprognostic> &x) override {
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

} // namespace pamc
