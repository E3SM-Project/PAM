#pragma once

namespace pamc {

struct CompressiblePressureLinearSystem : PressureLinearSystem {
  real3d linp_coeff;

  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  Equations &equations) override {
    PressureLinearSystem::initialize(params, primal_geom, dual_geom, equations);
    linp_coeff = real3d("linp coeff", VS::ndensity, primal_geometry.topology.ni,
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
          SArray<real, 1, VS::ndensity> p_coeff;
          varset.linear_pressure_coeffs(p_coeff, k, pks, n);
          for (int d = 0; d < VS::ndensity; ++d) {
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

          for (int dd = 0; dd < VS::ndensity_prognostic; ++dd) {
            for (int d = 0; d < ndims; ++d) {
              tri_d(k, j, i, n) -= alpha * alpha * linp_coeff(dd, k, n) *
                                   fHn1bar * fH1(d) * fD0Dbar(d) *
                                   q_pi(dd, k + pks, n);
            }
          }

          for (int d = 0; d < VS::ndensity_prognostic; ++d) {
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
          SArray<real, 1, VS::ndensity_prognostic> tend;
          SArray<real, 1, VS::ndensity_prognostic> vtend;

          compute_wDnm1bar<VS::ndensity_prognostic>(tend, q_pi, Fvar, dis, djs,
                                                    dks, i, j, k, n);
          compute_wDnm1bar_vert<VS::ndensity_prognostic>(
              vtend, q_di, FWvar, dis, djs, dks, i, j, k, n);

          for (int d = 0; d < VS::ndensity_prognostic; ++d) {
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
          compute_Hn1bar<VS::ndensity_prognostic, diff_ord, vert_diff_ord>(
              Bvar, mfvar, primal_geometry, dual_geometry, pis, pjs, pks, i, j,
              k, n);
        });

    parallel_for(
        "Linear solve rhs 4",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_LAMBDA(int k, int j, int i, int n) {
          p_transform(k, j, i, n) = 0;
          for (int d = 0; d < VS::ndensity_prognostic; ++d) {
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
} // namespace pamc
