#pragma once

namespace pamc {

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

} // namespace pamc
