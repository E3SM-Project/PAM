#pragma once

namespace pamc {

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
      throw std::runtime_error("semi-implicit not implemented in 3d yet");
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
} // namespace pamc
