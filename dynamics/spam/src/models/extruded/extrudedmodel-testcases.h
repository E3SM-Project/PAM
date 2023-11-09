#pragma once

namespace pamc {

// Universal

real YAKL_INLINE isentropic_T(real z, real theta0, real g,
                              const ThermoPotential &thermo) {
  return theta0 - z * g / thermo.cst.Cpd;
}

real YAKL_INLINE isentropic_p(real z, real theta0, real g,
                              const ThermoPotential &thermo) {
  return thermo.cst.pr * pow(isentropic_T(z, theta0, g, thermo) / theta0,
                             1. / thermo.cst.kappa_d);
}

real YAKL_INLINE isentropic_rho(real z, real theta0, real g,
                                const ThermoPotential &thermo) {
  real p = isentropic_p(z, theta0, g, thermo);
  real T = isentropic_T(z, theta0, g, thermo);
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

real YAKL_INLINE flat_geop(real z, real g) { return g * z; }

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
  using T::Lz;
  using T::xc;
  using T::yc;

  using T::h_f;
#ifdef PAMC_TSWE
  using T::S_f;
#endif
  using T::coriolis_f;
  using T::v_f;

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
    params.ylen = T::Ly;
    params.yc = T::yc;
  }

  std::array<real, 3> get_domain() const override { return {Lx, Ly, Lz}; }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return h_f(x, y, z); },
        progvars.fields_arr[DENSVAR], 0);
#ifdef PAMC_TSWE
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return S_f(x, y, z); },
        progvars.fields_arr[DENSVAR], 1);
#endif
    primal_geom.set_10form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[VVAR], 0);
    primal_geom.set_01form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[WVAR], 0);

    if (ndims == 1) {
      primal_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) { return coriolis_f(x, y, z).v; },
          constvars.fields_arr[CORIOLISHZVAR], 0);
    } else {
      primal_geom.set_nm11form_values(
          YAKL_LAMBDA(real x, real y, real z) { return coriolis_f(x, y, z); },
          constvars.fields_arr[CORIOLISHZVAR], 0);

      primal_geom.set_n0form_values(
          YAKL_LAMBDA(real x, real y, real z) { return coriolis_f(x, y, z); },
          constvars.fields_arr[CORIOLISXYVAR], 0);
    }

    // YAKL_SCOPE(tracer_f, this->tracer_f);
    // for (int i = 0; i < ntracers_dycore; i++) {
    //   dual_geom.set_11form_values(
    //       YAKL_LAMBDA(real x, real y) {
    //         return h_f(x, y) * tracer_f(i)->compute(x, y, Lx, Ly, xc, yc);
    //       },
    //       progvars.fields_arr[DENSVAR], i + ndensity_dycore);
    // }
    this->equations->Hs.set_parameters(g);
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = this->equations->reference_state;

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    const int pks = primal_topology.ks;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(varset, equations->varset);
    YAKL_SCOPE(Hs, equations->Hs);

    // real href = T::H0;
    real href = 0;

    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return 0; }, refstate.geop, 0);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return href; }, refstate.dens,
        varset.dens_id_mass);

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return href; }, refstate.rho_pi,
        varset.dens_id_mass);
    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return 1; }, refstate.q_pi, varset.dens_id_mass);

    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return href; }, refstate.rho_di,
        varset.dens_id_mass);
    dual_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return 1; }, refstate.q_di, varset.dens_id_mass);
  }
};

template <class T> class EulerTestCase : public TestCase, public T {
  using VS = VariableSet;

public:
  using T::g;
  using T::Lx;
  using T::Lz;
  using T::xc;
  using T::zc;

  using T::entropicvar_f;

  using T::refentropicdensity_f;
  using T::refnsq_f;
  using T::refrho_f;
  using T::v_f;

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
#ifdef PAMC_AN
    return refrho_f(z, thermo);
#else
    return T::rho_f(x, y, z, thermo);
#endif
  }

  std::array<real, 3> get_domain() const override {
    real Ly = 1;
    if constexpr (T::max_ndims > 1) {
      Ly = T::Ly;
    }
    return {Lx, Ly, Lz};
  }

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
    if constexpr (T::max_ndims > 1) {
      params.ylen = T::Ly;
      params.yc = T::yc;
    }
  }

  void add_diagnostics(
      std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {
    T::add_diagnostics(diagnostics);
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    equations->Hs.set_parameters(g);

    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(varset, equations->varset);
#ifndef PAMC_AN
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return rho_f(x, y, z, thermo); },
        progvars.fields_arr[DENSVAR], varset.dens_id_mass);
#endif
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) {
          return rho_f(x, y, z, thermo) * entropicvar_f(x, y, z, thermo);
        },
        progvars.fields_arr[DENSVAR], varset.dens_id_entr);

    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return flat_geop(z, g); },
        constvars.fields_arr[HSVAR], 0);

    primal_geom.set_10form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[VVAR], 0);
    primal_geom.set_01form_values(
        YAKL_LAMBDA(real x, real y, real z) { return v_f(x, y, z); },
        progvars.fields_arr[WVAR], 0);

    YAKL_SCOPE(tracers, this->tracers);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return rho_f(x, y, z, thermo) *
                   TracerFunctor{}(tracers(i), x, z, Lx, Lz, xc, zc);
          },
          progvars.fields_arr[DENSVAR], i + VS::ndensity_dycore_prognostic);
    }
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = this->equations->reference_state;

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    const int pks = primal_topology.ks;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(varset, equations->varset);
    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(Hs, equations->Hs);

    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return flat_geop(z, g); }, refstate.geop, 0);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.dens,
        varset.dens_id_mass);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.dens, varset.dens_id_entr);

    parallel_for(
        "compute rho_pi and unscaled q_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          SArray<real, 1, VS::ndensity> dens0;
          compute_Hn1bar<VS::ndensity, vert_diff_ord>(
              dens0, refstate.dens.data, primal_geom, dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k + pks, n) = dens0(varset.dens_id_mass);
          for (int d = 0; d < VS::ndensity; ++d) {
            refstate.q_pi.data(d, k + pks, n) = dens0(d);
          }
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refnsq_f(z, thermo); }, refstate.Nsq_pi,
        0);

    parallel_for(
        "compute rho_di and unscaled q_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, dks, n) = refstate.rho_pi.data(0, pks, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1 + dks, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1 + pks, n);
          } else {
            refstate.rho_di.data(0, k + dks, n) =
                0.5_fp * (refstate.rho_pi.data(0, k + pks, n) +
                          refstate.rho_pi.data(0, k - 1 + pks, n));
          }

          for (int d = 0; d < VS::ndensity; ++d) {
            if (k == 0) {
              refstate.q_di.data(d, dks, n) = refstate.q_pi.data(d, pks, n);
            } else if (k == dual_topology.ni - 1) {
              refstate.q_di.data(d, dual_topology.ni - 1 + dks, n) =
                  refstate.q_pi.data(d, primal_topology.ni - 1 + pks, n);
            } else {
              refstate.q_di.data(d, k + dks, n) =
                  0.5_fp * (refstate.q_pi.data(d, k + pks, n) +
                            refstate.q_pi.data(d, k - 1 + pks, n));
            }
          }
        });

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_pi.data(l, k + pks, n) /=
                refstate.rho_pi.data(0, k + pks, n);
          }
        });

    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_di.data(l, k + dks, n) /=
                refstate.rho_di.data(0, k + dks, n);
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
        "Compute refstate pres_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_pi.data(0, k + pks, n);
          const real entropicvar =
              refstate.q_pi.data(varset.dens_id_entr, k + pks, n);
          refstate.pres_pi.data(0, k + pks, n) =
              thermo.solve_p(rho, entropicvar, 1, 0, 0, 0);
        });

    parallel_for(
        "Compute refstate pres_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_di.data(0, k + dks, n);
          const real entropicvar =
              refstate.q_di.data(varset.dens_id_entr, k + dks, n);
          refstate.pres_di.data(0, k + dks, n) =
              thermo.solve_p(rho, entropicvar, 1, 0, 0, 0);
        });
  }
};

template <class T> class MoistEulerTestCase : public TestCase, public T {
  using VS = VariableSet;

public:
  using T::g;
  using T::Lx;
  using T::Lz;
  using T::xc;
  using T::zc;

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
#ifdef PAMC_MAN
    return T::refrho_f(z, thermo);
#elif defined PAMC_MCErhod || defined PAMC_MCErhodp
    return T::rhod_f(x, y, z, thermo);
#else
    return T::rho_f(x, y, z, thermo);
#endif
  }
  using T::entropicvar_f;
  using T::refentropicdensity_f;
  using T::refnsq_f;
  using T::refrho_f;
  using T::refrhov_f;

  std::array<real, 3> get_domain() const override {
    real Ly = 1;
    if constexpr (T::max_ndims > 1) {
      Ly = T::Ly;
    }
    return {Lx, Ly, Lz};
  }

  void set_domain(ModelParameters &params) override {
    params.xlen = Lx;
    params.xc = xc;
    if constexpr (T::max_ndims > 1) {
      params.ylen = T::Ly;
      params.yc = T::yc;
    }
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    equations->Hs.set_parameters(g);

    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(varset, equations->varset);

#ifndef PAMC_MAN
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return rho_f(x, y, z, thermo); },
        progvars.fields_arr[DENSVAR], varset.dens_id_mass);
#endif
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) {
          return rho_f(x, y, z, thermo) * entropicvar_f(x, y, z, thermo);
        },
        progvars.fields_arr[DENSVAR], varset.dens_id_entr);

    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return flat_geop(z, g); },
        constvars.fields_arr[HSVAR], 0);

    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) {
          return T::rhov_f(x, y, z, thermo);
        },
        progvars.fields_arr[DENSVAR], varset.dens_id_vap);

    YAKL_SCOPE(tracers, this->tracers);
    for (int i = 0; i < ntracers_dycore; i++) {
      dual_geom.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return rho_f(x, y, z, thermo) *
                   TracerFunctor{}(tracers(i), x, z, Lx, Lz, xc, zc);
          },
          progvars.fields_arr[DENSVAR], i + VS::ndensity_dycore_prognostic);
    }
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {
    auto &refstate = this->equations->reference_state;

    const auto primal_topology = primal_geom.topology;
    const auto dual_topology = dual_geom.topology;

    const int pks = primal_topology.ks;
    const int dks = dual_topology.ks;

    YAKL_SCOPE(thermo, equations->thermo);
    YAKL_SCOPE(Hs, equations->Hs);
    YAKL_SCOPE(varset, equations->varset);

    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return flat_geop(z, g); }, refstate.geop, 0);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return refrho_f(z, thermo); }, refstate.dens,
        varset.dens_id_mass);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return refentropicdensity_f(z, thermo); },
        refstate.dens, varset.dens_id_entr);
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return refrhov_f(z, thermo); }, refstate.dens,
        varset.dens_id_vap);

    parallel_for(
        "compute rho_pi and unscaled q_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const auto total_density_f = TotalDensityFunctor{varset};
          SArray<real, 1, 1> rho0;
          compute_Hn1bar<1, vert_diff_ord>(total_density_f, rho0,
                                           refstate.dens.data, primal_geom,
                                           dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k + pks, n) = rho0(0);

          SArray<real, 1, VS::ndensity> dens0;
          compute_Hn1bar<VS::ndensity, vert_diff_ord>(
              dens0, refstate.dens.data, primal_geom, dual_geom, pks, k, n);
          for (int d = 0; d < VS::ndensity; ++d) {
            refstate.q_pi.data(d, k + pks, n) = dens0(d);
          }
        });

    primal_geom.set_profile_00form_values(
        YAKL_LAMBDA(real z) { return refnsq_f(z, thermo); }, refstate.Nsq_pi,
        0);

    parallel_for(
        "compute rho_di and unscaled q_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, dks, n) = refstate.rho_pi.data(0, pks, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1 + dks, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1 + pks, n);
          } else {
            refstate.rho_di.data(0, k + dks, n) =
                0.5_fp * (refstate.rho_pi.data(0, k + pks, n) +
                          refstate.rho_pi.data(0, k - 1 + pks, n));
          }

          for (int d = 0; d < VS::ndensity; ++d) {
            if (k == 0) {
              refstate.q_di.data(d, dks, n) = refstate.q_pi.data(d, pks, n);
            } else if (k == dual_topology.ni - 1) {
              refstate.q_di.data(d, dual_topology.ni - 1 + dks, n) =
                  refstate.q_pi.data(d, primal_topology.ni - 1 + pks, n);
            } else {
              refstate.q_di.data(d, k + dks, n) =
                  0.5_fp * (refstate.q_pi.data(d, k + pks, n) +
                            refstate.q_pi.data(d, k - 1 + pks, n));
            }
          }
        });

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_pi.data(l, k + pks, n) /=
                refstate.rho_pi.data(0, k + pks, n);
          }
        });
    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_di.data(l, k + dks, n) /=
                refstate.rho_di.data(0, k + dks, n);
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
        "Compute refstate pres_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_pi.data(0, k + pks, n);
          const real entropicvar =
              refstate.q_pi.data(varset.dens_id_entr, k + pks, n);
          const real qv = refstate.q_pi.data(varset.dens_id_vap, k + pks, n);
          refstate.pres_pi.data(0, k + pks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });

    parallel_for(
        "Compute refstate pres_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_di.data(0, k + dks, n);
          const real entropicvar =
              refstate.q_di.data(varset.dens_id_entr, k + dks, n);
          const real qv = refstate.q_di.data(varset.dens_id_vap, k + dks, n);
          refstate.pres_di.data(0, k + dks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });
  }

  void add_diagnostics(
      std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {
    T::add_diagnostics(diagnostics);
  }
};

class CoupledTestCase : public TestCase {
  using VS = VariableSet;

public:
  PamCoupler &coupler;
  bool set_coupler_state = false;

  CoupledTestCase(PamCoupler &coupler) : coupler(coupler) {}

  std::array<real, 3> get_domain() const override { return {0, 1, 0}; }

  void set_domain(ModelParameters &params) override {
    params.xlen = coupler.get_xlen();
    params.xc = 0;
  }

  void set_initial_conditions(FieldSet<nprognostic> &progvars,
                              FieldSet<nconstant> &constvars,
                              const Geometry<Straight> &primal_geom,
                              const Geometry<Twisted> &dual_geom) override {

    const real g = coupler.get_option<real>("grav");
    equations->Hs.set_parameters(g);
    auto &varset = this->equations->varset;
    dual_geom.set_n1form_values(
        YAKL_LAMBDA(real x, real y, real z) { return flat_geop(z, g); },
        constvars.fields_arr[HSVAR], 0);

    convert_coupler_to_dynamics_densities(varset, coupler, progvars, constvars);
    convert_coupler_to_dynamics_wind(varset, coupler, progvars, constvars,
                                     false);
  }

  void set_reference_state(const Geometry<Straight> &primal_geom,
                           const Geometry<Twisted> &dual_geom) override {

    auto &refstate = this->equations->reference_state;
    const auto &varset = this->equations->varset;
    auto &thermo = this->equations->thermo;
    const auto &primal_topology = primal_geom.topology;
    const auto &dual_topology = dual_geom.topology;

    // Set thermo constants based on coupler values
    // Need a better way to separate fundamental vs derived constants
    // Also I think only P3 defines all of these

    thermo.cst.Rd = coupler.get_option<real>("R_d");
    thermo.cst.Rv = coupler.get_option<real>("R_v");
    thermo.cst.pr = coupler.get_option<real>("p0");
    thermo.cst.Cpd = coupler.get_option<real>("cp_d");
    thermo.cst.Cvd = thermo.cst.Cpd - thermo.cst.Rd;
    thermo.cst.Cpv = coupler.get_option<real>("cp_v");

    thermo.cst.Lvr = coupler.get_option<real>("latvap");
    thermo.cst.Lfr = coupler.get_option<real>("latice");

    thermo.cst.gamma_d = thermo.cst.Cpd / thermo.cst.Cvd;
    thermo.cst.kappa_d = thermo.cst.Rd / thermo.cst.Cpd;
    thermo.cst.delta_d = thermo.cst.Rd / thermo.cst.Cvd;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    const int pis = primal_topology.is;
    const int pjs = primal_topology.js;
    const int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readonly();

    auto dm_ref_dens_dry = dm.get<real const, 2>("ref_density_dry");
    auto dm_ref_dens_vap = dm.get<real const, 2>("ref_density_vapor");
    auto dm_ref_dens_liq = dm.get<real const, 2>("ref_density_liq");
    auto dm_ref_dens_ice = dm.get<real const, 2>("ref_density_ice");
    auto dm_ref_temp = dm.get<real const, 2>("ref_temp");

    const real grav = coupler.get_option<real>("grav");
    dual_geom.set_profile_n1form_values(
        YAKL_LAMBDA(real z) { return flat_geop(z, grav); }, refstate.geop, 0);

    // sets dens and unscaled q_pi
    parallel_for(
        "Coupled reference state 1",
        SimpleBounds<2>(dual_topology.nl, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int n) {
          const real temp = dm_ref_temp(k, n);
          const real dens_dry = dm_ref_dens_dry(k, n);
          const real dens_vap = dm_ref_dens_vap(k, n);
          const real dens_liq = dm_ref_dens_liq(k, n);
          const real dens_ice = dm_ref_dens_ice(k, n);
          const real dens = dens_dry + dens_vap;

          const real qd = dens_dry / dens;
          const real qv = dens_vap / dens;
          const real ql = dens_liq / dens;
          const real qi = dens_ice / dens;

          const real alpha = 1.0_fp / dens;
          const real entropic_var = thermo.compute_entropic_var_from_alpha_T(
              alpha, temp, qd, qv, ql, qi);

          const real dual_volume =
              dual_geom.get_area_n1entity(k + dks, djs, dis, n);
          refstate.dens.data(varset.dens_id_mass, k + dks, n) =
              dens * dual_volume;
          refstate.dens.data(varset.dens_id_entr, k + dks, n) =
              entropic_var * dens * dual_volume;
          refstate.dens.data(varset.dens_id_vap, k + dks, n) =
              dens_vap * dual_volume;

          refstate.q_pi.data(varset.dens_id_mass, k + dks, n) = dens;
          refstate.q_pi.data(varset.dens_id_entr, k + dks, n) =
              dens * entropic_var;
          refstate.q_pi.data(varset.dens_id_vap, k + dks, n) = dens_vap;
        });

    parallel_for(
        "compute unscaled q_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int d = 0; d < VS::ndensity; ++d) {
            if (k == 0) {
              refstate.q_di.data(d, dks, n) = refstate.q_pi.data(d, pks, n);
            } else if (k == dual_topology.ni - 1) {
              refstate.q_di.data(d, dual_topology.ni - 1 + dks, n) =
                  refstate.q_pi.data(d, primal_topology.ni - 1 + pks, n);
            } else {
              real q_km1 = refstate.q_pi.data(d, k - 1 + pks, n);
              real q_k = refstate.q_pi.data(d, k + pks, n);

              refstate.q_di.data(d, k + dks, n) =
                  q_km1 + (q_k - q_km1) *
                              (dual_geom.zint(k + pks, n) -
                               primal_geom.zint(k - 1 + pks, n)) /
                              primal_geom.dz(k - 1 + pks, n);
            }
          }
        });

    parallel_for(
        "compute rho_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const auto total_density_f = TotalDensityFunctor{varset};
          SArray<real, 1, 1> rho0;
          compute_Hn1bar<1, vert_diff_ord>(total_density_f, rho0,
                                           refstate.dens.data, primal_geom,
                                           dual_geom, pks, k, n);
          refstate.rho_pi.data(0, k + pks, n) = rho0(0);
        });

    parallel_for(
        "scale q_pi", SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_pi.data(l, k + pks, n) /=
                refstate.rho_pi.data(0, k + pks, n);
          }
#if defined PAMC_AN || defined PAMC_MAN
          refstate.q_pi.data(varset.dens_id_mass, k + pks, n) = 1;
#endif
        });

    parallel_for(
        "compute rho_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            refstate.rho_di.data(0, dks, n) = refstate.rho_pi.data(0, pks, n);
          } else if (k == dual_topology.ni - 1) {
            refstate.rho_di.data(0, dual_topology.ni - 1 + dks, n) =
                refstate.rho_pi.data(0, primal_topology.ni - 1 + pks, n);
          } else {
            real rho_k = refstate.rho_pi.data(0, k + pks, n);
            real rho_km1 = refstate.rho_pi.data(0, k - 1 + pks, n);

            refstate.rho_di.data(0, k + dks, n) =
                rho_km1 + (rho_k - rho_km1) *
                              (dual_geom.zint(k + pks, n) -
                               primal_geom.zint(k - 1 + pks, n)) /
                              primal_geom.dz(k - 1 + pks, n);
          }
        });
    parallel_for(
        "scale q_di", SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          for (int l = 0; l < VS::ndensity; ++l) {
            refstate.q_di.data(l, k + dks, n) /=
                refstate.rho_di.data(0, k + dks, n);
          }
#if defined PAMC_AN || defined PAMC_MAN
          refstate.q_di.data(varset.dens_id_mass, k + pks, n) = 1;
#endif
        });

    YAKL_SCOPE(Hs, equations->Hs);

    parallel_for(
        "Compute refstate B",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          Hs.compute_dHsdx(refstate.B.data, refstate.dens.data,
                           refstate.geop.data, pks, k, n, -1);
        });

    parallel_for(
        "compute Nsq",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real Rd = thermo.cst.Rd;
          const real Rv = thermo.cst.Rv;
          const real Cpd = thermo.cst.Cpd;
          const real Lvr = thermo.cst.Lvr;
          const real eta = Rv / Rd;

          real T_kp, T_km, T;
          real rv_kp, rv_km, rv;
          real dz;

          if (k == 0) {
            T_km = dm_ref_temp(k, n);
            T_kp = dm_ref_temp(k + 1, n);
            T = T_km;

            rv_kp = dm_ref_dens_vap(k + 1, n) / dm_ref_dens_dry(k + 1, n);
            rv_km = dm_ref_dens_vap(k, n) / dm_ref_dens_dry(k, n);
            rv = rv_km;

            dz = primal_geom.dz(k + pks, n);
          } else if (k == primal_topology.ni - 1) {
            T_km = dm_ref_temp(k - 1, n);
            T_kp = dm_ref_temp(k, n);
            T = dm_ref_temp(k, n);

            rv_kp = dm_ref_dens_vap(k, n) / dm_ref_dens_dry(k, n);
            rv_km = dm_ref_dens_vap(k - 1, n) / dm_ref_dens_dry(k - 1, n);
            rv = rv_kp;

            dz = primal_geom.dz(k - 1 + pks, n);
          } else {
            T = dm_ref_temp(k, n);
            T_km = dm_ref_temp(k - 1, n);
            T_kp = dm_ref_temp(k + 1, n);

            rv_km = dm_ref_dens_vap(k - 1, n) / dm_ref_dens_dry(k - 1, n);
            rv_kp = dm_ref_dens_vap(k + 1, n) / dm_ref_dens_dry(k + 1, n);
            rv = dm_ref_dens_vap(k, n) / dm_ref_dens_dry(k, n);

            dz = primal_geom.dz(k + pks, n) + primal_geom.dz(k - 1 + pks, n);
          }
          real dTdz = (T_kp - T_km) / dz;
          real drvdz = (rv_kp - rv_km) / dz;

          real Tv = T * (1 + eta * rv) / (1 + rv);
          real es = saturation_vapor_pressure(T);
          real rsw = (es / (Rd * T) - 1) * Rd / Rv;
          real qsw = rsw / (1 + rsw);

          real D1w = 1 + (1 + eta * rsw) * Lvr * qsw / (Rd * Tv);
          real D2w = 1 + (1 + eta * rsw) * Lvr * Lvr * qsw / (Cpd * Rv * T * T);
          real gamma_m = grav / Cpd * D1w / D2w;

          refstate.Nsq_pi.data(0, k + pks, n) =
              grav / T * D1w * (dTdz + gamma_m) - grav / (1 + rv) * drvdz;
        });

    parallel_for(
        "Compute refstate pres_pi",
        SimpleBounds<2>(primal_topology.ni, primal_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_pi.data(0, k + pks, n);
          const real entropicvar =
              refstate.q_pi.data(varset.dens_id_entr, k + pks, n);
          const real qv = refstate.q_pi.data(varset.dens_id_vap, k + pks, n);
          refstate.pres_pi.data(0, k + pks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });

    parallel_for(
        "Compute refstate pres_di",
        SimpleBounds<2>(dual_topology.ni, dual_topology.nens),
        YAKL_LAMBDA(int k, int n) {
          const real rho = refstate.rho_di.data(0, k + dks, n);
          const real entropicvar =
              refstate.q_di.data(varset.dens_id_entr, k + dks, n);
          const real qv = refstate.q_di.data(varset.dens_id_vap, k + dks, n);
          refstate.pres_di.data(0, k + dks, n) =
              thermo.solve_p(rho, entropicvar, 1 - qv, qv, 0, 0);
        });
  }

  void add_diagnostics(
      std::vector<std::unique_ptr<Diagnostic>> &diagnostics) override {}
};

struct DoubleVortex {
  enum class PLANE { XZ, YZ, XY };

  static PLANE constexpr plane = PLANE::XZ;
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 5000000._fp;
  static real constexpr Ly = 5000000._fp;
  static real constexpr Lz = 5000000._fp;
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
  static real constexpr zc = 0.5_fp * Lz;
  static real constexpr c = 0.05_fp;
  static real constexpr a = 1.0_fp / 3.0_fp;
  static real constexpr D = 0.5_fp * Lx;

  static std::pair<real, real> YAKL_INLINE get_plane_coords(real x, real y,
                                                            real z) {
    if (plane == PLANE::XY) {
      return {x, y};
    }
    if (plane == PLANE::XZ) {
      return {x, z};
    }
    if (plane == PLANE::YZ) {
      return {y, z};
    }
  }

  static VecXYZ YAKL_INLINE coriolis_f(real x, real y, real z) {
    VecXYZ v = {0, 0, 0};
    if (plane == PLANE::XY) {
      v.w = coriolis;
    }
    if (plane == PLANE::XZ) {
      v.v = ndims > 1 ? -coriolis : coriolis;
    }
    if (plane == PLANE::YZ) {
      v.u = coriolis;
    }

    return v;
  }

  static real YAKL_INLINE h_f(real xx, real yy, real zz) {
    auto [x, y] = get_plane_coords(xx, yy, zz);

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

  static VecXYZ YAKL_INLINE v_f(real xx, real yy, real zz) {
    auto [x, y] = get_plane_coords(xx, yy, zz);

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

    real u1 =
        -g * dh / coriolis / sigmay *
        (yprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         yprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));
    real u2 =
        g * dh / coriolis / sigmax *
        (xprimeprime1 * exp(-0.5_fp * (xprime1 * xprime1 + yprime1 * yprime1)) +
         xprimeprime2 * exp(-0.5_fp * (xprime2 * xprime2 + yprime2 * yprime2)));

    VecXYZ vvec = {0, 0, 0};

    if (plane == PLANE::XY) {
      vvec.u = u1;
      vvec.v = u2;
    }
    if (plane == PLANE::XZ) {
      vvec.u = u1;
      vvec.w = u2;
    }
    if (plane == PLANE::YZ) {
      vvec.v = u1;
      vvec.w = u2;
    }
    return vvec;
  }

  static real YAKL_INLINE S_f(real xx, real yy, real zz) {
    auto [x, y] = get_plane_coords(xx, yy, zz);
    real sval =
        g * (1._fp + c * exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) /
                             (a * a * D * D)));
    return sval * h_f(xx, yy, zz);
  }
};

template <bool acoustic_balance> struct RisingBubble {
  static int constexpr max_ndims = 2;
  static real constexpr g = 9.80616_fp;
  static real constexpr Lx = 1000._fp;
  static real constexpr Ly = 1000._fp;
  static real constexpr Lz = 1500._fp;
  static real constexpr xc = 0.5_fp * Lx;
  static real constexpr yc = 0.5_fp * Ly;
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

    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rho_b = isentropic_rho(z, theta0, g, thermo);
    if (acoustic_balance) {
      real theta = entropicvar_f(x, y, z, thermo);
      return rho_b * theta0 / theta;
    } else {
      return rho_b;
    }
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real r;
    if (ndims == 1) {
      r = sqrt((x - xc) * (x - xc) + (z - bzc) * (z - bzc));
    } else {
      r = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc) +
               (z - bzc) * (z - bzc));
    }
    real dtheta = (r < rc) ? dss * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.v = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct TwoBubbles {
  static int constexpr max_ndims = 1;
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
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    return isentropic_rho(z, theta0, g, thermo);
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);

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
    return thermo.compute_entropic_var_from_p_T(p, T + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct DensityCurrent {
  static int constexpr max_ndims = 1;
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
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rho_b = isentropic_rho(z, theta0, g, thermo);
    return rho_b;
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real r = sqrt((x - bxc) * (x - bxc) / (bxr * bxr) +
                  (z - bzc) * (z - bzc) / (bzr * bzr));
    real dtheta = (r < 1) ? dss * 0.5_fp * (1._fp + cos(pi * r)) : 0._fp;
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct MoistRisingBubble : public RisingBubble<false> {

  static real YAKL_INLINE rhov_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {
    real r = sqrt((x - xc) * (x - xc) + (z - bzc) * (z - bzc));
    real rh = (r < rc) ? rh0 * 0.5_fp * (1._fp + cos(pi * r / rc)) : 0._fp;
    real Th = isentropic_T(z, theta0, g, thermo);
    real svp = saturation_vapor_pressure(Th);
    real pv = svp * rh;
    return pv / (thermo.cst.Rv * Th);
  }

  static real YAKL_INLINE rhod_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rhod = rhod_f(x, y, z, thermo);
    real rhov = rhov_f(x, y, z, thermo);
    return rhod + rhov;
  }

  static real YAKL_INLINE refrhov_f(real z, const ThermoPotential &thermo) {
    return 0;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct LargeRisingBubble {
  static int constexpr max_ndims = 1;
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
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    return isentropic_rho(z, theta0, g, thermo);
  }
  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {

    real p = isentropic_p(z, theta0, g, thermo);
    real T0 = isentropic_T(z, theta0, g, thermo);
    real dtheta = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_theta);
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T0 + dT, 1, 0, 0, 0);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = 0;
    vvec.w = 0;
    return vvec;
  }

  static void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {}
};

struct MoistLargeRisingBubble : LargeRisingBubble {

  static real YAKL_INLINE rhod_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {
    real p = isentropic_p(z, theta0, g, thermo);
    real T = isentropic_T(z, theta0, g, thermo);
    real alpha = thermo.compute_alpha(p, T, 1, 0, 0, 0);
    return 1._fp / alpha;
  }

  static real YAKL_INLINE rhov_f(real x, real y, real z,
                                 const ThermoPotential &thermo) {

    real pert = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_vapor);
    real Th = isentropic_T(z, theta0, g, thermo);
    real svp = saturation_vapor_pressure(Th);
    real pv = svp * pert;
    return pv / (thermo.cst.Rv * Th);
  }

  static real YAKL_INLINE refrhov_f(real z, const ThermoPotential &thermo) {
    return 0;
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
    real rhod = rhod_f(x, y, z, thermo);
    real rhov = rhov_f(x, y, z, thermo);
    return rhod + rhov;
  }

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
                                        const ThermoPotential &thermo) {

    real p = isentropic_p(z, theta0, g, thermo);
    real T0 = isentropic_T(z, theta0, g, thermo);
    real dtheta = linear_ellipsoid(x, z, xc, bzc, xrad, zrad, amp_theta);
    real dT = dtheta * pow(p / thermo.cst.pr, thermo.cst.kappa_d);
    return thermo.compute_entropic_var_from_p_T(p, T0 + dT, 1, 0, 0, 0);
  }
};

template <bool add_perturbation> struct GravityWave {
  static int constexpr max_ndims = 1;
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
    return rho_ref *
           thermo.compute_entropic_var_from_p_T(p_ref, T_ref, 1, 0, 0, 0);
  }

  static real YAKL_INLINE rho_f(real x, real y, real z,
                                const ThermoPotential &thermo) {
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

  static real YAKL_INLINE entropicvar_f(real x, real y, real z,
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

    return thermo.compute_entropic_var_from_p_T(p, T, 1, 0, 0, 0);
  }

  static real YAKL_INLINE entropicdensity_f(real x, real y, real z,
                                            const ThermoPotential &thermo) {
    return rho_f(x, y, z, thermo) * entropicvar_f(x, y, z, thermo);
  }

  static VecXYZ YAKL_INLINE v_f(real x, real y, real z) {
    VecXYZ vvec;
    vvec.u = u_0;
    vvec.w = 0;
    return vvec;
  }

  static real YAKL_INLINE rhoexact_f(real x, real y, real z, real t,
                                     const ThermoPotential &thermo) {
    real rho = refrho_f(z, thermo);
    if (add_perturbation) {
      rho += sum_series(x, z, t, thermo).drho;
    }
    return rho;
  }

  static real YAKL_INLINE entropicdensityexact_f(
      real x, real y, real z, real t, const ThermoPotential &thermo) {

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

    return rho * thermo.compute_entropic_var_from_p_T(p, T, 1, 0, 0, 0);
  }

  static real YAKL_INLINE Texact_f(real x, real y, real z, real t,
                                   const ThermoPotential &thermo) {
    real T = T_ref;
    if (add_perturbation) {
      T += sum_series(x, z, t, thermo).dT;
    }
    return T;
  }

  static VecXYZ YAKL_INLINE vexact_f(real x, real y, real z, real t,
                                     const ThermoPotential &thermo) {
    VecXYZ v;

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

  struct ExactDensityDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "dense";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return rhoexact_f(x, y, z, time, thermo);
          },
          field, 0);

      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return entropicdensityexact_f(x, y, z, time, thermo);
          },
          field, 1);
    }
  };

  struct ExactTemperatureDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "Te";
      topology = pgeom.topology;
      dofs_arr = {0, 0, 1};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      dual_geometry.set_00form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return Texact_f(x, y, z, time, thermo);
          },
          field, 0);
    }
  };

  struct ExactWDiagnostic : public Diagnostic {
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "we";
      topology = pgeom.topology;
      dofs_arr = {0, 1, 1};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      primal_geometry.set_01form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return vexact_f(x, y, z, time, thermo);
          },
          field, 0);
    }
  };

  struct BackgroundDensityDiagnostic : public Diagnostic {

    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "densb";
      topology = dgeom.topology;
      dofs_arr = {1, 1, 2};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      YAKL_SCOPE(thermo, equations->thermo);
      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) { return refrho_f(z, thermo); },
          field, 0);

      dual_geometry.set_n1form_values(
          YAKL_LAMBDA(real x, real y, real z) {
            return refentropicdensity_f(z, thermo);
          },
          field, 1);
    }
  };

  struct TemperatureDiagnostic : public Diagnostic {
    static constexpr bool linear_T = true;
    void initialize(const Geometry<Straight> &pgeom,
                    const Geometry<Twisted> &dgeom,
                    Equations &equations) override {
      name = "T";
      topology = pgeom.topology;
      dofs_arr = {0, 0, 1};
      Diagnostic::initialize(pgeom, dgeom, equations);
    }

    void compute(real time, const FieldSet<nconstant> &const_vars,
                 const FieldSet<nprognostic> &x) override {

      const auto &primal_topology = primal_geometry.topology;

      int pis = primal_topology.is;
      int pjs = primal_topology.js;
      int pks = primal_topology.ks;

      const auto &densvar = x.fields_arr[DENSVAR].data;
      const auto &refstate = this->equations->reference_state;
      const auto &refdens = refstate.dens.data;

      YAKL_SCOPE(thermo, equations->thermo);
      YAKL_SCOPE(varset, equations->varset);
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
              T = thermo.compute_T_from_alpha(alpha, entropic_var, 1, 0, 0, 0);
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
                          bool acoustic_balance) {
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
  } else if (name == "doublevortex") {
    testcase = std::make_unique<SWETestCase<DoubleVortex>>();
  } else {
    throw std::runtime_error("unknown test case");
  }
}

#ifdef PAM_STANDALONE
void testcase_from_config(std::unique_ptr<TestCase> &testcase,
                          const YAML::Node &config) {
  const std::string name = config["init_data"].as<std::string>();
  const bool acoustic_balance =
      config["balance_initial_density"].as<bool>(false);
  testcase_from_string(testcase, name, acoustic_balance);
}
#endif

std::unique_ptr<TestCase> make_coupled_test_case(PamCoupler &coupler) {
  return std::make_unique<CoupledTestCase>(coupler);
}
} // namespace pamc
