#pragma once

#include "common.h"
#include "field_sets.h"
#include "functionals.h"
#include "geometry.h"
#include "hamiltonian.h"
#include "refstate.h"
#include "topology.h"
#include "variableset.h"
#include "weno_func_recon.h" // needed to set TransformMatrices related stuff
#include "weno_func_recon_variable.h" // needed to set TransformMatrices related stuff

class Equations {
public:
  Hamiltonian Hs;
#ifdef _LAYER
  Hamiltonian_Hk Hk;
  Functional_PVPE PVPE;
#elif _EXTRUDED
  Hamiltonian_Hk_extruded Hk;
  Functional_PVPE_extruded PVPE;
#endif
  VariableSet varset;
  ThermoPotential thermo;
  ReferenceState reference_state;
  bool is_initialized;

  Equations() { this->is_initialized = false; }
  void initialize(PamCoupler &coupler, ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {

    this->reference_state.initialize<VariableSet>(primal_geom.topology,
                                                  dual_geom.topology);
    this->varset.initialize(coupler, params, thermo, reference_state,
                            primal_geom, dual_geom);
    this->PVPE.initialize(varset);
    this->Hk.initialize(varset, primal_geom, dual_geom);
    this->Hs.initialize(thermo, varset, primal_geom, dual_geom);
    this->is_initialized = true;
  }
};

class Diagnostic {
public:
  std::string name;
  std::array<int, 3> dofs_arr;
  Topology topology;
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  Equations *equations;
  Field field;

  bool is_initialized;

  Diagnostic() { this->is_initialized = false; }

  virtual void compute(real time, const FieldSet<nconstant> &const_vars,
                       const FieldSet<nprognostic> &x) = 0;

  virtual void initialize(const Geometry<Straight> &pgeom,
                          const Geometry<Twisted> &dgeom, Equations &eqs) {

    this->primal_geometry = pgeom;
    this->dual_geometry = dgeom;
    this->equations = &eqs;
    field.initialize(topology, nullptr, name, dofs_arr[0], dofs_arr[1],
                     dofs_arr[2]);
    this->is_initialized = true;
  }
  virtual ~Diagnostic() = default;
};

enum class TRACER_TAG { CONSTANT, SQUARE, DOUBLESQUARE, GAUSSIAN };

struct TracerFunctor {
  YAKL_INLINE real operator()(TRACER_TAG tracer_tag, real x, real y, real Lx,
                              real Ly, real xc, real yc) {
    real tracer = -1;

    if (tracer_tag == TRACER_TAG::CONSTANT) {
      tracer = 1000;
    } else if (tracer_tag == TRACER_TAG::SQUARE) {
      tracer = (x > 0.35_fp * Lx && x < 0.65_fp * Lx && y > 0.35_fp * Ly &&
                y < 0.65_fp * Ly)
                   ? 0.005_fp
                   : 0.;
    } else if (tracer_tag == TRACER_TAG::DOUBLESQUARE) {
      const real ur = (x > 0.6_fp * Lx && x < 0.9_fp * Lx && y > 0.6_fp * Ly &&
                       y < 0.9_fp * Ly)
                          ? 0.005_fp
                          : 0.;
      const real ll = (x > 0.1_fp * Lx && x < 0.4_fp * Lx && y > 0.1_fp * Ly &&
                       y < 0.4_fp * Ly)
                          ? 0.005_fp
                          : 0.;

      tracer = ur + ll;
    } else if (tracer_tag == TRACER_TAG::GAUSSIAN) {
      const real a = 1.0_fp / 3.0_fp;
      const real D = 0.5_fp * Lx;
      tracer = 0.005_fp * exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) /
                              (a * a * D * D));
    }

    return tracer;
  }
};

class TestCase {
public:
  SArray<TRACER_TAG, 1, ntracers_dycore> tracers;
  Equations *equations;
  bool is_initialized;

  TestCase() { this->is_initialized = false; }

  virtual void initialize(Equations &eqs) {
    this->equations = &eqs;
    this->is_initialized = true;
  }

  virtual void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics){};

  void set_tracers(ModelParameters &params) {
    for (int i = 0; i < ntracers_dycore; i++) {
      if (params.tracerdataStr[i] == "gaussian") {
        tracers(i) = TRACER_TAG::GAUSSIAN;
      } else if (params.tracerdataStr[i] == "square") {
        tracers(i) = TRACER_TAG::SQUARE;
      } else if (params.tracerdataStr[i] == "doublesquare") {
        tracers(i) = TRACER_TAG::DOUBLESQUARE;
      } else if (params.tracerdataStr[i] == "constant") {
        tracers(i) = TRACER_TAG::CONSTANT;
      }
    }
  }

  virtual void set_domain(ModelParameters &params) = 0;
  virtual std::array<real, 3> get_domain() const = 0;
  virtual void set_initial_conditions(FieldSet<nprognostic> &progvars,
                                      FieldSet<nconstant> &constvars,
                                      const Geometry<Straight> &primal_geom,
                                      const Geometry<Twisted> &dual_geom) = 0;
  virtual void set_reference_state(const Geometry<Straight> &primal_geom,
                                   const Geometry<Twisted> &dual_geom){};
  virtual ~TestCase() = default;
};

class Tendencies {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  Equations *equations;

  SArray<real, 2, reconstruction_order, 2> primal_to_gll;
  SArray<real, 3, reconstruction_order, reconstruction_order,
         reconstruction_order>
      primal_wenoRecon;
  SArray<real, 1, (reconstruction_order - 1) / 2 + 2> primal_wenoIdl;
  real primal_wenoSigma;

  SArray<real, 2, dual_reconstruction_order, 2> dual_to_gll;
  SArray<real, 3, dual_reconstruction_order, dual_reconstruction_order,
         dual_reconstruction_order>
      dual_wenoRecon;
  SArray<real, 1, (dual_reconstruction_order - 1) / 2 + 2> dual_wenoIdl;
  real dual_wenoSigma;

  SArray<real, 2, coriolis_reconstruction_order, 2> coriolis_to_gll;
  SArray<real, 3, coriolis_reconstruction_order, coriolis_reconstruction_order,
         coriolis_reconstruction_order>
      coriolis_wenoRecon;
  SArray<real, 1, (coriolis_reconstruction_order - 1) / 2 + 2> coriolis_wenoIdl;
  real coriolis_wenoSigma;

  bool is_initialized;

  Tendencies() { this->is_initialized = false; }

  void initialize(ModelParameters &params, Equations &equations,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    this->equations = &equations;
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;

    TransformMatrices::coefs_to_gll_lower(primal_to_gll);
    TransformMatrices::weno_sten_to_coefs(primal_wenoRecon);
    wenoSetIdealSigma<reconstruction_order>(primal_wenoIdl, primal_wenoSigma);

    TransformMatrices::coefs_to_gll_lower(dual_to_gll);
    TransformMatrices::weno_sten_to_coefs(dual_wenoRecon);
    wenoSetIdealSigma<dual_reconstruction_order>(dual_wenoIdl, dual_wenoSigma);

    TransformMatrices::coefs_to_gll_lower(coriolis_to_gll);
    TransformMatrices::weno_sten_to_coefs(coriolis_wenoRecon);
    wenoSetIdealSigma<coriolis_reconstruction_order>(coriolis_wenoIdl,
                                                     coriolis_wenoSigma);

    this->is_initialized = true;
  }
  virtual void compute_constants(FieldSet<nconstant> &const_vars,
                                 FieldSet<nprognostic> &x) = 0;

  virtual void compute_functional_derivatives(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x,
      FieldSet<nauxiliary> &auxiliary_vars, real fac = 1,
      ADD_MODE addmode = ADD_MODE::REPLACE) = 0;

  virtual void compute_functional_derivatives_two_point(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x1,
      FieldSet<nprognostic> &x2, FieldSet<nauxiliary> &auxiliary_vars) {
    throw std::runtime_error("Not implemented");
  }

  virtual void apply_symplectic(real dt, FieldSet<nconstant> &const_vars,
                                FieldSet<nprognostic> &x,
                                FieldSet<nauxiliary> &auxiliary_vars,
                                FieldSet<nprognostic> &xtend,
                                ADD_MODE addmode = ADD_MODE::REPLACE) = 0;

  virtual void project_to_anelastic(FieldSet<nconstant> &const_vars,
                                    FieldSet<nprognostic> &x,
                                    FieldSet<nauxiliary> &auxiliary_vars) {}

  virtual void add_pressure_perturbation(real dt,
                                         FieldSet<nconstant> &const_vars,
                                         FieldSet<nprognostic> &x,
                                         FieldSet<nauxiliary> &auxiliary_vars,
                                         FieldSet<nprognostic> &xtend) {}

  virtual void compute_rhs(real dt, FieldSet<nconstant> &const_vars,
                           FieldSet<nprognostic> &x,
                           FieldSet<nauxiliary> &auxiliary_vars,
                           FieldSet<nprognostic> &xtend,
                           ADD_MODE addmode = ADD_MODE::REPLACE) {
    compute_functional_derivatives(dt, const_vars, x, auxiliary_vars);
    apply_symplectic(dt, const_vars, x, auxiliary_vars, xtend, addmode);
    add_pressure_perturbation(dt, const_vars, x, auxiliary_vars, xtend);
  }

  virtual void remove_negative_densities(FieldSet<nprognostic> &x) {}
};

class ExtrudedTendencies : public Tendencies {
public:
  SArray<real, 1, (vert_reconstruction_order - 1) / 2 + 2> primal_vert_wenoIdl;
  SArray<real, 2, vert_reconstruction_order, 2> primal_vert_to_gll;
  SArray<real, 3, vert_reconstruction_order, vert_reconstruction_order,
         vert_reconstruction_order>
      primal_vert_wenoRecon;
  real primal_vert_wenoSigma;
  real4d primal_vert_coefs_to_gll_arr;
  real4d primal_vert_sten_to_gll_arr;
  real4d primal_vert_sten_to_coefs_arr;
  real5d primal_vert_weno_recon_lower_arr;

  SArray<real, 1, (coriolis_vert_reconstruction_order - 1) / 2 + 2>
      coriolis_vert_wenoIdl;
  SArray<real, 2, coriolis_vert_reconstruction_order, 2> coriolis_vert_to_gll;
  SArray<real, 3, coriolis_vert_reconstruction_order,
         coriolis_vert_reconstruction_order, coriolis_vert_reconstruction_order>
      coriolis_vert_wenoRecon;
  real coriolis_vert_wenoSigma;
  real4d coriolis_vert_coefs_to_gll_arr;
  real4d coriolis_vert_sten_to_gll_arr;
  real4d coriolis_vert_sten_to_coefs_arr;
  real5d coriolis_vert_weno_recon_lower_arr;

  SArray<real, 1, (dual_vert_reconstruction_order - 1) / 2 + 2>
      dual_vert_wenoIdl;
  SArray<real, 3, dual_vert_reconstruction_order,
         dual_vert_reconstruction_order, dual_vert_reconstruction_order>
      dual_vert_wenoRecon;
  SArray<real, 2, dual_vert_reconstruction_order, 2> dual_vert_to_gll;
  real dual_vert_wenoSigma;
  real4d dual_vert_coefs_to_gll_arr;
  real4d dual_vert_sten_to_gll_arr;
  real4d dual_vert_sten_to_coefs_arr;
  real5d dual_vert_weno_recon_lower_arr;

  void initialize(ModelParameters &params, Equations &equations,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {

    Tendencies::initialize(params, equations, primal_geom, dual_geom);

    const auto &primal_topology = primal_geom.topology;
    const auto &dual_topology = dual_geom.topology;

    if (primal_geom.uniform_vertical) {
      TransformMatrices::coefs_to_gll_lower(primal_vert_to_gll);
      TransformMatrices::weno_sten_to_coefs(primal_vert_wenoRecon);
      wenoSetIdealSigma<vert_reconstruction_order>(primal_vert_wenoIdl,
                                                   primal_vert_wenoSigma);

      TransformMatrices::coefs_to_gll_lower(dual_vert_to_gll);
      TransformMatrices::weno_sten_to_coefs(dual_vert_wenoRecon);
      wenoSetIdealSigma<dual_vert_reconstruction_order>(dual_vert_wenoIdl,
                                                        dual_vert_wenoSigma);

      TransformMatrices::coefs_to_gll_lower(coriolis_vert_to_gll);
      TransformMatrices::weno_sten_to_coefs(coriolis_vert_wenoRecon);
      wenoSetIdealSigma<coriolis_vert_reconstruction_order>(
          coriolis_vert_wenoIdl, coriolis_vert_wenoSigma);
    } else {
      this->primal_vert_coefs_to_gll_arr =
          real4d("Primal Vert Coeffs to GLL Array", primal_topology.nl,
                 vert_reconstruction_order, 2, primal_topology.nens);

      this->primal_vert_sten_to_gll_arr =
          real4d("Primal Vert Sten to GLL Array", primal_topology.nl,
                 vert_reconstruction_order, 2, primal_topology.nens);

      this->primal_vert_sten_to_coefs_arr =
          real4d("Primal Vert Sten to Coeffs Array", primal_topology.nl,
                 vert_reconstruction_order, vert_reconstruction_order,
                 primal_topology.nens);

      this->primal_vert_weno_recon_lower_arr =
          real5d("Primal Vert wenoReconLower Array", primal_topology.nl,
                 (vert_reconstruction_order - 1) / 2 + 1,
                 (vert_reconstruction_order - 1) / 2 + 1,
                 (vert_reconstruction_order - 1) / 2 + 1, primal_topology.nens);

      this->coriolis_vert_coefs_to_gll_arr =
          real4d("Coriolis Vert Coeffs to GLL Array", primal_topology.nl,
                 coriolis_vert_reconstruction_order, 2, primal_topology.nens);

      this->coriolis_vert_sten_to_gll_arr =
          real4d("Coriolis Vert Sten to GLL Array", primal_topology.nl,
                 coriolis_vert_reconstruction_order, 2, primal_topology.nens);

      this->coriolis_vert_sten_to_coefs_arr =
          real4d("Coriolis Vert Sten to Coeffs Array", primal_topology.nl,
                 coriolis_vert_reconstruction_order,
                 coriolis_vert_reconstruction_order, primal_topology.nens);

      this->coriolis_vert_weno_recon_lower_arr =
          real5d("Coriolis Vert wenoReconLower Array", primal_topology.nl,
                 (coriolis_vert_reconstruction_order - 1) / 2 + 1,
                 (coriolis_vert_reconstruction_order - 1) / 2 + 1,
                 (coriolis_vert_reconstruction_order - 1) / 2 + 1,
                 primal_topology.nens);

      this->dual_vert_coefs_to_gll_arr =
          real4d("Dual Vert Coeffs to GLL Array", dual_topology.nl,
                 dual_vert_reconstruction_order, 2, dual_topology.nens);

      this->dual_vert_sten_to_gll_arr =
          real4d("Dual Vert Sten to GLL Array", dual_topology.nl,
                 dual_vert_reconstruction_order, 2, dual_topology.nens);

      this->dual_vert_sten_to_coefs_arr =
          real4d("Dual Vert Sten to Coeffs Array", dual_topology.nl,
                 dual_vert_reconstruction_order, dual_vert_reconstruction_order,
                 dual_topology.nens);

      this->dual_vert_weno_recon_lower_arr = real5d(
          "Dual Vert wenoReconLower Array", dual_topology.nl,
          (dual_vert_reconstruction_order - 1) / 2 + 1,
          (dual_vert_reconstruction_order - 1) / 2 + 1,
          (dual_vert_reconstruction_order - 1) / 2 + 1, dual_topology.nens);

      create_variable_WENO<Twisted, dual_vert_reconstruction_order>(
          dual_vert_coefs_to_gll_arr, dual_vert_sten_to_gll_arr,
          dual_vert_sten_to_coefs_arr, dual_vert_weno_recon_lower_arr,
          dual_vert_wenoSigma, dual_vert_wenoIdl, dual_geom);

      create_variable_WENO<Straight, vert_reconstruction_order>(
          primal_vert_coefs_to_gll_arr, primal_vert_sten_to_gll_arr,
          primal_vert_sten_to_coefs_arr, primal_vert_weno_recon_lower_arr,
          primal_vert_wenoSigma, primal_vert_wenoIdl, primal_geom);

      create_variable_WENO<Straight, coriolis_vert_reconstruction_order>(
          coriolis_vert_coefs_to_gll_arr, coriolis_vert_sten_to_gll_arr,
          coriolis_vert_sten_to_coefs_arr, coriolis_vert_weno_recon_lower_arr,
          coriolis_vert_wenoSigma, coriolis_vert_wenoIdl, primal_geom);
    }
  }

  virtual real
  compute_max_anelastic_constraint(FieldSet<nprognostic> &x,
                                   FieldSet<nauxiliary> &auxiliary_vars,
                                   bool has_f_and_fw = false) = 0;
};

class LinearSystem {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  Equations *equations;

  bool is_initialized = false;

  virtual void initialize(ModelParameters &params,
                          const Geometry<Straight> &primal_geom,
                          const Geometry<Twisted> &dual_geom,
                          Equations &equations) {

    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->equations = &equations;

    this->is_initialized = true;
  }
  virtual void compute_coefficients(real dt){};

  virtual void solve(real dt, FieldSet<nprognostic> &rhs,
                     FieldSet<nconstant> &const_vars,
                     FieldSet<nauxiliary> &auxiliary_vars,
                     FieldSet<nprognostic> &solution) = 0;
};
