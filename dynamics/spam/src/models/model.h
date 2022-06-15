#pragma once

#include "common.h"
#include "field_sets.h"
#include "geometry.h"
#include "topology.h"
#include "weno_func_recon.h" // needed to set TransformMatrices related stuff

class Diagnostic {
public:
  std::string name;
  std::array<int, 3> dofs_arr;
  Topology topology;
  Topology primal_topology;
  Geometry *primal_geometry;
  Topology dual_topology;
  Geometry *dual_geometry;
  Field field;

  bool is_initialized;

  Diagnostic() { this->is_initialized = false; }

  virtual void compute(real time, const FieldSet<nconstant> &const_vars,
                       const FieldSet<nprognostic> &x) = 0;

  virtual void initialize(const Topology &ptopo, const Topology &dtopo,
                          Geometry &pgeom, Geometry &dgeom) {
    this->primal_topology = ptopo;
    this->dual_topology = dtopo;
    this->primal_geometry = &pgeom;
    this->dual_geometry = &dgeom;
    field.initialize(topology, name, dofs_arr[0], dofs_arr[1], dofs_arr[2]);
    this->is_initialized = true;
  }
  virtual ~Diagnostic() = default;
};

real YAKL_INLINE tracer_constant(real x, real y, real Lx, real Ly, real xc,
                                 real yc) {
  return 1000;
}
real YAKL_INLINE tracer_square_cent(real x, real y, real Lx, real Ly, real xc,
                                    real yc) {
  return (x > 0.35_fp * Lx && x < 0.65_fp * Lx && y > 0.35_fp * Ly &&
          y < 0.65_fp * Ly)
             ? 0.005_fp
             : 0.;
}
real YAKL_INLINE tracer_square_ur(real x, real y, real Lx, real Ly, real xc,
                                  real yc) {
  return (x > 0.6_fp * Lx && x < 0.9_fp * Lx && y > 0.6_fp * Ly &&
          y < 0.9_fp * Ly)
             ? 0.005_fp
             : 0.;
}
real YAKL_INLINE tracer_square_ll(real x, real y, real Lx, real Ly, real xc,
                                  real yc) {
  return (x > 0.1_fp * Lx && x < 0.4_fp * Lx && y > 0.1_fp * Ly &&
          y < 0.4_fp * Ly)
             ? 0.005_fp
             : 0.;
}
real YAKL_INLINE tracer_square_urpll(real x, real y, real Lx, real Ly, real xc,
                                     real yc) {
  return tracer_square_ur(x, y, Lx, Ly, xc, yc) +
         tracer_square_ll(x, y, Lx, Ly, xc, yc);
}

real YAKL_INLINE tracer_gaussian(real x, real y, real Lx, real Ly, real xc,
                                 real yc) {
  real const a = 1.0_fp / 3.0_fp;
  real const D = 0.5_fp * Lx;
  return 0.005_fp *
         exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) / (a * a * D * D));
}

class TestCase {
public:
  using tracer_ptr = real (*)(real, real, real, real, real, real);
  std::array<tracer_ptr, ntracers_dycore> tracer_f;

  bool is_initialized;
  TestCase() { this->is_initialized = false; }

  void set_tracers(ModelParameters &params) {
    for (int i = 0; i < ntracers_dycore; i++) {
      if (params.tracerdataStr[i] == "gaussian") {
        tracer_f[i] = tracer_gaussian;
      } else if (params.tracerdataStr[i] == "square") {
        tracer_f[i] = tracer_square_cent;
      } else if (params.tracerdataStr[i] == "doublesquare") {
        tracer_f[i] = tracer_square_urpll;
      } else {
        // by default set tracers to constant
        tracer_f[i] = tracer_constant;
      }
    }
  }

  virtual void set_domain(ModelParameters &params) = 0;
  virtual void set_initial_conditions(FieldSet<nprognostic> &progvars,
                                      FieldSet<nconstant> &constvars,
                                      Geometry &primal_geom,
                                      Geometry &dual_geom) = 0;
  // virtual void add_diagnostics(real time, const FieldSet<nconstant>
  // &const_vars);
  virtual ~TestCase() = default;
};

class Tendencies {
public:
  Topology primal_topology;
  Topology dual_topology;
  ExchangeSet<nauxiliary> *aux_exchange;
  ExchangeSet<nconstant> *const_exchange;
  Geometry *primal_geometry;
  Geometry *dual_geometry;

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

  void initialize(ModelParameters &params, Topology &primal_topo,
                  Topology &dual_topo, Geometry &primal_geom,
                  Geometry &dual_geom, ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nconstant> &const_exchange) {
    this->primal_topology = primal_topo;
    this->dual_topology = dual_topo;
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;
    this->aux_exchange = &aux_exchange;
    this->const_exchange = &const_exchange;

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
                                 FieldSet<nprognostic> &x){};
  virtual void YAKL_INLINE compute_rhs(real dt, FieldSet<nconstant> &const_vars,
                                       FieldSet<nprognostic> &x,
                                       FieldSet<nauxiliary> &auxiliary_vars,
                                       FieldSet<nprognostic> &xtend){};
};

class ExtrudedTendencies : public Tendencies {
public:
  SArray<real, 2, vert_reconstruction_order, 2> primal_vert_to_gll;
  SArray<real, 3, vert_reconstruction_order, vert_reconstruction_order,
         vert_reconstruction_order>
      primal_vert_wenoRecon;
  SArray<real, 1, (vert_reconstruction_order - 1) / 2 + 2> primal_vert_wenoIdl;
  real primal_vert_wenoSigma;

  SArray<real, 2, dual_vert_reconstruction_order, 2> dual_vert_to_gll;
  SArray<real, 3, dual_vert_reconstruction_order,
         dual_vert_reconstruction_order, dual_vert_reconstruction_order>
      dual_vert_wenoRecon;
  SArray<real, 1, (dual_vert_reconstruction_order - 1) / 2 + 2>
      dual_vert_wenoIdl;
  real dual_vert_wenoSigma;

  SArray<real, 2, coriolis_vert_reconstruction_order, 2> coriolis_vert_to_gll;
  SArray<real, 3, coriolis_vert_reconstruction_order,
         coriolis_vert_reconstruction_order, coriolis_vert_reconstruction_order>
      coriolis_vert_wenoRecon;
  SArray<real, 1, (coriolis_vert_reconstruction_order - 1) / 2 + 2>
      coriolis_vert_wenoIdl;
  real coriolis_vert_wenoSigma;

  void initialize(ModelParameters &params, Topology &primal_topo,
                  Topology &dual_topo, Geometry &primal_geom,
                  Geometry &dual_geom, ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nconstant> &const_exchange) {
    Tendencies::initialize(params, primal_topo, dual_topo, primal_geom,
                           dual_geom, aux_exchange, const_exchange);

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

    this->is_initialized = true;
  }
};
