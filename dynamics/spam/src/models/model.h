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
  Geometry<Straight> primal_geometry;
  Topology dual_topology;
  Geometry<Twisted> dual_geometry;
  Field field;

  bool is_initialized;

  Diagnostic() { this->is_initialized = false; }

  virtual void compute(real time, const FieldSet<nconstant> &const_vars,
                       const FieldSet<nprognostic> &x) = 0;

  virtual void initialize(const Topology &ptopo, const Topology &dtopo,
                          const Geometry<Straight> &pgeom,
                          const Geometry<Twisted> &dgeom) {
    this->primal_topology = ptopo;
    this->dual_topology = dtopo;
    this->primal_geometry = pgeom;
    this->dual_geometry = dgeom;
    field.initialize(topology, name, dofs_arr[0], dofs_arr[1], dofs_arr[2]);
    this->is_initialized = true;
  }
  virtual ~Diagnostic() = default;
};

enum class TRACER_TAG { CONSTANT, SQUARE, DOUBLESQUARE, GAUSSIAN };
struct Tracer {
  YAKL_INLINE virtual real compute(real x, real y, real Lx, real Ly, real xc,
                                   real yc) = 0;
};

struct TracerConstant : Tracer {
  YAKL_INLINE real compute(real x, real y, real Lx, real Ly, real xc,
                           real yc) override {
    return 1000;
  }
};

struct TracerSquare : Tracer {
  YAKL_INLINE real compute(real x, real y, real Lx, real Ly, real xc,
                           real yc) override {
    return (x > 0.35_fp * Lx && x < 0.65_fp * Lx && y > 0.35_fp * Ly &&
            y < 0.65_fp * Ly)
               ? 0.005_fp
               : 0.;
  }
};

struct TracerDoubleSquare : Tracer {
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
  YAKL_INLINE real compute(real x, real y, real Lx, real Ly, real xc,
                           real yc) override {

    return tracer_square_ur(x, y, Lx, Ly, xc, yc) +
           tracer_square_ll(x, y, Lx, Ly, xc, yc);
  }
};

struct TracerGaussian : Tracer {
  YAKL_INLINE real compute(real x, real y, real Lx, real Ly, real xc,
                           real yc) override {
    real const a = 1.0_fp / 3.0_fp;
    real const D = 0.5_fp * Lx;
    return 0.005_fp *
           exp(-((x - xc) * (x - xc) + (y - yc) * (y - yc)) / (a * a * D * D));
  }
};

class TestCase {
public:
  using TracerArr = yakl::Array<Tracer *, 1, yakl::memDevice, yakl::styleC>;
  TracerArr tracer_f;

  TestCase() { this->tracer_f = TracerArr("tracer_f", ntracers_dycore); }

  virtual void
  add_diagnostics(std::vector<std::unique_ptr<Diagnostic>> &diagnostics) = 0;

  void set_tracers(ModelParameters &params) {
    SArray<TRACER_TAG, 1, ntracers_dycore> tracer_tag;
    for (int i = 0; i < ntracers_dycore; i++) {
      if (params.tracerdataStr[i] == "gaussian") {
        tracer_tag(i) = TRACER_TAG::GAUSSIAN;
      } else if (params.tracerdataStr[i] == "square") {
        tracer_tag(i) = TRACER_TAG::SQUARE;
      } else if (params.tracerdataStr[i] == "doublesquare") {
        tracer_tag(i) = TRACER_TAG::DOUBLESQUARE;
      } else {
        // by default set tracers to constant
        tracer_tag(i) = TRACER_TAG::CONSTANT;
      }
    }

    YAKL_SCOPE(tracer_f, this->tracer_f);
    parallel_for(
        ntracers_dycore, YAKL_LAMBDA(int i) {
          if (tracer_tag(i) == TRACER_TAG::GAUSSIAN) {
            tracer_f(i) = new TracerGaussian();
          }
          if (tracer_tag(i) == TRACER_TAG::SQUARE) {
            tracer_f(i) = new TracerSquare();
          }
          if (tracer_tag(i) == TRACER_TAG::DOUBLESQUARE) {
            tracer_f(i) = new TracerDoubleSquare();
          }
          if (tracer_tag(i) == TRACER_TAG::CONSTANT) {
            tracer_f(i) = new TracerConstant();
          }
        });
  }

  virtual void set_domain(ModelParameters &params) = 0;
  virtual void set_initial_conditions(FieldSet<nprognostic> &progvars,
                                      FieldSet<nconstant> &constvars,
                                      ExchangeSet<nconstant> &const_exchange,
                                      const Geometry<Straight> &primal_geom,
                                      const Geometry<Twisted> &dual_geom) = 0;
  virtual ~TestCase() = default;

  // why doesn't this work ? Tracers need to be deallocated !
  // virtual ~TestCase() {
  //   YAKL_SCOPE(tracer_f, this->tracer_f);
  //   parallel_for(ntracers_dycore, YAKL_LAMBDA(int i) {
  //     delete tracer_f(i);
  //   });
  // }
};

class Tendencies {
public:
  Topology primal_topology;
  Topology dual_topology;
  ExchangeSet<nauxiliary> *aux_exchange;
  ExchangeSet<nconstant> *const_exchange;
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

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
                  Topology &dual_topo, const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  ExchangeSet<nauxiliary> &aux_exchange,
                  ExchangeSet<nconstant> &const_exchange) {
    this->primal_topology = primal_topo;
    this->dual_topology = dual_topo;
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
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
                                 FieldSet<nprognostic> &x) = 0;

  virtual void compute_linrhs(real dt, FieldSet<nconstant> &const_vars,
                              FieldSet<nprognostic> &x,
                              FieldSet<nauxiliary> &auxiliary_vars,
                              FieldSet<nprognostic> &xtend) = 0;

  virtual void YAKL_INLINE compute_functional_derivatives(
      ADD_MODE addmode, real fac, real dt, FieldSet<nconstant> &const_vars,
      FieldSet<nprognostic> &x, FieldSet<nauxiliary> &auxiliary_vars) = 0;

  virtual void YAKL_INLINE apply_symplectic(
      real dt, FieldSet<nconstant> &const_vars, FieldSet<nprognostic> &x,
      FieldSet<nauxiliary> &auxiliary_vars, FieldSet<nprognostic> &xtend) = 0;

  virtual void YAKL_INLINE compute_rhs(real dt, FieldSet<nconstant> &const_vars,
                                       FieldSet<nprognostic> &x,
                                       FieldSet<nauxiliary> &auxiliary_vars,
                                       FieldSet<nprognostic> &xtend) {
    compute_functional_derivatives(ADD_MODE::REPLACE, 1._fp, dt, const_vars, x,
                                   auxiliary_vars);
    apply_symplectic(dt, const_vars, x, auxiliary_vars, xtend);
    // compute_linrhs(dt, const_vars,
    //                x,
    //                auxiliary_vars,
    //                xtend);
  }
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
                  Topology &dual_topo, const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  ExchangeSet<nauxiliary> &aux_exchange,
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

class LinearSystem {
public:
  Topology primal_topology;
  Topology dual_topology;
  ExchangeSet<nprognostic> *prog_exchange;
  ExchangeSet<nauxiliary> *aux_exchange;
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

  bool is_initialized = false;

  virtual void initialize(ModelParameters &params, Tendencies *tend,
                          FieldSet<nprognostic> &x,
                          FieldSet<nconstant> &const_vars,
                          FieldSet<nauxiliary> &auxiliary_vars,
                          ExchangeSet<nprognostic> &prog_exchange) {

    this->primal_topology = tend->primal_topology;
    this->dual_topology = tend->dual_topology;
    this->primal_geometry = tend->primal_geometry;
    this->dual_geometry = tend->dual_geometry;
    this->aux_exchange = tend->aux_exchange;
    this->prog_exchange = &prog_exchange;

    this->is_initialized = true;
  }
  virtual void compute_coefficients(real dt,FieldSet<nconstant> &constvars) = 0;

  virtual void YAKL_INLINE solve(real dt, FieldSet<nprognostic> &rhs,
                                 FieldSet<nconstant> &const_vars,
                                 FieldSet<nauxiliary> &auxiliary_vars,
                                 FieldSet<nprognostic> &solution) = 0;
};
