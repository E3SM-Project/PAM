#pragma once

namespace pamc {

struct PressureLinearSystem : public LinearSystem {
  yakl::RealFFT1D<real> fftp_x;
  yakl::RealFFT1D<real> fftp_y;

  bool is_initialized = false;

  int nxf, nyf;
  real4d p_transform;

  real4d tri_l;
  real4d tri_d;
  real4d tri_u;
  real4d tri_c;

  using VS = VariableSet;

public:
  void initialize(ModelParameters &params,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  Equations &equations) override {

    LinearSystem::initialize(params, primal_geom, dual_geom, equations);

    const auto &primal_topology = primal_geom.topology;

    auto pni = primal_topology.ni;
    auto pnl = primal_topology.nl;
    auto nx = primal_topology.n_cells_x;
    auto ny = primal_topology.n_cells_y;
    auto nens = primal_topology.nens;

    this->nxf = nx + 2 - nx % 2;
    this->nyf = ndims > 1 ? ny + 2 - ny % 2 : ny;

    p_transform = real4d("p transform", pni, nyf, nxf, nens);
    yakl::memset(p_transform, 0);

    fftp_x.init(p_transform, 2, nx);
    if (ndims > 1) {
      fftp_y.init(p_transform, 1, ny);
    }

    tri_d = real4d("tri d", pni, nyf, nxf, nens);
    tri_l = real4d("tri l", pni, nyf, nxf, nens);
    tri_u = real4d("tri u", pni, nyf, nxf, nens);
    tri_c = real4d("tri c", pni, nyf, nxf, nens);

    this->is_initialized = true;
  }

  virtual void prepare_pressure_rhs(real dt, FieldSet<nprognostic> &rhs,
                                    FieldSet<nconstant> &const_vars,
                                    FieldSet<nauxiliary> &auxiliary_vars) = 0;

  virtual void solve_for_pressure(real dt, FieldSet<nprognostic> &rhs,
                                  FieldSet<nconstant> &const_vars,
                                  FieldSet<nauxiliary> &auxiliary_vars) = 0;

  virtual void update_velocity(real dt, FieldSet<nprognostic> &rhs,
                               FieldSet<nconstant> &const_vars,
                               FieldSet<nauxiliary> &auxiliary_vars,
                               FieldSet<nprognostic> &solution) = 0;

  virtual void update_densities(real dt, FieldSet<nprognostic> &rhs,
                                FieldSet<nconstant> &const_vars,
                                FieldSet<nauxiliary> &auxiliary_vars,
                                FieldSet<nprognostic> &solution) = 0;

  void solve(real dt, FieldSet<nprognostic> &rhs,
             FieldSet<nconstant> &const_vars,
             FieldSet<nauxiliary> &auxiliary_vars,
             FieldSet<nprognostic> &solution) override {

    prepare_pressure_rhs(dt, rhs, const_vars, auxiliary_vars);
    solve_for_pressure(dt, rhs, const_vars, auxiliary_vars);
    update_velocity(dt, rhs, const_vars, auxiliary_vars, solution);
    update_densities(dt, rhs, const_vars, auxiliary_vars, solution);
  }
};
} // namespace pamc
