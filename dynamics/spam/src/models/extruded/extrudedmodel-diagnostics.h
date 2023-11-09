#pragma once

namespace pamc {

struct TotalDensityDiagnostic : public Diagnostic {
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    name = "total_dens";
    topology = dgeom.topology;
    dofs_arr = {1, 1, 1};
    Diagnostic::initialize(pgeom, dgeom, equations);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &dual_topology = dual_geometry.topology;

    const int dis = dual_topology.is;
    const int djs = dual_topology.js;
    const int dks = dual_topology.ks;

    const auto &densvar = x.fields_arr[DENSVAR].data;
    const auto &varset = this->equations->varset;

    parallel_for(
        "Compute total density diagnostic",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          const real total_dens =
              varset.get_total_density(densvar, k, j, i, dks, djs, dis, n);

          field.data(0, k + dks, j + djs, i + dis, n) = total_dens;
        });
  }
};

class Dens0Diagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    // concentration 0-forms for dens
    name = "densl";
    topology = pgeom.topology;
    dofs_arr = {
        0, 0,
        VariableSet::ndensity_prognostic}; // densldiag = straight (0,0)-form
    Diagnostic::initialize(pgeom, dgeom, equations);
  }

  void compute(real time, const FieldSet<nconstant> &const_vars,
               const FieldSet<nprognostic> &x) override {

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    parallel_for(
        "Compute DENS0 DIAG",
        SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                        primal_topology.n_cells_x, primal_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          compute_Hn1bar<VariableSet::ndensity_prognostic, diff_ord,
                         vert_diff_ord>(
              field.data, x.fields_arr[DENSVAR].data, this->primal_geometry,
              this->dual_geometry, pis, pjs, pks, i, j, k, n);
        });
  }
};

class QHZDiagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    name = "QHZl";
    topology = dgeom.topology;
    dofs_arr = {ndims - 1, 0, 1}; // // Qldiag = twisted (n-1,0)-form
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
        "Compute Q0 DIAG",
        SimpleBounds<4>(dual_topology.ni - 3, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_qhz(field.data, x.fields_arr[VVAR].data,
                           x.fields_arr[WVAR].data, x.fields_arr[DENSVAR].data,
                           const_vars.fields_arr[CORIOLISHZVAR].data, dis, djs,
                           dks, i, j, k + 2, n);
        });

    // Bottom is k=1 and top is k=dual_topology.ni-2
    parallel_for(
        "Compute Q0 DIAG TOP/BOTTOM",
        SimpleBounds<3>(dual_topology.n_cells_y, dual_topology.n_cells_x,
                        dual_topology.nens),
        YAKL_CLASS_LAMBDA(int j, int i, int n) {
          PVPE.compute_qhz_bottom(field.data, x.fields_arr[VVAR].data,
                                  x.fields_arr[WVAR].data,
                                  x.fields_arr[DENSVAR].data,
                                  const_vars.fields_arr[CORIOLISHZVAR].data,
                                  dis, djs, dks, i, j, 1, n);
          PVPE.compute_qhz_top(field.data, x.fields_arr[VVAR].data,
                               x.fields_arr[WVAR].data,
                               x.fields_arr[DENSVAR].data,
                               const_vars.fields_arr[CORIOLISHZVAR].data, dis,
                               djs, dks, i, j, dual_topology.ni - 2, n);
        });

    field.set_bnd(0.0);
  }
};

class QXYDiagnostic : public Diagnostic {
public:
  void initialize(const Geometry<Straight> &pgeom,
                  const Geometry<Twisted> &dgeom,
                  Equations &equations) override {
    name = "QXYl";
    topology = dgeom.topology;
    dofs_arr = {0, 1, 1}; // // Qxydiag = twisted (0,1)-form
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
        "Compute QXYl DIAG",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          PVPE.compute_qxy(field.data, x.fields_arr[VVAR].data,
                           x.fields_arr[DENSVAR].data,
                           const_vars.fields_arr[CORIOLISXYVAR].data, dis, djs,
                           dks, i, j, k, n);
          real dz = dual_geometry.get_dz(k + dks, n);
          field.data(0, k + dks, j + djs, i + dis, n) *= dz;
        });
  }
};

void add_model_diagnostics(
    std::vector<std::unique_ptr<Diagnostic>> &diagnostics) {
  diagnostics.emplace_back(std::make_unique<Dens0Diagnostic>());
  diagnostics.emplace_back(std::make_unique<QHZDiagnostic>());
  if (ndims > 1) {
    diagnostics.emplace_back(std::make_unique<QXYDiagnostic>());
  }
  diagnostics.emplace_back(std::make_unique<TotalDensityDiagnostic>());
}
} // namespace pamc
