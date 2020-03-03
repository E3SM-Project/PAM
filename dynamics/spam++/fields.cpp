



  #include "fields.h"


  void Field::initialize(Topology &topo, Exchange &exch, std::string fieldName, int ndof0, int ndof1, int ndof2 = 0, int ndof3 = 0) {
    topology = topo;
    topology = exch;
    geometry = geom;
    name = FieldName;
    ndofs0 = ndof0;
    ndofs1 = ndof1;
    ndofs2 = ndof2;
    ndofs3 = ndof3;

    if (ndims == 1) { total_dofs = ndofs0 + ndofs1; }
    if (ndims == 2) { total_dofs = ndofs0 + 2*ndofs1 + ndofs2; }
    if (ndims == 3) { total_dofs = ndofs0 + 3*ndofs1 + 3*ndofs2 + ndofs3; }

    data = topology.create_field(name, &data, ndofs0, ndofs1, ndofs2, ndofs3);

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;
    yakl::parallel_for("ZeroField", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (ndof=0, ndof<total_dofs; ndof++) {
        data(ndof, k+ks, j+js, i+is) = 0.0;
      }
    });

  }

    int YAKL_INLINE Field::get_offset(int ndof) {
      if (ndof == 0)                  { return 0; }
      if (ndof == 1)                  { return ndofs0; }
      if ((ndof == 2) && (ndim == 2)) { return ndofs0 + 2*ndofs1; }
      if ((ndof == 2) && (ndim == 3)) { return ndofs0 + 3*ndofs1; }
      if ((ndof == 3) && (ndim == 3)) { return ndofs0 + 3*ndofs1 + 3*ndofs2; }
    }


  // creates a new Field f that is a deep clone of self, without copying data over
  void Field::clone(Field &f, std::string fieldName) {
    f.initialize(topology, fieldName, ndofs0, ndofs1, ndofs2, ndofs3);
  }

  // copies data from f into self

// Ideally ndofs here is a compile time constant...

  void Field::copy(Field & f) {
    yakl::parallel_for("CopyField", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (ndof=0, ndof<total_dofs; ndof++) {
        f.data(ndof, k+ks, j+js, i+is) = data(ndof, k+ks, j+js, i+is);
      }
    });
  }

// Ideally ndofs here is a compile time constant...
  // Computes w (self) = alpha x + y
  void Field::waxpy(real alpha, Field &x, Field &y) {

    yakl::parallel_for("WAXPYField", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (ndof=0, ndof<total_dofs; ndof++) {
        data(ndof, k+ks, j+js, i+is) = alpha * x.data(ndof, k+ks, j+js, i+is) + y.data(ndof, k+ks, j+js, i+is);
      }
    });

  }
