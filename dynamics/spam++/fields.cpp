



  #include "fields.h"


  int YAKL_INLINE get_total_dofs(int ndims, int ndof0, int ndof1, int ndof2, int ndof3)
  {
  if (ndims == 1) { return ndof0 + ndof1; }
  if (ndims == 2) { return ndof0 + 2*ndof1 + ndof2; }
  if (ndims == 3) { return ndof0 + 3*ndof1 + 3*ndof2 + ndof3; }
  }


  // creates a new Field f with same parameters as self, without copying data over
  void Field::clone(Field &f, std::string fieldName) {
    initialize(f.topology, f.fieldName);
  }

  void Field::initialize(Topology &topo, std::string fieldName) {
    topology = topo;
    geometry = geom;
    name = FieldName;

    total_dofs = get_total_dofs(ndims, ndof0, ndof1, ndof2, ndof3);

    data = topology.create_arr(name, &data, ndof0, ndof1, ndof2, ndof3);

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;
    yakl::parallel_for("ZeroField", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int ndof=0, ndof<total_dofs; ndof++) {
        data(ndof, k+ks, j+js, i+is) = 0.0;
      }
    });

  }

    int YAKL_INLINE Field::get_offset(int ndof) {
      if (ndof == 0)                  { return 0; }
      if (ndof == 1)                  { return ndof0; }
      if ((ndof == 2) && (ndim == 2)) { return ndof0 + 2*ndof1; }
      if ((ndof == 2) && (ndim == 3)) { return ndof0 + 3*ndof1; }
      if ((ndof == 3) && (ndim == 3)) { return ndof0 + 3*ndof1 + 3*ndof2; }
    }





  // copies data from f into self

// Ideally ndofs here is a compile time constant...

  void Field::copy(Field & f) {
    yakl::parallel_for("CopyField", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int ndof=0, ndof<total_dofs; ndof++) {
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
      for (int ndof=0, ndof<total_dofs; ndof++) {
        data(ndof, k+ks, j+js, i+is) = alpha * x.data(ndof, k+ks, j+js, i+is) + y.data(ndof, k+ks, j+js, i+is);
      }
    });

  }
