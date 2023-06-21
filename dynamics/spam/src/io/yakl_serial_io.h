#pragma once

#include "YAKL_netcdf.h"
#include "common.h"
#include "field_sets.h"
#include "model.h"
#include "stats.h"

class FileIO {

public:
  bool is_initialized;
  int masterproc;
  yakl::SimpleNetCDF nc;
  int ulIndex = 0; // Unlimited dimension index to place this data at

  std::array<real5d, nprognostic> prog_temp_arr;
  std::array<real5d, nconstant> const_temp_arr;
  std::vector<real5d> diag_temp_arr;
  std::string outputName;

  const FieldSet<nprognostic> *prog_vars;
  const FieldSet<nconstant> *const_vars;
  const std::vector<std::unique_ptr<Diagnostic>> *diagnostics;
  Stats *statistics;

  FileIO();
  FileIO(const FileIO &fio) = delete;
  FileIO &operator=(const FileIO &fio) = delete;
  void initialize(std::string outputName, Topology &ptopo, Topology &dtopo,
                  Parallel &par, const FieldSet<nprognostic> &progvars,
                  const FieldSet<nconstant> &const_vars,
                  const std::vector<std::unique_ptr<Diagnostic>> &diag,
                  Stats &stats);
  void output(real time);
  void outputInit(real time, const Geometry<Straight> &primal_geometry,
                  const Geometry<Twisted> &dual_geometry);
  void outputStats(const Stats &stats);
};

FileIO::FileIO() { this->is_initialized = false; }

void FileIO::initialize(std::string outName, Topology &ptopo, Topology &dtopo,
                        Parallel &par, const FieldSet<nprognostic> &progvars,
                        const FieldSet<nconstant> &constvars,
                        const std::vector<std::unique_ptr<Diagnostic>> &diag,
                        Stats &stats) {

  this->outputName = outName + std::to_string(par.actualrank) + ".nc";
  this->prog_vars = &progvars;
  this->const_vars = &constvars;
  this->diagnostics = &diag;
  this->statistics = &stats;
  this->masterproc = par.masterproc;

  // int nranks;
  // int ierr = MPI_Comm_size(MPI_COMM_WORLD,&nranks);
  // if (nranks > 1) {endrun("spam++ cannot use serial IO in a parallel run");}

  // NEED TO OUTPUT A LOT MORE INFO HERE
  //  start/end indices, density names, etc.
  //  basically all the info that is output to terminal
  //  also coordinate values for various dofs
  //  and dof values for various fields
  //  some of these are field attributes, and some are global attributes

  // field attributes- dof values, density names
  // global attributes- start/end indices, initial cond, dt's, nsteps/out/stat,
  // nranks, nprocx/nprocy, crm_per_phys, etc. basically all the yaml file info;
  // plus some parallel decomp stuff new arrays- coordinate values for all the
  // various dofs

  nc.create(this->outputName);
  nc.createDim("t");
  // nc.createDim( "primal_ncells_x" ,  ptopo.nx_glob );
  // nc.createDim( "primal_ncells_y" ,  ptopo.ny_glob );
  nc.createDim("primal_ncells_x", ptopo.n_cells_x);
  nc.createDim("primal_ncells_y", ptopo.n_cells_y);
  nc.createDim("primal_nlayers", ptopo.nl);
  nc.createDim("primal_ninterfaces", ptopo.ni);
  // nc.createDim( "dual_ncells_x" ,  dtopo.nx_glob );
  // nc.createDim( "dual_ncells_y" ,  dtopo.ny_glob );
  nc.createDim("dual_ncells_x", dtopo.n_cells_x);
  nc.createDim("dual_ncells_y", dtopo.n_cells_y);
  nc.createDim("dual_nlayers", dtopo.nl);
  nc.createDim("dual_ninterfaces", dtopo.ni);
  nc.createDim("nens", ptopo.nens);

  if (this->masterproc) {
    nc.createDim("statsize", this->statistics->statsize);
    for (int l = 0; l < this->statistics->stats_arr.size(); l++) {
      nc.createDim(this->statistics->stats_arr[l].name + "_ndofs",
                   this->statistics->stats_arr[l].ndofs);
    }
  }

  for (int i = 0; i < this->const_vars->fields_arr.size(); i++) {
    nc.createDim(this->const_vars->fields_arr[i].name + "_ndofs",
                 this->const_vars->fields_arr[i].total_dofs);
    this->const_temp_arr[i] =
        real5d(this->const_vars->fields_arr[i].name.c_str(),
               this->const_vars->fields_arr[i].total_dofs,
               this->const_vars->fields_arr[i]._nz,
               this->const_vars->fields_arr[i].topology.n_cells_y,
               this->const_vars->fields_arr[i].topology.n_cells_x,
               this->const_vars->fields_arr[i].topology.nens);
  }

  for (int i = 0; i < this->prog_vars->fields_arr.size(); i++) {
    nc.createDim(this->prog_vars->fields_arr[i].name + "_ndofs",
                 this->prog_vars->fields_arr[i].total_dofs);
    this->prog_temp_arr[i] =
        real5d(this->prog_vars->fields_arr[i].name.c_str(),
               this->prog_vars->fields_arr[i].total_dofs,
               this->prog_vars->fields_arr[i]._nz,
               this->prog_vars->fields_arr[i].topology.n_cells_y,
               this->prog_vars->fields_arr[i].topology.n_cells_x,
               this->prog_vars->fields_arr[i].topology.nens);
  }

  for (auto &diag : *diagnostics) {
    auto &field = diag->field;
    nc.createDim(field.name + "_ndofs", field.total_dofs);
    diag_temp_arr.emplace_back(real5d(field.name.c_str(), field.total_dofs,
                                      field._nz, field.topology.n_cells_y,
                                      field.topology.n_cells_x,
                                      field.topology.nens));
  }

  nc.close();

  this->is_initialized = true;
}

void FileIO::output(real time) {
  nc.open(this->outputName, yakl::NETCDF_MODE_WRITE);
  ulIndex = nc.getDimSize("t");
  // Write the elapsed time
  nc.write1(time, "t", ulIndex, "t");

  for (int l = 0; l < this->prog_vars->fields_arr.size(); l++) {

    int is = this->prog_vars->fields_arr[l].topology.is;
    int js = this->prog_vars->fields_arr[l].topology.js;
    int ks = this->prog_vars->fields_arr[l].topology.ks;
    YAKL_SCOPE(prog_vars_fields, this->prog_vars->fields_arr);
    YAKL_SCOPE(prog_temp_arr, this->prog_temp_arr);
    parallel_for(
        "spam output prognostic",
        SimpleBounds<5>(this->prog_vars->fields_arr[l].total_dofs,
                        this->prog_vars->fields_arr[l]._nz,
                        this->prog_vars->fields_arr[l].topology.n_cells_y,
                        this->prog_vars->fields_arr[l].topology.n_cells_x,
                        this->prog_vars->fields_arr[l].topology.nens),
        YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
          prog_temp_arr[l](ndof, k, j, i, n) =
              prog_vars_fields[l].data(ndof, k + ks, j + js, i + is, n);
        });

    if (this->prog_vars->fields_arr[l].topology.primal) {
      if (this->prog_vars->fields_arr[l].extdof == 1) {
        nc.write1(this->prog_temp_arr[l].createHostCopy(),
                  this->prog_vars->fields_arr[l].name,
                  {this->prog_vars->fields_arr[l].name + "_ndofs",
                   "primal_nlayers", "primal_ncells_y", "primal_ncells_x",
                   "nens"},
                  ulIndex, "t");
      }
      if (this->prog_vars->fields_arr[l].extdof == 0) {
        nc.write1(this->prog_temp_arr[l].createHostCopy(),
                  this->prog_vars->fields_arr[l].name,
                  {this->prog_vars->fields_arr[l].name + "_ndofs",
                   "primal_ninterfaces", "primal_ncells_y", "primal_ncells_x",
                   "nens"},
                  ulIndex, "t");
      }
    } else {
      if (this->prog_vars->fields_arr[l].extdof == 1) {
        nc.write1(this->prog_temp_arr[l].createHostCopy(),
                  this->prog_vars->fields_arr[l].name,
                  {this->prog_vars->fields_arr[l].name + "_ndofs",
                   "dual_nlayers", "dual_ncells_y", "dual_ncells_x", "nens"},
                  ulIndex, "t");
      }
      if (this->prog_vars->fields_arr[l].extdof == 0) {
        nc.write1(this->prog_temp_arr[l].createHostCopy(),
                  this->prog_vars->fields_arr[l].name,
                  {this->prog_vars->fields_arr[l].name + "_ndofs",
                   "dual_ninterfaces", "dual_ncells_y", "dual_ncells_x",
                   "nens"},
                  ulIndex, "t");
      }
    }
  }

  for (int l = 0; l < diagnostics->size(); l++) {
    auto &field = (*diagnostics)[l]->field;

    int is = field.topology.is;
    int js = field.topology.js;
    int ks = field.topology.ks;
    YAKL_SCOPE(diag_temp_arr_l, this->diag_temp_arr[l]);
    parallel_for(
        "spam output diagnostics",
        SimpleBounds<5>(field.total_dofs, field._nz, field.topology.n_cells_y,
                        field.topology.n_cells_x, field.topology.nens),
        YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
          diag_temp_arr_l(ndof, k, j, i, n) =
              field.data(ndof, k + ks, j + js, i + is, n);
        });

    if (field.topology.primal) {
      if (field.extdof == 1) {
        nc.write1(this->diag_temp_arr[l].createHostCopy(), field.name,
                  {field.name + "_ndofs", "primal_nlayers", "primal_ncells_y",
                   "primal_ncells_x", "nens"},
                  ulIndex, "t");
      }
      if (field.extdof == 0) {
        nc.write1(this->diag_temp_arr[l].createHostCopy(), field.name,
                  {field.name + "_ndofs", "primal_ninterfaces",
                   "primal_ncells_y", "primal_ncells_x", "nens"},
                  ulIndex, "t");
      }
    } else {
      if (field.extdof == 1) {
        nc.write1(this->diag_temp_arr[l].createHostCopy(), field.name,
                  {field.name + "_ndofs", "dual_nlayers", "dual_ncells_y",
                   "dual_ncells_x", "nens"},
                  ulIndex, "t");
      }
      if (field.extdof == 0) {
        nc.write1(this->diag_temp_arr[l].createHostCopy(), field.name,
                  {field.name + "_ndofs", "dual_ninterfaces", "dual_ncells_y",
                   "dual_ncells_x", "nens"},
                  ulIndex, "t");
      }
    }
  }

  nc.close();
}

void FileIO::outputInit(real time, const Geometry<Straight> &primal_geometry,
                        const Geometry<Twisted> &dual_geometry) {
  nc.open(this->outputName, yakl::NETCDF_MODE_WRITE);

  nc.write(primal_geometry.dx, "dx");
  nc.write(primal_geometry.Lx, "Lx");
  nc.write(primal_geometry.xc, "xc");
  if (ndims > 1) {
    nc.write(primal_geometry.dy, "dy");
    nc.write(primal_geometry.Ly, "Ly");
    nc.write(primal_geometry.yc, "yc");
  }

#ifdef _EXTRUDED
  {
    const int pks = primal_geometry.topology.ks;

    // primal zint
    real2d primal_zint("primal_zint", primal_geometry.topology.ni,
                       primal_geometry.topology.nens);
    parallel_for(
        "output primal zint",
        SimpleBounds<2>(primal_geometry.topology.ni,
                        primal_geometry.topology.nens),
        YAKL_LAMBDA(int k, int n) {
          primal_zint(k, n) = primal_geometry.zint(k + pks, n);
        });
    nc.write(primal_zint.createHostCopy(), "primal_zint",
             {"primal_ninterfaces", "nens"});

    // primal dz
    real2d primal_dz("primal_dz", primal_geometry.topology.nl,
                     primal_geometry.topology.nens);
    parallel_for(
        "output primal dz",
        SimpleBounds<2>(primal_geometry.topology.nl,
                        primal_geometry.topology.nens),
        YAKL_LAMBDA(int k, int n) {
          primal_dz(k, n) = primal_geometry.dz(k + pks, n);
        });
    nc.write(primal_dz.createHostCopy(), "primal_dz",
             {"primal_nlayers", "nens"});

    const int dks = dual_geometry.topology.ks;

    // dual zint
    real2d dual_zint("dual_zint", dual_geometry.topology.ni,
                     dual_geometry.topology.nens);
    parallel_for(
        "output dual zint",
        SimpleBounds<2>(dual_geometry.topology.ni, dual_geometry.topology.nens),
        YAKL_LAMBDA(int k, int n) {
          dual_zint(k, n) = dual_geometry.zint(k + dks, n);
        });
    nc.write(dual_zint.createHostCopy(), "dual_zint",
             {"dual_ninterfaces", "nens"});

    // dual dz
    real2d dual_dz("dual_dz", dual_geometry.topology.nl,
                   dual_geometry.topology.nens);
    parallel_for(
        "output dual dz",
        SimpleBounds<2>(dual_geometry.topology.nl, dual_geometry.topology.nens),
        YAKL_LAMBDA(int k, int n) {
          dual_dz(k, n) = dual_geometry.dz(k + dks, n);
        });
    nc.write(dual_dz.createHostCopy(), "dual_dz", {"dual_nlayers", "nens"});
  }
#endif

  for (int l = 0; l < this->const_vars->fields_arr.size(); l++) {

    int is = this->const_vars->fields_arr[l].topology.is;
    int js = this->const_vars->fields_arr[l].topology.js;
    int ks = this->const_vars->fields_arr[l].topology.ks;

    YAKL_SCOPE(const_vars_fields, this->const_vars->fields_arr);
    YAKL_SCOPE(const_temp_arr, this->const_temp_arr);
    parallel_for(
        "spam output constant",
        SimpleBounds<5>(this->const_vars->fields_arr[l].total_dofs,
                        this->const_vars->fields_arr[l]._nz,
                        this->const_vars->fields_arr[l].topology.n_cells_y,
                        this->const_vars->fields_arr[l].topology.n_cells_x,
                        this->const_vars->fields_arr[l].topology.nens),
        YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
          const_temp_arr[l](ndof, k, j, i, n) =
              const_vars_fields[l].data(ndof, k + ks, j + js, i + is, n);
        });

    if (this->const_vars->fields_arr[l].topology.primal) {
      if (this->const_vars->fields_arr[l].extdof == 1) {
        nc.write(this->const_temp_arr[l].createHostCopy(),
                 this->const_vars->fields_arr[l].name,
                 {this->const_vars->fields_arr[l].name + "_ndofs",
                  "primal_nlayers", "primal_ncells_y", "primal_ncells_x",
                  "nens"});
      }
      if (this->const_vars->fields_arr[l].extdof == 0) {
        nc.write(this->const_temp_arr[l].createHostCopy(),
                 this->const_vars->fields_arr[l].name,
                 {this->const_vars->fields_arr[l].name + "_ndofs",
                  "primal_ninterfaces", "primal_ncells_y", "primal_ncells_x",
                  "nens"});
      }
    } else {
      if (this->const_vars->fields_arr[l].extdof == 1) {
        nc.write(this->const_temp_arr[l].createHostCopy(),
                 this->const_vars->fields_arr[l].name,
                 {this->const_vars->fields_arr[l].name + "_ndofs",
                  "dual_nlayers", "dual_ncells_y", "dual_ncells_x", "nens"});
      }
      if (this->const_vars->fields_arr[l].extdof == 0) {
        nc.write(this->const_temp_arr[l].createHostCopy(),
                 this->const_vars->fields_arr[l].name,
                 {this->const_vars->fields_arr[l].name + "_ndofs",
                  "dual_ninterfaces", "dual_ncells_y", "dual_ncells_x",
                  "nens"});
      }
    }
  }

  nc.close();

  output(time);
}

void FileIO::outputStats(const Stats &stats) {
  if (this->masterproc) {
    nc.open(this->outputName, yakl::NETCDF_MODE_WRITE);

    for (int l = 0; l < this->statistics->stats_arr.size(); l++) {
      nc.write(
          this->statistics->stats_arr[l].data,
          this->statistics->stats_arr[l].name,
          {this->statistics->stats_arr[l].name + "_ndofs", "statsize", "nens"});
    }

    nc.close();
  }
}
