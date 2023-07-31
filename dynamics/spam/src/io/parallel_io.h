
#pragma once

#include "common.h"
#include "field_sets.h"
#include "pnetcdf.h"

#include "common.h"
#include "stats.h"

namespace pamc {

// Error reporting routine for the PNetCDF I/O
void ncwrap(int ierr, int line) {
  if (ierr != NC_NOERR) {
    std::cout << "NetCDF Error at line: " << line << "\n";
    std::cout << ncmpi_strerror(ierr) << "\n";
    exit(-1);
  }
}

template <uint nprog, uint nconst, uint ndiag, uint nstats> class FileIO {

public:
  bool is_initialized;
  int masterproc;
  int ncid;
  int const_dim_ids[5], prog_dim_ids[6], diag_dim_ids[6];
  int tDim, pxDim, pyDim, plDim, piDim, dxDim, dyDim, dlDim, diDim, nensDim;
  int const_var_dim_ids[nconst], prog_var_dim_ids[nprog],
      diag_var_dim_ids[ndiag];
  int const_var_ids[nconst], prog_var_ids[nprog], diag_var_ids[ndiag];
  int tVar;
  std::string outputName;
  std::array<real5d, nprog> prog_temp_arr;
  std::array<real5d, nconst> const_temp_arr;
  std::array<real5d, ndiag> diag_temp_arr;
  const FieldSet<nprog> *prog_vars;
  const FieldSet<nconst> *const_vars;
  const FieldSet<ndiag> *diag_vars;
  Stats<nprog, nconst, nstats> *statistics;

  int stat_var_dim_ids[nstats];
  int statDimIds[3];
  int stat_ids[nstats];

  int numOut;

  MPI_Offset const_start[5], const_count[5];
  MPI_Offset prog_start[6], prog_count[6];
  MPI_Offset diag_start[6], diag_count[6];

  FileIO();
  FileIO(const FileIO<nprog, nconst, ndiag, nstats> &fio) = delete;
  FileIO &operator=(const FileIO<nprog, nconst, ndiag, nstats> &fio) = delete;
  void initialize(std::string outputName, Topology &ptopo, Topology &dtopo,
                  Parallel &par, const FieldSet<nprog> &progvars,
                  const FieldSet<nconst> &const_vars,
                  const FieldSet<ndiag> &diagvars,
                  Stats<nprog, nconst, nstats> &stats);
  void output(int nstep, real time);
  void outputInit(real time);
  void outputStats(const Stats<nprog, nconst, nstats> &stats);
  void close();
};

template <uint nprog, uint nconst, uint ndiag, uint nstats>
FileIO<nprog, nconst, ndiag, nstats>::FileIO() {
  this->is_initialized = false;
}

template <uint nprog, uint nconst, uint ndiag, uint nstats>
void FileIO<nprog, nconst, ndiag, nstats>::initialize(
    std::string outName, Topology &ptopo, Topology &dtopo, Parallel &par,
    const FieldSet<nprog> &progvars, const FieldSet<nconst> &constvars,
    const FieldSet<ndiag> &diagvars, Stats<nprog, nconst, nstats> &stats) {
  this->outputName = outName;
  this->prog_vars = &progvars;
  this->const_vars = &constvars;
  this->diag_vars = &diagvars;
  this->statistics = &stats;
  this->masterproc = par.masterproc;

  ncwrap(ncmpi_create(MPI_COMM_WORLD, this->outputName.c_str(), NC_CLOBBER,
                      MPI_INFO_NULL, &ncid),
         __LINE__);

  // Define time and cell dimensions
  ncwrap(ncmpi_def_dim(ncid, "t", (MPI_Offset)NC_UNLIMITED, &tDim), __LINE__);
  ncwrap(
      ncmpi_def_dim(ncid, "primal_ncells_x", (MPI_Offset)ptopo.nx_glob, &pxDim),
      __LINE__);
  ncwrap(
      ncmpi_def_dim(ncid, "primal_ncells_y", (MPI_Offset)ptopo.ny_glob, &pyDim),
      __LINE__);
  ncwrap(ncmpi_def_dim(ncid, "primal_nlayers", (MPI_Offset)ptopo.nl, &plDim),
         __LINE__);
  ncwrap(
      ncmpi_def_dim(ncid, "primal_ninterfaces", (MPI_Offset)ptopo.ni, &piDim),
      __LINE__);
  ncwrap(
      ncmpi_def_dim(ncid, "dual_ncells_x", (MPI_Offset)dtopo.nx_glob, &dxDim),
      __LINE__);
  ncwrap(
      ncmpi_def_dim(ncid, "dual_ncells_y", (MPI_Offset)dtopo.ny_glob, &dyDim),
      __LINE__);
  ncwrap(ncmpi_def_dim(ncid, "dual_nlayers", (MPI_Offset)dtopo.nl, &dlDim),
         __LINE__);
  ncwrap(ncmpi_def_dim(ncid, "dual_ninterfaces", (MPI_Offset)dtopo.ni, &diDim),
         __LINE__);
  ncwrap(ncmpi_def_dim(ncid, "nens", (MPI_Offset)ptopo.nens, &nensDim),
         __LINE__);

  // Create time dimension
  const_dim_ids[0] = tDim;
  ncwrap(ncmpi_def_var(ncid, "t", REAL_NC, 1, const_dim_ids, &tVar), __LINE__);
  ncwrap(ncmpi_put_att_text(ncid, tVar, "units", std::strlen("seconds"),
                            "seconds"),
         __LINE__);

  for (int i = 0; i < this->const_vars->fields_arr.size(); i++) {
    if (this->const_vars->fields_arr[i].total_dofs > 0) {
      ncwrap(ncmpi_def_dim(
                 ncid,
                 (this->const_vars->fields_arr[i].name + "_ndofs").c_str(),
                 (MPI_Offset)this->const_vars->fields_arr[i].total_dofs,
                 &const_var_dim_ids[i]),
             __LINE__);
      if (this->const_vars->fields_arr[i].topology->primal) {
        if (this->const_vars->fields_arr[i].extdof == 1) {
          const_dim_ids[0] = const_var_dim_ids[i];
          const_dim_ids[1] = plDim;
          const_dim_ids[2] = pyDim;
          const_dim_ids[3] = pxDim;
          const_dim_ids[4] = nensDim;
        }
        if (this->const_vars->fields_arr[i].extdof == 0) {
          const_dim_ids[0] = const_var_dim_ids[i];
          const_dim_ids[1] = piDim;
          const_dim_ids[2] = pyDim;
          const_dim_ids[3] = pxDim;
          const_dim_ids[4] = nensDim;
        }
      } else {
        if (this->const_vars->fields_arr[i].extdof == 1) {
          const_dim_ids[0] = const_var_dim_ids[i];
          const_dim_ids[1] = dlDim;
          const_dim_ids[2] = dyDim;
          const_dim_ids[3] = dxDim;
          const_dim_ids[4] = nensDim;
        }
        if (this->const_vars->fields_arr[i].extdof == 0) {
          const_dim_ids[0] = const_var_dim_ids[i];
          const_dim_ids[1] = diDim;
          const_dim_ids[2] = dyDim;
          const_dim_ids[3] = dxDim;
          const_dim_ids[4] = nensDim;
        }
      }
      ncwrap(ncmpi_def_var(ncid, this->const_vars->fields_arr[i].name.c_str(),
                           REAL_NC, 5, const_dim_ids, &const_var_ids[i]),
             __LINE__);
      this->const_temp_arr[i] =
          real5d(this->const_vars->fields_arr[i].name.c_str(),
                 this->const_vars->fields_arr[i].total_dofs,
                 this->const_vars->fields_arr[i]._nz,
                 this->const_vars->fields_arr[i].topology->n_cells_y,
                 this->const_vars->fields_arr[i].topology->n_cells_x,
                 this->const_vars->fields_arr[i].topology->nens);
    }
  }

  for (int i = 0; i < this->prog_vars->fields_arr.size(); i++) {
    if (this->prog_vars->fields_arr[i].total_dofs > 0) {
      ncwrap(ncmpi_def_dim(
                 ncid, (this->prog_vars->fields_arr[i].name + "_ndofs").c_str(),
                 (MPI_Offset)this->prog_vars->fields_arr[i].total_dofs,
                 &prog_var_dim_ids[i]),
             __LINE__);
      if (this->prog_vars->fields_arr[i].topology->primal) {
        if (this->prog_vars->fields_arr[i].extdof == 1) {
          prog_dim_ids[0] = tVar;
          prog_dim_ids[1] = prog_var_dim_ids[i];
          prog_dim_ids[2] = plDim;
          prog_dim_ids[3] = pyDim;
          prog_dim_ids[4] = pxDim;
          prog_dim_ids[5] = nensDim;
        }
        if (this->prog_vars->fields_arr[i].extdof == 0) {
          prog_dim_ids[0] = tVar;
          prog_dim_ids[1] = prog_var_dim_ids[i];
          prog_dim_ids[2] = piDim;
          prog_dim_ids[3] = pyDim;
          prog_dim_ids[4] = pxDim;
          prog_dim_ids[5] = nensDim;
        }
      } else {
        if (this->prog_vars->fields_arr[i].extdof == 1) {
          prog_dim_ids[0] = tVar;
          prog_dim_ids[1] = prog_var_dim_ids[i];
          prog_dim_ids[2] = dlDim;
          prog_dim_ids[3] = dyDim;
          prog_dim_ids[4] = dxDim;
          prog_dim_ids[5] = nensDim;
        }
        if (this->prog_vars->fields_arr[i].extdof == 0) {
          prog_dim_ids[0] = tVar;
          prog_dim_ids[1] = prog_var_dim_ids[i];
          prog_dim_ids[2] = diDim;
          prog_dim_ids[3] = dyDim;
          prog_dim_ids[4] = dxDim;
          prog_dim_ids[5] = nensDim;
        }
      }
      ncwrap(ncmpi_def_var(ncid, this->prog_vars->fields_arr[i].name.c_str(),
                           REAL_NC, 6, prog_dim_ids, &prog_var_ids[i]),
             __LINE__);
      this->prog_temp_arr[i] =
          real5d(this->prog_vars->fields_arr[i].name.c_str(),
                 this->prog_vars->fields_arr[i].total_dofs,
                 this->prog_vars->fields_arr[i]._nz,
                 this->prog_vars->fields_arr[i].topology->n_cells_y,
                 this->prog_vars->fields_arr[i].topology->n_cells_x,
                 this->prog_vars->fields_arr[i].topology->nens);
    }
  }

  for (int i = 0; i < this->diag_vars->fields_arr.size(); i++) {
    if (this->diag_vars->fields_arr[i].total_dofs > 0) {
      ncwrap(ncmpi_def_dim(
                 ncid, (this->diag_vars->fields_arr[i].name + "_ndofs").c_str(),
                 (MPI_Offset)this->diag_vars->fields_arr[i].total_dofs,
                 &diag_var_dim_ids[i]),
             __LINE__);
      if (this->diag_vars->fields_arr[i].topology->primal) {
        if (this->diag_vars->fields_arr[i].extdof == 1) {
          diag_dim_ids[0] = tVar;
          diag_dim_ids[1] = diag_var_dim_ids[i];
          diag_dim_ids[2] = plDim;
          diag_dim_ids[3] = pyDim;
          diag_dim_ids[4] = pxDim;
          diag_dim_ids[5] = nensDim;
        }
        if (this->diag_vars->fields_arr[i].extdof == 0) {
          diag_dim_ids[0] = tVar;
          diag_dim_ids[1] = diag_var_dim_ids[i];
          diag_dim_ids[2] = piDim;
          diag_dim_ids[3] = pyDim;
          diag_dim_ids[4] = pxDim;
          diag_dim_ids[5] = nensDim;
        }
      } else {
        if (this->diag_vars->fields_arr[i].extdof == 1) {
          diag_dim_ids[0] = tVar;
          diag_dim_ids[1] = diag_var_dim_ids[i];
          diag_dim_ids[2] = dlDim;
          diag_dim_ids[3] = dyDim;
          diag_dim_ids[4] = dxDim;
          diag_dim_ids[5] = nensDim;
        }
        if (this->diag_vars->fields_arr[i].extdof == 0) {
          diag_dim_ids[0] = tVar;
          diag_dim_ids[1] = diag_var_dim_ids[i];
          diag_dim_ids[2] = diDim;
          diag_dim_ids[3] = dyDim;
          diag_dim_ids[4] = dxDim;
          diag_dim_ids[5] = nensDim;
        }
      }
      ncwrap(ncmpi_def_var(ncid, this->diag_vars->fields_arr[i].name.c_str(),
                           REAL_NC, 6, diag_dim_ids, &diag_var_ids[i]),
             __LINE__);
      this->diag_temp_arr[i] =
          real5d(this->diag_vars->fields_arr[i].name.c_str(),
                 this->diag_vars->fields_arr[i].total_dofs,
                 this->diag_vars->fields_arr[i]._nz,
                 this->diag_vars->fields_arr[i].topology->n_cells_y,
                 this->diag_vars->fields_arr[i].topology->n_cells_x,
                 this->diag_vars->fields_arr[i].topology->nens);
    }
  }

  // define statistics dimension and variables

  ncwrap(ncmpi_def_dim(ncid, "nstat", (MPI_Offset)this->statistics->statsize,
                       &statDimIds[1]),
         __LINE__);
  statDimIds[2] = nensDim;

  for (int l = 0; l < this->statistics->stats_arr.size(); l++) {
    if (this->statistics->stats_arr[l].ndofs > 0) {
      ncwrap(ncmpi_def_dim(
                 ncid, (this->statistics->stats_arr[l].name + "_ndofs").c_str(),
                 (MPI_Offset)this->statistics->stats_arr[l].ndofs,
                 &stat_var_dim_ids[l]),
             __LINE__);
      statDimIds[0] = stat_var_dim_ids[l];
      ncwrap(ncmpi_def_var(ncid, this->statistics->stats_arr[l].name.c_str(),
                           REAL_NC, 3, statDimIds, &stat_ids[l]),
             __LINE__);
    }
  }

  ncwrap(ncmpi_enddef(ncid), __LINE__);
  ncwrap(ncmpi_close(ncid), __LINE__);
  this->is_initialized = true;
}

template <uint nprog, uint nconst, uint ndiag, uint nstats>
void FileIO<nprog, nconst, ndiag, nstats>::output(int nstep, real time) {

  ncwrap(ncmpi_open(MPI_COMM_WORLD, this->outputName.c_str(), NC_WRITE,
                    MPI_INFO_NULL, &ncid),
         __LINE__);

  for (int l = 0; l < this->prog_vars->fields_arr.size(); l++) {
    if (this->prog_vars->fields_arr[l].total_dofs > 0) {
      ncwrap(ncmpi_inq_varid(ncid, this->prog_vars->fields_arr[l].name.c_str(),
                             &prog_var_ids[l]),
             __LINE__);
      int is = this->prog_vars->fields_arr[l].topology->is;
      int js = this->prog_vars->fields_arr[l].topology->js;
      int ks = this->prog_vars->fields_arr[l].topology->ks;
      // yakl::parallel_for("CopyFieldToOutputBuffer",
      // this->prog_vars->fields_arr[l]._nloop, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, this->prog_vars->fields_arr[l]._nz,
      //   this->prog_vars->fields_arr[l].topology->n_cells_y,
      //   this->prog_vars->fields_arr[l].topology->n_cells_x, k, j, i); for
      //   (int ndof=0; ndof<this->prog_vars->fields_arr[l].total_dofs; ndof++)
      //   {
      parallel_for(
          SimpleBounds<5>(this->prog_vars->fields_arr[l].total_dofs,
                          this->prog_vars->fields_arr[l]._nz,
                          this->prog_vars->fields_arr[l].topology->n_cells_y,
                          this->prog_vars->fields_arr[l].topology->n_cells_x,
                          this->prog_vars->fields_arr[l].topology->nens),
          YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
            this->prog_temp_arr[l](ndof, k, j, i, n) =
                this->prog_vars->fields_arr[l].data(ndof, k + ks, j + js,
                                                    i + is, n);
            //  }
          });

      prog_start[0] = this->numOut;
      prog_start[1] = 0;
      prog_start[2] = 0;
      prog_start[3] = this->prog_vars->fields_arr[l].topology->j_beg;
      prog_start[4] = this->prog_vars->fields_arr[l].topology->i_beg;
      prog_start[5] = 0;
      prog_count[0] = 1;
      prog_count[1] = this->prog_vars->fields_arr[l].total_dofs;
      prog_count[2] = this->prog_vars->fields_arr[l]._nz;
      prog_count[3] = this->prog_vars->fields_arr[l].topology->n_cells_y;
      prog_count[4] = this->prog_vars->fields_arr[l].topology->n_cells_x;
      prog_count[4] = this->prog_vars->fields_arr[l].topology->nens;
      ncwrap(
          PNETCDF_PUT_VAR_ALL(ncid, prog_var_ids[l], prog_start, prog_count,
                              this->prog_temp_arr[l].createHostCopy().data()),
          __LINE__);
    }
  }

  for (int l = 0; l < this->diag_vars->fields_arr.size(); l++) {
    if (this->diag_vars->fields_arr[l].total_dofs > 0) {
      ncwrap(ncmpi_inq_varid(ncid, this->diag_vars->fields_arr[l].name.c_str(),
                             &diag_var_ids[l]),
             __LINE__);
      int is = this->diag_vars->fields_arr[l].topology->is;
      int js = this->diag_vars->fields_arr[l].topology->js;
      int ks = this->diag_vars->fields_arr[l].topology->ks;
      // yakl::parallel_for("CopyFieldToOutputBuffer",
      // this->diag_vars->fields_arr[l]._nloop, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, this->diag_vars->fields_arr[l]._nz,
      //   this->diag_vars->fields_arr[l].topology->n_cells_y,
      //   this->diag_vars->fields_arr[l].topology->n_cells_x, k, j, i); for
      //   (int ndof=0; ndof<this->diag_vars->fields_arr[l].total_dofs; ndof++)
      //   {
      parallel_for(
          SimpleBounds<5>(this->diag_vars->fields_arr[l].total_dofs,
                          this->diag_vars->fields_arr[l]._nz,
                          this->diag_vars->fields_arr[l].topology->n_cells_y,
                          this->diag_vars->fields_arr[l].topology->n_cells_x,
                          this->diag_vars->fields_arr[l].topology->nens),
          YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
            this->diag_temp_arr[l](ndof, k, j, i, n) =
                this->diag_vars->fields_arr[l].data(ndof, k + ks, j + js,
                                                    i + is, n);
            //}
          });

      diag_start[0] = this->numOut;
      diag_start[1] = 0;
      diag_start[2] = 0;
      diag_start[3] = this->diag_vars->fields_arr[l].topology->j_beg;
      diag_start[4] = this->diag_vars->fields_arr[l].topology->i_beg;
      diag_start[5] = 0;
      diag_count[0] = 1;
      diag_count[1] = this->diag_vars->fields_arr[l].total_dofs;
      diag_count[2] = this->diag_vars->fields_arr[l]._nz;
      diag_count[3] = this->diag_vars->fields_arr[l].topology->n_cells_y;
      diag_count[4] = this->diag_vars->fields_arr[l].topology->n_cells_x;
      diag_count[5] = this->diag_vars->fields_arr[l].topology->nens;
      ncwrap(
          PNETCDF_PUT_VAR_ALL(ncid, diag_var_ids[l], diag_start, diag_count,
                              this->diag_temp_arr[l].createHostCopy().data()),
          __LINE__);
    }
  }

  // ADD THIS
  // write elapsed time/time step

  ncwrap(ncmpi_close(ncid), __LINE__);
  this->numOut++;
}

template <uint nprog, uint nconst, uint ndiag, uint nstats>
void FileIO<nprog, nconst, ndiag, nstats>::outputInit(real time) {
  ncwrap(ncmpi_open(MPI_COMM_WORLD, this->outputName.c_str(), NC_WRITE,
                    MPI_INFO_NULL, &ncid),
         __LINE__);

  // std::cout << "writing initial fields\n"  <<std::flush;

  for (int l = 0; l < this->const_vars->fields_arr.size(); l++) {
    // std::cout << "writing initial field "  <<
    // this->const_vars->fields_arr[l].name.c_str() << "\n" <<std::flush;
    if (this->const_vars->fields_arr[l].total_dofs > 0) {
      ncwrap(ncmpi_inq_varid(ncid, this->const_vars->fields_arr[l].name.c_str(),
                             &const_var_ids[l]),
             __LINE__);
      int is = this->const_vars->fields_arr[l].topology->is;
      int js = this->const_vars->fields_arr[l].topology->js;
      int ks = this->const_vars->fields_arr[l].topology->ks;
      // yakl::parallel_for("CopyFieldToOutputBuffer",
      // this->const_vars->fields_arr[l]._nloop, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, this->const_vars->fields_arr[l]._nz,
      //   this->const_vars->fields_arr[l].topology->n_cells_y,
      //   this->const_vars->fields_arr[l].topology->n_cells_x, k, j, i); for
      //   (int ndof=0; ndof<this->const_vars->fields_arr[l].total_dofs; ndof++)
      //   {
      parallel_for(
          SimpleBounds<5>(this->const_vars->fields_arr[l].total_dofs,
                          this->const_vars->fields_arr[l]._nz,
                          this->const_vars->fields_arr[l].topology->n_cells_y,
                          this->const_vars->fields_arr[l].topology->n_cells_x,
                          this->const_vars->fields_arr[l].topology->nens),
          YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
            this->const_temp_arr[l](ndof, k, j, i, n) =
                this->const_vars->fields_arr[l].data(ndof, k + ks, j + js,
                                                     i + is, n);
            //  }
          });

      const_start[0] = 0;
      const_start[1] = 0;
      const_start[2] = this->const_vars->fields_arr[l].topology->j_beg;
      const_start[3] = this->const_vars->fields_arr[l].topology->i_beg;
      const_start[4] = 0;
      const_count[0] = this->const_vars->fields_arr[l].total_dofs;
      const_count[1] = this->const_vars->fields_arr[l]._nz;
      const_count[2] = this->const_vars->fields_arr[l].topology->n_cells_y;
      const_count[3] = this->const_vars->fields_arr[l].topology->n_cells_x;
      const_count[4] = this->const_vars->fields_arr[l].topology->nens;
      ncwrap(
          PNETCDF_PUT_VAR_ALL(ncid, const_var_ids[l], const_start, const_count,
                              this->const_temp_arr[l].createHostCopy().data()),
          __LINE__);
    }
  }

  for (int l = 0; l < this->prog_vars->fields_arr.size(); l++) {
    if (this->prog_vars->fields_arr[l].total_dofs > 0) {
      // std::cout << "writing initial field "  <<
      // this->prog_vars->fields_arr[l].name.c_str() << "\n" <<std::flush;
      ncwrap(ncmpi_inq_varid(ncid, this->prog_vars->fields_arr[l].name.c_str(),
                             &prog_var_ids[l]),
             __LINE__);
      int is = this->prog_vars->fields_arr[l].topology->is;
      int js = this->prog_vars->fields_arr[l].topology->js;
      int ks = this->prog_vars->fields_arr[l].topology->ks;
      // yakl::parallel_for("CopyFieldToOutputBuffer",
      // this->prog_vars->fields_arr[l]._nloop, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, this->prog_vars->fields_arr[l]._nz,
      //   this->prog_vars->fields_arr[l].topology->n_cells_y,
      //   this->prog_vars->fields_arr[l].topology->n_cells_x, k, j, i); for
      //   (int ndof=0; ndof<this->prog_vars->fields_arr[l].total_dofs; ndof++)
      //   {
      parallel_for(
          SimpleBounds<5>(this->prog_vars->fields_arr[l].total_dofs,
                          this->prog_vars->fields_arr[l]._nz,
                          this->prog_vars->fields_arr[l].topology->n_cells_y,
                          this->prog_vars->fields_arr[l].topology->n_cells_x,
                          this->prog_vars->fields_arr[l].topology->nens),
          YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
            this->prog_temp_arr[l](ndof, k, j, i, n) =
                this->prog_vars->fields_arr[l].data(ndof, k + ks, j + js,
                                                    i + is, n);
            //  }
          });

      prog_start[0] = this->numOut;
      prog_start[1] = 0;
      prog_start[2] = 0;
      prog_start[3] = this->prog_vars->fields_arr[l].topology->j_beg;
      prog_start[4] = this->prog_vars->fields_arr[l].topology->i_beg;
      prog_start[5] = 0;
      prog_count[0] = 1;
      prog_count[1] = this->prog_vars->fields_arr[l].total_dofs;
      prog_count[2] = this->prog_vars->fields_arr[l]._nz;
      prog_count[3] = this->prog_vars->fields_arr[l].topology->n_cells_y;
      prog_count[4] = this->prog_vars->fields_arr[l].topology->n_cells_x;
      prog_count[5] = this->prog_vars->fields_arr[l].topology->nens;
      ncwrap(
          PNETCDF_PUT_VAR_ALL(ncid, prog_var_ids[l], prog_start, prog_count,
                              this->prog_temp_arr[l].createHostCopy().data()),
          __LINE__);
    }
  }

  for (int l = 0; l < this->diag_vars->fields_arr.size(); l++) {
    if (this->diag_vars->fields_arr[l].total_dofs > 0) {
      // std::cout << "writing initial field "  <<
      // this->diag_vars->fields_arr[l].name.c_str() << "\n" <<std::flush;
      ncwrap(ncmpi_inq_varid(ncid, this->diag_vars->fields_arr[l].name.c_str(),
                             &diag_var_ids[l]),
             __LINE__);
      int is = this->diag_vars->fields_arr[l].topology->is;
      int js = this->diag_vars->fields_arr[l].topology->js;
      int ks = this->diag_vars->fields_arr[l].topology->ks;
      // yakl::parallel_for("CopyFieldToOutputBuffer",
      // this->diag_vars->fields_arr[l]._nloop, YAKL_LAMBDA (int iGlob) {
      //   int k, j, i;
      //   yakl::unpackIndices(iGlob, this->diag_vars->fields_arr[l]._nz,
      //   this->diag_vars->fields_arr[l].topology->n_cells_y,
      //   this->diag_vars->fields_arr[l].topology->n_cells_x, k, j, i); for
      //   (int ndof=0; ndof<this->diag_vars->fields_arr[l].total_dofs; ndof++)
      //   {
      parallel_for(
          SimpleBounds<5>(this->diag_vars->fields_arr[l].total_dofs,
                          this->diag_vars->fields_arr[l]._nz,
                          this->diag_vars->fields_arr[l].topology->n_cells_y,
                          this->diag_vars->fields_arr[l].topology->n_cells_x,
                          this->diag_vars->fields_arr[l].topology->nens),
          YAKL_LAMBDA(int ndof, int k, int j, int i, int n) {
            this->diag_temp_arr[l](ndof, k, j, i, n) =
                this->diag_vars->fields_arr[l].data(ndof, k + ks, j + js,
                                                    i + is, n);
            //  }
          });

      diag_start[0] = this->numOut;
      diag_start[1] = 0;
      diag_start[2] = 0;
      diag_start[3] = this->diag_vars->fields_arr[l].topology->j_beg;
      diag_start[4] = this->diag_vars->fields_arr[l].topology->i_beg;
      diag_start[5] = 0;
      diag_count[0] = 1;
      diag_count[1] = this->diag_vars->fields_arr[l].total_dofs;
      diag_count[2] = this->diag_vars->fields_arr[l]._nz;
      diag_count[3] = this->diag_vars->fields_arr[l].topology->n_cells_y;
      diag_count[4] = this->diag_vars->fields_arr[l].topology->n_cells_x;
      diag_count[5] = this->diag_vars->fields_arr[l].topology->nens;
      ncwrap(
          PNETCDF_PUT_VAR_ALL(ncid, diag_var_ids[l], diag_start, diag_count,
                              this->diag_temp_arr[l].createHostCopy().data()),
          __LINE__);
    }
  }

  // ADD THIS
  // write elapsed time/time step

  ncwrap(ncmpi_close(ncid), __LINE__);

  this->numOut++;
}

template <uint nprog, uint nconst, uint ndiag, uint nstats>
void FileIO<nprog, nconst, ndiag, nstats>::close() {}

template <uint nprog, uint nconst, uint ndiag, uint nstats>
void FileIO<nprog, nconst, ndiag, nstats>::outputStats(
    const Stats<nprog, nconst, nstats> &stats) {
  ncwrap(ncmpi_open(MPI_COMM_WORLD, this->outputName.c_str(), NC_WRITE,
                    MPI_INFO_NULL, &ncid),
         __LINE__);
  ncwrap(ncmpi_begin_indep_data(ncid), __LINE__);

  if (masterproc) {
    MPI_Offset statStart[3], statCount[3];
    statStart[0] = 0;
    statStart[1] = 0;
    statStart[2] = 0;
    statCount[1] = this->statistics->statsize;
    statCount[2] = this->statistics->nens;
    for (int l = 0; l < this->statistics->stats_arr.size(); l++) {
      if (this->statistics->stats_arr[l].ndofs > 0) {
        statCount[0] = this->statistics->stats_arr[l].ndofs;
        ncwrap(ncmpi_inq_varid(ncid,
                               this->statistics->stats_arr[l].name.c_str(),
                               &stat_ids[l]),
               __LINE__);
        ncwrap(PNETCDF_PUT_VAR(
                   ncid, stat_ids[l], statStart, statCount,
                   this->statistics->stats_arr[l].data.createHostCopy().data()),
               __LINE__);
      }
    }
  }

  ncwrap(ncmpi_end_indep_data(ncid), __LINE__);
  ncwrap(ncmpi_close(ncid), __LINE__);
}
} // namespace pamc
