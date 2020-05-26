
#ifndef _FILEIO_H_
#define _FILEIO_H_




#include "common.h"
#include <iostream>
#include "variable_sets.h"
#include "pnetcdf.h"
#include "mpi.h"
#include <cstring>
#include "model.h"

//Error reporting routine for the PNetCDF I/O
void ncwrap( int ierr , int line ) {
  if (ierr != NC_NOERR) {
    std::cout << "NetCDF Error at line: " << line << "\n";
    std::cout << ncmpi_strerror(ierr) << "\n";
    exit(-1);
  }
  }


template<uint nprog, uint nconst, uint ndiag, uint nstats> class FileIO {

public:
  bool is_initialized;
  int masterproc;
  int ncid;
  int const_dim_ids[4], prog_dim_ids[5], diag_dim_ids[5];
  int tDim, xDim, yDim, zDim;
  int const_var_dim_ids[nconst], prog_var_dim_ids[nprog], diag_var_dim_ids[ndiag];
  int const_var_ids[nconst], prog_var_ids[nprog], diag_var_ids[ndiag];
  int tVar;
  std::string outputName;
  std::array<realArr, nprog> prog_temp_arr;
  std::array<realArr, nconst> const_temp_arr;
  std::array<realArr, ndiag> diag_temp_arr;
  const VariableSet<nprog> *prog_vars;
  const VariableSet<nconst> *const_vars;
  const VariableSet<ndiag> *diag_vars;
  Stats<nprog, nconst, nstats> *statistics;

  int statDim;
  int stat_ids[nstats];

  int numOut;

  MPI_Offset const_start[4], const_count[4];
  MPI_Offset prog_start[5], prog_count[5];
  MPI_Offset diag_start[5], diag_count[5];

  FileIO();
  FileIO( const FileIO<nprog,nconst,ndiag,nstats> &fio) = delete;
  FileIO& operator=( const FileIO<nprog,nconst,ndiag,nstats> &fio) = delete;
  void initialize(std::string outputName, Topology &topo, Parallel &par, const VariableSet<nprog> &progvars, const VariableSet<nconst> &const_vars, const VariableSet<ndiag> &diagvars, Stats<nprog, nconst, nstats> &stats);
  void output(int nstep, real time);
  void outputInit(real time);
  void outputStats(const Stats<nprog, nconst, nstats> &stats);
  void close();


};




    template<uint nprog, uint nconst, uint ndiag, uint nstats> FileIO<nprog,nconst,ndiag,nstats>::FileIO()
    {
      this->is_initialized = false;
      std::cout << "CREATED FILEIO\n";
    }


  template<uint nprog, uint nconst, uint ndiag, uint nstats> void FileIO<nprog,nconst,ndiag,nstats>::initialize(std::string outName, Topology &topo, Parallel &par, const VariableSet<nprog> &progvars, const VariableSet<nconst> &constvars, const VariableSet<ndiag> &diagvars, Stats<nprog, nconst, nstats> &stats)
  {
       this->outputName = outName;
       this->prog_vars = &progvars;
       this->const_vars = &constvars;
       this->diag_vars = &diagvars;
       this->statistics = &stats;
       this->masterproc = par.masterproc;

       ncwrap( ncmpi_create( MPI_COMM_WORLD , this->outputName.c_str() , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

       //Define time and cell dimensions
       ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
       ncwrap( ncmpi_def_dim( ncid , "ncells_x" , (MPI_Offset) topo.nx_glob  , &xDim ) , __LINE__ );
       ncwrap( ncmpi_def_dim( ncid , "ncells_y" , (MPI_Offset) topo.ny_glob  , &yDim ) , __LINE__ );
       ncwrap( ncmpi_def_dim( ncid , "ncells_z" , (MPI_Offset) topo.nz_glob  , &zDim ) , __LINE__ );

       //Create time dimension
       const_dim_ids[0] = tDim;
       ncwrap( ncmpi_def_var( ncid , "t"      , REAL_NC , 1 , const_dim_ids , &tVar ) , __LINE__ );
       ncwrap( ncmpi_put_att_text (ncid, tVar, "units", std::strlen("seconds"), "seconds"), __LINE__ );


       for (int i=0; i<this->const_vars->fields_arr.size(); i++)
       {
         ncwrap( ncmpi_def_dim( ncid , (this->const_vars->fields_arr[i].name + "_ndofs").c_str() , (MPI_Offset) this->const_vars->fields_arr[i].total_dofs , &const_var_dim_ids[i] ) , __LINE__ );
         const_dim_ids[0] = const_var_dim_ids[i]; const_dim_ids[1] = zDim; const_dim_ids[2] = yDim; const_dim_ids[3] = xDim;
         ncwrap( ncmpi_def_var( ncid , this->const_vars->fields_arr[i].name.c_str() , REAL_NC , 4 , const_dim_ids , &const_var_ids[i]  ) , __LINE__ );
         this->const_temp_arr[i] = realArr(this->const_vars->fields_arr[i].name.c_str(), this->const_vars->fields_arr[i].total_dofs, this->const_vars->fields_arr[i].topology->n_cells_z, this->const_vars->fields_arr[i].topology->n_cells_y, this->const_vars->fields_arr[i].topology->n_cells_x);
       }

       for (int i=0; i<this->prog_vars->fields_arr.size(); i++)
       {
         ncwrap( ncmpi_def_dim( ncid , (this->prog_vars->fields_arr[i].name + "_ndofs").c_str() , (MPI_Offset) this->prog_vars->fields_arr[i].total_dofs , &prog_var_dim_ids[i] ) , __LINE__ );
         prog_dim_ids[0] = tVar; prog_dim_ids[1] = prog_var_dim_ids[i]; prog_dim_ids[2] = zDim; prog_dim_ids[3] = yDim; prog_dim_ids[4] = xDim;
         ncwrap( ncmpi_def_var( ncid , this->prog_vars->fields_arr[i].name.c_str() , REAL_NC , 5 , prog_dim_ids , &prog_var_ids[i]  ) , __LINE__ );
         this->prog_temp_arr[i] = realArr(this->prog_vars->fields_arr[i].name.c_str(), this->prog_vars->fields_arr[i].total_dofs, this->prog_vars->fields_arr[i].topology->n_cells_z, this->prog_vars->fields_arr[i].topology->n_cells_y, this->prog_vars->fields_arr[i].topology->n_cells_x);
       }

       for (int i=0; i<this->diag_vars->fields_arr.size(); i++)
       {
         ncwrap( ncmpi_def_dim( ncid , (this->diag_vars->fields_arr[i].name + "_ndofs").c_str() , (MPI_Offset) this->diag_vars->fields_arr[i].total_dofs , &diag_var_dim_ids[i] ) , __LINE__ );
         diag_dim_ids[0] = tVar; diag_dim_ids[1] = diag_var_dim_ids[i]; diag_dim_ids[2] = zDim; diag_dim_ids[3] = yDim; diag_dim_ids[4] = xDim;
         ncwrap( ncmpi_def_var( ncid , this->diag_vars->fields_arr[i].name.c_str() , REAL_NC , 5 , diag_dim_ids , &diag_var_ids[i]  ) , __LINE__ );
         this->diag_temp_arr[i] = realArr(this->diag_vars->fields_arr[i].name.c_str(), this->diag_vars->fields_arr[i].total_dofs, this->diag_vars->fields_arr[i].topology->n_cells_z, this->diag_vars->fields_arr[i].topology->n_cells_y, this->diag_vars->fields_arr[i].topology->n_cells_x);
       }

       // define statistics dimension and variables

       ncwrap( ncmpi_def_dim( ncid , "nstat" , (MPI_Offset) this->statistics->statsize , &statDim ) , __LINE__ );
       for (int l=0; l<this->statistics->stats_arr.size(); l++)
       {
       ncwrap( ncmpi_def_var( ncid , this->statistics->stats_arr[l].name.c_str() , REAL_NC , 1 , &statDim , &stat_ids[l]  ) , __LINE__ );
      }

       ncwrap( ncmpi_enddef( ncid ) , __LINE__ );


       ncwrap( ncmpi_close(ncid) , __LINE__ );
       this->is_initialized = true;
    }



  template<uint nprog, uint nconst, uint ndiag, uint nstats> void FileIO<nprog,nconst,ndiag,nstats>::output(int nstep, real time)
  {

    ncwrap( ncmpi_open( MPI_COMM_WORLD , this->outputName.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );

      for (int l=0; l<this->prog_vars->fields_arr.size(); l++)
      {
        ncwrap( ncmpi_inq_varid( ncid , this->prog_vars->fields_arr[l].name.c_str() , &prog_var_ids[l]  ) , __LINE__ );
        int is = this->prog_vars->fields_arr[l].topology->is;
        int js = this->prog_vars->fields_arr[l].topology->js;
        int ks = this->prog_vars->fields_arr[l].topology->ks;
        yakl::parallel_for("CopyFieldToOutputBuffer", this->prog_vars->fields_arr[l].topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, this->prog_vars->fields_arr[l].topology->n_cells_z, this->prog_vars->fields_arr[l].topology->n_cells_y, this->prog_vars->fields_arr[l].topology->n_cells_x, k, j, i);
          for (int ndof=0; ndof<this->prog_vars->fields_arr[l].total_dofs; ndof++) {
            this->prog_temp_arr[l](ndof, k, j, i) = this->prog_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
          }
        });

        prog_start[0] = this->numOut; prog_start[1] = 0; prog_start[2] = this->prog_vars->fields_arr[l].topology->k_beg; prog_start[3] = this->prog_vars->fields_arr[l].topology->j_beg; prog_start[4] = this->prog_vars->fields_arr[l].topology->i_beg;
        prog_count[0] = 1; prog_count[1] = this->prog_vars->fields_arr[l].total_dofs; prog_count[2] = this->prog_vars->fields_arr[l].topology->n_cells_z; prog_count[3] = this->prog_vars->fields_arr[l].topology->n_cells_y; prog_count[4] = this->prog_vars->fields_arr[l].topology->n_cells_x;
        ncwrap( PNETCDF_PUT_VAR_ALL( ncid , prog_var_ids[l] , prog_start , prog_count , this->prog_temp_arr[l].createHostCopy().data() ) , __LINE__ );
      }


      for (int l=0; l<this->diag_vars->fields_arr.size(); l++)
      {
        ncwrap( ncmpi_inq_varid( ncid , this->diag_vars->fields_arr[l].name.c_str() , &diag_var_ids[l]  ) , __LINE__ );
        int is = this->diag_vars->fields_arr[l].topology->is;
        int js = this->diag_vars->fields_arr[l].topology->js;
        int ks = this->diag_vars->fields_arr[l].topology->ks;
        yakl::parallel_for("CopyFieldToOutputBuffer", this->diag_vars->fields_arr[l].topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          yakl::unpackIndices(iGlob, this->diag_vars->fields_arr[l].topology->n_cells_z, this->diag_vars->fields_arr[l].topology->n_cells_y, this->diag_vars->fields_arr[l].topology->n_cells_x, k, j, i);
          for (int ndof=0; ndof<this->diag_vars->fields_arr[l].total_dofs; ndof++) {
            this->diag_temp_arr[l](ndof, k, j, i) = this->diag_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
          }
        });

        diag_start[0] = this->numOut; diag_start[1] = 0; diag_start[2] = this->diag_vars->fields_arr[l].topology->k_beg; diag_start[3] = this->diag_vars->fields_arr[l].topology->j_beg; diag_start[4] = this->diag_vars->fields_arr[l].topology->i_beg;
        diag_count[0] = 1; diag_count[1] = this->diag_vars->fields_arr[l].total_dofs; diag_count[2] = this->diag_vars->fields_arr[l].topology->n_cells_z; diag_count[3] = this->diag_vars->fields_arr[l].topology->n_cells_y; diag_count[4] = this->diag_vars->fields_arr[l].topology->n_cells_x;
        ncwrap( PNETCDF_PUT_VAR_ALL( ncid , diag_var_ids[l] , diag_start , diag_count , this->diag_temp_arr[l].createHostCopy().data() ) , __LINE__ );
      }

      // ADD THIS
            // write elapsed time/time step

      ncwrap( ncmpi_close(ncid) , __LINE__ );
      this->numOut++;
    }

  template<uint nprog, uint nconst, uint ndiag, uint nstats> void FileIO<nprog,nconst,ndiag,nstats>::outputInit(real time)
  {
    ncwrap( ncmpi_open( MPI_COMM_WORLD , this->outputName.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );

    for (int l=0; l<this->const_vars->fields_arr.size(); l++)
    {
      ncwrap( ncmpi_inq_varid( ncid , this->const_vars->fields_arr[l].name.c_str() , &const_var_ids[l]  ) , __LINE__ );
      int is = this->const_vars->fields_arr[l].topology->is;
      int js = this->const_vars->fields_arr[l].topology->js;
      int ks = this->const_vars->fields_arr[l].topology->ks;
      yakl::parallel_for("CopyFieldToOutputBuffer", this->const_vars->fields_arr[l].topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, this->const_vars->fields_arr[l].topology->n_cells_z, this->const_vars->fields_arr[l].topology->n_cells_y, this->const_vars->fields_arr[l].topology->n_cells_x, k, j, i);
        for (int ndof=0; ndof<this->const_vars->fields_arr[l].total_dofs; ndof++) {
          this->const_temp_arr[l](ndof, k, j, i) = this->const_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
        }
      });

      const_start[0] = 0; const_start[1] = this->const_vars->fields_arr[l].topology->k_beg; const_start[2] = this->const_vars->fields_arr[l].topology->j_beg; const_start[3] = this->const_vars->fields_arr[l].topology->i_beg;
      const_count[0] = this->const_vars->fields_arr[l].total_dofs; const_count[1] = this->const_vars->fields_arr[l].topology->n_cells_z; const_count[2] = this->const_vars->fields_arr[l].topology->n_cells_y; const_count[3] = this->const_vars->fields_arr[l].topology->n_cells_x;
      ncwrap( PNETCDF_PUT_VAR_ALL( ncid , const_var_ids[l] , const_start , const_count , this->const_temp_arr[l].createHostCopy().data() ) , __LINE__ );
    }

    for (int l=0; l<this->prog_vars->fields_arr.size(); l++)
    {
      ncwrap( ncmpi_inq_varid( ncid , this->prog_vars->fields_arr[l].name.c_str() , &prog_var_ids[l]  ) , __LINE__ );
      int is = this->prog_vars->fields_arr[l].topology->is;
      int js = this->prog_vars->fields_arr[l].topology->js;
      int ks = this->prog_vars->fields_arr[l].topology->ks;
      yakl::parallel_for("CopyFieldToOutputBuffer", this->prog_vars->fields_arr[l].topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, this->prog_vars->fields_arr[l].topology->n_cells_z, this->prog_vars->fields_arr[l].topology->n_cells_y, this->prog_vars->fields_arr[l].topology->n_cells_x, k, j, i);
        for (int ndof=0; ndof<this->prog_vars->fields_arr[l].total_dofs; ndof++) {
          this->prog_temp_arr[l](ndof, k, j, i) = this->prog_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
        }
      });

      prog_start[0] = this->numOut; prog_start[1] = 0; prog_start[2] = this->prog_vars->fields_arr[l].topology->k_beg; prog_start[3] = this->prog_vars->fields_arr[l].topology->j_beg; prog_start[4] = this->prog_vars->fields_arr[l].topology->i_beg;
      prog_count[0] = 1; prog_count[1] = this->prog_vars->fields_arr[l].total_dofs; prog_count[2] = this->prog_vars->fields_arr[l].topology->n_cells_z; prog_count[3] = this->prog_vars->fields_arr[l].topology->n_cells_y; prog_count[4] = this->prog_vars->fields_arr[l].topology->n_cells_x;
      ncwrap( PNETCDF_PUT_VAR_ALL( ncid , prog_var_ids[l] , prog_start , prog_count , this->prog_temp_arr[l].createHostCopy().data() ) , __LINE__ );
    }

    for (int l=0; l<this->diag_vars->fields_arr.size(); l++)
    {
      ncwrap( ncmpi_inq_varid( ncid , this->diag_vars->fields_arr[l].name.c_str() , &diag_var_ids[l]  ) , __LINE__ );
      int is = this->diag_vars->fields_arr[l].topology->is;
      int js = this->diag_vars->fields_arr[l].topology->js;
      int ks = this->diag_vars->fields_arr[l].topology->ks;
      yakl::parallel_for("CopyFieldToOutputBuffer", this->diag_vars->fields_arr[l].topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, this->diag_vars->fields_arr[l].topology->n_cells_z, this->diag_vars->fields_arr[l].topology->n_cells_y, this->diag_vars->fields_arr[l].topology->n_cells_x, k, j, i);
        for (int ndof=0; ndof<this->diag_vars->fields_arr[l].total_dofs; ndof++) {
          this->diag_temp_arr[l](ndof, k, j, i) = this->diag_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
        }
      });

      diag_start[0] = this->numOut; diag_start[1] = 0; diag_start[2] = this->diag_vars->fields_arr[l].topology->k_beg; diag_start[3] = this->diag_vars->fields_arr[l].topology->j_beg; diag_start[4] = this->diag_vars->fields_arr[l].topology->i_beg;
      diag_count[0] = 1; diag_count[1] = this->diag_vars->fields_arr[l].total_dofs; diag_count[2] = this->diag_vars->fields_arr[l].topology->n_cells_z; diag_count[3] = this->diag_vars->fields_arr[l].topology->n_cells_y; diag_count[4] = this->diag_vars->fields_arr[l].topology->n_cells_x;
      ncwrap( PNETCDF_PUT_VAR_ALL( ncid , diag_var_ids[l] , diag_start , diag_count , this->diag_temp_arr[l].createHostCopy().data() ) , __LINE__ );
    }

// ADD THIS
      // write elapsed time/time step

      ncwrap( ncmpi_close(ncid) , __LINE__ );

      this->numOut++;

    }

  template<uint nprog, uint nconst, uint ndiag, uint nstats> void FileIO<nprog,nconst,ndiag,nstats>::close() { }

  template<uint nprog, uint nconst, uint ndiag, uint nstats>  void FileIO<nprog,nconst,ndiag,nstats>::outputStats(const Stats<nprog, nconst, nstats> &stats)
  {
    ncwrap( ncmpi_open( MPI_COMM_WORLD , this->outputName.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );
    ncwrap( ncmpi_begin_indep_data(ncid) , __LINE__ );

    if (masterproc)
    {
    MPI_Offset statStart[1], statCount[1];
    statStart[0] = 0; statCount[0] = this->statistics->statsize;
    for (int l=0; l<this->statistics->stats_arr.size(); l++)
    {
        ncwrap( ncmpi_inq_varid( ncid , this->statistics->stats_arr[l].name.c_str() , &stat_ids[l]  ) , __LINE__ );
      ncwrap( PNETCDF_PUT_VAR(ncid, stat_ids[l], statStart, statCount, this->statistics->stats_arr[l].data.createHostCopy().data()) , __LINE__ );
    }
    }

    ncwrap( ncmpi_end_indep_data(ncid) , __LINE__ );
    ncwrap( ncmpi_close(ncid) , __LINE__ );

  }




#endif
