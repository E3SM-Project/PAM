
#ifndef _FILEIO_H_
#define _FILEIO_H_




#include "common.h"
#include <iostream>
#include <fstream>
#include "variable_sets.h"
#include "pnetcdf.h"
#include "mpi.h"
#include <cstring>

//Error reporting routine for the PNetCDF I/O
void ncwrap( int ierr , int line ) {
  if (ierr != NC_NOERR) {
    std::cout << "NetCDF Error at line: " << line << "\n";
    std::cout << ncmpi_strerror(ierr) << "\n";
    exit(-1);
  }
  }


template<uint ndims, uint nprog, uint nconst, uint ndiag> class FileIO {
  std::ofstream file;

public:
  bool is_initialized;
  int ncid;
  int const_dim_ids[4], prog_dim_ids[5];
  int tDim, xDim, yDim, zDim;
  int const_var_dim_ids[nconst], prog_var_dim_ids[nprog];
  int const_var_ids[nconst], prog_var_ids[nprog];
  int tVar;
  std::string outputName;
  std::array<realArr, nprog> prog_temp_arr;
  std::array<realArr, nconst> const_temp_arr;
  const VariableSet<ndims, nprog> *prog_vars;
  const VariableSet<ndims, nconst> *const_vars;

  int numOut;

  MPI_Offset const_start[4], const_count[4];
  MPI_Offset prog_start[5], prog_count[5];

  FileIO();
  FileIO( const FileIO<ndims,nprog,nconst,ndiag> &fio) = delete;
  FileIO& operator=( const FileIO<ndims,nprog,nconst,ndiag> &fio) = delete;
  void initialize(std::string outputName, Topology<ndims> &topo, const VariableSet<ndims, nprog> &progvars, const VariableSet<ndims, nconst> &const_vars);
  void output(int nstep, real time);
  void outputInit(real time);
  void close();


};


    template<uint ndims, uint nprog, uint nconst, uint ndiag> FileIO<ndims,nprog,nconst,ndiag>::FileIO()
    {
      this->is_initialized = false;
      std::cout << "CREATED FILEIO\n";
    }


  template<uint ndims, uint nprog, uint nconst, uint ndiag> void FileIO<ndims,nprog,nconst,ndiag>::initialize(std::string outName, Topology<ndims> &topo, const VariableSet<ndims, nprog> &progvars, const VariableSet<ndims, nconst> &constvars)
  {
       this->file.open ("output.txt");
       this->outputName = outName;
       this->prog_vars = &progvars;
       this->const_vars = &constvars;

       ncwrap( ncmpi_create( MPI_COMM_WORLD , this->outputName.c_str() , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

       //Define time and cell dimensions
       ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
       ncwrap( ncmpi_def_dim( ncid , "ncells_x" , (MPI_Offset) topo.n_cells_x  , &xDim ) , __LINE__ );
       ncwrap( ncmpi_def_dim( ncid , "ncells_y" , (MPI_Offset) topo.n_cells_y  , &yDim ) , __LINE__ );
       ncwrap( ncmpi_def_dim( ncid , "ncells_z" , (MPI_Offset) topo.n_cells_z  , &zDim ) , __LINE__ );

       //Create time dimension
       const_dim_ids[0] = tDim;
       ncwrap( ncmpi_def_var( ncid , "t"      , NC_FLOAT , 1 , const_dim_ids , &tVar ) , __LINE__ );
       ncwrap( ncmpi_put_att_text (ncid, tVar, "units", std::strlen("seconds"), "seconds"), __LINE__ );


       for (int i=0; i<this->const_vars->fields_arr.size(); i++)
       {
         ncwrap( ncmpi_def_dim( ncid , (this->const_vars->fields_arr[i].name + "_ndofs").c_str() , (MPI_Offset) this->const_vars->fields_arr[i].total_dofs , &const_var_dim_ids[i] ) , __LINE__ );
         const_dim_ids[0] = const_var_dim_ids[i]; const_dim_ids[1] = zDim; const_dim_ids[2] = yDim; const_dim_ids[3] = xDim;
         ncwrap( ncmpi_def_var( ncid , this->const_vars->fields_arr[i].name.c_str() , NC_FLOAT , 4 , const_dim_ids , &const_var_ids[i]  ) , __LINE__ );
         this->const_temp_arr[i] = realArr(this->const_vars->fields_arr[i].name.c_str(), this->const_vars->fields_arr[i].total_dofs, this->const_vars->fields_arr[i].topology->n_cells_z, this->const_vars->fields_arr[i].topology->n_cells_y, this->const_vars->fields_arr[i].topology->n_cells_x);
       }

       for (int i=0; i<this->prog_vars->fields_arr.size(); i++)
       {
         ncwrap( ncmpi_def_dim( ncid , (this->prog_vars->fields_arr[i].name + "_ndofs").c_str() , (MPI_Offset) this->prog_vars->fields_arr[i].total_dofs , &prog_var_dim_ids[i] ) , __LINE__ );
         prog_dim_ids[0] = tVar; prog_dim_ids[1] = prog_var_dim_ids[i]; prog_dim_ids[2] = zDim; prog_dim_ids[3] = yDim; prog_dim_ids[4] = xDim;
         ncwrap( ncmpi_def_var( ncid , this->prog_vars->fields_arr[i].name.c_str() , NC_FLOAT , 5 , prog_dim_ids , &prog_var_ids[i]  ) , __LINE__ );
         this->prog_temp_arr[i] = realArr(this->prog_vars->fields_arr[i].name.c_str(), this->prog_vars->fields_arr[i].total_dofs, this->prog_vars->fields_arr[i].topology->n_cells_z, this->prog_vars->fields_arr[i].topology->n_cells_y, this->prog_vars->fields_arr[i].topology->n_cells_x);
       }


       ncwrap( ncmpi_enddef( ncid ) , __LINE__ );


       ncwrap( ncmpi_close(ncid) , __LINE__ );
       this->is_initialized = true;
    }

// have an internal write field routine
// that writes a field out
// how big should temporary arrays be? this is the question...

  template<uint ndims, uint nprog, uint nconst, uint ndiag> void FileIO<ndims,nprog,nconst,ndiag>::output(int nstep, real time)
  {
    ncwrap( ncmpi_open( MPI_COMM_WORLD , this->outputName.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );

      // open netcdf file
      // open variables for constants and vars
      // write variables for constants and vars
      // write elapsed time/time step

      for (int l=0; l<this->prog_vars->fields_arr.size(); l++)
      {
        ncwrap( ncmpi_inq_varid( ncid , this->prog_vars->fields_arr[l].name.c_str() , &prog_var_ids[l]  ) , __LINE__ );
        int is = this->prog_vars->fields_arr[l].topology->is;
        int js = this->prog_vars->fields_arr[l].topology->js;
        int ks = this->prog_vars->fields_arr[l].topology->ks;
        yakl::parallel_for("CopyFieldToOutputBuffer", this->prog_vars->fields_arr[l].topology->n_cells, YAKL_LAMBDA (int iGlob) {
          int k, j, i;
          //std::cout <<"copy field to output buffer " << l << " " << this->const_vars->fields_arr[l].name << "\n" << std::flush;
          yakl::unpackIndices(iGlob, this->prog_vars->fields_arr[l].topology->n_cells_z, this->prog_vars->fields_arr[l].topology->n_cells_y, this->prog_vars->fields_arr[l].topology->n_cells_x, k, j, i);
          for (int ndof=0; ndof<this->prog_vars->fields_arr[l].total_dofs; ndof++) {
            //std::cout << i << " " << j << " " << k << " " << ndof << "\n" << std::flush;
            this->prog_temp_arr[l](ndof, k, j, i) = this->prog_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
          }
        });

        prog_start[0] = this->numOut; prog_start[1] = 0; prog_start[2] = 0; prog_start[3] = 0; prog_start[4] = 0;
        prog_count[0] = 1; prog_count[1] = this->prog_vars->fields_arr[l].total_dofs; prog_count[2] = this->prog_vars->fields_arr[l].topology->n_cells_z; prog_count[3] = this->prog_vars->fields_arr[l].topology->n_cells_y; prog_count[4] = this->prog_vars->fields_arr[l].topology->n_cells_x;
        ncwrap( ncmpi_put_vara_float_all( ncid , prog_var_ids[l] , prog_start , prog_count , this->prog_temp_arr[l].createHostCopy().data() ) , __LINE__ );
      }


      this->file << "Variables at step " << nstep << " and t=" << time << "\n";
      for (int i=0; i<this->prog_vars->fields_arr.size(); i++)
      {
        this->file << this->prog_vars->fields_arr[i].name << "\n" << std::flush;
      }

      // ADD THIS
            // write elapsed time/time step

      ncwrap( ncmpi_close(ncid) , __LINE__ );
      this->numOut++;
    }

  template<uint ndims, uint nprog, uint nconst, uint ndiag> void FileIO<ndims,nprog,nconst,ndiag>::outputInit(real time)
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
        //std::cout <<"copy field to output buffer " << l << " " << this->const_vars->fields_arr[l].name << "\n" << std::flush;
        yakl::unpackIndices(iGlob, this->const_vars->fields_arr[l].topology->n_cells_z, this->const_vars->fields_arr[l].topology->n_cells_y, this->const_vars->fields_arr[l].topology->n_cells_x, k, j, i);
        for (int ndof=0; ndof<this->const_vars->fields_arr[l].total_dofs; ndof++) {
          //std::cout << i << " " << j << " " << k << " " << ndof << "\n" << std::flush;
          this->const_temp_arr[l](ndof, k, j, i) = this->const_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
        }
      });

      const_start[0] = 0; const_start[1] = 0; const_start[2] = 0; const_start[3] = 0;
      const_count[0] = this->const_vars->fields_arr[l].total_dofs; const_count[1] = this->const_vars->fields_arr[l].topology->n_cells_z; const_count[2] = this->const_vars->fields_arr[l].topology->n_cells_y; const_count[3] = this->const_vars->fields_arr[l].topology->n_cells_x;
      ncwrap( ncmpi_put_vara_float_all( ncid , const_var_ids[l] , const_start , const_count , this->const_temp_arr[l].createHostCopy().data() ) , __LINE__ );
    }

    for (int l=0; l<this->prog_vars->fields_arr.size(); l++)
    {
      ncwrap( ncmpi_inq_varid( ncid , this->prog_vars->fields_arr[l].name.c_str() , &prog_var_ids[l]  ) , __LINE__ );
      int is = this->prog_vars->fields_arr[l].topology->is;
      int js = this->prog_vars->fields_arr[l].topology->js;
      int ks = this->prog_vars->fields_arr[l].topology->ks;
      yakl::parallel_for("CopyFieldToOutputBuffer", this->prog_vars->fields_arr[l].topology->n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        //std::cout <<"copy field to output buffer " << l << " " << this->const_vars->fields_arr[l].name << "\n" << std::flush;
        yakl::unpackIndices(iGlob, this->prog_vars->fields_arr[l].topology->n_cells_z, this->prog_vars->fields_arr[l].topology->n_cells_y, this->prog_vars->fields_arr[l].topology->n_cells_x, k, j, i);
        for (int ndof=0; ndof<this->prog_vars->fields_arr[l].total_dofs; ndof++) {
          //std::cout << i << " " << j << " " << k << " " << ndof << "\n" << std::flush;
          this->prog_temp_arr[l](ndof, k, j, i) = this->prog_vars->fields_arr[l].data(ndof, k+ks, j+js, i+is);
        }
      });

      prog_start[0] = this->numOut; prog_start[1] = 0; prog_start[2] = 0; prog_start[3] = 0; prog_start[4] = 0;
      prog_count[0] = 1; prog_count[1] = this->prog_vars->fields_arr[l].total_dofs; prog_count[2] = this->prog_vars->fields_arr[l].topology->n_cells_z; prog_count[3] = this->prog_vars->fields_arr[l].topology->n_cells_y; prog_count[4] = this->prog_vars->fields_arr[l].topology->n_cells_x;
      ncwrap( ncmpi_put_vara_float_all( ncid , prog_var_ids[l] , prog_start , prog_count , this->prog_temp_arr[l].createHostCopy().data() ) , __LINE__ );
    }





      this->file << "Constants at step " << 0 << " and t=" << time << "\n";
      for (int i=0; i<this->const_vars->fields_arr.size(); i++)
      {
        this->file << this->const_vars->fields_arr[i].name << "\n" << std::flush;
        //  ncwrap( ncmpi_put_vara_float_all( ncid , uVar , start , count , data.createHostCopy().data() ) , __LINE__ );

      }

      file << "Variables at step " << 0 << " and t=" << time << "\n";
      for (int i=0; i<this->prog_vars->fields_arr.size(); i++)
      {
        this->file << this->prog_vars->fields_arr[i].name << "\n" << std::flush;
      }

// ADD THIS
      // write elapsed time/time step
      ncwrap( ncmpi_close(ncid) , __LINE__ );

      this->numOut++;

    }

  template<uint ndims, uint nprog, uint nconst, uint ndiag> void FileIO<ndims,nprog,nconst,ndiag>::close()
  {
      this->file.close();
    }





#endif
