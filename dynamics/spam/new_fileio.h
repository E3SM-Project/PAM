#pragma once

#include "common.h"
#include "variable_sets.h"
#include "stats.h"
#include "YAKL_pnetcdf.h"

template<uint nprog, uint nconst, uint ndiag, uint nstats> class FileIO {

public:
  
  bool is_initialized;
  int masterproc;
  yakl::SimplePNetCDF nc;
  //MPI_Offset ulIndex = 0; // Unlimited dimension index to place this data at
  
  std::array<real5d, nprog> prog_temp_arr;
  std::array<real5d, nconst> const_temp_arr;
  std::array<real5d, ndiag> diag_temp_arr;
  std::string outputName;

  const VariableSet<nprog> *prog_vars;
  const VariableSet<nconst> *const_vars;
  const VariableSet<ndiag> *diag_vars;
  Stats<nprog, nconst, nstats> *statistics;
  
  FileIO();
  FileIO( const FileIO<nprog,nconst,ndiag,nstats> &fio) = delete;
  FileIO& operator=( const FileIO<nprog,nconst,ndiag,nstats> &fio) = delete;
  void initialize(std::string outputName, Topology &ptopo, Topology &dtopo, Parallel &par, const VariableSet<nprog> &progvars, const VariableSet<nconst> &const_vars, const VariableSet<ndiag> &diagvars, Stats<nprog, nconst, nstats> &stats);
  void output(int nstep, real time);
  void outputInit(real time);
  void outputStats(const Stats<nprog, nconst, nstats> &stats);
};


template<uint nprog, uint nconst, uint ndiag, uint nstats> FileIO<nprog,nconst,ndiag,nstats>::FileIO()
{
  this->is_initialized = false;
  std::cout << "CREATED FILEIO\n";
}

template<uint nprog, uint nconst, uint ndiag, uint nstats> void FileIO<nprog,nconst,ndiag,nstats>::initialize(std::string outName, Topology &ptopo, Topology &dtopo, Parallel &par, const VariableSet<nprog> &progvars, const VariableSet<nconst> &constvars, const VariableSet<ndiag> &diagvars, Stats<nprog, nconst, nstats> &stats)
{
     this->outputName = outName;
     this->prog_vars = &progvars;
     this->const_vars = &constvars;
     this->diag_vars = &diagvars;
     this->statistics = &stats;
     this->masterproc = par.masterproc;
     std::cout << outName << "\n";
     
     nc.create(this->outputName);
     std::cout << outName << "\n";

     // nc.create_unlim_dim( "t" );
     // nc.create_dim( "primal_ncells_x" ,  ptopo.nx_glob );
     // nc.create_dim( "primal_ncells_y" ,  ptopo.ny_glob );
     // nc.create_dim( "primal_nlayers" ,  ptopo.nl );
     // nc.create_dim( "primal_ninterfaces" ,  ptopo.ni );
     // nc.create_dim( "dual_ncells_x" ,  dtopo.nx_glob );
     // nc.create_dim( "dual_ncells_y" ,  dtopo.ny_glob );
     // nc.create_dim( "dual_nlayers" ,  dtopo.nl );
     // nc.create_dim( "dual_ninterfaces" ,  dtopo.ni );
     // nc.create_dim( "nens" ,  ptopo.nens );
     // 
     // nc.create_var<real>( "t" , {"t"} );
     //EVENTUALLY CREATE X,Y,Z COORDINATES AS WELL...
     //BUT THESE WILL BE STAGGERED....!
     
    //  for (int i=0; i<this->const_vars->fields_arr.size(); i++)
    //  {
    //    if (this->const_vars->fields_arr[i].topology->primal) { 
    //      if (this->const_vars->fields_arr[i].extdof == 1)
    //      {nc.create_var<real>( this->const_vars->fields_arr[i].name , {"primal_nlayers","primal_ncells_y","primal_ncells_x", "nens"});}
    //     if (this->const_vars->fields_arr[i].extdof == 0)
    //     {nc.create_var<real>( this->const_vars->fields_arr[i].name , {"primal_ninterfaces","primal_ncells_y","primal_ncells_x", "nens"});}
    //     }
    //    else {
    //     if (this->const_vars->fields_arr[i].extdof == 1)
    //     {nc.create_var<real>( this->const_vars->fields_arr[i].name , {"dual_nlayers","dual_ncells_y","dual_ncells_x", "nens"});}
    //     if (this->const_vars->fields_arr[i].extdof == 0)
    //     {nc.create_var<real>( this->const_vars->fields_arr[i].name , {"dual_ninterfaces","dual_ncells_y","dual_ncells_x", "nens"});}
    //     }
    // this->const_temp_arr[i] = real5d(this->const_vars->fields_arr[i].name.c_str(), this->const_vars->fields_arr[i].total_dofs, this->const_vars->fields_arr[i]._nz, this->const_vars->fields_arr[i].topology->n_cells_y, this->const_vars->fields_arr[i].topology->n_cells_x, this->const_vars->fields_arr[i].topology->nens);
    // }

   //  for (int i=0; i<this->prog_vars->fields_arr.size(); i++)
   //  {
   //    if (this->prog_vars->fields_arr[i].topology->primal) { 
   //      if (this->prog_vars->fields_arr[i].extdof == 1)
   //      {nc.create_var<real>( this->prog_vars->fields_arr[i].name , {"t", "primal_nlayers","primal_ncells_y","primal_ncells_x", "nens"});}
   //     if (this->prog_vars->fields_arr[i].extdof == 0)
   //     {nc.create_var<real>( this->prog_vars->fields_arr[i].name , {"t", "primal_ninterfaces","primal_ncells_y","primal_ncells_x", "nens"});}
   //     }
   //    else {
   //     if (this->prog_vars->fields_arr[i].extdof == 1)
   //     {nc.create_var<real>( this->prog_vars->fields_arr[i].name , {"t", "dual_nlayers","dual_ncells_y","dual_ncells_x", "nens"});}
   //     if (this->prog_vars->fields_arr[i].extdof == 0)
   //     {nc.create_var<real>( this->prog_vars->fields_arr[i].name , {"t", "dual_ninterfaces","dual_ncells_y","dual_ncells_x", "nens"});}
   //     }
   // this->prog_temp_arr[i] = real5d(this->prog_vars->fields_arr[i].name.c_str(), this->prog_vars->fields_arr[i].total_dofs, this->prog_vars->fields_arr[i]._nz, this->prog_vars->fields_arr[i].topology->n_cells_y, this->prog_vars->fields_arr[i].topology->n_cells_x, this->prog_vars->fields_arr[i].topology->nens);
   // }

  //  for (int i=0; i<this->diag_vars->fields_arr.size(); i++)
  //  {
  //    if (this->diag_vars->fields_arr[i].topology->primal) { 
  //      if (this->diag_vars->fields_arr[i].extdof == 1)
  //      {nc.create_var<real>( this->diag_vars->fields_arr[i].name , {"t", "primal_nlayers","primal_ncells_y","primal_ncells_x", "nens"});}
  //     if (this->diag_vars->fields_arr[i].extdof == 0)
  //     {nc.create_var<real>( this->diag_vars->fields_arr[i].name , {"t", "primal_ninterfaces","primal_ncells_y","primal_ncells_x", "nens"});}
  //     }
  //    else {
  //     if (this->diag_vars->fields_arr[i].extdof == 1)
  //     {nc.create_var<real>( this->diag_vars->fields_arr[i].name , {"t", "dual_nlayers","dual_ncells_y","dual_ncells_x", "nens"});}
  //     if (this->diag_vars->fields_arr[i].extdof == 0)
  //     {nc.create_var<real>( this->diag_vars->fields_arr[i].name , {"t", "dual_ninterfaces","dual_ncells_y","dual_ncells_x", "nens"});}
  //     }
  // this->diag_temp_arr[i] = real5d(this->diag_vars->fields_arr[i].name.c_str(), this->diag_vars->fields_arr[i].total_dofs, this->diag_vars->fields_arr[i]._nz, this->diag_vars->fields_arr[i].topology->n_cells_y, this->diag_vars->fields_arr[i].topology->n_cells_x, this->diag_vars->fields_arr[i].topology->nens);
  // }

//   nc.create_dim( "statsize" ,  this->statistics->statsize );
//   for (int l=0; l<this->statistics->stats_arr.size(); l++)
//   {
//     nc.create_dim( this->statistics->stats_arr[l].name + "_ndofs" ,  this->statistics->stats_arr[l].ndofs );
//     nc.create_var<real>( this->statistics->stats_arr[l].name , {this->statistics->stats_arr[l].name + "_ndofs", "statsize", "nens"});
// }
      //nc.enddef();
      nc.close();
  this->is_initialized = true;

   }