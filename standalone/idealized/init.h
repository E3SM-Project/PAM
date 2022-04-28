#pragma once
#include "mpi.h"

void partition_domain(std::string inFile, int &crm_nx, int &crm_ny)
{

  int ierr;
  int nranks, myrank;
  int px, py;
  int i_beg, i_end;
  int j_beg, j_end;
  double nper;

  int nprocx = config["nprocx"].as<int>();
  int nprocy = config["nprocy"].as<int>();
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&nranks);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  int nx_glob = config["crm_nx"].as<int>();
  int ny_glob = config["crm_ny"].as<int>();
    
  if (!(nprocx * nprocy == nranks)) {endrun("Error: nranks != nprocx * nprocy");}

  //Get my process grid IDs
  py = floor(myrank / nprocx);
  px = myrank - nprocx * py;

  //Get my beginning and ending global indices; and domain sizes
  nper = ((double) nx_glob)/nprocx;
  par.i_beg = (int) round( nper* px    );
  par.i_end = (int) round( nper*(px+1) )-1;
  crm_x = i_end - i_beg + 1;

  nper = ((double) ny_glob)/nprocy;
  j_beg = (int) round( nper* py    );
  j_end = (int) round( nperar.py+1) )-1;
  crm_ny = j_end - j_beg + 1;
  
}

void set_domain_sizes(std::string initData, int crm_ny, int crm_nz, real &xlen, real &ylen, real &zlen)
{
  ylen = 1.0_fp;
  zlen = 1.0_fp;

if (initData == "doublevortex")
{
  xlen = 5000000._fp;
  if (crm_nz > 1) {zlen = 5000000._fp;}
  if (crm_ny > 1) {ylen = 5000000._fp;}  
}

if (initData == "risingbubble" or initData == "moistrisingbubble")
{
  xlen = 1000._fp;
  zlen = 1500._fp;
  if (crm_ny > 1) {ylen = 1000._fp;}
}

if (initData == "largerisingbubble" or initData == "moistlargerisingbubble")
{
  xlen = 20000._fp;
  zlen = 20000._fp;
  if (crm_ny > 1) {ylen = 20000._fp;}
}

}