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

  //Read config file
  YAML::Node config = YAML::LoadFile(inFile);
  
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
  i_beg = (int) round( nper* px    );
  i_end = (int) round( nper*(px+1) )-1;
  crm_nx = i_end - i_beg + 1;

  nper = ((double) ny_glob)/nprocy;
  j_beg = (int) round( nper* py    );
  j_end = (int) round( nper*(py+1) )-1;
  crm_ny = j_end - j_beg + 1;
  
}

void set_domain_sizes(const YAML::Node &config, real &xlen, real &ylen, real &zlen)
{
#ifdef PAMC_DYCORE
  std::unique_ptr<pamc::TestCase> testcase;
  testcase_from_config(testcase, config);
  const auto [Lx, Ly, Lz] = testcase->get_domain();
  xlen = Lx;
  ylen = Ly;
  zlen = Lz;
#endif
}
