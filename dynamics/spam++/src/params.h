#ifndef _PARAMS_H_
#define _PARAMS_H_

#include <fstream>
#include <sstream>
#include "string.h"
#include "mpi.h"
#include "parallel.h"
#include "common.h"









template<uint ndims> void readParamsFile(std::string inFile, Parameters &params, Parallel &par) {

  int ierr;

  // Determine if I'm the master process
  if (par.myrank == 0) { par.masterproc = 1;}
  else { par.masterproc = 0; }

  // Read in equals-separated key = value file line by line
  std::ifstream fInStream(inFile);
  std::string line;
  while (std::getline(fInStream, line)) {
    // Remove spaces and tabs from the line
    line.erase (std::remove(line.begin(), line.end(), ' '), line.end());
    line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());

    // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
    if (!line.empty() && line.find("//",0) != 0) {
      // Find the colon
      uint splitloc = line.find('=',0);
      // Store the key and value strings
      std::string key   = line.substr(0,splitloc);
      std::string value = line.substr(splitloc+1,line.length()-splitloc);

      // Transform the value into a string stream for convenience
      std::stringstream ssVal(value);

      // Match the key, and store the value
      if      ( !strcmp( "nx"         , key.c_str() ) ) { ssVal >> params.nx_glob    ; }
      else if ( !strcmp( "ny"         , key.c_str() ) ) { ssVal >> params.ny_glob    ; }
      else if ( !strcmp( "nz"         , key.c_str() ) ) { ssVal >> params.nz_glob    ; }

      else if ( !strcmp( "dt"         , key.c_str() ) ) { ssVal >> params.dt         ; }
      else if ( !strcmp( "Nsteps"     , key.c_str() ) ) { ssVal >> params.Nsteps     ; }
      else if ( !strcmp( "Nout"       , key.c_str() ) ) { ssVal >> params.Nout       ; }
      else if ( !strcmp( "Nstat"      , key.c_str() ) ) { ssVal >> params.Nstat      ; }
      else if ( !strcmp( "TStype"     , key.c_str() ) ) { ssVal >> params.TStype     ; }
      else if ( !strcmp( "cfl"        , key.c_str() ) ) { ssVal >> params.cfl        ; }
      else if ( !strcmp( "outputName" , key.c_str() ) ) { ssVal >> params.outputName ; }

      else if ( !strcmp( "nprocx"     , key.c_str() ) ) { ssVal >> par.nprocx        ; }
      else if ( !strcmp( "nprocy"     , key.c_str() ) ) { ssVal >> par.nprocy        ; }
      else if ( !strcmp( "nprocz"     , key.c_str() ) ) { ssVal >> par.nprocz        ; }



      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

// ADD TSType fixes here

  // Test to make sure all required values were initialized
  if (params.nx_glob == -1) { std::cout << "Error: key " << "nx"          << " not set.\n"; exit(-1); }
  if (params.ny_glob == -1) { std::cout << "Error: key " << "ny"          << " not set.\n"; exit(-1); }
  if (params.nz_glob == -1) { std::cout << "Error: key " << "nz"          << " not set.\n"; exit(-1); }
  if (params.dt      == -1) { std::cout << "Error: key " << "dt"          << " not set.\n"; exit(-1); }
  if (params.Nsteps  == -1) { std::cout << "Error: key " << "Nsteps"      << " not set.\n"; exit(-1); }
  if (params.Nout    == -1) { std::cout << "Error: key " << "Nout"        << " not set.\n"; exit(-1); }
  if (par.nprocx     == -1) { std::cout << "Error: key " << "nprocx"      << " not set.\n"; exit(-1); }
  if (par.nprocy     == -1) { std::cout << "Error: key " << "nprocy"      << " not set.\n"; exit(-1); }
  if (par.nprocz     == -1) { std::cout << "Error: key " << "nprocz"      << " not set.\n"; exit(-1); }

  if (!(par.nprocx * par.nprocy * par.nprocz == par.nranks)) { std::cout << "Error: nranks != nprocx * nprocy * nprocz\n"; exit(-1); }

  //Get my process grid IDs
  par.pz = floor(par.myrank / (par.nprocx * par.nprocy));
  par.py = floor((par.myrank - par.nprocx * par.nprocy * par.pz) / par.nprocx);
  par.px = par.myrank - par.nprocx * par.nprocy * par.pz - par.nprocx * par.py;


    //Get my beginning and ending global indices; and domain sizes
    double nper;
    nper = floor(((double) params.nx_glob)/par.nprocx);
    par.i_beg = (int) round( nper* par.px    );
    par.i_end = (int) round( nper*(par.px+1) )-1;
    par.nx = par.i_end - par.i_beg + 1;
    par.nx_glob = params.nx_glob;

    if (ndims >=2)
    {
    nper = ((double) params.ny_glob)/par.nprocy;
    par.j_beg = (int) round( nper* par.py    );
    par.j_end = (int) round( nper*(par.py+1) )-1;
    par.ny = par.j_end - par.j_beg + 1;
    par.ny_glob = params.ny_glob;
    }

    if (ndims ==3)
    {
    nper = ((double) params.nz_glob)/par.nprocz;
    par.k_beg = (int) round( nper* par.pz    );
    par.k_end = (int) round( nper*(par.pz+1) )-1;
    par.nz = par.k_end - par.k_beg + 1;
    par.nz_glob = params.nz_glob;
    }

    // Determine my neighbors
    // x-dir
    int pxloc_neg = par.px-1;
    if (pxloc_neg < 0            ) pxloc_neg = pxloc_neg + par.nprocx;
    if (pxloc_neg > par.nprocx-1) pxloc_neg = pxloc_neg - par.nprocx;
    par.x_neigh(0) = (par.py + par.pz * par.nprocy) * par.nprocx + pxloc_neg;
    int pxloc_pos = par.px+1;
    if (pxloc_pos < 0            ) pxloc_pos = pxloc_pos + par.nprocx;
    if (pxloc_pos > par.nprocx-1) pxloc_pos = pxloc_pos - par.nprocx;
    par.x_neigh(1) = (par.py + par.pz * par.nprocy) * par.nprocx + pxloc_pos;

    // y-dir
    if (ndims>=2)
    {
    int pyloc_neg = par.py-1;
    if (pyloc_neg < 0            ) pyloc_neg = pyloc_neg + par.nprocy;
    if (pyloc_neg > par.nprocy-1) pyloc_neg = pyloc_neg - par.nprocy;
    par.y_neigh(0) = (pyloc_neg + par.pz * par.nprocy) * par.nprocx + par.px;
    int pyloc_pos = par.py+1;
    if (pyloc_pos < 0            ) pyloc_pos = pyloc_pos + par.nprocy;
    if (pyloc_pos > par.nprocy-1) pyloc_pos = pyloc_pos - par.nprocy;
    par.y_neigh(1) = (pyloc_pos + par.pz * par.nprocy) * par.nprocx + par.px;
    }

    // z-dir
    if (ndims == 3)
    {
    int pzloc_neg = par.pz-1;
    if (pzloc_neg < 0            ) pzloc_neg = pzloc_neg + par.nprocz;
    if (pzloc_neg > par.nprocz-1) pzloc_neg = pzloc_neg - par.nprocz;
    par.z_neigh(0) = (par.py + pzloc_neg * par.nprocy) * par.nprocx + par.px;
    int pzloc_pos = par.pz+1;
    if (pzloc_pos < 0            ) pzloc_pos = pzloc_pos + par.nprocz;
    if (pzloc_pos > par.nprocz-1) pzloc_pos = pzloc_pos - par.nprocz;
    par.z_neigh(1) = (par.py + pzloc_pos * par.nprocy) * par.nprocx + par.px;
    }

    // set halos
    par.halox = maxhalosize;
    par.haloy = maxhalosize;
    par.haloz = maxhalosize;

    // Debug output for the parallel decomposition
    if (1) {
      ierr = MPI_Barrier(MPI_COMM_WORLD);
      for (int rr=0; rr < par.nranks; rr++) {
        if (rr == par.myrank) {
          std::cout << "Hello! My rank is: " << par.myrank << "\n";
          std::cout << "My proc grid ID is: " << par.px << " , " << par.py << " , " << par.pz << "\n";
          std::cout << "I have: " << par.nx << " x " << par.ny << " x " << par.nz << " grid points." << "\n";
          std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg << " x " << par.k_beg << "\n";
          std::cout << "I end at index: " << par.i_end << " x " << par.j_end << " x " << par.k_end << "\n";
          std::cout << "My x neighbors are: " << par.x_neigh(0) << " " << par.x_neigh(1) << "\n";
          std::cout << "My y neighbors are: " << par.y_neigh(0) << " " << par.y_neigh(1) << "\n";
          std::cout << "My z neighbors are: " << par.z_neigh(0) << " " << par.z_neigh(1) << "\n";
        }
        ierr = MPI_Barrier(MPI_COMM_WORLD);
      }
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    }



  // Print out the values
  if (par.masterproc) {
    std::cout << "nx:         " << params.nx_glob    << "\n";
    std::cout << "ny:         " << params.nx_glob    << "\n";
    std::cout << "nz:         " << params.nx_glob    << "\n";

    std::cout << "halox:      " << par.halox    << "\n";
    std::cout << "haloy:      " << par.haloy    << "\n";
    std::cout << "haloz:      " << par.haloz    << "\n";

    std::cout << "dt:         " << params.dt         << "\n";
    std::cout << "Nsteps:     " << params.Nsteps     << "\n";
    std::cout << "Nout:       " << params.Nout       << "\n";
    std::cout << "outputName: " << params.outputName << "\n";

    std::cout << "nranks:     " << par.nranks        << "\n";
    std::cout << "nprocx:     " << par.nprocx        << "\n";
    std::cout << "nprocy:     " << par.nprocy        << "\n";
    std::cout << "nprocz:     " << par.nprocz        << "\n";

    std::cout << "xlen:       " << params.xlen       << "\n";
    std::cout << "ylen:       " << params.ylen       << "\n";
    std::cout << "zlen:       " << params.zlen       << "\n";
    std::cout << "xc:         " << params.xc         << "\n";
    std::cout << "yc:         " << params.yc         << "\n";
    std::cout << "zc:         " << params.zc         << "\n";
    std::cout << "etime:      " << params.etime      << "\n";
  }


};





#endif
