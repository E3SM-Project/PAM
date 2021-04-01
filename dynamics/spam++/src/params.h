#ifndef _PARAMS_H_
#define _PARAMS_H_

#include <fstream>
#include <sstream>
#include "string.h"
#include "mpi.h"
#include "parallel.h"
#include "common.h"





int ij_to_l(int i, int j, int nx)
{
  return j * nx + i;
}

int wrap(int i, int nx)
{
  if (i < 0) {return i + nx;}
  else if (i > nx-1) {return i - nx;}
  else {return i;}
}




void readParamsFile(std::string inFile, Parameters &params, Parallel &par) {

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
      else if ( !strcmp( "nz"         , key.c_str() ) ) { ssVal >> params.nz    ; }

      else if ( !strcmp( "dt"         , key.c_str() ) ) { ssVal >> params.dt         ; }
      else if ( !strcmp( "Nsteps"     , key.c_str() ) ) { ssVal >> params.Nsteps     ; }
      else if ( !strcmp( "Nout"       , key.c_str() ) ) { ssVal >> params.Nout       ; }
      else if ( !strcmp( "Nstat"      , key.c_str() ) ) { ssVal >> params.Nstat      ; }
      else if ( !strcmp( "TStype"     , key.c_str() ) ) { ssVal >> params.TStype     ; }
      else if ( !strcmp( "cfl"        , key.c_str() ) ) { ssVal >> params.cfl        ; }
      else if ( !strcmp( "outputName" , key.c_str() ) ) { ssVal >> params.outputName ; }

      else if ( !strcmp( "nprocx"     , key.c_str() ) ) { ssVal >> par.nprocx        ; }
      else if ( !strcmp( "nprocy"     , key.c_str() ) ) { ssVal >> par.nprocy        ; }

      else if ( !strcmp( "xbnd"         , key.c_str() ) ) { ssVal >> params.xbnd    ; }
      else if ( !strcmp( "ybnd"         , key.c_str() ) ) { ssVal >> params.ybnd    ; }


      //else {
      //  std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      //}
    }
  }

// ADD TSType fixes here

  // Test to make sure all required values were initialized
  if (params.nx_glob == -1) { std::cout << "Error: key " << "nx"          << " not set.\n"; exit(-1); }
  if (params.ny_glob == -1) { std::cout << "Error: key " << "ny"          << " not set.\n"; exit(-1); }
  if (params.nz == -1) { std::cout << "Error: key " << "nz"          << " not set.\n"; exit(-1); }
  if (params.dt      == -1) { std::cout << "Error: key " << "dt"          << " not set.\n"; exit(-1); }
  if (params.Nsteps  == -1) { std::cout << "Error: key " << "Nsteps"      << " not set.\n"; exit(-1); }
  if (params.Nout    == -1) { std::cout << "Error: key " << "Nout"        << " not set.\n"; exit(-1); }
  if (par.nprocx     == -1) { std::cout << "Error: key " << "nprocx"      << " not set.\n"; exit(-1); }
  if (par.nprocy     == -1) { std::cout << "Error: key " << "nprocy"      << " not set.\n"; exit(-1); }

  if (!(par.nprocx * par.nprocy == par.nranks)) { std::cout << "Error: nranks != nprocx * nprocy\n"; exit(-1); }

//FIX THIS....
  //Get my process grid IDs
//  par.pz = floor(par.myrank / (par.nprocx * par.nprocy));
//  par.py = floor((par.myrank - par.nprocx * par.nprocy * par.pz) / par.nprocx);
//  par.px = par.myrank - par.nprocx * par.nprocy * par.pz - par.nprocx * par.py;

  par.py = floor(par.myrank / par.nprocx);
  par.px = par.myrank - par.nprocx * par.py;


    //Get my beginning and ending global indices; and domain sizes
    double nper;
    nper = ((double) params.nx_glob)/par.nprocx;
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
    
    par.nz = params.nz;

    // Determine my neighbors
    // x-dir
    par.x_neigh(0) = ij_to_l(wrap(par.px-1,par.nprocx), par.py, par.nprocx);
    par.x_neigh(1) = ij_to_l(wrap(par.px+1,par.nprocx), par.py, par.nprocx);

    // y-dir
    if (ndims==2)
    {
      par.y_neigh(0) = ij_to_l(par.px, wrap(par.py-1,par.nprocy), par.nprocx);
      par.y_neigh(1) = ij_to_l(par.px, wrap(par.py+1,par.nprocy), par.nprocx);
      
    par.ll_neigh = ij_to_l(wrap(par.px-1,par.nprocx), wrap(par.py-1,par.nprocy), par.nprocx); 
    par.lr_neigh = ij_to_l(wrap(par.px+1,par.nprocx), wrap(par.py-1,par.nprocy), par.nprocx); 
    par.ur_neigh = ij_to_l(wrap(par.px+1,par.nprocx), wrap(par.py+1,par.nprocy), par.nprocx); 
    par.ul_neigh = ij_to_l(wrap(par.px-1,par.nprocx), wrap(par.py+1,par.nprocy), par.nprocx);
    }
    
    // set boundaries
    if (!strcmp(params.xbnd.c_str(),"periodic")) {par.xbnd = BND_TYPE::PERIODIC;}
    if (!strcmp(params.xbnd.c_str(),"none")) {par.xbnd = BND_TYPE::NONE;}
    if (!strcmp(params.ybnd.c_str(),"periodic")) {par.ybnd = BND_TYPE::PERIODIC;}
    if (!strcmp(params.ybnd.c_str(),"none")) {par.ybnd = BND_TYPE::NONE;}
        
    // set halos
    par.halox = maxhalosize;
    par.haloy = maxhalosize;
        
    // Debug output for the parallel decomposition
    if (1) {
      ierr = MPI_Barrier(MPI_COMM_WORLD);
      for (int rr=0; rr < par.nranks; rr++) {
        if (rr == par.myrank) {
          std::cout << "Hello! My rank is: " << par.myrank << "\n";
          std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
          std::cout << "I have: " << par.nx << " x " << par.ny << " grid points." << "\n";
          std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg << "\n";
          std::cout << "I end at index: " << par.i_end << " x " << par.j_end << "\n";
          std::cout << "My x neighbors are: " << par.x_neigh(0) << " " << par.x_neigh(1) << "\n";
          std::cout << "My y neighbors are: " << par.y_neigh(0) << " " << par.y_neigh(1) << "\n";
          std::cout << "My corner neighbors are: " << par.ll_neigh << " " << par.ul_neigh << " " << par.ur_neigh << " " << par.lr_neigh << "\n";
        }
        ierr = MPI_Barrier(MPI_COMM_WORLD);
      }
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    }



  // Print out the values
  if (par.masterproc) {
    std::cout << "nx:         " << params.nx_glob    << "\n";
    std::cout << "ny:         " << params.ny_glob    << "\n";
    std::cout << "nl:         " << params.nz    << "\n";
    std::cout << "ni:         " << params.nz+1    << "\n";

    std::cout << "halox:      " << par.halox    << "\n";
    std::cout << "haloy:      " << par.haloy    << "\n";

    std::cout << "xbnd:      " << params.xbnd    << "\n";
    std::cout << "ybnd:      " << params.ybnd    << "\n";
    
    std::cout << "dt:         " << params.dt         << "\n";
    std::cout << "Nsteps:     " << params.Nsteps     << "\n";
    std::cout << "Nout:       " << params.Nout       << "\n";
    std::cout << "outputName: " << params.outputName << "\n";

    std::cout << "nranks:     " << par.nranks        << "\n";
    std::cout << "nprocx:     " << par.nprocx        << "\n";
    std::cout << "nprocy:     " << par.nprocy        << "\n";

    std::cout << "etime:      " << params.etime      << "\n";

//THESE ARE UNIFORM GEOMETRY SPECIFIC....
//SO WRAP UP INTO A PRINT ROUTINE FOR GEOMETRY?
    std::cout << "xlen:       " << params.xlen       << "\n";
    std::cout << "ylen:       " << params.ylen       << "\n";
    std::cout << "zlen:       " << params.zlen       << "\n";
    std::cout << "xc:         " << params.xc         << "\n";
    std::cout << "yc:         " << params.yc         << "\n";
    std::cout << "zc:         " << params.zc         << "\n";
  }


};





#endif
