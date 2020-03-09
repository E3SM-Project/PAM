#ifndef _PARSER_H_
#define _PARSER_H_

#include <fstream>
#include <sstream>
#include "string.h"

class Parameters
{
public:
  int nx = -1;
  int ny = -1;
  int nz = -1;

  int Nsteps = -1;
  int Nout = -1;
  real dt = -1.;
  std::string outputName = "output.nc";

  int nranks = -1;
  int myrank = -1;
  int masterproc = -1;

  real etime;
//THESE ARE REALLY SPECIFIC TO UNIFORM RECT GEOM...
  real xlen, ylen, zlen;
  real xc, yc, zc;
};




void readParamsFile(std::string inFile, Parameters &params) {


// EVENTUALLY THIS SHOULD BE SET BY THE INITIAL CONDITION CHOICE!


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
      if      ( !strcmp( "nx"         , key.c_str() ) ) { ssVal >> params.nx         ; }
      else if ( !strcmp( "ny"         , key.c_str() ) ) { ssVal >> params.ny         ; }
      else if ( !strcmp( "nz"         , key.c_str() ) ) { ssVal >> params.nz         ; }

      else if ( !strcmp( "dt"         , key.c_str() ) ) { ssVal >> params.dt         ; }
      else if ( !strcmp( "Nsteps"     , key.c_str() ) ) { ssVal >> params.Nsteps     ; }
      else if ( !strcmp( "Nout"       , key.c_str() ) ) { ssVal >> params.Nout       ; }
      else if ( !strcmp( "outputName" , key.c_str() ) ) { ssVal >> params.outputName ; }
      else {
        std::cout << "Error: key " << key << " not understood in file " << inFile << "\n";
      }
    }
  }

  // Test to make sure all required values were initialized
  if (params.nx     == -1) { std::cout << "Error: key " << "nx"        << " not set.\n"; exit(-1); }
  if (params.ny     == -1) { std::cout << "Error: key " << "ny"        << " not set.\n"; exit(-1); }
  if (params.nz     == -1) { std::cout << "Error: key " << "nz"        << " not set.\n"; exit(-1); }
  if (params.dt     == -1) { std::cout << "Error: key " << "dt"        << " not set.\n"; exit(-1); }
  if (params.Nsteps == -1) { std::cout << "Error: key " << "Nsteps"    << " not set.\n"; exit(-1); }
  if (params.Nout   == -1) { std::cout << "Error: key " << "Nout"      << " not set.\n"; exit(-1); }

  // Print out the values
  if (params.masterproc) {
    std::cout << "nx:         " << params.nx         << "\n";
    std::cout << "ny:         " << params.ny         << "\n";
    std::cout << "nz:         " << params.nz         << "\n";

    std::cout << "dt:         " << params.dt         << "\n";
    std::cout << "Nsteps:     " << params.Nsteps     << "\n";
    std::cout << "Nout:       " << params.Nout       << "\n";
    std::cout << "outputName: " << params.outputName << "\n";


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
