#pragma once

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


void readParamsFile(std::string inFile, ModelParameters &params, Parallel &par, int nz) {
  
  //Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  // Determine if I'm the master process
  if (par.myrank == 0) { par.masterproc = 1;}
  else { par.masterproc = 0; }
  
  int ierr;

 params.nx_glob = config["crm_nx"].as<int>();
 params.ny_glob = config["crm_ny"].as<int>();
 params.nens = config["nens"].as<int>();
 params.nz = nz;

 par.nprocx = config["nprocx"].as<int>();
 par.nprocy = config["nprocy"].as<int>();

 if (!(par.nprocx * par.nprocy == par.nranks)) { 
   std::cout << par.nprocx << " " << par.nprocy << " " << par.nranks << "\n"; exit(-1); 

   std::cout << "Error: nranks != nprocx * nprocy\n"; exit(-1); 
 }



  params.dtphys = config["dtphys"].as<real>();
  params.crm_per_phys = config["crm_per_phys"].as<int>(0);
  params.Nout = config["outSteps"].as<int>(0);
  params.Nstat = config["statSteps"].as<int>(0);
  int simSteps = config["simSteps"].as<int>(0);
  if (not (params.dtphys > 0.0_fp) or  not (simSteps > 0) or not (params.crm_per_phys > 0)  or not (params.Nout > 0) or not (params.Nstat > 0)) 
  {endrun("spam++ must use step-based time control logic ie set simSteps >0, dtphys>0, crm_per_phys >0, outSteps >0, statSteps >0");}
  
  params.dtcrm =  params.dtphys / params.crm_per_phys;
  params.Nsteps = simSteps * params.crm_per_phys;

  params.tstype = config["tstype"].as<std::string>();

  params.outputName = config["out_prefix"].as<std::string>("output") + "_dycore" + ".nc";

  //Get my process grid IDs
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
    par.nens = params.nens;

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
    // if (!strcmp(params.xbnd.c_str(),"periodic")) {par.xbnd = BND_TYPE::PERIODIC;}
    // if (!strcmp(params.xbnd.c_str(),"none")) {par.xbnd = BND_TYPE::NONE;}
    // if (!strcmp(params.ybnd.c_str(),"periodic")) {par.ybnd = BND_TYPE::PERIODIC;}
    // if (!strcmp(params.ybnd.c_str(),"none")) {par.ybnd = BND_TYPE::NONE;}
    par.xbnd = BND_TYPE::PERIODIC;
    par.ybnd = BND_TYPE::PERIODIC;
    
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
    std::cout << "nens:         " << params.nens    << "\n";

    std::cout << "halox:      " << par.halox    << "\n";
    std::cout << "haloy:      " << par.haloy    << "\n";
    
    std::cout << "dtcrm:         " << params.dtcrm         << "\n";
    std::cout << "dtphys:         " << params.dtphys         << "\n";
    std::cout << "Nsteps:     " << params.Nsteps     << "\n";
    std::cout << "Nout:       " << params.Nout       << "\n";
    std::cout << "outputName: " << params.outputName << "\n";

    std::cout << "nranks:     " << par.nranks        << "\n";
    std::cout << "nprocx:     " << par.nprocx        << "\n";
    std::cout << "nprocy:     " << par.nprocy        << "\n";

    std::cout << "xlen:       " << params.xlen       << "\n";
    std::cout << "ylen:       " << params.ylen       << "\n";
    std::cout << "xc:         " << params.xc         << "\n";
    std::cout << "yc:         " << params.yc         << "\n";
  }


};