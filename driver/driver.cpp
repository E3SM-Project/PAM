
#include "stdlib.h"
#include <iostream>
#include <string>
#include "const.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"
#include "Initializer.h"
#include "TimeIntegrator.h"
#include "FileIO.h"
#include "Exchange.h"

int main(int argc, char** argv) {

  yakl::init();

  {
    // Create the model objects
    realArr         state;
    Domain         dom;
    Parallel       par;
    FileIO         io;
    Exchange       exch;
    TimeIntegrator tint;

    int nstep = 0;

    // Initialize MPI and read the input file
    initialize_mpi( &argc , &argv , par );

    // Default input file is "input.txt" unless the user passes in another file
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    readParamsFile(inFile, dom, par, io);

    // Initialize the model
    initialize(state, dom, par, exch, tint);

    // Output the initial model state
    io.outputInit(state, dom, par);

    // while (dom.etime < dom.simLength) {
    //   if (dom.etime + dom.dt > dom.simLength) { dom.dt = dom.simLength - dom.etime; }

    //   yakl::fence();
    //   tint.stepForward(state, dom, exch, par);
    //   yakl::fence();

    //   dom.etime += dom.dt;
    //   if (par.masterproc && nstep%100 == 0) {std::cout << dom.etime << "\n";}

    //   io.output(state, dom, par);

    //   nstep += 1;
    // }
  }

  int ierr = MPI_Finalize();

}
