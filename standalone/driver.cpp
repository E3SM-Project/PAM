
#include "pam_const.h"
#include "coupler_state.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "DataManager.h"


int main(int argc, char** argv) {
  yakl::init();
  {

    DataManager dm;

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real simTime = config["simTime"].as<real>();
    real outFreq = config["outFreq"].as<real>();
    int  nx      = config["nx"     ].as<int>();
    int  ny      = config["ny"     ].as<int>();
    int  nens    = config["nens"   ].as<int>();
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int nz = nc.getDimSize("num_interfaces") - 1;
    nc.close();


    int numOut = 0;

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    allocate_coupler_state( nz , ny , nx , nens , dm );

    // Initialize the dycore and the microphysics
    dycore.init( inFile , micro.num_tracers , dm );
    micro .init( inFile , dycore            , dm );

    // Initialize the dry state
    dycore.init_state_and_tracers( dm , micro );

    real etime = 0;

    dycore.output( dm , micro , etime );

    while (etime < simTime) {
      real dt = dycore.compute_time_step( 0.8 , dm , micro );
      if (etime + dt > simTime) { dt = simTime - etime; }
      dycore.timeStep( dm , micro , dt );
      dycore.convert_dynamics_to_coupler_state( dm , micro );
      micro.timeStep( dm , dt );
      dycore.convert_coupler_state_to_dynamics( dm , micro );
      etime += dt;
      if (etime / outFreq >= numOut+1) {
        std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
        dycore.output( dm , micro , etime );
        numOut++;
      }
    }

    dycore.timeStep( dm , micro , 1. );

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( dm );

  }
  yakl::finalize();
}


