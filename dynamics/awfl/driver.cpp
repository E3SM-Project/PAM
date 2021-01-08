
#include "const.h"
#include "Spatial_euler3d_cons_expl_cart_fv_Agrid.h"
#include "Temporal_ader.h"
#include "Profiles.h"
#include "Microphysics_saturation_adjustment.h"
#include "DataManager.h"

// Define the Spatial operator based on constants from the Temporal operator
typedef Spatial_euler3d_cons_expl_cart_fv_Agrid<nTimeDerivs,timeAvg,nAder> Spatial;

// Define the Temporal operator based on the Spatial operator
typedef Temporal_ader<Spatial> Dycore;

typedef Microphysics_saturation_adjustment Microphysics;

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
    int numOut = 0;

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    // Initialize the dycore and the microphysics
    dycore.init( inFile , micro.num_tracers+1 , dm );
    micro .init( inFile , dycore              , dm );

    // Initialize the dry state
    dycore.init_state_and_tracers( dm , micro );

    real etime = 0;

    dycore.space_op.output( dm , micro , etime );
    
    while (etime < simTime) {
      real dt = dycore.compute_time_step( 0.8 , dm , micro );
      if (etime + dt > simTime) { dt = simTime - etime; }
      dycore.timeStep( dm , micro , dt );
      micro.timeStep( dm , dt );
      etime += dt;
      if (etime / outFreq >= numOut+1) {
        std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
        dycore.space_op.output( dm , micro , etime );
        numOut++;
      }
    }

    dycore.timeStep( dm , micro , 1. );

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( dm );

  }
  yakl::finalize();
}



