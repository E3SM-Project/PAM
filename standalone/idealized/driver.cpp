
#include "pam_const.h"
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "DataManager.h"


int main(int argc, char** argv) {
  yakl::init();
  {
    yakl::timer_start("main");

    pam::PamCoupler coupler;

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real simTime   = config["simTime"].as<real>();
    real outFreq   = config["outFreq"].as<real>();
    int  nx        = config["nx"     ].as<int>();
    int  ny        = config["ny"     ].as<int>();
    int  nens      = config["nens"   ].as<int>();
    real xlen      = config["xlen"   ].as<real>();
    real ylen      = config["ylen"   ].as<real>();
    real dtphys_in = config["dtphys" ].as<real>();

    // Store vertical coordinates
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in("zint_in",nz+1);
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    // Use microphysics gas constants values in the coupler
    coupler.set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Allocate coupler state
    coupler.allocate_coupler_state( nz , ny , nx , nens );

    // Set the horizontal domain lengths and the vertical grid in the coupler
    coupler.set_grid( xlen , ylen , zint_in );

    int numOut = 0;

    // This is for the dycore to pull out to determine how to do idealized test cases
    coupler.add_note( "standalone_input_file" , inFile );

    // Initialize the dycore and the microphysics
    dycore.init( micro.get_num_tracers() , coupler );
    micro .init( dycore , coupler );

    #ifdef PAM_STANDALONE
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name() << std::endl;
    #endif

    // Only for the idealized standalone driver; clearly not going to be used for the MMF driver
    dycore.init_idealized_state_and_tracers( coupler );

    // Now that we have an initial state, define hydrostasis for each ensemble member
    coupler.update_hydrostasis( coupler.compute_pressure_array() );

    real etime = 0;

    if (outFreq >= 0) dycore.output( coupler , etime );

    real dtphys = dtphys_in;
    while (etime < simTime) {
      if (dtphys_in == 0.) { dtphys = dycore.compute_time_step(coupler); }
      if (etime + dtphys > simTime) { dtphys = simTime - etime; }

      yakl::timer_start("micro");
      micro.timeStep( coupler , dtphys );
      yakl::timer_stop("micro");

      yakl::timer_start("dycore");
      dycore.timeStep( coupler , dtphys );
      yakl::timer_stop("dycore");

      etime += dtphys;
      if (outFreq >= 0. && etime / outFreq >= numOut+1) {
        std::cout << "Etime , dtphys: " << etime << " , " << dtphys << "\n";
        yakl::timer_start("output");
        dycore.output( coupler , etime );
        yakl::timer_stop("output");
        numOut++;
      }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( coupler );

    yakl::timer_stop("main");
  }
  yakl::finalize();
}
