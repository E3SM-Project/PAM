
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
    coupler.set_gas_constants(micro.constants.R_d , micro.constants.R_v);

    // Allocate coupler state
    coupler.allocate_coupler_state( nz , ny , nx , nens );

    // Set the vertical grid in the coupler
    coupler.set_vertical_grid( zint_in );

    int numOut = 0;

    // Initialize the dycore and the microphysics
    dycore.init( inFile , ny , nx , nens , xlen , ylen , micro.get_num_tracers() , coupler.dm );
    micro .init( inFile , ny , nx , nens , dycore , coupler.dm );

    #ifdef PAM_STANDALONE
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name() << std::endl;
    #endif

    // Initialize the dry state
    dycore.init_state_and_tracers( coupler.dm , micro );

    // Now that we have an initial state, define hydrostasis for each ensemble member
    auto press = pam::compute_pressure_array( coupler.dm , micro.constants.R_d , micro.constants.R_v );
    coupler.update_hydrostasis( press );

    real etime = 0;

    if (outFreq >= 0) dycore.output( coupler.dm , micro , etime );

    real dtphys = dtphys_in;
    while (etime < simTime) {
      if (dtphys_in == 0.) { dtphys = dycore.compute_time_step(coupler.dm, micro); }
      if (etime + dtphys > simTime) { dtphys = simTime - etime; }

      yakl::timer_start("micro");
      micro.timeStep( coupler.dm , dtphys );
      yakl::timer_stop("micro");

      yakl::timer_start("dycore");
      dycore.timeStep( coupler.dm , micro , dtphys );
      yakl::timer_stop("dycore");

      etime += dtphys;
      if (outFreq >= 0. && etime / outFreq >= numOut+1) {
        std::cout << "Etime , dtphys: " << etime << " , " << dtphys << "\n";
        yakl::timer_start("output");
        dycore.output( coupler.dm , micro , etime );
        yakl::timer_stop("output");
        numOut++;
      }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( coupler.dm );

    yakl::timer_stop("main");
  }
  yakl::finalize();
}
