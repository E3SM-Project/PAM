
#include "pam_const.h"
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "DataManager.h"


int main(int argc, char** argv) {
  yakl::init();
  {
    yakl::timer_start("main");

    // Creates one PamCoupler object in static scope (no parameters assumes 1 thread)
    mmf_interface_init();

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real simTime     = config["simTime"    ].as<real>();
    int  crm_nx      = config["crm_nx"     ].as<int>();
    int  crm_ny      = config["crm_ny"     ].as<int>();
    int  nens        = config["nens"       ].as<int>();
    real xlen        = config["xlen"       ].as<real>();
    real ylen        = config["ylen"       ].as<real>();
    real dt_gcm      = config["dt_gcm"     ].as<real>();
    real dt_crm_phys = config["dt_crm_phys"].as<real>();
    real outFreq     = config["outFreq"    ].as<real>();

    // Store vertical coordinates
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in("zint_in",nz+1);
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    // Allocates the coupler state (density_dry, uvel, vvel, wvel, temp, vert grid, hydro background) for thread 0
    mmf_allocate_coupler_state( nz , crm_ny , crm_nx , nens );

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    // Set physical constants for coupler at thread 0 using microphysics data
    mmf_set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Set the vertical grid in the coupler
    mmf_set_vertical_grid( zint_in );

    int numOut = 0;

    // This is for the dycore to pull out to determine how to do idealized test cases
    coupler.add_note( "standalone_input_file" , inFile );

    // Initialize the microphysics (registers the microphysics's tracers)
    micro .init( coupler );
    // Before calling dycore.init(coupler), be sure ALL TRACERS HAVE BEEN REGISTERED ALREADY in the coupler
    dycore.init( coupler );

    #ifdef PAM_STANDALONE
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name() << std::endl;
    #endif

    // TODO: Initialize CRM data using supercell background state

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
