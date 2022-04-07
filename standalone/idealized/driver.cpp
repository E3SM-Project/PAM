
#include "pam_const.h"
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "DataManager.h"


int main(int argc, char** argv) {
  yakl::init();
  {
    using yakl::intrinsics::abs;
    using yakl::intrinsics::maxval;
    yakl::timer_start("main");

    pam::PamCoupler coupler;

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real        simTime      = config["simTime" ].as<real>(0.0);
    real        simSteps     = config["simSteps"].as<int>(0);
    int         crm_nx       = config["crm_nx"  ].as<int>();
    int         crm_ny       = config["crm_ny"  ].as<int>();
    int         nens         = config["nens"    ].as<int>();
    real        xlen         = config["xlen"    ].as<real>();
    real        ylen         = config["ylen"    ].as<real>();
    real        dtphys_in    = config["dtphys"  ].as<real>();
    std::string vcoords_file = config["vcoords" ].as<std::string>("vcoords_none.nc");
    bool        use_coupler_hydrostasis = config["use_coupler_hydrostasis"].as<bool>(false);

    // Read vertical coordinates
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int crm_nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in("zint_in",crm_nz+1);
    nc.read(zint_in,"vertical_interfaces");
    nc.close();      

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;
    SGS          sgs;

    // Use microphysics gas constants values in the coupler
    coupler.set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Allocate coupler state
    coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );
    
    // Set the horizontal domain lengths and the vertical grid in the coupler
    coupler.set_grid( xlen , ylen , zint_in );

    // This is for the dycore to pull out to determine how to do idealized test cases
    coupler.set_option<std::string>( "standalone_input_file" , inFile );

    micro .init( coupler );
    sgs   .init( coupler );
    dycore.init( coupler ); // Dycore should initialize its own state here

    #ifdef PAM_STANDALONE
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name () << std::endl;
      std::cout << "SGS   : " << sgs   .sgs_name   () << std::endl;
      std::cout << "\n";
    #endif

    // Now that we have an initial state, define hydrostasis for each ensemble member
    if (use_coupler_hydrostasis) coupler.update_hydrostasis( coupler.compute_pressure_array() );

    real etime = 0;
    if (simTime == 0.0) {  simTime = simSteps * dtphys_in; }
    
    real dtphys = dtphys_in;
    while (etime < simTime) {
      if (dtphys_in == 0.) { dtphys = dycore.compute_time_step(coupler); }
      if (etime + dtphys > simTime) { dtphys = simTime - etime; }

      yakl::timer_start("sgs");
      sgs.timeStep( coupler , dtphys );
      yakl::timer_stop("sgs");

      yakl::timer_start("micro");
      micro.timeStep( coupler , dtphys );
      yakl::timer_stop("micro");

      yakl::timer_start("dycore");
      dycore.timeStep( coupler , dtphys );
      yakl::timer_stop("dycore");

      etime += dtphys;
      real maxw = maxval(abs(coupler.dm.get_collapsed<real const>("wvel")));
      std::cout << "Etime , dtphys, maxw: " << etime  << " , " 
                                            << dtphys << " , "
                                            << std::setw(10) << maxw << "\n";
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    dycore.finalize( coupler );

    yakl::timer_stop("main");
  }
  yakl::finalize();
}
