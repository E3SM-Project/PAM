
#include "pam_const.h"
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "DataManager.h"
#include "mmf_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "gcm_density_forcing.h"


int main(int argc, char** argv) {
  yakl::init();
  {
    yakl::timer_start("main");

    // Creates one PamCoupler object in static scope (no parameters assumes 1 thread)
    mmf_interface::init();

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real simTime        = config["simTime"    ].as<real>();
    int  crm_nx         = config["crm_nx"     ].as<int>();
    int  crm_ny         = config["crm_ny"     ].as<int>();
    int  nens           = config["nens"       ].as<int>();
    real xlen           = config["xlen"       ].as<real>();
    real ylen           = config["ylen"       ].as<real>();
    real dt_gcm         = config["dt_gcm"     ].as<real>();
    real dt_crm_phys    = config["dt_crm_phys"].as<real>();
    real outFreq        = config["outFreq"    ].as<real>();
    std::string coldata = config["column_data"].as<std::string>();

    // How to apply the density forcing: loose (force it like sam does); strict (enforce it strictly every time step)
    if (config["density_forcing"]) {
      mmf_interface::set_option<std::string>("density_forcing",config["density_forcing"].as<std::string>());
    } else {
      mmf_interface::set_option<std::string>("density_forcing","loose");
    }

    // Apply forcing at dynamics time step?  If false, it's applied at the CRM physics time step
    if (config["forcing_at_dycore_time_step"]) {
      mmf_interface::set_option<bool>( "forcing_at_dycore_time_step" , config["forcing_at_dycore_time_step"].as<bool>() );
    } else {
      mmf_interface::set_option<bool>( "forcing_at_dycore_time_step" , false );
    }

    // Store vertical coordinates
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int crm_nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in;
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    // Allocates the coupler state (density_dry, uvel, vvel, wvel, temp, vert grid, hydro background) for thread 0
    mmf_interface::allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;

    // Set physical constants for coupler at thread 0 using microphysics data
    mmf_interface::set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Set the vertical grid in the coupler
    mmf_interface::set_grid( xlen , ylen , zint_in );

    // This is for the dycore to pull out to determine how to do idealized test cases
    mmf_interface::set_option<std::string>( "standalone_input_file" , inFile );

    mmf_interface::register_and_allocate_array<real>("gcm_density_dry","GCM column dry density"     ,{crm_nz,nens});
    mmf_interface::register_and_allocate_array<real>("gcm_uvel"       ,"GCM column u-velocity"      ,{crm_nz,nens});
    mmf_interface::register_and_allocate_array<real>("gcm_vvel"       ,"GCM column v-velocity"      ,{crm_nz,nens});
    mmf_interface::register_and_allocate_array<real>("gcm_wvel"       ,"GCM column w-velocity"      ,{crm_nz,nens});
    mmf_interface::register_and_allocate_array<real>("gcm_temp"       ,"GCM column temperature"     ,{crm_nz,nens});
    mmf_interface::register_and_allocate_array<real>("gcm_water_vapor","GCM column water vapor mass",{crm_nz,nens});

    ///////////////////////////////////////////////////////////////////////////////
    // This is the end of the mmf_interface code. Using coupler directly from here
    ///////////////////////////////////////////////////////////////////////////////
    auto &coupler = mmf_interface::get_coupler();

    // Initialize the microphysics (registers the microphysics's tracers)
    micro .init( coupler );
    // Before calling dycore.init(coupler), be sure ALL TRACERS HAVE BEEN REGISTERED ALREADY in the coupler
    dycore.init( coupler );

    #ifdef PAM_STANDALONE
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name() << std::endl;
    #endif

    // Only for the idealized standalone driver; clearly not going to be used for the MMF driver
    dycore.init_idealized_state_and_tracers( coupler );

    auto gcm_rho_d = coupler.dm.get<real,2>("gcm_density_dry");
    auto gcm_uvel  = coupler.dm.get<real,2>("gcm_uvel"       );
    auto gcm_vvel  = coupler.dm.get<real,2>("gcm_vvel"       );
    auto gcm_wvel  = coupler.dm.get<real,2>("gcm_wvel"       );
    auto gcm_temp  = coupler.dm.get<real,2>("gcm_temp"       );
    auto gcm_rho_v = coupler.dm.get<real,2>("gcm_water_vapor");

    auto rho_d = coupler.dm.get<real const,4>("density_dry" );
    auto uvel  = coupler.dm.get<real const,4>("uvel"        );
    auto vvel  = coupler.dm.get<real const,4>("vvel"        );
    auto wvel  = coupler.dm.get<real const,4>("wvel"        );
    auto temp  = coupler.dm.get<real const,4>("temp"        );
    auto rho_v = coupler.dm.get<real const,4>("water_vapor" );

    // Compute a column to force the model with by averaging the columns at init
    parallel_for( Bounds<4>(crm_nz,crm_ny,crm_nx,nens) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int iens) {
      gcm_rho_d(k,iens) = 0;
      gcm_uvel (k,iens) = 0;
      gcm_vvel (k,iens) = 0;
      gcm_wvel (k,iens) = 0;
      gcm_temp (k,iens) = 0;
      gcm_rho_v(k,iens) = 0;
    });
    real r_nx_ny = 1._fp / (crm_nx*crm_ny);  // Avoid costly divisions
    parallel_for( Bounds<4>(crm_nz,crm_ny,crm_nx,nens) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int iens) {
      yakl::atomicAdd( gcm_rho_d(k,iens) , rho_d(k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_uvel (k,iens) , uvel (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_vvel (k,iens) , vvel (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_wvel (k,iens) , wvel (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_temp (k,iens) , temp (k,j,i,iens) * r_nx_ny );
      yakl::atomicAdd( gcm_rho_v(k,iens) , rho_v(k,j,i,iens) * r_nx_ny );
    });

    perturb_temperature( coupler , 0 );

    // Now that we have an initial state, define hydrostasis for each ensemble member
    coupler.update_hydrostasis( coupler.compute_pressure_array() );

    bool forcing_at_dycore_time_step = coupler.get_option<bool>("forcing_at_dycore_time_step");

    if (forcing_at_dycore_time_step) {
      coupler.add_dycore_function( gcm_density_forcing          );
      coupler.add_dycore_function( apply_gcm_forcing_tendencies );
    }

    real etime_gcm = 0;
    int numOut = 0;

    if (outFreq >= 0) dycore.output( coupler , etime_gcm );

    while (etime_gcm < simTime) {
      if (etime_gcm + dt_gcm > simTime) { dt_gcm = simTime - etime_gcm; }

      compute_gcm_forcing_tendencies( coupler , dt_gcm );

      real etime_crm = 0;
      real simTime_crm = dt_gcm;
      real dt_crm = dt_crm_phys;
      while (etime_crm < simTime_crm) {
        if (dt_crm == 0.) { dt_crm = dycore.compute_time_step(coupler); }
        if (etime_crm + dt_crm > simTime_crm) { dt_crm = simTime_crm - etime_crm; }

        yakl::timer_start("micro");
        micro.timeStep( coupler , dt_crm , etime_crm );
        yakl::timer_stop("micro");

        yakl::timer_start("dycore");
        dycore.timeStep( coupler , dt_crm , etime_crm );
        yakl::timer_stop("dycore");

        if (! forcing_at_dycore_time_step) {
          gcm_density_forcing         (coupler , dt_crm);
          apply_gcm_forcing_tendencies(coupler , dt_crm);
        }

        etime_crm += dt_crm;
        etime_gcm += dt_crm;
        if (outFreq >= 0. && etime_gcm / outFreq >= numOut+1) {
          std::cout << "Etime , dt_crm , maxw: " << etime_gcm << " , "
                                                 << dt_crm    << " , "
                                                 << yakl::intrinsics::maxval(wvel) << "\n";
          yakl::timer_start("output");
          dycore.output( coupler , etime_gcm );
          yakl::timer_stop("output");
          numOut++;
        }
      }
    }

    std::cout << "Elapsed Time: " << etime_gcm << "\n";

    dycore.finalize( coupler );

    mmf_interface::finalize();

    yakl::timer_stop("main");
  }
  yakl::finalize();
}


