
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "pam_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"
#include "broadcast_initial_gcm_column.h"
#include "output.h"
#include "supercell_init.h"
#include <iostream>
#include <chrono>
#include "scream_cxx_interface_finalize.h"
#include "pamc_init.h"


int main(int argc, char** argv) {
  MPI_Init( &argc , &argv );
  yakl::init();
  {

    using yakl::intrinsics::abs;
    using yakl::intrinsics::maxval;
    yakl::timer_start("main");

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }

    auto        nsteps_gcm     = config["nsteps_gcm"].as<int>();
    auto        crm_nx       = config["crm_nx"  ].as<int>();
    auto        crm_ny       = config["crm_ny"  ].as<int>();
    int         nens         = config["nens"    ].as<int>();
    auto dt_gcm         = config["dt_gcm"     ].as<real>();
    auto dt_crm_phys    = config["dt_crm_phys"].as<real>();
    std::string vcoords_file = config["vcoords" ].as<std::string>();
    auto        use_coupler_hydrostasis = config["use_coupler_hydrostasis"].as<bool>(false);
    auto out_freq                = config["out_freq"   ].as<real>(-1.);
    auto out_prefix              = config["out_prefix" ].as<std::string>("test");
    auto        inner_mpi = config["inner_mpi"].as<bool>(false);

//ADD CHECK THAT dt_gcm is divisible by both dt_crm_phys and out_freq
    int nsteps_crm_phys = dt_gcm / dt_crm_phys;

    //set xlen, ylen, zlen based on init cond
//THIS NEEDS TO BE FIXED FOR PAM-A
    real xlen, ylen, zlen;
    set_domain_sizes(config, xlen, ylen, zlen);


      // Read vertical coordinates
//THIS IS BROKEN FOR PARALLEL IO CASE- maybe this is okay ie switch entirely to standard netcdf?
    real1d zint_in;
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    int crm_nz = nc.getDimSize("num_interfaces") - 1;
    zint_in = real1d("zint_in",crm_nz+1);
    nc.read(zint_in,"vertical_interfaces");
    nc.close();
//TODO: Eventually add support for ensemble dependent vertical grids

    int nranks;
    int myrank;
    MPI_Comm_size( MPI_COMM_WORLD , &nranks );
    MPI_Comm_rank( MPI_COMM_WORLD , &myrank );
    bool mainproc = (myrank == 0);

    auto &coupler = pam_interface::get_coupler();

    coupler.set_option<real>("gcm_physics_dt",dt_gcm);
    coupler.set_option<real>("crm_dt",dt_crm_phys);

    //this partitions the domain if INNER_MPI is set, otherwise it does nothing
    if (inner_mpi)
    {partition_domain(inFile, crm_nx, crm_ny);}

    // Allocates the coupler state (density_dry, uvel, vvel, wvel, temp, vert grid, hydro background) for thread 0
    coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

    // Set the vertical grid in the coupler
    coupler.set_grid( xlen , ylen , zint_in );

    // NORMALLY THIS WOULD BE DONE INSIDE THE CRM, BUT WE'RE USING CONSTANTS DEFINED BY THE CRM MICRO SCHEME
    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;
    SGS          sgs;

    micro .init( coupler );
    sgs   .init( coupler );
    yakl::timer_start("dycore");
    dycore.init( coupler ); // Dycore should initialize its own state here
    yakl::timer_stop("dycore");

#if PAM_SCREAM_USE_CXX
    if ( (micro.micro_name() == "p3" || sgs.sgs_name() == "shoc") && crm_nz != PAM_NLEV ) {
      endrun("ERROR: Running with a different number of vertical levels than compiled for");
    }
#endif

//ELIMINATE THIS LOGIC!
    // // Compute a supercell initial column
    // real1d rho_d_col("rho_d_col",crm_nz);
    // real1d uvel_col ("uvel_col" ,crm_nz);
    // real1d vvel_col ("vvel_col" ,crm_nz);
    // real1d wvel_col ("wvel_col" ,crm_nz);
    // real1d temp_col ("temp_col" ,crm_nz);
    // real1d rho_v_col("rho_v_col",crm_nz);
    //
    // auto R_d  = coupler.get_option<real>("R_d" );
    // auto R_v  = coupler.get_option<real>("R_v" );
    // auto grav = coupler.get_option<real>("grav");
    // supercell_init( zint_in, rho_d_col, uvel_col, vvel_col, wvel_col, temp_col, rho_v_col, R_d , R_v, grav );
    //
    // // Set the GCM column data for each ensemble to the supercell initial state
    // auto &dm = coupler.get_data_manager_device_readwrite();
    //
    // auto gcm_rho_d = dm.get<real,2>("gcm_density_dry");
    // auto gcm_uvel  = dm.get<real,2>("gcm_uvel"       );
    // auto gcm_vvel  = dm.get<real,2>("gcm_vvel"       );
    // auto gcm_wvel  = dm.get<real,2>("gcm_wvel"       );
    // auto gcm_temp  = dm.get<real,2>("gcm_temp"       );
    // auto gcm_rho_v = dm.get<real,2>("gcm_water_vapor");
    //
    // parallel_for( Bounds<2>(crm_nz,nens) , YAKL_LAMBDA (int k, int iens) {
    //   gcm_rho_d(k,iens) = rho_d_col(k);
    //   gcm_uvel (k,iens) = uvel_col (k);
    //   gcm_vvel (k,iens) = vvel_col (k);
    //   gcm_wvel (k,iens) = wvel_col (k);
    //   gcm_temp (k,iens) = temp_col (k);
    //   gcm_rho_v(k,iens) = rho_v_col(k);
    // });
///

    yakl::timer_start("dycore");
    // this sets up initial conditons and reference state
    dycore.pre_time_loop(coupler);
    yakl::timer_stop("dycore");

    if (mainproc) {
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name () << std::endl;
      std::cout << "SGS   : " << sgs   .sgs_name   () << std::endl;
      std::cout << "\n";
      std::cout << "crm_nx:   " << crm_nx << "\n";
      std::cout << "crm_ny:   " << crm_ny << "\n";
      std::cout << "crm_nz:   " << crm_nz << "\n";
      std::cout << "xlen (m): " << xlen << "\n";
      std::cout << "ylen (m): " << ylen << "\n";
      std::cout << "Vertical interface heights: ";
      auto zint_host = zint_in.createHostCopy();
      for (int k=0; k < crm_nz+1; k++) {
        std::cout << zint_host(k) << "  ";
      }
      std::cout << "\n\n";
    }

//FIX THIS UP A LITTLE
    // Initialize the CRM internal state from the initial GCM column and random temperature perturbations
  //  modules::broadcast_initial_gcm_column( coupler );
///

    // Now that we have an initial state, define hydrostasis for each ensemble member
    if (use_coupler_hydrostasis) coupler.update_hydrostasis( );

//FIX THIS UP A LITTLE
  //  int1d seeds("seeds",nens);
  //  seeds = 0;
  //  modules::perturb_temperature( coupler , seeds );
///

    real etime_gcm = 0;
    int  num_out = 0;

    // Output the initial state
    yakl::timer_start("output");
    if (out_freq >= 0. ) output( coupler , out_prefix , etime_gcm );
    yakl::timer_stop("output");

    yakl::fence();
    auto ts = std::chrono::steady_clock::now();
    yakl::timer_start("main_loop");
    for (int step_gcm = 0; step_gcm < nsteps_gcm; ++step_gcm) {

//EVENTUALLY WE NEED TO CREATE/DESTROY THINGS FOR EACH GCM TIME STEP!
//ACTUALLY CREATE TWO DIFFERENT DRIVERS HERE....

//    modules::compute_gcm_forcing_tendencies( coupler );

      for (int step_crm_phys = 0; step_crm_phys < nsteps_crm_phys; ++step_crm_phys) {

//SAME ORDER AS MMF? IT SHOULD BE!
        //yakl::timer_start("apply_gcm_forcing_tendencies");
        //coupler.run_module( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies                        );
        //yakl::timer_stop("apply_gcm_forcing_tendencies");

        yakl::timer_start("dycore");
        coupler.run_module( "dycore"                       , [&] (pam::PamCoupler &coupler) { dycore.timeStep(coupler); } );
        yakl::timer_stop("dycore");

        //yakl::timer_start("sponge_layer");
        //coupler.run_module( "sponge_layer"                 , modules::sponge_layer                                        );
        //yakl::timer_stop("sponge_layer");

        yakl::timer_start("sgs");
        coupler.run_module( "sgs"                          , [&] (pam::PamCoupler &coupler) { sgs   .timeStep(coupler); } );
        yakl::timer_stop("sgs");

        yakl::timer_start("micro");
        coupler.run_module( "micro"                        , [&] (pam::PamCoupler &coupler) { micro .timeStep(coupler); } );
        yakl::timer_stop("micro");

        etime_gcm = step_gcm * dt_gcm + (step_crm_phys + 1) * dt_crm_phys;

        if (out_freq >= 0. && etime_gcm / out_freq >= num_out+1) {
          yakl::timer_start("output");
          output( coupler , out_prefix , etime_gcm);
          yakl::timer_stop("output");

          auto &dm = coupler.get_data_manager_device_readonly();
          real maxw = maxval(abs(dm.get_collapsed<real const>("wvel")));
          if (mainproc) {
            std::cout << "Etime , dtphys, maxw: " << etime_gcm << " , "
                                                  << dt_crm_phys    << " , "
                                                  << std::setw(10) << maxw << std::endl;
          }
          num_out++;
        }
      }
    }

    yakl::timer_stop("main_loop");
    yakl::fence();
    auto te = std::chrono::steady_clock::now();

    auto runtime = std::chrono::duration<double>(te - ts).count();

    if (mainproc) {
      std::cout << "Simulation Time: " << etime_gcm << "\n";
      std::cout << "Run Time: " << runtime << "\n";
    }

    yakl::timer_start("micro");
    micro .finalize(coupler);
    yakl::timer_stop("micro");

    yakl::timer_start("sgs");
    sgs   .finalize(coupler);
    yakl::timer_stop("sgs");

    yakl::timer_start("dycore");
    dycore.finalize( coupler );
    yakl::timer_stop("dycore");

    pam_interface::finalize();

    yakl::timer_stop("main");
  }

  #if defined(P3_CXX) || defined(SHOC_CXX)
    pam::deallocate_scream_cxx_globals();
    pam::call_kokkos_finalize();
  #endif
  yakl::finalize();
  MPI_Finalize();
}
