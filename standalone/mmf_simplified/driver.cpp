
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "mmf_interface.h"
#include "perturb_temperature.h"
#include "gcm_forcing.h"
#include "sponge_layer.h"
#include "saturation_adjustment.h"
#include "broadcast_initial_gcm_column.h"
#include "output.h"
#include "supercell_init.h"
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
//ADD VERTICAL GRID STUFF HERE
    auto simTime = config["simTime" ].as<real>();
    auto gcm_physics_dt = config["gcm_physics_dt" ].as<real>();
    auto crm_nx         = config["crm_nx"     ].as<int>();
    auto crm_ny         = config["crm_ny"     ].as<int>();
    auto crm_nz         = config["crm_nz"  ].as<int>(0);
    auto nens           = config["nens"       ].as<int>();
    auto xlen           = config["xlen"       ].as<real>(-1.0_fp);
    auto ylen           = config["ylen"       ].as<real>(-1.0_fp);
    auto crm_dt_in      = config["crm_dt"  ].as<real>();
    auto out_freq       = config["out_freq"   ].as<real>();
    auto out_prefix     = config["out_prefix" ].as<std::string>();
    real zlen;
    
    int nranks;
    int myrank;
    MPI_Comm_size( MPI_COMM_WORLD , &nranks );
    MPI_Comm_rank( MPI_COMM_WORLD , &myrank );
    bool mainproc = (myrank == 0);

    auto &coupler = mmf_interface::get_coupler();

    //set xlen, ylen based on init cond if needed
    if (xlen < 0 || ylen < 0) { set_domain_sizes(config, xlen, ylen, zlen); }

//FIX UP A LITTLE
    int crm_nz = -1;
    real1d zint_in;
    if (vcoords_file == "uniform") {
      if (config["crm_nz"]) {
        crm_nz = config["crm_nz"].as<int>();
      } else {
        endrun("To use uniform vertical grid you need to specify crm_nz");
      }
      zint_in = real1d("zint_in",crm_nz+1);
      const real dz = zlen / (crm_nz - 1);
      parallel_for("uniform zint", crm_nz+1, YAKL_LAMBDA(int k) {
          if (k == 0) {
            zint_in(k) = 0;
          } else if (k == crm_nz) {
            zint_in(k) = zlen;
          } else {
            zint_in(k) = k * dz - dz / 2;
          }
      });
    } else {
      // Read vertical coordinates
      //THIS IS BROKEN FOR PARALLEL IO CASE- maybe this is okay ie switch entirely to standard netcdf?
      yakl::SimpleNetCDF nc;
      nc.open(vcoords_file);
      crm_nz = nc.getDimSize("num_interfaces") - 1;
      zint_in = real1d("zint_in",crm_nz+1);
      nc.read(zint_in,"vertical_interfaces");
      nc.close();
      //TODO: Coupler needs to eventually support ensemble dependent vertical grids
      //real2d zint_expanded = real2d("zint_expanded",crm_nz+1,nens);
      //parallel_for( "Set zint expanded" , SimpleBounds<2>(crm_nz+1,nens) , YAKL_LAMBDA (int k, int n) {
      //  zint_expanded(k,n) = zint_in(k);
      //});
    }

    // Allocates the coupler state (density_dry, uvel, vvel, wvel, temp, vert grid, hydro background) for thread 0
    coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

    if !(config["initData" ].as<std::string>() == 'coupler')
    {endrun("mmf_simplified only supports setting the initial condition through the coupler\n");}

    real1d rho_d_col("rho_d_col",crm_nz);
    real1d uvel_col ("uvel_col" ,crm_nz);
    real1d vvel_col ("vvel_col" ,crm_nz);
    real1d wvel_col ("wvel_col" ,crm_nz);
    real1d temp_col ("temp_col" ,crm_nz);
    real1d rho_v_col("rho_v_col",crm_nz);

    auto initCouplerData = config["initCouplerData" ].as<std::string>();

    // Compute a supercell initial column
    if (initCouplerData == "supercell") {
    supercell_init( zint_in, rho_d_col, uvel_col, vvel_col, wvel_col, temp_col, rho_v_col, coupler.get_R_d() ,
                    coupler.get_R_v(), coupler.get_grav() );
    }

    // Set the GCM column data for each ensemble to the initial state
    auto &dm = coupler.get_data_manager_readwrite();

    auto gcm_rho_d = dm.get<real,2>("gcm_density_dry");
    auto gcm_uvel  = dm.get<real,2>("gcm_uvel"       );
    auto gcm_vvel  = dm.get<real,2>("gcm_vvel"       );
    auto gcm_wvel  = dm.get<real,2>("gcm_wvel"       );
    auto gcm_temp  = dm.get<real,2>("gcm_temp"       );
    auto gcm_rho_v = dm.get<real,2>("gcm_water_vapor");

    parallel_for( Bounds<2>(crm_nz,nens) , YAKL_LAMBDA (int k, int iens) {
      gcm_rho_d(k,iens) = rho_d_col(k);
      gcm_uvel (k,iens) = uvel_col (k);
      gcm_vvel (k,iens) = vvel_col (k);
      gcm_wvel (k,iens) = wvel_col (k);
      gcm_temp (k,iens) = temp_col (k);
      gcm_rho_v(k,iens) = rho_v_col(k);
    });

    // Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;
    SGS          sgs;

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

    // Set physical constants for coupler at thread 0 using microphysics data
    coupler.set_phys_constants( micro.R_d , micro.R_v , micro.cp_d , micro.cp_v , micro.grav , micro.p0 );

    // Set the vertical grid in the coupler
    coupler.set_grid( xlen , ylen , zint_in );

    micro .init( coupler );
    sgs   .init( coupler );
    dycore.init( coupler );

    // Initialize the CRM internal state from the initial GCM column and random temperature perturbations
    modules::broadcast_initial_gcm_column( coupler );

    // Now that we have an initial state, define hydrostasis for each ensemble member
    coupler.update_hydrostasis();

    modules::perturb_temperature( coupler , 0 );

    coupler.add_mmf_function( "apply_gcm_forcing_tendencies" , modules::apply_gcm_forcing_tendencies );
    coupler.add_mmf_function( "dycore" , [&] (pam::PamCoupler &coupler, real dt) { dycore.timeStep(coupler,dt); } );
    coupler.add_mmf_function( "sgs"    , [&] (pam::PamCoupler &coupler, real dt) { sgs   .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "micro"  , [&] (pam::PamCoupler &coupler, real dt) { micro .timeStep(coupler,dt); } );
    coupler.add_mmf_function( "sponge_layer"                 , modules::sponge_layer                 );

    real etime_gcm = 0;
    int  num_out = 0;

    // Output the initial state
    if (out_freq >= 0. ) output( coupler , out_prefix , etime_gcm );

    yakl::timer_start("main_loop");

//REMOVE ONE LAYER OF LOOPING HERE
//ACTUALLY, NO, KEEP IT!
//maybe rework this a little though?
//etime_crm can probably be dropped? no, it is needed to clip time step at the end of a gcm step if needed

    while (etime_gcm < simTime) {
      if (etime_gcm + dt_gcm > simTime) { dt_gcm = simTime - etime_gcm; }

      modules::compute_gcm_forcing_tendencies( coupler , dt_gcm );

      real etime_crm = 0;
      real simTime_crm = dt_gcm;
      real dt_crm = dt_crm_phys;
      while (etime_crm < simTime_crm) {
        if (dt_crm == 0.) { dt_crm = dycore.compute_time_step(coupler); }
        if (etime_crm + dt_crm > simTime_crm) { dt_crm = simTime_crm - etime_crm; }

        coupler.run_mmf_function( "apply_gcm_forcing_tendencies" , dt_crm );
        coupler.run_mmf_function( "dycore"                       , dt_crm );
        coupler.run_mmf_function( "sponge_layer"                 , dt_crm );
        coupler.run_mmf_function( "sgs"                          , dt_crm );
        coupler.run_mmf_function( "micro"                        , dt_crm );

        etime_crm += dt_crm;
        etime_gcm += dt_crm;
        if (out_freq >= 0. && etime_gcm / out_freq >= num_out+1) {
          yakl::timer_start("output");
          output( coupler , out_prefix , etime_gcm );
          yakl::timer_stop("output");
          real maxw = maxval(abs(dm.get_collapsed<real const>("wvel")));
          if (mainproc) {
            std::cout << "Etime , dtphys, maxw: " << etime_gcm << " , "
                                                  << dt_crm    << " , "
                                                  << std::setw(10) << maxw << std::endl;
          }
          num_out++;
        }
      }
    }
    yakl::timer_stop("main_loop");

    if (mainproc) {
      std::cout << "Elapsed Time: " << etime_gcm << "\n";
    }

    dycore.finalize( coupler );

    mmf_interface::finalize();

    yakl::timer_stop("main");
  }
  yakl::finalize();
  MPI_Finalize();
}
