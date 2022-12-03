
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "mpi.h"
#include "output.h"
#include "pamc_init.h"

int main(int argc, char** argv) {
  int ierr = MPI_Init( &argc , &argv );
  yakl::init();

  {

    int myrank;
    bool masterproc;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    // Determine if I'm the master process
    masterproc = myrank == 0;

    using yakl::intrinsics::abs;
    using yakl::intrinsics::maxval;
    yakl::timer_start("main");

    pam::PamCoupler coupler;

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    real        simTime = config["simTime" ].as<real>();
    real        gcm_physics_dt_in = config["gcm_physics_dt" ].as<real>();
    int         crm_nx       = config["crm_nx"  ].as<int>();
    int         crm_ny       = config["crm_ny"  ].as<int>();
    int         nens         = config["nens"    ].as<int>();
    real        xlen         = config["xlen"    ].as<real>(-1.0_fp);
    real        ylen         = config["ylen"    ].as<real>(-1.0_fp);
    real        crm_dt_in    = config["crm_dt"  ].as<real>();
    std::string vcoords_file = config["vcoords" ].as<std::string>();

//WHAT DOES THIS DO? IS IT EVEN NEEDED?
    bool        use_coupler_hydrostasis = config["use_coupler_hydrostasis"].as<bool>(false);
    std::string out_freq                = config["out_freq"   ].as<real>();
    std::string out_prefix              = config["out_prefix" ].as<std::string>("test");
    bool        inner_mpi = config["inner_mpi"].as<bool>(false);

//ELMINATE?
    real zlen;

//BROKEN FOR ENSEMBLES!
//do we even really need zlen?
//used ONLY for setting uniform vertical grids, which we will fix up anyways...
    if (xlen < 0 || ylen < 0) { set_domain_sizes(config, xlen, ylen, zlen); }

//FIX THIS UP A LITTLE I THINK?
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

    //// Create the dycore and the microphysics
    Dycore       dycore;
    Microphysics micro;
    SGS          sgs;

    //this partitions the domain if INNER_MPI is set, otherwise it does nothing
    if (inner_mpi)
    {partition_domain(inFile, crm_nx, crm_ny);}

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
    yakl::timer_start("dycore");
    dycore.init( coupler ); // Dycore should initialize its own state here
    yakl::timer_stop("dycore");

    if (masterproc) {
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name () << std::endl;
      std::cout << "SGS   : " << sgs   .sgs_name   () << std::endl;
      std::cout << "crm_dt_in:         " << crm_dt_in << "\n";
      std::cout << "gcm_physics_dt:         " << gcm_physics_dt << "\n";
      std::cout << "\n";
    }

    // Now that we have an initial state, define hydrostasis for each ensemble member
    if (use_coupler_hydrostasis) coupler.update_hydrostasis( );

    real etime = 0;

    int  num_out = 0;
    // Output the initial state
    yakl::timer_start("output");
    if (out_freq >= 0. ) output( coupler , out_prefix , etime );
    yakl::timer_stop("output");

    real gcm_dt = gcm_physics_dt_in;
    real crm_dt = crm_dt_in;
    while (etime < simTime) {
      if (etime + gcm_dt > simTime) { gcm_dt = simTime - gcm_dt; }

      //      modules::compute_gcm_forcing_tendencies( coupler , dt_gcm );

//ALL BROKEN TIME STEPPING CURRENTLY
      real etime_crm = 0;
      real simTime_crm = dt_gcm;
      real dt_crm = dt_crm_phys;

      while (etime_crm < simTime_crm) {
        if (dt_crm == 0.) { dt_crm = dycore.compute_time_step(coupler); }
        if (etime_crm + dt_crm > simTime_crm) { dt_crm = simTime_crm - etime_crm; }

      yakl::timer_start("sgs");
      sgs.timeStep( coupler , crm_dt );
      yakl::timer_stop("sgs");

      yakl::timer_start("micro");
      micro.timeStep( coupler , crm_dt );
      yakl::timer_stop("micro");

      yakl::timer_start("dycore");
      dycore.timeStep( coupler , crm_dt );
      yakl::timer_stop("dycore");

      etime += crm_dt;

      auto &dm = coupler.get_data_manager_readonly();
      real maxw = maxval(abs(dm.get_collapsed<real const>("wvel")));
      if (masterproc) {
      std::cout << "Etime , crm_dt, maxw: " << etime  << " , "
                                            << crm_dt << " , "
                                            << std::setw(10) << maxw << "\n";
      }
        if (out_freq >= 0. && etime / out_freq >= num_out+1) {
          yakl::timer_start("output");
          output( coupler , out_prefix , etime );
          yakl::timer_stop("output");
          num_out++;
        }

    }

    if (masterproc) {std::cout << "Elapsed Time: " << etime << "\n";}

    yakl::timer_start("dycore");
    dycore.finalize( coupler );
    yakl::timer_stop("dycore");

    yakl::timer_stop("main");
  }

  yakl::finalize();

  ierr = MPI_Finalize();

}
