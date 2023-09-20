
#include "pam_coupler.h"
#include "Dycore.h"
#include "Microphysics.h"
#include "SGS.h"
#include "mpi.h"
#include "output.h"
#include "pamc_init.h"
#include <chrono>

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

    real        simTime      = config["simTime" ].as<real>(0.0_fp);
    real        simSteps     = config["simSteps"].as<int>(0);
    int         crm_nx       = config["crm_nx"  ].as<int>();
    int         crm_ny       = config["crm_ny"  ].as<int>();
    int         nens         = config["nens"    ].as<int>();
    real        xlen         = config["xlen"    ].as<real>(-1.0_fp);
    real        ylen         = config["ylen"    ].as<real>(-1.0_fp);
    real        zlen         = config["zlen"    ].as<real>(-1.0_fp);
    real        dtphys_in    = config["dtphys"  ].as<real>();
    std::string vcoords_file = config["vcoords" ].as<std::string>();
    bool        use_coupler_hydrostasis = config["use_coupler_hydrostasis"].as<bool>(false);
    auto out_freq                = config["out_freq"   ].as<real>(0);
    auto out_prefix              = config["out_prefix" ].as<std::string>("test");
    bool        inner_mpi = config["inner_mpi"].as<bool>(false);

    //set xlen, ylen, zlen based on init cond if needed
    if (xlen < 0 || ylen < 0 || zlen < 0) { set_domain_sizes(config, xlen, ylen, zlen); }

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

      using yakl::c::parallel_for;
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

    // Allocate coupler state
    coupler.allocate_coupler_state( crm_nz , crm_ny , crm_nx , nens );

    // Set the horizontal domain lengths and the vertical grid in the coupler
    coupler.set_grid( xlen , ylen , zint_in );

    // This is for the dycore to pull out to determine how to do idealized test cases
    coupler.set_option<std::string>( "standalone_input_file" , inFile );
    
    coupler.set_option<real>("crm_dt", dtphys_in);

    micro .init( coupler );
    sgs   .init( coupler );
    yakl::timer_start("dycore");
    dycore.init( coupler ); // Dycore should initialize its own state here
    yakl::timer_stop("dycore");

    #ifdef PAM_STANDALONE
    if (masterproc) {
      std::cout << "Dycore: " << dycore.dycore_name() << std::endl;
      std::cout << "Micro : " << micro .micro_name () << std::endl;
      std::cout << "SGS   : " << sgs   .sgs_name   () << std::endl;
      std::cout << "\n";
    }
    #endif

    yakl::timer_start("dycore");
    // this sets up initial conditons and reference state 
    dycore.pre_time_loop(coupler);
    yakl::timer_stop("dycore");

    real etime = 0;
    // There are two ways of time control- setting total simulation time (simTime) or setting number of physics time steps (simSteps)
    if (simTime == 0.0) {  simTime = simSteps * dtphys_in; }

    int  num_out = 0;
    // Output the initial state
    yakl::timer_start("output");
    if (out_freq >= 0. ) output( coupler , out_prefix , etime );
    yakl::timer_stop("output");

    real dtphys = dtphys_in;

    yakl::fence();
    auto ts = std::chrono::steady_clock::now();
    while (etime < simTime) {
      yakl::timer_start("dycore");
      if (dtphys_in == 0.) { dtphys = dycore.compute_time_step(coupler); }
      yakl::timer_stop("dycore");
      if (etime + dtphys > simTime) { dtphys = simTime - etime; }

      yakl::timer_start("sgs");
      sgs.timeStep( coupler);
      yakl::timer_stop("sgs");

      yakl::timer_start("micro");
      micro.timeStep( coupler);
      yakl::timer_stop("micro");

      yakl::timer_start("dycore");
      dycore.timeStep( coupler);
      yakl::timer_stop("dycore");

      etime += dtphys;

      auto &dm = coupler.get_data_manager_device_readonly();
      real maxw = maxval(abs(dm.get_collapsed<real const>("wvel")));
      if (masterproc) {
      std::cout << "Etime , dtphys, maxw: " << etime  << " , "
                                            << dtphys << " , "
                                            << std::setw(10) << maxw << "\n";
      }
        if (out_freq >= 0. && etime / out_freq >= num_out+1) {
          yakl::timer_start("output");
          output( coupler , out_prefix , etime );
          yakl::timer_stop("output");
          num_out++;
        }

    }
    yakl::fence();
    auto te = std::chrono::steady_clock::now();
    auto runtime = std::chrono::duration<double>(te - ts).count();

    if (masterproc) {
      std::cout << "Simulation Time: " << etime << "\n";
      std::cout << "Run Time: " <<  runtime << "\n";
    }

    yakl::timer_start("dycore");
    dycore.finalize( coupler );
    yakl::timer_stop("dycore");

    yakl::timer_stop("main");
  }

  yakl::finalize();

  ierr = MPI_Finalize();

}
