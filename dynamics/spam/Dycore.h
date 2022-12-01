#pragma once

#include "DataManager.h"
#include "MultipleFields.h"
#include "common.h"
#include "field_sets.h"
#include "fileio.h"
#include "geometry.h"
#include "pam_coupler.h" //Has DataManager and pam_const
#include "parallel.h"
#include "params.h"
#include "topology.h"
#if defined _HAMILTONIAN && defined _LAYER
#include "layermodel.h"
#elif defined _HAMILTONIAN && defined _EXTRUDED
#include "extrudedmodel.h"
#elif defined _ADVECTION && defined _LAYER
#include "layeradvection.h"
#elif defined _ADVECTION && defined _EXTRUDED
#include "extrudedadvection.h"
#endif
#include "RKSimple.h"
#include "SI.h"
#include "SSPRK.h"
#include <memory>

using pam::PamCoupler;

class Dycore {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

  std::unique_ptr<TestCase> testcase;
  ModelStats stats;
  FieldSet<nprognostic> prognostic_vars;
  FieldSet<nconstant> constant_vars;
  FieldSet<nauxiliary> auxiliary_vars;
  ExchangeSet<nprognostic> prog_exchange;
  ExchangeSet<nconstant> const_exchange;
  ExchangeSet<nauxiliary> aux_exchange;
  FileIO io;
  ModelTendencies tendencies;
  std::vector<std::unique_ptr<Diagnostic>> diagnostics;
  Topology primal_topology;
  Topology dual_topology;
  ModelParameters params;
  ModelReferenceState reference_state;
  Parallel par;
#if _TIME_TYPE == 0
  RKSimpleTimeIntegrator tint;
#endif
#if _TIME_TYPE == 1
  SSPKKTimeIntegrator tint;
#endif
#if _TIME_TYPE == 2
  ModelLinearSystem linear_system;
  SITimeIntegrator<si_quad_pts> tint;
#endif

  int ierr;
  real etime = 0.0;
  real dt;
  int niters;

  void init(PamCoupler &coupler) {

    serial_print("setting up dycore", par.masterproc);

    // Set parameters

    debug_print(
        "reading parameters and partitioning domain/setting domain sizes",
        par.masterproc);
    std::string inFile =
        coupler.get_option<std::string>("standalone_input_file");
    readModelParamsFile(inFile, params, par, coupler, testcase);
    debug_print("read parameters and partitioned domain/setting domain sizes",
                par.masterproc);
    set_timestepping_variables(inFile):
    
    // Initialize the grid
    debug_print("start init topo/geom", par.masterproc);
    primal_topology.initialize(par, true);
    primal_geometry.initialize(primal_topology, params);
    dual_topology.initialize(par, false);
    dual_geometry.initialize(dual_topology, params);
    debug_print("finish init topo/geom", par.masterproc);

    debug_print("start init reference state", par.masterproc);
    reference_state.initialize(primal_topology, dual_topology);
    debug_print("finish init reference state", par.masterproc);

    // Allocate the variables
    debug_print("start init field/exchange sets", par.masterproc);
    // this gives basedof, extdof and ndofs
    std::array<FieldDescription, nprognostic> prog_desc_arr;
    std::array<FieldDescription, nconstant> const_desc_arr;
    std::array<FieldDescription, nauxiliary> aux_desc_arr;
    initialize_variables(primal_topology, dual_topology, prog_desc_arr,
                         const_desc_arr, aux_desc_arr);

    prog_exchange.initialize(prog_desc_arr);
    const_exchange.initialize(const_desc_arr);
    aux_exchange.initialize(aux_desc_arr);

    prognostic_vars.initialize("x", prog_desc_arr, prog_exchange);
    constant_vars.initialize("cons", const_desc_arr, const_exchange);
    auxiliary_vars.initialize("aux", aux_desc_arr, aux_exchange);
    debug_print("finish init field/exchange sets", par.masterproc);

    debug_print("start diagnostics init", par.masterproc);
    add_model_diagnostics(diagnostics);
    testcase->add_diagnostics(diagnostics);
    for (auto &diag : diagnostics) {
      diag->initialize(primal_geometry, dual_geometry);
    }
    debug_print("end diagnostics init", par.masterproc);

    // Initialize the statistics
    debug_print("start stats init", par.masterproc);
    stats.initialize(params, par, primal_geometry, dual_geometry);
    debug_print("end stats init", par.masterproc);

    // Create the outputter
    debug_print("start io init", par.masterproc);
    io.initialize(params.outputName, primal_topology, dual_topology, par,
                  prognostic_vars, constant_vars, diagnostics, stats);
    debug_print("finish io init", par.masterproc);

    // // Initialize the tendencies and diagnostics
    debug_print("start tendencies init", par.masterproc);
    tendencies.initialize(coupler, params, primal_geometry, dual_geometry,
                          reference_state);
    debug_print("end tendencies init", par.masterproc);

    //  set the initial conditions and compute initial stats
    debug_print("start ic setting", par.masterproc);
//FIX THESE FOR COUPLERDATA
    testcase->set_reference_state(reference_state, primal_geometry,
                                  dual_geometry);
    testcase->set_initial_conditions(prognostic_vars, constant_vars,
                                     primal_geometry, dual_geometry);
    prognostic_vars.exchange();
    constant_vars.exchange();

    tendencies.compute_constants(constant_vars, prognostic_vars);
    constant_vars.exchange();
    stats.compute(prognostic_vars, constant_vars, 0);
    debug_print("end ic setting", par.masterproc);

    // // Initialize the time stepper
    debug_print("start ts init", par.masterproc);
#if _TIME_TYPE == 2
    linear_system.initialize(params, primal_geometry, dual_geometry,
                             reference_state);
    linear_system.compute_coefficients(params.dtcrm);
    tint.initialize(params, tendencies, linear_system, prognostic_vars,
                    constant_vars, auxiliary_vars);
#else
    tint.initialize(params, tendencies, prognostic_vars, constant_vars,
                    auxiliary_vars);
#endif
    debug_print("end ts init", par.masterproc);

    // convert dynamics state to Coupler state
    tendencies.convert_dynamics_to_coupler_state(coupler, prognostic_vars,
                                                 constant_vars);

    // Output the initial model state
    debug_print("start initial io", par.masterproc);
    for (auto &diag : diagnostics) {
      diag->compute(0, reference_state, constant_vars, prognostic_vars);
    }
    stats.compute(prognostic_vars, constant_vars, 0);
    io.outputInit(etime);
    io.outputStats(stats);
    debug_print("end initial io", par.masterproc);

    prevstep = 1;
  };

  // Given the model data and CFL value, compute the maximum stable time step
  real compute_time_step(PamCoupler const &coupler, real cfl_in = -1) {
    
//GET CFL FROM TIMESTEPPER!
    real cfl = cfl_in;
    if (cfl < 0) cfl = 0.75;

    // If we've already computed the time step, then don't compute it again
    if (dtInit <= 0) {

//GET FROM THERMO/PARAMS/GEOM AS APPROPRIATE
      YAKL_SCOPE( dx                   , this->dx                  );
      YAKL_SCOPE( dy                   , this->dy                  );
      YAKL_SCOPE( dz                   , this->dz                  );
      YAKL_SCOPE( gamma                , this->gamma               );
      YAKL_SCOPE( Rd                   , this->Rd                  );
      YAKL_SCOPE( Rv                   , this->Rv                  );

      auto &dm = coupler.get_data_manager_readonly();

      // Convert data from DataManager to state and tracers array for convenience
      auto dm_dens_dry = dm.get<real const,4>( "density_dry" );
      auto dm_uvel     = dm.get<real const,4>( "uvel"        );
      auto dm_vvel     = dm.get<real const,4>( "vvel"        );
      auto dm_wvel     = dm.get<real const,4>( "wvel"        );
      auto dm_temp     = dm.get<real const,4>( "temp"        );
      auto dm_dens_vap = dm.get<real const,4>( "water_vapor" );

      // Allocate a 3-D array for the max stable time steps (we'll use this for a reduction later)
      real4d dt3d("dt3d",nz,ny,nx,nens);

      // Loop through the cells, calculate the max stable time step for each cell
      parallel_for( "Spatial.h compute_time_step" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        real rho_d = dm_dens_dry(k,j,i,iens);
        real u     = dm_uvel    (k,j,i,iens);
        real v     = dm_vvel    (k,j,i,iens);
        real w     = dm_wvel    (k,j,i,iens);
        real temp  = dm_temp    (k,j,i,iens);
        real rho_v = dm_dens_vap(k,j,i,iens);
        real p = Rd * rho_d * temp + Rv * rho_v * temp;
        // This neglects non-wv mass-adding tracers, but these are small, and their lack only increases cs
        // Thus the resulting time step is conservative w/r to these missing masses, which is more stable
        real cs = sqrt(gamma*p/(rho_v+rho_d));

        // Compute the maximum stable time step in each direction
        real udt = cfl * dx         / max( abs(u-cs) , abs(u+cs) );
        real vdt = cfl * dy         / max( abs(v-cs) , abs(v+cs) );
        if (sim2d) vdt = std::numeric_limits<real>::max();
        real wdt = cfl * dz(k,iens) / max( abs(w-cs) , abs(w+cs) );

        // Compute the min of the max stable time steps
        dt3d(k,j,i,iens) = min( min(udt,vdt) , wdt );
      });
      // Store to dtInit so we don't have to compute this again
      dtInit = yakl::intrinsics::minval( dt3d );
    }

    return dtInit;
    
    return 0._fp;
  };

  void timeStep(PamCoupler &coupler, real dtphys) {

    // convert Coupler state to dynamics state
    tendencies.convert_coupler_to_dynamics_state(coupler, prognostic_vars,
                                                 constant_vars);

    // Time stepping loop
    debug_print("start time stepping loop", par.masterproc);
//SWAP THIS CAREFULLY
//Actually this is fine -> either set crm_per_phys based on 

real dt;
if (dycore_per_phys_step > 0)
{
dt = compute_time_step( coupler );
int n_iter = ceil( dtphys / dt );
dt = dtphys / n_iter;
}
else
{
  dt = dtphys / dycore_per_phys_step;
}

//PRINT OUT ACTUAL DT USED
serial_print("taking a crm_dt step", par.masterproc);

    for (uint nstep = 0; nstep < n_iter; nstep++) {
      yakl::fence();
      tint.stepForward(dt);
      yakl::fence();

      etime += dt;

//ADD NUM_OUT STUFF
      if (dycore_out_freq >= 0. && etime / dycore_out_freq >= num_out+1) {
        
        serial_print("dycore ouput at " + std::to_string(etime),par.masterproc);
        for (auto &diag : diagnostics) {
          diag->compute(etime, reference_state, constant_vars, prognostic_vars);
        }
        io.output(etime);
        io.outputStats(stats);
      }

//ADD NUM_STAT STUFF
      if (dycore_stat_freq >= 0. && etime / dycore_stat_freq >= num_stat+1) {
        stats.compute(prognostic_vars, constant_vars,
                      (nstep + prevstep) / params.Nstat);
        //HOW DO WE HANDLE THIS LAST ARGUMENT?
      }
    }

    if (remove_negative_densities) {
      tendencies.remove_negative_densities(prognostic_vars);
    }

    // convert dynamics state to Coupler state
    tendencies.convert_dynamics_to_coupler_state(coupler, prognostic_vars,
                                                 constant_vars);

    debug_print("end time stepping loop", par.masterproc);
  };

//CAN ACTUALLY MOSTLY MOVE THIS ALL BACK INTO THE PARAMS
//IE ALL THAT REALLY MATTERS IS OUTPUT/STAT FREQUENCY!!!
//AND dycore_per_phys, etc.
void set_timestepping()
{
  
  YAML::Node config = YAML::LoadFile(inFile);


  //two types of time stepping control: total time based and step based
  simSteps = config["simSteps"].as<int>(0);
  simTime = config["simTime"].as<real>(0.0_fp);


  std::cout << "dtcrm:         " << dtcrm << "\n";
  std::cout << "dtphys:         " << dtphys << "\n";
  std::cout << "Nsteps:     " << Nsteps << "\n";
  std::cout << "simSteps:     " << simSteps << "\n";
  std::cout << "crm per phys:     " << crm_per_phys << "\n";
  std::cout << "Nout:       " << Nout << "\n";


  dtphys = config["dtphys"].as<real>();
  crm_per_phys = config["crm_per_phys"].as<int>(0);
  Nout = config["outSteps"].as<int>(0);
  Nstat = config["statSteps"].as<int>(0);
  if (not(params.dtphys > 0.0_fp) or not(simSteps > 0) or
      not(params.crm_per_phys > 0) or not(params.Nout > 0) or
      not(params.Nstat > 0)) {
    endrun("spam++ must use step-based time control logic ie set simSteps >0, "
           "dtphys>0, crm_per_phys >0, outSteps >0, statSteps >0");
  }

//FIX THIS A LITTLE I THINK
  params.dtcrm = params.dtphys / params.crm_per_phys;
  params.Nsteps = simSteps * params.crm_per_phys;
}
  void finalize(PamCoupler &coupler) { io.outputStats(stats); }

  const char *dycore_name() const { return "SPAM++"; }
};
