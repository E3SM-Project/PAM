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
#include "KGRK.h"
#include "LSRK.h"
#include "SI.h"
#include "SSPRK.h"
#include "time_integrator.h"
#include <memory>

using pam::PamCoupler;

class Dycore {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

  Equations equations;
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
  Parallel par;
  ModelLinearSystem linear_system;
  std::unique_ptr<TimeIntegrator> time_integrator;

  int ierr;
  real etime = 0.0;
  uint prevstep = 0;

  std::unique_ptr<TimeIntegrator>
  choose_time_integrator(const std::string &tstype) {
    if (tstype.substr(0, 5) == "ssprk") {
      return std::make_unique<SSPRKTimeIntegrator>(tstype);
    } else if (tstype.substr(0, 4) == "lsrk") {
      return std::make_unique<LSRKTimeIntegrator>(tstype);
    } else if (tstype.substr(0, 2) == "si") {
      return std::make_unique<SITimeIntegrator<si_quad_pts>>(tstype);
    } else if (tstype.substr(0, 4) == "kgrk") {
      return std::make_unique<KGRKTimeIntegrator>(tstype);
    } else {
      throw std::runtime_error("unknown time integrator type");
    }
  }

  void init(PamCoupler &coupler) {
    serial_print("setting up dycore", par.masterproc);

    // Set parameters
    debug_print(
        "reading parameters and partitioning domain/setting domain sizes",
        par.masterproc);

    if (coupler.option_exists("standalone_input_file")) {
      std::string inFile =
          coupler.get_option<std::string>("standalone_input_file");
      read_model_params_file(inFile, params, par, coupler, testcase);
    } else {
      read_model_params_coupler(params, par, coupler, testcase);
    }

    finalize_parallel(params, par);

    testcase->set_tracers(params);
    testcase->set_domain(params);

    check_and_print_model_parameters(params, par);

    debug_print("read parameters and partitioned domain/setting domain sizes",
                par.masterproc);

    time_integrator = choose_time_integrator(params.tstype);

    // Initialize the grid
    debug_print("start init topo/geom", par.masterproc);
    primal_topology.initialize(par, true);
    primal_geometry.initialize(primal_topology, params);
    dual_topology.initialize(par, false);
    dual_geometry.initialize(dual_topology, params);
    debug_print("finish init topo/geom", par.masterproc);

    debug_print("start init equations", par.masterproc);
    equations.initialize(coupler, params, primal_geometry, dual_geometry);
    debug_print("finish init equations", par.masterproc);

    debug_print("start init testcase", par.masterproc);
    testcase->initialize(equations);
    debug_print("finish init testcase", par.masterproc);

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
      diag->initialize(primal_geometry, dual_geometry, equations);
    }
    debug_print("end diagnostics init", par.masterproc);

    // Initialize the statistics
    debug_print("start stats init", par.masterproc);
    stats.initialize(params, par, primal_geometry, dual_geometry, equations);
    debug_print("end stats init", par.masterproc);

    // Create the outputter
    debug_print("start io init", par.masterproc);
    io.initialize(params.outputName, primal_topology, dual_topology, par,
                  prognostic_vars, constant_vars, diagnostics, stats);
    debug_print("finish io init", par.masterproc);
  };

  void pre_time_loop(PamCoupler &coupler) {
    // // Set the reference state and initialize the tendencies
    debug_print("start tendencies init", par.masterproc);
    // the reference state has to be set before the tendencies are initialized
    // for anelastic
    testcase->set_reference_state(primal_geometry, dual_geometry);
    tendencies.initialize(params, equations, primal_geometry, dual_geometry);
    debug_print("end tendencies init", par.masterproc);

    // EVENTUALLY THIS NEEDS TO BE MORE CLEVER IE POSSIBLY DO NOTHING BASED ON
    // THE IC STRING?
    //  set the initial conditions and compute initial stats
    debug_print("start ic setting", par.masterproc);
    testcase->set_initial_conditions(prognostic_vars, constant_vars,
                                     primal_geometry, dual_geometry);
    prognostic_vars.exchange();
    constant_vars.exchange();

    tendencies.compute_constants(constant_vars, prognostic_vars);
    constant_vars.exchange();
    debug_print("end ic setting", par.masterproc);

    // // Initialize the time integrator
    debug_print("start ts init", par.masterproc);
    time_integrator->initialize(params, tendencies, linear_system,
                                prognostic_vars, constant_vars, auxiliary_vars);
    if (time_integrator->is_semi_implicit) {
      linear_system.initialize(params, primal_geometry, dual_geometry,
                               equations);
      linear_system.compute_coefficients(params.dtcrm);
    }
    debug_print("end ts init", par.masterproc);

    // convert dynamics state to Coupler state
    tendencies.convert_dynamics_to_coupler_state(coupler, prognostic_vars,
                                                 constant_vars);

    // Output the initial model state
#ifndef _NOIO
    debug_print("start initial io", par.masterproc);
    for (auto &diag : diagnostics) {
      diag->compute(0, constant_vars, prognostic_vars);
    }
    stats.compute(prognostic_vars, constant_vars, 0);
    io.outputInit(etime);
    io.outputStats(stats);
    debug_print("end initial io", par.masterproc);
#endif
    prevstep = 1;
  }

  // Given the model data and CFL value, compute the maximum stable time step
  real compute_time_step(PamCoupler const &coupler, real cfl_in = -1) {
    return 0._fp;
  };

  void timeStep(PamCoupler &coupler) {

    serial_print("taking a dycore dtphys step", par.masterproc);

    // convert Coupler state to dynamics state
    tendencies.convert_coupler_to_dynamics_state(coupler, prognostic_vars,
                                                 auxiliary_vars, constant_vars);

    // Time stepping loop
    debug_print("start time stepping loop", par.masterproc);
    for (uint nstep = 0; nstep < params.crm_per_phys; nstep++) {
      yakl::fence();
      time_integrator->step_forward(params.dtcrm);

#ifdef CHECK_ANELASTIC_CONSTRAINT
      //real max_div = tendencies.compute_max_anelastic_constraint(
      //    prognostic_vars, auxiliary_vars);
      //std::cout << "Anelastic constraint: " << max_div << std::endl;
#endif//

      yakl::fence();

      etime += params.dtcrm;
#ifndef _NOIO
      if ((nstep + prevstep) % params.Nout == 0) {
        serial_print("dycore step " + std::to_string((nstep + prevstep)) +
                         " time " + std::to_string(etime),
                     par.masterproc);
        for (auto &diag : diagnostics) {
          diag->compute(etime, constant_vars, prognostic_vars);
        }
        io.output(etime);
        io.outputStats(stats);
      }

      if ((nstep + prevstep) % params.Nstat == 0) {
        stats.compute(prognostic_vars, constant_vars,
                      (nstep + prevstep) / params.Nstat);
      }
#endif
    }
    prevstep += params.crm_per_phys;

    if (!time_integrator->is_ssp) {
      tendencies.remove_negative_densities(prognostic_vars);
    }

    // convert dynamics state to Coupler state
    tendencies.convert_dynamics_to_coupler_state(coupler, prognostic_vars,
                                                 constant_vars);

    // ADD COUPLER DYCORE FUNCTIONS HERE AS WELL

    debug_print("end time stepping loop", par.masterproc);
  };

  void finalize(PamCoupler &coupler) { io.outputStats(stats); }

  const char *dycore_name() const { return "SPAM++"; }
};
