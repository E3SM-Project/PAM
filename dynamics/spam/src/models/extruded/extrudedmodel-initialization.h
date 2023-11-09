#pragma once

namespace pamc {

std::unique_ptr<LinearSystem> model_linear_system() {
  if (VariableSet::compressible) {
    // return std::make_unique<CompressibleVelocityLinearSystem>();
    return std::make_unique<CompressiblePressureLinearSystem>();
  } else {
    return std::make_unique<AnelasticLinearSystem>();
  }
}

// *******   FieldSet Initialization   ***********//
void initialize_variables(
    const Topology &ptopo, const Topology &dtopo,
    std::array<FieldDescription, nprognostic> &prog_desc_arr,
    std::array<FieldDescription, nconstant> &const_desc_arr,
    std::array<FieldDescription, nauxiliary> &aux_desc_arr) {

  using VS = VariableSet;

  // primal grid represents straight quantities, dual grid twisted quantities
  // ndims is the BASEDIM size!

  // v, w, dens
  prog_desc_arr[VVAR] = {"v", ptopo, 1, 0, 1}; // v = straight (1,0)-form
  prog_desc_arr[WVAR] = {"w", ptopo, 0, 1, 1}; // w = straight (0,1)-form
  prog_desc_arr[DENSVAR] = {
      "dens", dtopo, ndims, 1,
      VS::ndensity_prognostic}; // dens = twisted (n,1)-form

  // hs
  const_desc_arr[HSVAR] = {"hs", dtopo, ndims, 1, 1}; // hs = twisted (n,1)-form
  const_desc_arr[CORIOLISHZVAR] = {"coriolishz", ptopo, 1, 1,
                                   1}; // f = straight (1,1)-form

  // functional derivatives = F, B, K, he, U
  aux_desc_arr[FVAR] = {"F", dtopo, ndims - 1, 1, 1}; // F = twisted
                                                      // (n-1,1)-form

  aux_desc_arr[BVAR] = {"B", ptopo, 0, 0,
                        VS::ndensity_active}; // B = straight (0,0)-form

  aux_desc_arr[KVAR] = {"K", dtopo, ndims, 1, 1}; // K = twisted (n,1)-form

  aux_desc_arr[UVAR] = {"U", dtopo, ndims - 1, 1, 1}; // U = twisted
                                                      // (n-1,1)-form

  aux_desc_arr[FWVAR] = {"Fw", dtopo, ndims, 0, 1}; // Fw = twisted (n,0)-form

  aux_desc_arr[UWVAR] = {"Uw", dtopo, ndims, 0, 1}; // Uw = twisted (n,0)-form

  // second copy on F/Fw, needed for the semi-implicit solver because
  // the symplecitic operator gets evaluated at a different state
  // and needs its own F and Fw
  aux_desc_arr[F2VAR] = {"F2", dtopo, ndims - 1, 1, 1}; // F2 = twisted
                                                        // (n-1,1)-form
  aux_desc_arr[FW2VAR] = {"Fw2", dtopo, ndims, 0,
                          1}; // Fw2 = twisted (n,0)-form

  aux_desc_arr[FTVAR] = {"FT", ptopo, 0, 1, ndims}; // FT
  aux_desc_arr[FTWVAR] = {"FTW", ptopo, 1, 0, 1};   // FTW

  // dens primal grid reconstruction stuff- dens0, edgerecon, recon
  aux_desc_arr[DENS0VAR] = {
      "dens0", ptopo, 0, 0,
      VS::ndensity_prognostic}; // dens0 = straight (0,0)-form
  aux_desc_arr[DENSEDGERECONVAR] = {
      "densedgerecon", dtopo, ndims, 1,
      2 * ndims * VS::ndensity_prognostic}; // densedgerecon lives on dual
                                            // cells, associated with F
  aux_desc_arr[DENSRECONVAR] = {"densrecon", dtopo, ndims - 1, 1,
                                VS::ndensity}; // densrecon lives on horiz dual
                                               // edges, associated with F
  aux_desc_arr[DENSVERTEDGERECONVAR] = {
      "densvertedgerecon", dtopo, ndims, 1,
      2 * VS::ndensity_prognostic}; // densedgerecon lives on dual cells,
                                    // associated with Fw
  aux_desc_arr[DENSVERTRECONVAR] = {
      "densvertrecon", dtopo, ndims, 0,
      VS::ndensity}; // densvertrecon lives on vert dual
                     // edges, associated with Fw

  // fct stuff- Phi, Mf, edgeflux
  aux_desc_arr[MFVAR] = {"Mf", dtopo, ndims, 1, VS::ndensity_prognostic};
  aux_desc_arr[EDGEFLUXVAR] = {"edgeflux", dtopo, ndims - 1, 1,
                               VS::ndensity_prognostic};
  aux_desc_arr[VERTEDGEFLUXVAR] = {"vertedgeflux", dtopo, ndims, 0,
                                   VS::ndensity_prognostic};

  // Q stuff
  aux_desc_arr[QHZVAR] = {"QHZ", dtopo, ndims - 1, 0,
                          1}; // in 2d QHZ=QXZ00 = twisted (0,0)-form,
                              // in 3d QHZ=(QXZ10,QYZ10) lives on twisted edges
  aux_desc_arr[QHZRECONVAR] = {
      "qhzrecon", ptopo, 0, 1,
      ndims}; // qhzrecon lives on vert primal edges, associated with w
  aux_desc_arr[QHZEDGERECONVAR] = {
      "qhzedgerecon", ptopo, 1, 1,
      2 * 1}; // qhzedgerecon lives on primal cells in 2d or primal faces in 3d,
              // associated with Fw/w
  aux_desc_arr[QHZVERTRECONVAR] = {
      "qhzvertrecon", ptopo, 1, 0,
      1}; // qhzsvertrecon lives on horiz primal edges, associated with v
  aux_desc_arr[QHZVERTEDGERECONVAR] = {
      "qhzvertedgerecon", ptopo, 1, 1,
      2 * 1}; // qhzvertedgerecon lives on primal cells 2d or primal faces in
              // 3d, associated with F/v
  aux_desc_arr[QHZFLUXVAR] = {
      "qhzflux", ptopo, 0, 1,
      1}; // qhzflux lives on vert primal edges, associated with w
  aux_desc_arr[QHZVERTFLUXVAR] = {
      "qhzvertflux", ptopo, 1, 0,
      1}; // qhzvertflux lives on horiz primal edges, associated with v

  aux_desc_arr[FHZVAR] = {"fhz", dtopo, ndims - 1, 0,
                          1}; // in 2d fhz=fxz00 is a twisted 0,0 form
                              // in 3d fhz=(fxz10,fyz10) lives on twisted edges
  aux_desc_arr[CORIOLISHZRECONVAR] = {
      "coriolishzrecon", ptopo, 0, 1,
      ndims}; // coriolishzrecon lives on vert primal edges, associated with w
  aux_desc_arr[CORIOLISHZEDGERECONVAR] = {
      "coriolishzedgerecon", ptopo, 1, 1,
      2 * 1}; // coriolishzedgerecon lives on primal cells in 2d or primal faces
              // in 3d, associated with Fw/w
  aux_desc_arr[CORIOLISHZVERTRECONVAR] = {
      "coriolishzvertrecon", ptopo, 1, 0,
      1}; // coriolishzsvertrecon lives on horiz primal edges, associated with v
  aux_desc_arr[CORIOLISHZVERTEDGERECONVAR] = {
      "coriolishzvertedgerecon", ptopo, 1, 1,
      2 * 1}; // coriolishzvertedgerecon lives on primal cells in 2d or primal
              // faces in 3d, associated with F/v

  // diffusion variables
  aux_desc_arr[FDIFFVAR] = {
      "F", dtopo, ndims - 1, 1,
      VS::ndensity_diffused}; // diffusive F = twisted (n-1,1)-form
  aux_desc_arr[FWDIFFVAR] = {
      "Fw", dtopo, ndims, 0,
      VS::ndensity_diffused}; // diffusive Fw = twisted (n,0)-form

  // additonal Q and coriolis variables for 3d
  if (ndims > 1) {
    aux_desc_arr[QXYVAR] = {"QXY", dtopo, 0, ndims - 1,
                            1}; // QXY lives on twisted edges
    aux_desc_arr[QXYRECONVAR] = {"qxyrecon", ptopo, 1, 0,
                                 2}; // qxyrecon lives on primal edges
    aux_desc_arr[QXYEDGERECONVAR] = {
        "qxyedgerecon", ptopo, 2, 0,
        4}; // qxyedgerecon lives on primal vertical faces

    aux_desc_arr[FXYVAR] = {"FXY", dtopo, 0, ndims - 1,
                            1}; // FXY lives on twisted edges
    aux_desc_arr[CORIOLISXYRECONVAR] = {
        "coriolisxyrecon", ptopo, 1, 0,
        2}; // coriolisxyrecon lives on primal edges
    aux_desc_arr[CORIOLISXYEDGERECONVAR] = {
        "corliolisxyedgerecon", ptopo, 2, 0,
        4}; // coriolisxyedgerecon lives on primal vertical faces

    aux_desc_arr[FTXYVAR] = {"FTXY", ptopo, 1, 0, 1}; // FTXY
    const_desc_arr[CORIOLISXYVAR] = {"coriolisxy", ptopo, 2, 0, 1};
  }
}

std::unique_ptr<TestCase> make_coupled_test_case(PamCoupler &coupler);
void testcase_from_string(std::unique_ptr<TestCase> &testcase, std::string name,
                          bool acoustic_balance);

void read_model_params_file(std::string inFile, ModelParameters &params,
                            Parallel &par, PamCoupler &coupler,
                            std::unique_ptr<TestCase> &testcase) {
#ifdef PAM_STANDALONE
  // Read common parameters
  int nz = coupler.get_nz();
  readParamsFile(inFile, params, par, nz);

  // Read config file
  YAML::Node config = YAML::LoadFile(inFile);

  params.acoustic_balance = config["balance_initial_density"].as<bool>(false);
  params.uniform_vertical = (config["vcoords"].as<std::string>() == "uniform");
  // Read diffusion coefficients
  params.scalar_diffusion_coeff = config["scalar_diffusion_coeff"].as<real>(0);
  params.scalar_diffusion_subtract_refstate =
      config["scalar_diffusion_subtract_refstate"].as<bool>(true);
  params.velocity_diffusion_coeff =
      config["velocity_diffusion_coeff"].as<real>(0);
  // Read the data initialization options
  params.init_data = config["init_data"].as<std::string>();
  params.force_refstate_hydrostatic_balance =
      config["force_refstate_hydrostatic_balance"].as<bool>(false);
  params.check_anelastic_constraint =
      config["check_anelastic_constraint"].as<bool>(false);

  for (int i = 0; i < ntracers_dycore; i++) {
    params.init_dycore_tracer[i] =
        config["init_dycore_tracer" + std::to_string(i)].as<std::string>(
            "constant");
    params.dycore_tracer_pos[i] =
        config["dycore_tracer" + std::to_string(i) + "_pos"].as<bool>(false);
  }

  // Store vertical cell interface heights in the data manager
  auto &dm = coupler.get_data_manager_device_readonly();
  params.zint = dm.get<real const, 2>("vertical_interface_height");

  params.ylen = 1.0;
  params.yc = 0.5;

  testcase_from_string(testcase, params.init_data, params.acoustic_balance);
#endif
}

void read_model_params_coupler(ModelParameters &params, Parallel &par,
                               PamCoupler &coupler,
                               std::unique_ptr<TestCase> &testcase) {

  // Read common parameters
  read_params_coupler(params, par, coupler);

  params.acoustic_balance = false;
  params.uniform_vertical = false;
  params.scalar_diffusion_coeff = 0;
  params.scalar_diffusion_subtract_refstate = true;
  params.velocity_diffusion_coeff = 0;
  params.init_data = "coupler";
  params.force_refstate_hydrostatic_balance = true;
  params.check_anelastic_constraint = false;

  // Store vertical cell interface heights in the data manager
  auto &dm = coupler.get_data_manager_device_readonly();
  params.zint = dm.get<real const, 2>("vertical_interface_height");

  params.ylen = 1.0;
  params.yc = 0.5;

  testcase = make_coupled_test_case(coupler);
}

void check_and_print_model_parameters(const ModelParameters &params,
                                      const Parallel &par,
                                      bool verbose = false) {

  check_and_print_parameters(params, par, verbose);

  if (verbose) {
    serial_print("IC: " + params.init_data, par.masterproc);
    serial_print("acoustically balanced: " +
                     std::to_string(params.acoustic_balance),
                 par.masterproc);
    serial_print("scalar_diffusion_coeff: " +
                     std::to_string(params.scalar_diffusion_coeff),
                 par.masterproc);
    serial_print("scalar_diffusion_subtract_refstate: " +
                     std::to_string(params.scalar_diffusion_subtract_refstate),
                 par.masterproc);
    serial_print("velocity_diffusion_coeff: " +
                     std::to_string(params.velocity_diffusion_coeff),
                 par.masterproc);

    for (int i = 0; i < ntracers_dycore; i++) {
      serial_print("Dycore Tracer" + std::to_string(i) +
                       " IC: " + params.init_dycore_tracer[i],
                   par.masterproc);
    }
  }
}
} // namespace pamc
