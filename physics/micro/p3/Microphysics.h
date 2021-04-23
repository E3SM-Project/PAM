
#pragma once

#include "awfl_const.h"
#include "DataManager.h"

class Microphysics {
public:
  // Doesn't actually have to be static or constexpr. Could be assigned in the constructor
  int static constexpr num_tracers = 9;

  // You should set these int he constructor
  struct Constants {
    real R_d    ;
    real cp_d   ;
    real cv_d   ;
    real gamma_d;
    real kappa_d;
    real R_v    ;
    real cp_v   ;
    real cv_v   ;
    real p0     ;
  };

  // This must be set during init() so we can return it in the get_water_vapor_index function
  int tracer_index_vapor;

  Constants constants;

  // TODO: Change this to type int instead of real
  SArray<real,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  // Indices for all of your tracer quantities
  int static constexpr ID_QC = 0;  // Local index for water vapor
  int static constexpr ID_NC = 1;  // Local index for cloud liquid
  int static constexpr ID_QR = 2;  // Local index for precipitated liquid (rain)
  int static constexpr ID_NR = 3;  // Local index for precipitated liquid (rain)
  int static constexpr ID_QI = 4;  // Local index for precipitated liquid (rain)
  int static constexpr ID_QM = 5;  // Local index for precipitated liquid (rain)
  int static constexpr ID_NI = 6;  // Local index for precipitated liquid (rain)
  int static constexpr ID_BM = 7;  // Local index for precipitated liquid (rain)
  int static constexpr ID_QV = 8;  // Local index for precipitated liquid (rain)



  // TODO: Make sure the constants vibe with P3
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  Microphysics() {
    constants.R_d         = 287.;
    constants.cp_d        = 1003.;
    constants.cv_d        = constants.cp_d - constants.R_d;
    constants.gamma_d     = constants.cp_d / constants.cv_d;
    constants.kappa_d     = constants.R_d  / constants.cp_d;
    constants.R_v         = 461.;
    constants.cp_v        = 1859;
    constants.cv_v        = constants.R_v - constants.cp_v;
    constants.p0          = 1.e5;
  }



  // This must return the correct # of tracers **BEFORE** init(...) is called
  int get_num_tracers() const {
    return num_tracers;
  }



  // This must return the correct index of water vapor **AFTER** init(...) is called
  int get_water_vapor_index() const {
    return tracer_index_vapor;
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  // and storing the water vapor tracer index
  template <class DC>
  void init(std::string infile , int ny, int nx, int nens , DC &dycore , DataManager &dm) {
    int nz = dm.get_dimension_size("z");

    // Register tracers in the dycore
    //                                        name                 description            positive   adds mass
    tracer_IDs(ID_QC) = dycore.add_tracer(dm , "cloud_water"     , "Cloud Water Mass"   , true     , true);
    tracer_IDs(ID_NC) = dycore.add_tracer(dm , "cloud_water_num" , "Cloud Water Number" , true     , true);
    tracer_IDs(ID_QR) = dycore.add_tracer(dm , "rain"            , "Rain Water Mass"    , true     , true);
    tracer_IDs(ID_NR) = dycore.add_tracer(dm , "rain_num"        , "Rain Water Number"  , true     , true);
    tracer_IDs(ID_QI) = dycore.add_tracer(dm , "ice"             , "Ice Mass"           , true     , true);
    tracer_IDs(ID_NI) = dycore.add_tracer(dm , "ice_num"         , "Ice Number"         , true     , true);
    tracer_IDs(ID_QM) = dycore.add_tracer(dm , "ice_rime"        , "Ice-Rime Mass"      , true     , true);
    tracer_IDs(ID_BM) = dycore.add_tracer(dm , "ice_rime_vol"    , "Ice-Rime Volume"    , true     , true);
    tracer_IDs(ID_QV) = dycore.add_tracer(dm , "water_vapor"     , "Water Vapor"        , true     , true);

    // Register and allocate the tracers in the DataManager
    dm.register_and_allocate<real>( "cloud_water"     , "Cloud Water Mass"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cloud_water_num" , "Cloud Water Number" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "rain"            , "Rain Water Mass"    , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "rain_num"        , "Rain Water Number"  , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice"             , "Ice Mass"           , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice_num"         , "Ice Number"         , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice_rime"        , "Ice-Rime Mass"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "ice_rime_vol"    , "Ice-Rime Volume"    , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "water_vapor"     , "Water Vapor"        , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    tracer_index_vapor = tracer_IDs(ID_QV);

    // Register and allocation non-tracer quantities used by the microphysics
    dm.register_and_allocate<real>( "precip_liq_surf"    , "precipitation rate, liquid       m s-1"              , {        ny,nx,nens} , {           "y","x","nens"} );
    dm.register_and_allocate<real>( "precip_ice_surf"    , "precipitation rate, solid        m s-1"              , {        ny,nx,nens} , {           "y","x","nens"} );
    dm.register_and_allocate<real>( "diag_eff_radius_qc" , "effective radius, cloud          m"                  , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "diag_eff_radius_qi" , "effective radius, ice            m"                  , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "rho_qi"             , "bulk density of ice              kg m-3"             , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "mu_c"               , "Size distribution shape parameter for radiation"     , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "lamc"               , "Size distribution slope parameter for radiation"     , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qv2qi_depos_tend"   , "qitend due to deposition/sublimation"                , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "precip_total_tend"  , "Total precipitation (rain + snow)"                   , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "nevapr"             , "evaporation of total precipitation (rain + snow)"    , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qr_evap_tend"       , "evaporation of rain"                                 , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "precip_liq_flux"    , "grid-box average rain flux (kg m^-2 s^-1) pverp"     , {   nz+1,ny,nx,nens} , {     "zp1","y","x","nens"} );
    dm.register_and_allocate<real>( "precip_ice_flux"    , "grid-box average ice/snow flux (kg m^-2 s^-1) pverp" , {   nz+1,ny,nx,nens} , {     "zp1","y","x","nens"} );
    dm.register_and_allocate<real>( "liq_ice_exchange"   , "sum of liq-ice phase change tendenices"              , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "vap_liq_exchange"   , "sum of vap-liq phase change tendenices"              , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "vap_ice_exchange"   , "sum of vap-ice phase change tendenices"              , {   nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "p3_tend_out"        , "micro physics tendencies"                            , {49,nz  ,ny,nx,nens} , {"49","z"  ,"y","x","nens"} );
  }



  void timeStep( DataManager &dm , real dt ) {
    // Get the dimensions sizes
    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    // Get tracers dimensioned as (nz,ny*nx*nens)
    auto rho_qc = dm.get_lev_col<real>("cloud_water"    );
    auto rho_nc = dm.get_lev_col<real>("cloud_water_num");
    auto rho_qr = dm.get_lev_col<real>("rain"           );
    auto rho_nr = dm.get_lev_col<real>("rain_num"       );
    auto rho_qi = dm.get_lev_col<real>("ice"            );
    auto rho_ni = dm.get_lev_col<real>("ice_num"        );
    auto rho_qm = dm.get_lev_col<real>("ice_rime"       );
    auto rho_bm = dm.get_lev_col<real>("ice_rime_vol"   );
    auto rho_qv = dm.get_lev_col<real>("water_vapor"    );

    // Get coupler state
    auto rho_dry      = dm.get_lev_col<real>("density_dry");
    auto temp         = dm.get_lev_col<real>("temp");
    auto pressure_dry = dm.get_lev_col<real>("pressure_dry");

    // Calculate the grid spacing
    auto zint_in = dm.get<real,2>("vertical_interface_height");
    real2d dz("dz",nz,ny*nx*nens);
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      dz(k,j*nx*nens + i*nens + iens) = zint_in(k+1,iens) - zint_in(k,iens);
    });

    // Get everything from the DataManager that's not a tracer but is persistent across multiple micro calls
    auto precip_liq_surf    = dm.get_collapsed<real>("precip_liq_surf"    );
    auto precip_ice_surf    = dm.get_collapsed<real>("precip_ice_surf"    );
    auto diag_eff_radius_qc = dm.get_lev_col  <real>("diag_eff_radius_qc" );
    auto diag_eff_radius_qi = dm.get_lev_col  <real>("diag_eff_radius_qi" );
    auto rho_qi             = dm.get_lev_col  <real>("rho_qi"             );
    auto mu_c               = dm.get_lev_col  <real>("mu_c"               );
    auto lamc               = dm.get_lev_col  <real>("lamc"               );
    auto qv2qi_depos_tend   = dm.get_lev_col  <real>("qv2qi_depos_tend"   );
    auto precip_total_tend  = dm.get_lev_col  <real>("precip_total_tend"  );
    auto nevapr             = dm.get_lev_col  <real>("nevapr"             );
    auto qr_evap_tend       = dm.get_lev_col  <real>("qr_evap_tend"       );
    auto liq_ice_exchange   = dm.get_lev_col  <real>("liq_ice_exchange"   );
    auto vap_liq_exchange   = dm.get_lev_col  <real>("vap_liq_exchange"   );
    auto vap_ice_exchange   = dm.get_lev_col  <real>("vap_ice_exchange"   );

    auto precip_liq_flux_dm = dm.get<real,4>("precip_liq_flux");
    auto precip_ice_flux_dm = dm.get<real,4>("precip_ice_flux");
    real2d precip_liq_flux("precip_liq_flux",nz+1,ny*nx*nens);
    real2d precip_ice_flux("precip_ice_flux",nz+1,ny*nx*nens);
    parallel_for( Bounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      precip_liq_flux(k,nx*nens*j + nens*i + iens) = precip_liq_flux_dm(k,j,i,iens);
      precip_ice_flux(k,nx*nens*j + nens*i + iens) = precip_ice_flux_dm(k,j,i,iens);
    });

    auto p3_tend_out_dm = dm.get<real,4>("p3_tend_out");
    real3d p3_tend_out("p3_tend_out",49,nz,ny*nx*nens);
    parallel_for( Bounds<5>(49,nz+1,ny,nx,nens) , YAKL_LAMBDA (int l, int k, int j, int i, int iens) {
      p3_tend_out(l,k,nx*nens*j + nens*i + iens) = p3_tend_out_dm(l,k,j,i,iens);
    });

    // These are inputs to p3
    real2d qc          ("qc"          ,nz,ncol);
    real2d nc          ("nc"          ,nz,ncol);
    real2d qr          ("qr"          ,nz,ncol);
    real2d nr          ("nr"          ,nz,ncol);
    real2d qi          ("qi"          ,nz,ncol);
    real2d ni          ("ni"          ,nz,ncol);
    real2d qm          ("qm"          ,nz,ncol);
    real2d bm          ("bm"          ,nz,ncol);
    real2d qv          ("qv"          ,nz,ncol);
    real2d theta_dry   ("theta_dry"   ,nz,ncol);
    real2d exner_dry   ("exner_dry"   ,nz,ncol);

    //////////////////////////////////////////////////////////////////////////////
    // Compute quantities needed for inputs to P3
    //////////////////////////////////////////////////////////////////////////////
    // Force constants into local scope
    real gamma_d = this->constants.gamma_d;
    real R_d     = this->constants.R_d;
    real R_v     = this->constants.R_v;
    real cp_d    = this->constants.cp_d;
    real p0      = this->constants.p0;
    // Save initial state, and compute inputs for kessler(...)
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      qc          (k,i) = rho_qc(k,i) / rho_dry(k,i);
      nc          (k,i) = rho_nc(k,i) / rho_dry(k,i);
      qr          (k,i) = rho_qr(k,i) / rho_dry(k,i);
      nr          (k,i) = rho_nr(k,i) / rho_dry(k,i);
      qi          (k,i) = rho_qi(k,i) / rho_dry(k,i);
      ni          (k,i) = rho_ni(k,i) / rho_dry(k,i);
      qm          (k,i) = rho_qm(k,i) / rho_dry(k,i);
      bm          (k,i) = rho_bm(k,i) / rho_dry(k,i);
      qv          (k,i) = rho_qv(k,i) / rho_dry(k,i);
      exner_dry   (k,i) = pow( pressure_dry(k,i) / p0 , R_d / cp_d );
      theta_dry   (k,i) = temp(k,i) / exner_dry(k,i);
    });




  // TODO: Find out if we need full or dry pressure and exner and theta
  SUBROUTINE p3_main(qc,nc,qr,nr,th_atm,qv,dt,qi,qm,ni,bm,   &
       pres,dz,nc_nuceat_tend,nccn_prescribed,ni_activated,inv_qc_relvar,it,precip_liq_surf,precip_ice_surf,its,ite,kts,kte,diag_eff_radius_qc,     &
       diag_eff_radius_qi,rho_qi,do_predict_nc, do_prescribed_CCN, &
       dpres,exner,qv2qi_depos_tend,precip_total_tend,nevapr,qr_evap_tend,precip_liq_flux,precip_ice_flux,cld_frac_r,cld_frac_l,cld_frac_i,  &
       p3_tend_out,mu_c,lamc,liq_ice_exchange,vap_liq_exchange, &
       vap_ice_exchange,qv_prev,t_prev,col_location &
       ,elapsed_s &
      )
    implicit none
    !----- Input/ouput arguments:  ----------------------------------------------------------!
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: nc_nuceat_tend      ! IN ccn activated number tendency kg-1 s-1
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: nccn_prescribed
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: ni_activated       ! IN actived ice nuclei concentration  1/kg
    real(rtype), intent(in)                                     :: dt         ! model time step                  s
    integer, intent(in)                                  :: its,ite    ! array bounds (horizontal)
    integer, intent(in)                                  :: kts,kte    ! array bounds (vertical)
    integer, intent(in)                                  :: it         ! time step counter NOTE: starts at 1 for first time step
    logical(btype), intent(in)                           :: do_predict_nc ! .T. (.F.) for prediction (specification) of Nc
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: dpres       ! pressure thickness               Pa
    logical(btype), intent(in)                                  :: do_prescribed_CCN
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: cld_frac_i, cld_frac_l, cld_frac_r ! Ice, Liquid and Rain cloud fraction
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: qv_prev, t_prev                    ! qv and t from previous p3_main call
    real(rtype), intent(in),    dimension(its:ite,3)            :: col_location
    real(rtype), intent(in),    dimension(its:ite,kts:kte)      :: inv_qc_relvar
    real(rtype), intent(out) :: elapsed_s ! duration of main loop in seconds





    ///////////////////////////////////////////////////////////////////////////////
    // Convert P3 outputs into dynamics coupler state and tracer masses
    ///////////////////////////////////////////////////////////////////////////////
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_qc(k,i) = qc(k,i)*rho_dry(k,i);
      rho_nc(k,i) = nc(k,i)*rho_dry(k,i);
      rho_qr(k,i) = qr(k,i)*rho_dry(k,i);
      rho_nr(k,i) = nr(k,i)*rho_dry(k,i);
      rho_qi(k,i) = qi(k,i)*rho_dry(k,i);
      rho_ni(k,i) = ni(k,i)*rho_dry(k,i);
      rho_qm(k,i) = qm(k,i)*rho_dry(k,i);
      rho_bm(k,i) = bm(k,i)*rho_dry(k,i);
      rho_qv(k,i) = qv(k,i)*rho_dry(k,i);
      temp (k,i) = theta_dry(k,i) * exner_dry(k,i);
    });

  }



  // These are outputs that are not tracer mass. Tracer mass is handled by the dycore instead
  // This assumes the NetCDF handler "nc" is already open and will be closed later
  // This is for a single ensemble index
  void output(DataManager &dm, yakl::SimpleNetCDF &nc, int ulIndex, int iens) const {
  }



  std::string micro_name() const {
    return "p3";
  }



};



