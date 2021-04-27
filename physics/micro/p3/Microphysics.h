
#pragma once

#include "awfl_const.h"
#include "DataManager.h"


extern "C"
void p3_main_fortran(double *qc , double *nc , double *qr , double *nr , double *th_atm , double *qv ,
                     double &dt , double *qi , double *qm , double *ni , double *bm , double *pres ,
                     double *dz , double *nc_nuceat_tend , double *nccn_prescribed , double *ni_activated ,
                     double *inv_qc_relvar , int &it , double *precip_liq_surf , double *precip_ice_surf ,
                     int &its , int &ite , int &kts , int &kte , double *diag_eff_radius_qc ,
                     double *diag_eff_radius_qi , double *rho_qi , bool &do_predict_nc , 
                     bool &do_prescribed_CCN ,double *dpres , double *exner , double *qv2qi_depos_tend ,
                     double *precip_total_tend , double *nevapr , double *qr_evap_tend ,
                     double *precip_liq_flux , double *precip_ice_flux , double *cld_frac_r ,
                     double *cld_frac_l , double *cld_frac_i , double *p3_tend_out , double *mu_c ,
                     double *lamc , double *liq_ice_exchange , double *vap_liq_exchange , 
                     double *vap_ice_exchange , double *qv_prev , double *t_prev , double *col_location ,
                     double *elapsed_s );



extern "C"
void micro_p3_utils_init_fortran(real &cpair , real &rair , real &rh2o , real &rhoh2o , real &mwh2o ,
                                 real &mwdry , real &gravit , real &latvap , real &latice , real &cpliq ,
                                 real &tmelt , real &pi , int &iulog , bool &masterproc );



extern "C"
void p3_init_fortran(char const *lookup_file_dir , int &dir_len , char const *version_p3 , int &ver_len );



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

  real grav;

  // This must be set during init() so we can return it in the get_water_vapor_index function
  int tracer_index_vapor;

  Constants constants;

  SArray<int,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  // Indices for all of your tracer quantities
  int static constexpr ID_C  = 0;  // Local index for Cloud Water Mass  
  int static constexpr ID_NC = 1;  // Local index for Cloud Water Number
  int static constexpr ID_R  = 2;  // Local index for Rain Water Mass   
  int static constexpr ID_NR = 3;  // Local index for Rain Water Number 
  int static constexpr ID_I  = 4;  // Local index for Ice Mass          
  int static constexpr ID_M  = 5;  // Local index for Ice Number        
  int static constexpr ID_NI = 6;  // Local index for Ice-Rime Mass     
  int static constexpr ID_BM = 7;  // Local index for Ice-Rime Volume   
  int static constexpr ID_V  = 8;  // Local index for Water Vapor       



  // TODO: Make sure the constants vibe with P3
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  Microphysics() {
    constants.R_d     = 287.;
    constants.cp_d    = 1003.;
    constants.cv_d    = constants.cp_d - constants.R_d;
    constants.gamma_d = constants.cp_d / constants.cv_d;
    constants.kappa_d = constants.R_d  / constants.cp_d;
    constants.R_v     = 461.;
    constants.cp_v    = 1859;
    constants.cv_v    = constants.R_v - constants.cp_v;
    constants.p0      = 1.e5;
    grav              = 9.81;
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
    tracer_IDs(ID_C ) = dycore.add_tracer(dm , "cloud_water"     , "Cloud Water Mass"   , true     , true);
    tracer_IDs(ID_NC) = dycore.add_tracer(dm , "cloud_water_num" , "Cloud Water Number" , true     , true);
    tracer_IDs(ID_R ) = dycore.add_tracer(dm , "rain"            , "Rain Water Mass"    , true     , true);
    tracer_IDs(ID_NR) = dycore.add_tracer(dm , "rain_num"        , "Rain Water Number"  , true     , true);
    tracer_IDs(ID_I ) = dycore.add_tracer(dm , "ice"             , "Ice Mass"           , true     , true);
    tracer_IDs(ID_NI) = dycore.add_tracer(dm , "ice_num"         , "Ice Number"         , true     , true);
    tracer_IDs(ID_M ) = dycore.add_tracer(dm , "ice_rime"        , "Ice-Rime Mass"      , true     , true);
    tracer_IDs(ID_BM) = dycore.add_tracer(dm , "ice_rime_vol"    , "Ice-Rime Volume"    , true     , true);
    tracer_IDs(ID_V ) = dycore.add_tracer(dm , "water_vapor"     , "Water Vapor"        , true     , true);

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

    tracer_index_vapor = tracer_IDs(ID_V);

    // Register and allocation non-tracer quantities used by the microphysics
    int p3_out = 49;
    dm.register_and_allocate<real>( "precip_liq_surf"    , "precipitation rate, liquid       m s-1"              , {            ny,nx,nens} , {           "y","x","nens"} );
    dm.register_and_allocate<real>( "precip_ice_surf"    , "precipitation rate, solid        m s-1"              , {            ny,nx,nens} , {           "y","x","nens"} );
    dm.register_and_allocate<real>( "diag_eff_radius_qc" , "effective radius, cloud          m"                  , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "diag_eff_radius_qi" , "effective radius, ice            m"                  , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "bulk_qi"            , "bulk density of ice              kg m-3"             , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "mu_c"               , "Size distribution shape parameter for radiation"     , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "lamc"               , "Size distribution slope parameter for radiation"     , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qv2qi_depos_tend"   , "qitend due to deposition/sublimation"                , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "precip_total_tend"  , "Total precipitation (rain + snow)"                   , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "nevapr"             , "evaporation of total precipitation (rain + snow)"    , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qr_evap_tend"       , "evaporation of rain"                                 , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "precip_liq_flux"    , "grid-box average rain flux (kg m^-2 s^-1) pverp"     , {       nz+1,ny,nx,nens} , {     "zp1","y","x","nens"} );
    dm.register_and_allocate<real>( "precip_ice_flux"    , "grid-box average ice/snow flux (kg m^-2 s^-1) pverp" , {       nz+1,ny,nx,nens} , {     "zp1","y","x","nens"} );
    dm.register_and_allocate<real>( "liq_ice_exchange"   , "sum of liq-ice phase change tendenices"              , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "vap_liq_exchange"   , "sum of vap-liq phase change tendenices"              , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "vap_ice_exchange"   , "sum of vap-ice phase change tendenices"              , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "p3_tend_out"        , "micro physics tendencies"                            , {p3_out,nz  ,ny,nx,nens} , {"p3","z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "qv_prev"            , "qv from the previous step"                           , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );
    dm.register_and_allocate<real>( "t_prev"             , "Temperature from the previous step"                  , {       nz  ,ny,nx,nens} , {     "z"  ,"y","x","nens"} );

    real rhoh2o = 1000.;
    real mwdry  = 28.966;
    real mwh2o  = 0.622 * mwdry;
    real latvap = 2.5E6;
    real latice = 3.50E5;
    real cpliq  = 4218.;
    real tmelt  = 273.15;
    real pi     = M_PI;
    int  iulog  = 1;
    bool masterproc = true;
    micro_p3_utils_init_fortran(constants.cp_d , constants.R_d , constants.R_v ,rhoh2o , mwh2o , mwdry , grav ,
                                latvap , latice, cpliq , tmelt , pi , iulog , masterproc);

    std::string dir = "../../physics/micro/p3";
    std::string ver = "4";
    int dir_len = dir.length();
    int ver_len = ver.length();
    p3_init_fortran( dir.c_str() , dir_len , ver.c_str() , ver_len );
  }



  void timeStep( DataManager &dm , real dt ) {
    // Get the dimensions sizes
    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    // Get tracers dimensioned as (nz,ny*nx*nens)
    auto rho_c  = dm.get_lev_col<real>("cloud_water"    );
    auto rho_nc = dm.get_lev_col<real>("cloud_water_num");
    auto rho_r  = dm.get_lev_col<real>("rain"           );
    auto rho_nr = dm.get_lev_col<real>("rain_num"       );
    auto rho_i  = dm.get_lev_col<real>("ice"            );
    auto rho_ni = dm.get_lev_col<real>("ice_num"        );
    auto rho_m  = dm.get_lev_col<real>("ice_rime"       );
    auto rho_bm = dm.get_lev_col<real>("ice_rime_vol"   );
    auto rho_v  = dm.get_lev_col<real>("water_vapor"    );

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
    auto precip_liq_surf    = dm.get_collapsed<real>("precip_liq_surf"   );
    auto precip_ice_surf    = dm.get_collapsed<real>("precip_ice_surf"   );
    auto diag_eff_radius_qc = dm.get_lev_col  <real>("diag_eff_radius_qc");
    auto diag_eff_radius_qi = dm.get_lev_col  <real>("diag_eff_radius_qi");
    auto bulk_qi            = dm.get_lev_col  <real>("bulk_qi"           );
    auto mu_c               = dm.get_lev_col  <real>("mu_c"              );
    auto lamc               = dm.get_lev_col  <real>("lamc"              );
    auto qv2qi_depos_tend   = dm.get_lev_col  <real>("qv2qi_depos_tend"  );
    auto precip_total_tend  = dm.get_lev_col  <real>("precip_total_tend" );
    auto nevapr             = dm.get_lev_col  <real>("nevapr"            );
    auto qr_evap_tend       = dm.get_lev_col  <real>("qr_evap_tend"      );
    auto liq_ice_exchange   = dm.get_lev_col  <real>("liq_ice_exchange"  );
    auto vap_liq_exchange   = dm.get_lev_col  <real>("vap_liq_exchange"  );
    auto vap_ice_exchange   = dm.get_lev_col  <real>("vap_ice_exchange"  );
    auto qv_prev            = dm.get_lev_col  <real>("qv_prev"           );
    auto t_prev             = dm.get_lev_col  <real>("t_prev"            );

    auto precip_liq_flux_dm = dm.get<real,4>("precip_liq_flux");
    auto precip_ice_flux_dm = dm.get<real,4>("precip_ice_flux");
    real2d precip_liq_flux("precip_liq_flux",nz+1,ny*nx*nens);
    real2d precip_ice_flux("precip_ice_flux",nz+1,ny*nx*nens);
    parallel_for( Bounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      precip_liq_flux(k,nx*nens*j + nens*i + iens) = precip_liq_flux_dm(k,j,i,iens);
      precip_ice_flux(k,nx*nens*j + nens*i + iens) = precip_ice_flux_dm(k,j,i,iens);
    });

    auto p3_tend_out_dm = dm.get<real,5>("p3_tend_out");
    real3d p3_tend_out("p3_tend_out",49,nz,ny*nx*nens);
    parallel_for( Bounds<5>(49,nz,ny,nx,nens) , YAKL_LAMBDA (int l, int k, int j, int i, int iens) {
      p3_tend_out(l,k,nx*nens*j + nens*i + iens) = p3_tend_out_dm(l,k,j,i,iens);
    });

    // These are inputs to p3
    real2d qc             ("qc"             ,nz,ncol);
    real2d nc             ("nc"             ,nz,ncol);
    real2d qr             ("qr"             ,nz,ncol);
    real2d nr             ("nr"             ,nz,ncol);
    real2d qi             ("qi"             ,nz,ncol);
    real2d ni             ("ni"             ,nz,ncol);
    real2d qm             ("qm"             ,nz,ncol);
    real2d bm             ("bm"             ,nz,ncol);
    real2d qv             ("qv"             ,nz,ncol);
    real2d pressure       ("pressure"       ,nz,ncol);
    real2d theta          ("theta"          ,nz,ncol);
    real2d exner          ("exner"          ,nz,ncol);
    real2d inv_exner      ("inv_exner"      ,nz,ncol);
    real2d dpres          ("dpres"          ,nz,ncol);
    real2d nc_nuceat_tend ("nc_nuceat_tend" ,nz,ncol);
    real2d nccn_prescribed("nccn_prescribed",nz,ncol);
    real2d ni_activated   ("ni_activated"   ,nz,ncol);
    real2d cld_frac_i     ("cld_frac_i"     ,nz,ncol);
    real2d cld_frac_l     ("cld_frac_l"     ,nz,ncol);
    real2d cld_frac_r     ("cld_frac_r"     ,nz,ncol);
    real2d inv_qc_relvar  ("inv_qc_relvar"  ,nz,ncol);
    real2d col_location   ("col_location"   ,3 ,ncol);

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
      qc           (k,i) = rho_c (k,i) / rho_dry(k,i);
      nc           (k,i) = rho_nc(k,i) / rho_dry(k,i);
      qr           (k,i) = rho_r (k,i) / rho_dry(k,i);
      nr           (k,i) = rho_nr(k,i) / rho_dry(k,i);
      qi           (k,i) = rho_i (k,i) / rho_dry(k,i);
      ni           (k,i) = rho_ni(k,i) / rho_dry(k,i);
      qm           (k,i) = rho_m (k,i) / rho_dry(k,i);
      bm           (k,i) = rho_bm(k,i) / rho_dry(k,i);
      qv           (k,i) = rho_v (k,i) / rho_dry(k,i);
      pressure     (k,i) = pressure_dry(k,i) + R_v * rho_v(k,i) * temp(k,i);
      exner        (k,i) = pow( pressure(k,i) / p0 , R_d / cp_d );
      inv_exner    (k,i) = 1. / exner(k,i);
      theta        (k,i) = temp(k,i) / exner(k,i);
      // P3 uses dpres to calculate density via the hydrostatic assumption.
      // So we just reverse this to compute dpres to give true density
      real rho = rho_dry(k,i) + rho_c(k,i) + rho_r(k,i) + rho_i(k,i) + rho_m(k,i) + rho_v(k,i);
      dpres        (k,i) = rho * grav * dz(k,i);
      // nc_nuceat_tend, nccn_prescribed, and ni_activated are not used
      // cld_frac_[lir] are set to one if there is mass; otherwise zero
      cld_frac_l   (k,i) = rho_c(k,i) > 0 ? 1 : 0;
      cld_frac_i   (k,i) = rho_i(k,i) > 0 ? 1 : 0;
      // cld_frac_r   (k,i) = rho_r(k,i) > 0 ? 1 : 0;
      // Peter Caldwell recommended setting cld_frac_r to 1 all the time
      cld_frac_r   (k,i) = 1;
      // inv_qc_relvar is always set to one
      inv_qc_relvar(k,i) = 1;
      // col_location is for debugging only, and it will be ignored for now
      if (k < 3) { col_location(k,i) = 1; }
    });

    int it = 1;
    int its = 1;
    int ite = ncol;
    int kts = 1;
    int kte = ncol;
    bool do_predict_nc = false;
    bool do_prescribed_CCN = false;
    double elapsed_s;

    p3_main_fortran(qc.data() , nc.data() , qr.data() , nr.data() , theta.data() , qv.data() , dt , qi.data() ,
                    qm.data() , ni.data() , bm.data() , pressure.data() , dz.data() , nc_nuceat_tend.data() ,
                    nccn_prescribed.data() , ni_activated.data() , inv_qc_relvar.data() , it ,
                    precip_liq_surf.data() , precip_ice_surf.data() , its , ite , kts , kte ,
                    diag_eff_radius_qc.data() , diag_eff_radius_qi.data() , bulk_qi.data() , do_predict_nc ,
                    do_prescribed_CCN , dpres.data() , inv_exner.data() , qv2qi_depos_tend.data() ,
                    precip_total_tend.data() , nevapr.data() , qr_evap_tend.data() , precip_liq_flux.data() ,
                    precip_ice_flux.data() , cld_frac_r.data() , cld_frac_l.data() , cld_frac_i.data() , p3_tend_out.data() , mu_c.data() ,
                    lamc.data() , liq_ice_exchange.data() , vap_liq_exchange.data() , vap_ice_exchange.data() , 
                    qv_prev.data() , t_prev.data() , col_location.data() , &elapsed_s );
                    
    ///////////////////////////////////////////////////////////////////////////////
    // Convert P3 outputs into dynamics coupler state and tracer masses
    ///////////////////////////////////////////////////////////////////////////////
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_c  (k,i) = qc(k,i)*rho_dry(k,i);
      rho_nc (k,i) = nc(k,i)*rho_dry(k,i);
      rho_r  (k,i) = qr(k,i)*rho_dry(k,i);
      rho_nr (k,i) = nr(k,i)*rho_dry(k,i);
      rho_i  (k,i) = qi(k,i)*rho_dry(k,i);
      rho_ni (k,i) = ni(k,i)*rho_dry(k,i);
      rho_m  (k,i) = qm(k,i)*rho_dry(k,i);
      rho_bm (k,i) = bm(k,i)*rho_dry(k,i);
      rho_v  (k,i) = qv(k,i)*rho_dry(k,i);
      temp   (k,i) = theta(k,i) * exner(k,i);
      // Save qv and temperature for the next call to p3_main
      qv_prev(k,i) = qv  (k,i);
      t_prev (k,i) = temp(k,i);
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



