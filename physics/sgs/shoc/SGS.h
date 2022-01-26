
#pragma once

#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"

using pam::PamCoupler;


extern "C"
void shoc_init(int &nlev, double &gravit, double &rair, double &rh2o, double &cpair, double &zvir, double &latvap,
               double &latice, double &karman, double *pref_mid, int &nbot_shoc, int &ntop_shoc);


extern "C"
void shoc_main(int &shcol, int &nlev, int &nlevi, double *dtime, int &nadv, real *host_dx, double *host_dy,
               double *thv, double *zt_grid, double *zi_grid, double *pres, double *presi, double *pdel,
               double *wthl_sfc, double *wqw_sfc, double *uw_sfc, double *vw_sfc, double *wtracer_sfc,
               int &num_qtracers, double *w_field, double *exner, double *phis, double *host_dse, double *tke,
               double *thetal, double *qw, double *u_wind, double *v_wind, double *qtracers, double *wthv_sec,
               double *tkh, double *tk, double *shoc_ql, double *shoc_cldfrac, double *pblh, double *shoc_mix,
               double *isotropy, double *w_sec, double *thl_sec, double *qw_sec, double *qwthl_sec,
               double *wthl_sec, double *wqw_sec, double *wtke_sec, double *uw_sec, double *vw_sec, double *w3,
               double *wqls_sec, double *brunt, double *shoc_ql2 );


class SGS {
public:
  // Doesn't actually have to be static or constexpr. Could be assigned in the constructor
  int static constexpr num_tracers = 1;

  // You should set these in the constructor
  real R_d    ;
  real cp_d   ;
  real cv_d   ;
  real gamma_d;
  real kappa_d;
  real R_v    ;
  real cp_v   ;
  real cv_v   ;
  real p0     ;

  real grav;
  real cp_l;

  bool first_step;

  // Indices for all of your tracer quantities
  int static constexpr ID_TKE  = 0;  // Local index for Turbulent Kinetic Energy (m^2/s^2)



  // TODO: Make sure the constants vibe with P3
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  SGS() {
    R_d        = 287.;
    cp_d       = 1003.;
    cv_d       = cp_d - R_d;
    gamma_d    = cp_d / cv_d;
    kappa_d    = R_d  / cp_d;
    R_v        = 461.;
    cp_v       = 1859;
    cv_v       = R_v - cp_v;
    p0         = 1.e5;
    grav       = 9.81;
    first_step = true;
    cp_l       = 4218.;
  }



  // This must return the correct # of tracers **BEFORE** init(...) is called
  YAKL_INLINE static int get_num_tracers() {
    return num_tracers;
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  void init(PamCoupler &coupler) {
    int nx   = coupler.get_nx  ();
    int ny   = coupler.get_ny  ();
    int nz   = coupler.get_nz  ();
    int nens = coupler.get_nens();

    // Register tracers in the coupler
    //                 name    description                              positive   adds mass
    coupler.add_tracer("tke" , "Turbulent Kinetic Energy (m^2/s^2)"   , true     , false );

    // Register and allocation non-tracer quantities used by the microphysics
    int p3_nout = 49;
    coupler.dm.register_and_allocate<real>( "wthv_sec"     , "Buoyancy flux [K m/s]"                , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "tk"           , "Eddy coefficient for momentum [m2/s]" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "tkh"          , "Eddy coefficent for heat [m2/s]"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    coupler.dm.register_and_allocate<real>( "shoc_cldfrac" , "Cloud fraction [-]"                   , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    auto tke          = coupler.dm.get<real,4>( "tke"          );
    auto wthv_sec     = coupler.dm.get<real,4>( "wthv_sec"     );
    auto tk           = coupler.dm.get<real,4>( "tk"           );
    auto tkh          = coupler.dm.get<real,4>( "tkh"          );
    auto shoc_cldfrac = coupler.dm.get<real,4>( "shoc_cldfrac" );

    parallel_for( "sgs zero" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tke         (k,j,i,iens) = 0;
      wthv_sec    (k,j,i,iens) = 0;
      tk          (k,j,i,iens) = 0;
      tkh         (k,j,i,iens) = 0;
      shoc_cldfrac(k,j,i,iens) = 0;
    });
  }



  void timeStep( PamCoupler &coupler , real dt , real etime ) {
    // Get the dimensions sizes
    int nz   = coupler.get_nz();
    int ny   = coupler.get_ny();
    int nx   = coupler.get_nx();
    int nens = coupler.get_nens();
    int ncol = ny*nx*nens;

    auto pres     = coupler.compute_pressure_array(); 
    auto pres_int = coupler.interp_pressure_interfaces( pres );

    // SHOC init requires reference pressure, which we do not have available for the init() call
    if (first_step) {
      real constexpr latvap = 2.5E6 ;
      real constexpr latice = 3.50E5;
      real constexpr karman = 0.4;
      // TODO: See if a more appropriate reference pressure should be used
      // pref is only used to limit PBL height to 400mb, so this should be OK
      // Invert the first column in x, z, and ensemble to use as reference pressure for shoc
      real1d pref_shoc("pref_shoc",nz);
      parallel_for( Bounds<2>(nz) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        pref_shoc(k,j,i,iens) = pres(nz-1-k,0,0,0);
      });
      shoc_init( nz , grav , R_d , R_v , cp_d , R_v / R_d -1 , latvap , latice , karman ,
                 pref_shoc.createHostCopy().data() , nlev , 1 );
      first_step = false;
    }

    // Get saved SHOC-related variables
    auto tke          = coupler.dm.get_lev_col<real>( "tke"          ); // PAM Tracer
    auto wthv_sec     = coupler.dm.get_lev_col<real>( "wthv_sec"     ); // Reuse from last SHOC output
    auto tk           = coupler.dm.get_lev_col<real>( "tk"           ); // Reuse from last SHOC output
    auto tkh          = coupler.dm.get_lev_col<real>( "tkh"          ); // Reuse from last SHOC output
    auto shoc_cldfrac = coupler.dm.get_lev_col<real>( "shoc_cldfrac" ); // Reuse from last SHOC output
    // Get coupler state
    auto rho_d        = coupler.dm.get_lev_col<real>( "density_dry"  );
    auto uvel         = coupler.dm.get_lev_col<real>( "uvel"         );
    auto vvel         = coupler.dm.get_lev_col<real>( "vvel"         );
    auto wvel         = coupler.dm.get_lev_col<real>( "wvel"         );
    auto temp         = coupler.dm.get_lev_col<real>( "temp"         );
    auto rho_v        = coupler.dm.get_lev_col<real>( "water_vapor"  );

    // Calculate the grid spacing
    auto zint_in = coupler.dm.get<real,2>("vertical_interface_height");

    int num_qtracers = ???;

    // Create variables for SHOC call (these are all inverted for vertical indices
    real1d host_dx     ("host_dx     ",                  ncol); // grid spacing of host model in x direction [m]
    real1d host_dy     ("host_dy     ",                  ncol); // grid spacing of host model in y direction [m]
    real2d zt_grid     ("zt_grid     ",             nz  ,ncol); // heights, for thermo grid [m]
    real2d zi_grid     ("zi_grid     ",             nz+1,ncol); // heights, for interface grid [m]
    real2d pres        ("pres        ",             nz  ,ncol); // pressure levels on thermo grid [Pa]
    real2d presi       ("presi       ",             nz+1,ncol); // pressure levels on interface grid [Pa]
    real2d pdel        ("pdel        ",             nz  ,ncol); // Differences in pressure levels [Pa]
    real2d thv         ("thv         ",             nz  ,ncol); // virtual potential temperature [K]
    real2d w_field     ("w_field     ",             nz  ,ncol); // large scale vertical velocity [m/s]
    real1d wthl_sfc    ("wthl_sfc    ",                  ncol); // Surface sensible heat flux [K m/s]
    real1d wqw_sfc     ("wqw_sfc     ",                  ncol); // Surface latent heat flux [kg/kg m/s]
    real1d uw_sfc      ("uw_sfc      ",                  ncol); // Surface momentum flux (u-direction) [m2/s2]
    real1d vw_sfc      ("vw_sfc      ",                  ncol); // Surface momentum flux (v-direction) [m2/s2]
    real2d wtracer_sfc ("wtracer_sfc ",num_qtracers,     ncol); // Surface flux for tracers [varies]
    real2d exner       ("exner       ",             nz  ,ncol); // Exner function [-]
    real1d phis        ("phis        ",                  ncol); // Host model surface geopotential height
    real2d host_dse    ("host_dse    ",             nz  ,ncol); // dry static energy [J/kg];  dse = Cp*T + g*z + phis
    real2d tke         ("tke         ",             nz  ,ncol); // turbulent kinetic energy [m2/s2]
    real2d thetal      ("thetal      ",             nz  ,ncol); // liquid water potential temperature [K]
    real2d qw          ("qw          ",             nz  ,ncol); // total water mixing ratio [kg/kg]
    real2d u_wind      ("u_wind      ",             nz  ,ncol); // u wind component [m/s]
    real2d v_wind      ("v_wind      ",             nz  ,ncol); // v wind component [m/s]
    real2d wthv_sec    ("wthv_sec    ",             nz  ,ncol); // buoyancy flux [K m/s]
    real3d qtracers    ("qtracers    ",num_qtracers,nz  ,ncol); // tracers [varies]
    real2d tk          ("tk          ",             nz  ,ncol); // eddy coefficient for momentum [m2/s]
    real2d tkh         ("tkh         ",             nz  ,ncol); // eddy coefficent for heat [m2/s]
    real2d shoc_cldfrac("shoc_cldfrac",             nz  ,ncol); // Cloud fraction [-]
    real2d shoc_ql     ("shoc_ql     ",             nz  ,ncol); // cloud liquid mixing ratio [kg/kg]
    real1d pblh        ("pblh        ",                  ncol); // planetary boundary layer depth [m]
    real2d shoc_ql2    ("shoc_ql2    ",             nz  ,ncol); // cloud liquid mixing ratio variance [kg^2/kg^2]
    real2d shoc_mix    ("shoc_mix    ",             nz  ,ncol); // Turbulent length scale [m]
    real2d w_sec       ("w_sec       ",             nz  ,ncol); // vertical velocity variance [m2/s2]
    real2d thl_sec     ("thl_sec     ",             nz+1,ncol); // temperature variance [K^2]
    real2d qw_sec      ("qw_sec      ",             nz+1,ncol); // moisture variance [kg2/kg2]
    real2d qwthl_sec   ("qwthl_sec   ",             nz+1,ncol); // temp moisture covariance [K kg/kg]
    real2d wthl_sec    ("wthl_sec    ",             nz+1,ncol); // vertical heat flux [K m/s]
    real2d wqw_sec     ("wqw_sec     ",             nz+1,ncol); // vertical moisture flux [K m/s]
    real2d wtke_sec    ("wtke_sec    ",             nz+1,ncol); // vertical tke flux [m3/s3]
    real2d uw_sec      ("uw_sec      ",             nz+1,ncol); // vertical zonal momentum flux [m2/s2]
    real2d vw_sec      ("vw_sec      ",             nz+1,ncol); // vertical meridional momentum flux [m2/s2]
    real2d w3          ("w3          ",             nz+1,ncol); // third moment vertical velocity [m3/s3]
    real2d wqls_sec    ("wqls_sec    ",             nz  ,ncol); // liquid water flux [kg/kg m/s]
    real2d brunt       ("brunt       ",             nz  ,ncol); // brunt vaisala frequency [s-1]
    real2d isotropy    ("isotropy    ",             nz  ,ncol); // return to isotropic timescale [s]





//      subroutine shoc_main ( &
//           shcol, nlev, nlevi, dtime, nadv, &   ! Input
//           host_dx, host_dy,thv, &              ! Input
//           zt_grid,zi_grid,pres,presi,pdel,&    ! Input
//           wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, & ! Input
//           wtracer_sfc,num_qtracers,w_field, &  ! Input
//           exner,phis, &                        ! Input
//           host_dse, tke, thetal, qw, &         ! Input/Output
//           u_wind, v_wind,qtracers,&            ! Input/Output
//           wthv_sec,tkh,tk,&                    ! Input/Output
//           shoc_ql,shoc_cldfrac,&               ! Input/Output
//           pblh,&                               ! Output
//           shoc_mix, isotropy,&                 ! Output (diagnostic)
//           w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
//           wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
//           uw_sec, vw_sec, w3,&                 ! Output (diagnostic)
//           wqls_sec, brunt, shoc_ql2 &          ! Output (diagnostic)
//        integer    , intent(in   ) :: shcol                                  !UNDERSTOOD   number of SHOC columns in the array
//        integer    , intent(in   ) :: nlev                                   !UNDERSTOOD   number of levels [-]
//        integer    , intent(in   ) :: nlevi                                  !UNDERSTOOD   number of levels on interface grid [-]
//        integer    , intent(in   ) :: num_qtracers                           !UNDERSTOOD   number of tracers [-]
//        integer    , intent(in   ) :: nadv                                   !UNDERSTOOD   THIS WILL MOST LIKELY BE 1 ;  number of times to loop SHOC
//        real(rtype), intent(in   ) :: dtime                                  !UNDERSTOOD   NAMELISTS HAVE DTIME = 150s;  SHOC timestep [s]
//        real(rtype), intent(in   ) :: host_dx     (shcol                   ) !UNDERSTOOD   grid spacing of host model in x direction [m]
//        real(rtype), intent(in   ) :: host_dy     (shcol                   ) !UNDERSTOOD   grid spacing of host model in y direction [m]
//        real(rtype), intent(in   ) :: zt_grid     (shcol,nlev              ) !UNDERSTOOD   heights, for thermo grid [m]
//        real(rtype), intent(in   ) :: zi_grid     (shcol,nlevi             ) !UNDERSTOOD   heights, for interface grid [m]
//        real(rtype), intent(in   ) :: pres        (shcol,nlev              ) !UNDERSTOOD   pressure levels on thermo grid [Pa]
//        real(rtype), intent(in   ) :: presi       (shcol,nlevi             ) !UNDERSTOOD   pressure levels on interface grid [Pa]
//        real(rtype), intent(in   ) :: pdel        (shcol,nlev              ) !UNDERSTOOD   Differences in pressure levels [Pa]
//        real(rtype), intent(in   ) :: thv         (shcol,nlev              ) !UNDERSTOOD   https://glossary.ametsoc.org/wiki/Virtual_potential_temperature ??  virtual potential temperature [K]
//        real(rtype), intent(in   ) :: w_field     (shcol,nlev              ) !UNDERSTOOD   large scale vertical velocity [m/s]
//        real(rtype), intent(in   ) :: wthl_sfc    (shcol                   ) !UNDERSTOOD   Surface sensible heat flux [K m/s]
//        real(rtype), intent(in   ) :: wqw_sfc     (shcol                   ) !UNDERSTOOD   SURFACE FLUX OF TOTAL WATER DRY MIXING RATIO (VAPOR + CLOUD LIQUID);  Surface latent heat flux [kg/kg m/s]
//        real(rtype), intent(in   ) :: uw_sfc      (shcol                   ) !UNDERSTOOD   Surface momentum flux (u-direction) [m2/s2]
//        real(rtype), intent(in   ) :: vw_sfc      (shcol                   ) !UNDERSTOOD   Surface momentum flux (v-direction) [m2/s2]
//        real(rtype), intent(in   ) :: wtracer_sfc (shcol      ,num_qtracers) !UNDERSTOOD   DRY MIXING RATIOS, SAME DEFINITIONS AS QTRACERS;  Surface flux for tracers [varies]
//        real(rtype), intent(in   ) :: exner       (shcol,nlev              ) !UNDERSTOOD   Exner function [-]
//        real(rtype), intent(in   ) :: phis        (shcol                   ) !UNDERSTOOD   Host model surface geopotential height
//        real(rtype), intent(inout) :: host_dse    (shcol,nlev              ) !UNDERSTOOD   prognostic temp variable of host model;  dry static energy [J/kg];  dse = Cp*T + g*z + phis
//        real(rtype), intent(inout) :: tke         (shcol,nlev              ) !UNDERSTOOD   turbulent kinetic energy [m2/s2]
//        real(rtype), intent(inout) :: thetal      (shcol,nlev              ) !UNDERSTOOD   https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature ??  liquid water potential temperature [K]
//        real(rtype), intent(inout) :: qw          (shcol,nlev              ) !UNDERSTOOD   VAPOR + CLOUD LIQUID;   total water mixing ratio [kg/kg]
//        real(rtype), intent(inout) :: u_wind      (shcol,nlev              ) !UNDERSTOOD   u wind component [m/s]
//        real(rtype), intent(inout) :: v_wind      (shcol,nlev              ) !UNDERSTOOD   v wind component [m/s]
//        real(rtype), intent(inout) :: wthv_sec    (shcol,nlev              ) !UNDERSTOOD   https://glossary.ametsoc.org/wiki/Buoyancy_flux ??   buoyancy flux [K m/s]
//        real(rtype), intent(inout) :: qtracers    (shcol,nlev ,num_qtracers) !UNDERSTOOD   THIS IS DRY MIXING RATIOS. EXCLUDE QV AND QC         tracers [varies]
//        real(rtype), intent(inout) :: tk          (shcol,nlev              ) !UNDERSTOOD   INIT TO ZERO, THEN LET SHOC HANDLE IT FROM THERE     eddy coefficient for momentum [m2/s]
//        real(rtype), intent(inout) :: tkh         (shcol,nlev              ) !UNDERSTOOD   INIT TO ZERO, THEN LET SHOC HANDLE IT FROM THERE     eddy coefficent for heat [m2/s]
//        real(rtype), intent(inout) :: shoc_cldfrac(shcol,nlev              ) !UNDERSTOOD   INIT TO ZERO, THEN LET SHOC HANDLE IT FROM THERE     Cloud fraction [-]
//        real(rtype), intent(inout) :: shoc_ql     (shcol,nlev              ) !UNDERSTOOD   cloud liquid mixing ratio [kg/kg]
//        real(rtype), intent(  out) :: pblh        (shcol                   ) ! planetary boundary layer depth [m]
//        real(rtype), intent(  out) :: shoc_ql2    (shcol,nlev              ) ! cloud liquid mixing ratio variance [kg^2/kg^2]
//        real(rtype), intent(  out) :: shoc_mix    (shcol,nlev              ) ! Turbulent length scale [m]
//        real(rtype), intent(  out) :: w_sec       (shcol,nlev              ) ! vertical velocity variance [m2/s2]
//        real(rtype), intent(  out) :: thl_sec     (shcol,nlevi             ) ! temperature variance [K^2]
//        real(rtype), intent(  out) :: qw_sec      (shcol,nlevi             ) ! moisture variance [kg2/kg2]
//        real(rtype), intent(  out) :: qwthl_sec   (shcol,nlevi             ) ! temp moisture covariance [K kg/kg]
//        real(rtype), intent(  out) :: wthl_sec    (shcol,nlevi             ) ! vertical heat flux [K m/s]
//        real(rtype), intent(  out) :: wqw_sec     (shcol,nlevi             ) ! vertical moisture flux [K m/s]
//        real(rtype), intent(  out) :: wtke_sec    (shcol,nlevi             ) ! vertical tke flux [m3/s3]
//        real(rtype), intent(  out) :: uw_sec      (shcol,nlevi             ) ! vertical zonal momentum flux [m2/s2]
//        real(rtype), intent(  out) :: vw_sec      (shcol,nlevi             ) ! vertical meridional momentum flux [m2/s2]
//        real(rtype), intent(  out) :: w3          (shcol,nlevi             ) ! third moment vertical velocity [m3/s3]
//        real(rtype), intent(  out) :: wqls_sec    (shcol,nlev              ) ! liquid water flux [kg/kg m/s]
//        real(rtype), intent(  out) :: brunt       (shcol,nlev              ) ! brunt vaisala frequency [s-1]
//        real(rtype), intent(  out) :: isotropy    (shcol,nlev              ) ! return to isotropic timescale [s]





  }



  // Returns saturation vapor pressure
  // TODO: Make this the same as P3's
  YAKL_INLINE static real saturation_vapor_pressure(real temp) {
    real tc = temp - 273.15;
    return 610.94 * exp( 17.625*tc / (243.04+tc) );
  }



  YAKL_INLINE static real latent_heat_condensation(real temp) {
    real tc = temp - 273.15;
    return (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*1000;
  }



  YAKL_INLINE static real cp_moist(real rho_d, real rho_v, real rho_c, real cp_d, real cp_v, real cp_l) {
    // For the moist specific heat, ignore other species than water vapor and cloud droplets
    real rho = rho_d + rho_v + rho_c;
    return rho_d / rho * cp_d  +  rho_v / rho * cp_v  +  rho_c / rho * cp_l;
  }



  // Compute an instantaneous adjustment of sub or super saturation
  YAKL_INLINE static void compute_adjusted_state(real rho, real rho_d , real &rho_v , real &rho_c , real &temp,
                                                 real R_v , real cp_d , real cp_v , real cp_l) {
    // Define a tolerance for convergence
    real tol = 1.e-6;

    // Saturation vapor pressure at this temperature
    real svp = saturation_vapor_pressure( temp );

    // Vapor pressure at this temperature
    real pv = rho_v * R_v * temp;

    // If we're super-saturated, we need to condense until saturation is reached
    if        (pv > svp) {
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to condense
      real cond1  = 0;     // Minimum amount we can condense out
      real cond2 = rho_v;  // Maximum amount we can condense out

      bool keep_iterating = true;
      while (keep_iterating) {
        real rho_cond = (cond1 + cond2) / 2;                    // How much water vapor to condense for this iteration
        real rv_loc = max( 0._fp , rho_v - rho_cond );          // New vapor density
        real rc_loc = max( 0._fp , rho_c + rho_cond );          // New cloud liquid density
        real Lv = latent_heat_condensation(temp);               // Compute latent heat of condensation
        real cp = cp_moist(rho_d,rv_loc,rc_loc,cp_d,cp_v,cp_l); // New moist specific heat at constant pressure
        real temp_loc = temp + rho_cond*Lv/(rho*cp);            // New temperature after condensation
        real svp_loc = saturation_vapor_pressure(temp_loc);     // New saturation vapor pressure after condensation
        real pv_loc = rv_loc * R_v * temp_loc;                  // New vapor pressure after condensation
        // If we're supersaturated still, we need to condense out more water vapor
        // otherwise, we need to condense out less water vapor
        if (pv_loc > svp_loc) {
          cond1 = rho_cond;
        } else {
          cond2 = rho_cond;
        }
        // If we've converged, then we can stop iterating
        if (abs(cond2-cond1) <= tol) {
          rho_v = rv_loc;
          rho_c = rc_loc;
          temp  = temp_loc;
          keep_iterating = false;
        }
      }

    // If we are unsaturated and have cloud liquid
    } else if (pv < svp && rho_c > 0) {
      // If there's cloud, evaporate enough to achieve saturation
      // or all of it if there isn't enough to reach saturation
      ////////////////////////////////////////////////////////
      // Bisection method
      ////////////////////////////////////////////////////////
      // Set bounds on how much mass to evaporate
      real evap1 = 0;     // minimum amount we can evaporate
      real evap2 = rho_c; // maximum amount we can evaporate

      bool keep_iterating = true;
      while (keep_iterating) {
        real rho_evap = (evap1 + evap2) / 2;                    // How much water vapor to evapense
        real rv_loc = max( 0._fp , rho_v + rho_evap );          // New vapor density
        real rc_loc = max( 0._fp , rho_c - rho_evap );          // New cloud liquid density
        real Lv = latent_heat_condensation(temp);               // Compute latent heat of condensation for water
        real cp = cp_moist(rho_d,rv_loc,rc_loc,cp_d,cp_v,cp_l); // New moist specific heat
        real temp_loc = temp - rho_evap*Lv/(rho*cp);            // New temperature after evaporation
        real svp_loc = saturation_vapor_pressure(temp_loc);     // New saturation vapor pressure after evaporation
        real pv_loc = rv_loc * R_v * temp_loc;                  // New vapor pressure after evaporation
        // If we're unsaturated still, we need to evaporate out more water vapor
        // otherwise, we need to evaporate out less water vapor
        if (pv_loc < svp_loc) {
          evap1 = rho_evap;
        } else {
          evap2 = rho_evap;
        }
        // If we've converged, then we can stop iterating
        if (abs(evap2-evap1) <= tol) {
          rho_v = rv_loc;
          rho_c = rc_loc;
          temp  = temp_loc;
          keep_iterating = false;
        }
      }
    }
  }



  std::string micro_name() const {
    return "p3";
  }



};



