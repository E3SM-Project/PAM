
! MAIN QUESTIONS:
! * Does SHOC expect ground to be at nlev-1 or zero? (I'm prett
! * How do I back out water vapor and cloud liquid from SHOC's outputs?
! * Are there tunable parameters to increase or decrease overall dissipation?


subroutine shoc_init( &
         nlev, gravit, rair, rh2o, cpair, &
         zvir, latvap, latice, karman, &
         pref_mid, nbot_shoc, ntop_shoc)
  integer    , intent(in) :: nlev           !UNDERSTOOD number of levels
  real(rtype), intent(in) :: gravit         !UNDERSTOOD gravity
  real(rtype), intent(in) :: rair           !UNDERSTOOD dry air gas constant
  real(rtype), intent(in) :: rh2o           !UNDERSTOOD water vapor gas constant
  real(rtype), intent(in) :: cpair          !UNDERSTOOD specific heat of dry air
  real(rtype), intent(in) :: zvir           !UNDERSTOOD rh2o/rair - 1
  real(rtype), intent(in) :: latvap         !UNDERSTOOD latent heat of vaporization
  real(rtype), intent(in) :: latice         !UNDERSTOOD latent heat of fusion
  real(rtype), intent(in) :: karman         !UNDERSTOOD (STEAL FROM PHYSCONST) Von Karman's constant
  real(rtype), intent(in) :: pref_mid(nlev) !UNDERSTOOD reference pressures at midpoints
  integer    , intent(in) :: nbot_shoc      !UNDERSTOOD FROM SHOC_INTR ntop_shoc = 1   ;   Bottom level to which SHOC is applied
  integer    , intent(in) :: ntop_shoc      !UNDERSTOOD FROM SHOC_INTR nbot_shoc = pver;   Top level to which SHOC is applied


subroutine shoc_main ( &
     shcol, nlev, nlevi, dtime, nadv, &   ! Input
     host_dx, host_dy,thv, &              ! Input
     zt_grid,zi_grid,pres,presi,pdel,&    ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, & ! Input
     wtracer_sfc,num_qtracers,w_field, &  ! Input
     exner,phis, &                        ! Input
     host_dse, tke, thetal, qw, &         ! Input/Output
     u_wind, v_wind,qtracers,&            ! Input/Output
     wthv_sec,tkh,tk,&                    ! Input/Output
     shoc_ql,shoc_cldfrac,&               ! Input/Output
     pblh,&                               ! Output
     shoc_mix, isotropy,&                 ! Output (diagnostic)
     w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
     wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
     uw_sec, vw_sec, w3,&                 ! Output (diagnostic)
     wqls_sec, brunt, shoc_ql2 &          ! Output (diagnostic)
  integer    , intent(in   ) :: shcol                                  !UNDERSTOOD   number of SHOC columns in the array
  integer    , intent(in   ) :: nlev                                   !UNDERSTOOD   number of levels [-]
  integer    , intent(in   ) :: nlevi                                  !UNDERSTOOD   number of levels on interface grid [-]
  integer    , intent(in   ) :: num_qtracers                           !UNDERSTOOD   number of tracers [-]
  integer    , intent(in   ) :: nadv                                   !UNDERSTOOD   THIS WILL MOST LIKELY BE 1 ;  number of times to loop SHOC
  real(rtype), intent(in   ) :: dtime                                  !UNDERSTOOD   NAMELISTS HAVE DTIME = 150s;  SHOC timestep [s]
  real(rtype), intent(in   ) :: host_dx     (shcol                   ) !UNDERSTOOD   grid spacing of host model in x direction [m]
  real(rtype), intent(in   ) :: host_dy     (shcol                   ) !UNDERSTOOD   grid spacing of host model in y direction [m]
  real(rtype), intent(in   ) :: zt_grid     (shcol,nlev              ) !UNDERSTOOD   heights, for thermo grid [m]
  real(rtype), intent(in   ) :: zi_grid     (shcol,nlevi             ) !UNDERSTOOD   heights, for interface grid [m]
  real(rtype), intent(in   ) :: pres        (shcol,nlev              ) !UNDERSTOOD   pressure levels on thermo grid [Pa]
  real(rtype), intent(in   ) :: presi       (shcol,nlevi             ) !UNDERSTOOD   pressure levels on interface grid [Pa]
  real(rtype), intent(in   ) :: pdel        (shcol,nlev              ) !UNDERSTOOD   Differences in pressure levels [Pa]
  real(rtype), intent(in   ) :: thv         (shcol,nlev              ) !UNDERSTOOD   https://glossary.ametsoc.org/wiki/Virtual_potential_temperature ??  virtual potential temperature [K]
  real(rtype), intent(in   ) :: w_field     (shcol,nlev              ) !UNDERSTOOD   large scale vertical velocity [m/s]
  real(rtype), intent(in   ) :: wthl_sfc    (shcol                   ) !UNDERSTOOD   Surface sensible heat flux [K m/s]
  real(rtype), intent(in   ) :: wqw_sfc     (shcol                   ) !UNDERSTOOD   Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in   ) :: uw_sfc      (shcol                   ) !UNDERSTOOD   Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in   ) :: vw_sfc      (shcol                   ) !UNDERSTOOD   Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in   ) :: wtracer_sfc (shcol      ,num_qtracers) !QUESTION     MASS OR DRY MIXING RATIO? Surface flux for tracers [varies]
  real(rtype), intent(in   ) :: exner       (shcol,nlev              ) !UNDERSTOOD   Exner function [-]
  real(rtype), intent(in   ) :: phis        (shcol                   ) !UNDERSTOOD   Host model surface geopotential height
  real(rtype), intent(inout) :: host_dse    (shcol,nlev              ) !UNDERSTOOD   prognostic temp variable of host model;  dry static energy [J/kg];  dse = Cp*T + g*z + phis
  real(rtype), intent(inout) :: tke         (shcol,nlev              ) !UNDERSTOOD   turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: thetal      (shcol,nlev              ) !UNDERSTOOD   https://glossary.ametsoc.org/wiki/Liquid_water_potential_temperature ??  liquid water potential temperature [K]
  real(rtype), intent(inout) :: qw          (shcol,nlev              ) !QUESTION     INCLUDE CLOUD ICE AND / OR PRECIPITANTS?  total water mixing ratio [kg/kg]
  real(rtype), intent(inout) :: u_wind      (shcol,nlev              ) !UNDERSTOOD   u wind component [m/s]
  real(rtype), intent(inout) :: v_wind      (shcol,nlev              ) !UNDERSTOOD   v wind component [m/s]
  real(rtype), intent(inout) :: wthv_sec    (shcol,nlev              ) !UNDERSTOOD   https://glossary.ametsoc.org/wiki/Buoyancy_flux ??  buoyancy flux [K m/s]
  real(rtype), intent(inout) :: qtracers    (shcol,nlev ,num_qtracers) !QUESTION     MASS OR DRY MIXING RATIOS? tracers [varies]
  real(rtype), intent(inout) :: tk          (shcol,nlev              ) !QUESTION     HOW TO INITIALIZE?   eddy coefficient for momentum [m2/s]
  real(rtype), intent(inout) :: tkh         (shcol,nlev              ) !QUESTION     HOW TO INITIALIZE?   eddy coefficent for heat [m2/s]
  real(rtype), intent(inout) :: shoc_cldfrac(shcol,nlev              ) !QUESTION     HOW TO INITIALIZE?   Cloud fraction [-]
  real(rtype), intent(inout) :: shoc_ql     (shcol,nlev              ) !UNDERSTOOD   cloud liquid mixing ratio [kg/kg]
  real(rtype), intent(  out) :: pblh        (shcol                   ) ! planetary boundary layer depth [m]
  real(rtype), intent(  out) :: shoc_ql2    (shcol,nlev              ) ! cloud liquid mixing ratio variance [kg^2/kg^2]
  real(rtype), intent(  out) :: shoc_mix    (shcol,nlev              ) ! Turbulent length scale [m]
  real(rtype), intent(  out) :: w_sec       (shcol,nlev              ) ! vertical velocity variance [m2/s2]
  real(rtype), intent(  out) :: thl_sec     (shcol,nlevi             ) ! temperature variance [K^2]
  real(rtype), intent(  out) :: qw_sec      (shcol,nlevi             ) ! moisture variance [kg2/kg2]
  real(rtype), intent(  out) :: qwthl_sec   (shcol,nlevi             ) ! temp moisture covariance [K kg/kg]
  real(rtype), intent(  out) :: wthl_sec    (shcol,nlevi             ) ! vertical heat flux [K m/s]
  real(rtype), intent(  out) :: wqw_sec     (shcol,nlevi             ) ! vertical moisture flux [K m/s]
  real(rtype), intent(  out) :: wtke_sec    (shcol,nlevi             ) ! vertical tke flux [m3/s3]
  real(rtype), intent(  out) :: uw_sec      (shcol,nlevi             ) ! vertical zonal momentum flux [m2/s2]
  real(rtype), intent(  out) :: vw_sec      (shcol,nlevi             ) ! vertical meridional momentum flux [m2/s2]
  real(rtype), intent(  out) :: w3          (shcol,nlevi             ) ! third moment vertical velocity [m3/s3]
  real(rtype), intent(  out) :: wqls_sec    (shcol,nlev              ) ! liquid water flux [kg/kg m/s]
  real(rtype), intent(  out) :: brunt       (shcol,nlev              ) ! brunt vaisala frequency [s-1]
  real(rtype), intent(  out) :: isotropy    (shcol,nlev              ) ! return to isotropic timescale [s]

