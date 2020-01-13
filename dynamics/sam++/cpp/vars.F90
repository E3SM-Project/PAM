module vars
  use grid
  use params, only: crm_rknd, crm_iknd
  implicit none
  !--------------------------------------------------------------------
  ! prognostic variables:

  real(crm_rknd), bind(C) :: u(ncrms,dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm) ! x-wind
  real(crm_rknd), bind(C) :: v(ncrms,dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm) ! y-wind
  real(crm_rknd), bind(C) :: w(ncrms,dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz ) ! z-wind
  real(crm_rknd), bind(C) :: t(ncrms,dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:
  real(crm_rknd), bind(C) :: p   (ncrms,0:nx, (1-YES3D):ny, nzm)               ! perturbation pressure (from Poison eq)
  real(crm_rknd), bind(C) :: tabs(ncrms,nx, ny, nzm)                           ! temperature
  real(crm_rknd), bind(C) :: qv  (ncrms,nx, ny, nzm)                           ! water vapor
  real(crm_rknd), bind(C) :: qcl (ncrms,nx, ny, nzm)                           ! liquid water  (condensate)
  real(crm_rknd), bind(C) :: qpl (ncrms,nx, ny, nzm)                           ! liquid water  (precipitation)
  real(crm_rknd), bind(C) :: qci (ncrms,nx, ny, nzm)                           ! ice water  (condensate)
  real(crm_rknd), bind(C) :: qpi (ncrms,nx, ny, nzm)                           ! ice water  (precipitation)
  real(crm_rknd), bind(C) :: tke2(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! SGS TKE
  real(crm_rknd), bind(C) :: tk2 (ncrms,0:nxp1, (1-YES3D):nyp1, nzm)           ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables
  real(crm_rknd), bind(C) :: dudt(ncrms,nxp1, ny, nzm, 3)
  real(crm_rknd), bind(C) :: dvdt(ncrms,nx, nyp1, nzm, 3)
  real(crm_rknd), bind(C) :: dwdt(ncrms,nx, ny  , nz,  3)

  !----------------------------------------------------------------
  ! Temporary storage array:

  real(crm_rknd), bind(C) :: misc(ncrms,nx, ny, nz)
  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:
  real(crm_rknd), bind(C) :: fluxbu  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fluxbv  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fluxbt  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fluxbq  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fluxtu  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fluxtv  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fluxtt  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fluxtq  (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: fzero   (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: precsfc (ncrms,nx,ny) ! surface precip. rate
  real(crm_rknd), bind(C) :: precssfc(ncrms,nx,ny) ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles
  real(crm_rknd), bind(C) :: t0   (ncrms,nzm)
  real(crm_rknd), bind(C) :: q0   (ncrms,nzm)
  real(crm_rknd), bind(C) :: qv0  (ncrms,nzm)
  real(crm_rknd), bind(C) :: tabs0(ncrms,nzm)
  real(crm_rknd), bind(C) :: tv0  (ncrms,nzm)
  real(crm_rknd), bind(C) :: u0   (ncrms,nzm)
  real(crm_rknd), bind(C) :: v0   (ncrms,nzm)
  real(crm_rknd), bind(C) :: tg0  (ncrms,nzm)
  real(crm_rknd), bind(C) :: qg0  (ncrms,nzm)
  real(crm_rknd), bind(C) :: ug0  (ncrms,nzm)
  real(crm_rknd), bind(C) :: vg0  (ncrms,nzm)
  real(crm_rknd), bind(C) :: p0   (ncrms,nzm)
  real(crm_rknd), bind(C) :: tke0 (ncrms,nzm)
  real(crm_rknd), bind(C) :: t01  (ncrms,nzm)
  real(crm_rknd), bind(C) :: q01  (ncrms,nzm)
  real(crm_rknd), bind(C) :: qp0  (ncrms,nzm)
  real(crm_rknd), bind(C) :: qn0  (ncrms,nzm)

  !-----------------------------------------------------------------
  ! reference vertical profiles:
  real(crm_rknd), bind(C) :: prespot(ncrms,nzm)  ! (1000./pres)**R/cp
  real(crm_rknd), bind(C) :: rho    (ncrms,nzm)   ! air density at pressure levels,kg/m3
  real(crm_rknd), bind(C) :: rhow   (ncrms,nz )   ! air density at vertical velocity levels,kg/m3
  real(crm_rknd), bind(C) :: bet    (ncrms,nzm)   ! = ggr/tv0
  real(crm_rknd), bind(C) :: gamaz  (ncrms,nzm) ! ggr/cp*z
  real(crm_rknd), bind(C) :: wsub   (ncrms,nz )   ! Large-scale subsidence velocity,m/s
  real(crm_rknd), bind(C) :: qtend  (ncrms,nzm) ! Large-scale tendency for total water
  real(crm_rknd), bind(C) :: ttend  (ncrms,nzm) ! Large-scale tendency for temp.
  real(crm_rknd), bind(C) :: utend  (ncrms,nzm) ! Large-scale tendency for u
  real(crm_rknd), bind(C) :: vtend  (ncrms,nzm) ! Large-scale tendency for v

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------
  !  Horizontally varying stuff (as a function of xy)
  real(crm_rknd), bind(C) :: sstxy    (ncrms,0:nx,(1-YES3D):ny) !  surface temperature xy-distribution
  real(crm_rknd), bind(C) :: fcory    (ncrms,0:ny)                 !  Coriolis parameter xy-distribution
  real(crm_rknd), bind(C) :: fcorzy   (ncrms,ny)                   !  z-Coriolis parameter xy-distribution
  real(crm_rknd), bind(C) :: latitude (ncrms,nx,ny)                  ! latitude (degrees,:)
  real(crm_rknd), bind(C) :: longitude(ncrms,nx,ny)                  ! longitude(degrees,:)
  real(crm_rknd), bind(C) :: prec_xy  (ncrms,nx,ny)             ! mean precip. rate for outout
  real(crm_rknd), bind(C) :: pw_xy    (ncrms,nx,ny)               ! precipitable water
  real(crm_rknd), bind(C) :: cw_xy    (ncrms,nx,ny)               ! cloud water path
  real(crm_rknd), bind(C) :: iw_xy    (ncrms,nx,ny)               ! ice water path
  real(crm_rknd), bind(C) :: cld_xy   (ncrms,nx,ny)               ! cloud frequency
  real(crm_rknd), bind(C) :: u200_xy  (ncrms,nx,ny)             ! u-wind at 200 mb
  real(crm_rknd), bind(C) :: usfc_xy  (ncrms,nx,ny)             ! u-wind at at the surface
  real(crm_rknd), bind(C) :: v200_xy  (ncrms,nx,ny)             ! v-wind at 200 mb
  real(crm_rknd), bind(C) :: vsfc_xy  (ncrms,nx,ny)             ! v-wind at the surface
  real(crm_rknd), bind(C) :: w500_xy  (ncrms,nx,ny)             ! w at 500 mb

  !----------------------------------------------------------------------
  ! Vertical profiles of quantities sampled for statitistics purposes:
  real(crm_rknd), bind(C) :: w_max    (ncrms)
  real(crm_rknd), bind(C) :: u_max    (ncrms)
  real(crm_rknd), bind(C) :: twsb     (ncrms,nz)
  real(crm_rknd), bind(C) :: precflux (ncrms,nz)
  real(crm_rknd), bind(C) :: uwle     (ncrms,nz)
  real(crm_rknd), bind(C) :: uwsb     (ncrms,nz)
  real(crm_rknd), bind(C) :: vwle     (ncrms,nz)
  real(crm_rknd), bind(C) :: vwsb     (ncrms,nz)
  real(crm_rknd), bind(C) :: tkelediss(ncrms,nz)
  real(crm_rknd), bind(C) :: tdiff    (ncrms,nz)
  real(crm_rknd), bind(C) :: tlat     (ncrms,nz)
  real(crm_rknd), bind(C) :: tlatqi   (ncrms,nz)
  real(crm_rknd), bind(C) :: qifall   (ncrms,nz)
  real(crm_rknd), bind(C) :: qpfall   (ncrms,nz)

  ! energy conservation diagnostics:
  real(8), bind(C) :: total_water_evap(ncrms)
  real(8), bind(C) :: total_water_prec(ncrms)

  real(crm_rknd), bind(C) :: CF3D(ncrms,1:nx, 1:ny, 1:nzm)  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd), bind(C) :: u850_xy(ncrms,nx,ny) ! zonal velocity at 850 mb
  real(crm_rknd), bind(C) :: v850_xy(ncrms,nx,ny) ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd), bind(C) :: psfc_xy(ncrms,nx,ny) ! pressure (in millibar) at lowest grid point

  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd), bind(C) :: swvp_xy(ncrms,nx,ny)  ! saturated water vapor path (wrt water)

  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd), bind(C) :: cloudtopheight(ncrms,nx,ny)
  real(crm_rknd), bind(C) :: echotopheight (ncrms,nx,ny)
  real(crm_rknd), bind(C) :: cloudtoptemp  (ncrms,nx,ny)

  ! END UW ADDITIONS
  !===========================================================================


contains


  subroutine allocate_vars()
    implicit none
    real(crm_rknd) :: zero

    zero = 0

    u = zero
    v = zero
    w = zero
    t = zero
    p = zero
    tabs = zero
    qv = zero
    qcl = zero
    qpl = zero
    qci = zero
    qpi = zero
    tke2 = zero
    tk2 = zero
    dudt = zero
    dvdt = zero
    dwdt = zero
    misc = zero
    fluxbu = zero
    fluxbv = zero
    fluxbt = zero
    fluxbq = zero
    fluxtu = zero
    fluxtv = zero
    fluxtt = zero
    fluxtq = zero
    fzero = zero
    precsfc = zero
    precssfc = zero
    t0 = zero
    q0 = zero
    qv0 = zero
    tabs0 = zero
    tv0 = zero
    u0 = zero
    v0 = zero
    tg0 = zero
    qg0 = zero
    ug0 = zero
    vg0 = zero
    p0 = zero
    tke0 = zero
    t01 = zero
    q01 = zero
    qp0 = zero
    qn0 = zero
    prespot = zero
    rho = zero
    rhow = zero
    bet = zero
    gamaz = zero
    wsub = zero
    qtend = zero
    ttend = zero
    utend = zero
    vtend = zero
    sstxy = zero
    fcory = zero
    fcorzy = zero
    latitude = zero
    longitude = zero
    prec_xy = zero
    pw_xy = zero
    cw_xy = zero
    iw_xy = zero
    cld_xy = zero
    u200_xy = zero
    usfc_xy = zero
    v200_xy = zero
    vsfc_xy = zero
    w500_xy = zero
    twsb = zero
    precflux = zero
    uwle = zero
    uwsb = zero
    vwle = zero
    vwsb = zero
    tkelediss = zero
    tdiff = zero
    tlat = zero
    tlatqi = zero
    qifall = zero
    qpfall = zero
    CF3D = 1.
    u850_xy = zero
    v850_xy = zero
    psfc_xy = zero
    swvp_xy = zero
    cloudtopheight = zero
    echotopheight = zero
    cloudtoptemp = zero
    u_max = zero
    w_max = zero
    total_water_evap = zero
    total_water_prec = zero
  end subroutine allocate_vars


  subroutine deallocate_vars()
    implicit none
end subroutine deallocate_vars


end module vars
