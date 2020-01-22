module vars
  use grid
  use params, only: crm_rknd
  use gator_mod, only: gator_allocate, gator_deallocate

  implicit none
  !--------------------------------------------------------------------
  ! prognostic variables:

  real(crm_rknd), pointer, contiguous :: u   (:,:,:,:) ! x-wind
  real(crm_rknd), pointer, contiguous :: v   (:,:,:,:) ! y-wind
  real(crm_rknd), pointer, contiguous :: w   (:,:,:,:) ! z-wind
  real(crm_rknd), pointer, contiguous :: t   (:,:,:,:) ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:

  real(crm_rknd), pointer, contiguous :: p       (:,:,:,:)     ! perturbation pressure (from Poison eq)
  real(crm_rknd), pointer, contiguous :: tabs    (:,:,:,:)                 ! temperature
  real(crm_rknd), pointer, contiguous :: qv      (:,:,:,:)                ! water vapor
  real(crm_rknd), pointer, contiguous :: qcl     (:,:,:,:)                ! liquid water  (condensate)
  real(crm_rknd), pointer, contiguous :: qpl     (:,:,:,:)                ! liquid water  (precipitation)
  real(crm_rknd), pointer, contiguous :: qci     (:,:,:,:)                ! ice water  (condensate)
  real(crm_rknd), pointer, contiguous :: qpi     (:,:,:,:)                ! ice water  (precipitation)
  real(crm_rknd), pointer, contiguous :: tke2    (:,:,:,:)   ! SGS TKE
  real(crm_rknd), pointer, contiguous :: tk2     (:,:,:,:) ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables

  real(crm_rknd), pointer, contiguous :: dudt   (:,:,:,:,:)
  real(crm_rknd), pointer, contiguous :: dvdt   (:,:,:,:,:)
  real(crm_rknd), pointer, contiguous :: dwdt   (:,:,:,:,:)

  !----------------------------------------------------------------
  ! Temporary storage array:

  real(crm_rknd), pointer, contiguous :: misc(:,:,:,:)
  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:

  real(crm_rknd), pointer, contiguous :: fluxbu  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fluxbv  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fluxbt  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fluxbq  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fluxtu  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fluxtv  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fluxtt  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fluxtq  (:,:,:)
  real(crm_rknd), pointer, contiguous :: fzero   (:,:,:)
  real(crm_rknd), pointer, contiguous :: precsfc (:,:,:) ! surface precip. rate
  real(crm_rknd), pointer, contiguous :: precssfc(:,:,:) ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles

  real(crm_rknd), pointer, contiguous :: t0   (:,:)
  real(crm_rknd), pointer, contiguous :: q0   (:,:)
  real(crm_rknd), pointer, contiguous :: qv0  (:,:)
  real(crm_rknd), pointer, contiguous :: tabs0(:,:)
  real(crm_rknd), pointer, contiguous :: tv0  (:,:)
  real(crm_rknd), pointer, contiguous :: u0   (:,:)
  real(crm_rknd), pointer, contiguous :: v0   (:,:)
  real(crm_rknd), pointer, contiguous :: tg0  (:,:)
  real(crm_rknd), pointer, contiguous :: qg0  (:,:)
  real(crm_rknd), pointer, contiguous :: ug0  (:,:)
  real(crm_rknd), pointer, contiguous :: vg0  (:,:)
  real(crm_rknd), pointer, contiguous :: p0   (:,:)
  real(crm_rknd), pointer, contiguous :: tke0 (:,:)
  real(crm_rknd), pointer, contiguous :: t01  (:,:)
  real(crm_rknd), pointer, contiguous :: q01  (:,:)
  real(crm_rknd), pointer, contiguous :: qp0  (:,:)
  real(crm_rknd), pointer, contiguous :: qn0  (:,:)

  !-----------------------------------------------------------------
  ! reference vertical profiles:
  real(crm_rknd), pointer, contiguous :: prespot(:,:)  ! (1000./pres)**R/cp
  real(crm_rknd), pointer, contiguous :: rho    (:,:)   ! air density at pressure levels,kg/m3
  real(crm_rknd), pointer, contiguous :: rhow   (:,:)   ! air density at vertical velocity levels,kg/m3
  real(crm_rknd), pointer, contiguous :: bet    (:,:)   ! = ggr/tv0
  real(crm_rknd), pointer, contiguous :: gamaz  (:,:) ! ggr/cp*z
  real(crm_rknd), pointer, contiguous :: wsub   (:,:)   ! Large-scale subsidence velocity,m/s
  real(crm_rknd), pointer, contiguous :: qtend  (:,:) ! Large-scale tendency for total water
  real(crm_rknd), pointer, contiguous :: ttend  (:,:) ! Large-scale tendency for temp.
  real(crm_rknd), pointer, contiguous :: utend  (:,:) ! Large-scale tendency for u
  real(crm_rknd), pointer, contiguous :: vtend  (:,:) ! Large-scale tendency for v

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------
  !  Horizontally varying stuff (as a function of xy)
  !
  real(crm_rknd), pointer, contiguous :: sstxy    (:,:,:) !  surface temperature xy-distribution
  real(crm_rknd), pointer, contiguous :: fcory    (:,:)      !  Coriolis parameter xy-distribution
  real(crm_rknd), pointer, contiguous :: fcorzy   (:,:)      !  z-Coriolis parameter xy-distribution
  real(crm_rknd), pointer, contiguous :: latitude (:,:,:)      ! latitude (degrees,:)
  real(crm_rknd), pointer, contiguous :: longitude(:,:,:)      ! longitude(degrees,:)
  real(crm_rknd), pointer, contiguous :: prec_xy  (:,:,:) ! mean precip. rate for outout
  real(crm_rknd), pointer, contiguous :: pw_xy    (:,:,:)   ! precipitable water
  real(crm_rknd), pointer, contiguous :: cw_xy    (:,:,:)   ! cloud water path
  real(crm_rknd), pointer, contiguous :: iw_xy    (:,:,:)   ! ice water path
  real(crm_rknd), pointer, contiguous :: cld_xy   (:,:,:)   ! cloud frequency
  real(crm_rknd), pointer, contiguous :: u200_xy  (:,:,:) ! u-wind at 200 mb
  real(crm_rknd), pointer, contiguous :: usfc_xy  (:,:,:) ! u-wind at at the surface
  real(crm_rknd), pointer, contiguous :: v200_xy  (:,:,:) ! v-wind at 200 mb
  real(crm_rknd), pointer, contiguous :: vsfc_xy  (:,:,:) ! v-wind at the surface
  real(crm_rknd), pointer, contiguous :: w500_xy  (:,:,:) ! w at 500 mb

  !----------------------------------------------------------------------
  ! Vertical profiles of quantities sampled for statitistics purposes:

  real(crm_rknd), pointer, contiguous :: w_max(:)
  real(crm_rknd), pointer, contiguous :: u_max(:)

  real(crm_rknd), pointer, contiguous :: twsb(:,:)
  real(crm_rknd), pointer, contiguous :: precflux(:,:)
  real(crm_rknd), pointer, contiguous :: uwle(:,:)
  real(crm_rknd), pointer, contiguous :: uwsb(:,:)
  real(crm_rknd), pointer, contiguous :: vwle(:,:)
  real(crm_rknd), pointer, contiguous :: vwsb(:,:)
  real(crm_rknd), pointer, contiguous :: tkelediss(:,:)
  real(crm_rknd), pointer, contiguous :: tdiff(:,:)
  real(crm_rknd), pointer, contiguous :: tlat(:,:)
  real(crm_rknd), pointer, contiguous :: tlatqi(:,:)
  real(crm_rknd), pointer, contiguous :: qifall(:,:)
  real(crm_rknd), pointer, contiguous :: qpfall(:,:)

  ! energy conservation diagnostics:
  real(8), pointer, contiguous :: total_water_evap(:)
  real(8), pointer, contiguous :: total_water_prec(:)

  real(crm_rknd), pointer, contiguous :: CF3D(:,:,:,:)  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd), pointer, contiguous :: u850_xy(:,:,:) ! zonal velocity at 850 mb
  real(crm_rknd), pointer, contiguous :: v850_xy(:,:,:) ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd), pointer, contiguous :: psfc_xy(:,:,:) ! pressure (in millibar) at lowest grid point

  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd), pointer, contiguous :: swvp_xy(:,:,:)  ! saturated water vapor path (wrt water)

  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd), pointer, contiguous :: cloudtopheight(:,:,:)
  real(crm_rknd), pointer, contiguous :: echotopheight (:,:,:)
  real(crm_rknd), pointer, contiguous :: cloudtoptemp  (:,:,:)

  ! END UW ADDITIONS
  !===========================================================================


contains


  subroutine allocate_vars()
    implicit none
    real(crm_rknd) :: zero
    call gator_allocate( u                 , (/ncrms,dimx2_u-dimx1_u+1,dimy2_u-dimy1_u+1,nzm/) , (/1,dimx1_u,dimy1_u,1/) )
    call gator_allocate( v                 , (/ncrms,dimx2_v-dimx1_v+1,dimy2_v-dimy1_v+1,nzm/) , (/1,dimx1_v,dimy1_v,1/) )
    call gator_allocate( w                 , (/ncrms,dimx2_w-dimx1_w+1,dimy2_w-dimy1_w+1,nz /) , (/1,dimx1_w,dimy1_w,1/) )
    call gator_allocate( t                 , (/ncrms,dimx2_s-dimx1_s+1,dimy2_s-dimy1_s+1,nzm/) , (/1,dimx1_s,dimy1_s,1/) )
    call gator_allocate( tke2              , (/ncrms,dimx2_s-dimx1_s+1,dimy2_s-dimy1_s+1,nzm/) , (/1,dimx1_s,dimy1_s,1/) )
    call gator_allocate( p                 , (/ncrms,nxp1,ny  -(1-YES3D)+1,nzm/)               , (/1,0,1-YES3D,1/) )
    call gator_allocate( tk2               , (/ncrms,nxp2,nyp1-(1-YES3D)+1,nzm/)               , (/1,0,1-YES3D,1/) )
    call gator_allocate( sstxy             , (/ncrms,nxp1,ny-(1-YES3D)+1/)                     , (/1,0,1-YES3D/) )
    call gator_allocate( fcory             , (/ncrms,nyp1/)                                    , (/1,0/) )

    call gator_allocate( tabs              , (/ncrms,nx, ny, nzm/)      )         
    call gator_allocate( qv                , (/ncrms,nx, ny, nzm/)      )
    call gator_allocate( qcl               , (/ncrms,nx, ny, nzm/)      )
    call gator_allocate( qpl               , (/ncrms,nx, ny, nzm/)      )
    call gator_allocate( qci               , (/ncrms,nx, ny, nzm/)      )
    call gator_allocate( qpi               , (/ncrms,nx, ny, nzm/)      )
    call gator_allocate( dudt              , (/ncrms,nxp1, ny, nzm, 3/) )
    call gator_allocate( dvdt              , (/ncrms,nx, nyp1, nzm, 3/) )
    call gator_allocate( dwdt              , (/ncrms,nx, ny  , nz,  3/) )
    call gator_allocate( misc              , (/ncrms,nx, ny, nz/)       )
    call gator_allocate( fluxbu            , (/ncrms,nx,ny/)            )
    call gator_allocate( fluxbv            , (/ncrms,nx,ny/)            )
    call gator_allocate( fluxbt            , (/ncrms,nx,ny/)            )
    call gator_allocate( fluxbq            , (/ncrms,nx,ny/)            )
    call gator_allocate( fluxtu            , (/ncrms,nx,ny/)            )
    call gator_allocate( fluxtv            , (/ncrms,nx,ny/)            )
    call gator_allocate( fluxtt            , (/ncrms,nx,ny/)            )
    call gator_allocate( fluxtq            , (/ncrms,nx,ny/)            )
    call gator_allocate( fzero             , (/ncrms,nx,ny/)            )
    call gator_allocate( precsfc           , (/ncrms,nx,ny/)            )
    call gator_allocate( precssfc          , (/ncrms,nx,ny/)            )
    call gator_allocate( t0                , (/ncrms,nzm/)              )
    call gator_allocate( q0                , (/ncrms,nzm/)              )
    call gator_allocate( qv0               , (/ncrms,nzm/)              )
    call gator_allocate( tabs0             , (/ncrms,nzm/)              )
    call gator_allocate( tv0               , (/ncrms,nzm/)              )
    call gator_allocate( u0                , (/ncrms,nzm/)              )
    call gator_allocate( v0                , (/ncrms,nzm/)              )
    call gator_allocate( tg0               , (/ncrms,nzm/)              )
    call gator_allocate( qg0               , (/ncrms,nzm/)              )
    call gator_allocate( ug0               , (/ncrms,nzm/)              )
    call gator_allocate( vg0               , (/ncrms,nzm/)              )
    call gator_allocate( p0                , (/ncrms,nzm/)              )
    call gator_allocate( tke0              , (/ncrms,nzm/)              )
    call gator_allocate( t01               , (/ncrms,nzm/)              )
    call gator_allocate( q01               , (/ncrms,nzm/)              )
    call gator_allocate( qp0               , (/ncrms,nzm/)              )
    call gator_allocate( qn0               , (/ncrms,nzm/)              )
    call gator_allocate( prespot           , (/ncrms,nzm/)              )
    call gator_allocate( rho               , (/ncrms,nzm/)              )
    call gator_allocate( rhow              , (/ncrms,nz /)              )
    call gator_allocate( bet               , (/ncrms,nzm/)              )
    call gator_allocate( gamaz             , (/ncrms,nzm/)              )
    call gator_allocate( wsub              , (/ncrms,nz /)              )
    call gator_allocate( qtend             , (/ncrms,nzm/)              )
    call gator_allocate( ttend             , (/ncrms,nzm/)              )
    call gator_allocate( utend             , (/ncrms,nzm/)              )
    call gator_allocate( vtend             , (/ncrms,nzm/)              )
    call gator_allocate( fcorzy            , (/ncrms,ny/)               )
    call gator_allocate( latitude          , (/ncrms,nx,ny/)            )
    call gator_allocate( longitude         , (/ncrms,nx,ny/)            )
    call gator_allocate( prec_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( pw_xy             , (/ncrms,nx,ny/)            )
    call gator_allocate( cw_xy             , (/ncrms,nx,ny/)            )
    call gator_allocate( iw_xy             , (/ncrms,nx,ny/)            )
    call gator_allocate( cld_xy            , (/ncrms,nx,ny/)            )
    call gator_allocate( u200_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( usfc_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( v200_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( vsfc_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( w500_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( twsb              , (/ncrms,nz/)               )
    call gator_allocate( precflux          , (/ncrms,nz/)               )
    call gator_allocate( uwle              , (/ncrms,nz/)               )
    call gator_allocate( uwsb              , (/ncrms,nz/)               )
    call gator_allocate( vwle              , (/ncrms,nz/)               )
    call gator_allocate( vwsb              , (/ncrms,nz/)               )
    call gator_allocate( tkelediss         , (/ncrms,nz/)               )
    call gator_allocate( tdiff             , (/ncrms,nz/)               )
    call gator_allocate( tlat              , (/ncrms,nz/)               )
    call gator_allocate( tlatqi            , (/ncrms,nz/)               )
    call gator_allocate( qifall            , (/ncrms,nz/)               )
    call gator_allocate( qpfall            , (/ncrms,nz/)               )
    call gator_allocate( cf3d              , (/ncrms,nx,ny,nzm/)        )
    call gator_allocate( u850_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( v850_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( psfc_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( swvp_xy           , (/ncrms,nx,ny/)            )
    call gator_allocate( cloudtopheight    , (/ncrms,nx,ny/)            )
    call gator_allocate( echotopheight     , (/ncrms,nx,ny/)            )
    call gator_allocate( cloudtoptemp      , (/ncrms,nx,ny/)            )
    call gator_allocate( u_max             , (/ncrms/)                  )
    call gator_allocate( w_max             , (/ncrms/)                  )
    call gator_allocate( total_water_evap  , (/ncrms/)                  )
    call gator_allocate( total_water_prec  , (/ncrms/)                  )

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
    call gator_deallocate( u )
    call gator_deallocate( v )
    call gator_deallocate( w )
    call gator_deallocate( t )
    call gator_deallocate( p )
    call gator_deallocate( tabs )
    call gator_deallocate( qv )
    call gator_deallocate( qcl )
    call gator_deallocate( qpl )
    call gator_deallocate( qci )
    call gator_deallocate( qpi )
    call gator_deallocate( tke2 )
    call gator_deallocate( tk2 )
    call gator_deallocate( dudt )
    call gator_deallocate( dvdt )
    call gator_deallocate( dwdt )
    call gator_deallocate( misc )
    call gator_deallocate( fluxbu )
    call gator_deallocate( fluxbv )
    call gator_deallocate( fluxbt )
    call gator_deallocate( fluxbq )
    call gator_deallocate( fluxtu )
    call gator_deallocate( fluxtv )
    call gator_deallocate( fluxtt )
    call gator_deallocate( fluxtq )
    call gator_deallocate( fzero )
    call gator_deallocate( precsfc )
    call gator_deallocate( precssfc )
    call gator_deallocate( t0 )
    call gator_deallocate( q0 )
    call gator_deallocate( qv0 )
    call gator_deallocate( tabs0 )
    call gator_deallocate( tv0 )
    call gator_deallocate( u0 )
    call gator_deallocate( v0 )
    call gator_deallocate( tg0 )
    call gator_deallocate( qg0 )
    call gator_deallocate( ug0 )
    call gator_deallocate( vg0 )
    call gator_deallocate( p0 )
    call gator_deallocate( tke0 )
    call gator_deallocate( t01 )
    call gator_deallocate( q01 )
    call gator_deallocate( qp0 )
    call gator_deallocate( qn0 )
    call gator_deallocate( prespot )
    call gator_deallocate( rho )
    call gator_deallocate( rhow )
    call gator_deallocate( bet )
    call gator_deallocate( gamaz )
    call gator_deallocate( wsub )
    call gator_deallocate( qtend )
    call gator_deallocate( ttend )
    call gator_deallocate( utend )
    call gator_deallocate( vtend )
    call gator_deallocate( sstxy )
    call gator_deallocate( fcory )
    call gator_deallocate( fcorzy )
    call gator_deallocate( latitude )
    call gator_deallocate( longitude )
    call gator_deallocate( prec_xy )
    call gator_deallocate( pw_xy )
    call gator_deallocate( cw_xy )
    call gator_deallocate( iw_xy )
    call gator_deallocate( cld_xy )
    call gator_deallocate( u200_xy )
    call gator_deallocate( usfc_xy )
    call gator_deallocate( v200_xy )
    call gator_deallocate( vsfc_xy )
    call gator_deallocate( w500_xy )
    call gator_deallocate( twsb )
    call gator_deallocate( precflux )
    call gator_deallocate( uwle )
    call gator_deallocate( uwsb )
    call gator_deallocate( vwle )
    call gator_deallocate( vwsb )
    call gator_deallocate( tkelediss )
    call gator_deallocate( tdiff )
    call gator_deallocate( tlat )
    call gator_deallocate( tlatqi )
    call gator_deallocate( qifall )
    call gator_deallocate( qpfall )
    call gator_deallocate( CF3D )
    call gator_deallocate( u850_xy )
    call gator_deallocate( v850_xy )
    call gator_deallocate( psfc_xy )
    call gator_deallocate( swvp_xy )
    call gator_deallocate( cloudtopheight )
    call gator_deallocate( echotopheight )
    call gator_deallocate( cloudtoptemp )
    call gator_deallocate( u_max )
    call gator_deallocate( w_max )
    call gator_deallocate( total_water_evap )
    call gator_deallocate( total_water_prec )
end subroutine deallocate_vars


end module vars
