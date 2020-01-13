module params
  use iso_c_binding
  ! use grid, only: nzm

!#ifdef CRM
!  use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
!  shr_const_latice, shr_const_latsub, shr_const_rgas, &
!  shr_const_mwwv, shr_const_stebol, shr_const_tkfrz, &
!  shr_const_mwdair, shr_const_g, shr_const_karman, &
!  shr_const_rhofw
!#endif /*CRM*/

  implicit none
  integer, parameter :: r8 = c_double
  integer, parameter :: crm_rknd = c_double
  integer, parameter :: crm_iknd = c_int
  integer, parameter :: crm_lknd = c_bool

  !   Constants:

! #ifndef CRM
  real(crm_rknd), parameter :: cp    = 1004.          ! Specific heat of air, J/kg/K
  real(crm_rknd), parameter :: ggr   = 9.81           ! Gravity acceleration, m/s2
  real(crm_rknd), parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
  real(crm_rknd), parameter :: lfus  = 0.3336e+06   ! Latent heat of fusion, J/kg
  real(crm_rknd), parameter :: lsub  = 2.8440e+06     ! Latent heat of sublimation, J/kg
  real(crm_rknd), parameter :: rv    = 461.           ! Gas constant for water vapor, J/kg/K
  real(crm_rknd), parameter :: rgas  = 287.           ! Gas constant for dry air, J/kg/K
!#else
!  real(crm_rknd), parameter :: cp    = real( shr_const_cpdair ,crm_rknd)
!  real(crm_rknd), parameter :: ggr   = real( shr_const_g      ,crm_rknd)
!  real(crm_rknd), parameter :: lcond = real( shr_const_latvap ,crm_rknd)
!  real(crm_rknd), parameter :: lfus  = real( shr_const_latice ,crm_rknd)
!  real(crm_rknd), parameter :: lsub  = real( lcond + lfus     ,crm_rknd)
!  real(crm_rknd), parameter :: rgas  = real( shr_const_rdair  ,crm_rknd)
!  real(crm_rknd), parameter :: rv    = real( shr_const_rgas/shr_const_mwwv ,crm_rknd)
!#endif
  real(crm_rknd), parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
  real(crm_rknd), parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
  real(crm_rknd), parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

  real(crm_rknd), parameter :: fac_cond = lcond/cp
  real(crm_rknd), parameter :: fac_fus  = lfus/cp
  real(crm_rknd), parameter :: fac_sub  = lsub/cp

  real(crm_rknd), parameter ::  pi = 3.141592653589793

  integer(crm_iknd), parameter :: asyncid = 1

  !
  ! internally set parameters:

  real(crm_rknd)    , bind(C) :: epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
  logical(crm_lknd) , bind(C) :: dosubsidence = .false.

  !----------------------------------------------
  ! Parameters set by PARAMETERS namelist:
  ! Initialized to default values.
  !----------------------------------------------

  real(crm_rknd)    , bind(C) :: ug = 0.        ! Velocity of the Domain's drift in x direction
  real(crm_rknd)    , bind(C) :: vg = 0.        ! Velocity of the Domain's drift in y direction

  logical(crm_lknd) , bind(C) :: les =.false.    ! flag for Large-Eddy Simulation
  logical(crm_lknd) , bind(C) :: sfc_flx_fxd =.false. ! surface sensible flux is fixed
  logical(crm_lknd) , bind(C) :: sfc_tau_fxd =.false.! surface drag is fixed

  logical(crm_lknd) , bind(C) :: dodamping = .false.
  logical(crm_lknd) , bind(C) :: docloud = .false.
  logical(crm_lknd) , bind(C) :: docam_sfc_fluxes = .false.   ! Apply the surface fluxes within CAM
  logical(crm_lknd) , bind(C) :: doprecip = .false.
  logical(crm_lknd) , bind(C) :: dosgs = .false.
  logical(crm_lknd) , bind(C) :: docoriolis = .false.
  logical(crm_lknd) , bind(C) :: dosurface = .false.
  logical(crm_lknd) , bind(C) :: dowallx = .false.
  logical(crm_lknd) , bind(C) :: dowally = .false.
  logical(crm_lknd) , bind(C) :: docolumn = .false.
  logical(crm_lknd) , bind(C) :: dosmoke = .false.

  real(crm_rknd)   , bind(C) :: fcor      (NCRMS)       ! Coriolis parameter
  real(crm_rknd)   , bind(C) :: fcorz     (NCRMS)      ! Vertical Coriolis parameter
  real(crm_rknd)   , bind(C) :: longitude0(NCRMS)    ! latitude of the domain's center
  real(crm_rknd)   , bind(C) :: latitude0 (NCRMS)    ! longitude of the domain's center
  real(crm_rknd)   , bind(C) :: z0        (NCRMS)            ! roughness length
  logical(crm_lknd), bind(C) :: ocean     (NCRMS)           ! flag indicating that surface is water
  logical(crm_lknd), bind(C) :: land      (NCRMS)            ! flag indicating that surface is land
  real(crm_rknd)   , bind(C) :: uhl       (NCRMS)    ! current large-scale velocity in x near sfc
  real(crm_rknd)   , bind(C) :: vhl       (NCRMS)    ! current large-scale velocity in y near sfc
  real(crm_rknd)   , bind(C) :: taux0     (NCRMS)    ! surface stress in x, m2/s2
  real(crm_rknd)   , bind(C) :: tauy0     (NCRMS)    ! surface stress in y, m2/s2


contains

  
  subroutine allocate_params()
    implicit none
    fcor  = 0
    fcorz = 0
    longitude0 = 0
    latitude0  = 0
    z0 = 0.035
    ocean = .false.
    land = .false.
    uhl = 0
    vhl = 0
    taux0 = 0
    tauy0 = 0
  end subroutine allocate_params

  
  subroutine deallocate_params()
    implicit none
  end subroutine deallocate_params


end module params
