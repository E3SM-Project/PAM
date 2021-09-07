
module micro_mod
  implicit none
  real(8), parameter :: rhor = 1000.D0           ! Density of water, kg/m3
  real(8), parameter :: rhos = 100.D0            ! Density of snow, kg/m3
  real(8), parameter :: rhog = 400.D0            ! Density of graupel, kg/m3
  real(8), parameter :: tbgmin = 253.16D0        ! Minimum temperature for cloud water., K
  real(8), parameter :: tbgmax = 273.16D0        ! Maximum temperature for cloud ice, K
  real(8), parameter :: tprmin = 268.16D0        ! Minimum temperature for rain, K
  real(8), parameter :: tprmax = 283.16D0        ! Maximum temperature for snow+graupel, K
  real(8), parameter :: tgrmin = 223.16D0        ! Minimum temperature for snow, K
  real(8), parameter :: tgrmax = 283.16D0        ! Maximum temperature for graupel, K
  real(8), parameter :: a_rain = 842.D0          ! Coeff.for rain term vel
  real(8), parameter :: b_rain = 0.8D0           ! Fall speed exponent for rain
  real(8), parameter :: a_snow = 4.84D0          ! Coeff.for snow term vel
  real(8), parameter :: b_snow = 0.25D0          ! Fall speed exponent for snow
  real(8), parameter :: a_grau = 94.5D0          ! Lin (1983) (rhog=400)
  real(8), parameter :: b_grau = 0.5D0           ! Fall speed exponent for graupel
  real(8), parameter :: qcw0 = 1.D-3             ! Threshold for water autoconversion, g/g
  real(8), parameter :: qci0 = 1.D-4             ! Threshold for ice autoconversion, g/g
  real(8), parameter :: alphaelq = 1.D-3         ! autoconversion of cloud water rate coef
  real(8), parameter :: betaelq = 1.D-3          ! autoconversion of cloud ice rate coef
  real(8), parameter :: erccoef = 1.0D0          ! Rain/Cloud water collection efficiency
  real(8), parameter :: esccoef = 1.0D0          ! Snow/Cloud water collection efficiency
  real(8), parameter :: esicoef = 0.1D0          ! Snow/cloud ice collection efficiency
  real(8), parameter :: egccoef = 1.0D0          ! Graupel/Cloud water collection efficiency
  real(8), parameter :: egicoef = 0.1D0          ! Graupel/Cloud ice collection efficiency
  real(8), parameter :: nzeror = 8.D6            ! Intercept coeff. for rain
  real(8), parameter :: nzeros = 3.D6            ! Intersept coeff. for snow
  real(8), parameter :: nzerog = 4.D6            ! Intersept coeff. for graupel
  real(8), parameter :: qp_threshold = 1.D-8     ! minimal rain/snow water content
  real(8), parameter :: cp = 1004.D0             ! Specific heat of air, J/kg/K
  real(8), parameter :: ggr = 9.81D0             ! Gravity acceleration, m/s2
  real(8), parameter :: lcond = 2.5104D+06       ! Latent heat of condensation, J/kg
  real(8), parameter :: lfus = 0.3336D+06        ! Latent heat of fusion, J/kg
  real(8), parameter :: lsub = 2.8440D+06        ! Latent heat of sublimation, J/kg
  real(8), parameter :: rv = 461.D0              ! Gas constant for water vapor, J/kg/K
  real(8), parameter :: rgas = 287.D0            ! Gas constant for dry air, J/kg/K
  real(8), parameter :: diffelq = 2.21D-05       ! Diffusivity of water vapor, m2/s
  real(8), parameter :: therco = 2.40D-02        ! Thermal conductivity of air, J/m/s/K
  real(8), parameter :: muelq = 1.717D-05        ! Dynamic viscosity of air
  real(8), parameter :: fac_cond = lcond/cp      
  real(8), parameter :: fac_fus = lfus/cp        
  real(8), parameter :: fac_sub = lsub/cp        
  real(8), parameter :: pi = 3.141592653589793D0
  real(8), parameter :: a_bg = 1.D0/(tbgmax-tbgmin)
  real(8), parameter :: a_pr = 1.D0/(tprmax-tprmin)
  real(8), parameter :: a_gr = 1.D0/(tgrmax-tgrmin)


contains

  subroutine micro(tabs, pres, rho, q, t, gamaz, qp, qn, qpsrc, qpevp, qv, qcl, qci, qpl, qpi, dt, ncol, nz)
    use precip_init_mod
    use cloud_mod
    use precip_proc_mod
    use micro_diagnose_mod
    implicit none
    real(8), intent(inout) :: tabs   (ncol,nz) ! temperature
    real(8), intent(in   ) :: pres   (ncol,nz) ! pressure,mb at scalar levels
    real(8), intent(in   ) :: rho    (ncol,nz) ! air density at pressure levels,kg/m3 
    real(8), intent(inout) :: q      (ncol,nz) ! total nonprecipitating water
    real(8), intent(in   ) :: t      (ncol,nz) ! liquid/ice water static energy 
    real(8), intent(in   ) :: gamaz  (ncol,nz) ! grav/cp*z
    real(8), intent(inout) :: qp     (ncol,nz) ! total precipitating water
    real(8), intent(  out) :: qn     (ncol,nz) ! cloud condensate (liquid + ice)
    real(8), intent(  out) :: qpsrc  (ncol,nz) ! source of precipitation microphysical processes
    real(8), intent(  out) :: qpevp  (ncol,nz) ! sink of precipitating water due to evaporation
    real(8), intent(  out) :: qv     (ncol,nz) ! water vapor
    real(8), intent(  out) :: qcl    (ncol,nz) ! liquid water  (condensate)
    real(8), intent(  out) :: qci    (ncol,nz) ! ice water  (condensate)
    real(8), intent(  out) :: qpl    (ncol,nz) ! liquid water  (precipitation)
    real(8), intent(  out) :: qpi    (ncol,nz) ! ice water  (precipitation)
    real(8), intent(in   ) :: dt
    integer, intent(in   ) :: ncol, nz

    ! The following are computed by precip_init
    real(8) :: accrsi (ncol,nz) ! Undocumented
    real(8) :: accrsc (ncol,nz) ! Undocumented
    real(8) :: coefice(ncol,nz) ! Undocumented
    real(8) :: evaps1 (ncol,nz) ! Undocumented
    real(8) :: evaps2 (ncol,nz) ! Undocumented
    real(8) :: accrgi (ncol,nz) ! Undocumented
    real(8) :: accrgc (ncol,nz) ! Undocumented
    real(8) :: evapg1 (ncol,nz) ! Undocumented
    real(8) :: evapg2 (ncol,nz) ! Undocumented
    real(8) :: accrrc (ncol,nz) ! Undocumented
    real(8) :: evapr1 (ncol,nz) ! Undocumented
    real(8) :: evapr2 (ncol,nz) ! Undocumented
    real(8) :: gam3, gamr1, gamr2, gamr3, gams1, gams2, gams3, gamg1, gamg2, gamg3

    call precip_init(tabs, pres, rho, accrsi, accrsc, coefice, evaps1, evaps2, accrgi, accrgc, evapg1, evapg2,  &
                     accrrc, evapr1, evapr2, b_rain, b_snow, b_grau, a_grau, a_rain, a_snow, diffelq, egccoef,  &
                     egicoef, erccoef, esccoef, esicoef, lcond, lsub, muelq, nzerog, nzeror, nzeros, rhog,      &
                     rhor, rhos, rv, therco, gam3, gamr1, gamr2, gamr3, gams1, gams2, gams3, gamg1, gamg2,      &
                     gamg3, ncol, nz)
    ! real(8), intent(in   ) :: tabs0  (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: pres   (ncol,nz) ! pressure,mb at scalar levels
    ! real(8), intent(in   ) :: rho    (ncol,nz) ! air density at pressure levels,kg/m3 
    ! real(8), intent(  out) :: accrsi (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: accrsc (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: coefice(ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: evaps1 (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: evaps2 (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: accrgi (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: accrgc (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: evapg1 (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: evapg2 (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: accrrc (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: evapr1 (ncol,nz) ! Undocumented
    ! real(8), intent(  out) :: evapr2 (ncol,nz) ! Undocumented


    call cloud(q, tabs, t, gamaz, qp, pres, qn, tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, &
               tgrmax, tgrmin, ncol, nz)
    ! real(8), intent(inout) :: q    (ncol,nz)  ! total nonprecipitating water
    ! real(8), intent(  out) :: tabs (ncol,nz)  ! temperature
    ! real(8), intent(in   ) :: t    (ncol,nz)  ! liquid/ice water static energy 
    ! real(8), intent(in   ) :: gamaz(ncol,nz)  ! grav/cp*z
    ! real(8), intent(inout) :: qp   (ncol,nz)  ! total precipitating water
    ! real(8), intent(in   ) :: pres (ncol,nz)  ! pressure,mb at scalar levels
    ! real(8), intent(  out) :: qn   (ncol,nz)  ! cloud condensate (liquid + ice)

    call precip_proc(qpsrc, qpevp, qn, qp, tabs, coefice, accrrc, accrsc, accrsi, accrgc, accrgi, q, pres, &
                     evapr1, evapr2, evaps1, evaps2, evapg1, evapg2, a_bg, a_gr, a_pr, alphaelq, b_grau,   &
                     b_rain, b_snow, betaelq, dt, qci0, qcw0, qp_threshold, tbgmin, tgrmin, tprmin,       &
                     ncol, nz)
    ! real(8), intent(  out) :: qpsrc  (ncol,nz) ! source of precipitation microphysical processes
    ! real(8), intent(  out) :: qpevp  (ncol,nz) ! sink of precipitating water due to evaporation
    ! real(8), intent(inout) :: qn     (ncol,nz) ! cloud condensate (liquid + ice)
    ! real(8), intent(inout) :: qp     (ncol,nz) ! total precipitating water
    ! real(8), intent(in   ) :: tabs   (ncol,nz) ! temperature
    ! real(8), intent(in   ) :: coefice(ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: accrrc (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: accrsc (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: accrsi (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: accrgc (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: accrgi (ncol,nz) ! Undocumented
    ! real(8), intent(inout) :: q      (ncol,nz) ! total nonprecipitating water
    ! real(8), intent(in   ) :: pres   (ncol,nz) ! pressure,mb at scalar levels
    ! real(8), intent(in   ) :: evapr1 (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: evapr2 (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: evaps1 (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: evaps2 (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: evapg1 (ncol,nz) ! Undocumented
    ! real(8), intent(in   ) :: evapg2 (ncol,nz) ! Undocumented

    call micro_diagnose(qv, q, qn, tabs, qcl, qci, qpl, qpi, qp, a_bg, a_pr, tbgmin, tprmin, ncol, nz)
    ! real(8), intent(  out) :: qv  (ncol,nz) ! water vapor
    ! real(8), intent(in   ) :: q   (ncol,nz) ! total nonprecipitating water
    ! real(8), intent(in   ) :: qn  (ncol,nz) ! cloud condensate (liquid + ice)
    ! real(8), intent(in   ) :: tabs(ncol,nz) ! temperature
    ! real(8), intent(  out) :: qcl (ncol,nz) ! liquid water  (condensate)
    ! real(8), intent(  out) :: qci (ncol,nz) ! ice water  (condensate)
    ! real(8), intent(  out) :: qpl (ncol,nz) ! liquid water  (precipitation)
    ! real(8), intent(  out) :: qpi (ncol,nz) ! ice water  (precipitation)
    ! real(8), intent(in   ) :: qp  (ncol,nz) ! total precipitating water

  endsubroutine micro

endmodule micro_mod


