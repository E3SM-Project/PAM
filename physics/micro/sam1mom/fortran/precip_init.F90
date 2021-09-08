module precip_init_mod
  implicit none

contains

  subroutine precip_init(tabs0, pres, rho, accrsi, accrsc, coefice, evaps1, evaps2, accrgi, accrgc, evapg1, evapg2, &
                         accrrc, evapr1, evapr2, b_rain, b_snow, b_grau, a_grau, a_rain, a_snow, diffelq, egccoef,  &
                         egicoef, erccoef, esccoef, esicoef, lcond, lsub, muelq, nzerog, nzeror, nzeros, rhog,      &
                         rhor, rhos, rv, therco, gam3, gamr1, gamr2, gamr3, gams1, gams2, gams3, gamg1, gamg2,      &
                         gamg3, ncol, nz)
    use sat_mod
    implicit none
    real(8), intent(in   ) :: tabs0  (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: pres   (ncol,nz) ! pressure,mb at scalar levels
    real(8), intent(in   ) :: rho    (ncol,nz) ! air density at pressure levels,kg/m3 
    real(8), intent(  out) :: accrsi (ncol,nz) ! Undocumented
    real(8), intent(  out) :: accrsc (ncol,nz) ! Undocumented
    real(8), intent(  out) :: coefice(ncol,nz) ! Undocumented
    real(8), intent(  out) :: evaps1 (ncol,nz) ! Undocumented
    real(8), intent(  out) :: evaps2 (ncol,nz) ! Undocumented
    real(8), intent(  out) :: accrgi (ncol,nz) ! Undocumented
    real(8), intent(  out) :: accrgc (ncol,nz) ! Undocumented
    real(8), intent(  out) :: evapg1 (ncol,nz) ! Undocumented
    real(8), intent(  out) :: evapg2 (ncol,nz) ! Undocumented
    real(8), intent(  out) :: accrrc (ncol,nz) ! Undocumented
    real(8), intent(  out) :: evapr1 (ncol,nz) ! Undocumented
    real(8), intent(  out) :: evapr2 (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: b_rain, b_snow, b_grau, a_grau, a_rain, a_snow, diffelq, egccoef, egicoef, erccoef, &
                              esccoef, esicoef, lcond, lsub, muelq, nzerog, nzeror, nzeros, rhog, rhor, rhos, rv, therco
    real(8), intent(  out) :: gam3, gamr1, gamr2, gamr3, gams1, gams2, gams3, gamg1, gamg2, gamg3
    integer, intent(in   ) :: ncol, nz



    real(8) :: pratio, coef1, coef2,estw,esti,rrr1,rrr2
    integer :: k,icol
    real(8), parameter :: pi = 3.14159265358979323846D0

    gam3 = 3.D0
    gamr1 = 3.D0+b_rain
    gamr2 = (5.D0+b_rain)/2.D0
    gamr3 = 4.D0+b_rain
    gams1 = 3.D0+b_snow
    gams2 = (5.D0+b_snow)/2.D0
    gams3 = 4.D0+b_snow
    gamg1 = 3.D0+b_grau
    gamg2 = (5.D0+b_grau)/2.D0
    gamg3 = 4.D0+b_grau
    gam3  = gamma(gam3 )
    gamr1 = gamma(gamr1)
    gamr2 = gamma(gamr2)
    gamr3 = gamma(gamr3)
    gams1 = gamma(gams1)
    gams2 = gamma(gams2)
    gams3 = gamma(gams3)
    gamg1 = gamma(gamg1)
    gamg2 = gamma(gamg2)
    gamg3 = gamma(gamg3)

    if(nint(gam3).ne.2) then
      print*,'cannot compute gamma-function in precip_init. Exiting...'
      stop
    endif

    !$acc parallel loop collapse(2)  async(asyncid)
    do k=1,nz
      do icol = 1 , ncol
        pratio = sqrt(1.29D0 / rho(icol,k))
        rrr1=393.D0/(tabs0(icol,k)+120.D0)*(tabs0(icol,k)/273.D0)**1.5D0
        rrr2=(tabs0(icol,k)/273.D0)**1.94D0*(1000.D0/pres(icol,k))
        estw = 100.D0*esatw_crm(tabs0(icol,k))
        esti = 100.D0*esati_crm(tabs0(icol,k))

        ! accretion by snow:
        coef1 = 0.25D0 * pi * nzeros * a_snow * gams1 * pratio/(pi * rhos * nzeros/rho(icol,k) ) ** ((3.D0+b_snow)/4.D0)
        coef2 = exp(0.025D0*(tabs0(icol,k) - 273.15D0))
        accrsi(icol,k) =  coef1 * coef2 * esicoef
        accrsc(icol,k) =  coef1 * esccoef
        coefice(icol,k) =  coef2

        ! evaporation of snow:
        coef1  =(lsub/(tabs0(icol,k)*rv)-1.)*lsub/(therco*rrr1*tabs0(icol,k))
        coef2  = rv*tabs0(icol,k)/(diffelq*rrr2*esti)
        evaps1(icol,k)  =  0.65D0*4.D0*nzeros/sqrt(pi*rhos*nzeros)/(coef1+coef2)/sqrt(rho(icol,k))
        evaps2(icol,k)  =  0.49D0*4.D0*nzeros*gams2*sqrt(a_snow/(muelq*rrr1))/(pi*rhos*nzeros)**((5.D0+b_snow)/8.D0) / (coef1+coef2) * rho(icol,k)**((1.D0+b_snow)/8.D0)*sqrt(pratio)

        ! accretion by graupel:
        coef1 = 0.25D0*pi*nzerog*a_grau*gamg1*pratio/(pi*rhog*nzerog/rho(icol,k))**((3.D0+b_grau)/4.D0)
        coef2 = exp(0.025D0*(tabs0(icol,k) - 273.15D0))
        accrgi(icol,k) =  coef1 * coef2 * egicoef
        accrgc(icol,k) =  coef1 * egccoef

        ! evaporation of graupel:
        coef1  =(lsub/(tabs0(icol,k)*rv)-1.D0)*lsub/(therco*rrr1*tabs0(icol,k))
        coef2  = rv*tabs0(icol,k)/(diffelq*rrr2*esti)
        evapg1(icol,k)  = 0.65D0*4.D0*nzerog/sqrt(pi*rhog*nzerog)/(coef1+coef2)/sqrt(rho(icol,k))
        evapg2(icol,k)  = 0.49D0*4.D0*nzerog*gamg2*sqrt(a_grau/(muelq*rrr1))/(pi * rhog * nzerog)**((5D0+b_grau)/8.D0) / (coef1+coef2) * rho(icol,k)**((1.D0+b_grau)/8.D0)*sqrt(pratio)

        ! accretion by rain:
        accrrc(icol,k)=  0.25D0 * pi * nzeror * a_rain * gamr1 * pratio/(pi * rhor * nzeror / rho(icol,k)) ** ((3.D0+b_rain)/4.D0)* erccoef

        ! evaporation of rain:
        coef1  =(lcond/(tabs0(icol,k)*rv)-1.D0)*lcond/(therco*rrr1*tabs0(icol,k))
        coef2  = rv*tabs0(icol,k)/(diffelq * rrr2 * estw)
        evapr1(icol,k)  =  0.78D0 * 2.D0 * pi * nzeror / &
        sqrt(pi * rhor * nzeror) / (coef1+coef2) / sqrt(rho(icol,k))
        evapr2(icol,k)  =  0.31D0 * 2.D0 * pi  * nzeror * gamr2 * 0.89D0 * sqrt(a_rain/(muelq*rrr1))/(pi * rhor * nzeror)**((5.D0+b_rain)/8.D0) / (coef1+coef2) * rho(icol,k)**((1.D0+b_rain)/8.D0)*sqrt(pratio)
      enddo
    enddo

  end subroutine precip_init


end module precip_init_mod
