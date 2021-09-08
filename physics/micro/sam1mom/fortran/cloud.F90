module cloud_mod
  implicit none

contains

  ! Condensation of cloud water/cloud ice.
  subroutine cloud(q, tabs, t, gamaz, qp, pres, qn, tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, &
                   tgrmax, tgrmin, ncol, nz)
    use sat_mod
    implicit none
    real(8), intent(inout) :: q    (ncol,nz)  ! total nonprecipitating water
    real(8), intent(  out) :: tabs (ncol,nz)  ! temperature
    real(8), intent(in   ) :: t    (ncol,nz)  ! liquid/ice water static energy 
    real(8), intent(in   ) :: gamaz(ncol,nz)  ! grav/cp*z
    real(8), intent(inout) :: qp   (ncol,nz)  ! total precipitating water
    real(8), intent(in   ) :: pres (ncol,nz)  ! pressure,mb at scalar levels
    real(8), intent(  out) :: qn   (ncol,nz)  ! cloud condensate (liquid + ice)
    real(8), intent(in   ) :: tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, tgrmax, tgrmin
    integer, intent(in   ) :: ncol, nz

    integer :: icol, k
    real(8) :: dtabs, tabs1, an, bn, ap, bp, om, ag, omp
    real(8) :: fac1, fac2
    real(8) :: fff, dfff, qsatt, dqsat
    real(8) :: lstarn, dlstarn, lstarp, dlstarp
    integer :: niter

    an = 1./(tbgmax-tbgmin)
    bn = tbgmin * an
    ap = 1./(tprmax-tprmin)
    bp = tprmin * ap
    fac1 = fac_cond+(1+bp)*fac_fus
    fac2 = fac_fus*ap
    ag = 1./(tgrmax-tgrmin)

    !$acc parallel loop collapse(2)
    do k = 1, nz
      do icol = 1, ncol
        q(icol,k)=max(real(0.,8),q(icol,k))
        ! Initial guess for temperature assuming no cloud water/ice:
        tabs(icol,k) = t(icol,k)-gamaz(icol,k)
        tabs1=(tabs(icol,k)+fac1*qp(icol,k))/(1.+fac2*qp(icol,k))
        ! Warm cloud:
        if(tabs1.ge.tbgmax) then
          tabs1=tabs(icol,k)+fac_cond*qp(icol,k)
          qsatt = qsatw_crm(tabs1,pres(icol,k))
          ! Ice cloud:
        elseif(tabs1.le.tbgmin) then
          tabs1=tabs(icol,k)+fac_sub*qp(icol,k)
          qsatt = qsati_crm(tabs1,pres(icol,k))
          ! Mixed-phase cloud:
        else
          om = an*tabs1-bn
          qsatt = om*qsatw_crm(tabs1,pres(icol,k))+(1.-om)*qsati_crm(tabs1,pres(icol,k))
        endif
        !  Test if condensation is possible:
        if(q(icol,k).gt.qsatt) then
          niter=0
          dtabs = 100.
          do while(abs(dtabs).gt.0.01D0.and.niter.lt.10)
            if(tabs1.ge.tbgmax) then
              om=1.
              lstarn=fac_cond
              dlstarn=0.
              qsatt=qsatw_crm(tabs1,pres(icol,k))
              dqsat=dtqsatw_crm(tabs1,pres(icol,k))
            else if(tabs1.le.tbgmin) then
              om=0.
              lstarn=fac_sub
              dlstarn=0.
              qsatt=qsati_crm(tabs1,pres(icol,k))
              dqsat=dtqsati_crm(tabs1,pres(icol,k))
            else
              om=an*tabs1-bn
              lstarn=fac_cond+(1.-om)*fac_fus
              dlstarn=an*fac_fus
              qsatt=om*qsatw_crm(tabs1,pres(icol,k))+(1.-om)*qsati_crm(tabs1,pres(icol,k))
              dqsat=om*dtqsatw_crm(tabs1,pres(icol,k))+(1.-om)*dtqsati_crm(tabs1,pres(icol,k))
            endif
            if(tabs1.ge.tprmax) then
              omp=1.
              lstarp=fac_cond
              dlstarp=0.
            else if(tabs1.le.tprmin) then
              omp=0.
              lstarp=fac_sub
              dlstarp=0.
            else
              omp=ap*tabs1-bp
              lstarp=fac_cond+(1.-omp)*fac_fus
              dlstarp=ap*fac_fus
            endif
            fff = tabs(icol,k)-tabs1+lstarn*(q(icol,k)-qsatt)+lstarp*qp(icol,k)
            dfff=dlstarn*(q(icol,k)-qsatt)+dlstarp*qp(icol,k)-lstarn*dqsat-1.
            dtabs=-fff/dfff
            niter=niter+1
            tabs1=tabs1+dtabs
          enddo
          qsatt = qsatt + dqsat * dtabs
          qn(icol,k) = max(real(0.,8),q(icol,k)-qsatt)
        else
          qn(icol,k) = 0.
        endif
        tabs(icol,k) = tabs1
        qp(icol,k) = max(real(0.,8),qp(icol,k)) ! just in case
      enddo
    enddo

  end subroutine cloud

end module cloud_mod
