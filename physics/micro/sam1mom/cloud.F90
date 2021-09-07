module cloud_mod
  implicit none

contains

  ! Condensation of cloud water/cloud ice.
  subroutine cloud(q, tabs, t, gamaz, qp, pres, qn, tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, &
                   tgrmax, tgrmin, ncrms, nx, ny, nzm)
    use sat_mod
    implicit none
    real(8), intent(inout) :: q    (ncrms,nx,ny,nzm)
    real(8), intent(  out) :: tabs (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: t    (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: gamaz(ncrms      ,nzm)
    real(8), intent(inout) :: qp   (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: pres (ncrms      ,nzm)
    real(8), intent(  out) :: qn   (ncrms,nx,ny,nzm)  ! cloud condensate (liquid + ice)
    real(8), intent(in   ) :: tbgmax, tbgmin, tprmax, tprmin, fac_cond, fac_fus, fac_sub, tgrmax, tgrmin
    integer, intent(in   ) :: ncrms, nx, ny, nzm

    integer :: i, j, k, kb, kc,icrm
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

    !$acc parallel loop collapse(4)
    do k = 1, nzm
      do j = 1, ny
        do i = 1, nx
          do icrm = 1 , ncrms
            q(icrm,i,j,k)=max(real(0.,8),q(icrm,i,j,k))
            ! Initial guess for temperature assuming no cloud water/ice:
            tabs(icrm,i,j,k) = t(icrm,i,j,k)-gamaz(icrm,k)
            tabs1=(tabs(icrm,i,j,k)+fac1*qp(icrm,i,j,k))/(1.+fac2*qp(icrm,i,j,k))
            ! Warm cloud:
            if(tabs1.ge.tbgmax) then
              tabs1=tabs(icrm,i,j,k)+fac_cond*qp(icrm,i,j,k)
              qsatt = qsatw_crm(tabs1,pres(icrm,k))
              ! Ice cloud:
            elseif(tabs1.le.tbgmin) then
              tabs1=tabs(icrm,i,j,k)+fac_sub*qp(icrm,i,j,k)
              qsatt = qsati_crm(tabs1,pres(icrm,k))
              ! Mixed-phase cloud:
            else
              om = an*tabs1-bn
              qsatt = om*qsatw_crm(tabs1,pres(icrm,k))+(1.-om)*qsati_crm(tabs1,pres(icrm,k))
            endif
            !  Test if condensation is possible:
            if(q(icrm,i,j,k).gt.qsatt) then
              niter=0
              dtabs = 100.
              do while(abs(dtabs).gt.0.01D0.and.niter.lt.10)
                if(tabs1.ge.tbgmax) then
                  om=1.
                  lstarn=fac_cond
                  dlstarn=0.
                  qsatt=qsatw_crm(tabs1,pres(icrm,k))
                  dqsat=dtqsatw_crm(tabs1,pres(icrm,k))
                else if(tabs1.le.tbgmin) then
                  om=0.
                  lstarn=fac_sub
                  dlstarn=0.
                  qsatt=qsati_crm(tabs1,pres(icrm,k))
                  dqsat=dtqsati_crm(tabs1,pres(icrm,k))
                else
                  om=an*tabs1-bn
                  lstarn=fac_cond+(1.-om)*fac_fus
                  dlstarn=an*fac_fus
                  qsatt=om*qsatw_crm(tabs1,pres(icrm,k))+(1.-om)*qsati_crm(tabs1,pres(icrm,k))
                  dqsat=om*dtqsatw_crm(tabs1,pres(icrm,k))+(1.-om)*dtqsati_crm(tabs1,pres(icrm,k))
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
                fff = tabs(icrm,i,j,k)-tabs1+lstarn*(q(icrm,i,j,k)-qsatt)+lstarp*qp(icrm,i,j,k)
                dfff=dlstarn*(q(icrm,i,j,k)-qsatt)+dlstarp*qp(icrm,i,j,k)-lstarn*dqsat-1.
                dtabs=-fff/dfff
                niter=niter+1
                tabs1=tabs1+dtabs
              enddo
              qsatt = qsatt + dqsat * dtabs
              qn(icrm,i,j,k) = max(real(0.,8),q(icrm,i,j,k)-qsatt)
            else
              qn(icrm,i,j,k) = 0.
            endif
            tabs(icrm,i,j,k) = tabs1
            qp(icrm,i,j,k) = max(real(0.,8),qp(icrm,i,j,k)) ! just in case
          enddo
        enddo
      enddo
    enddo

  end subroutine cloud

end module cloud_mod
