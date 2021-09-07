module precip_proc_mod
  implicit none

contains

  subroutine precip_proc(qpsrc, qpevp, qn, qp, tabs, coefice, accrrc, accrsc, accrsi, accrgc, accrgi, q, pres, &
                         evapr1, evapr2, evaps1, evaps2, evapg1, evapg2, a_bg, a_gr, a_pr, alphaelq, b_grau,   &
                         b_rain, b_snow, betaelq, dtn, qci0, qcw0, qp_threshold, tbgmin, tgrmin, tprmin,       &
                         ncol, nz)
    use sat_mod
    implicit none
    real(8), intent(  out) :: qpsrc  (ncol,nz) ! source of precipitation microphysical processes
    real(8), intent(  out) :: qpevp  (ncol,nz) ! sink of precipitating water due to evaporation
    real(8), intent(inout) :: qn     (ncol,nz) ! cloud condensate (liquid + ice)
    real(8), intent(inout) :: qp     (ncol,nz) ! total precipitating water
    real(8), intent(in   ) :: tabs   (ncol,nz) ! temperature
    real(8), intent(in   ) :: coefice(ncol,nz) ! Undocumented
    real(8), intent(in   ) :: accrrc (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: accrsc (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: accrsi (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: accrgc (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: accrgi (ncol,nz) ! Undocumented
    real(8), intent(inout) :: q      (ncol,nz) ! total nonprecipitating water
    real(8), intent(in   ) :: pres   (ncol,nz) ! pressure,mb at scalar levels
    real(8), intent(in   ) :: evapr1 (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: evapr2 (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: evaps1 (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: evaps2 (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: evapg1 (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: evapg2 (ncol,nz) ! Undocumented
    real(8), intent(in   ) :: a_bg, a_gr, a_pr, alphaelq, b_grau, b_rain, b_snow, betaelq, dtn, qci0, qcw0, &
                              qp_threshold, tbgmin, tgrmin, tprmin
    integer, intent(in   ) :: ncol, nz

    integer :: icol, k
    real(8) :: autor, autos, accrr, accris, accrcs, accrig, accrcg
    real(8) :: dq, omn, omp, omg, qsatt
    real(8) :: pows1, pows2, powg1, powg2, powr1, powr2, tmp
    real(8) :: qii, qcc, qrr, qss, qgg

    powr1 = (3 + b_rain) / 4.D0
    powr2 = (5 + b_rain) / 8.D0
    pows1 = (3 + b_snow) / 4.D0
    pows2 = (5 + b_snow) / 8.D0
    powg1 = (3 + b_grau) / 4.D0
    powg2 = (5 + b_grau) / 8.D0

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nz
      do icol = 1 , ncol
        qpsrc(icol,k)=0.
        qpevp(icol,k)=0.
      enddo
    enddo

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nz
      do icol = 1 , ncol
        !-------     Autoconversion/accretion
        if(qn(icol,k)+qp(icol,k).gt.0.) then
          omn = max(real(0.,8),min(real(1.,8),(tabs(icol,k)-tbgmin)*a_bg))
          omp = max(real(0.,8),min(real(1.,8),(tabs(icol,k)-tprmin)*a_pr))
          omg = max(real(0.,8),min(real(1.,8),(tabs(icol,k)-tgrmin)*a_gr))

          if(qn(icol,k).gt.0.) then

            qcc = qn(icol,k) * omn
            qii = qn(icol,k) * (1.-omn)

            if(qcc .gt. qcw0) then
              autor = alphaelq
            else
              autor = 0.
            endif

            if(qii .gt. qci0) then
              autos = betaelq*coefice(icol,k)
            else
              autos = 0.
            endif

            accrr = 0.
            if(omp.gt.0.001D0) then
              qrr = qp(icol,k) * omp
              accrr = accrrc(icol,k) * qrr ** powr1
            endif
            accrcs = 0.
            accris = 0.
            if(omp.lt.0.999D0.and.omg.lt.0.999D0) then
              qss = qp(icol,k) * (1.-omp)*(1.-omg)
              tmp = qss ** pows1
              accrcs = accrsc(icol,k) * tmp
              accris = accrsi(icol,k) * tmp
            endif
            accrcg = 0.
            accrig = 0.
            if(omp.lt.0.999D0.and.omg.gt.0.001D0) then
              qgg = qp(icol,k) * (1.-omp)*omg
              tmp = qgg ** powg1
              accrcg = accrgc(icol,k) * tmp
              accrig = accrgi(icol,k) * tmp
            endif
            qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
            qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
            dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+(accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))
            dq = min(dq,qn(icol,k))
            qp(icol,k) = qp(icol,k) + dq
            q(icol,k) = q(icol,k) - dq
            qn(icol,k) = qn(icol,k) - dq
            qpsrc(icol,k) = qpsrc(icol,k) + dq

          elseif(qp(icol,k).gt.qp_threshold.and.qn(icol,k).eq.0.) then

            qsatt = 0.
            if(omn.gt.0.001D0) qsatt = qsatt + omn*qsatw_crm(tabs(icol,k),pres(icol,k))
            if(omn.lt.0.999D0) qsatt = qsatt + (1.-omn)*qsati_crm(tabs(icol,k),pres(icol,k))
            dq = 0.
            if(omp.gt.0.001D0) then
              qrr = qp(icol,k) * omp
              dq = dq + evapr1(icol,k)*sqrt(qrr) + evapr2(icol,k)*qrr**powr2
            endif
            if(omp.lt.0.999D0.and.omg.lt.0.999D0) then
              qss = qp(icol,k) * (1.-omp)*(1.-omg)
              dq = dq + evaps1(icol,k)*sqrt(qss) + evaps2(icol,k)*qss**pows2
            endif
            if(omp.lt.0.999D0.and.omg.gt.0.001D0) then
              qgg = qp(icol,k) * (1.-omp)*omg
              dq = dq + evapg1(icol,k)*sqrt(qgg) + evapg2(icol,k)*qgg**powg2
            endif
            dq = dq * dtn * (q(icol,k) /qsatt-1.)
            dq = max(-0.5D0*qp(icol,k),dq)
            qp(icol,k) = qp(icol,k) + dq
            q(icol,k) = q(icol,k) - dq
            qpevp(icol,k) = qpevp(icol,k) + dq

          else

            q(icol,k) = q(icol,k) + qp(icol,k)
            qpevp(icol,k) = qpevp(icol,k) - qp(icol,k)
            qp(icol,k) = 0.

          endif

        endif

        dq = qp(icol,k)
        qp(icol,k)=max(real(0.,8),qp(icol,k))
        q(icol,k) = q(icol,k) + (dq-qp(icol,k))
      enddo
    enddo

  end subroutine precip_proc

end module precip_proc_mod
