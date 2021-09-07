
module precip_fall_mod
  implicit none

contains


  real(8) function term_vel_qp(qploc, rho, tabs, qp_threshold, tprmin, &
                               a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, csnow)
    implicit none
    real(8), intent(in   ) :: qploc, rho, tabs, qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, &
                              cgrau, vsnow, csnow
    real(8) :: wmax, omp, omg, qrr, qss, qgg

    term_vel_qp = 0.
    if(qploc.gt.qp_threshold) then
      omp = max(real(0.,8),min(real(1.,8),(tabs-tprmin)*a_pr))
      if(omp.eq.1.) then
        term_vel_qp = vrain*(rho*qploc)**crain
      elseif(omp.eq.0.) then
        omg = max(real(0.,8),min(real(1.,8),(tabs-tgrmin)*a_gr))
        qgg=omg*qploc
        qss=qploc-qgg
        term_vel_qp = (omg*vgrau*(rho*qgg)**cgrau &
        +(1.-omg)*vsnow*(rho*qss)**csnow)
      else
        omg = max(real(0.,8),min(real(1.,8),(tabs-tgrmin)*a_gr))
        qrr=omp*qploc
        qss=qploc-qrr
        qgg=omg*qss
        qss=qss-qgg
        term_vel_qp = (omp*vrain*(rho*qrr)**crain + (1.-omp)*(omg*vgrau*(rho*qgg)**cgrau + &
                       (1.-omg)*vsnow*(rho*qss)**csnow))
      endif
    endif
  end function term_vel_qp



  subroutine micro_precip_fall(rho, adz, dz, rhow, qp, t, tabs, qpfall, tlat, precflux, precsfc, precssfc,       &
                               prec_xy, qp_threshold, tprmin, a_pr, tgrmin, a_gr, vsnow, dtn, fac_cond, fac_fus, &
                               b_rain, b_snow, b_grau, a_rain, a_snow, a_grau, gamr3, gams3, gamg3, rhor, rhos,  &
                               rhog, nzeror, nzeros, nzerog, ncrms, nx, ny, nzm)
    implicit none
    real(8), intent(in   ) :: rho     (ncrms      ,nzm)
    real(8), intent(in   ) :: adz     (ncrms      ,nzm)
    real(8), intent(in   ) :: dz      (ncrms          )
    real(8), intent(in   ) :: rhow    (ncrms      ,nzm)
    real(8), intent(inout) :: qp      (ncrms,nx,ny,nzm)
    real(8), intent(inout) :: t       (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: tabs    (ncrms,nx,ny,nzm)
    real(8), intent(inout) :: qpfall  (ncrms      ,nzm)
    real(8), intent(inout) :: tlat    (ncrms      ,nzm)
    real(8), intent(inout) :: precflux(ncrms      ,nzm)
    real(8), intent(inout) :: precsfc (ncrms,nx,ny    )
    real(8), intent(inout) :: precssfc(ncrms,nx,ny    )
    real(8), intent(inout) :: prec_xy (ncrms,nx,ny    )
    real(8), intent(in   ) :: qp_threshold, tprmin, a_pr, tgrmin, a_gr, dtn, fac_cond, fac_fus, b_rain, b_snow, &
                              b_grau, a_rain, a_snow, a_grau, gamr3, gams3, gamg3, rhor, rhos, rhog, nzeror,    &
                              nzeros, nzerog
    integer, intent(in   ) :: ncrms, nx, ny, nzm

    real(8), allocatable :: omega(:,:,:,:)
    real(8)              :: crain, csnow, cgrau, vrain, vsnow, vgrau
    integer              :: i,j,k,icrm
    real(8), parameter :: pi = 3.14159265358979323846D0

    allocate(omega(ncrms,nx,ny,nzm))

    crain = b_rain / 4.D0
    csnow = b_snow / 4.D0
    cgrau = b_grau / 4.D0
    vrain = a_rain * gamr3 / 6.D0 / (pi * rhor * nzeror) ** crain
    vsnow = a_snow * gams3 / 6.D0 / (pi * rhos * nzeros) ** csnow
    vgrau = a_grau * gamg3 / 6.D0 / (pi * rhog * nzerog) ** cgrau

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            omega(icrm,i,j,k) = max(real(0.,8),min(real(1.,8),(tabs(icrm,i,j,k)-tprmin)*a_pr))
          enddo
        enddo
      enddo
    enddo

    call precip_fall(rho, adz, dz, omega, rhow, qp, t, tabs, qpfall, tlat, precflux, precsfc, precssfc,    &
                     prec_xy, qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, &
                     csnow, dtn, fac_cond, fac_fus, ncrms, nx, ny, nzm)

    deallocate(omega)

  end subroutine micro_precip_fall



  ! Positively definite monotonic advection with non-oscillatory option
  ! and gravitational sedimentation
  subroutine precip_fall(rho, adz, dz, omega, rhow, qp, t, tabs, qpfall, tlat, precflux, precsfc, precssfc,    &
                         prec_xy, qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, &
                         csnow, dtn, fac_cond, fac_fus, ncrms, nx, ny, nzm)
    implicit none
    real(8), intent(in   ) :: rho     (ncrms      ,nzm)
    real(8), intent(in   ) :: adz     (ncrms      ,nzm)
    real(8), intent(in   ) :: dz      (ncrms          )
    real(8), intent(in   ) :: omega   (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: rhow    (ncrms      ,nzm)
    real(8), intent(inout) :: qp      (ncrms,nx,ny,nzm)
    real(8), intent(inout) :: t       (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: tabs    (ncrms,nx,ny,nzm)
    real(8), intent(inout) :: qpfall  (ncrms      ,nzm)
    real(8), intent(inout) :: tlat    (ncrms      ,nzm)
    real(8), intent(inout) :: precflux(ncrms      ,nzm)
    real(8), intent(inout) :: precsfc (ncrms,nx,ny    )
    real(8), intent(inout) :: precssfc(ncrms,nx,ny    )
    real(8), intent(inout) :: prec_xy (ncrms,nx,ny    )
    real(8), intent(in   ) :: qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, csnow, dtn, &
                              fac_cond, fac_fus
    integer, intent(in   ) :: ncrms, nx, ny, nzm

    ! Local:
    real(8), allocatable :: mx     (:,:,:,:)
    real(8), allocatable :: mn     (:,:,:,:)
    real(8), allocatable :: lfac   (:,:,:,:)
    real(8), allocatable :: www    (:,:,:,:)
    real(8), allocatable :: fz     (:,:,:,:)
    real(8), allocatable :: wp     (:,:,:,:)
    real(8), allocatable :: tmp_qp (:,:,:,:)
    real(8), allocatable :: irhoadz(:,:)
    real(8), allocatable :: iwmax  (:,:)
    real(8), allocatable :: rhofac (:,:)
    real(8) :: prec_cfl
    real(8) :: eps
    integer :: i,j,k,kc,kb,icrm
    real(8) :: y,pp,pn
    real(8) :: lat_heat, wmax
    integer nprec, iprec
    real(8) :: flagstat, tmp, tvel

    !Statement functions
    pp(y)= max(real(0.,8),y)
    pn(y)=-min(real(0.,8),y)

    eps = 1.D-10

    allocate( mx     (ncrms,nx,ny,nzm  ) )
    allocate( mn     (ncrms,nx,ny,nzm  ) )
    allocate( lfac   (ncrms,nx,ny,nzm+1) )
    allocate( www    (ncrms,nx,ny,nzm+1) )
    allocate( fz     (ncrms,nx,ny,nzm+1) )
    allocate( wp     (ncrms,nx,ny,nzm  ) )
    allocate( tmp_qp (ncrms,nx,ny,nzm  ) )
    allocate( irhoadz(ncrms,nzm) )
    allocate( iwmax  (ncrms,nzm) )
    allocate( rhofac (ncrms,nzm) )

    !$acc parallel loop gang vector collapse(2) async(asyncid)
    do k = 1,nzm
      do icrm = 1 , ncrms
        rhofac(icrm,k) = sqrt(1.29D0/rho(icrm,k))
        irhoadz(icrm,k) = 1.D0/(rho(icrm,k)*adz(icrm,k)) ! Useful factor
        kb = max(1,k-1)
        wmax       = dz(icrm)*adz(icrm,kb)/dtn   ! Velocity equivalent to a cfl of 1.0.
        iwmax(icrm,k)   = 1.D0/wmax
      enddo
    enddo

    !   Add sedimentation of precipitation field to the vert. vel.
    prec_cfl = 0.D0
    !$acc parallel loop gang vector collapse(4) reduction(max:prec_cfl) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            lfac(icrm,i,j,k) = fac_cond + (1-omega(icrm,i,j,k))*fac_fus
            flagstat = 1.D0
            tvel = term_vel_qp(qp(icrm,i,j,k),rho(icrm,k),tabs(icrm,i,j,k),qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,&
                               a_gr,vgrau,cgrau,vsnow,csnow)
            wp(icrm,i,j,k)=rhofac(icrm,k)*tvel
            tmp = wp(icrm,i,j,k)*iwmax(icrm,k)
            prec_cfl = max(prec_cfl,tmp) ! Keep column maximum CFL
            wp(icrm,i,j,k) = -wp(icrm,i,j,k)*rhow(icrm,k)*dtn/dz(icrm)
            if (k == 1) then
              fz(icrm,i,j,nzm+1)=0.D0
              www(icrm,i,j,nzm+1)=0.D0
              lfac(icrm,i,j,nzm+1)=0D0
            endif
          enddo  ! k
        enddo
      enddo
    enddo

    ! If maximum CFL due to precipitation velocity is greater than 0.9,
    ! take more than one advection step to maintain stability.
    if (prec_cfl.gt.0.9D0) then
      nprec = CEILING(prec_cfl/0.9D0)
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k = 1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              ! wp already includes factor of dt, so reduce it by a
              ! factor equal to the number of precipitation steps.
              wp(icrm,i,j,k) = wp(icrm,i,j,k)/real(nprec,8)
            enddo
          enddo
        enddo
      enddo
    else
      nprec = 1
    endif

    !  loop over iterations
    do iprec = 1,nprec
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k = 1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              tmp_qp(icrm,i,j,k) = qp(icrm,i,j,k) ! Temporary array for qp in this column
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
                kc=min(nzm,k+1)
                kb=max(1,k-1)
                mx(icrm,i,j,k)=max(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k))
                mn(icrm,i,j,k)=min(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k))
              ! Define upwind precipitation flux
              fz(icrm,i,j,k)=tmp_qp(icrm,i,j,k)*wp(icrm,i,j,k)
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc=k+1
              tmp_qp(icrm,i,j,k)=tmp_qp(icrm,i,j,k)-(fz(icrm,i,j,kc)-fz(icrm,i,j,k))*irhoadz(icrm,k) !Update temporary qp
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              ! Also, compute anti-diffusive correction to previous
              ! (upwind) approximation to the flux
              kb=max(1,k-1)
              ! The precipitation velocity is a cell-centered quantity,
              ! since it is computed from the cell-centered
              ! precipitation mass fraction.  Therefore, a reformulated
              ! anti-diffusive flux is used here which accounts for
              ! this and results in reduced numerical diffusion.
              www(icrm,i,j,k) = 0.5D0*(1.+wp(icrm,i,j,k)*irhoadz(icrm,k))*(tmp_qp(icrm,i,j,kb)*wp(icrm,i,j,kb) - &
                                     tmp_qp(icrm,i,j,k)*wp(icrm,i,j,k)) ! works for wp(k)<0
            enddo
          enddo
        enddo
      enddo

      !---------- non-osscilatory option ---------------
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc=min(nzm,k+1)
              kb=max(1,k-1)
              mx(icrm,i,j,k)=max(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k),mx(icrm,i,j,k))
              mn(icrm,i,j,k)=min(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k),mn(icrm,i,j,k))
              kc=min(nzm,k+1)
              mx(icrm,i,j,k)=rho(icrm,k)*adz(icrm,k)*(mx(icrm,i,j,k)-tmp_qp(icrm,i,j,k))/(pn(www(icrm,i,j,kc)) + pp(www(icrm,i,j,k))+eps)
              mn(icrm,i,j,k)=rho(icrm,k)*adz(icrm,k)*(tmp_qp(icrm,i,j,k)-mn(icrm,i,j,k))/(pp(www(icrm,i,j,kc)) + pn(www(icrm,i,j,k))+eps)
            enddo
          enddo
        enddo
      enddo
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kb=max(1,k-1)
              ! Add limited flux correction to fz(k).
              fz(icrm,i,j,k) = fz(icrm,i,j,k) + pp(www(icrm,i,j,k))*min(real(1.,8),mx(icrm,i,j,k), mn(icrm,i,j,kb)) - &
                                                pn(www(icrm,i,j,k))*min(real(1.,8),mx(icrm,i,j,kb),mn(icrm,i,j,k)) ! Anti-diffusive flux
            enddo
          enddo
        enddo
      enddo

      ! Update precipitation mass fraction and liquid-ice static
      ! energy using precipitation fluxes computed in this column.
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do j=1,ny
        do i=1,nx
          do k=1,nzm
            do icrm = 1 , ncrms
              kc=k+1
              ! Update precipitation mass fraction.
              ! Note that fz is the total flux, including both the
              ! upwind flux and the anti-diffusive correction.
              qp(icrm,i,j,k)=qp(icrm,i,j,k)-(fz(icrm,i,j,kc)-fz(icrm,i,j,k))*irhoadz(icrm,k)
              tmp = -(fz(icrm,i,j,kc)-fz(icrm,i,j,k))*irhoadz(icrm,k)*flagstat  ! For qp budget
              !$acc atomic update
              qpfall(icrm,k)=qpfall(icrm,k)+tmp
              lat_heat = -(lfac(icrm,i,j,kc)*fz(icrm,i,j,kc)-lfac(icrm,i,j,k)*fz(icrm,i,j,k))*irhoadz(icrm,k)
              t(icrm,i,j,k)=t(icrm,i,j,k)-lat_heat
              !$acc atomic update
              tlat(icrm,k)=tlat(icrm,k)-lat_heat            ! For energy budget
              tmp = fz(icrm,i,j,k)*flagstat
              !$acc atomic update
              precflux(icrm,k) = precflux(icrm,k) - tmp   ! For statistics
              if (k == 1) then
                precsfc(icrm,i,j) = precsfc(icrm,i,j) - fz(icrm,i,j,1)*flagstat ! For statistics
                precssfc(icrm,i,j) = precssfc(icrm,i,j) - fz(icrm,i,j,1)*(1.-omega(icrm,i,j,1))*flagstat ! For statistics
                prec_xy(icrm,i,j) = prec_xy(icrm,i,j) - fz(icrm,i,j,1)*flagstat ! For 2D output
              endif
            enddo
          enddo
        enddo
      enddo

      if (iprec.lt.nprec) then
        ! Re-compute precipitation velocity using new value of qp.
        !$acc parallel loop gang vector collapse(4) async(asyncid)
        do j=1,ny
          do i=1,nx
            do k=1,nzm
              do icrm = 1 , ncrms
                !Passing variables via first index because of PGI bug with pointers
                tvel = term_vel_qp(qp(icrm,i,j,k),rho(icrm,k),tabs(icrm,i,j,k),qp_threshold,tprmin,a_pr,vrain,crain,&
                                   tgrmin,a_gr,vgrau,cgrau,vsnow,csnow)
                wp(icrm,i,j,k) = rhofac(icrm,k)*tvel
                ! Decrease precipitation velocity by factor of nprec
                wp(icrm,i,j,k) = -wp(icrm,i,j,k)*rhow(icrm,k)*dtn/dz(icrm)/real(nprec,8)
                ! Note: Don't bother checking CFL condition at each
                ! substep since it's unlikely that the CFL will
                ! increase very much between substeps when using
                ! monotonic advection schemes.
                if (k == 1) then
                  fz(icrm,i,j,nzm+1)=0.D0
                  www(icrm,i,j,nzm+1)=0.D0
                  lfac(icrm,i,j,nzm+1)=0.D0
                endif
              enddo
            enddo
          enddo
        enddo
      endif

    enddo
    
    deallocate( mx      )
    deallocate( mn      )
    deallocate( lfac    )
    deallocate( www     )
    deallocate( fz      )
    deallocate( wp      )
    deallocate( tmp_qp  )
    deallocate( irhoadz )
    deallocate( iwmax   )
    deallocate( rhofac  )

  end subroutine precip_fall

endmodule precip_fall_mod
