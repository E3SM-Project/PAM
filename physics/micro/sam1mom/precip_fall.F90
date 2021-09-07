
module precip_fall_mod
  implicit none

contains


  real(8) function term_vel_qp(qploc, rho, tabs, qp_threshold, tprmin, &
                               a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, csnow)
    implicit none
    real(8), intent(in   ) :: qploc, rho, tabs, qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, &
                              cgrau, vsnow, csnow
    real(8) :: omp, omg, qrr, qss, qgg

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
                               prec_xy, qp_threshold, tprmin, a_pr, tgrmin, a_gr, dtn, fac_cond, fac_fus, &
                               b_rain, b_snow, b_grau, a_rain, a_snow, a_grau, gamr3, gams3, gamg3, rhor, rhos,  &
                               rhog, nzeror, nzeros, nzerog, ncol, nz)
    implicit none
    real(8), intent(in   ) :: rho     (ncol,nz  ) ! air density at pressure levels,kg/m3 
    real(8), intent(in   ) :: adz     (ncol,nz  ) ! ratio of the thickness of scalar levels to dz 
    real(8), intent(in   ) :: dz      (ncol     ) ! constant grid spacing in z direction (when dz_constant=.true.)
    real(8), intent(in   ) :: rhow    (ncol,nz+1) ! air density at vertical velocity levels,kg/m3
    real(8), intent(inout) :: qp      (ncol,nz  ) ! total precipitating water
    real(8), intent(inout) :: t       (ncol,nz  ) ! liquid/ice water static energy 
    real(8), intent(in   ) :: tabs    (ncol,nz  ) ! temperature
    real(8), intent(inout) :: qpfall  (ncol,nz  ) ! for statistics
    real(8), intent(inout) :: tlat    (ncol,nz  ) ! for statistics
    real(8), intent(inout) :: precflux(ncol,nz  ) ! for statistics
    real(8), intent(inout) :: precsfc (ncol     ) ! surface precip. rate
    real(8), intent(inout) :: precssfc(ncol     ) ! surface ice precip. rate
    real(8), intent(inout) :: prec_xy (ncol     ) ! mean precip. rate for outout
    real(8), intent(in   ) :: qp_threshold, tprmin, a_pr, tgrmin, a_gr, dtn, fac_cond, fac_fus, b_rain, b_snow, &
                              b_grau, a_rain, a_snow, a_grau, gamr3, gams3, gamg3, rhor, rhos, rhog, nzeror,    &
                              nzeros, nzerog
    integer, intent(in   ) :: ncol, nz

    real(8), allocatable :: omega(:,:)
    real(8)              :: crain, csnow, cgrau, vrain, vsnow, vgrau
    integer              :: icol, k
    real(8), parameter :: pi = 3.14159265358979323846D0

    allocate(omega(ncol,nz))

    crain = b_rain / 4.D0
    csnow = b_snow / 4.D0
    cgrau = b_grau / 4.D0
    vrain = a_rain * gamr3 / 6.D0 / (pi * rhor * nzeror) ** crain
    vsnow = a_snow * gams3 / 6.D0 / (pi * rhos * nzeros) ** csnow
    vgrau = a_grau * gamg3 / 6.D0 / (pi * rhog * nzerog) ** cgrau

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nz
      do icol = 1 , ncol
        omega(icol,k) = max(real(0.,8),min(real(1.,8),(tabs(icol,k)-tprmin)*a_pr))
      enddo
    enddo

    call precip_fall(rho, adz, dz, omega, rhow, qp, t, tabs, qpfall, tlat, precflux, precsfc, precssfc,    &
                     prec_xy, qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, &
                     csnow, dtn, fac_cond, fac_fus, ncol, nz)

    deallocate(omega)

  end subroutine micro_precip_fall



  ! Positively definite monotonic advection with non-oscillatory option
  ! and gravitational sedimentation
  subroutine precip_fall(rho, adz, dz, omega, rhow, qp, t, tabs, qpfall, tlat, precflux, precsfc, precssfc,    &
                         prec_xy, qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, &
                         csnow, dtn, fac_cond, fac_fus, ncol, nz)
    implicit none
    real(8), intent(in   ) :: rho     (ncol,nz  ) ! air density at pressure levels,kg/m3 
    real(8), intent(in   ) :: adz     (ncol,nz  ) ! ratio of the thickness of scalar levels to dz 
    real(8), intent(in   ) :: dz      (ncol      ) ! constant grid spacing in z direction (when dz_constant=.true.)
    real(8), intent(in   ) :: omega   (ncol,nz  ) ! Undocumented
    real(8), intent(in   ) :: rhow    (ncol,nz+1) ! air density at vertical velocity levels,kg/m3
    real(8), intent(inout) :: qp      (ncol,nz  ) ! total precipitating water
    real(8), intent(inout) :: t       (ncol,nz  ) ! liquid/ice water static energy 
    real(8), intent(in   ) :: tabs    (ncol,nz  ) ! temperature
    real(8), intent(inout) :: qpfall  (ncol,nz  ) ! for statistics
    real(8), intent(inout) :: tlat    (ncol,nz  ) ! for statistics
    real(8), intent(inout) :: precflux(ncol,nz  ) ! for statistics
    real(8), intent(inout) :: precsfc (ncol      ) ! surface precip. rate
    real(8), intent(inout) :: precssfc(ncol      ) ! surface ice precip. rate
    real(8), intent(inout) :: prec_xy (ncol      ) ! mean precip. rate for outout
    real(8), intent(in   ) :: qp_threshold, tprmin, a_pr, vrain, crain, tgrmin, a_gr, vgrau, cgrau, vsnow, csnow, dtn, &
                              fac_cond, fac_fus
    integer, intent(in   ) :: ncol, nz

    ! Local:
    real(8), allocatable :: mx     (:,:)
    real(8), allocatable :: mn     (:,:)
    real(8), allocatable :: lfac   (:,:)
    real(8), allocatable :: www    (:,:)
    real(8), allocatable :: fz     (:,:)
    real(8), allocatable :: wp     (:,:)
    real(8), allocatable :: tmp_qp (:,:)
    real(8), allocatable :: irhoadz(:,:)
    real(8), allocatable :: iwmax  (:,:)
    real(8), allocatable :: rhofac (:,:)
    real(8) :: prec_cfl
    real(8) :: eps
    integer :: icol,k,kc,kb
    real(8) :: y,pp,pn
    real(8) :: lat_heat, wmax
    integer nprec, iprec
    real(8) :: flagstat, tmp, tvel

    !Statement functions
    pp(y)= max(real(0.,8),y)
    pn(y)=-min(real(0.,8),y)

    eps = 1.D-10

    allocate( mx     (ncol,nz  ) )
    allocate( mn     (ncol,nz  ) )
    allocate( lfac   (ncol,nz+1) )
    allocate( www    (ncol,nz+1) )
    allocate( fz     (ncol,nz+1) )
    allocate( wp     (ncol,nz  ) )
    allocate( tmp_qp (ncol,nz  ) )
    allocate( irhoadz(ncol,nz) )
    allocate( iwmax  (ncol,nz) )
    allocate( rhofac (ncol,nz) )

    !$acc parallel loop gang vector collapse(2) async(asyncid)
    do k = 1,nz
      do icol = 1 , ncol
        rhofac(icol,k) = sqrt(1.29D0/rho(icol,k))
        irhoadz(icol,k) = 1.D0/(rho(icol,k)*adz(icol,k)) ! Useful factor
        kb = max(1,k-1)
        wmax       = dz(icol)*adz(icol,kb)/dtn   ! Velocity equivalent to a cfl of 1.0.
        iwmax(icol,k)   = 1.D0/wmax
      enddo
    enddo

    !   Add sedimentation of precipitation field to the vert. vel.
    prec_cfl = 0.D0
    !$acc parallel loop gang vector collapse(2) reduction(max:prec_cfl) async(asyncid)
    do k=1,nz
      do icol = 1 , ncol
        lfac(icol,k) = fac_cond + (1-omega(icol,k))*fac_fus
        flagstat = 1.D0
        tvel = term_vel_qp(qp(icol,k),rho(icol,k),tabs(icol,k),qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,&
                           a_gr,vgrau,cgrau,vsnow,csnow)
        wp(icol,k)=rhofac(icol,k)*tvel
        tmp = wp(icol,k)*iwmax(icol,k)
        prec_cfl = max(prec_cfl,tmp) ! Keep column maximum CFL
        wp(icol,k) = -wp(icol,k)*rhow(icol,k)*dtn/dz(icol)
        if (k == 1) then
          fz(icol,nz+1)=0.D0
          www(icol,nz+1)=0.D0
          lfac(icol,nz+1)=0D0
        endif
      enddo
    enddo

    ! If maximum CFL due to precipitation velocity is greater than 0.9,
    ! take more than one advection step to maintain stability.
    if (prec_cfl.gt.0.9D0) then
      nprec = CEILING(prec_cfl/0.9D0)
      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k = 1,nz
        do icol = 1 , ncol
          ! wp already includes factor of dt, so reduce it by a
          ! factor equal to the number of precipitation steps.
          wp(icol,k) = wp(icol,k)/real(nprec,8)
        enddo
      enddo
    else
      nprec = 1
    endif

    !  loop over iterations
    do iprec = 1,nprec
      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k = 1,nz
        do icol = 1 , ncol
          tmp_qp(icol,k) = qp(icol,k) ! Temporary array for qp in this column
        enddo
      enddo

      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k=1,nz
        do icol = 1 , ncol
          kc=min(nz,k+1)
          kb=max(1,k-1)
          mx(icol,k)=max(tmp_qp(icol,kb),tmp_qp(icol,kc),tmp_qp(icol,k))
          mn(icol,k)=min(tmp_qp(icol,kb),tmp_qp(icol,kc),tmp_qp(icol,k))
          ! Define upwind precipitation flux
          fz(icol,k)=tmp_qp(icol,k)*wp(icol,k)
        enddo
      enddo

      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k=1,nz
        do icol = 1 , ncol
          kc=k+1
          tmp_qp(icol,k)=tmp_qp(icol,k)-(fz(icol,kc)-fz(icol,k))*irhoadz(icol,k) !Update temporary qp
        enddo
      enddo

      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k=1,nz
        do icol = 1 , ncol
          ! Also, compute anti-diffusive correction to previous
          ! (upwind) approximation to the flux
          kb=max(1,k-1)
          ! The precipitation velocity is a cell-centered quantity,
          ! since it is computed from the cell-centered
          ! precipitation mass fraction.  Therefore, a reformulated
          ! anti-diffusive flux is used here which accounts for
          ! this and results in reduced numerical diffusion.
          www(icol,k) = 0.5D0*(1.+wp(icol,k)*irhoadz(icol,k))*(tmp_qp(icol,kb)*wp(icol,kb) - &
                                 tmp_qp(icol,k)*wp(icol,k)) ! works for wp(k)<0
        enddo
      enddo

      !---------- non-osscilatory option ---------------
      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k=1,nz
        do icol = 1 , ncol
          kc=min(nz,k+1)
          kb=max(1,k-1)
          mx(icol,k)=max(tmp_qp(icol,kb),tmp_qp(icol,kc),tmp_qp(icol,k),mx(icol,k))
          mn(icol,k)=min(tmp_qp(icol,kb),tmp_qp(icol,kc),tmp_qp(icol,k),mn(icol,k))
          kc=min(nz,k+1)
          mx(icol,k)=rho(icol,k)*adz(icol,k)*(mx(icol,k)-tmp_qp(icol,k))/(pn(www(icol,kc)) + pp(www(icol,k))+eps)
          mn(icol,k)=rho(icol,k)*adz(icol,k)*(tmp_qp(icol,k)-mn(icol,k))/(pp(www(icol,kc)) + pn(www(icol,k))+eps)
        enddo
      enddo
      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k=1,nz
        do icol = 1 , ncol
          kb=max(1,k-1)
          ! Add limited flux correction to fz(k).
          fz(icol,k) = fz(icol,k) + pp(www(icol,k))*min(real(1.,8),mx(icol,k), mn(icol,kb)) - &
                                            pn(www(icol,k))*min(real(1.,8),mx(icol,kb),mn(icol,k)) ! Anti-diffusive flux
        enddo
      enddo

      ! Update precipitation mass fraction and liquid-ice static
      ! energy using precipitation fluxes computed in this column.
      !$acc parallel loop gang vector collapse(2) async(asyncid)
      do k=1,nz
        do icol = 1 , ncol
          kc=k+1
          ! Update precipitation mass fraction.
          ! Note that fz is the total flux, including both the
          ! upwind flux and the anti-diffusive correction.
          qp(icol,k)=qp(icol,k)-(fz(icol,kc)-fz(icol,k))*irhoadz(icol,k)
          tmp = -(fz(icol,kc)-fz(icol,k))*irhoadz(icol,k)*flagstat  ! For qp budget
          qpfall(icol,k)=qpfall(icol,k)+tmp
          lat_heat = -(lfac(icol,kc)*fz(icol,kc)-lfac(icol,k)*fz(icol,k))*irhoadz(icol,k)
          t(icol,k)=t(icol,k)-lat_heat
          tlat(icol,k)=tlat(icol,k)-lat_heat            ! For energy budget
          tmp = fz(icol,k)*flagstat
          precflux(icol,k) = precflux(icol,k) - tmp   ! For statistics
          if (k == 1) then
            precsfc(icol) = precsfc(icol) - fz(icol,1)*flagstat ! For statistics
            precssfc(icol) = precssfc(icol) - fz(icol,1)*(1.-omega(icol,1))*flagstat ! For statistics
            prec_xy(icol) = prec_xy(icol) - fz(icol,1)*flagstat ! For 2D output
          endif
        enddo
      enddo

      if (iprec.lt.nprec) then
        ! Re-compute precipitation velocity using new value of qp.
        !$acc parallel loop gang vector collapse(2) async(asyncid)
        do k=1,nz
          do icol = 1 , ncol
            !Passing variables via first index because of PGI bug with pointers
            tvel = term_vel_qp(qp(icol,k),rho(icol,k),tabs(icol,k),qp_threshold,tprmin,a_pr,vrain,crain,&
                               tgrmin,a_gr,vgrau,cgrau,vsnow,csnow)
            wp(icol,k) = rhofac(icol,k)*tvel
            ! Decrease precipitation velocity by factor of nprec
            wp(icol,k) = -wp(icol,k)*rhow(icol,k)*dtn/dz(icol)/real(nprec,8)
            ! Note: Don't bother checking CFL condition at each
            ! substep since it's unlikely that the CFL will
            ! increase very much between substeps when using
            ! monotonic advection schemes.
            if (k == 1) then
              fz(icol,nz+1)=0.D0
              www(icol,nz+1)=0.D0
              lfac(icol,nz+1)=0.D0
            endif
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
