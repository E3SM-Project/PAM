
module ice_fall_mod
  implicit none

contains

  ! Sedimentation of ice:
  subroutine ice_fall( qcl, qci, tabs, adz, dz, rho, q, t, precsfc, precssfc, dtn, fac_cond, fac_fus, ncol, nz )
    implicit none
    real(8), intent(in   ) :: qcl     (ncol,nz) ! liquid water  (condensate)
    real(8), intent(in   ) :: qci     (ncol,nz) ! ice water  (condensate)
    real(8), intent(in   ) :: tabs    (ncol,nz) ! temperature
    real(8), intent(in   ) :: adz     (ncol,nz) ! ratio of the thickness of scalar levels to dz 
    real(8), intent(in   ) :: dz      (ncol    ) ! constant grid spacing in z direction (when dz_constant=.true.)
    real(8), intent(in   ) :: rho     (ncol,nz) ! air density at pressure levels,kg/m3 
    real(8), intent(inout) :: q       (ncol,nz) ! total nonprecipitating water
    real(8), intent(inout) :: t       (ncol,nz) ! liquid/ice water static energy 
    real(8), intent(inout) :: precsfc (ncol    ) ! surface precip. rate
    real(8), intent(inout) :: precssfc(ncol    ) ! surface ice precip. rate
    real(8), intent(in   ) :: dtn                       ! current dynamical timestep (can be smaller than dt)
    real(8), intent(in   ) :: fac_cond, fac_fus
    integer, intent(in   ) :: ncol, nz

    integer, allocatable :: kmax(:)
    integer, allocatable :: kmin(:)
    real(8), allocatable :: fz(:,:)
    integer :: icol, k, kb, kc
    real(8) coef, dqi, lat_heat, vt_ice
    real(8) qiu, qic, qid, tmp_theta, tmp_phi

    allocate( kmax(ncol) )
    allocate( kmin(ncol) )
    allocate( fz(ncol,nz+1) )

    !$acc parallel loop async(asyncid)
    do icol = 1 , ncol
      kmax(icol) = 0
      kmin(icol) = nz+1
      do k = 1,nz
        if(qcl(icol,k)+qci(icol,k).gt.0.D0.and. tabs(icol,k).lt.273.15D0) then
          kmin(icol) = min(kmin(icol),k)
          kmax(icol) = max(kmax(icol),k)
        endif
      enddo
    enddo

    !$acc parallel loop collapse(2) async(asyncid)
    do k = 1,nz+1
      do icol = 1 , ncol
        fz(icol,k) = 0.D0
      enddo
    enddo

    ! Compute cloud ice flux (using flux limited advection scheme, as in
    ! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
    ! LeVeque, Cambridge University Press, 2002).
    !$acc parallel loop collapse(2) async(asyncid)
    do k = 1 , nz+1
      do icol = 1 , ncol
        if (k >= max(1,kmin(icol)-1) .and. k <= kmax(icol) ) then
          ! Set up indices for x-y planes above and below current plane.
          kc = min(nz,k+1)
          kb = max(1,k-1)
          ! CFL number based on grid spacing interpolated to interface i,j,k-1/2
          coef = dtn/(0.5D0*(adz(icol,kb)+adz(icol,k))*dz(icol))

          ! Compute cloud ice density in this cell and the ones above/below.
          ! Since cloud ice is falling, the above cell is u (upwind),
          ! this cell is c (center) and the one below is d (downwind).
          qiu = rho(icol,kc)*qci(icol,kc)
          qic = rho(icol,k) *qci(icol,k)
          qid = rho(icol,kb)*qci(icol,kb)

          ! Ice sedimentation velocity depends on ice content. The fiting is
          ! based on the data by Heymsfield (JAS,2003). -Marat
          vt_ice = min(real(0.4D0,8),8.66D0*(max(real(0.,8),qic)+1.D-10)**0.24D0)   ! Heymsfield (JAS, 2003, p.2607)

          ! Use MC flux limiter in computation of flux correction.
          ! (MC = monotonized centered difference).
          !         if (qic.eq.qid) then
          if (abs(qic-qid).lt.1.0D-25) then  ! when qic, and qid is very small, qic_qid can still be zero
            ! even if qic is not equal to qid. so add a fix here +++mhwang
            tmp_phi = 0.
          else
            tmp_theta = (qiu-qic)/(qic-qid)
            tmp_phi = max(real(0.,8),min(0.5D0*(1.+tmp_theta),real(2.D0,8),2.D0*tmp_theta))
          endif

          ! Compute limited flux.
          ! Since falling cloud ice is a 1D advection problem, this
          ! flux-limited advection scheme is monotonic.
          fz(icol,k) = -vt_ice*(qic - 0.5D0*(1.D0-coef*vt_ice)*tmp_phi*(qic-qid))
        endif
      enddo
    enddo

    !$acc parallel loop async(asyncid)
    do icol = 1 , ncol
      fz(icol,nz+1) = 0.
    enddo

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1, nz+1
      do icol = 1 , ncol
        if ( k >= max(1,kmin(icol)-2) .and. k <= kmax(icol) ) then
          coef=dtn/(dz(icol)*adz(icol,k)*rho(icol,k))
          ! The cloud ice increment is the difference of the fluxes.
          dqi=coef*(fz(icol,k)-fz(icol,k+1))
          ! Add this increment to both non-precipitating and total water.
          q(icol,k) = q(icol,k) + dqi

          ! The latent heat flux induced by the falling cloud ice enters
          ! the liquid-ice static energy budget in the same way as the
          ! precipitation.  Note: use latent heat of sublimation.
          lat_heat = (fac_cond+fac_fus)*dqi
          ! Add divergence of latent heat flux to liquid-ice static energy.
          t(icol,k)  = t(icol,k)  - lat_heat
        endif
      enddo
    enddo

    !$acc parallel loop async(asyncid)
    do icol = 1 , ncol
      coef=dtn/dz(icol)
      dqi=-coef*fz(icol,1)
      precsfc (icol) = precsfc (icol)+dqi
      precssfc(icol) = precssfc(icol)+dqi
    enddo

    deallocate( kmax )
    deallocate( kmin )
    deallocate( fz )

  end subroutine ice_fall

end module ice_fall_mod
