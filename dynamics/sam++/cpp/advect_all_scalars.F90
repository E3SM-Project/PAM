module advect_all_scalars_mod
  use advect_scalar_mod
  use params, only: crm_iknd
  implicit none

contains

  subroutine advect_all_scalars()

    use vars
    use microphysics
    use sgs
    implicit none
    integer(crm_iknd) k,icrm, i, j, kk
    real(crm_rknd), allocatable :: dummy(:,:)

    allocate( dummy(ncrms,nz) )

    !      advection of scalars :
    call advect_scalar(t,dummy,dummy)

    !    Advection of microphysics prognostics:
    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
      .or. doprecip.and.flag_precip(k).eq.1 ) then
        call advect_scalar(micro_field(1,dimx1_s,dimy1_s,1,k),mkadv(1,1,k),mkwle(1,1,k))
      endif
    end do

    !    Advection of sgs prognostics:
    if(dosgs.and.advect_sgs) then
      do k = 1,nsgs_fields
        call advect_scalar(sgs_field(1,dimx1_s,dimy1_s,1,k),dummy,dummy)
      end do
    end if

    !   Precipitation fallout:
    if(doprecip) then
      !do icrm = 1 , ncrms
      !  total_water_prec(icrm) = total_water_prec(icrm) + total_water(ncrms,icrm)
      !enddo
      call micro_precip_fall()
      !do icrm = 1 , ncrms
      !  total_water_prec(icrm) = total_water_prec(icrm) - total_water(ncrms,icrm)
      !enddo
    end if
    deallocate( dummy )

  end subroutine advect_all_scalars

end module advect_all_scalars_mod
