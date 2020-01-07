module advect_mom_mod
  use advect2_mom_xy_mod
  use advect2_mom_z_mod
  implicit none

contains

  subroutine advect_mom(ncrms)
    use vars
    use params, only: docolumn
    implicit none
    integer(crm_iknd), intent(in) :: ncrms
    integer(crm_iknd) i,j,k

    if(docolumn) return

    call advect2_mom_xy(ncrms)
    call advect2_mom_z(ncrms)

  end subroutine advect_mom

end module advect_mom_mod
