
module abcoefs_mod
  implicit none

contains

  subroutine abcoefs(dt3, nstep, na, nb, nc, at, bt, ct) bind(C, name="abcoefs_f90")
    ! Coefficients for the Adams-Bashforth scheme
    use params, only: crm_rknd, crm_iknd
    implicit none
    real(crm_rknd)          , intent(in   ) :: dt3(3)
    integer(crm_iknd), value, intent(in   ) :: na, nb, nc, nstep
    real(crm_rknd)          , intent(  out) :: at, bt, ct
    real(crm_rknd) :: alpha, beta

    if(nstep.ge.3) then
      alpha = dt3(nb) / dt3(na)
      beta = dt3(nc) / dt3(na)
      ct = (2.+3.* alpha) / (6.* (alpha + beta) * beta)
      bt = -(1.+2.*(alpha + beta) * ct)/(2. * alpha)
      at = 1. - bt - ct
    else if(nstep.ge.2) then
      at = 3./2.
      bt = -1./2.
      ct = 0.
    else
      at = 1.
      bt = 0.
      ct = 0.
    end if

  end subroutine abcoefs

end module abcoefs_mod

