
module micro_diagnose_mod
  implicit none

contains

  subroutine micro_diagnose(qv, q, qn, tabs, qcl, qci, qpl, qpi, qp, a_bg, a_pr, tbgmin, tprmin, ncol, nzm)
    implicit none
    real(8), intent(  out) :: qv  (ncol,nzm) ! water vapor
    real(8), intent(in   ) :: q   (ncol,nzm) ! total nonprecipitating water
    real(8), intent(in   ) :: qn  (ncol,nzm) ! cloud condensate (liquid + ice)
    real(8), intent(in   ) :: tabs(ncol,nzm) ! temperature
    real(8), intent(  out) :: qcl (ncol,nzm) ! liquid water  (condensate)
    real(8), intent(  out) :: qci (ncol,nzm) ! ice water  (condensate)
    real(8), intent(  out) :: qpl (ncol,nzm) ! liquid water  (precipitation)
    real(8), intent(  out) :: qpi (ncol,nzm) ! ice water  (precipitation)
    real(8), intent(in   ) :: qp  (ncol,nzm) ! total precipitating water
    real(8), intent(in   ) :: a_bg, a_pr, tbgmin, tprmin
    integer, intent(in   ) :: ncol, nzm

    real(8) :: omn, omp
    integer :: icol,k

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icol = 1 , ncol
        qv(icol,k) = q(icol,k) - qn(icol,k)
        omn = max(real(0.,8),min(real(1.,8),(tabs(icol,k)-tbgmin)*a_bg))
        qcl(icol,k) = qn(icol,k)*omn
        qci(icol,k) = qn(icol,k)*(1.-omn)
        omp = max(real(0.,8),min(real(1.,8),(tabs(icol,k)-tprmin)*a_pr))
        qpl(icol,k) = qp(icol,k)*omp
        qpi(icol,k) = qp(icol,k)*(1.-omp)
      enddo
    enddo
  end subroutine micro_diagnose

endmodule micro_diagnose_mod
