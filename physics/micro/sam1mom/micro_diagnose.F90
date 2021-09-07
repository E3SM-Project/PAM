
module micro_diagnose_mod
  implicit none

contains

  subroutine micro_diagnose(qv, q, qn, tabs, qcl, qci, qpl, qpi, qp, a_bg, a_pr, tbgmin, tprmin, ncrms, nx, ny, nzm)
    implicit none
    real(8), intent(  out) :: qv  (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: q   (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: qn  (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: tabs(ncrms,nx,ny,nzm)
    real(8), intent(  out) :: qcl (ncrms,nx,ny,nzm)
    real(8), intent(  out) :: qci (ncrms,nx,ny,nzm)
    real(8), intent(  out) :: qpl (ncrms,nx,ny,nzm)
    real(8), intent(  out) :: qpi (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: qp  (ncrms,nx,ny,nzm)
    real(8), intent(in   ) :: a_bg, a_pr, tbgmin, tprmin
    integer, intent(in   ) :: ncrms, nx, ny, nzm

    real(8) :: omn, omp
    integer :: i,j,k,icrm

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            qv(icrm,i,j,k) = q(icrm,i,j,k) - qn(icrm,i,j,k)
            omn = max(real(0.,8),min(real(1.,8),(tabs(icrm,i,j,k)-tbgmin)*a_bg))
            qcl(icrm,i,j,k) = qn(icrm,i,j,k)*omn
            qci(icrm,i,j,k) = qn(icrm,i,j,k)*(1.-omn)
            omp = max(real(0.,8),min(real(1.,8),(tabs(icrm,i,j,k)-tprmin)*a_pr))
            qpl(icrm,i,j,k) = qp(icrm,i,j,k)*omp
            qpi(icrm,i,j,k) = qp(icrm,i,j,k)*(1.-omp)
          enddo
        enddo
      enddo
    enddo
  end subroutine micro_diagnose

endmodule micro_diagnose_mod
