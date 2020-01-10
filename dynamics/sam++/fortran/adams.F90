module adams_mod
  use params, only: asyncid
  implicit none

contains

  subroutine adams(ncrms, dtn, dx, dy, rho, rhow, dudt, dvdt, dwdt, misc, u, v, w, dt3)
    ! Adams-Bashforth scheme
    use vars
    use params, only: crm_rknd
    implicit none
    integer(crm_iknd), intent(in) :: ncrms
    real(crm_rknd) :: dtn, dx, dy
    real(crm_rknd) :: u(ncrms,dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm)
    real(crm_rknd) :: v(ncrms,dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm)
    real(crm_rknd) :: w(ncrms,dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz )
    real(crm_rknd) :: dudt(ncrms,nxp1, ny, nzm, 3)
    real(crm_rknd) :: dvdt(ncrms,nx, nyp1, nzm, 3)
    real(crm_rknd) :: dwdt(ncrms,nx, ny  , nz,  3)
    real(crm_rknd) :: misc(ncrms,nx, ny, nz)
    real(crm_rknd) :: dt3(3)
    real(crm_rknd) :: rho (ncrms,nzm)
    real(crm_rknd) :: rhow(ncrms,nz )
    real(crm_rknd) dtdx, dtdy, dtdz, rhox, rhoy, rhoz
    integer(crm_iknd) i,j,k,icrm

    dtdx = dtn/dx
    dtdy = dtn/dy

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dtdz = dtn/dz(icrm)
            rhox = rho(icrm,k)*dtdx
            rhoy = rho(icrm,k)*dtdy
            rhoz = rhow(icrm,k)*dtdz
            dudt(icrm,i,j,k,nc) = u(icrm,i,j,k) + dt3(na) *(at*dudt(icrm,i,j,k,na)+bt*dudt(icrm,i,j,k,nb)+ct*dudt(icrm,i,j,k,nc))
            dvdt(icrm,i,j,k,nc) = v(icrm,i,j,k) + dt3(na) *(at*dvdt(icrm,i,j,k,na)+bt*dvdt(icrm,i,j,k,nb)+ct*dvdt(icrm,i,j,k,nc))
            dwdt(icrm,i,j,k,nc) = w(icrm,i,j,k) + dt3(na) *(at*dwdt(icrm,i,j,k,na)+bt*dwdt(icrm,i,j,k,nb)+ct*dwdt(icrm,i,j,k,nc))
            u(icrm,i,j,k) = 0.5*(u(icrm,i,j,k)+dudt(icrm,i,j,k,nc)) * rhox
            v(icrm,i,j,k) = 0.5*(v(icrm,i,j,k)+dvdt(icrm,i,j,k,nc)) * rhoy
            misc(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc))
            w(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc)) * rhoz
          end do
        end do
      end do
    end do

  end subroutine adams

end module adams_mod
