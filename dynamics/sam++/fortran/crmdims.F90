
#ifndef CRM_NX
#define CRM_NX 32
#endif

#ifndef CRM_NY
#define CRM_NY 1
#endif

#ifndef CRM_NZ
#define CRM_NZ 28
#endif

#ifndef CRM_NX_RAD
#define CRM_NX_RAD 2
#endif

#ifndef CRM_NY_RAD
#define CRM_NY_RAD 2
#endif

#ifndef CRM_DX
#define CRM_DX 1000
#endif

#ifndef CRM_DT
#define CRM_DT 5
#endif

module crmdims
! whannah - non-SP not compiling without this #ifdef - says it can't find params.mod - not sure how to fix
    use params, only: crm_rknd
    implicit none
    integer, parameter ::  crm_nx=CRM_NX
    integer, parameter ::  crm_ny=CRM_NY
    integer, parameter ::  crm_nz=CRM_NZ

    integer, parameter ::  crm_nx_rad=CRM_NX_RAD
    integer, parameter ::  crm_ny_rad=CRM_NY_RAD

    real(crm_rknd), parameter :: crm_dx=CRM_DX
    real(crm_rknd), parameter :: crm_dy=crm_dx
    real(crm_rknd), parameter :: crm_dt=CRM_DT

end module crmdims
