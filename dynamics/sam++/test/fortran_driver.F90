
program driver
  use crmdims
  use params, only: crm_rknd, crm_iknd, crm_lknd, r8
  use crm_input_module
  use crm_output_module
  use crm_state_module
  use crm_rad_module
  use dmdf
  use crm_module
  implicit none
  type(crm_input_type)         :: crm_input
  type(crm_output_type)        :: crm_output
  type(crm_state_type)         :: crm_state
  type(crm_rad_type)           :: crm_rad
  integer          , parameter :: plev   = PLEV
  character(len=64), parameter :: prefix = INPUT_FILE
  real(crm_rknd), allocatable  :: lat0  (:)
  real(crm_rknd), allocatable  :: long0 (:)
  real(crm_rknd), allocatable  :: dt_gl (:)
  integer                      :: ncrms, icrm
  real(crm_rknd), allocatable :: read_crm_input_zmid       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_zint       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_tl         (:,:)
  real(crm_rknd), allocatable :: read_crm_input_ql         (:,:)
  real(crm_rknd), allocatable :: read_crm_input_qccl       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_qiil       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_pmid       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_pint       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_pdel       (:,:)
  real(crm_rknd), allocatable :: read_crm_input_ul         (:,:)
  real(crm_rknd), allocatable :: read_crm_input_vl         (:,:)
  real(crm_rknd), allocatable :: read_crm_state_u_wind     (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_v_wind     (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_w_wind     (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_temperature(:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_qt         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_qp         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_state_qn         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qrad         (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_temperature  (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qv           (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qc           (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_qi           (:,:,:,:)
  real(crm_rknd), allocatable :: read_crm_rad_cld          (:,:,:,:)
  character(len=64) :: fprefix = 'fortran_output'

  call dmdf_num_records(prefix,ncrms)
#ifdef NCRMS
  ncrms = NCRMS
#endif

  write(*,*) "File   : ", trim(prefix)
  write(*,*) "Samples: ", ncrms
  write(*,*) "crm_nx : ", crm_nx
  write(*,*) "crm_ny : ", crm_ny
  write(*,*) "crm_dx : ", crm_dx
  write(*,*) "crm_dt : ", crm_dt
  write(*,*) "plev   : ", plev 

  ! Allocate model data
  call crm_input%initialize (           ncrms,plev)
  call crm_output_initialize(crm_output,ncrms,plev)
  ! These are normally allocated by pbuf, so we have to do it explicitly
  allocate( crm_state%u_wind     (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%v_wind     (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%w_wind     (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%temperature(ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%qt         (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%qp         (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_state%qn         (ncrms,crm_nx,crm_ny,crm_nz) )
  allocate( crm_rad%qrad         (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%temperature  (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%qv           (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%qc           (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%qi           (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( crm_rad%cld          (ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
  allocate( lat0                 (ncrms) )
  allocate( long0                (ncrms) )
  allocate( dt_gl                (ncrms) )
  
  ! Allocate transposed arrays because this is the storage format in netcdf
  allocate( read_crm_input_zmid       (plev  ,ncrms))
  allocate( read_crm_input_zint       (plev+1,ncrms))
  allocate( read_crm_input_tl         (plev  ,ncrms))
  allocate( read_crm_input_ql         (plev  ,ncrms))
  allocate( read_crm_input_qccl       (plev  ,ncrms))
  allocate( read_crm_input_qiil       (plev  ,ncrms))
  allocate( read_crm_input_pmid       (plev  ,ncrms))
  allocate( read_crm_input_pint       (plev+1,ncrms))
  allocate( read_crm_input_pdel       (plev  ,ncrms))
  allocate( read_crm_input_ul         (plev  ,ncrms))
  allocate( read_crm_input_vl         (plev  ,ncrms))
  allocate( read_crm_state_u_wind     (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_v_wind     (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_w_wind     (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_temperature(crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_qt         (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_qp         (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_state_qn         (crm_nx    ,crm_ny    ,crm_nz,ncrms) )
  allocate( read_crm_rad_qrad         (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_temperature  (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_qv           (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_qc           (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_qi           (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )
  allocate( read_crm_rad_cld          (crm_nx_rad,crm_ny_rad,crm_nz,ncrms) )

  ! Read in the samples to drive the code
  write(*,*) 'Reading the data'
  call dmdf_read( dt_gl                      , prefix , trim("dt_gl            ") , 1 , ncrms , .true.  , .false. )
  call dmdf_read( lat0                       , prefix , trim("latitude0        ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( long0                      , prefix , trim("longitude0       ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_zmid        , prefix , trim("in_zmid          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_zint        , prefix , trim("in_zint          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_tl          , prefix , trim("in_tl            ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_ql          , prefix , trim("in_ql            ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_qccl        , prefix , trim("in_qccl          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_qiil        , prefix , trim("in_qiil          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_pmid        , prefix , trim("in_pmid          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_pint        , prefix , trim("in_pint          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_pdel        , prefix , trim("in_pdel          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_ul          , prefix , trim("in_ul            ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_input_vl          , prefix , trim("in_vl            ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_state_u_wind      , prefix , trim("state_u_wind     ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_state_v_wind      , prefix , trim("state_v_wind     ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_state_w_wind      , prefix , trim("state_w_wind     ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_state_temperature , prefix , trim("state_temperature") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_state_qt          , prefix , trim("state_qt         ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_state_qp          , prefix , trim("state_qp         ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_state_qn          , prefix , trim("state_qn         ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_rad_qrad          , prefix , trim("rad_qrad         ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_rad_temperature   , prefix , trim("rad_temperature  ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_rad_qv            , prefix , trim("rad_qv           ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_rad_qc            , prefix , trim("rad_qc           ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_rad_qi            , prefix , trim("rad_qi           ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( read_crm_rad_cld           , prefix , trim("rad_cld          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%ps               , prefix , trim("in_ps            ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%phis             , prefix , trim("in_phis          ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%ocnfrac          , prefix , trim("in_ocnfrac       ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%tau00            , prefix , trim("in_tau00         ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%wndls            , prefix , trim("in_wndls         ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%bflxls           , prefix , trim("in_bflxls        ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%fluxu00          , prefix , trim("in_fluxu00       ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%fluxv00          , prefix , trim("in_fluxv00       ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%fluxt00          , prefix , trim("in_fluxt00       ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_input%fluxq00          , prefix , trim("in_fluxq00       ") , 1 , ncrms , .false. , .false. )
  call dmdf_read( crm_output%timing_factor   , prefix , trim("out_timing_factor") , 1 , ncrms , .false. , .true.  )

  do icrm = 1 , ncrms
    crm_input%zmid       (icrm,:)     = read_crm_input_zmid       (:    ,icrm)                       
    crm_input%zint       (icrm,:)     = read_crm_input_zint       (:    ,icrm)                       
    crm_input%tl         (icrm,:)     = read_crm_input_tl         (:    ,icrm)                       
    crm_input%ql         (icrm,:)     = read_crm_input_ql         (:    ,icrm)                       
    crm_input%qccl       (icrm,:)     = read_crm_input_qccl       (:    ,icrm)                       
    crm_input%qiil       (icrm,:)     = read_crm_input_qiil       (:    ,icrm)                       
    crm_input%pmid       (icrm,:)     = read_crm_input_pmid       (:    ,icrm)                       
    crm_input%pint       (icrm,:)     = read_crm_input_pint       (:    ,icrm)                       
    crm_input%pdel       (icrm,:)     = read_crm_input_pdel       (:    ,icrm)                       
    crm_input%ul         (icrm,:)     = read_crm_input_ul         (:    ,icrm)                       
    crm_input%vl         (icrm,:)     = read_crm_input_vl         (:    ,icrm)                       
    crm_state%u_wind     (icrm,:,:,:) = read_crm_state_u_wind     (:,:,:,icrm) 
    crm_state%v_wind     (icrm,:,:,:) = read_crm_state_v_wind     (:,:,:,icrm) 
    crm_state%w_wind     (icrm,:,:,:) = read_crm_state_w_wind     (:,:,:,icrm) 
    crm_state%temperature(icrm,:,:,:) = read_crm_state_temperature(:,:,:,icrm) 
    crm_state%qt         (icrm,:,:,:) = read_crm_state_qt         (:,:,:,icrm) 
    crm_state%qp         (icrm,:,:,:) = read_crm_state_qp         (:,:,:,icrm) 
    crm_state%qn         (icrm,:,:,:) = read_crm_state_qn         (:,:,:,icrm) 
    crm_rad%qrad         (icrm,:,:,:) = read_crm_rad_qrad         (:,:,:,icrm) 
    crm_rad%temperature  (icrm,:,:,:) = read_crm_rad_temperature  (:,:,:,icrm) 
    crm_rad%qv           (icrm,:,:,:) = read_crm_rad_qv           (:,:,:,icrm) 
    crm_rad%qc           (icrm,:,:,:) = read_crm_rad_qc           (:,:,:,icrm) 
    crm_rad%qi           (icrm,:,:,:) = read_crm_rad_qi           (:,:,:,icrm) 
    crm_rad%cld          (icrm,:,:,:) = read_crm_rad_cld          (:,:,:,icrm) 
  enddo

  write(*,*) 'Running the CRM'

  ! Run the code
  call crm( ncrms, dt_gl(1), plev, crm_input, crm_state, crm_rad, crm_output, lat0, long0 )

  write(*,*) 'Writing output data'
  ! dmdf_write(dat,rank,fprefix,vname       ,first,last) !For scalar values
  ! dmdf_write(dat,rank,fprefix,vname,dnames,first,last) !For array values
  do icrm = 1 , ncrms
    call dmdf_write( crm_state%u_wind         (icrm,:,:,:) , 1 , fprefix , trim('state_u_wind        ') , (/'crm_nx','crm_ny','crm_nz'/)             , .true.  , .false. )
    call dmdf_write( crm_state%v_wind         (icrm,:,:,:) , 1 , fprefix , trim('state_v_wind        ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_state%w_wind         (icrm,:,:,:) , 1 , fprefix , trim('state_w_wind        ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_state%temperature    (icrm,:,:,:) , 1 , fprefix , trim('state_temperature   ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_state%qt             (icrm,:,:,:) , 1 , fprefix , trim('state_qt            ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_state%qp             (icrm,:,:,:) , 1 , fprefix , trim('state_qp            ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_state%qn             (icrm,:,:,:) , 1 , fprefix , trim('state_qn            ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%qcl           (icrm,:,:,:) , 1 , fprefix , trim('output_qcl          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%qci           (icrm,:,:,:) , 1 , fprefix , trim('output_qci          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%qpl           (icrm,:,:,:) , 1 , fprefix , trim('output_qpl          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%qpi           (icrm,:,:,:) , 1 , fprefix , trim('output_qpi          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%tk            (icrm,:,:,:) , 1 , fprefix , trim('output_tk           ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%tkh           (icrm,:,:,:) , 1 , fprefix , trim('output_tkh          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%prec_crm      (icrm,:,:)   , 1 , fprefix , trim('output_prec_crm     ') , (/'crm_nx','crm_ny'/)                      , .false. , .false. )
    call dmdf_write( crm_output%wvar          (icrm,:,:,:) , 1 , fprefix , trim('output_wvar         ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%aut           (icrm,:,:,:) , 1 , fprefix , trim('output_aut          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%acc           (icrm,:,:,:) , 1 , fprefix , trim('output_acc          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%evpc          (icrm,:,:,:) , 1 , fprefix , trim('output_evpc         ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%evpr          (icrm,:,:,:) , 1 , fprefix , trim('output_evpr         ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%mlt           (icrm,:,:,:) , 1 , fprefix , trim('output_mlt          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%sub           (icrm,:,:,:) , 1 , fprefix , trim('output_sub          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%dep           (icrm,:,:,:) , 1 , fprefix , trim('output_dep          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%con           (icrm,:,:,:) , 1 , fprefix , trim('output_con          ') , (/'crm_nx','crm_ny','crm_nz'/)             , .false. , .false. )
    call dmdf_write( crm_output%cltot         (icrm)       , 1 , fprefix , trim('output_cltot        ')                                              , .false. , .false. )
    call dmdf_write( crm_output%cllow         (icrm)       , 1 , fprefix , trim('output_cllow        ')                                              , .false. , .false. )
    call dmdf_write( crm_output%clmed         (icrm)       , 1 , fprefix , trim('output_clmed        ')                                              , .false. , .false. )
    call dmdf_write( crm_output%clhgh         (icrm)       , 1 , fprefix , trim('output_clhgh        ')                                              , .false. , .false. )
    call dmdf_write( crm_output%precc         (icrm)       , 1 , fprefix , trim('output_precc        ')                                              , .false. , .false. )
    call dmdf_write( crm_output%precl         (icrm)       , 1 , fprefix , trim('output_precl        ')                                              , .false. , .false. )
    call dmdf_write( crm_output%precsc        (icrm)       , 1 , fprefix , trim('output_precsc       ')                                              , .false. , .false. )
    call dmdf_write( crm_output%precsl        (icrm)       , 1 , fprefix , trim('output_precsl       ')                                              , .false. , .false. )
    call dmdf_write( crm_output%cldtop        (icrm,:)     , 1 , fprefix , trim('output_cldtop       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qc_mean       (icrm,:)     , 1 , fprefix , trim('output_qc_mean      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qi_mean       (icrm,:)     , 1 , fprefix , trim('output_qi_mean      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qs_mean       (icrm,:)     , 1 , fprefix , trim('output_qs_mean      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qg_mean       (icrm,:)     , 1 , fprefix , trim('output_qg_mean      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qr_mean       (icrm,:)     , 1 , fprefix , trim('output_qr_mean      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%sltend        (icrm,:)     , 1 , fprefix , trim('output_sltend       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qltend        (icrm,:)     , 1 , fprefix , trim('output_qltend       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qcltend       (icrm,:)     , 1 , fprefix , trim('output_qcltend      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qiltend       (icrm,:)     , 1 , fprefix , trim('output_qiltend      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%cld           (icrm,:)     , 1 , fprefix , trim('output_cld          ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%gicewp        (icrm,:)     , 1 , fprefix , trim('output_gicewp       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%gliqwp        (icrm,:)     , 1 , fprefix , trim('output_gliqwp       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%mctot         (icrm,:)     , 1 , fprefix , trim('output_mctot        ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%mcup          (icrm,:)     , 1 , fprefix , trim('output_mcup         ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%mcdn          (icrm,:)     , 1 , fprefix , trim('output_mcdn         ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%mcuup         (icrm,:)     , 1 , fprefix , trim('output_mcuup        ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%mcudn         (icrm,:)     , 1 , fprefix , trim('output_mcudn        ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%mu_crm        (icrm,:)     , 1 , fprefix , trim('output_mu_crm       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%md_crm        (icrm,:)     , 1 , fprefix , trim('output_md_crm       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%du_crm        (icrm,:)     , 1 , fprefix , trim('output_du_crm       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%eu_crm        (icrm,:)     , 1 , fprefix , trim('output_eu_crm       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%ed_crm        (icrm,:)     , 1 , fprefix , trim('output_ed_crm       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%jt_crm        (icrm)       , 1 , fprefix , trim('output_jt_crm       ')                                              , .false. , .false. )
    call dmdf_write( crm_output%mx_crm        (icrm)       , 1 , fprefix , trim('output_mx_crm       ')                                              , .false. , .false. )
    call dmdf_write( crm_output%flux_qt       (icrm,:)     , 1 , fprefix , trim('output_flux_qt      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%fluxsgs_qt    (icrm,:)     , 1 , fprefix , trim('output_fluxsgs_qt   ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%tkez          (icrm,:)     , 1 , fprefix , trim('output_tkez         ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%tkesgsz       (icrm,:)     , 1 , fprefix , trim('output_tkesgsz      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%tkz           (icrm,:)     , 1 , fprefix , trim('output_tkz          ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%flux_u        (icrm,:)     , 1 , fprefix , trim('output_flux_u       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%flux_v        (icrm,:)     , 1 , fprefix , trim('output_flux_v       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%flux_qp       (icrm,:)     , 1 , fprefix , trim('output_flux_qp      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%precflux      (icrm,:)     , 1 , fprefix , trim('output_precflux     ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qt_ls         (icrm,:)     , 1 , fprefix , trim('output_qt_ls        ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qt_trans      (icrm,:)     , 1 , fprefix , trim('output_qt_trans     ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qp_trans      (icrm,:)     , 1 , fprefix , trim('output_qp_trans     ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qp_fall       (icrm,:)     , 1 , fprefix , trim('output_qp_fall      ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qp_src        (icrm,:)     , 1 , fprefix , trim('output_qp_src       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%qp_evp        (icrm,:)     , 1 , fprefix , trim('output_qp_evp       ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%t_ls          (icrm,:)     , 1 , fprefix , trim('output_t_ls         ') , (/'nlev'/)                                 , .false. , .false. )
    call dmdf_write( crm_output%prectend      (icrm)       , 1 , fprefix , trim('output_prectend     ')                                              , .false. , .false. )
    call dmdf_write( crm_output%precstend     (icrm)       , 1 , fprefix , trim('output_precstend    ')                                              , .false. , .false. )
    call dmdf_write( crm_output%taux          (icrm)       , 1 , fprefix , trim('output_taux         ')                                              , .false. , .false. )
    call dmdf_write( crm_output%tauy          (icrm)       , 1 , fprefix , trim('output_tauy         ')                                              , .false. , .false. )
    call dmdf_write( crm_output%z0m           (icrm)       , 1 , fprefix , trim('output_z0m          ')                                              , .false. , .false. )
    call dmdf_write( crm_output%timing_factor (icrm)       , 1 , fprefix , trim('output_timing_factor')                                              , .false. , .false. )
    call dmdf_write( crm_rad%qrad             (icrm,:,:,:) , 1 , fprefix , trim('rad_qrad            ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
    call dmdf_write( crm_rad%temperature      (icrm,:,:,:) , 1 , fprefix , trim('rad_temperature     ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
    call dmdf_write( crm_rad%qv               (icrm,:,:,:) , 1 , fprefix , trim('rad_qv              ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
    call dmdf_write( crm_rad%qc               (icrm,:,:,:) , 1 , fprefix , trim('rad_qc              ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
    call dmdf_write( crm_rad%qi               (icrm,:,:,:) , 1 , fprefix , trim('rad_qi              ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .false. )
    call dmdf_write( crm_rad%cld              (icrm,:,:,:) , 1 , fprefix , trim('rad_cld             ') , (/'crm_nx_rad','crm_ny_rad','crm_nz    '/) , .false. , .true.  )
  enddo

end program driver


