module grid

  use domain
  use params, only: crm_rknd, crm_iknd, crm_lknd

  implicit none

  character(6), parameter :: version = '6.10.4'
  character(8), parameter :: version_date = 'Feb 2013'

  integer(crm_iknd), parameter :: nx = nx_gl/nsubdomains_x
  integer(crm_iknd), parameter :: ny = ny_gl/nsubdomains_y
  integer(crm_iknd), parameter :: nz = nz_gl+1        ! note that nz_gl = crm_nz
  integer(crm_iknd), parameter :: nzm = nz-1          ! note that nzm   = crm_nz

  integer(crm_iknd), parameter :: nsubdomains = nsubdomains_x * nsubdomains_y

  logical(crm_lknd), parameter :: RUN3D = ny_gl.gt.1
  logical(crm_lknd), parameter :: RUN2D = .not.RUN3D

  integer(crm_iknd), parameter :: nxp1 = nx + 1
  integer(crm_iknd), parameter :: nyp1 = ny + 1 * YES3D
  integer(crm_iknd), parameter :: nxp2 = nx + 2
  integer(crm_iknd), parameter :: nyp2 = ny + 2 * YES3D
  integer(crm_iknd), parameter :: nxp3 = nx + 3
  integer(crm_iknd), parameter :: nyp3 = ny + 3 * YES3D
  integer(crm_iknd), parameter :: nxp4 = nx + 4
  integer(crm_iknd), parameter :: nyp4 = ny + 4 * YES3D

  integer(crm_iknd), parameter :: dimx1_u = -1                !!-1        -1        -1        -1
  integer(crm_iknd), parameter :: dimx2_u = nxp3              !!nxp3      nxp3      nxp3      nxp3
  integer(crm_iknd), parameter :: dimy1_u = 1-(2)*YES3D  !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-2*YES3D
  integer(crm_iknd), parameter :: dimy2_u = nyp2         !!nyp5      nyp4      nyp3      nyp2
  integer(crm_iknd), parameter :: dimx1_v = -1           !!-4        -3        -2        -1
  integer(crm_iknd), parameter :: dimx2_v = nxp2         !!nxp5      nxp4      nxp3      nxp2
  integer(crm_iknd), parameter :: dimy1_v = 1-2*YES3D         !!1-2*YES3D 1-2*YES3D 1-2*YES3D 1-2*YES3D
  integer(crm_iknd), parameter :: dimy2_v = nyp3              !!nyp3       nyp3      nyp3      nyp3
  integer(crm_iknd), parameter :: dimx1_w = -1           !!-4        -3        -2        -1
  integer(crm_iknd), parameter :: dimx2_w = nxp2         !!nxp5      nxp4      nxp3      nxp2
  integer(crm_iknd), parameter :: dimy1_w = 1-(2)*YES3D  !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-2*YES3D
  integer(crm_iknd), parameter :: dimy2_w = nyp2         !!nyp5      nyp4      nyp3      nyp2
  integer(crm_iknd), parameter :: dimx1_s = -2          !!-4        -3        -2        -2
  integer(crm_iknd), parameter :: dimx2_s = nxp3        !!nxp5      nxp4      nxp3      nxp3
  integer(crm_iknd), parameter :: dimy1_s = 1-(3)*YES3D !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-3*YES3D
  integer(crm_iknd), parameter :: dimy2_s = nyp3        !!nyp5      nyp4      nyp3      nyp3

  integer(crm_iknd), parameter :: ncols = nx*ny
  integer(crm_iknd), parameter :: nadams = 3

  ! Vertical grid parameters:
  ! real(crm_rknd) pres0      ! Reference surface pressure, Pa

  integer(crm_iknd) :: nstep = 0! current number of performed time steps
  integer(crm_iknd)  ncycle  ! number of subcycles over the dynamical timestep
  integer(crm_iknd) icycle  ! current subcycle
  integer(crm_iknd) :: na, nb, nc ! indices for swapping the rhs arrays for AB scheme
  real(crm_rknd) at, bt, ct ! coefficients for the Adams-Bashforth scheme
  real(crm_rknd) dtn  ! current dynamical timestep (can be smaller than dt)
  real(crm_rknd) dtfactor   ! dtn/dt

  !  MPI staff:
  integer(crm_iknd) rank   ! rank of the current subdomain task (default 0)
  integer(crm_iknd) ranknn ! rank of the "northern" subdomain task
  integer(crm_iknd) rankss ! rank of the "southern" subdomain task
  integer(crm_iknd) rankee ! rank of the "eastern"  subdomain task
  integer(crm_iknd) rankww ! rank of the "western"  subdomain task
  integer(crm_iknd) rankne ! rank of the "north-eastern" subdomain task
  integer(crm_iknd) ranknw ! rank of the "north-western" subdomain task
  integer(crm_iknd) rankse ! rank of the "south-eastern" subdomain task
  integer(crm_iknd) ranksw ! rank of the "south-western" subdomain task
  logical(crm_lknd) dompi  ! logical(crm_lknd) switch to do multitasking
  logical(crm_lknd) masterproc ! .true. if rank.eq.0

  character(80) case   ! id-string to identify a case-name(set in CaseName file)

  logical(crm_lknd) dostatis     ! flag to permit the gathering of statistics
  logical(crm_lknd) dostatisrad  ! flag to permit the gathering of radiation statistics
  integer(crm_iknd) nstatis ! the interval between substeps to compute statistics

  logical(crm_lknd) :: compute_reffc = .false.
  logical(crm_lknd) :: compute_reffi = .false.

  logical(crm_lknd) notopened2D  ! flag to see if the 2D output datafile is opened
  logical(crm_lknd) notopened3D  ! flag to see if the 3D output datafile is opened
  logical(crm_lknd) notopenedmom ! flag to see if the statistical moment file is opened

  !-----------------------------------------
  ! Parameters controled by namelist PARAMETERS

  real(crm_rknd) :: dx =0.   ! grid spacing in x direction
  real(crm_rknd) :: dy =0.   ! grid spacing in y direction
  real(crm_rknd), allocatable :: dz(:)    ! constant grid spacing in z direction (when dz_constant=.true.)
  logical(crm_lknd):: doconstdz = .false.  ! do constant vertical grid spacing set by dz

  integer(crm_iknd):: nstop =0   ! time step number to stop the integration
  integer(crm_iknd):: nelapse =999999999! time step number to elapse before stoping

  real(crm_rknd):: dt=0.  ! dynamical timestep
  real(crm_rknd):: day0=0.  ! starting day (including fraction)

  integer(crm_iknd):: nrad =1  ! frequency of calling the radiation routines
  integer(crm_iknd):: nrestart =0 ! switch to control starting/restarting of the model
  integer(crm_iknd):: nstat =1000 ! the interval in time steps to compute statistics
  integer(crm_iknd):: nstatfrq =50 ! frequency of computing statistics

  logical(crm_lknd):: restart_sep =.false.  ! write separate restart files for sub-domains
  integer(crm_iknd):: nrestart_skip =0 ! number of skips of writing restart (default 0)
  logical(crm_lknd):: output_sep =.false.   ! write separate 3D and 2D files for sub-domains

  character(80):: caseid =''! id-string to identify a run
  character(80):: caseid_restart =''! id-string for branch restart file
  character(80):: case_restart =''! id-string for branch restart file

  logical(crm_lknd):: doisccp = .false.
  logical(crm_lknd):: domodis = .false.
  logical(crm_lknd):: domisr = .false.
  logical(crm_lknd):: dosimfilesout = .false.

  logical(crm_lknd):: doSAMconditionals = .false. !core updraft,downdraft conditional statistics
  logical(crm_lknd):: dosatupdnconditionals = .false.!cloudy updrafts,downdrafts and cloud-free
  logical(crm_lknd):: doscamiopdata = .false.! initialize the case from a SCAM IOP netcdf input file
  logical(crm_lknd):: dozero_out_day0 = .false.
  character(len=120):: iopfile=''
  character(256):: rundatadir ='./RUNDATA' ! path to data directory

  integer(crm_iknd):: nsave3Dstart =99999999! timestep to start writting 3D fields
  integer(crm_iknd):: nsave3Dend  =99999999 ! timestep to end writting 3D fields
  logical(crm_lknd):: save3Dbin =.false.   ! save 3D data in binary format(no 2-byte compression)
  logical(crm_lknd):: save3Dsep =.false.   ! use separate file for each time point for2-model
  real(crm_rknd)   :: qnsave3D =0.    !threshold manimum cloud water(kg/kg) to save 3D fields
  logical(crm_lknd):: dogzip3D =.false.    ! gzip compress a 3D output file
  logical(crm_lknd):: rad3Dout = .false. ! output additional 3D radiation foelds (like reff)

  integer(crm_iknd):: nsave2D =1000     ! frequency of writting 2D fields (steps)
  integer(crm_iknd):: nsave2Dstart =99999999! timestep to start writting 2D fields
  integer(crm_iknd):: nsave2Dend =99999999  ! timestep to end writting 2D fields
  logical(crm_lknd):: save2Dbin =.false.   ! save 2D data in binary format, rather than compressed
  logical(crm_lknd):: save2Dsep =.false.   ! write separate file for each time point for 2D output
  logical(crm_lknd):: save2Davg =.false.   ! flag to time-average 2D output fields (default .false.)
  logical(crm_lknd):: dogzip2D =.false.    ! gzip compress a 2D output file if save2Dsep=.true.

  integer(crm_iknd):: nstatmom =1000! frequency of writting statistical moment fields (steps)
  integer(crm_iknd):: nstatmomstart =99999999! timestep to start writting statistical moment fields
  integer(crm_iknd):: nstatmomend =99999999  ! timestep to end writting statistical moment fields
  logical(crm_lknd):: savemomsep =.false.! use one file with stat moments for each time point
  logical(crm_lknd):: savemombin =.false.! save statistical moment data in binary format

  integer(crm_iknd):: nmovie =1000! frequency of writting movie fields (steps)
  integer(crm_iknd):: nmoviestart =99999999! timestep to start writting statistical moment fields
  integer(crm_iknd):: nmovieend =99999999  ! timestep to end writting statistical moment fields

  logical(crm_lknd) :: isInitialized_scamiopdata = .false.
  logical(crm_lknd) :: wgls_holds_omega = .false.

  real(crm_rknd), allocatable :: z    (:,:)      ! height of the pressure levels above surface,m
  real(crm_rknd), allocatable :: pres (:,:)  ! pressure,mb at scalar levels
  real(crm_rknd), allocatable :: zi   (:,:)     ! height of the interface levels
  real(crm_rknd), allocatable :: presi(:,:)  ! pressure,mb at interface levels
  real(crm_rknd), allocatable :: adz  (:,:)   ! ratio of the thickness of scalar levels to dz
  real(crm_rknd), allocatable :: adzw (:,:)  ! ratio of the thinckness of w levels to dz
  real(crm_rknd), allocatable :: dt3  (:)   ! dynamical timesteps for three most recent time steps

  !-----------------------------------------

contains



  subroutine allocate_grid()
    implicit none
    real(crm_rknd) :: zero

    allocate( z(ncrms,nz)       )
    allocate( pres(ncrms,nzm)   )
    allocate( zi(ncrms,nz)      )
    allocate( presi(ncrms,nz)   )
    allocate( adz(ncrms,nzm)    )
    allocate( adzw(ncrms,nz)    )
    allocate( dt3(3)      )
    allocate( dz(ncrms)         )

    zero = 0

    na=1
    nb=2
    nc=3
    z = zero
    pres = zero
    zi = zero
    presi = zero
    adz = zero
    adzw = zero
    dt3 = zero
    dz = zero
  end subroutine allocate_grid


  subroutine deallocate_grid()
    implicit none
    deallocate( z )
    deallocate( pres )
    deallocate( zi )
    deallocate( presi )
    deallocate( adz )
    deallocate( adzw )
    deallocate( dt3 )
    deallocate( dz )
  end subroutine deallocate_grid


end module grid
