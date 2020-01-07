
program unit_abcoefs
  use abcoefs
  use params, only: crm_rknd
  implicit none
  real(crm_rknd) :: dt3(3), at, bt, ct
  integer        :: na, nb, nc

  dt3(1) = 4.2
  dt3(2) = 8.4
  dt3(3) = 5.1

  ! Test 1
  na = 1 ; nb = 2 ; nc = 3
  call abcoefs(dt3, na, nb, nc, at, bt, ct)
  write(*,*) 'abcoefs Fortran Test 1'
  write(*,*) 'at: ', at
  write(*,*) 'bt: ', bt
  write(*,*) 'ct: ', ct
  write(*,*)

  ! Test 2
  na = 3 ; nb = 1 ; nc = 2
  call abcoefs(dt3, na, nb, nc, at, bt, ct)
  write(*,*) 'abcoefs Fortran Test 2'
  write(*,*) 'at: ', at
  write(*,*) 'bt: ', bt
  write(*,*) 'ct: ', ct
  write(*,*)

  ! Test 3
  na = 2 ; nb = 3 ; nc = 1
  call abcoefs(dt3, na, nb, nc, at, bt, ct)
  write(*,*) 'abcoefs Fortran Test 3'
  write(*,*) 'at: ', at
  write(*,*) 'bt: ', bt
  write(*,*) 'ct: ', ct
  write(*,*)

end program unit_abcoefs


