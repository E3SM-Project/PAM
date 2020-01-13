module boundaries_mod
  use periodic_mod
  use task_util_mod
  implicit none

contains

  subroutine boundaries(flag)
    use grid, only: dompi
    implicit none
    integer(crm_iknd), intent(in) :: flag

    call periodic(flag)

  end subroutine boundaries
end module boundaries_mod
