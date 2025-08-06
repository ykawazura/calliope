!-----------------------------------------------!
!> @author  YK
!! @date    26 Feb 2024
!! @brief   Model specific subroutines
!-----------------------------------------------!
module model_specific
  implicit none

  public init_model_specific

contains


!-----------------------------------------------!
!> @author  YK
!! @date    26 Feb 2024
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine init_model_specific
    use mp, only: proc0
    use grid, only: kx, ky
    implicit none

    ! if(proc0 .and. maxval(kx**2) /= maxval(ky**2)) then
      ! print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      ! print *, '!              Error!              !'
      ! print *, '!  kx_max must be equal to ky_max  !'
      ! print *, '!  i.e., lx/nlx == ly/nly          !'
      ! print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      ! stop
    ! endif


  end subroutine init_model_specific
end module model_specific

