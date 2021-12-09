!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2021
!! @brief   Forcing setting specific to RMHD
!!          'force_common' is inherited
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public get_force

  complex(r8), allocatable, dimension(:,:,:) :: fmx, fmx_old1
  complex(r8), allocatable, dimension(:,:,:) :: fmy, fmy_old1
  complex(r8), allocatable, dimension(:,:,:) :: fmz, fmz_old1

contains


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2021
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine alloc_force
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none

    allocate(fmx     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmx      = 0.d0
    allocate(fmx_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmx_old1 = 0.d0
    allocate(fmy     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmy      = 0.d0
    allocate(fmy_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmy_old1 = 0.d0
    allocate(fmz     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmz      = 0.d0
    allocate(fmz_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmz_old1 = 0.d0
  end subroutine alloc_force

end module force

