!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2021
!! @brief   Forcing setting specific to MHD_COMP_ISOTH
!!          'force_common' is inherited
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public get_force

  complex(r8), allocatable, dimension(:,:,:) :: fmx, fmx_old
  complex(r8), allocatable, dimension(:,:,:) :: fmy, fmy_old
  complex(r8), allocatable, dimension(:,:,:) :: fmz, fmz_old

contains


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2021
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine alloc_force
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    allocate(fmx    , source=src)
    allocate(fmx_old, source=src)
    allocate(fmy    , source=src)
    allocate(fmy_old, source=src)
    allocate(fmz    , source=src)
    allocate(fmz_old, source=src)
    deallocate(src)
  end subroutine alloc_force

end module force

