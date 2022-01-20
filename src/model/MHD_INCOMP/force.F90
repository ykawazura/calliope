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

  complex(r8), allocatable, dimension(:,:,:) :: fux, fux_old1
  complex(r8), allocatable, dimension(:,:,:) :: fuy, fuy_old1
  complex(r8), allocatable, dimension(:,:,:) :: fuz, fuz_old1

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
    allocate(fux     , source=src)
    allocate(fux_old1, source=src)
    allocate(fuy     , source=src)
    allocate(fuy_old1, source=src)
    allocate(fuz     , source=src)
    allocate(fuz_old1, source=src)
    deallocate(src)
  end subroutine alloc_force

end module force

