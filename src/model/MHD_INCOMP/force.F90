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

    allocate(fux     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fux      = 0.d0
    allocate(fux_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fux_old1 = 0.d0
    allocate(fuy     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fuy      = 0.d0
    allocate(fuy_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fuy_old1 = 0.d0
    allocate(fuz     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fuz      = 0.d0
    allocate(fuz_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fuz_old1 = 0.d0
  end subroutine alloc_force

end module force
