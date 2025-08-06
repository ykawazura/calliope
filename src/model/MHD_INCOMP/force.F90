!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../force_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2021
!! @brief   Forcing setting specific to MHD_INCOMP
!!          'force_common' is inherited
!-----------------------------------------------!
module force
  use force_common
  implicit none

  public get_force

  complex(r8), allocatable, dimension(:,:,:) :: fux, fux_old
  complex(r8), allocatable, dimension(:,:,:) :: fuy, fuy_old
  complex(r8), allocatable, dimension(:,:,:) :: fuz, fuz_old
  complex(r8), allocatable, dimension(:,:,:) :: fbx, fbx_old
  complex(r8), allocatable, dimension(:,:,:) :: fby, fby_old
  complex(r8), allocatable, dimension(:,:,:) :: fbz, fbz_old
  complex(r8), allocatable, dimension(:,:,:) :: fzpx
  complex(r8), allocatable, dimension(:,:,:) :: fzpy
  complex(r8), allocatable, dimension(:,:,:) :: fzpz
  complex(r8), allocatable, dimension(:,:,:) :: fzmx
  complex(r8), allocatable, dimension(:,:,:) :: fzmy
  complex(r8), allocatable, dimension(:,:,:) :: fzmz

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
    allocate(fux    , source=src)
    allocate(fux_old, source=src)
    allocate(fuy    , source=src)
    allocate(fuy_old, source=src)
    allocate(fuz    , source=src)
    allocate(fuz_old, source=src)
    allocate(fbx    , source=src)
    allocate(fbx_old, source=src)
    allocate(fby    , source=src)
    allocate(fby_old, source=src)
    allocate(fbz    , source=src)
    allocate(fbz_old, source=src)
    if(elsasser) then
      allocate(fzpx, source=src)
      allocate(fzpy, source=src)
      allocate(fzpz, source=src)
      allocate(fzmx, source=src)
      allocate(fzmy, source=src)
      allocate(fzmz, source=src)
    endif
    deallocate(src)
  end subroutine alloc_force


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2021
!! @brief   Compute the power of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force(fux, fuy, fuz, fbx, fby, fbz)
    use grid, only: nkx, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: fux, fuy, fuz, fbx, fby, fbz
    real(r8) :: p_ext_ene_sum
    real(r8), allocatable, dimension(:,:,:) :: p_ext_ene

    integer :: i, j, k, istir

    if (proc0) call put_time_stamp(timer_force)

    allocate(p_ext_ene(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    
    if (fix_power) then
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            p_ext_ene(i, k, j) = 0.5d0*( &
                                    (fux(i, k, j)*conjg(ux(i, k, j)) + conjg(fux(i, k, j))*ux(i, k, j)) &
                                  + (fuy(i, k, j)*conjg(uy(i, k, j)) + conjg(fuy(i, k, j))*uy(i, k, j)) &
                                  + (fuz(i, k, j)*conjg(uz(i, k, j)) + conjg(fuz(i, k, j))*uz(i, k, j)) &
                                  + (fbx(i, k, j)*conjg(bx(i, k, j)) + conjg(fbx(i, k, j))*bx(i, k, j)) &
                                  + (fby(i, k, j)*conjg(by(i, k, j)) + conjg(fby(i, k, j))*by(i, k, j)) &
                                  + (fbz(i, k, j)*conjg(bz(i, k, j)) + conjg(fbz(i, k, j))*bz(i, k, j)) &
                                )
            ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
            !-----------------------------------------------------------------------------------------------------------------------------------
            ! The volume integral of a quadratic function is
            ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
            ! |f_{k_x, k_y}|^2 = \left( \sum_{k_y = -n_{k_y}/2}^{-1}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 0}
            ! \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 1}^{n_{k_y}/2}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} \right) |f_{k_x, k_y}|^2\]
            ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
            !-----------------------------------------------------------------------------------------------------------------------------------
            if (j /= 1) then
              p_ext_ene(i, k, j) = 2.0d0*p_ext_ene(i, k, j)
            endif

          end do
        end do
      end do

      p_ext_ene_sum = sum(p_ext_ene); call sum_allreduce(p_ext_ene_sum)

      if(p_ext_ene_sum > 1.0d-6) then
        fux = ene_inj*fux/p_ext_ene_sum
        fuy = ene_inj*fuy/p_ext_ene_sum
        fuz = ene_inj*fuz/p_ext_ene_sum
        fbx = ene_inj*fbx/p_ext_ene_sum
        fby = ene_inj*fby/p_ext_ene_sum
        fbz = ene_inj*fbz/p_ext_ene_sum
      endif
    endif

    deallocate(p_ext_ene)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2021
!! @brief   Compute the power and helicity of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
    use grid, only: nkx, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: fzpx, fzpy, fzpz, fzmx, fzmy, fzmz
    complex(r8) :: zpx, zpy, zpz, zmx, zmy, zmz
    real(r8)    :: zp_dot_fzp_sum, zm_dot_fzm_sum
    real(r8), allocatable, dimension(:,:,:) :: zp_dot_fzp, zm_dot_fzm

    integer :: i, j, k, istir

    if (proc0) call put_time_stamp(timer_force)

    allocate(zp_dot_fzp(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(zm_dot_fzm(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    
    if (fix_power) then
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            zpx = ux(i, k, j) + bx(i, k, j)
            zpy = uy(i, k, j) + by(i, k, j)
            zpz = uz(i, k, j) + bz(i, k, j)
            zmx = ux(i, k, j) - bx(i, k, j)
            zmy = uy(i, k, j) - by(i, k, j)
            zmz = uz(i, k, j) - bz(i, k, j)

            zp_dot_fzp(i, k, j) = 0.5d0*( &
                                    (fzpx(i, k, j)*conjg(zpx) + conjg(fzpx(i, k, j))*zpx) &
                                  + (fzpy(i, k, j)*conjg(zpy) + conjg(fzpy(i, k, j))*zpy) &
                                  + (fzpz(i, k, j)*conjg(zpz) + conjg(fzpz(i, k, j))*zpz) &
                                )
            zm_dot_fzm(i, k, j) = 0.5d0*( &
                                    (fzmx(i, k, j)*conjg(zmx) + conjg(fzmx(i, k, j))*zmx) &
                                  + (fzmy(i, k, j)*conjg(zmy) + conjg(fzmy(i, k, j))*zmy) &
                                  + (fzmz(i, k, j)*conjg(zmz) + conjg(fzmz(i, k, j))*zmz) &
                                )
            ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
            !-----------------------------------------------------------------------------------------------------------------------------------
            ! The volume integral of a quadratic function is
            ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
            ! |f_{k_x, k_y}|^2 = \left( \sum_{k_y = -n_{k_y}/2}^{-1}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 0}
            ! \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 1}^{n_{k_y}/2}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} \right) |f_{k_x, k_y}|^2\]
            ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
            !-----------------------------------------------------------------------------------------------------------------------------------
            if (j /= 1) then
              zp_dot_fzp(i, k, j) = 2.0d0*zp_dot_fzp(i, k, j)
              zm_dot_fzm(i, k, j) = 2.0d0*zm_dot_fzm(i, k, j)
            endif

          end do
        end do
      end do

      zp_dot_fzp_sum = sum(zp_dot_fzp); call sum_allreduce(zp_dot_fzp_sum)
      zm_dot_fzm_sum = sum(zm_dot_fzm); call sum_allreduce(zm_dot_fzm_sum)

      if(abs(zp_dot_fzp_sum) > 1.0d-6 .and. abs(zm_dot_fzm_sum) > 1.0d-6) then
        fzpx = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzpx/abs(zp_dot_fzp_sum)
        fzpy = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzpy/abs(zp_dot_fzp_sum)
        fzpz = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzpz/abs(zp_dot_fzp_sum)
        fzmx = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmx/abs(zm_dot_fzm_sum)
        fzmy = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmy/abs(zm_dot_fzm_sum)
        fzmz = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmz/abs(zm_dot_fzm_sum)
      endif
    endif

    deallocate(zp_dot_fzp)
    deallocate(zm_dot_fzm)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force_els

end module force

