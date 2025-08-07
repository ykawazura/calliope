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

  complex(r8), allocatable, dimension(:,:,:) :: fphi, fphi_old
  complex(r8), allocatable, dimension(:,:,:) :: fpsi, fpsi_old
  complex(r8), allocatable, dimension(:,:,:) :: fzppe
  complex(r8), allocatable, dimension(:,:,:) :: fzmpe

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
    allocate(fphi    , source=src)
    allocate(fphi_old, source=src)
    allocate(fpsi    , source=src)
    allocate(fpsi_old, source=src)
    if(elsasser) then
      allocate(fzppe, source=src)
      allocate(fzmpe, source=src)
    endif
    deallocate(src)
  end subroutine alloc_force


!-----------------------------------------------!
!> @author  YK
!! @date    19 Feb 2025
!! @brief   Compute the power of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force(fphi, fpsi)
    use grid, only: nkx, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: kprp2
    use fields, only: phi
    use fields, only: psi
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: fphi, fpsi

    real(r8) :: phi_dot_nbl2_fphi_sum, psi_dot_nbl2_fpsi_sum
    real(r8), allocatable, dimension(:,:,:) :: phi_dot_nbl2_fphi, psi_dot_nbl2_fpsi

    integer :: i, j, k, istir

    if (proc0) call put_time_stamp(timer_force)

    allocate(phi_dot_nbl2_fphi(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(psi_dot_nbl2_fpsi(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    
    if (fix_power) then
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            phi_dot_nbl2_fphi(i, k, j) = - 0.5d0*( &
                                    phi(i, k, j)*conjg(kprp2(i, k, j)*fphi(i, k, j)) &
                                  + conjg(phi(i, k, j))*kprp2(i, k, j)*fphi(i, k, j) &
                                ) 
            psi_dot_nbl2_fpsi(i, k, j) = - 0.5d0*( &
                                    psi(i, k, j)*conjg(kprp2(i, k, j)*fpsi(i, k, j)) &
                                  + conjg(psi(i, k, j))*kprp2(i, k, j)*fpsi(i, k, j) &
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
              phi_dot_nbl2_fphi(i, k, j) = 2.0d0*phi_dot_nbl2_fphi(i, k, j)
              psi_dot_nbl2_fpsi(i, k, j) = 2.0d0*psi_dot_nbl2_fpsi(i, k, j)
            endif

          end do
        end do
      end do

      phi_dot_nbl2_fphi_sum = sum(phi_dot_nbl2_fphi); call sum_allreduce(phi_dot_nbl2_fphi_sum)
      psi_dot_nbl2_fpsi_sum = sum(psi_dot_nbl2_fpsi); call sum_allreduce(psi_dot_nbl2_fpsi_sum)

      if(abs(phi_dot_nbl2_fphi_sum) > 1.0d-6 .and. abs(psi_dot_nbl2_fpsi_sum) > 1.0d-6) then
        fphi = 0.5d0*ene_inj*(1.d0 + res_inj)*fphi/abs(phi_dot_nbl2_fphi_sum)*(-sign(1.d0, phi_dot_nbl2_fphi_sum))
        fpsi = 0.5d0*ene_inj*(1.d0 - res_inj)*fpsi/abs(psi_dot_nbl2_fpsi_sum)*(-sign(1.d0, psi_dot_nbl2_fpsi_sum))
      endif
    endif

    deallocate(phi_dot_nbl2_fphi)
    deallocate(psi_dot_nbl2_fpsi)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force


!-----------------------------------------------!
!> @author  YK
!! @date    19 Feb 2025
!! @brief   Compute the power and helicity of forcing 
!!          and normalize with it 
!-----------------------------------------------!
  subroutine normalize_force_els(fzppe, fzmpe)
    use grid, only: nkx, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: kprp2
    use fields, only: phi
    use fields, only: psi
    use mp, only: sum_allreduce
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: fzppe, fzmpe
    complex(r8) :: zppe, zmpe
    real(r8) :: zppe_dot_nbl2_fzppe_sum, zmpe_dot_nbl2_fzmpe_sum
    real(r8), allocatable, dimension(:,:,:) :: zppe_dot_nbl2_fzppe, zmpe_dot_nbl2_fzmpe

    integer :: i, j, k, istir

    if (proc0) call put_time_stamp(timer_force)

    allocate(zppe_dot_nbl2_fzppe(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(zmpe_dot_nbl2_fzmpe(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    
    if (fix_power) then
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            zppe = phi(i, k, j) + psi(i, k, j)
            zmpe = phi(i, k, j) - psi(i, k, j)
            zppe_dot_nbl2_fzppe(i, k, j) = - 0.5d0*( &
                                    zppe*conjg(kprp2(i, k, j)*fzppe(i, k, j)) &
                                  + conjg(zppe)*kprp2(i, k, j)*fzppe(i, k, j) &
                                ) 
            zmpe_dot_nbl2_fzmpe(i, k, j) = - 0.5d0*( &
                                    zmpe*conjg(kprp2(i, k, j)*fzmpe(i, k, j)) &
                                  + conjg(zmpe)*kprp2(i, k, j)*fzmpe(i, k, j) &
                                ) 
            ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
            ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
            !-----------------------------------------------------------------------------------------------------------------------------------
            ! The volume integral of a quadratic function is
            ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
            ! |f_{k_x, k_y}|^2 = \left( \sum_{k_y = -n_{k_y}/2}^{-1}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 0}
            ! \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 1}^{n_{k_y}/2}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} \right) |f_{k_x, k_y}|^2\]
            ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
            !-----------------------------------------------------------------------------------------------------------------------------------
            if (j /= 1) then
              zppe_dot_nbl2_fzppe(i, k, j) = 2.0d0*zppe_dot_nbl2_fzppe(i, k, j)
              zmpe_dot_nbl2_fzmpe(i, k, j) = 2.0d0*zmpe_dot_nbl2_fzmpe(i, k, j)
            endif

          end do
        end do
      end do

      zppe_dot_nbl2_fzppe_sum = sum(zppe_dot_nbl2_fzppe); call sum_allreduce(zppe_dot_nbl2_fzppe_sum)
      zmpe_dot_nbl2_fzmpe_sum = sum(zmpe_dot_nbl2_fzmpe); call sum_allreduce(zmpe_dot_nbl2_fzmpe_sum)

      if(abs(zppe_dot_nbl2_fzppe_sum) > 1.0d-6 .and. abs(zmpe_dot_nbl2_fzmpe_sum) > 1.0d-6) then
        fzppe = 0.5d0*ene_inj*(1.d0 + xhl_inj)*fzppe/abs(zppe_dot_nbl2_fzppe_sum)*(-sign(1.d0, zppe_dot_nbl2_fzppe_sum))
        fzmpe = 0.5d0*ene_inj*(1.d0 - xhl_inj)*fzmpe/abs(zmpe_dot_nbl2_fzmpe_sum)*(-sign(1.d0, zmpe_dot_nbl2_fzmpe_sum))
      endif
    endif

    deallocate(zppe_dot_nbl2_fzppe)
    deallocate(zmpe_dot_nbl2_fzmpe)

    if (proc0) call put_time_stamp(timer_force)

  end subroutine normalize_force_els


end module force
