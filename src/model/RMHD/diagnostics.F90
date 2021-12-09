!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Diagnostics for RMHD
!-----------------------------------------------!
module diagnostics
  use diagnostics_common
  use p3dfft
  implicit none

  public :: init_diagnostics, finish_diagnostics
  public :: loop_diagnostics, loop_diagnostics_fields_secion, loop_diagnostics_SF2

  private
contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialization of diagnostics
!-----------------------------------------------!
  subroutine init_diagnostics
    use diagnostics_common, only: init_polar_spectrum_2d
    use io, only: init_io 
    implicit none

    call init_polar_spectrum_2d
    call init_io(nkpolar, kpbin)
  end subroutine init_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Diagnostics in loop
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io
    use mp, only: proc0
    use grid, only: kperp2, kz2, kperp2_max, kz2_max
    use grid, only: nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use fields, only: phi, psi
    use fields, only: phi_old1, psi_old1
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics
    use params, only: zi, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp
    use force, only: fomg, fpsi, fomg_old1, fpsi_old1
    implicit none
    integer :: i, j, k

    real(r8)   , allocatable, dimension(:,:,:) :: upe2, bpe2
    real(r8)   , allocatable, dimension(:,:,:) :: upe2old1, bpe2old1
    real(r8)   , allocatable, dimension(:,:,:) :: upe2dissip_x, upe2dissip_z
    real(r8)   , allocatable, dimension(:,:,:) :: bpe2dissip_x, bpe2dissip_z
    real(r8)   , allocatable, dimension(:,:,:) :: p_omg, p_psi

    real(r8) :: upe2_sum, bpe2_sum
    real(r8) :: upe2dot_sum, bpe2dot_sum
    real(r8) :: upe2dissip_sum, bpe2dissip_sum
    real(r8) :: p_omg_sum, p_psi_sum

    real(r8), dimension(:, :), allocatable :: upe2_bin, bpe2_bin      ! [kperp, kz]
    complex(r8) :: phi_mid, psi_mid, jpa_mid, fomg_mid, fpsi_mid

    if (proc0) call put_time_stamp(timer_diagnostics)

    allocate(upe2        (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); upe2     = 0.d0
    allocate(bpe2        (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bpe2     = 0.d0
    allocate(upe2old1    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); upe2old1 = 0.d0
    allocate(bpe2old1    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bpe2old1 = 0.d0
    allocate(upe2dissip_x(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); upe2dissip_x = 0.d0
    allocate(upe2dissip_z(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); upe2dissip_z = 0.d0
    allocate(bpe2dissip_x(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bpe2dissip_x = 0.d0
    allocate(bpe2dissip_z(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bpe2dissip_z = 0.d0
    allocate(p_omg       (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); p_omg        = 0.d0
    allocate(p_psi       (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); p_psi        = 0.d0

    allocate (upe2_bin(1:nkpolar, nkz)); upe2_bin = 0.d0
    allocate (bpe2_bin(1:nkpolar, nkz)); bpe2_bin = 0.d0

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          upe2    (i, k, j) = 0.5d0*abs(phi(i, k, j))**2*kperp2(i, k, j)
          bpe2    (i, k, j) = 0.5d0*abs(psi(i, k, j))**2*kperp2(i, k, j)

          upe2old1(i, k, j) = 0.5d0*abs(phi_old1(i, k, j))**2*kperp2(i, k, j)
          bpe2old1(i, k, j) = 0.5d0*abs(psi_old1(i, k, j))**2*kperp2(i, k, j)

          phi_mid  = 0.5d0*(phi (i, k, j) + phi_old1 (i, k, j))
          psi_mid  = 0.5d0*(psi (i, k, j) + psi_old1 (i, k, j))
          jpa_mid  = -kperp2(i, k, j)*psi_mid
          fomg_mid = 0.5d0*(fomg(i, k, j) + fomg_old1(i, k, j))
          fpsi_mid = 0.5d0*(fpsi(i, k, j) + fpsi_old1(i, k, j))

          upe2dissip_x(i, k, j) =  nupe_x*(kperp2(i, k, j)/kperp2_max)** nupe_x_exp*abs(phi_mid)**2*kperp2(i, k, j)
          upe2dissip_z(i, k, j) =  nupe_z*(kz2   (k)      /kz2_max   )** nupe_z_exp*abs(phi_mid)**2*kperp2(i, k, j)
          bpe2dissip_x(i, k, j) = etape_x*(kperp2(i, k, j)/kperp2_max)**etape_x_exp*abs(psi_mid)**2*kperp2(i, k, j)
          bpe2dissip_z(i, k, j) = etape_z*(kz2   (k)      /kz2_max   )**etape_z_exp*abs(psi_mid)**2*kperp2(i, k, j)
          p_omg      (i, k, j) = - 0.5d0*(fomg_mid*conjg(phi_mid) + conjg(fomg_mid)*phi_mid) 
          p_psi      (i, k, j) = - 0.5d0*(fpsi_mid*conjg(jpa_mid) + conjg(fpsi_mid)*jpa_mid) 

          ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
          ! |f_{k_x, k_y}|^2 = \left( \sum_{k_y = -n_{k_y}/2}^{-1}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 0}
          ! \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 1}^{n_{k_y}/2}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} \right) |f_{k_x, k_y}|^2\]
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (j /= 1) then
            upe2    (i, k, j) = 2.0d0*upe2    (i, k, j)
            bpe2    (i, k, j) = 2.0d0*bpe2    (i, k, j)

            upe2old1(i, k, j) = 2.0d0*upe2old1(i, k, j)
            bpe2old1(i, k, j) = 2.0d0*bpe2old1(i, k, j)

            upe2dissip_x(i, k, j) = 2.0d0*upe2dissip_x(i, k, j)
            upe2dissip_z(i, k, j) = 2.0d0*upe2dissip_z(i, k, j)
            bpe2dissip_x(i, k, j) = 2.0d0*bpe2dissip_x(i, k, j)
            bpe2dissip_z(i, k, j) = 2.0d0*bpe2dissip_z(i, k, j)
            p_omg       (i, k, j) = 2.0d0*p_omg       (i, k, j)
            p_psi       (i, k, j) = 2.0d0*p_psi       (i, k, j)
          endif

        end do
      end do
    end do

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    upe2_sum = sum(upe2); call sum_reduce(upe2_sum, 0)
    bpe2_sum = sum(bpe2); call sum_reduce(bpe2_sum, 0)

    upe2dot_sum = sum((upe2 - upe2old1)/dt); call sum_reduce(upe2dot_sum, 0)
    bpe2dot_sum = sum((bpe2 - bpe2old1)/dt); call sum_reduce(bpe2dot_sum, 0)

    upe2dissip_sum = sum(upe2dissip_x + upe2dissip_z); call sum_reduce(upe2dissip_sum, 0)
    bpe2dissip_sum = sum(bpe2dissip_x + bpe2dissip_z); call sum_reduce(bpe2dissip_sum, 0)

    p_omg_sum      = sum(p_omg); call sum_reduce(p_omg_sum, 0)
    p_psi_sum      = sum(p_psi); call sum_reduce(p_psi_sum, 0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv          bin over kperp           vvvvvvvvvvvvvvvvvv!
    call get_polar_spectrum_2d(upe2, upe2_bin)
    call get_polar_spectrum_2d(bpe2, bpe2_bin)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    call loop_io( &
                  upe2_sum, bpe2_sum, &
                  upe2dot_sum, bpe2dot_sum, &
                  upe2dissip_sum, bpe2dissip_sum, &
                  p_omg_sum, p_psi_sum, &
                  !
                  nkpolar, &
                  upe2_bin, bpe2_bin &
                )

    deallocate(upe2)
    deallocate(bpe2)
    deallocate(upe2old1)
    deallocate(bpe2old1)
    deallocate(upe2dissip_x)
    deallocate(upe2dissip_z)
    deallocate(bpe2dissip_x)
    deallocate(bpe2dissip_z)
    deallocate(p_omg)
    deallocate(p_psi)

    deallocate (upe2_bin)
    deallocate (bpe2_bin)

    if (proc0) call put_time_stamp(timer_diagnostics)
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    15 Apr 2020
!! @brief   Diagnostics in loop
!-----------------------------------------------!
  subroutine loop_diagnostics_fields_secion
    use io, only: loop_io, loop_io_fields_section
    use mp, only: proc0
    use grid, only: nlx, nly, nlz
    use grid, only: kx, ky, kperp2
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use fields, only: phi, psi
    use time_stamp, only: put_time_stamp, timer_diagnostics
    use params, only: zi
    implicit none
    integer :: i, j, k

    complex(r8), allocatable, dimension(:,:,:) :: omg, jpa, ux, uy, bx, by
    complex(r8), allocatable, dimension(:,:,:) :: f
    real(r8)   , allocatable, dimension(:,:,:) :: fr
    real(r8)   , allocatable, dimension(:,:)   :: phi_r_z0, phi_r_x0, phi_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: psi_r_z0, psi_r_x0, psi_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: omg_r_z0, omg_r_x0, omg_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: jpa_r_z0, jpa_r_x0, jpa_r_y0
    real(r8)   , allocatable, dimension(:,:)   ::  ux_r_z0,  ux_r_x0,  ux_r_y0
    real(r8)   , allocatable, dimension(:,:)   ::  uy_r_z0,  uy_r_x0,  uy_r_y0
    real(r8)   , allocatable, dimension(:,:)   ::  bx_r_z0,  bx_r_x0,  bx_r_y0
    real(r8)   , allocatable, dimension(:,:)   ::  by_r_z0,  by_r_x0,  by_r_y0

    if (proc0) call put_time_stamp(timer_diagnostics)

    allocate(omg(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); omg = 0.d0
    allocate(jpa(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); jpa = 0.d0
    allocate(ux (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); ux  = 0.d0
    allocate(uy (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); uy  = 0.d0
    allocate(bx (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bx  = 0.d0
    allocate(by (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); by  = 0.d0
    allocate(f  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f   = 0.d0
    allocate(fr (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); fr  = 0.d0

    allocate(phi_r_z0(nlx, nly)); phi_r_z0 = 0.d0
    allocate(phi_r_x0(nly, nlz)); phi_r_x0 = 0.d0
    allocate(phi_r_y0(nlx, nlz)); phi_r_y0 = 0.d0

    allocate(psi_r_z0(nlx, nly)); psi_r_z0 = 0.d0
    allocate(psi_r_x0(nly, nlz)); psi_r_x0 = 0.d0
    allocate(psi_r_y0(nlx, nlz)); psi_r_y0 = 0.d0

    allocate(omg_r_z0(nlx, nly)); omg_r_z0 = 0.d0
    allocate(omg_r_x0(nly, nlz)); omg_r_x0 = 0.d0
    allocate(omg_r_y0(nlx, nlz)); omg_r_y0 = 0.d0

    allocate(jpa_r_z0(nlx, nly)); jpa_r_z0 = 0.d0
    allocate(jpa_r_x0(nly, nlz)); jpa_r_x0 = 0.d0
    allocate(jpa_r_y0(nlx, nlz)); jpa_r_y0 = 0.d0

    allocate( ux_r_z0(nlx, nly));  ux_r_z0 = 0.d0
    allocate( ux_r_x0(nly, nlz));  ux_r_x0 = 0.d0
    allocate( ux_r_y0(nlx, nlz));  ux_r_y0 = 0.d0

    allocate( uy_r_z0(nlx, nly));  uy_r_z0 = 0.d0
    allocate( uy_r_x0(nly, nlz));  uy_r_x0 = 0.d0
    allocate( uy_r_y0(nlx, nlz));  uy_r_y0 = 0.d0

    allocate( bx_r_z0(nlx, nly));  bx_r_z0 = 0.d0
    allocate( bx_r_x0(nly, nlz));  bx_r_x0 = 0.d0
    allocate( bx_r_y0(nlx, nlz));  bx_r_y0 = 0.d0

    allocate( by_r_z0(nlx, nly));  by_r_z0 = 0.d0
    allocate( by_r_x0(nly, nlz));  by_r_x0 = 0.d0
    allocate( by_r_y0(nlx, nlz));  by_r_y0 = 0.d0


    !vvvvvvvvvvvvvvvvvv         2D cut of fields          vvvvvvvvvvvvvvvvvv!
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          omg(i, k, j) = -kperp2(i, k, j)*phi(i, k, j)
          jpa(i, k, j) = -kperp2(i, k, j)*psi(i, k, j)
          ux (i, k, j) = -zi*ky(j)*phi(i, k, j)
          uy (i, k, j) =  zi*kx(i)*phi(i, k, j)
          bx (i, k, j) = -zi*ky(j)*psi(i, k, j)
          by (i, k, j) =  zi*kx(i)*psi(i, k, j)
        enddo
      enddo
    enddo
    !$omp end parallel do
    f = phi; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr, phi_r_z0, phi_r_x0, phi_r_y0)
    f = psi; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr, psi_r_z0, psi_r_x0, psi_r_y0)
    f = omg; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr, omg_r_z0, omg_r_x0, omg_r_y0)
    f = jpa; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr, jpa_r_z0, jpa_r_x0, jpa_r_y0)
    f = ux ; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr,  ux_r_z0,  ux_r_x0,  ux_r_y0)
    f = uy ; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr,  uy_r_z0,  uy_r_x0,  uy_r_y0)
    f = bx ; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr,  bx_r_z0,  bx_r_x0,  bx_r_y0)
    f = by ; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr,  by_r_z0,  by_r_x0,  by_r_y0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call loop_io_fields_section( &
                  phi_r_z0, phi_r_x0, phi_r_y0, &
                  psi_r_z0, psi_r_x0, psi_r_y0, &
                  omg_r_z0, omg_r_x0, omg_r_y0, &
                  jpa_r_z0, jpa_r_x0, jpa_r_y0, &
                   ux_r_z0,  ux_r_x0,  ux_r_y0, &
                   uy_r_z0,  uy_r_x0,  uy_r_y0, &
                   bx_r_z0,  bx_r_x0,  bx_r_y0, &
                   by_r_z0,  by_r_x0,  by_r_y0  &
                )

    deallocate (omg)
    deallocate (jpa)
    deallocate (ux)
    deallocate (uy)
    deallocate (bx)
    deallocate (by)
    deallocate (f)
    deallocate (fr)

    deallocate (phi_r_z0)
    deallocate (phi_r_x0)
    deallocate (phi_r_y0)

    deallocate (psi_r_z0)
    deallocate (psi_r_x0)
    deallocate (psi_r_y0)

    deallocate (omg_r_z0)
    deallocate (omg_r_x0)
    deallocate (omg_r_y0)

    deallocate (jpa_r_z0)
    deallocate (jpa_r_x0)
    deallocate (jpa_r_y0)

    deallocate ( ux_r_z0)
    deallocate ( ux_r_x0)
    deallocate ( ux_r_y0)

    deallocate ( uy_r_z0)
    deallocate ( uy_r_x0)
    deallocate ( uy_r_y0)

    deallocate ( bx_r_z0)
    deallocate ( bx_r_x0)
    deallocate ( bx_r_y0)

    deallocate ( by_r_z0)
    deallocate ( by_r_x0)
    deallocate ( by_r_y0)

    if (proc0) call put_time_stamp(timer_diagnostics)
  end subroutine loop_diagnostics_fields_secion


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Second order structure function
!-----------------------------------------------!
  subroutine loop_diagnostics_SF2
  end subroutine loop_diagnostics_SF2

end module diagnostics

