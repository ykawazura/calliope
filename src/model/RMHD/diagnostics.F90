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
  public :: loop_diagnostics, loop_diagnostics_2D, loop_diagnostics_kpar, loop_diagnostics_SF2
  public :: loop_diagnostics_nltrans

  private
contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialization of diagnostics
!-----------------------------------------------!
  subroutine init_diagnostics
    use params, only: inputfile
    use diagnostics_common, only: read_parameters
    use diagnostics_common, only: init_polar_spectrum_2d
    use io, only: init_io 
    implicit none

    call read_parameters(inputfile)

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
    use grid, only: kprp2, kz2, kprp2_max, kz2_max
    use grid, only: nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use fields, only: phi, psi
    use fields, only: phi_old, psi_old
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp
    use force, only: fphi, fpsi, fphi_old, fpsi_old
    use utils, only: cabs2
    implicit none
    integer :: i, j, k

    real(r8), allocatable, dimension(:,:,:) :: upe2, bpe2
    real(r8), allocatable, dimension(:,:,:) :: upe2old, bpe2old
    real(r8), allocatable, dimension(:,:,:) :: upe2dissip_x, upe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpe2dissip_x, bpe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: p_phi, p_psi, p_xhl
    real(r8), allocatable, dimension(:,:,:) :: zppe2, zmpe2
    real(r8), allocatable, dimension(:,:,:) :: hp, hm ! two-time power spectra (10.1103/PhysRevResearch.2.023357)
    real(r8), allocatable, dimension(:,:,:) :: src

    real(r8) :: upe2_sum, bpe2_sum
    real(r8) :: upe2dot_sum, bpe2dot_sum
    real(r8) :: upe2dissip_sum, bpe2dissip_sum
    real(r8) :: zppe2_sum, zmpe2_sum
    real(r8) :: p_phi_sum, p_psi_sum, p_xhl_sum

    real(r8), dimension(:, :), allocatable :: upe2_bin , bpe2_bin      ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: zppe2_bin, zmpe2_bin
    real(r8), dimension(:, :), allocatable :: hp_bin , hm_bin     ! [nt, kprp]
    complex(r8) :: phi_mid, psi_mid, jpa_mid, fphi_mid, fpsi_mid
    complex(r8) :: zppe_mid, zmpe_mid, fzppe_mid, fzmpe_mid
    complex(r8), allocatable, dimension(:,:,:), save :: phi0, psi0
    logical :: is_first = .true.

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(upe2        , source=src)
    allocate(bpe2        , source=src)
    allocate(upe2old     , source=src)
    allocate(bpe2old     , source=src)
    allocate(upe2dissip_x, source=src)
    allocate(upe2dissip_z, source=src)
    allocate(bpe2dissip_x, source=src)
    allocate(bpe2dissip_z, source=src)
    allocate(p_phi       , source=src)
    allocate(p_psi       , source=src)
    allocate(p_xhl       , source=src)
    allocate(zppe2       , source=src)
    allocate(zmpe2       , source=src)
    allocate(hp          , source=src)
    allocate(hm          , source=src)
    deallocate(src)

    allocate (upe2_bin (1:nkpolar, nkz)); upe2_bin  = 0.d0
    allocate (bpe2_bin (1:nkpolar, nkz)); bpe2_bin  = 0.d0
    allocate (zppe2_bin(1:nkpolar, nkz)); zppe2_bin = 0.d0
    allocate (zmpe2_bin(1:nkpolar, nkz)); zmpe2_bin = 0.d0
    allocate (hp_bin   (1:nkpolar, nkz)); hp_bin    = 0.d0
    allocate (hm_bin   (1:nkpolar, nkz)); hm_bin    = 0.d0

    if (is_first) then
      allocate(phi0(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
      allocate(psi0(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
      phi0 = phi
      psi0 = psi
      is_first = .false.
    endif

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          upe2   (i, k, j) = 0.5d0*cabs2(phi(i, k, j))*kprp2(i, k, j)
          bpe2   (i, k, j) = 0.5d0*cabs2(psi(i, k, j))*kprp2(i, k, j)

          upe2old(i, k, j) = 0.5d0*cabs2(phi_old(i, k, j))*kprp2(i, k, j)
          bpe2old(i, k, j) = 0.5d0*cabs2(psi_old(i, k, j))*kprp2(i, k, j)

          phi_mid  = 0.5d0*(phi (i, k, j) + phi_old (i, k, j))
          psi_mid  = 0.5d0*(psi (i, k, j) + psi_old (i, k, j))
          jpa_mid  = -kprp2(i, k, j)*psi_mid
          fphi_mid = 0.5d0*(fphi(i, k, j) + fphi_old(i, k, j))
          fpsi_mid = 0.5d0*(fpsi(i, k, j) + fpsi_old(i, k, j))

          upe2dissip_x(i, k, j) =  nupe_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*cabs2(phi_mid)*kprp2(i, k, j)
          upe2dissip_z(i, k, j) =  nupe_z*(kz2  (k)      /kz2_max  )** nupe_z_exp*cabs2(phi_mid)*kprp2(i, k, j)
          bpe2dissip_x(i, k, j) = etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*cabs2(psi_mid)*kprp2(i, k, j)
          bpe2dissip_z(i, k, j) = etape_z*(kz2  (k)      /kz2_max  )**etape_z_exp*cabs2(psi_mid)*kprp2(i, k, j)
          p_phi       (i, k, j) = - 0.5d0*(-kprp2(i, k, j)*fphi_mid*conjg(phi_mid) + conjg(-kprp2(i, k, j)*fphi_mid)*phi_mid) 
          p_psi       (i, k, j) = - 0.5d0*(fpsi_mid*conjg(jpa_mid) + conjg(fpsi_mid)*jpa_mid) 
          zppe_mid  =  phi_mid +  psi_mid
          zmpe_mid  =  phi_mid -  psi_mid
          fzppe_mid = fphi_mid + fpsi_mid
          fzmpe_mid = fphi_mid - fpsi_mid
          zppe2       (i, k, j) = 0.5d0*cabs2(phi_mid + psi_mid)*kprp2(i, k, j)
          zmpe2       (i, k, j) = 0.5d0*cabs2(phi_mid - psi_mid)*kprp2(i, k, j)
          hp          (i, k, j) = 0.25d0*( &
                                            conjg(phi0(i, k, j) + psi0(i, k, j))*(phi_mid + psi_mid) &
                                          + (phi0(i, k, j) + psi0(i, k, j))*conjg(phi_mid + psi_mid) &
                                         )*kprp2(i, k, j)
          hm          (i, k, j) = 0.25d0*( &
                                            conjg(phi0(i, k, j) - psi0(i, k, j))*(phi_mid - psi_mid) &
                                          + (phi0(i, k, j) - psi0(i, k, j))*conjg(phi_mid - psi_mid) &
                                         )*kprp2(i, k, j)

          p_xhl       (i, k, j) = 0.25d0*( &
                                            zppe_mid*conjg(kprp2(i, k, j)*fzppe_mid) + conjg(zppe_mid)*kprp2(i, k, j)*fzppe_mid &
                                          - zmpe_mid*conjg(kprp2(i, k, j)*fzmpe_mid) - conjg(zmpe_mid)*kprp2(i, k, j)*fzmpe_mid &
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
            upe2   (i, k, j) = 2.0d0*upe2   (i, k, j)
            bpe2   (i, k, j) = 2.0d0*bpe2   (i, k, j)

            upe2old(i, k, j) = 2.0d0*upe2old(i, k, j)
            bpe2old(i, k, j) = 2.0d0*bpe2old(i, k, j)

            upe2dissip_x(i, k, j) = 2.0d0*upe2dissip_x(i, k, j)
            upe2dissip_z(i, k, j) = 2.0d0*upe2dissip_z(i, k, j)
            bpe2dissip_x(i, k, j) = 2.0d0*bpe2dissip_x(i, k, j)
            bpe2dissip_z(i, k, j) = 2.0d0*bpe2dissip_z(i, k, j)
            p_phi       (i, k, j) = 2.0d0*p_phi       (i, k, j)
            p_psi       (i, k, j) = 2.0d0*p_psi       (i, k, j)
            zppe2       (i, k, j) = 2.0d0*zppe2       (i, k, j)
            zmpe2       (i, k, j) = 2.0d0*zmpe2       (i, k, j)
            hp          (i, k, j) = 2.0d0*hp          (i, k, j)
            hm          (i, k, j) = 2.0d0*hm          (i, k, j)
            p_xhl       (i, k, j) = 2.0d0*p_xhl       (i, k, j)
          endif

        end do
      end do
    end do

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    upe2_sum = sum(upe2); call sum_reduce(upe2_sum, 0)
    bpe2_sum = sum(bpe2); call sum_reduce(bpe2_sum, 0)

    upe2dot_sum = sum((upe2 - upe2old)/dt); call sum_reduce(upe2dot_sum, 0)
    bpe2dot_sum = sum((bpe2 - bpe2old)/dt); call sum_reduce(bpe2dot_sum, 0)

    upe2dissip_sum = sum(upe2dissip_x + upe2dissip_z); call sum_reduce(upe2dissip_sum, 0)
    bpe2dissip_sum = sum(bpe2dissip_x + bpe2dissip_z); call sum_reduce(bpe2dissip_sum, 0)

    p_phi_sum      = sum(p_phi); call sum_reduce(p_phi_sum, 0)
    p_psi_sum      = sum(p_psi); call sum_reduce(p_psi_sum, 0)
    p_xhl_sum      = sum(p_xhl); call sum_reduce(p_xhl_sum, 0)

    zppe2_sum = sum(zppe2); call sum_reduce(zppe2_sum, 0)
    zmpe2_sum = sum(zmpe2); call sum_reduce(zmpe2_sum, 0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv          bin over kprp           vvvvvvvvvvvvvvvvvv!
    call get_polar_spectrum_2d(upe2 , upe2_bin )
    call get_polar_spectrum_2d(bpe2 , bpe2_bin )
    call get_polar_spectrum_2d(zppe2, zppe2_bin)
    call get_polar_spectrum_2d(zmpe2, zmpe2_bin)
    call get_polar_spectrum_2d(hp   , hp_bin   )
    call get_polar_spectrum_2d(hm   , hm_bin   )
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  upe2_sum, bpe2_sum, &
                  upe2dot_sum, bpe2dot_sum, &
                  upe2dissip_sum, bpe2dissip_sum, &
                  p_phi_sum, p_psi_sum, p_xhl_sum, &
                  zppe2_sum, zmpe2_sum, &
                  !
                  nkpolar, &
                  upe2_bin , bpe2_bin , &
                  zppe2_bin, zmpe2_bin, &
                  hp_bin   , hm_bin     &
                )

    deallocate(upe2)
    deallocate(bpe2)
    deallocate(upe2old)
    deallocate(bpe2old)
    deallocate(upe2dissip_x)
    deallocate(upe2dissip_z)
    deallocate(bpe2dissip_x)
    deallocate(bpe2dissip_z)
    deallocate(p_phi)
    deallocate(p_psi)
    deallocate(p_xhl)
    deallocate(zppe2)
    deallocate(zmpe2)
    deallocate(hp)
    deallocate(hm)

    deallocate (upe2_bin)
    deallocate (bpe2_bin)
    deallocate (zppe2_bin)
    deallocate (zmpe2_bin)
    deallocate (hp_bin)
    deallocate (hm_bin)
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    15 Apr 2020
!! @brief   Diagnostics in loop
!-----------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io, loop_io_2D
    use mp, only: proc0
    use grid, only: kx, ky, kprp2
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use fields, only: phi, psi
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi
    implicit none
    integer :: i, j, k

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
    complex(r8), allocatable, dimension(:,:,:) :: omg, jpa, ux, uy, bx, by

    real(r8)   , allocatable, dimension(:,:)   :: src1, src2, src3
    complex(r8), allocatable, dimension(:,:,:) :: src4

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(f  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f   = 0.d0
    allocate(fr (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); fr  = 0.d0

    allocate(src1(ilx_st:ilx_en, ily_st:ily_en), source=0.d0) 
    allocate(src2(ily_st:ily_en, ilz_st:ilz_en), source=0.d0) 
    allocate(src3(ilx_st:ilx_en, ilz_st:ilz_en), source=0.d0)
    allocate(phi_r_z0, source=src1)
    allocate(phi_r_x0, source=src2)
    allocate(phi_r_y0, source=src3)

    allocate(psi_r_z0, source=src1)
    allocate(psi_r_x0, source=src2)
    allocate(psi_r_y0, source=src3)

    allocate(omg_r_z0, source=src1)
    allocate(omg_r_x0, source=src2)
    allocate(omg_r_y0, source=src3)

    allocate(jpa_r_z0, source=src1)
    allocate(jpa_r_x0, source=src2)
    allocate(jpa_r_y0, source=src3)

    allocate( ux_r_z0, source=src1)
    allocate( ux_r_x0, source=src2)
    allocate( ux_r_y0, source=src3)

    allocate( uy_r_z0, source=src1)
    allocate( uy_r_x0, source=src2)
    allocate( uy_r_y0, source=src3)

    allocate( bx_r_z0, source=src1)
    allocate( bx_r_x0, source=src2)
    allocate( bx_r_y0, source=src3)

    allocate( by_r_z0, source=src1)
    allocate( by_r_x0, source=src2)
    allocate( by_r_y0, source=src3)
    deallocate(src1, src2, src3)

    allocate(src4(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    allocate(omg, source=src4)
    allocate(jpa, source=src4)
    allocate(ux , source=src4)
    allocate(uy , source=src4)
    allocate(bx , source=src4)
    allocate(by , source=src4)
    deallocate(src4)

    !vvvvvvvvvvvvvvvvvv         2D cut of fields          vvvvvvvvvvvvvvvvvv!
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          omg(i, k, j) = -kprp2(i, k, j)*phi(i, k, j)
          jpa(i, k, j) = -kprp2(i, k, j)*psi(i, k, j)
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

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io_2D( &
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
  end subroutine loop_diagnostics_2D


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Second order structure function
!-----------------------------------------------!
  subroutine loop_diagnostics_SF2
  ! under development...
  end subroutine loop_diagnostics_SF2


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Return kpar(k) & delta b/b0
!           k_\|(k) = \left(\frac{\langle|\mathbf{b}_{0,k} \cdot\nabla \delta\mathbf{b}_k|^2\rangle}
!           {\langle b_{0,k}^2\rangle\langle \delta b_k^2\rangle}\right)^{1/2}  \\ 
!           \mathbf{b}_{0,k}(\mathbf{x}) = \calF^{-1}\sum_{|\bm{k}|' \le k/2} \mathbf{b}_{\mathbf{k}'} \\  
!           \delta\mathbf{b}_{k}(\mathbf{x}) = \calF^{-1}\sum_{k/2 \le |\bm{k}|' \le 2k} \mathbf{b}_{\mathbf{k}'}
!-----------------------------------------------!
  subroutine loop_diagnostics_kpar
  ! under development...
  end subroutine loop_diagnostics_kpar


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate shell-to-shell transfer
!-----------------------------------------------!
  subroutine loop_diagnostics_nltrans
  ! under development...
  end subroutine loop_diagnostics_nltrans


end module diagnostics

