!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Diagnostics for RRMHD
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
    use grid, only: kprp2, kz2, kprp2_max, kz2_max, kx, ky
    use grid, only: nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use fields, only: phi, psi, aaa
    use fields, only: phi_old, psi_old, aaa_old
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, nonlinear,  nupe , nupe_exp , nuz , nuz_exp , &
                                      etape, etape_exp, etaz, etaz_exp, &
                                      chipe, chipe_exp, chiz, chiz_exp   
                                     
    implicit none
    integer :: i, j, k

    real(r8), allocatable, dimension(:,:,:) :: upe2, ux2, uy2, bpe2, bx2, by2, aaa2
    real(r8), allocatable, dimension(:,:,:) :: upe2old, bpe2old, aaa2old
    real(r8), allocatable, dimension(:,:,:) :: upe2dissip_x, upe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpe2dissip_x, bpe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: aaa2dissip_x, aaa2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: p
    real(r8), allocatable, dimension(:,:,:) :: zpep2, zpem2
    real(r8), allocatable, dimension(:,:,:) :: src

    real(r8) :: upe2_sum, bpe2_sum, aaa2_sum
    real(r8) :: upe2dot_sum, bpe2dot_sum, aaa2dot_sum
    real(r8) :: upe2dissip_sum, bpe2dissip_sum, aaa2dissip_sum
    real(r8) :: p_sum
    real(r8) :: zpep2_sum, zpem2_sum, zpap2_sum, zpam2_sum

    real(r8), dimension(:, :), allocatable :: upe2_bin, bpe2_bin, aaa2_bin                    ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: ux2_bin , uy2_bin , bx2_bin , by2_bin           ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: p_bin                                           ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: zpep2_bin, zpem2_bin                            ! [kprp, kz]
    complex(r8) :: phi_mid, psi_mid, jpa_mid, aaa_mid

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(upe2        , source=src)
    allocate(bpe2        , source=src)
    allocate(aaa2        , source=src)
    allocate(ux2         , source=src)
    allocate(bx2         , source=src)
    allocate(uy2         , source=src)
    allocate(by2         , source=src)
    allocate(upe2old     , source=src)
    allocate(bpe2old     , source=src)
    allocate(aaa2old     , source=src)
    allocate(upe2dissip_x, source=src)
    allocate(upe2dissip_z, source=src)
    allocate(bpe2dissip_x, source=src)
    allocate(bpe2dissip_z, source=src)
    allocate(aaa2dissip_x, source=src)
    allocate(aaa2dissip_z, source=src)
    allocate(p           , source=src)
    allocate(zpep2       , source=src)
    allocate(zpem2       , source=src)
    deallocate(src)

    allocate (upe2_bin           (1:nkpolar, nkz)); upe2_bin              = 0.d0
    allocate (bpe2_bin           (1:nkpolar, nkz)); bpe2_bin              = 0.d0
    allocate (aaa2_bin           (1:nkpolar, nkz)); aaa2_bin              = 0.d0
    allocate (ux2_bin            (1:nkpolar, nkz)); ux2_bin               = 0.d0
    allocate (uy2_bin            (1:nkpolar, nkz)); uy2_bin               = 0.d0
    allocate (bx2_bin            (1:nkpolar, nkz)); bx2_bin               = 0.d0
    allocate (by2_bin            (1:nkpolar, nkz)); by2_bin               = 0.d0
    allocate (p_bin              (1:nkpolar, nkz)); p_bin                 = 0.d0
    allocate (zpep2_bin          (1:nkpolar, nkz)); zpep2_bin             = 0.d0
    allocate (zpem2_bin          (1:nkpolar, nkz)); zpem2_bin             = 0.d0

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          upe2   (i, k, j) = 0.5d0*abs(phi(i, k, j))**2*kprp2(i, k, j)
          bpe2   (i, k, j) = 0.5d0*abs(psi(i, k, j))**2*kprp2(i, k, j)
          aaa2   (i, k, j) = 0.5d0*abs(aaa(i, k, j))**2

          ux2    (i, k, j) = 0.5d0*abs(ky(j)*phi(i, k, j))**2
          uy2    (i, k, j) = 0.5d0*abs(kx(i)*phi(i, k, j))**2
          bx2    (i, k, j) = 0.5d0*abs(ky(j)*psi(i, k, j))**2
          by2    (i, k, j) = 0.5d0*abs(kx(i)*psi(i, k, j))**2

          upe2old(i, k, j) = 0.5d0*abs(phi_old(i, k, j))**2*kprp2(i, k, j)
          bpe2old(i, k, j) = 0.5d0*abs(psi_old(i, k, j))**2*kprp2(i, k, j)
          aaa2old(i, k, j) = 0.5d0*abs(aaa_old(i, k, j))**2

          phi_mid = 0.5d0*(phi(i, k, j) + phi_old(i, k, j))
          psi_mid = 0.5d0*(psi(i, k, j) + psi_old(i, k, j))
          jpa_mid = -kprp2(i, k, j)*psi_mid
          aaa_mid = 0.5d0*(aaa(i, k, j) + aaa_old(i, k, j))

          upe2dissip_x(i, k, j) =  nupe*(kprp2(i, k, j)/kprp2_max)** nupe_exp*abs(phi_mid)**2*kprp2(i, k, j)
          upe2dissip_z(i, k, j) =  nuz *(kz2  (k)      /kz2_max  )** nuz_exp *abs(phi_mid)**2*kprp2(i, k, j)
          bpe2dissip_x(i, k, j) = etape*(kprp2(i, k, j)/kprp2_max)**etape_exp*abs(psi_mid)**2*kprp2(i, k, j)
          bpe2dissip_z(i, k, j) = etaz *(kz2  (k)      /kz2_max  )**etaz_exp *abs(psi_mid)**2*kprp2(i, k, j)
          aaa2dissip_x(i, k, j) = chipe*(kprp2(i, k, j)/kprp2_max)**chipe_exp*abs(aaa_mid)**2
          aaa2dissip_z(i, k, j) = chiz *(kz2(k)        /kz2_max  )**chiz_exp *abs(aaa_mid)**2
          p           (i, k, j) = - zi*ky(j)*0.5d0*(phi_mid*conjg(aaa_mid) - aaa_mid*conjg(phi_mid))

          zpep2       (i, k, j) = abs(phi(i, k, j) + psi(i, k, j))**2*kprp2(i, k, j)
          zpem2       (i, k, j) = abs(phi(i, k, j) - psi(i, k, j))**2*kprp2(i, k, j)

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
            aaa2   (i, k, j) = 2.0d0*aaa2   (i, k, j)

            ux2    (i, k, j) = 2.0d0*ux2    (i, k, j)
            uy2    (i, k, j) = 2.0d0*uy2    (i, k, j)
            bx2    (i, k, j) = 2.0d0*bx2    (i, k, j)
            by2    (i, k, j) = 2.0d0*by2    (i, k, j)

            upe2old(i, k, j) = 2.0d0*upe2old(i, k, j)
            bpe2old(i, k, j) = 2.0d0*bpe2old(i, k, j)
            aaa2old(i, k, j) = 2.0d0*aaa2old(i, k, j)

            upe2dissip_x(i, k, j) = 2.0d0*upe2dissip_x(i, k, j)
            upe2dissip_z(i, k, j) = 2.0d0*upe2dissip_z(i, k, j)
            bpe2dissip_x(i, k, j) = 2.0d0*bpe2dissip_x(i, k, j)
            bpe2dissip_z(i, k, j) = 2.0d0*bpe2dissip_z(i, k, j)
            aaa2dissip_x(i, k, j) = 2.0d0*aaa2dissip_x(i, k, j)
            aaa2dissip_z(i, k, j) = 2.0d0*aaa2dissip_z(i, k, j)

            p           (i, k, j) = 2.0d0*p           (i, k, j)

            zpep2    (i, k, j) = 2.0d0*zpep2(i, k, j)
            zpem2    (i, k, j) = 2.0d0*zpem2(i, k, j)
          endif

        end do
      end do
    end do

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    upe2_sum = sum(upe2); call sum_reduce(upe2_sum, 0)
    bpe2_sum = sum(bpe2); call sum_reduce(bpe2_sum, 0)
    aaa2_sum = sum(aaa2); call sum_reduce(aaa2_sum, 0)

    upe2dot_sum = sum((upe2 - upe2old)/dt); call sum_reduce(upe2dot_sum, 0)
    bpe2dot_sum = sum((bpe2 - bpe2old)/dt); call sum_reduce(bpe2dot_sum, 0)
    aaa2dot_sum = sum((aaa2 - aaa2old)/dt); call sum_reduce(aaa2dot_sum, 0)

    upe2dissip_sum = sum(upe2dissip_x + upe2dissip_z); call sum_reduce(upe2dissip_sum, 0)
    bpe2dissip_sum = sum(bpe2dissip_x + bpe2dissip_z); call sum_reduce(bpe2dissip_sum, 0)
    aaa2dissip_sum = sum(aaa2dissip_x + aaa2dissip_z); call sum_reduce(aaa2dissip_sum, 0)
    p_sum          = sum(p); call sum_reduce(p_sum   , 0)

    zpep2_sum = sum(zpep2); call sum_reduce(zpep2_sum, 0)
    zpem2_sum = sum(zpem2); call sum_reduce(zpem2_sum, 0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv          bin over kprp           vvvvvvvvvvvvvvvvvv!
    call get_polar_spectrum_2d(upe2, upe2_bin)
    call get_polar_spectrum_2d(bpe2, bpe2_bin)
    call get_polar_spectrum_2d(aaa2, aaa2_bin)
    call get_polar_spectrum_2d(ux2 , ux2_bin)
    call get_polar_spectrum_2d(uy2 , uy2_bin)
    call get_polar_spectrum_2d(bx2 , bx2_bin)
    call get_polar_spectrum_2d(by2 , by2_bin)
    call get_polar_spectrum_2d(dble(p), p_bin)
    call get_polar_spectrum_2d(zpep2, zpep2_bin)
    call get_polar_spectrum_2d(zpem2, zpem2_bin)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    ! call get_nonlinear_transfer(ntrans_aw_l_bin, ntrans_aw_g_bin, ntrans_compr_l_bin, ntrans_compr_g_bin)

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  upe2_sum, bpe2_sum, aaa2_sum, &
                  upe2dot_sum, bpe2dot_sum, aaa2dot_sum, &
                  upe2dissip_sum, bpe2dissip_sum, aaa2dissip_sum, &
                  p_sum, &
                  zpep2_sum, zpem2_sum, &
                  !
                  nkpolar, &
                  upe2_bin, bpe2_bin, aaa2_bin, &
                  ux2_bin , uy2_bin , bx2_bin , by2_bin , &
                  zpep2_bin, zpem2_bin, &
                  p_bin &
                )

    deallocate(upe2)
    deallocate(bpe2)
    deallocate(aaa2)
    deallocate(upe2old)
    deallocate(bpe2old)
    deallocate(aaa2old)
    deallocate(upe2dissip_x)
    deallocate(upe2dissip_z)
    deallocate(bpe2dissip_x)
    deallocate(bpe2dissip_z)
    deallocate(aaa2dissip_x)
    deallocate(aaa2dissip_z)
    deallocate(p)

    deallocate (upe2_bin)
    deallocate (bpe2_bin)
    deallocate (aaa2_bin)
    deallocate (p_bin)
    deallocate (zpep2_bin)
    deallocate (zpem2_bin)
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
    use fields, only: phi, psi, aaa
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi
    implicit none
    integer :: i, j, k

    complex(r8), allocatable, dimension(:,:,:) :: f
    real(r8)   , allocatable, dimension(:,:,:) :: fr
    real(r8)   , allocatable, dimension(:,:)   :: phi_r_z0, phi_r_x0, phi_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: psi_r_z0, psi_r_x0, psi_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: aaa_r_z0, aaa_r_x0, aaa_r_y0
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

    allocate(aaa_r_z0, source=src1)
    allocate(aaa_r_x0, source=src2)
    allocate(aaa_r_y0, source=src3)

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
    f = aaa; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr, aaa_r_z0, aaa_r_x0, aaa_r_y0)
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
                  aaa_r_z0, aaa_r_x0, aaa_r_y0, &
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

    deallocate (aaa_r_z0)
    deallocate (aaa_r_x0)
    deallocate (aaa_r_y0)

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

