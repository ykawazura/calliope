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
    use fields, only: phi, psi, upa, bpa
    use fields, only: phi_old, psi_old, upa_old, bpa_old
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, nonlinear, q, &
                      va2cs2_plus_1, nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                                     nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                                     etape_x, etape_x_exp, etape_z, etape_z_exp, &
                                     etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
    implicit none
    integer :: i, j, k

    real(r8), allocatable, dimension(:,:,:) :: upe2, ux2, uy2, bpe2, bx2, by2, upa2, bpa2
    real(r8), allocatable, dimension(:,:,:) :: upe2old, bpe2old, upa2old, bpa2old
    real(r8), allocatable, dimension(:,:,:) :: upe2dissip_x, upe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpe2dissip_x, bpe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: upa2dissip_x, upa2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpa2dissip_x, bpa2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: p_aw, p_compr
    real(r8), allocatable, dimension(:,:,:) :: zpep2, zpem2, zpap2, zpam2
    real(r8), allocatable, dimension(:,:,:) :: src

    real(r8) :: upe2_sum, bpe2_sum, upa2_sum, bpa2_sum
    real(r8) :: upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum
    real(r8) :: upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum
    real(r8) :: p_aw_sum, p_compr_sum
    real(r8) :: zpep2_sum, zpem2_sum, zpap2_sum, zpam2_sum

    real(r8), dimension(:, :), allocatable :: upe2_bin, bpe2_bin, upa2_bin, bpa2_bin          ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: ux2_bin , uy2_bin , bx2_bin , by2_bin           ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: p_aw_bin, p_compr_bin                           ! [kprp, kz]
    real(r8), dimension(:, :, :), allocatable :: ntrans_aw_l_bin   , ntrans_aw_g_bin          ! [4, kprp, kz]
                                                                                              !   1 : -upe.(upe.grad upe) 
                                                                                              !   2 : +upe.(bpe.grad bpe) 
                                                                                              !   3 : -bpe.(upe.grad bpe) 
                                                                                              !   4 : +bpe.(bpe.grad upe) 
    real(r8), dimension(:, :, :), allocatable :: ntrans_compr_l_bin, ntrans_compr_g_bin   ! [4, kprp, kz]
                                                                                              !   1 : -upa (upe.grad upa) 
                                                                                              !   2 : +upa (bpe.grad bpa) 
                                                                                              !   3 : -bpa (upe.grad bpa) 
                                                                                              !   4 : +bpa (bpe.grad upa) 
    real(r8), dimension(:, :), allocatable :: dissip_aw_bin, dissip_compr_bin             ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: zpep2_bin, zpem2_bin, zpap2_bin, zpam2_bin  ! [kprp, kz]
    complex(r8) :: phi_mid, psi_mid, jpa_mid, upa_mid, bpa_mid

    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(upe2        , source=src)
    allocate(bpe2        , source=src)
    allocate(upa2        , source=src)
    allocate(bpa2        , source=src)
    allocate(ux2         , source=src)
    allocate(bx2         , source=src)
    allocate(uy2         , source=src)
    allocate(by2         , source=src)
    allocate(upe2old     , source=src)
    allocate(bpe2old     , source=src)
    allocate(upa2old     , source=src)
    allocate(bpa2old     , source=src)
    allocate(upe2dissip_x, source=src)
    allocate(upe2dissip_z, source=src)
    allocate(bpe2dissip_x, source=src)
    allocate(bpe2dissip_z, source=src)
    allocate(upa2dissip_x, source=src)
    allocate(upa2dissip_z, source=src)
    allocate(bpa2dissip_x, source=src)
    allocate(bpa2dissip_z, source=src)
    allocate(p_aw        , source=src)
    allocate(p_compr     , source=src)
    allocate(zpep2       , source=src)
    allocate(zpem2       , source=src)
    allocate(zpap2       , source=src)
    allocate(zpam2       , source=src)
    deallocate(src)

    allocate (upe2_bin           (1:nkpolar, nkz)); upe2_bin              = 0.d0
    allocate (bpe2_bin           (1:nkpolar, nkz)); bpe2_bin              = 0.d0
    allocate (upa2_bin           (1:nkpolar, nkz)); upa2_bin              = 0.d0
    allocate (bpa2_bin           (1:nkpolar, nkz)); bpa2_bin              = 0.d0
    allocate (ux2_bin            (1:nkpolar, nkz)); ux2_bin               = 0.d0
    allocate (uy2_bin            (1:nkpolar, nkz)); uy2_bin               = 0.d0
    allocate (bx2_bin            (1:nkpolar, nkz)); bx2_bin               = 0.d0
    allocate (by2_bin            (1:nkpolar, nkz)); by2_bin               = 0.d0
    allocate (p_aw_bin           (1:nkpolar, nkz)); p_aw_bin              = 0.d0
    allocate (p_compr_bin        (1:nkpolar, nkz)); p_compr_bin           = 0.d0
    allocate (ntrans_aw_l_bin    (4, 1:nkpolar, nkz)); ntrans_aw_l_bin    = 0.d0
    allocate (ntrans_aw_g_bin    (4, 1:nkpolar, nkz)); ntrans_aw_g_bin    = 0.d0
    allocate (ntrans_compr_l_bin (4, 1:nkpolar, nkz)); ntrans_compr_l_bin = 0.d0
    allocate (ntrans_compr_g_bin (4, 1:nkpolar, nkz)); ntrans_compr_g_bin = 0.d0
    allocate (dissip_aw_bin      (1:nkpolar, nkz)); dissip_aw_bin         = 0.d0
    allocate (dissip_compr_bin   (1:nkpolar, nkz)); dissip_compr_bin      = 0.d0
    allocate (zpep2_bin          (1:nkpolar, nkz)); zpep2_bin             = 0.d0
    allocate (zpem2_bin          (1:nkpolar, nkz)); zpem2_bin             = 0.d0
    allocate (zpap2_bin          (1:nkpolar, nkz)); zpap2_bin             = 0.d0
    allocate (zpam2_bin          (1:nkpolar, nkz)); zpam2_bin             = 0.d0

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          upe2   (i, k, j) = 0.5d0*abs(phi(i, k, j))**2*kprp2(i, k, j)
          bpe2   (i, k, j) = 0.5d0*abs(psi(i, k, j))**2*kprp2(i, k, j)
          upa2   (i, k, j) = 0.5d0*abs(upa    (i, k, j))**2
          bpa2   (i, k, j) = 0.5d0*abs(bpa    (i, k, j))**2*va2cs2_plus_1

          ux2    (i, k, j) = 0.5d0*abs(ky(j)*phi(i, k, j))**2
          uy2    (i, k, j) = 0.5d0*abs(kx(i)*phi(i, k, j))**2
          bx2    (i, k, j) = 0.5d0*abs(ky(j)*psi(i, k, j))**2
          by2    (i, k, j) = 0.5d0*abs(kx(i)*psi(i, k, j))**2

          upe2old(i, k, j) = 0.5d0*abs(phi_old(i, k, j))**2*kprp2(i, k, j)
          bpe2old(i, k, j) = 0.5d0*abs(psi_old(i, k, j))**2*kprp2(i, k, j)
          upa2old(i, k, j) = 0.5d0*abs(upa_old(i, k, j))**2
          bpa2old(i, k, j) = 0.5d0*abs(bpa_old(i, k, j))**2*va2cs2_plus_1

          phi_mid = 0.5d0*(phi(i, k, j) + phi_old(i, k, j))
          psi_mid = 0.5d0*(psi(i, k, j) + psi_old(i, k, j))
          jpa_mid = -kprp2(i, k, j)*psi_mid
          upa_mid = 0.5d0*(upa(i, k, j) + upa_old(i, k, j))
          bpa_mid = 0.5d0*(bpa(i, k, j) + bpa_old(i, k, j))

          upe2dissip_x(i, k, j) =  nupe_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*abs(phi_mid)**2*kprp2(i, k, j)
          upe2dissip_z(i, k, j) =  nupe_z*(kz2  (k)      /kz2_max  )** nupe_z_exp*abs(phi_mid)**2*kprp2(i, k, j)
          bpe2dissip_x(i, k, j) = etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*abs(psi_mid)**2*kprp2(i, k, j)
          bpe2dissip_z(i, k, j) = etape_z*(kz2  (k)      /kz2_max  )**etape_z_exp*abs(psi_mid)**2*kprp2(i, k, j)
          upa2dissip_x(i, k, j) =  nupa_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*abs(upa_mid)**2
          upa2dissip_z(i, k, j) =  nupa_z*(kz2(k)        /kz2_max  )** nupe_z_exp*abs(upa_mid)**2
          bpa2dissip_x(i, k, j) = etapa_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*abs(bpa_mid)**2
          bpa2dissip_z(i, k, j) = etapa_z*(kz2(k)        /kz2_max  )**etape_z_exp*abs(bpa_mid)**2
          p_aw        (i, k, j) = - 2.d0*zi*ky(j)*0.5d0*(phi_mid*conjg(upa_mid) - upa_mid*conjg(phi_mid))
          p_compr     (i, k, j) = zi*ky(j)*( &
                                                       q *0.5d0*(psi_mid*conjg(bpa_mid) - bpa_mid*conjg(psi_mid)) &
                                             + (2.d0 - q)*0.5d0*(phi_mid*conjg(upa_mid) - upa_mid*conjg(phi_mid)) &
                                           )

          zpep2      (i, k, j) = abs(phi(i, k, j) + psi(i, k, j))**2*kprp2(i, k, j)
          zpem2      (i, k, j) = abs(phi(i, k, j) - psi(i, k, j))**2*kprp2(i, k, j)
          zpap2      (i, k, j) = abs(upa(i, k, j) + bpa(i, k, j)*sqrt(va2cs2_plus_1))**2
          zpam2      (i, k, j) = abs(upa(i, k, j) - bpa(i, k, j)*sqrt(va2cs2_plus_1))**2

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
            upa2   (i, k, j) = 2.0d0*upa2   (i, k, j)
            bpa2   (i, k, j) = 2.0d0*bpa2   (i, k, j)

            ux2    (i, k, j) = 2.0d0*ux2    (i, k, j)
            uy2    (i, k, j) = 2.0d0*uy2    (i, k, j)
            bx2    (i, k, j) = 2.0d0*bx2    (i, k, j)
            by2    (i, k, j) = 2.0d0*by2    (i, k, j)

            upe2old(i, k, j) = 2.0d0*upe2old(i, k, j)
            bpe2old(i, k, j) = 2.0d0*bpe2old(i, k, j)
            upa2old(i, k, j) = 2.0d0*upa2old(i, k, j)
            bpa2old(i, k, j) = 2.0d0*bpa2old(i, k, j)

            upe2dissip_x(i, k, j) = 2.0d0*upe2dissip_x(i, k, j)
            upe2dissip_z(i, k, j) = 2.0d0*upe2dissip_z(i, k, j)
            bpe2dissip_x(i, k, j) = 2.0d0*bpe2dissip_x(i, k, j)
            bpe2dissip_z(i, k, j) = 2.0d0*bpe2dissip_z(i, k, j)
            upa2dissip_x(i, k, j) = 2.0d0*upa2dissip_x(i, k, j)
            upa2dissip_z(i, k, j) = 2.0d0*upa2dissip_z(i, k, j)
            bpa2dissip_x(i, k, j) = 2.0d0*bpa2dissip_x(i, k, j)
            bpa2dissip_z(i, k, j) = 2.0d0*bpa2dissip_z(i, k, j)

            p_aw        (i, k, j) = 2.0d0*p_aw        (i, k, j)
            p_compr     (i, k, j) = 2.0d0*p_compr     (i, k, j)

            zpep2    (i, k, j) = 2.0d0*zpep2(i, k, j)
            zpem2    (i, k, j) = 2.0d0*zpem2(i, k, j)
            zpap2    (i, k, j) = 2.0d0*zpap2(i, k, j)
            zpam2    (i, k, j) = 2.0d0*zpam2(i, k, j)
          endif

        end do
      end do
    end do

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    upe2_sum = sum(upe2); call sum_reduce(upe2_sum, 0)
    bpe2_sum = sum(bpe2); call sum_reduce(bpe2_sum, 0)
    upa2_sum = sum(upa2); call sum_reduce(upa2_sum, 0)
    bpa2_sum = sum(bpa2); call sum_reduce(bpa2_sum, 0)

    upe2dot_sum = sum((upe2 - upe2old)/dt); call sum_reduce(upe2dot_sum, 0)
    bpe2dot_sum = sum((bpe2 - bpe2old)/dt); call sum_reduce(bpe2dot_sum, 0)
    upa2dot_sum = sum((upa2 - upa2old)/dt); call sum_reduce(upa2dot_sum, 0)
    bpa2dot_sum = sum((bpa2 - bpa2old)/dt); call sum_reduce(bpa2dot_sum, 0)

    upe2dissip_sum = sum(upe2dissip_x + upe2dissip_z); call sum_reduce(upe2dissip_sum, 0)
    bpe2dissip_sum = sum(bpe2dissip_x + bpe2dissip_z); call sum_reduce(bpe2dissip_sum, 0)
    upa2dissip_sum = sum(upa2dissip_x + upa2dissip_z); call sum_reduce(upa2dissip_sum, 0)
    bpa2dissip_sum = sum(bpa2dissip_x + bpa2dissip_z); call sum_reduce(bpa2dissip_sum, 0)
    p_aw_sum       = sum(p_aw   ); call sum_reduce(p_aw_sum   , 0)
    p_compr_sum    = sum(p_compr); call sum_reduce(p_compr_sum, 0)

    zpep2_sum = sum(zpep2); call sum_reduce(zpep2_sum, 0)
    zpem2_sum = sum(zpem2); call sum_reduce(zpem2_sum, 0)
    zpap2_sum = sum(zpap2); call sum_reduce(zpap2_sum, 0)
    zpam2_sum = sum(zpam2); call sum_reduce(zpam2_sum, 0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv          bin over kprp           vvvvvvvvvvvvvvvvvv!
    call get_polar_spectrum_2d(upe2, upe2_bin)
    call get_polar_spectrum_2d(bpe2, bpe2_bin)
    call get_polar_spectrum_2d(upa2, upa2_bin)
    call get_polar_spectrum_2d(bpa2, bpa2_bin)
    call get_polar_spectrum_2d(ux2 , ux2_bin)
    call get_polar_spectrum_2d(uy2 , uy2_bin)
    call get_polar_spectrum_2d(bx2 , bx2_bin)
    call get_polar_spectrum_2d(by2 , by2_bin)
    call get_polar_spectrum_2d(dble(p_aw   ), p_aw_bin   )
    call get_polar_spectrum_2d(dble(p_compr), p_compr_bin)
    call get_polar_spectrum_2d(upe2dissip_x + upe2dissip_z + bpe2dissip_x + bpe2dissip_z, dissip_aw_bin)
    call get_polar_spectrum_2d(upa2dissip_x + upa2dissip_z + bpa2dissip_x + bpa2dissip_z, dissip_compr_bin)
    call get_polar_spectrum_2d(zpep2, zpep2_bin)
    call get_polar_spectrum_2d(zpem2, zpem2_bin)
    call get_polar_spectrum_2d(zpap2, zpap2_bin)
    call get_polar_spectrum_2d(zpam2, zpam2_bin)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    ! call get_nonlinear_transfer(ntrans_aw_l_bin, ntrans_aw_g_bin, ntrans_compr_l_bin, ntrans_compr_g_bin)

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  upe2_sum, bpe2_sum, upa2_sum, bpa2_sum, &
                  upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum, &
                  upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum, &
                  p_aw_sum, p_compr_sum, &
                  zpep2_sum, zpem2_sum, zpap2_sum, zpam2_sum, &
                  !
                  nkpolar, &
                  upe2_bin, bpe2_bin, upa2_bin, bpa2_bin, &
                  ux2_bin , uy2_bin , bx2_bin , by2_bin , &
                  zpep2_bin, zpem2_bin, zpap2_bin, zpam2_bin, &
                  p_aw_bin, p_compr_bin, &
                  dissip_aw_bin, dissip_compr_bin, &
                  ntrans_aw_l_bin, ntrans_compr_l_bin, &
                  ntrans_aw_g_bin, ntrans_compr_g_bin  &
                )

    deallocate(upe2)
    deallocate(bpe2)
    deallocate(upa2)
    deallocate(bpa2)
    deallocate(upe2old)
    deallocate(bpe2old)
    deallocate(upa2old)
    deallocate(bpa2old)
    deallocate(upe2dissip_x)
    deallocate(upe2dissip_z)
    deallocate(bpe2dissip_x)
    deallocate(bpe2dissip_z)
    deallocate(upa2dissip_x)
    deallocate(upa2dissip_z)
    deallocate(bpa2dissip_x)
    deallocate(bpa2dissip_z)
    deallocate(p_aw)
    deallocate(p_compr)

    deallocate (upe2_bin)
    deallocate (bpe2_bin)
    deallocate (upa2_bin)
    deallocate (bpa2_bin)
    deallocate (p_aw_bin)
    deallocate (p_compr_bin)
    deallocate (ntrans_aw_l_bin)
    deallocate (ntrans_aw_g_bin)
    deallocate (ntrans_compr_l_bin)
    deallocate (ntrans_compr_g_bin)
    deallocate (dissip_aw_bin)
    deallocate (dissip_compr_bin)
    deallocate (zpep2_bin)
    deallocate (zpem2_bin)
    deallocate (zpap2_bin)
    deallocate (zpam2_bin)
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
    use fields, only: phi, psi, upa, bpa
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi
    implicit none
    integer :: i, j, k

    complex(r8), allocatable, dimension(:,:,:) :: f
    real(r8)   , allocatable, dimension(:,:,:) :: fr
    real(r8)   , allocatable, dimension(:,:)   :: phi_r_z0, phi_r_x0, phi_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: psi_r_z0, psi_r_x0, psi_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: upa_r_z0, upa_r_x0, upa_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: bpa_r_z0, bpa_r_x0, bpa_r_y0
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

    allocate(upa_r_z0, source=src1)
    allocate(upa_r_x0, source=src2)
    allocate(upa_r_y0, source=src3)

    allocate(bpa_r_z0, source=src1)
    allocate(bpa_r_x0, source=src2)
    allocate(bpa_r_y0, source=src3)

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
    f = upa; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr, upa_r_z0, upa_r_x0, upa_r_y0)
    f = bpa; call p3dfft_btran_c2r(f, fr, 'tff'); call cut_2d_r(fr, bpa_r_z0, bpa_r_x0, bpa_r_y0)
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
                  upa_r_z0, upa_r_x0, upa_r_y0, &
                  bpa_r_z0, bpa_r_x0, bpa_r_y0, &
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

    deallocate (upa_r_z0)
    deallocate (upa_r_x0)
    deallocate (upa_r_y0)

    deallocate (bpa_r_z0)
    deallocate (bpa_r_x0)
    deallocate (bpa_r_y0)

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


  ! subroutine get_nonlinear_transfer(ntrans_aw_l_bin, ntrans_aw_g_bin, ntrans_compr_l_bin, ntrans_compr_g_bin)
    ! use fields, only: phi, psi, upa, bpa
    ! use grid, only: kx, ky, nkz
    ! use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    ! use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    ! use params, only: zi, va2cs2_plus_1
    ! use mp, only: sum_reduce
    ! implicit none
    ! real(r8), dimension (4, 1:nkpolar, nkz), intent(inout) :: &
                                            ! ntrans_aw_l_bin, ntrans_compr_l_bin, &
                                            ! ntrans_aw_g_bin, ntrans_compr_g_bin   

    ! complex(r8), allocatable, dimension(:,:,:) :: phi_l, psi_l, upa_l, bpa_l
    ! complex(r8), allocatable, dimension(:,:,:) :: phi_g, psi_g, upa_g, bpa_g
    ! complex(r8), allocatable, dimension(:,:,:) :: nlin_upe_ux_l, nlin_upe_uy_l, nlin_upe_upa_l, &
                                                  ! nlin_upe_bx_l, nlin_upe_by_l, nlin_upe_bpa_l, &
                                                  ! nlin_bpe_ux_l, nlin_bpe_uy_l, nlin_bpe_upa_l, &
                                                  ! nlin_bpe_bx_l, nlin_bpe_by_l, nlin_bpe_bpa_l
    ! complex(r8), allocatable, dimension(:,:,:) :: nlin_upe_ux_g, nlin_upe_uy_g, nlin_upe_upa_g, &
                                                  ! nlin_upe_bx_g, nlin_upe_by_g, nlin_upe_bpa_g, &
                                                  ! nlin_bpe_ux_g, nlin_bpe_uy_g, nlin_bpe_upa_g, &
                                                  ! nlin_bpe_bx_g, nlin_bpe_by_g, nlin_bpe_bpa_g
    ! complex(r8), allocatable, dimension(:,:,:) :: uxk , uyk , bxk , byk
    ! complex(r8), allocatable, dimension(:,:,:) :: uxk_, uyk_, bxk_, byk_
    ! real   (r8), allocatable, dimension(:,:,:) :: ux , uy , bx , by 
    ! complex(r8), allocatable, dimension(:,:,:,:) :: ntrans_aw_l, ntrans_compr_l
    ! complex(r8), allocatable, dimension(:,:,:) :: src_c
    ! real   (r8), allocatable, dimension(:,:,:) :: src_r
    ! complex(r8), allocatable, dimension(:,:,:,:) :: ntrans_aw_g, ntrans_compr_g
    ! integer    , dimension(:), allocatable  :: num  !count of each bin

    ! integer :: idx, i, j, k

    ! allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    ! allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    ! allocate(phi_l, source=src_c)
    ! allocate(psi_l, source=src_c)
    ! allocate(upa_l, source=src_c)
    ! allocate(bpa_l, source=src_c)
    ! allocate(phi_g, source=src_c)
    ! allocate(psi_g, source=src_c)
    ! allocate(upa_g, source=src_c)
    ! allocate(bpa_g, source=src_c)

    ! allocate(nlin_upe_ux_l , source=src_c)
    ! allocate(nlin_upe_uy_l , source=src_c)
    ! allocate(nlin_upe_upa_l, source=src_c)
    ! allocate(nlin_upe_bx_l , source=src_c)
    ! allocate(nlin_upe_by_l , source=src_c)
    ! allocate(nlin_upe_bpa_l, source=src_c)
    ! allocate(nlin_bpe_ux_l , source=src_c)
    ! allocate(nlin_bpe_uy_l , source=src_c)
    ! allocate(nlin_bpe_upa_l, source=src_c)
    ! allocate(nlin_bpe_bx_l , source=src_c)
    ! allocate(nlin_bpe_by_l , source=src_c)
    ! allocate(nlin_bpe_bpa_l, source=src_c)

    ! allocate(nlin_upe_ux_g , source=src_c)
    ! allocate(nlin_upe_uy_g , source=src_c)
    ! allocate(nlin_upe_upa_g, source=src_c)
    ! allocate(nlin_upe_bx_g , source=src_c)
    ! allocate(nlin_upe_by_g , source=src_c)
    ! allocate(nlin_upe_bpa_g, source=src_c)
    ! allocate(nlin_bpe_ux_g , source=src_c)
    ! allocate(nlin_bpe_uy_g , source=src_c)
    ! allocate(nlin_bpe_upa_g, source=src_c)
    ! allocate(nlin_bpe_bx_g , source=src_c)
    ! allocate(nlin_bpe_by_g , source=src_c)
    ! allocate(nlin_bpe_bpa_g, source=src_c)

    ! allocate(uxk , source=src_c)
    ! allocate(uyk , source=src_c)
    ! allocate(uxk_, source=src_c)
    ! allocate(uyk_, source=src_c)
    ! allocate(ux  , source=src_r)
    ! allocate(uy  , source=src_r)

    ! allocate(bxk , source=src_c)
    ! allocate(byk , source=src_c)
    ! allocate(bxk_, source=src_c)
    ! allocate(byk_, source=src_c)
    ! allocate(bx  , source=src_r)
    ! allocate(by  , source=src_r)
    ! deallocate(src_c, src_r)

    ! allocate(ntrans_aw_l   (4, ikx_st:ikx_en, iky_st:iky_en, ikz_st:ikz_en), source=(0.d0,0.d0))
    ! allocate(ntrans_aw_g   (4, ikx_st:ikx_en, iky_st:iky_en, ikz_st:ikz_en), source=(0.d0,0.d0))
    ! allocate(ntrans_compr_l(4, ikx_st:ikx_en, iky_st:iky_en, ikz_st:ikz_en), source=(0.d0,0.d0))
    ! allocate(ntrans_compr_g(4, ikx_st:ikx_en, iky_st:iky_en, ikz_st:ikz_en), source=(0.d0,0.d0))

    ! allocate(num (1:nkpolar))
    ! num = 0


    ! ! 1. First get un-filtered ux, uy, bx, and by
    ! do j = iky_st, iky_en
      ! do k = ikz_st, ikz_en
        ! do i = ikx_st, ikx_en
          ! uyk(i, k, j) = -zi*kx(i)*phi(i, k, j)
          ! uxk(i, k, j) =  zi*ky(j)*phi(i, k, j)
          ! byk(i, k, j) = -zi*kx(i)*psi(i, k, j)
          ! bxk(i, k, j) =  zi*ky(j)*psi(i, k, j)
        ! enddo
      ! enddo
    ! enddo
    ! uxk_ = uxk
    ! uyk_ = uyk
    ! bxk_ = bxk
    ! byk_ = byk
    ! call p3dfft_btran_c2r(uxk_, ux, 'tff')
    ! call p3dfft_btran_c2r(uyk_, uy, 'tff')
    ! call p3dfft_btran_c2r(bxk_, bx, 'tff')
    ! call p3dfft_btran_c2r(byk_, by, 'tff')

    ! ! 2. Loop for kpbin
    ! do idx = 1, nkpolar - 1
      ! ! 3. Filter phi and psi for |k| < k and |k| > k
      ! do j = iky_st, iky_en
        ! do k = ikz_st, ikz_en
          ! do i = ikx_st, ikx_en
            ! if (kx(i)**2 + ky(j)**2 < kpbin(idx)**2) then
              ! phi_l(i, k, j) = phi(i, k, j)
              ! psi_l(i, k, j) = psi(i, k, j)
              ! upa_l(i, k, j) = upa(i, k, j)
              ! bpa_l(i, k, j) = bpa(i, k, j)
            ! endif
            ! if (kx(i)**2 + ky(j)**2 > kpbin(idx)**2) then
              ! phi_g(i, k, j) = phi(i, k, j)
              ! psi_g(i, k, j) = psi(i, k, j)
              ! upa_g(i, k, j) = upa(i, k, j)
              ! bpa_g(i, k, j) = bpa(i, k, j)
            ! endif
          ! enddo
        ! enddo
      ! enddo
      ! ! 4. get  low filtered nonlinear terms e.g., {phi, omg^<k}
      ! call get_nonlinear_terms(ux, uy, bx, by, &
                               ! phi_l, psi_l, upa_l, bpa_l, &
                               ! nlin_upe_ux_l, nlin_upe_uy_l, nlin_upe_upa_l, nlin_upe_bx_l, nlin_upe_by_l, nlin_upe_bpa_l, &
                               ! nlin_bpe_ux_l, nlin_bpe_uy_l, nlin_bpe_upa_l, nlin_bpe_bx_l, nlin_bpe_by_l, nlin_bpe_bpa_l)
      ! ! 4. get high filtered nonlinear terms e.g., {phi, omg^>k}
      ! call get_nonlinear_terms(ux, uy, bx, by, &
                               ! phi_g, psi_g, upa_g, bpa_g, &
                               ! nlin_upe_ux_g, nlin_upe_uy_g, nlin_upe_upa_g, nlin_upe_bx_g, nlin_upe_by_g, nlin_upe_bpa_g, &
                               ! nlin_bpe_ux_g, nlin_bpe_uy_g, nlin_bpe_upa_g, nlin_bpe_bx_g, nlin_bpe_by_g, nlin_bpe_bpa_g)
      ! ! 5. convolution with un-filtered fields and fileted nonlinear terms, e.g., phi^* F[{phi, omg^>k}]
      ! do j = iky_st, iky_en
        ! do k = ikz_st, ikz_en
          ! do i = ikx_st, ikx_en
            ! ntrans_aw_l   (1, i, k, j) = 0.5d0*( &
                                          ! - (conjg(uxk(i,k,j))*nlin_upe_ux_l(i,k,j) + uxk(i,k,j)*conjg(nlin_upe_ux_l(i,k,j))) & 
                                          ! - (conjg(uyk(i,k,j))*nlin_upe_uy_l(i,k,j) + uyk(i,k,j)*conjg(nlin_upe_uy_l(i,k,j))) & 
                                        ! )
            ! ntrans_aw_l   (2, i, k, j) = 0.5d0*( &
                                          ! + (conjg(uxk(i,k,j))*nlin_bpe_bx_l(i,k,j) + uxk(i,k,j)*conjg(nlin_bpe_bx_l(i,k,j))) & 
                                          ! + (conjg(uyk(i,k,j))*nlin_bpe_by_l(i,k,j) + uyk(i,k,j)*conjg(nlin_bpe_by_l(i,k,j))) & 
                                        ! )
            ! ntrans_aw_l   (3, i, k, j) = 0.5d0*( &
                                          ! - (conjg(bxk(i,k,j))*nlin_upe_bx_l(i,k,j) + bxk(i,k,j)*conjg(nlin_upe_bx_l(i,k,j))) & 
                                          ! - (conjg(byk(i,k,j))*nlin_upe_by_l(i,k,j) + byk(i,k,j)*conjg(nlin_upe_by_l(i,k,j))) & 
                                        ! )
            ! ntrans_aw_l   (4, i, k, j) = 0.5d0*( &
                                          ! + (conjg(bxk(i,k,j))*nlin_bpe_ux_l(i,k,j) + bxk(i,k,j)*conjg(nlin_bpe_ux_l(i,k,j))) & 
                                          ! + (conjg(byk(i,k,j))*nlin_bpe_uy_l(i,k,j) + byk(i,k,j)*conjg(nlin_bpe_uy_l(i,k,j))) & 
                                        ! )
            ! ntrans_aw_g   (1, i, k, j) = 0.5d0*( &
                                          ! - (conjg(uxk(i,k,j))*nlin_upe_ux_g(i,k,j) + uxk(i,k,j)*conjg(nlin_upe_ux_g(i,k,j))) & 
                                          ! - (conjg(uyk(i,k,j))*nlin_upe_uy_g(i,k,j) + uyk(i,k,j)*conjg(nlin_upe_uy_g(i,k,j))) & 
                                        ! )
            ! ntrans_aw_g   (2, i, k, j) = 0.5d0*( &
                                          ! + (conjg(uxk(i,k,j))*nlin_bpe_bx_g(i,k,j) + uxk(i,k,j)*conjg(nlin_bpe_bx_g(i,k,j))) & 
                                          ! + (conjg(uyk(i,k,j))*nlin_bpe_by_g(i,k,j) + uyk(i,k,j)*conjg(nlin_bpe_by_g(i,k,j))) & 
                                        ! )
            ! ntrans_aw_g   (3, i, k, j) = 0.5d0*( &
                                          ! - (conjg(bxk(i,k,j))*nlin_upe_bx_g(i,k,j) + bxk(i,k,j)*conjg(nlin_upe_bx_g(i,k,j))) & 
                                          ! - (conjg(byk(i,k,j))*nlin_upe_by_g(i,k,j) + byk(i,k,j)*conjg(nlin_upe_by_g(i,k,j))) & 
                                        ! )
            ! ntrans_aw_g   (4, i, k, j) = 0.5d0*( &
                                          ! + (conjg(bxk(i,k,j))*nlin_bpe_ux_g(i,k,j) + bxk(i,k,j)*conjg(nlin_bpe_ux_g(i,k,j))) & 
                                          ! + (conjg(byk(i,k,j))*nlin_bpe_uy_g(i,k,j) + byk(i,k,j)*conjg(nlin_bpe_uy_g(i,k,j))) & 
                                        ! )
            ! ntrans_compr_l(1, i, k, j) = 0.5d0*( &
                                          ! - (conjg(upa(i,k,j))*nlin_upe_upa_l(i,k,j) + upa(i,k,j)*conjg(nlin_upe_upa_l(i,k,j))) & 
                                        ! )
            ! ntrans_compr_l(2, i, k, j) = 0.5d0*( &
                                          ! + (conjg(upa(i,k,j))*nlin_bpe_bpa_l(i,k,j) + upa(i,k,j)*conjg(nlin_bpe_bpa_l(i,k,j))) & 
                                        ! )
            ! ntrans_compr_l(3, i, k, j) = 0.5d0*va2cs2_plus_1*( &
                                          ! - (conjg(bpa(i,k,j))*nlin_upe_bpa_l(i,k,j) + bpa(i,k,j)*conjg(nlin_upe_bpa_l(i,k,j))) & 
                                        ! )
            ! ntrans_compr_l(4, i, k, j) = 0.5d0*va2cs2_plus_1*( &
                                          ! + (conjg(bpa(i,k,j))*nlin_bpe_upa_l(i,k,j) + bpa(i,k,j)*conjg(nlin_bpe_upa_l(i,k,j))) & 
                                        ! )
            ! ntrans_compr_g(1, i, k, j) = 0.5d0*( &
                                          ! - (conjg(upa(i,k,j))*nlin_upe_upa_g(i,k,j) + upa(i,k,j)*conjg(nlin_upe_upa_g(i,k,j))) & 
                                        ! )
            ! ntrans_compr_g(2, i, k, j) = 0.5d0*( &
                                          ! + (conjg(upa(i,k,j))*nlin_bpe_bpa_g(i,k,j) + upa(i,k,j)*conjg(nlin_bpe_bpa_g(i,k,j))) & 
                                        ! )
            ! ntrans_compr_g(3, i, k, j) = 0.5d0*va2cs2_plus_1*( &
                                          ! - (conjg(bpa(i,k,j))*nlin_upe_bpa_g(i,k,j) + bpa(i,k,j)*conjg(nlin_upe_bpa_g(i,k,j))) & 
                                        ! )
            ! ntrans_compr_g(4, i, k, j) = 0.5d0*va2cs2_plus_1*( &
                                          ! + (conjg(bpa(i,k,j))*nlin_bpe_upa_g(i,k,j) + bpa(i,k,j)*conjg(nlin_bpe_upa_g(i,k,j))) & 
                                        ! )

            ! ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
            ! !-----------------------------------------------------------------------------------------------------------------------------------
            ! ! The volume integral of a quadratic function is
            ! ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
            ! ! |f_{k_x, k_y}|^2 = \left( \sum_{k_x = -n_{k_x}/2}^{-1}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} + \sum_{k_x = 0}
            ! ! \sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} + \sum_{k_x = 1}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} \right) |f_{k_x, k_y}|^2\]
            ! ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
            ! !-----------------------------------------------------------------------------------------------------------------------------------
            ! if (j /= 1) then
              ! ntrans_aw_l   (:, i, k, j) = 2.0d0*ntrans_aw_l   (:, i, k, j)
              ! ntrans_aw_g   (:, i, k, j) = 2.0d0*ntrans_aw_g   (:, i, k, j)
              ! ntrans_compr_l(:, i, k, j) = 2.0d0*ntrans_compr_l(:, i, k, j)
              ! ntrans_compr_g(:, i, k, j) = 2.0d0*ntrans_compr_g(:, i, k, j)
            ! endif

          ! end do
        ! end do
      ! end do

      ! ! 6. take kperp-bin
      ! do j = iky_st, iky_en
        ! do k = ikz_st, ikz_en
          ! do i = ikx_st, ikx_en
            ! if (kx(i)**2 + ky(j)**2 >= kpbin(idx)**2 .and. kx(i)**2 + ky(j)**2 < kpbin(idx+1)**2) then
              ! ntrans_aw_l_bin   (:, idx, k) = ntrans_aw_l_bin   (:, idx, k) + real(ntrans_aw_l   (:, i, j, k))
              ! ntrans_aw_g_bin   (:, idx, k) = ntrans_aw_g_bin   (:, idx, k) + real(ntrans_aw_g   (:, i, j, k))
              ! ntrans_compr_l_bin(:, idx, k) = ntrans_compr_l_bin(:, idx, k) + real(ntrans_compr_l(:, i, j, k))
              ! ntrans_compr_g_bin(:, idx, k) = ntrans_compr_g_bin(:, idx, k) + real(ntrans_compr_g(:, i, j, k))
              ! num (idx) = num (idx) + 1
            ! endif
          ! enddo
        ! enddo
      ! enddo
    ! enddo

    ! ! 7. combine all procs
    ! call sum_reduce(ntrans_aw_l_bin   , 0)
    ! call sum_reduce(ntrans_aw_g_bin   , 0)
    ! call sum_reduce(ntrans_compr_l_bin, 0)
    ! call sum_reduce(ntrans_compr_g_bin, 0)
    ! call sum_reduce(num , 0)

    ! deallocate(phi_l)
    ! deallocate(psi_l)
    ! deallocate(upa_l)
    ! deallocate(bpa_l)
    ! deallocate(phi_g)
    ! deallocate(psi_g)
    ! deallocate(upa_g)
    ! deallocate(bpa_g)
    ! deallocate(nlin_upe_ux_l, nlin_upe_uy_l, nlin_upe_upa_l, &
               ! nlin_upe_bx_l, nlin_upe_by_l, nlin_upe_bpa_l, &
               ! nlin_bpe_ux_l, nlin_bpe_uy_l, nlin_bpe_upa_l, &
               ! nlin_bpe_bx_l, nlin_bpe_by_l, nlin_bpe_bpa_l, &
               ! nlin_upe_ux_g, nlin_upe_uy_g, nlin_upe_upa_g, &
               ! nlin_upe_bx_g, nlin_upe_by_g, nlin_upe_bpa_g, &
               ! nlin_bpe_ux_g, nlin_bpe_uy_g, nlin_bpe_upa_g, &
               ! nlin_bpe_bx_g, nlin_bpe_by_g, nlin_bpe_bpa_g)
    ! deallocate(uxk , uyk , bxk , byk )
    ! deallocate(uxk_, uyk_, bxk_, byk_)
    ! deallocate(ux , uy , bx , by )
    ! deallocate(ntrans_aw_l)
    ! deallocate(ntrans_aw_g)
    ! deallocate(ntrans_compr_l)
    ! deallocate(ntrans_compr_g)
  ! end subroutine get_nonlinear_transfer

  ! subroutine get_nonlinear_terms(ux, uy, bx, by, &
                                 ! phi_, psi_, upa_, bpa_, &
                                 ! nlin_upe_ux_, nlin_upe_uy_, nlin_upe_upa_, nlin_upe_bx_, nlin_upe_by_, nlin_upe_bpa_, &
                                 ! nlin_bpe_ux_, nlin_bpe_uy_, nlin_bpe_upa_, nlin_bpe_bx_, nlin_bpe_by_, nlin_bpe_bpa_)
    ! use p3dfft
    ! use grid, only: nlx, nly, nlz, kx, ky
    ! use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    ! use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    ! use params, only: zi, va2cs2_plus_1, q
    ! use mpienv, only: ph, sp
    ! implicit none
    ! complex(r8), allocatable, dimension(:,:,:) :: dux_dxk_ , dux_dyk_
    ! real   (r8), allocatable, dimension(:,:,:) :: dux_dx_  , dux_dy_ 
    ! complex(r8), allocatable, dimension(:,:,:) :: duy_dxk_ , duy_dyk_
    ! real   (r8), allocatable, dimension(:,:,:) :: duy_dx_  , duy_dy_ 
    ! complex(r8), allocatable, dimension(:,:,:) :: dupa_dxk_, dupa_dyk_
    ! real   (r8), allocatable, dimension(:,:,:) :: dupa_dx_ , dupa_dy_ 

    ! complex(r8), allocatable, dimension(:,:,:) :: dbx_dxk_ , dbx_dyk_
    ! real   (r8), allocatable, dimension(:,:,:) :: dbx_dx_  , dbx_dy_ 
    ! complex(r8), allocatable, dimension(:,:,:) :: dby_dxk_ , dby_dyk_
    ! real   (r8), allocatable, dimension(:,:,:) :: dby_dx_  , dby_dy_ 
    ! complex(r8), allocatable, dimension(:,:,:) :: dbpa_dxk_, dbpa_dyk_
    ! real   (r8), allocatable, dimension(:,:,:) :: dbpa_dx_ , dbpa_dy_ 

    ! real   (r8), allocatable, dimension(:,:,:) :: nlin_upe_ux_r, nlin_upe_uy_r, nlin_upe_upa_r 
    ! real   (r8), allocatable, dimension(:,:,:) :: nlin_upe_bx_r, nlin_upe_by_r, nlin_upe_bpa_r 
    ! real   (r8), allocatable, dimension(:,:,:) :: nlin_bpe_ux_r, nlin_bpe_uy_r, nlin_bpe_upa_r 
    ! real   (r8), allocatable, dimension(:,:,:) :: nlin_bpe_bx_r, nlin_bpe_by_r, nlin_bpe_bpa_r 

    ! real   (r8), dimension (ily_st:ily_en, &
                            ! ilz_st:ilz_en, &
                            ! ilx_st:ilx_en), intent(in) :: &
                                    ! ux, uy, bx, by
    ! complex(r8), dimension (ikx_st:ikx_en, &
                            ! ikz_st:ikz_en, &
                            ! iky_st:iky_en), intent(in) :: &
                                    ! phi_, psi_, upa_, bpa_
    ! complex(r8), dimension (ikx_st:ikx_en, &
                            ! ikz_st:ikz_en, &
                            ! iky_st:iky_en), intent(inout) :: &
                                    ! nlin_upe_ux_, nlin_upe_uy_, nlin_upe_upa_, nlin_upe_bx_, nlin_upe_by_, nlin_upe_bpa_, &
                                    ! nlin_bpe_ux_, nlin_bpe_uy_, nlin_bpe_upa_, nlin_bpe_bx_, nlin_bpe_by_, nlin_bpe_bpa_

    ! complex(r8) :: ux_, uy_, bx_, by_
    ! integer :: i, j, k

    ! call alloc_z(dux_dxk_ , sp, opt_global=.true.); dux_dxk_  = 0._mytype
    ! call alloc_z(dux_dyk_ , sp, opt_global=.true.); dux_dyk_  = 0._mytype
    ! call alloc_x(dux_dx_  , ph, opt_global=.true.); dux_dx_   = 0._mytype
    ! call alloc_x(dux_dy_  , ph, opt_global=.true.); dux_dy_   = 0._mytype
                                                              
    ! call alloc_z(duy_dxk_ , sp, opt_global=.true.); duy_dxk_  = 0._mytype
    ! call alloc_z(duy_dyk_ , sp, opt_global=.true.); duy_dyk_  = 0._mytype
    ! call alloc_x(duy_dx_  , ph, opt_global=.true.); duy_dx_   = 0._mytype
    ! call alloc_x(duy_dy_  , ph, opt_global=.true.); duy_dy_   = 0._mytype

    ! call alloc_z(dupa_dxk_, sp, opt_global=.true.); dupa_dxk_ = 0._mytype
    ! call alloc_z(dupa_dyk_, sp, opt_global=.true.); dupa_dyk_ = 0._mytype
    ! call alloc_x(dupa_dx_ , ph, opt_global=.true.); dupa_dx_  = 0._mytype
    ! call alloc_x(dupa_dy_ , ph, opt_global=.true.); dupa_dy_  = 0._mytype

    ! call alloc_z(dbx_dxk_ , sp, opt_global=.true.); dbx_dxk_  = 0._mytype
    ! call alloc_z(dbx_dyk_ , sp, opt_global=.true.); dbx_dyk_  = 0._mytype
    ! call alloc_x(dbx_dx_  , ph, opt_global=.true.); dbx_dx_   = 0._mytype
    ! call alloc_x(dbx_dy_  , ph, opt_global=.true.); dbx_dy_   = 0._mytype
                                                              
    ! call alloc_z(dby_dxk_ , sp, opt_global=.true.); dby_dxk_  = 0._mytype
    ! call alloc_z(dby_dyk_ , sp, opt_global=.true.); dby_dyk_  = 0._mytype
    ! call alloc_x(dby_dx_  , ph, opt_global=.true.); dby_dx_   = 0._mytype
    ! call alloc_x(dby_dy_  , ph, opt_global=.true.); dby_dy_   = 0._mytype

    ! call alloc_z(dbpa_dxk_, sp, opt_global=.true.); dbpa_dxk_ = 0._mytype
    ! call alloc_z(dbpa_dyk_, sp, opt_global=.true.); dbpa_dyk_ = 0._mytype
    ! call alloc_x(dbpa_dx_ , ph, opt_global=.true.); dbpa_dx_  = 0._mytype
    ! call alloc_x(dbpa_dy_ , ph, opt_global=.true.); dbpa_dy_  = 0._mytype

    ! call alloc_x(nlin_upe_ux_r , ph, opt_global=.true.); nlin_upe_ux_r  = 0._mytype
    ! call alloc_x(nlin_upe_uy_r , ph, opt_global=.true.); nlin_upe_uy_r  = 0._mytype
    ! call alloc_x(nlin_upe_upa_r, ph, opt_global=.true.); nlin_upe_upa_r = 0._mytype
    ! call alloc_x(nlin_upe_bx_r , ph, opt_global=.true.); nlin_upe_bx_r  = 0._mytype
    ! call alloc_x(nlin_upe_by_r , ph, opt_global=.true.); nlin_upe_by_r  = 0._mytype
    ! call alloc_x(nlin_upe_bpa_r, ph, opt_global=.true.); nlin_upe_bpa_r = 0._mytype
    ! call alloc_x(nlin_bpe_ux_r , ph, opt_global=.true.); nlin_bpe_ux_r  = 0._mytype
    ! call alloc_x(nlin_bpe_uy_r , ph, opt_global=.true.); nlin_bpe_uy_r  = 0._mytype
    ! call alloc_x(nlin_bpe_upa_r, ph, opt_global=.true.); nlin_bpe_upa_r = 0._mytype
    ! call alloc_x(nlin_bpe_bx_r , ph, opt_global=.true.); nlin_bpe_bx_r  = 0._mytype
    ! call alloc_x(nlin_bpe_by_r , ph, opt_global=.true.); nlin_bpe_by_r  = 0._mytype
    ! call alloc_x(nlin_bpe_bpa_r, ph, opt_global=.true.); nlin_bpe_bpa_r = 0._mytype

    ! ! 1. Calculate grad in Fourier space
    ! do k = ikz_st, ikz_en
      ! do j = iky_st, iky_en
        ! do i = ikx_st, ikx_en
          ! uy_ = -zi*kx(i)*phi_(i, j, k)
          ! ux_ =  zi*ky(j)*phi_(i, j, k)
          ! by_ = -zi*kx(i)*psi_(i, j, k)
          ! bx_ =  zi*ky(j)*psi_(i, j, k)

          ! dux_dxk_ (i, j, k) = zi*kx(i)*ux_
          ! dux_dyk_ (i, j, k) = zi*ky(j)*ux_
          ! duy_dxk_ (i, j, k) = zi*kx(i)*uy_
          ! duy_dyk_ (i, j, k) = zi*ky(j)*uy_
          ! dupa_dxk_(i, j, k) = zi*kx(i)*upa_(i, j, k)
          ! dupa_dyk_(i, j, k) = zi*ky(j)*upa_(i, j, k)

          ! dbx_dxk_ (i, j, k) = zi*kx(i)*bx_
          ! dbx_dyk_ (i, j, k) = zi*ky(j)*bx_
          ! dby_dxk_ (i, j, k) = zi*kx(i)*by_
          ! dby_dyk_ (i, j, k) = zi*ky(j)*by_
          ! dbpa_dxk_(i, j, k) = zi*kx(i)*bpa_(i, j, k)
          ! dbpa_dyk_(i, j, k) = zi*ky(j)*bpa_(i, j, k)
        ! enddo
      ! enddo
    ! enddo

    ! ! 2. Inverse FFT
    ! call decomp_2d_fft_3d(dux_dxk_ , dux_dx_ )
    ! call decomp_2d_fft_3d(dux_dyk_ , dux_dy_ )
    ! call decomp_2d_fft_3d(duy_dxk_ , duy_dx_ )
    ! call decomp_2d_fft_3d(duy_dyk_ , duy_dy_ )
    ! call decomp_2d_fft_3d(dupa_dxk_, dupa_dx_)
    ! call decomp_2d_fft_3d(dupa_dyk_, dupa_dy_)

    ! call decomp_2d_fft_3d(dbx_dxk_ , dbx_dx_ )
    ! call decomp_2d_fft_3d(dbx_dyk_ , dbx_dy_ )
    ! call decomp_2d_fft_3d(dby_dxk_ , dby_dx_ )
    ! call decomp_2d_fft_3d(dby_dyk_ , dby_dy_ )
    ! call decomp_2d_fft_3d(dbpa_dxk_, dbpa_dx_)
    ! call decomp_2d_fft_3d(dbpa_dyk_, dbpa_dy_)

    ! ! 3. Calculate poisson brackets in real space
    ! do k = ilz_st, ilz_en
      ! do j = ily_st, ily_en
        ! do i = ilx_st, ilx_en
          ! nlin_upe_ux_r (i, j, k) =  ux(i, j, k)* dux_dx_(i, j, k) &
                                   ! + uy(i, j, k)* dux_dy_(i, j, k)
          ! nlin_upe_uy_r (i, j, k) =  ux(i, j, k)* duy_dx_(i, j, k) &
                                   ! + uy(i, j, k)* duy_dy_(i, j, k)
          ! nlin_upe_upa_r(i, j, k) =  ux(i, j, k)*dupa_dx_(i, j, k) &
                                   ! + uy(i, j, k)*dupa_dy_(i, j, k)

          ! nlin_upe_bx_r (i, j, k) =  ux(i, j, k)* dbx_dx_(i, j, k) &
                                   ! + uy(i, j, k)* dbx_dy_(i, j, k)
          ! nlin_upe_by_r (i, j, k) =  ux(i, j, k)* dby_dx_(i, j, k) &
                                   ! + uy(i, j, k)* dby_dy_(i, j, k)
          ! nlin_upe_bpa_r(i, j, k) =  ux(i, j, k)*dbpa_dx_(i, j, k) &
                                   ! + uy(i, j, k)*dbpa_dy_(i, j, k)

          ! nlin_bpe_ux_r (i, j, k) =  bx(i, j, k)* dux_dx_(i, j, k) &
                                   ! + by(i, j, k)* dux_dy_(i, j, k)
          ! nlin_bpe_uy_r (i, j, k) =  bx(i, j, k)* duy_dx_(i, j, k) &
                                   ! + by(i, j, k)* duy_dy_(i, j, k)
          ! nlin_bpe_upa_r(i, j, k) =  bx(i, j, k)*dupa_dx_(i, j, k) &
                                   ! + by(i, j, k)*dupa_dy_(i, j, k)

          ! nlin_bpe_bx_r (i, j, k) =  bx(i, j, k)* dbx_dx_(i, j, k) &
                                   ! + by(i, j, k)* dbx_dy_(i, j, k)
          ! nlin_bpe_by_r (i, j, k) =  bx(i, j, k)* dby_dx_(i, j, k) &
                                   ! + by(i, j, k)* dby_dy_(i, j, k)
          ! nlin_bpe_bpa_r(i, j, k) =  bx(i, j, k)*dbpa_dx_(i, j, k) &
                                   ! + by(i, j, k)*dbpa_dy_(i, j, k)
        ! enddo
      ! enddo
    ! enddo

    ! ! 4. Forward FFT
    ! call decomp_2d_fft_3d(nlin_upe_ux_r , nlin_upe_ux_ ); nlin_upe_ux_  = nlin_upe_ux_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_upe_uy_r , nlin_upe_uy_ ); nlin_upe_uy_  = nlin_upe_uy_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_upe_upa_r, nlin_upe_upa_); nlin_upe_upa_ = nlin_upe_upa_/nlx/nly/nlz

    ! call decomp_2d_fft_3d(nlin_upe_bx_r , nlin_upe_bx_ ); nlin_upe_bx_  = nlin_upe_bx_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_upe_by_r , nlin_upe_by_ ); nlin_upe_by_  = nlin_upe_by_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_upe_bpa_r, nlin_upe_bpa_); nlin_upe_bpa_ = nlin_upe_bpa_/nlx/nly/nlz

    ! call decomp_2d_fft_3d(nlin_bpe_ux_r , nlin_bpe_ux_ ); nlin_bpe_ux_  = nlin_bpe_ux_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_bpe_uy_r , nlin_bpe_uy_ ); nlin_bpe_uy_  = nlin_bpe_uy_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_bpe_upa_r, nlin_bpe_upa_); nlin_bpe_upa_ = nlin_bpe_upa_/nlx/nly/nlz

    ! call decomp_2d_fft_3d(nlin_bpe_bx_r , nlin_bpe_bx_ ); nlin_bpe_bx_  = nlin_bpe_bx_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_bpe_by_r , nlin_bpe_by_ ); nlin_bpe_by_  = nlin_bpe_by_ /nlx/nly/nlz
    ! call decomp_2d_fft_3d(nlin_bpe_bpa_r, nlin_bpe_bpa_); nlin_bpe_bpa_ = nlin_bpe_bpa_/nlx/nly/nlz

    ! deallocate(dux_dxk_ , dux_dyk_)
    ! deallocate(dux_dx_  , dux_dy_ )
    ! deallocate(duy_dxk_ , duy_dyk_)
    ! deallocate(duy_dx_  , duy_dy_ )
    ! deallocate(dupa_dxk_, dupa_dyk_)
    ! deallocate(dupa_dx_ , dupa_dy_ )

    ! deallocate(dbx_dxk_ , dbx_dyk_)
    ! deallocate(dbx_dx_  , dbx_dy_ )
    ! deallocate(dby_dxk_ , dby_dyk_)
    ! deallocate(dby_dx_  , dby_dy_ )
    ! deallocate(dbpa_dxk_, dbpa_dyk_)
    ! deallocate(dbpa_dx_ , dbpa_dy_ )

    ! deallocate(nlin_upe_ux_r, nlin_upe_uy_r, nlin_upe_upa_r)
    ! deallocate(nlin_upe_bx_r, nlin_upe_by_r, nlin_upe_bpa_r)
    ! deallocate(nlin_bpe_ux_r, nlin_bpe_uy_r, nlin_bpe_upa_r)
    ! deallocate(nlin_bpe_bx_r, nlin_bpe_by_r, nlin_bpe_bpa_r)

  ! end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate shell-to-shell transfer
!-----------------------------------------------!
  subroutine loop_diagnostics_nltrans
  ! under development...
  end subroutine loop_diagnostics_nltrans


end module diagnostics

