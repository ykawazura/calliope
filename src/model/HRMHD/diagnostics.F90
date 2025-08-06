!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Diagnostics for HRMHD
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
    use params, only: zi, nonlinear, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
    use force, only: fphi, fpsi, fupa, fbpa, fphi_old, fpsi_old, fupa_old, fbpa_old
    use utils, only: cabs2
    implicit none
    integer :: i, j, k

    real(r8), allocatable, dimension(:,:,:) :: upe2, ux2, uy2, bpe2, bx2, by2, upa2, bpa2
    real(r8), allocatable, dimension(:,:,:) :: upe2old, bpe2old, upa2old, bpa2old
    real(r8), allocatable, dimension(:,:,:) :: upe2dissip_x, upe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpe2dissip_x, bpe2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: upa2dissip_x, upa2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpa2dissip_x, bpa2dissip_z
    real(r8), allocatable, dimension(:,:,:) :: upe2_KAW_dissip_x, upe2_KAW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpe2_KAW_dissip_x, bpe2_KAW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: upa2_KAW_dissip_x, upa2_KAW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpa2_KAW_dissip_x, bpa2_KAW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: upe2_ICW_dissip_x, upe2_ICW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpe2_ICW_dissip_x, bpe2_ICW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: upa2_ICW_dissip_x, upa2_ICW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: bpa2_ICW_dissip_x, bpa2_ICW_dissip_z
    real(r8), allocatable, dimension(:,:,:) :: p_aw, p_compr, p_xhl
    real(r8), allocatable, dimension(:,:,:) :: zppe2, zmpe2, zppa2, zmpa2
    real(r8), allocatable, dimension(:,:,:) :: upe2_KAW, bpe2_KAW, upa2_KAW, bpa2_KAW
    real(r8), allocatable, dimension(:,:,:) :: upe2_ICW, bpe2_ICW, upa2_ICW, bpa2_ICW
    real(r8), allocatable, dimension(:,:,:) :: src

    complex(r8), allocatable, dimension(:,:,:) :: src_c
    complex(r8), allocatable, dimension(:,:,:) :: phi_KAW, psi_KAW, upa_KAW, bpa_KAW
    complex(r8), allocatable, dimension(:,:,:) :: phi_ICW, psi_ICW, upa_ICW, bpa_ICW

    real(r8) :: upe2_sum, bpe2_sum, upa2_sum, bpa2_sum
    real(r8) :: upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum
    real(r8) :: upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum
    real(r8) :: upe2_KAW_dissip_sum, bpe2_KAW_dissip_sum, upa2_KAW_dissip_sum, bpa2_KAW_dissip_sum
    real(r8) :: upe2_ICW_dissip_sum, bpe2_ICW_dissip_sum, upa2_ICW_dissip_sum, bpa2_ICW_dissip_sum
    real(r8) :: p_aw_sum, p_compr_sum, p_xhl_sum
    real(r8) :: zppe2_sum, zmpe2_sum, zppa2_sum, zmpa2_sum

    real(r8), dimension(:, :), allocatable :: upe2_bin, bpe2_bin, upa2_bin, bpa2_bin     ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: ux2_bin , uy2_bin , bx2_bin , by2_bin      ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: p_aw_bin, p_compr_bin                      ! [kprp, kz]
    real(r8), dimension(:, :, :), allocatable :: ntrans_aw_l_bin   , ntrans_aw_g_bin     ! [4, kprp, kz]
                                                                                         !   1 : -upe.(upe.grad upe) 
                                                                                         !   2 : +upe.(bpe.grad bpe) 
                                                                                         !   3 : -bpe.(upe.grad bpe) 
                                                                                         !   4 : +bpe.(bpe.grad upe) 
    real(r8), dimension(:, :, :), allocatable :: ntrans_compr_l_bin, ntrans_compr_g_bin  ! [4, kprp, kz]
                                                                                         !   1 : -upa (upe.grad upa) 
                                                                                         !   2 : +upa (bpe.grad bpa) 
                                                                                         !   3 : -bpa (upe.grad bpa) 
                                                                                         !   4 : +bpa (bpe.grad upa) 
    real(r8), dimension(:, :), allocatable :: dissip_aw_bin, dissip_compr_bin            ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: dissip_KAW_bin, dissip_ICW_bin             ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: zppe2_bin, zmpe2_bin, zppa2_bin, zmpa2_bin ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: upe2_KAW_bin, bpe2_KAW_bin                 ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: upa2_KAW_bin, bpa2_KAW_bin                 ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: upe2_ICW_bin, bpe2_ICW_bin                 ! [kprp, kz]
    real(r8), dimension(:, :), allocatable :: upa2_ICW_bin, bpa2_ICW_bin                 ! [kprp, kz]
    complex(r8) :: phi_mid, psi_mid, jpa_mid, upa_mid, bpa_mid, fphi_mid, fpsi_mid, fupa_mid, fbpa_mid
    complex(r8) :: zppe_mid, zmpe_mid, fzppe_mid, fzmpe_mid

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
    allocate(p_xhl       , source=src)
    allocate(zppe2       , source=src)
    allocate(zmpe2       , source=src)
    allocate(zppa2       , source=src)
    allocate(zmpa2       , source=src)
    allocate(upe2_KAW    , source=src)
    allocate(bpe2_KAW    , source=src)
    allocate(upa2_KAW    , source=src)
    allocate(bpa2_KAW    , source=src)
    allocate(upe2_ICW    , source=src)
    allocate(bpe2_ICW    , source=src)
    allocate(upa2_ICW    , source=src)
    allocate(bpa2_ICW    , source=src)
    allocate(upe2_KAW_dissip_x, source=src)
    allocate(upe2_KAW_dissip_z, source=src)
    allocate(bpe2_KAW_dissip_x, source=src)
    allocate(bpe2_KAW_dissip_z, source=src)
    allocate(upa2_KAW_dissip_x, source=src)
    allocate(upa2_KAW_dissip_z, source=src)
    allocate(bpa2_KAW_dissip_x, source=src)
    allocate(bpa2_KAW_dissip_z, source=src)
    allocate(upe2_ICW_dissip_x, source=src)
    allocate(upe2_ICW_dissip_z, source=src)
    allocate(bpe2_ICW_dissip_x, source=src)
    allocate(bpe2_ICW_dissip_z, source=src)
    allocate(upa2_ICW_dissip_x, source=src)
    allocate(upa2_ICW_dissip_z, source=src)
    allocate(bpa2_ICW_dissip_x, source=src)
    allocate(bpa2_ICW_dissip_z, source=src)
    deallocate(src)

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(phi_KAW, source=src_c)
    allocate(psi_KAW, source=src_c)
    allocate(upa_KAW, source=src_c)
    allocate(bpa_KAW, source=src_c)
    allocate(phi_ICW, source=src_c)
    allocate(psi_ICW, source=src_c)
    allocate(upa_ICW, source=src_c)
    allocate(bpa_ICW, source=src_c)
    deallocate(src_c)

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
    allocate (dissip_KAW_bin     (1:nkpolar, nkz)); dissip_KAW_bin        = 0.d0
    allocate (dissip_ICW_bin     (1:nkpolar, nkz)); dissip_ICW_bin        = 0.d0
    allocate (zppe2_bin          (1:nkpolar, nkz)); zppe2_bin             = 0.d0
    allocate (zmpe2_bin          (1:nkpolar, nkz)); zmpe2_bin             = 0.d0
    allocate (zppa2_bin          (1:nkpolar, nkz)); zppa2_bin             = 0.d0
    allocate (zmpa2_bin          (1:nkpolar, nkz)); zmpa2_bin             = 0.d0
    allocate (upe2_KAW_bin       (1:nkpolar, nkz)); upe2_KAW_bin          = 0.d0
    allocate (bpe2_KAW_bin       (1:nkpolar, nkz)); bpe2_KAW_bin          = 0.d0
    allocate (upa2_KAW_bin       (1:nkpolar, nkz)); upa2_KAW_bin          = 0.d0
    allocate (bpa2_KAW_bin       (1:nkpolar, nkz)); bpa2_KAW_bin          = 0.d0
    allocate (upe2_ICW_bin       (1:nkpolar, nkz)); upe2_ICW_bin          = 0.d0
    allocate (bpe2_ICW_bin       (1:nkpolar, nkz)); bpe2_ICW_bin          = 0.d0
    allocate (upa2_ICW_bin       (1:nkpolar, nkz)); upa2_ICW_bin          = 0.d0
    allocate (bpa2_ICW_bin       (1:nkpolar, nkz)); bpa2_ICW_bin          = 0.d0

    call split_KAW_ICW( &
                        phi    , psi    , upa    , bpa    , &
                        phi_KAW, psi_KAW, upa_KAW, bpa_KAW, &
                        phi_ICW, psi_ICW, upa_ICW, bpa_ICW  &
                      )

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          upe2   (i, k, j) = 0.5d0*cabs2(phi(i, k, j))*kprp2(i, k, j)
          bpe2   (i, k, j) = 0.5d0*cabs2(psi(i, k, j))*kprp2(i, k, j)
          upa2   (i, k, j) = 0.5d0*cabs2(upa(i, k, j))
          bpa2   (i, k, j) = 0.5d0*cabs2(bpa(i, k, j))

          ux2    (i, k, j) = 0.5d0*cabs2(ky(j)*phi(i, k, j))
          uy2    (i, k, j) = 0.5d0*cabs2(kx(i)*phi(i, k, j))
          bx2    (i, k, j) = 0.5d0*cabs2(ky(j)*psi(i, k, j))
          by2    (i, k, j) = 0.5d0*cabs2(kx(i)*psi(i, k, j))

          upe2old(i, k, j) = 0.5d0*cabs2(phi_old(i, k, j))*kprp2(i, k, j)
          bpe2old(i, k, j) = 0.5d0*cabs2(psi_old(i, k, j))*kprp2(i, k, j)
          upa2old(i, k, j) = 0.5d0*cabs2(upa_old(i, k, j))
          bpa2old(i, k, j) = 0.5d0*cabs2(bpa_old(i, k, j))

          phi_mid  = 0.5d0*(phi(i, k, j) + phi_old(i, k, j))
          psi_mid  = 0.5d0*(psi(i, k, j) + psi_old(i, k, j))
          jpa_mid  = -kprp2(i, k, j)*psi_mid
          upa_mid  = 0.5d0*(upa(i, k, j) + upa_old(i, k, j))
          bpa_mid  = 0.5d0*(bpa(i, k, j) + bpa_old(i, k, j))
          fphi_mid = 0.5d0*(fphi(i, k, j) + fphi_old(i, k, j))
          fpsi_mid = 0.5d0*(fpsi(i, k, j) + fpsi_old(i, k, j))
          fupa_mid = 0.5d0*(fupa(i, k, j) + fupa_old(i, k, j))
          fbpa_mid = 0.5d0*(fbpa(i, k, j) + fbpa_old(i, k, j))

          upe2dissip_x(i, k, j) =  nupe_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*cabs2(phi_mid)*kprp2(i, k, j)
          upe2dissip_z(i, k, j) =  nupe_z*(kz2  (k)      /kz2_max  )** nupe_z_exp*cabs2(phi_mid)*kprp2(i, k, j)
          bpe2dissip_x(i, k, j) = etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*cabs2(psi_mid)*kprp2(i, k, j)
          bpe2dissip_z(i, k, j) = etape_z*(kz2  (k)      /kz2_max  )**etape_z_exp*cabs2(psi_mid)*kprp2(i, k, j)
          upa2dissip_x(i, k, j) =  nupa_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*cabs2(upa_mid)
          upa2dissip_z(i, k, j) =  nupa_z*(kz2(k)        /kz2_max  )** nupe_z_exp*cabs2(upa_mid)
          bpa2dissip_x(i, k, j) = etapa_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*cabs2(bpa_mid)
          bpa2dissip_z(i, k, j) = etapa_z*(kz2(k)        /kz2_max  )**etape_z_exp*cabs2(bpa_mid)
          p_aw        (i, k, j) = - 0.5d0*(  -kprp2(i, k, j)*fphi_mid*conjg(phi_mid) + conjg(-kprp2(i, k, j)*fphi_mid)*phi_mid &
                                           + fpsi_mid*conjg(jpa_mid) + conjg(fpsi_mid)*jpa_mid) 
          p_compr     (i, k, j) =   0.5d0*(  fupa_mid*conjg(upa_mid) + conjg(fupa_mid)*upa_mid &
                                           + fbpa_mid*conjg(bpa_mid) + conjg(fbpa_mid)*bpa_mid) 

          zppe_mid  =  phi_mid +  psi_mid
          zmpe_mid  =  phi_mid -  psi_mid
          fzppe_mid = fphi_mid + fpsi_mid
          fzmpe_mid = fphi_mid - fpsi_mid
          zppe2      (i, k, j) = cabs2(phi_mid + psi_mid)*kprp2(i, k, j)
          zmpe2      (i, k, j) = cabs2(phi_mid - psi_mid)*kprp2(i, k, j)
          zppa2      (i, k, j) = cabs2(upa_mid + bpa_mid)
          zmpa2      (i, k, j) = cabs2(upa_mid - bpa_mid)

          p_xhl       (i, k, j) = 0.25d0*( &
                                            zppe_mid*conjg(kprp2(i, k, j)*fzppe_mid) + conjg(zppe_mid)*kprp2(i, k, j)*fzppe_mid &
                                          - zmpe_mid*conjg(kprp2(i, k, j)*fzmpe_mid) - conjg(zmpe_mid)*kprp2(i, k, j)*fzmpe_mid &
                                        ) 

          upe2_KAW   (i, k, j) = 0.5d0*cabs2(phi_KAW(i, k, j))*kprp2(i, k, j)
          bpe2_KAW   (i, k, j) = 0.5d0*cabs2(psi_KAW(i, k, j))*kprp2(i, k, j)
          upa2_KAW   (i, k, j) = 0.5d0*cabs2(upa_KAW(i, k, j))
          bpa2_KAW   (i, k, j) = 0.5d0*cabs2(bpa_KAW(i, k, j))
          upe2_ICW   (i, k, j) = 0.5d0*cabs2(phi_ICW(i, k, j))*kprp2(i, k, j)
          bpe2_ICW   (i, k, j) = 0.5d0*cabs2(psi_ICW(i, k, j))*kprp2(i, k, j)
          upa2_ICW   (i, k, j) = 0.5d0*cabs2(upa_ICW(i, k, j))
          bpa2_ICW   (i, k, j) = 0.5d0*cabs2(bpa_ICW(i, k, j))

          upe2_KAW_dissip_x(i, k, j) =  nupe_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*cabs2(phi_KAW(i, k, j))*kprp2(i, k, j)
          upe2_KAW_dissip_z(i, k, j) =  nupe_z*(kz2  (k)      /kz2_max  )** nupe_z_exp*cabs2(phi_KAW(i, k, j))*kprp2(i, k, j)
          bpe2_KAW_dissip_x(i, k, j) = etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*cabs2(psi_KAW(i, k, j))*kprp2(i, k, j)
          bpe2_KAW_dissip_z(i, k, j) = etape_z*(kz2  (k)      /kz2_max  )**etape_z_exp*cabs2(psi_KAW(i, k, j))*kprp2(i, k, j)
          upa2_KAW_dissip_x(i, k, j) =  nupa_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*cabs2(upa_KAW(i, k, j))
          upa2_KAW_dissip_z(i, k, j) =  nupa_z*(kz2(k)        /kz2_max  )** nupe_z_exp*cabs2(upa_KAW(i, k, j))
          bpa2_KAW_dissip_x(i, k, j) = etapa_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*cabs2(bpa_KAW(i, k, j))
          bpa2_KAW_dissip_z(i, k, j) = etapa_z*(kz2(k)        /kz2_max  )**etape_z_exp*cabs2(bpa_KAW(i, k, j))

          upe2_ICW_dissip_x(i, k, j) =  nupe_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*cabs2(phi_ICW(i, k, j))*kprp2(i, k, j)
          upe2_ICW_dissip_z(i, k, j) =  nupe_z*(kz2  (k)      /kz2_max  )** nupe_z_exp*cabs2(phi_ICW(i, k, j))*kprp2(i, k, j)
          bpe2_ICW_dissip_x(i, k, j) = etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*cabs2(psi_ICW(i, k, j))*kprp2(i, k, j)
          bpe2_ICW_dissip_z(i, k, j) = etape_z*(kz2  (k)      /kz2_max  )**etape_z_exp*cabs2(psi_ICW(i, k, j))*kprp2(i, k, j)
          upa2_ICW_dissip_x(i, k, j) =  nupa_x*(kprp2(i, k, j)/kprp2_max)** nupe_x_exp*cabs2(upa_ICW(i, k, j))
          upa2_ICW_dissip_z(i, k, j) =  nupa_z*(kz2(k)        /kz2_max  )** nupe_z_exp*cabs2(upa_ICW(i, k, j))
          bpa2_ICW_dissip_x(i, k, j) = etapa_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp*cabs2(bpa_ICW(i, k, j))
          bpa2_ICW_dissip_z(i, k, j) = etapa_z*(kz2(k)        /kz2_max  )**etape_z_exp*cabs2(bpa_ICW(i, k, j))

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
            p_xhl       (i, k, j) = 2.0d0*p_xhl       (i, k, j)

            zppe2    (i, k, j) = 2.0d0*zppe2   (i, k, j)
            zmpe2    (i, k, j) = 2.0d0*zmpe2   (i, k, j)
            zppa2    (i, k, j) = 2.0d0*zppa2   (i, k, j)
            zmpa2    (i, k, j) = 2.0d0*zmpa2   (i, k, j)

            upe2_KAW (i, k, j) = 2.0d0*upe2_KAW(i, k, j)
            bpe2_KAW (i, k, j) = 2.0d0*bpe2_KAW(i, k, j)
            upa2_KAW (i, k, j) = 2.0d0*upa2_KAW(i, k, j)
            bpa2_KAW (i, k, j) = 2.0d0*bpa2_KAW(i, k, j)
            upe2_ICW (i, k, j) = 2.0d0*upe2_ICW(i, k, j)
            bpe2_ICW (i, k, j) = 2.0d0*bpe2_ICW(i, k, j)
            upa2_ICW (i, k, j) = 2.0d0*upa2_ICW(i, k, j)
            bpa2_ICW (i, k, j) = 2.0d0*bpa2_ICW(i, k, j)

            upe2_KAW_dissip_x(i, k, j) = 2.0d0*upe2_KAW_dissip_x(i, k, j)
            upe2_KAW_dissip_z(i, k, j) = 2.0d0*upe2_KAW_dissip_z(i, k, j)
            bpe2_KAW_dissip_x(i, k, j) = 2.0d0*bpe2_KAW_dissip_x(i, k, j)
            bpe2_KAW_dissip_z(i, k, j) = 2.0d0*bpe2_KAW_dissip_z(i, k, j)
            upa2_KAW_dissip_x(i, k, j) = 2.0d0*upa2_KAW_dissip_x(i, k, j)
            upa2_KAW_dissip_z(i, k, j) = 2.0d0*upa2_KAW_dissip_z(i, k, j)
            bpa2_KAW_dissip_x(i, k, j) = 2.0d0*bpa2_KAW_dissip_x(i, k, j)
            bpa2_KAW_dissip_z(i, k, j) = 2.0d0*bpa2_KAW_dissip_z(i, k, j)

            upe2_ICW_dissip_x(i, k, j) = 2.0d0*upe2_ICW_dissip_x(i, k, j)
            upe2_ICW_dissip_z(i, k, j) = 2.0d0*upe2_ICW_dissip_z(i, k, j)
            bpe2_ICW_dissip_x(i, k, j) = 2.0d0*bpe2_ICW_dissip_x(i, k, j)
            bpe2_ICW_dissip_z(i, k, j) = 2.0d0*bpe2_ICW_dissip_z(i, k, j)
            upa2_ICW_dissip_x(i, k, j) = 2.0d0*upa2_ICW_dissip_x(i, k, j)
            upa2_ICW_dissip_z(i, k, j) = 2.0d0*upa2_ICW_dissip_z(i, k, j)
            bpa2_ICW_dissip_x(i, k, j) = 2.0d0*bpa2_ICW_dissip_x(i, k, j)
            bpa2_ICW_dissip_z(i, k, j) = 2.0d0*bpa2_ICW_dissip_z(i, k, j)
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
    p_xhl_sum      = sum(p_xhl); call sum_reduce(p_xhl_sum, 0)

    zppe2_sum = sum(zppe2); call sum_reduce(zppe2_sum, 0)
    zmpe2_sum = sum(zmpe2); call sum_reduce(zmpe2_sum, 0)
    zppa2_sum = sum(zppa2); call sum_reduce(zppa2_sum, 0)
    zmpa2_sum = sum(zmpa2); call sum_reduce(zmpa2_sum, 0)

    upe2_KAW_dissip_sum = sum(upe2_KAW_dissip_x + upe2_KAW_dissip_z); call sum_reduce(upe2_KAW_dissip_sum, 0)
    bpe2_KAW_dissip_sum = sum(bpe2_KAW_dissip_x + bpe2_KAW_dissip_z); call sum_reduce(bpe2_KAW_dissip_sum, 0)
    upa2_KAW_dissip_sum = sum(upa2_KAW_dissip_x + upa2_KAW_dissip_z); call sum_reduce(upa2_KAW_dissip_sum, 0)
    bpa2_KAW_dissip_sum = sum(bpa2_KAW_dissip_x + bpa2_KAW_dissip_z); call sum_reduce(bpa2_KAW_dissip_sum, 0)

    upe2_ICW_dissip_sum = sum(upe2_ICW_dissip_x + upe2_ICW_dissip_z); call sum_reduce(upe2_ICW_dissip_sum, 0)
    bpe2_ICW_dissip_sum = sum(bpe2_ICW_dissip_x + bpe2_ICW_dissip_z); call sum_reduce(bpe2_ICW_dissip_sum, 0)
    upa2_ICW_dissip_sum = sum(upa2_ICW_dissip_x + upa2_ICW_dissip_z); call sum_reduce(upa2_ICW_dissip_sum, 0)
    bpa2_ICW_dissip_sum = sum(bpa2_ICW_dissip_x + bpa2_ICW_dissip_z); call sum_reduce(bpa2_ICW_dissip_sum, 0)
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
    call get_polar_spectrum_2d(zppe2, zppe2_bin)
    call get_polar_spectrum_2d(zmpe2, zmpe2_bin)
    call get_polar_spectrum_2d(zppa2, zppa2_bin)
    call get_polar_spectrum_2d(zmpa2, zmpa2_bin)
    call get_polar_spectrum_2d(upe2_KAW, upe2_KAW_bin)
    call get_polar_spectrum_2d(bpe2_KAW, bpe2_KAW_bin)
    call get_polar_spectrum_2d(upa2_KAW, upa2_KAW_bin)
    call get_polar_spectrum_2d(bpa2_KAW, bpa2_KAW_bin)
    call get_polar_spectrum_2d(upe2_ICW, upe2_ICW_bin)
    call get_polar_spectrum_2d(bpe2_ICW, bpe2_ICW_bin)
    call get_polar_spectrum_2d(upa2_ICW, upa2_ICW_bin)
    call get_polar_spectrum_2d(bpa2_ICW, bpa2_ICW_bin)
    call get_polar_spectrum_2d(upe2_KAW_dissip_x + upe2_KAW_dissip_z + &
                               bpe2_KAW_dissip_x + bpe2_KAW_dissip_z + &
                               upa2_KAW_dissip_x + upa2_KAW_dissip_z + &
                               bpa2_KAW_dissip_x + bpa2_KAW_dissip_z, dissip_KAW_bin)
    call get_polar_spectrum_2d(upe2_ICW_dissip_x + upe2_ICW_dissip_z + &
                               bpe2_ICW_dissip_x + bpe2_ICW_dissip_z + &
                               upa2_ICW_dissip_x + upa2_ICW_dissip_z + &
                               bpa2_ICW_dissip_x + bpa2_ICW_dissip_z, dissip_ICW_bin)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    ! call get_nonlinear_transfer(ntrans_aw_l_bin, ntrans_aw_g_bin, ntrans_compr_l_bin, ntrans_compr_g_bin)

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  upe2_sum, bpe2_sum, upa2_sum, bpa2_sum, &
                  upe2dot_sum, bpe2dot_sum, upa2dot_sum, bpa2dot_sum, &
                  upe2dissip_sum, bpe2dissip_sum, upa2dissip_sum, bpa2dissip_sum, &
                  p_aw_sum, p_compr_sum, p_xhl_sum, &
                  zppe2_sum, zmpe2_sum, zppa2_sum, zmpa2_sum, &
                  upe2_KAW_dissip_sum, bpe2_KAW_dissip_sum, upa2_KAW_dissip_sum, bpa2_KAW_dissip_sum, &
                  upe2_ICW_dissip_sum, bpe2_ICW_dissip_sum, upa2_ICW_dissip_sum, bpa2_ICW_dissip_sum, &
                  !
                  nkpolar, &
                  upe2_bin, bpe2_bin, upa2_bin, bpa2_bin, &
                  ux2_bin , uy2_bin , bx2_bin , by2_bin , &
                  zppe2_bin, zmpe2_bin, zppa2_bin, zmpa2_bin, &
                  upe2_KAW_bin, bpe2_KAW_bin, upa2_KAW_bin, bpa2_KAW_bin, &
                  upe2_ICW_bin, bpe2_ICW_bin, upa2_ICW_bin, bpa2_ICW_bin, &
                  p_aw_bin, p_compr_bin, &
                  dissip_aw_bin, dissip_compr_bin, &
                  dissip_KAW_bin, dissip_ICW_bin, &
                  ntrans_aw_l_bin, ntrans_compr_l_bin, &
                  ntrans_aw_g_bin, ntrans_compr_g_bin  &
                )

    deallocate(upe2)
    deallocate(bpe2)
    deallocate(upa2)
    deallocate(bpa2)
    deallocate(ux2)
    deallocate(bx2)
    deallocate(uy2)
    deallocate(by2)
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
    deallocate(p_xhl)
    deallocate(zppe2   )
    deallocate(zmpe2   )
    deallocate(zppa2   )
    deallocate(zmpa2   )
    deallocate(upe2_KAW)
    deallocate(bpe2_KAW)
    deallocate(upa2_KAW)
    deallocate(bpa2_KAW)
    deallocate(upe2_ICW)
    deallocate(bpe2_ICW)
    deallocate(upa2_ICW)
    deallocate(bpa2_ICW)
    deallocate(upe2_KAW_dissip_x)
    deallocate(upe2_KAW_dissip_z)
    deallocate(bpe2_KAW_dissip_x)
    deallocate(bpe2_KAW_dissip_z)
    deallocate(upa2_KAW_dissip_x)
    deallocate(upa2_KAW_dissip_z)
    deallocate(bpa2_KAW_dissip_x)
    deallocate(bpa2_KAW_dissip_z)
    deallocate(upe2_ICW_dissip_x)
    deallocate(upe2_ICW_dissip_z)
    deallocate(bpe2_ICW_dissip_x)
    deallocate(bpe2_ICW_dissip_z)
    deallocate(upa2_ICW_dissip_x)
    deallocate(upa2_ICW_dissip_z)
    deallocate(bpa2_ICW_dissip_x)
    deallocate(bpa2_ICW_dissip_z)

    deallocate(phi_KAW)
    deallocate(psi_KAW)
    deallocate(upa_KAW)
    deallocate(bpa_KAW)
    deallocate(phi_ICW)
    deallocate(psi_ICW)
    deallocate(upa_ICW)
    deallocate(bpa_ICW)

    deallocate (upe2_bin)
    deallocate (bpe2_bin)
    deallocate (upa2_bin)
    deallocate (bpa2_bin)
    deallocate (ux2_bin)
    deallocate (uy2_bin)
    deallocate (bx2_bin)
    deallocate (by2_bin)
    deallocate (p_aw_bin)
    deallocate (p_compr_bin)
    deallocate (ntrans_aw_l_bin)
    deallocate (ntrans_aw_g_bin)
    deallocate (ntrans_compr_l_bin)
    deallocate (ntrans_compr_g_bin)
    deallocate (dissip_aw_bin)
    deallocate (dissip_compr_bin)
    deallocate (dissip_KAW_bin)
    deallocate (dissip_ICW_bin)
    deallocate (zppe2_bin)
    deallocate (zmpe2_bin)
    deallocate (zppa2_bin)
    deallocate (zmpa2_bin)
    deallocate (upe2_KAW_bin)
    deallocate (bpe2_KAW_bin)
    deallocate (upa2_KAW_bin)
    deallocate (bpa2_KAW_bin)
    deallocate (upe2_ICW_bin)
    deallocate (bpe2_ICW_bin)
    deallocate (upa2_ICW_bin)
    deallocate (bpa2_ICW_bin)
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    28 Mar 2025
!! @brief   Split KAWs and ICWs
!-----------------------------------------------!
  subroutine split_KAW_ICW( &
                            phi    , psi    , upa    , bpa    , &
                            phi_KAW, psi_KAW, upa_KAW, bpa_KAW, &
                            phi_ICW, psi_ICW, upa_ICW, bpa_ICW  &
                          )
    use grid, only: nkx, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: kprp2, kprp2inv
    use params, only: zi, sgm
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in   ) :: phi    , psi    , upa    , bpa
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: phi_KAW, psi_KAW, upa_KAW, bpa_KAW, &
                                                             phi_ICW, psi_ICW, upa_ICW, bpa_ICW

    integer :: i, j, k, istir

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
           bpa_KAW(i, k, j) = kprp2(i, k, j)*(bpa(i, k, j) + phi(i, k, j))/(1.d0 + kprp2(i, k, j))
           psi_KAW(i, k, j) = (kprp2(i, k, j)*psi(i, k, j) + sgm*upa(i, k, j))/(sgm**2 + kprp2(i, k, j))
           phi_KAW(i, k, j) = kprp2inv(i, k, j)*bpa_KAW(i, k, j)
           upa_KAW(i, k, j) = sgm*psi_KAW(i, k, j)

           bpa_ICW(i, k, j) = bpa(i, k, j) - bpa_KAW(i, k, j)
           psi_ICW(i, k, j) = psi(i, k, j) - psi_KAW(i, k, j)
           phi_ICW(i, k, j) = phi(i, k, j) - phi_KAW(i, k, j)
           upa_ICW(i, k, j) = upa(i, k, j) - upa_KAW(i, k, j)
        end do
      end do
    end do
  end subroutine split_KAW_ICW


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


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate shell-to-shell transfer
!-----------------------------------------------!
  subroutine loop_diagnostics_nltrans
  ! under development...
  end subroutine loop_diagnostics_nltrans


end module diagnostics

