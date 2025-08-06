!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Diagnostics for MHD_INCOMP
!-----------------------------------------------!
module diagnostics
  use diagnostics_common
  use p3dfft
  implicit none

  public :: init_diagnostics, finish_diagnostics
  public :: loop_diagnostics, loop_diagnostics_2D, loop_diagnostics_kpar, loop_diagnostics_SF2
  public :: loop_diagnostics_nltrans

  private

  integer  :: nl, nsamp_r0, nsamp_ang
  real(r8), allocatable :: ll(:), lpar(:), lper(:)
  integer  :: unit_u2kxy_vs_kz, unit_u2kxz_vs_ky
  integer  :: unit_b2kxy_vs_kz, unit_b2kxz_vs_ky
contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of diagnostics
!-----------------------------------------------!
  subroutine init_diagnostics
    use params, only: inputfile
    use diagnostics_common, only: read_parameters
    use diagnostics_common, only: init_polar_spectrum_2d, init_polar_spectrum_3d
    use diagnostics_common, only: init_series_modes
    use io, only: init_io 
    use grid, only: nlz
    implicit none

    call read_parameters(inputfile)

    if(nlz == 2) then
      call init_polar_spectrum_2d
    else
      call init_polar_spectrum_3d
    endif
    call init_SF2
    call init_io(nkpolar, kpbin, nkpolar_log, kpbin_log, nl, lpar, lper)
    call init_series_modes
    call set_unit_for_polar_spectrum_3d_in_2d
  end subroutine init_diagnostics

  subroutine set_unit_for_polar_spectrum_3d_in_2d
    use file, only: open_output_file_binary
    use mp, only: proc0
    implicit none

    if(proc0) then
      call open_output_file_binary (unit_u2kxy_vs_kz, 'out2d/u2_kxy_vs_kz.dat')
      call open_output_file_binary (unit_u2kxz_vs_ky, 'out2d/u2_kxz_vs_ky.dat')
      call open_output_file_binary (unit_b2kxy_vs_kz, 'out2d/b2_kxy_vs_kz.dat')
      call open_output_file_binary (unit_b2kxz_vs_ky, 'out2d/b2_kxz_vs_ky.dat')
    endif
  end subroutine set_unit_for_polar_spectrum_3d_in_2d


!-----------------------------------------------!
!> @author  YK
!! @date    15 Apr 2020
!! @brief   Diagnostics in loop
!-----------------------------------------------!
  subroutine loop_diagnostics
    use io, only: loop_io
    use mp, only: proc0
    use grid, only: kx, ky, kz, k2_max
    use grid, only: nlx, nly, nlz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: dlx, dly, dlz
    use grid, only:  lx,  ly,  lz
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use fields, only: ux_old, uy_old, uz_old
    use fields, only: bx_old, by_old, bz_old
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, shear, q
    use force, only: fux, fuy, fuz, fux_old, fuy_old, fuz_old
    use force, only: fbx, fby, fbz, fbx_old, fby_old, fbz_old
    use shearing_box, only: shear_flg, tsc, nremap, k2t
    use utils, only: cabs2
    implicit none
    integer :: i, j, k

    real(r8), allocatable, dimension(:,:,:) :: u2, ux2, uy2, uz2
    real(r8), allocatable, dimension(:,:,:) :: b2, bx2, by2, bz2
    real(r8), allocatable, dimension(:,:,:) :: u2old, b2old
    real(r8), allocatable, dimension(:,:,:) :: u2dissip, b2dissip
    real(r8), allocatable, dimension(:,:,:) :: p_ext_ene, p_ext_xhl, p_re, p_ma
    real(r8), allocatable, dimension(:,:,:) :: zp2, zm2
    real(r8), allocatable, dimension(:,:,:) :: src

    real(r8) :: u2_sum, b2_sum
    real(r8) :: u2dot_sum, b2dot_sum
    real(r8) :: u2dissip_sum, b2dissip_sum
    real(r8) :: p_ext_ene_sum, p_ext_xhl_sum, p_re_sum, p_ma_sum
    real(r8) :: zp2_sum, zm2_sum
    real(r8) :: bx0, by0, bz0 ! mean magnetic field

    real(r8), dimension(:), allocatable :: u2_bin, ux2_bin, uy2_bin, uz2_bin
    real(r8), dimension(:), allocatable :: b2_bin, bx2_bin, by2_bin, bz2_bin
    real(r8), dimension(:), allocatable :: zp2_bin, zm2_bin
    real(r8), dimension(:), allocatable :: u2dissip_bin, b2dissip_bin
    real(r8), dimension(:), allocatable :: p_re_bin, p_ma_bin
    complex(r8) ::  ux_mid,  uy_mid,  uz_mid
    complex(r8) ::  bx_mid,  by_mid,  bz_mid
    complex(r8) :: fux_mid, fuy_mid, fuz_mid
    complex(r8) :: fbx_mid, fby_mid, fbz_mid

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(u2       , source=src)
    allocate(ux2      , source=src)
    allocate(uy2      , source=src)
    allocate(uz2      , source=src)
    allocate(b2       , source=src)
    allocate(bx2      , source=src)
    allocate(by2      , source=src)
    allocate(bz2      , source=src)
    allocate(u2old    , source=src)
    allocate(b2old    , source=src)
    allocate(u2dissip , source=src)
    allocate(b2dissip , source=src)
    allocate(p_ext_ene, source=src)
    allocate(p_ext_xhl, source=src)
    allocate(p_re     , source=src)
    allocate(p_ma     , source=src)
    allocate(zp2      , source=src)
    allocate(zm2      , source=src)
    deallocate(src)

    allocate(  u2_bin    (1:nkpolar), source=0.d0)
    allocate( ux2_bin    (1:nkpolar), source=0.d0)
    allocate( uy2_bin    (1:nkpolar), source=0.d0)
    allocate( uz2_bin    (1:nkpolar), source=0.d0)
    allocate(  b2_bin    (1:nkpolar), source=0.d0)
    allocate( bx2_bin    (1:nkpolar), source=0.d0)
    allocate( by2_bin    (1:nkpolar), source=0.d0)
    allocate( bz2_bin    (1:nkpolar), source=0.d0)
    allocate( zp2_bin    (1:nkpolar), source=0.d0)
    allocate( zm2_bin    (1:nkpolar), source=0.d0)
    allocate(u2dissip_bin(1:nkpolar), source=0.d0)
    allocate(b2dissip_bin(1:nkpolar), source=0.d0)
    allocate(p_re_bin    (1:nkpolar), source=0.d0)
    allocate(p_ma_bin    (1:nkpolar), source=0.d0)

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
            bx0 = sqrt(cabs2(bx(i,k,j)))
            by0 = sqrt(cabs2(by(i,k,j)))
            bz0 = sqrt(cabs2(bz(i,k,j)))
          endif
          u2   (i, k, j) = 0.5d0*(cabs2(ux    (i, k, j)) + cabs2(uy    (i, k, j)) + cabs2(uz    (i, k, j)))
          ux2  (i, k, j) = 0.5d0*(cabs2(ux    (i, k, j)))
          uy2  (i, k, j) = 0.5d0*(cabs2(uy    (i, k, j)))
          uz2  (i, k, j) = 0.5d0*(cabs2(uz    (i, k, j)))
          b2   (i, k, j) = 0.5d0*(cabs2(bx    (i, k, j)) + cabs2(by    (i, k, j)) + cabs2(bz    (i, k, j)))
          bx2  (i, k, j) = 0.5d0*(cabs2(bx    (i, k, j)))
          by2  (i, k, j) = 0.5d0*(cabs2(by    (i, k, j)))
          bz2  (i, k, j) = 0.5d0*(cabs2(bz    (i, k, j)))
          u2old(i, k, j) = 0.5d0*(cabs2(ux_old(i, k, j)) + cabs2(uy_old(i, k, j)) + cabs2(uz_old(i, k, j)))
          b2old(i, k, j) = 0.5d0*(cabs2(bx_old(i, k, j)) + cabs2(by_old(i, k, j)) + cabs2(bz_old(i, k, j)))
          zp2  (i, k, j) = cabs2(ux(i,k,j) + bx(i,k,j)) + cabs2(uy(i,k,j) + by(i,k,j)) + cabs2(uz(i,k,j) + bz(i,k,j))
          zm2  (i, k, j) = cabs2(ux(i,k,j) - bx(i,k,j)) + cabs2(uy(i,k,j) - by(i,k,j)) + cabs2(uz(i,k,j) - bz(i,k,j))

           ux_mid  = 0.5d0*( ux(i, k, j) +  ux_old(i, k, j))
           uy_mid  = 0.5d0*( uy(i, k, j) +  uy_old(i, k, j))
           uz_mid  = 0.5d0*( uz(i, k, j) +  uz_old(i, k, j))
           bx_mid  = 0.5d0*( bx(i, k, j) +  bx_old(i, k, j))
           by_mid  = 0.5d0*( by(i, k, j) +  by_old(i, k, j))
           bz_mid  = 0.5d0*( bz(i, k, j) +  bz_old(i, k, j))
          fux_mid  = 0.5d0*(fux(i, k, j) + fux_old(i, k, j))
          fuy_mid  = 0.5d0*(fuy(i, k, j) + fuy_old(i, k, j))
          fuz_mid  = 0.5d0*(fuz(i, k, j) + fuz_old(i, k, j))
          fbx_mid  = 0.5d0*(fbx(i, k, j) + fbx_old(i, k, j))
          fby_mid  = 0.5d0*(fby(i, k, j) + fby_old(i, k, j))
          fbz_mid  = 0.5d0*(fbz(i, k, j) + fbz_old(i, k, j))

          u2dissip (i, k, j) = (nu *(k2t(i, k, j)/k2_max) + nu_h *(k2t(i, k, j)/k2_max)**nu_h_exp ) &
                                *(cabs2(ux_mid) + cabs2(uy_mid) + cabs2(uz_mid))
          b2dissip (i, k, j) = (eta*(k2t(i, k, j)/k2_max) + eta_h*(k2t(i, k, j)/k2_max)**eta_h_exp) &
                                *(cabs2(bx_mid) + cabs2(by_mid) + cabs2(bz_mid))
          p_ext_ene(i, k, j) = 0.5d0*( &
                                  (fux_mid*conjg(ux_mid) + conjg(fux_mid)*ux_mid) &
                                + (fuy_mid*conjg(uy_mid) + conjg(fuy_mid)*uy_mid) &
                                + (fuz_mid*conjg(uz_mid) + conjg(fuz_mid)*uz_mid) &
                                + (fbx_mid*conjg(bx_mid) + conjg(fbx_mid)*bx_mid) &
                                + (fby_mid*conjg(by_mid) + conjg(fby_mid)*by_mid) &
                                + (fbz_mid*conjg(bz_mid) + conjg(fbz_mid)*bz_mid) &
                              )
          p_ext_xhl(i, k, j) = 0.5d0*( &
                                  (fux_mid*conjg(bx_mid) + conjg(fux_mid)*bx_mid) &
                                + (fuy_mid*conjg(by_mid) + conjg(fuy_mid)*by_mid) &
                                + (fuz_mid*conjg(bz_mid) + conjg(fuz_mid)*bz_mid) &
                                + (fbx_mid*conjg(ux_mid) + conjg(fbx_mid)*ux_mid) &
                                + (fby_mid*conjg(uy_mid) + conjg(fby_mid)*uy_mid) &
                                + (fbz_mid*conjg(uz_mid) + conjg(fbz_mid)*uz_mid) &
                              )
          p_re    (i, k, j) = + 0.5d0*q*shear_flg*(ux_mid*conjg(uy_mid) + conjg(ux_mid)*uy_mid)
          p_ma    (i, k, j) = - 0.5d0*q*shear_flg*(bx_mid*conjg(by_mid) + conjg(bx_mid)*by_mid)

          ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
          ! |f_{k_x, k_y}|^2 = \left( \sum_{k_y = -n_{k_y}/2}^{-1}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 0}
          ! \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 1}^{n_{k_y}/2}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} \right) |f_{k_x, k_y}|^2\]
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (j /= 1) then
            u2       (i, k, j) = 2.0d0*u2       (i, k, j)
            ux2      (i, k, j) = 2.0d0*ux2      (i, k, j)
            uy2      (i, k, j) = 2.0d0*uy2      (i, k, j)
            uz2      (i, k, j) = 2.0d0*uz2      (i, k, j)
            b2       (i, k, j) = 2.0d0*b2       (i, k, j)
            bx2      (i, k, j) = 2.0d0*bx2      (i, k, j)
            by2      (i, k, j) = 2.0d0*by2      (i, k, j)
            bz2      (i, k, j) = 2.0d0*bz2      (i, k, j)
            u2old    (i, k, j) = 2.0d0*u2old    (i, k, j)
            b2old    (i, k, j) = 2.0d0*b2old    (i, k, j)
            u2dissip (i, k, j) = 2.0d0*u2dissip (i, k, j)
            b2dissip (i, k, j) = 2.0d0*b2dissip (i, k, j)
            p_ext_ene(i, k, j) = 2.0d0*p_ext_ene(i, k, j)
            p_ext_xhl(i, k, j) = 2.0d0*p_ext_xhl(i, k, j)
            p_re     (i, k, j) = 2.0d0*p_re      (i, k, j)
            p_ma     (i, k, j) = 2.0d0*p_ma      (i, k, j)
            zp2      (i, k, j) = 2.0d0*zp2       (i, k, j)
            zm2      (i, k, j) = 2.0d0*zm2       (i, k, j)
          endif

        end do
      end do
    end do

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    u2_sum        = sum(u2); call sum_reduce(u2_sum, 0)
    b2_sum        = sum(b2); call sum_reduce(b2_sum, 0)
                  
    u2dot_sum     = sum((u2 - u2old)/dt); call sum_reduce(u2dot_sum, 0)
    b2dot_sum     = sum((b2 - b2old)/dt); call sum_reduce(b2dot_sum, 0)
                  
    u2dissip_sum  = sum(u2dissip); call sum_reduce(u2dissip_sum, 0)
    b2dissip_sum  = sum(b2dissip); call sum_reduce(b2dissip_sum, 0)
                  
    p_ext_ene_sum = sum(p_ext_ene); call sum_reduce(p_ext_ene_sum, 0)
    p_ext_xhl_sum = sum(p_ext_xhl); call sum_reduce(p_ext_xhl_sum, 0)
    p_re_sum      = sum(p_re ); call sum_reduce(p_re_sum , 0)
    p_ma_sum      = sum(p_ma ); call sum_reduce(p_ma_sum , 0)
                  
    zp2_sum       = sum(zp2); call sum_reduce(zp2_sum, 0)
    zm2_sum       = sum(zm2); call sum_reduce(zm2_sum, 0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv       bin over (kx, ky, kz)       vvvvvvvvvvvvvvvvvv!
    call get_polar_spectrum_3d(  u2    ,   u2_bin    )
    call get_polar_spectrum_3d( ux2    ,  ux2_bin    )
    call get_polar_spectrum_3d( uy2    ,  uy2_bin    )
    call get_polar_spectrum_3d( uz2    ,  uz2_bin    )
    call get_polar_spectrum_3d(  b2    ,   b2_bin    )
    call get_polar_spectrum_3d( bx2    ,  bx2_bin    )
    call get_polar_spectrum_3d( by2    ,  by2_bin    )
    call get_polar_spectrum_3d( bz2    ,  bz2_bin    )
    call get_polar_spectrum_3d( zp2    ,  zp2_bin    )
    call get_polar_spectrum_3d( zm2    ,  zm2_bin    )
    call get_polar_spectrum_3d(u2dissip, u2dissip_bin)
    call get_polar_spectrum_3d(b2dissip, b2dissip_bin)
    call get_polar_spectrum_3d(p_re    , p_re_bin    )
    call get_polar_spectrum_3d(p_ma    , p_ma_bin    )
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  u2_sum, b2_sum, &
                  u2dot_sum, b2dot_sum, &
                  u2dissip_sum, b2dissip_sum, &
                  p_ext_ene_sum, p_ext_xhl_sum, &
                  p_re_sum, p_ma_sum, &
                  zp2_sum, zm2_sum, &
                  bx0, by0, bz0, &
                  !
                  nkpolar, &
                  u2_bin, ux2_bin, uy2_bin, uz2_bin, &
                  b2_bin, bx2_bin, by2_bin, bz2_bin, &
                  zp2_bin, zm2_bin, &
                  u2dissip_bin, b2dissip_bin, &
                  p_re_bin, p_ma_bin &
                )

    deallocate(u2)
    deallocate(ux2)
    deallocate(uy2)
    deallocate(uz2)
    deallocate(b2)
    deallocate(bx2)
    deallocate(by2)
    deallocate(bz2)
    deallocate(u2old)
    deallocate(b2old)
    deallocate(u2dissip)
    deallocate(b2dissip)
    deallocate(p_ext_ene)
    deallocate(p_ext_xhl)
    deallocate(p_re)
    deallocate(p_ma)
    deallocate(zp2)
    deallocate(zm2)

    deallocate(  u2_bin    )
    deallocate( ux2_bin    )
    deallocate( uy2_bin    )
    deallocate( uz2_bin    )
    deallocate(  b2_bin    )
    deallocate( bx2_bin    )
    deallocate( by2_bin    )
    deallocate( bz2_bin    )
    deallocate( zp2_bin    )
    deallocate( zm2_bin    )
    deallocate(u2dissip_bin)
    deallocate(b2dissip_bin)
    deallocate(p_re_bin    )
    deallocate(p_ma_bin    )
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Diagnostics for cross section of fileds
!-----------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io_2D
    use utils, only: curl
    use mp, only: proc0
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: shear
    use shearing_box, only: to_non_shearing_coordinate, tsc, nremap
    use utils, only: cabs2
    implicit none

    complex(r8), allocatable, dimension(:,:,:) :: f
    real(r8)   , allocatable, dimension(:,:,:) :: fr
    real(r8)   , allocatable, dimension(:,:)   :: ux_r_z0, ux_r_x0, ux_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: uy_r_z0, uy_r_x0, uy_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: uz_r_z0, uz_r_x0, uz_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: wx_r_z0, wx_r_x0, wx_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: wy_r_z0, wy_r_x0, wy_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: wz_r_z0, wz_r_x0, wz_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: bx_r_z0, bx_r_x0, bx_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: by_r_z0, by_r_x0, by_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: bz_r_z0, bz_r_x0, bz_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: jx_r_z0, jx_r_x0, jx_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: jy_r_z0, jy_r_x0, jy_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: jz_r_z0, jz_r_x0, jz_r_y0
    complex(r8), allocatable, dimension(:,:,:) :: wx, wy, wz, jx, jy, jz

    real(r8)   , allocatable, dimension(:,:,:) :: u2, b2
    real(r8)   , allocatable, dimension(:,:)   :: u2_kxy, u2_kyz, u2_kxz
    real(r8)   , allocatable, dimension(:,:)   :: b2_kxy, b2_kyz, b2_kxz

    real(r8)   , allocatable, dimension(:,:)   :: src1, src2, src3
    complex(r8), allocatable, dimension(:,:,:) :: src4

    integer :: i, j, k

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(f (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f   = 0.d0
    allocate(fr(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); fr  = 0.d0

    allocate(src1(ilx_st:ilx_en, ily_st:ily_en), source=0.d0) 
    allocate(src2(ily_st:ily_en, ilz_st:ilz_en), source=0.d0) 
    allocate(src3(ilx_st:ilx_en, ilz_st:ilz_en), source=0.d0)
    allocate(ux_r_z0, source=src1)
    allocate(ux_r_x0, source=src2)
    allocate(ux_r_y0, source=src3)

    allocate(uy_r_z0, source=src1)
    allocate(uy_r_x0, source=src2)
    allocate(uy_r_y0, source=src3)

    allocate(uz_r_z0, source=src1)
    allocate(uz_r_x0, source=src2)
    allocate(uz_r_y0, source=src3)

    allocate(wx_r_z0, source=src1)
    allocate(wx_r_x0, source=src2)
    allocate(wx_r_y0, source=src3)

    allocate(wy_r_z0, source=src1)
    allocate(wy_r_x0, source=src2)
    allocate(wy_r_y0, source=src3)

    allocate(wz_r_z0, source=src1)
    allocate(wz_r_x0, source=src2)
    allocate(wz_r_y0, source=src3)

    allocate(bx_r_z0, source=src1)
    allocate(bx_r_x0, source=src2)
    allocate(bx_r_y0, source=src3)

    allocate(by_r_z0, source=src1)
    allocate(by_r_x0, source=src2)
    allocate(by_r_y0, source=src3)

    allocate(bz_r_z0, source=src1)
    allocate(bz_r_x0, source=src2)
    allocate(bz_r_y0, source=src3)

    allocate(jx_r_z0, source=src1)
    allocate(jx_r_x0, source=src2)
    allocate(jx_r_y0, source=src3)

    allocate(jy_r_z0, source=src1)
    allocate(jy_r_x0, source=src2)
    allocate(jy_r_y0, source=src3)

    allocate(jz_r_z0, source=src1)
    allocate(jz_r_x0, source=src2)
    allocate(jz_r_y0, source=src3)
    deallocate(src1, src2, src3)

    allocate(src4(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    allocate(wx, source=src4)
    allocate(wy, source=src4)
    allocate(wz, source=src4)
    allocate(jx, source=src4)
    allocate(jy, source=src4)
    allocate(jz, source=src4)
    deallocate(src4)

    allocate(u2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(b2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)

    allocate(src1(ikx_st:ikx_en, iky_st:iky_en), source=0.d0); 
    allocate(src2(iky_st:iky_en, ikz_st:ikz_en), source=0.d0); 
    allocate(src3(ikx_st:ikx_en, ikz_st:ikz_en), source=0.d0)
    allocate(u2_kxy, source=src1)
    allocate(u2_kyz, source=src2)
    allocate(u2_kxz, source=src3)
    allocate(b2_kxy, source=src1)
    allocate(b2_kyz, source=src2)
    allocate(b2_kxz, source=src3)
    deallocate(src1, src2, src3)

    !vvvvvvvvvvvvvvvvvv         2D cut of fields          vvvvvvvvvvvvvvvvvv!
    call curl(ux, uy, uz, wx, wy, wz)
    call curl(bx, by, bz, jx, jy, jz)

    f = ux; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, ux_r_z0, ux_r_x0, ux_r_y0)
    f = uy; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, uy_r_z0, uy_r_x0, uy_r_y0)
    f = uz; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, uz_r_z0, uz_r_x0, uz_r_y0)
    f = bx; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, bx_r_z0, bx_r_x0, bx_r_y0)
    f = by; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, by_r_z0, by_r_x0, by_r_y0)
    f = bz; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, bz_r_z0, bz_r_x0, bz_r_y0)
    f = wx; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, wx_r_z0, wx_r_x0, wx_r_y0)
    f = wy; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, wy_r_z0, wy_r_x0, wy_r_y0)
    f = wz; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, wz_r_z0, wz_r_x0, wz_r_y0)
    f = jx; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, jx_r_z0, jx_r_x0, jx_r_y0)
    f = jy; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, jy_r_z0, jy_r_x0, jy_r_y0)
    f = jz; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, jz_r_z0, jz_r_x0, jz_r_y0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv      2D spectra (avg over 1D)     vvvvvvvvvvvvvvvvvv!
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          u2(i, k, j) = 0.5d0*(cabs2(ux(i, k, j)) + cabs2(uy(i, k, j)) + cabs2(uz(i, k, j)))
          b2(i, k, j) = 0.5d0*(cabs2(bx(i, k, j)) + cabs2(by(i, k, j)) + cabs2(bz(i, k, j)))
          ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
          !-----------------------------------------------------------------------------------------------------------------------------------
          ! The volume integral of a quadratic function is
          ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
          ! |f_{k_x, k_y}|^2 = \left( \sum_{k_y = -n_{k_y}/2}^{-1}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 0}
          ! \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 1}^{n_{k_y}/2}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} \right) |f_{k_x, k_y}|^2\]
          ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
          !-----------------------------------------------------------------------------------------------------------------------------------
          if (j /= 1) then
            u2(i, k, j) = 2.0d0*u2(i, k, j)
            b2(i, k, j) = 2.0d0*b2(i, k, j)
          endif
        end do
      end do
    end do
    call sum_2d_k(u2, u2_kxy, u2_kyz, u2_kxz)
    call sum_2d_k(b2, b2_kxy, b2_kyz, b2_kxz)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv       bin over (kx, ky) vs kz     vvvvvvvvvvvvvvvvvv!
    call write_polar_spectrum_2d_in_3d(u2, 'z', unit_u2kxy_vs_kz)
    call write_polar_spectrum_2d_in_3d(b2, 'z', unit_b2kxy_vs_kz)
    !vvvvvvvvvvvvvvvvvv       bin over (kx, kz) vs ky     vvvvvvvvvvvvvvvvvv!
    call write_polar_spectrum_2d_in_3d(u2, 'y', unit_u2kxz_vs_ky)
    call write_polar_spectrum_2d_in_3d(b2, 'y', unit_b2kxz_vs_ky)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io_2D( &
                  ux_r_z0, ux_r_x0, ux_r_y0, &
                  uy_r_z0, uy_r_x0, uy_r_y0, &
                  uz_r_z0, uz_r_x0, uz_r_y0, &
                  !
                  wx_r_z0, wx_r_x0, wx_r_y0, &
                  wy_r_z0, wy_r_x0, wy_r_y0, &
                  wz_r_z0, wz_r_x0, wz_r_y0, &
                  !
                  bx_r_z0, bx_r_x0, bx_r_y0, &
                  by_r_z0, by_r_x0, by_r_y0, &
                  bz_r_z0, bz_r_x0, bz_r_y0, &
                  !
                  jx_r_z0, jx_r_x0, jx_r_y0, &
                  jy_r_z0, jy_r_x0, jy_r_y0, &
                  jz_r_z0, jz_r_x0, jz_r_y0, &
                  !
                  u2_kxy, u2_kyz, u2_kxz, &
                  b2_kxy, b2_kyz, b2_kxz  &
                )

    deallocate(f)
    deallocate(fr)

    deallocate(ux_r_z0)
    deallocate(ux_r_x0)
    deallocate(ux_r_y0)

    deallocate(uy_r_z0)
    deallocate(uy_r_x0)
    deallocate(uy_r_y0)

    deallocate(uz_r_z0)
    deallocate(uz_r_x0)
    deallocate(uz_r_y0)

    deallocate(wx_r_z0)
    deallocate(wx_r_x0)
    deallocate(wx_r_y0)

    deallocate(wy_r_z0)
    deallocate(wy_r_x0)
    deallocate(wy_r_y0)

    deallocate(wz_r_z0)
    deallocate(wz_r_x0)
    deallocate(wz_r_y0)

    deallocate(bx_r_z0)
    deallocate(bx_r_x0)
    deallocate(bx_r_y0)

    deallocate(by_r_z0)
    deallocate(by_r_x0)
    deallocate(by_r_y0)

    deallocate(bz_r_z0)
    deallocate(bz_r_x0)
    deallocate(bz_r_y0)

    deallocate(jx_r_z0)
    deallocate(jx_r_x0)
    deallocate(jx_r_y0)

    deallocate(jy_r_z0)
    deallocate(jy_r_x0)
    deallocate(jy_r_y0)

    deallocate(jz_r_z0)
    deallocate(jz_r_x0)
    deallocate(jz_r_y0)

    deallocate(wx)
    deallocate(wy)
    deallocate(wz)
    deallocate(jx)
    deallocate(jy)
    deallocate(jz)

    deallocate(u2)
    deallocate(b2)

    deallocate(u2_kxy)
    deallocate(u2_kyz)
    deallocate(u2_kxz)
    deallocate(b2_kxy)
    deallocate(b2_kyz)
    deallocate(b2_kxz)
  end subroutine loop_diagnostics_2D


!-----------------------------------------------!
!> @author  YK
!! @date    4 Jul 2021
!! @brief   Initialize order structure function
!-----------------------------------------------!
  subroutine init_SF2
    use grid, only: lx, ly, lz, dlx, dly, dlz
    use mp, only: nproc
    implicit none
    real(r8) :: l, d
    integer  :: il

    l = 0.5d0*min(lx, ly, lz)
    d = max(dlx, dly, dlz)

    nl = int(l/d)

    allocate(ll(nl))

    do il = 1, nl
      ! ll(il) = 10.d0**((dlog10(l) - dlog10(d))/(nl - 1)*(il - 1) + dlog10(d))
      ll(il) = (l - d)/(nl - 1)*(il - 1) + d
    enddo

    allocate(lpar(nl), source=ll)
    allocate(lper(nl), source=ll)

    nsamp_r0  = int(sqrt(real(SF2_nsample)/nproc))
    nsamp_ang = int(sqrt(real(SF2_nsample)/nproc))

  end subroutine init_SF2


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Second order structure function
!-----------------------------------------------!
  subroutine loop_diagnostics_SF2
    use io, only: loop_io_SF2
    use mp, only: proc0, proc_id, nproc
    use mp, only: sum_allreduce
    use grid, only: xx, yy, zz, nlx, nly, nlz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use fields, only: bx, by, bz, ux, uy, uz
    use params, only: pi
    use utils, only: ranf
    use time, only: microsleep
    use time_stamp, only: put_time_stamp, timer_diagnostics_total, timer_diagnostics_SF2
    implicit none
    include 'mpif.h'
    real(r8) :: ll_x, ll_y, ll_z, theta, phi
    real(r8) :: x0, y0, z0, x1, y1, z1
    integer  :: i0, j0, k0, i1, j1, k1
    integer  :: il, isamp_r0, isamp_ang , iproc
    integer  :: demand(nproc, 3)
    real(r8) :: supply_b(nproc, 3), supply_u(nproc, 3)
    real(r8) :: b0x, b0y, b0z, b1x, b1y, b1z
    real(r8) :: u0x, u0y, u0z, u1x, u1y, u1z
    real(r8) :: bloc_x, bloc_y, bloc_z, lpar_, lper_
    real(r8) :: sf2b(nl, nl), sf2u(nl, nl), count(nl, nl)
    integer  :: ilpar, ilper
    complex(r8), allocatable, dimension(:,:,:) :: f
    real(r8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r, ux_r, uy_r, uz_r
    real(r8), allocatable, dimension(:,:,:) :: src

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_SF2)

    allocate(f   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f    = 0.d0

    allocate(src(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(bx_r, source=src)
    allocate(by_r, source=src)
    allocate(bz_r, source=src)
    allocate(ux_r, source=src)
    allocate(uy_r, source=src)
    allocate(uz_r, source=src)
    deallocate(src)

    f = bx ; call p3dfft_btran_c2r(f, bx_r, 'tff')
    f = by ; call p3dfft_btran_c2r(f, by_r, 'tff')
    f = bz ; call p3dfft_btran_c2r(f, bz_r, 'tff')
    f = ux ; call p3dfft_btran_c2r(f, ux_r, 'tff')
    f = uy ; call p3dfft_btran_c2r(f, uy_r, 'tff')
    f = uz ; call p3dfft_btran_c2r(f, uz_r, 'tff')

    sf2b (:, :) = 0.d0
    sf2u (:, :) = 0.d0
    count(:, :) = 0
    do il = 1, nl
      ! if(proc0) print *, il, nl
      do isamp_r0 = 1, nsamp_r0
        call microsleep(1000)

        ! Pick a random initial point
        i0 = ilx_st + floor((ilx_en + 1 - ilx_st)*ranf()); x0 = xx(i0)
        j0 = ily_st + floor((ily_en + 1 - ily_st)*ranf()); y0 = yy(j0)
        k0 = ilz_st + floor((ilz_en + 1 - ilz_st)*ranf()); z0 = zz(k0)

        b0x = bx_r(j0, k0, i0)
        b0y = by_r(j0, k0, i0)
        b0z = bz_r(j0, k0, i0)
        u0x = ux_r(j0, k0, i0)
        u0y = uy_r(j0, k0, i0)
        u0z = uz_r(j0, k0, i0)

        do isamp_ang = 1, nsamp_ang
          theta = pi*(ranf() - 0.5d0)
          phi   = 2.d0*pi*ranf() 

          ll_x = ll(il)*dsin(theta)*dcos(phi)
          ll_y = ll(il)*dsin(theta)*dsin(phi)
          ll_z = ll(il)*dcos(theta)

          x1 = x0 + ll_x
          y1 = y0 + ll_y
          z1 = z0 + ll_z

          ! When the separation vector crosses the boundary, flip the direction of the separation
          ! so that the separation stays inside the domain.
          ! if(x1 > maxval(xx) .or. x1 < minval(xx)) phi   = pi - phi
          ! if(y1 > maxval(yy) .or. y1 < minval(yy)) phi   =    - phi
          ! if(z1 > maxval(zz) .or. z1 < minval(zz)) theta = pi - theta
          ! x1 = x0 + ll(il)*dsin(theta)*dcos(phi)
          ! y1 = y0 + ll(il)*dsin(theta)*dsin(phi)
          ! z1 = z0 + ll(il)*dcos(theta)
          ! Ignore if the separation goes beyond the boundary
          ! if(x1 > maxval(xx) .or. x1 < minval(xx) .or. &
             ! y1 > maxval(yy) .or. y1 < minval(yy) .or. &
             ! z1 > maxval(zz) .or. z1 < minval(zz)      &
            ! ) then

            ! x1 = x0
            ! y1 = y0
            ! z1 = z0
          ! endif

          i1 = minloc(abs(xx - x1), 1); x1 = xx(i1)
          j1 = minloc(abs(yy - y1), 1); y1 = yy(j1)
          k1 = minloc(abs(zz - z1), 1); z1 = zz(k1)

          ll_x = (x1 - x0)
          ll_y = (y1 - y0)
          ll_z = (z1 - z0)

          ! Periodic boundary
          if(x1 > maxval(xx)) x1 = minval(xx) + x1 - maxval(xx)
          if(y1 > maxval(yy)) y1 = minval(yy) + y1 - maxval(yy) 
          if(z1 > maxval(zz)) z1 = minval(zz) + z1 - maxval(zz)
          if(x1 < minval(xx)) x1 = maxval(xx) + x1 - minval(xx)
          if(y1 < minval(yy)) y1 = maxval(yy) + y1 - minval(yy) 
          if(z1 < minval(zz)) z1 = maxval(zz) + z1 - minval(zz)
          i1 = minloc(abs(xx - x1), 1)
          j1 = minloc(abs(yy - y1), 1)
          k1 = minloc(abs(zz - z1), 1)

          ! Set demand array; "proc_id" process needs a value at (i1, j1, k1) = demand(proc_id+1, :) 
          demand(:, :) = 0
          demand(proc_id+1, :) = [i1, j1, k1]

          call sum_allreduce(demand)

          ! Get B and u vector at the demanded location
          supply_b(:, :) = 0
          supply_u(:, :) = 0
          do iproc = 1, nproc
            i1 = demand(iproc, 1)
            j1 = demand(iproc, 2)
            k1 = demand(iproc, 3)

            if(i1 >= ilx_st .and. i1 <= ilx_en .and. &
               j1 >= ily_st .and. j1 <= ily_en .and. &
               k1 >= ilz_st .and. k1 <= ilz_en) then

               supply_b(iproc, 1) = bx_r(j1, k1, i1)
               supply_b(iproc, 2) = by_r(j1, k1, i1)
               supply_b(iproc, 3) = bz_r(j1, k1, i1)
               supply_u(iproc, 1) = ux_r(j1, k1, i1)
               supply_u(iproc, 2) = uy_r(j1, k1, i1)
               supply_u(iproc, 3) = uz_r(j1, k1, i1)
            endif
          enddo

          call sum_allreduce(supply_b)
          call sum_allreduce(supply_u)

          ! Get B vector at the demanding location
          b1x = supply_b(proc_id+1, 1)
          b1y = supply_b(proc_id+1, 2)
          b1z = supply_b(proc_id+1, 3)
          u1x = supply_u(proc_id+1, 1)
          u1y = supply_u(proc_id+1, 2)
          u1z = supply_u(proc_id+1, 3)

          ! Local mean magnetic field
          bloc_x = (b0x + b1x)/2.d0
          bloc_y = (b0y + b1y)/2.d0
          bloc_z = (b0z + b1z)/2.d0

          ! Parallel and perpendicular components of the separation vector
          if(ll_x**2 + ll_y**2 + ll_z**2 > 0.d0) then
            lpar_ = dabs((ll_x*bloc_x + ll_y*bloc_y + ll_z*bloc_z))/dsqrt((bloc_x**2 + bloc_y**2 + bloc_z**2))
            lper_ = dsqrt((   (bloc_y*ll_z - bloc_z*ll_y)**2 &
                            + (bloc_z*ll_x - bloc_x*ll_z)**2 &
                            + (bloc_x*ll_y - bloc_y*ll_x)**2)/(bloc_x**2 + bloc_y**2 + bloc_z**2))
            ilpar = minloc(abs(lpar_ - lpar), 1)
            ilper = minloc(abs(lper_ - lper), 1)

            sf2b(ilpar, ilper) = sf2b(ilpar, ilper) &
                                + (b1x - b0x)**2 + (b1y - b0y)**2 + (b1z - b0z)**2
            sf2u(ilpar, ilper) = sf2u(ilpar, ilper) &
                                + (u1x - u0x)**2 + (u1y - u0y)**2 + (u1z - u0z)**2
            count(ilpar, ilper) = count(ilpar, ilper) + 1
          endif
        enddo
      enddo
    enddo

    call sum_allreduce(sf2b)
    call sum_allreduce(sf2u)
    call sum_allreduce(count)

    ! Average the structure function
    do ilper = 1, nl
      do ilpar = 1, nl
        if(count(ilpar, ilper) == 0) then
          sf2b(ilpar, ilper) = 0.d0
          sf2u(ilpar, ilper) = 0.d0
        else
          sf2b(ilpar, ilper) = sf2b(ilpar, ilper)/count(ilpar, ilper)
          sf2u(ilpar, ilper) = sf2u(ilpar, ilper)/count(ilpar, ilper)
        endif
      enddo
    enddo

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_SF2)

    call loop_io_SF2(nl, sf2b, sf2u)

    deallocate(f)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
  end subroutine loop_diagnostics_SF2


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate kpar(k) & delta b/b0
!           k_\|(k) = \left(\frac{\langle|\mathbf{b}_{0,k} \cdot\nabla \delta\mathbf{b}_k|^2\rangle}
!           {\langle b_{0,k}^2\rangle\langle \delta b_k^2\rangle}\right)^{1/2}  \\ 
!           \mathbf{b}_{0,k}(\mathbf{x}) = \calF^{-1}\sum_{|\bm{k}|' \le k/2} \mathbf{b}_{\mathbf{k}'} \\  
!           \delta\mathbf{b}_{k}(\mathbf{x}) = \calF^{-1}\sum_{k/2 \le |\bm{k}|' \le 2k} \mathbf{b}_{\mathbf{k}'}
!-----------------------------------------------!
  subroutine loop_diagnostics_kpar
    use fields, only: bx, by, bz, ux, uy, uz
    use mp, only: proc0, sum_allreduce
    use grid, only: ky, kz, nlx, nly, nlz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: zi
    use shearing_box, only: k2t, kxt
    use io, only: loop_io_kpar
    use time_stamp, only: put_time_stamp, timer_diagnostics_total, timer_diagnostics_kpar
    use utils, only: cabs2
    implicit none
    integer :: ii, i, j, k

    real   (r8), dimension(:), allocatable :: kpar_b, kpar_u, b1_ovr_b0

    complex(r8), allocatable, dimension(:,:,:) :: bx0, by0, bz0                ! local mean field
    complex(r8), allocatable, dimension(:,:,:) :: bx1, by1, bz1                ! local fluctuating field
    complex(r8), allocatable, dimension(:,:,:) :: ux1, uy1, uz1                ! local fluctuating field
    complex(r8), allocatable, dimension(:,:,:) :: dbx1_dx, dby1_dx, dbz1_dx    
    complex(r8), allocatable, dimension(:,:,:) :: dbx1_dy, dby1_dy, dbz1_dy    
    complex(r8), allocatable, dimension(:,:,:) :: dbx1_dz, dby1_dz, dbz1_dz    
    complex(r8), allocatable, dimension(:,:,:) :: dux1_dx, duy1_dx, duz1_dx    
    complex(r8), allocatable, dimension(:,:,:) :: dux1_dy, duy1_dy, duz1_dy    
    complex(r8), allocatable, dimension(:,:,:) :: dux1_dz, duy1_dz, duz1_dz    
    real   (r8), allocatable, dimension(:,:,:) :: bx0r, by0r, bz0r             
    real   (r8), allocatable, dimension(:,:,:) :: bx1r, by1r, bz1r             
    real   (r8), allocatable, dimension(:,:,:) :: ux1r, uy1r, uz1r             
    real   (r8), allocatable, dimension(:,:,:) :: dbx1_dxr, dby1_dxr, dbz1_dxr 
    real   (r8), allocatable, dimension(:,:,:) :: dbx1_dyr, dby1_dyr, dbz1_dyr 
    real   (r8), allocatable, dimension(:,:,:) :: dbx1_dzr, dby1_dzr, dbz1_dzr 
    real   (r8), allocatable, dimension(:,:,:) :: dux1_dxr, duy1_dxr, duz1_dxr 
    real   (r8), allocatable, dimension(:,:,:) :: dux1_dyr, duy1_dyr, duz1_dyr 
    real   (r8), allocatable, dimension(:,:,:) :: dux1_dzr, duy1_dzr, duz1_dzr 

    real   (r8), allocatable, dimension(:,:,:) :: b0_gradb1_sq, b0_gradu1_sq, b0sq, b1sq, u1sq
    real   (r8) :: b0_gradb1_sq_avg, b0_gradu1_sq_avg, b0sq_avg, b1sq_avg, u1sq_avg

    real   (r8), allocatable, dimension(:,:,:) :: bx0hat, by0hat, bz0hat ! local mean field unit vector
    real   (r8), allocatable, dimension(:,:,:) :: b1par, b1prpx, b1prpy  ! b1par : projection of b1 to b0 => Pseudo AW
                                                                         ! b1per : b1 - b1par*b0hat       => Shear AW
    real   (r8), allocatable, dimension(:,:,:) :: u1par, u1prpx, u1prpy  ! u1par : projection of u1 to b0 => Pseudo AW
                                                                         ! u1per : u1 - u1par*b0hat       => Shear AW
    real   (r8), allocatable, dimension(:) :: b1par2, b1prp2
    real   (r8), allocatable, dimension(:) :: u1par2, u1prp2

    complex(r8), allocatable, dimension(:,:,:) :: src_c
    real   (r8), allocatable, dimension(:,:,:) :: src_r

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_kpar)

    allocate(kpar_b   (nkpolar_log), source=0.d0)
    allocate(kpar_u   (nkpolar_log), source=0.d0)
    allocate(b1_ovr_b0(nkpolar_log), source=0.d0)
    allocate(b1par2   (nkpolar_log), source=0.d0)
    allocate(b1prp2   (nkpolar_log), source=0.d0)
    allocate(u1par2   (nkpolar_log), source=0.d0)
    allocate(u1prp2   (nkpolar_log), source=0.d0)

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    allocate(bx0    , source=src_c)
    allocate(by0    , source=src_c)
    allocate(bz0    , source=src_c)
    allocate(bx1    , source=src_c)
    allocate(by1    , source=src_c)
    allocate(bz1    , source=src_c)
    allocate(ux1    , source=src_c)
    allocate(uy1    , source=src_c)
    allocate(uz1    , source=src_c)
    allocate(dbx1_dx, source=src_c)
    allocate(dby1_dx, source=src_c)
    allocate(dbz1_dx, source=src_c)
    allocate(dbx1_dy, source=src_c)
    allocate(dby1_dy, source=src_c)
    allocate(dbz1_dy, source=src_c)
    allocate(dbx1_dz, source=src_c)
    allocate(dby1_dz, source=src_c)
    allocate(dbz1_dz, source=src_c)
    allocate(dux1_dx, source=src_c)
    allocate(duy1_dx, source=src_c)
    allocate(duz1_dx, source=src_c)
    allocate(dux1_dy, source=src_c)
    allocate(duy1_dy, source=src_c)
    allocate(duz1_dy, source=src_c)
    allocate(dux1_dz, source=src_c)
    allocate(duy1_dz, source=src_c)
    allocate(duz1_dz, source=src_c)
    deallocate(src_c)

    allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(bx0r    , source=src_r)
    allocate(by0r    , source=src_r)
    allocate(bz0r    , source=src_r)
    allocate(bx1r    , source=src_r)
    allocate(by1r    , source=src_r)
    allocate(bz1r    , source=src_r)
    allocate(ux1r    , source=src_r)
    allocate(uy1r    , source=src_r)
    allocate(uz1r    , source=src_r)
    allocate(dbx1_dxr, source=src_r)
    allocate(dby1_dxr, source=src_r)
    allocate(dbz1_dxr, source=src_r)
    allocate(dbx1_dyr, source=src_r)
    allocate(dby1_dyr, source=src_r)
    allocate(dbz1_dyr, source=src_r)
    allocate(dbx1_dzr, source=src_r)
    allocate(dby1_dzr, source=src_r)
    allocate(dbz1_dzr, source=src_r)
    allocate(dux1_dxr, source=src_r)
    allocate(duy1_dxr, source=src_r)
    allocate(duz1_dxr, source=src_r)
    allocate(dux1_dyr, source=src_r)
    allocate(duy1_dyr, source=src_r)
    allocate(duz1_dyr, source=src_r)
    allocate(dux1_dzr, source=src_r)
    allocate(duy1_dzr, source=src_r)
    allocate(duz1_dzr, source=src_r)

    allocate(b0_gradb1_sq, source=src_r)
    allocate(b0_gradu1_sq, source=src_r)
    allocate(b0sq        , source=src_r)
    allocate(b1sq        , source=src_r)
    allocate(u1sq        , source=src_r)

    allocate(bx0hat      , source=src_r)
    allocate(by0hat      , source=src_r)
    allocate(bz0hat      , source=src_r)
    allocate(b1par       , source=src_r)
    allocate(b1prpx      , source=src_r)
    allocate(b1prpy      , source=src_r)
    allocate(u1par       , source=src_r)
    allocate(u1prpx      , source=src_r)
    allocate(u1prpy      , source=src_r)
    deallocate(src_r)

    ! get kpar for each kprp_log(ii)
    do ii = 1, nkpolar_log

      ! filter out to get local mean field and local fluctuating field
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            ! smaller than kprp/2
            if(k2t(i, k, j) < (0.5d0*kpbin_log(ii))**2) then
              bx0(i, k, j) = bx(i, k, j)
              by0(i, k, j) = by(i, k, j)
              bz0(i, k, j) = bz(i, k, j)
            else
              bx0(i, k, j) = 0.d0
              by0(i, k, j) = 0.d0
              bz0(i, k, j) = 0.d0
            endif

            ! larger than kprp/2 and smaller than 2*kprp
            if(k2t(i, k, j) >= (0.5d0*kpbin_log(ii))**2 .and. k2t(i, k, j) < (2.0d0*kpbin_log(ii))**2) then
              bx1(i, k, j) = bx(i, k, j)
              by1(i, k, j) = by(i, k, j)
              bz1(i, k, j) = bz(i, k, j)
              ux1(i, k, j) = ux(i, k, j)
              uy1(i, k, j) = uy(i, k, j)
              uz1(i, k, j) = uz(i, k, j)
            else
              bx1(i, k, j) = 0.d0
              by1(i, k, j) = 0.d0
              bz1(i, k, j) = 0.d0
              ux1(i, k, j) = 0.d0
              uy1(i, k, j) = 0.d0
              uz1(i, k, j) = 0.d0
            endif

            if(cabs2(bx0(i,k,j)) < epsilon(1.d0) .or. cabs2(bx0(i,k,j)) > 1.d0/epsilon(1.d0)) bx0(i,k,j) = 0.d0
            if(cabs2(by0(i,k,j)) < epsilon(1.d0) .or. cabs2(by0(i,k,j)) > 1.d0/epsilon(1.d0)) by0(i,k,j) = 0.d0
            if(cabs2(bz0(i,k,j)) < epsilon(1.d0) .or. cabs2(bz0(i,k,j)) > 1.d0/epsilon(1.d0)) bz0(i,k,j) = 0.d0
            if(cabs2(bx1(i,k,j)) < epsilon(1.d0) .or. cabs2(bx1(i,k,j)) > 1.d0/epsilon(1.d0)) bx1(i,k,j) = 0.d0
            if(cabs2(by1(i,k,j)) < epsilon(1.d0) .or. cabs2(by1(i,k,j)) > 1.d0/epsilon(1.d0)) by1(i,k,j) = 0.d0
            if(cabs2(bz1(i,k,j)) < epsilon(1.d0) .or. cabs2(bz1(i,k,j)) > 1.d0/epsilon(1.d0)) bz1(i,k,j) = 0.d0
            if(cabs2(ux1(i,k,j)) < epsilon(1.d0) .or. cabs2(ux1(i,k,j)) > 1.d0/epsilon(1.d0)) ux1(i,k,j) = 0.d0
            if(cabs2(uy1(i,k,j)) < epsilon(1.d0) .or. cabs2(uy1(i,k,j)) > 1.d0/epsilon(1.d0)) uy1(i,k,j) = 0.d0
            if(cabs2(uz1(i,k,j)) < epsilon(1.d0) .or. cabs2(uz1(i,k,j)) > 1.d0/epsilon(1.d0)) uz1(i,k,j) = 0.d0

            dbx1_dx(i, k, j) = zi*kxt(i,j)*bx1(i, k, j)
            dby1_dx(i, k, j) = zi*kxt(i,j)*by1(i, k, j)
            dbz1_dx(i, k, j) = zi*kxt(i,j)*bz1(i, k, j)

            dbx1_dy(i, k, j) = zi*ky(j)   *bx1(i, k, j)
            dby1_dy(i, k, j) = zi*ky(j)   *by1(i, k, j)
            dbz1_dy(i, k, j) = zi*ky(j)   *bz1(i, k, j)
                                          
            dbx1_dz(i, k, j) = zi*kz(k)   *bx1(i, k, j)
            dby1_dz(i, k, j) = zi*kz(k)   *by1(i, k, j)
            dbz1_dz(i, k, j) = zi*kz(k)   *bz1(i, k, j)
                                          
            dux1_dx(i, k, j) = zi*kxt(i,j)*ux1(i, k, j)
            duy1_dx(i, k, j) = zi*kxt(i,j)*uy1(i, k, j)
            duz1_dx(i, k, j) = zi*kxt(i,j)*uz1(i, k, j)
                                          
            dux1_dy(i, k, j) = zi*ky(j)   *ux1(i, k, j)
            duy1_dy(i, k, j) = zi*ky(j)   *uy1(i, k, j)
            duz1_dy(i, k, j) = zi*ky(j)   *uz1(i, k, j)
                                          
            dux1_dz(i, k, j) = zi*kz(k)   *ux1(i, k, j)
            duy1_dz(i, k, j) = zi*kz(k)   *uy1(i, k, j)
            duz1_dz(i, k, j) = zi*kz(k)   *uz1(i, k, j)
          end do
        end do
      end do

      call p3dfft_btran_c2r(bx0, bx0r, 'tff')
      call p3dfft_btran_c2r(by0, by0r, 'tff')
      call p3dfft_btran_c2r(bz0, bz0r, 'tff')
      call p3dfft_btran_c2r(bx1, bx1r, 'tff')
      call p3dfft_btran_c2r(by1, by1r, 'tff')
      call p3dfft_btran_c2r(bz1, bz1r, 'tff')
      call p3dfft_btran_c2r(ux1, ux1r, 'tff')
      call p3dfft_btran_c2r(uy1, uy1r, 'tff')
      call p3dfft_btran_c2r(uz1, uz1r, 'tff')

      call p3dfft_btran_c2r(dbx1_dx, dbx1_dxr, 'tff')
      call p3dfft_btran_c2r(dby1_dx, dby1_dxr, 'tff')
      call p3dfft_btran_c2r(dbz1_dx, dbz1_dxr, 'tff')
      call p3dfft_btran_c2r(dbx1_dy, dbx1_dyr, 'tff')
      call p3dfft_btran_c2r(dby1_dy, dby1_dyr, 'tff')
      call p3dfft_btran_c2r(dbz1_dy, dbz1_dyr, 'tff')
      call p3dfft_btran_c2r(dbx1_dz, dbx1_dzr, 'tff')
      call p3dfft_btran_c2r(dby1_dz, dby1_dzr, 'tff')
      call p3dfft_btran_c2r(dbz1_dz, dbz1_dzr, 'tff')

      call p3dfft_btran_c2r(dux1_dx, dux1_dxr, 'tff')
      call p3dfft_btran_c2r(duy1_dx, duy1_dxr, 'tff')
      call p3dfft_btran_c2r(duz1_dx, duz1_dxr, 'tff')
      call p3dfft_btran_c2r(dux1_dy, dux1_dyr, 'tff')
      call p3dfft_btran_c2r(duy1_dy, duy1_dyr, 'tff')
      call p3dfft_btran_c2r(duz1_dy, duz1_dyr, 'tff')
      call p3dfft_btran_c2r(dux1_dz, dux1_dzr, 'tff')
      call p3dfft_btran_c2r(duy1_dz, duy1_dzr, 'tff')
      call p3dfft_btran_c2r(duz1_dz, duz1_dzr, 'tff')

      ! Get kpar & delta b/b0
      b0_gradb1_sq =   (bx0r*dbx1_dxr + by0r*dbx1_dyr + bz0r*dbx1_dzr)**2 &
                     + (bx0r*dby1_dxr + by0r*dby1_dyr + bz0r*dby1_dzr)**2 &
                     + (bx0r*dbz1_dxr + by0r*dbz1_dyr + bz0r*dbz1_dzr)**2 
      b0_gradu1_sq =   (bx0r*dux1_dxr + by0r*dux1_dyr + bz0r*dux1_dzr)**2 &
                     + (bx0r*duy1_dxr + by0r*duy1_dyr + bz0r*duy1_dzr)**2 &
                     + (bx0r*duz1_dxr + by0r*duz1_dyr + bz0r*duz1_dzr)**2 
      b0sq = bx0r**2 + by0r**2 + bz0r**2
      b1sq = bx1r**2 + by1r**2 + bz1r**2
      u1sq = ux1r**2 + uy1r**2 + uz1r**2

      b0_gradb1_sq_avg = sum(b0_gradb1_sq); call sum_allreduce(b0_gradb1_sq_avg); b0_gradb1_sq_avg = b0_gradb1_sq_avg/nlx/nly/nlz
      b0_gradu1_sq_avg = sum(b0_gradu1_sq); call sum_allreduce(b0_gradu1_sq_avg); b0_gradu1_sq_avg = b0_gradu1_sq_avg/nlx/nly/nlz
      b0sq_avg         = sum(b0sq)        ; call sum_allreduce(b0sq_avg        ); b0sq_avg         = b0sq_avg        /nlx/nly/nlz
      b1sq_avg         = sum(b1sq)        ; call sum_allreduce(b1sq_avg        ); b1sq_avg         = b1sq_avg        /nlx/nly/nlz
      u1sq_avg         = sum(u1sq)        ; call sum_allreduce(u1sq_avg        ); u1sq_avg         = u1sq_avg        /nlx/nly/nlz

      if (b0sq_avg /= 0.d0 .and. b1sq_avg /= 0.d0 .and. u1sq_avg /= 0.d0) then
        kpar_b   (ii) = dsqrt( b0_gradb1_sq_avg/(b1sq_avg*b0sq_avg) )
        kpar_u   (ii) = dsqrt( b0_gradu1_sq_avg/(u1sq_avg*b0sq_avg) )
        b1_ovr_b0(ii) = dsqrt( b1sq_avg/b0sq_avg )
      else
        kpar_b   (ii) = 0.d0
        kpar_u   (ii) = 0.d0
        b1_ovr_b0(ii) = 0.d0
      endif

      ! Get Shear AWs and pseudo AWs
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            if(bx0r(j,k,i)**2 + by0r(j,k,i)**2 + bz0r(j,k,i)**2 /= 0.d0) then 
              bx0hat(j,k,i) = bx0r(j,k,i)/dsqrt(bx0r(j,k,i)**2 + by0r(j,k,i)**2 + bz0r(j,k,i)**2)
              by0hat(j,k,i) = by0r(j,k,i)/dsqrt(bx0r(j,k,i)**2 + by0r(j,k,i)**2 + bz0r(j,k,i)**2)
              bz0hat(j,k,i) = bz0r(j,k,i)/dsqrt(bx0r(j,k,i)**2 + by0r(j,k,i)**2 + bz0r(j,k,i)**2)
            else 
              bx0hat(j,k,i) = 0.d0 
              by0hat(j,k,i) = 0.d0 
              bz0hat(j,k,i) = 0.d0 
            endif

            b1par (j,k,i) = bx1r(j,k,i)*bx0hat(j,k,i) + by1r(j,k,i)*by0hat(j,k,i) + bz1r(j,k,i)*bz0hat(j,k,i)
            b1prpx(j,k,i) = bx1r(j,k,i) - b1par(j,k,i)*bx0hat(j,k,i)
            b1prpy(j,k,i) = by1r(j,k,i) - b1par(j,k,i)*by0hat(j,k,i)

            u1par (j,k,i) = ux1r(j,k,i)*bx0hat(j,k,i) + uy1r(j,k,i)*by0hat(j,k,i) + uz1r(j,k,i)*bz0hat(j,k,i)
            u1prpx(j,k,i) = ux1r(j,k,i) - u1par(j,k,i)*bx0hat(j,k,i)
            u1prpy(j,k,i) = uy1r(j,k,i) - u1par(j,k,i)*by0hat(j,k,i)
          enddo
        enddo
      enddo

      b1par2(ii) = sum(b1par**2)             ; call sum_allreduce(b1par2(ii)); b1par2(ii) = b1par2(ii)/nlx/nly/nlz
      b1prp2(ii) = sum(b1prpx**2 + b1prpy**2); call sum_allreduce(b1prp2(ii)); b1prp2(ii) = b1prp2(ii)/nlx/nly/nlz
      u1par2(ii) = sum(u1par**2)             ; call sum_allreduce(u1par2(ii)); u1par2(ii) = u1par2(ii)/nlx/nly/nlz
      u1prp2(ii) = sum(u1prpx**2 + u1prpy**2); call sum_allreduce(u1prp2(ii)); u1prp2(ii) = u1prp2(ii)/nlx/nly/nlz
    enddo

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_kpar)

    call loop_io_kpar(nkpolar_log, kpar_b, kpar_u, b1_ovr_b0, &
                      b1par2, b1prp2, u1par2, u1prp2)

    deallocate(kpar_b)
    deallocate(kpar_u)
    deallocate(b1_ovr_b0)
    deallocate(bx0, by0, bz0)
    deallocate(bx1, by1, bz1)
    deallocate(ux1, uy1, uz1)
    deallocate(dbx1_dx, dby1_dx, dbz1_dx)
    deallocate(dbx1_dy, dby1_dy, dbz1_dy)
    deallocate(dbx1_dz, dby1_dz, dbz1_dz)
    deallocate(dux1_dx, duy1_dx, duz1_dx)
    deallocate(dux1_dy, duy1_dy, duz1_dy)
    deallocate(dux1_dz, duy1_dz, duz1_dz)
    deallocate(bx0r, by0r, bz0r)
    deallocate(bx1r, by1r, bz1r)
    deallocate(ux1r, uy1r, uz1r)
    deallocate(dbx1_dxr, dby1_dxr, dbz1_dxr)
    deallocate(dbx1_dyr, dby1_dyr, dbz1_dyr)
    deallocate(dbx1_dzr, dby1_dzr, dbz1_dzr)
    deallocate(dux1_dxr, duy1_dxr, duz1_dxr)
    deallocate(dux1_dyr, duy1_dyr, duz1_dyr)
    deallocate(dux1_dzr, duy1_dzr, duz1_dzr)
    deallocate(b0_gradb1_sq, b0_gradu1_sq, b0sq, b1sq, u1sq)
    deallocate(bx0hat)
    deallocate(by0hat)
    deallocate(bz0hat)
    deallocate(b1par )
    deallocate(b1prpx)
    deallocate(b1prpy)
    deallocate(u1par )
    deallocate(u1prpx)
    deallocate(u1prpy)
    deallocate(b1par2)
    deallocate(b1prp2)
    deallocate(u1par2)
    deallocate(u1prp2)

  end subroutine loop_diagnostics_kpar


!-----------------------------------------------!
!> @author  YK
!! @date    26 Jul 2019
!! @brief   Calculate shell-to-shell transfer
!!          trans_uu(K, Q) : u(Q) -> u(K)
!!          trans_bb(K, Q) : b(Q) -> b(K)
!!          trans_ub(K, Q) : u(Q) -> b(K)
!!          trans_bu(K, Q) : b(Q) -> u(K)
!-----------------------------------------------!
  subroutine loop_diagnostics_nltrans
    use fields, only: nfields
    use fields, only: iux, iuy, iuz
    use fields, only: ibx, iby, ibz
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use params, only: zi
    use grid, only: ky, kz
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use mp, only: proc0, sum_reduce
    use shearing_box, only: k2t, kxt
    use io, only: loop_io_nltrans
    use time_stamp, only: put_time_stamp, timer_diagnostics_total, timer_diagnostics_nltrans
    implicit none


    ! Forward FFT variables
    integer :: nftran = 36
    integer :: iflx_uu_xx = 1 , iflx_uu_xy = 2 , iflx_uu_xz = 3  !  
    integer :: iflx_uu_yx = 4 , iflx_uu_yy = 5 , iflx_uu_yz = 6  ! uu^Q (u^Q is filtered on |k| = Q) 
    integer :: iflx_uu_zx = 7 , iflx_uu_zy = 8 , iflx_uu_zz = 9  !  

    integer :: iflx_bb_xx = 10, iflx_bb_xy = 11, iflx_bb_xz = 12 !  
    integer :: iflx_bb_yx = 13, iflx_bb_yy = 14, iflx_bb_yz = 15 ! bb^Q (b^Q is filtered on |k| = Q) 
    integer :: iflx_bb_zx = 16, iflx_bb_zy = 17, iflx_bb_zz = 18 !  

    integer :: iflx_ub_xx = 19, iflx_ub_xy = 20, iflx_ub_xz = 21 !  
    integer :: iflx_ub_yx = 22, iflx_ub_yy = 23, iflx_ub_yz = 24 ! ub^Q (b^Q is filtered on |k| = Q) 
    integer :: iflx_ub_zx = 25, iflx_ub_zy = 26, iflx_ub_zz = 27 !  

    integer :: iflx_bu_xx = 28, iflx_bu_xy = 29, iflx_bu_xz = 30 !  
    integer :: iflx_bu_yx = 31, iflx_bu_yy = 32, iflx_bu_yz = 33 ! bu^Q (u^Q is filtered on |k| = Q) 
    integer :: iflx_bu_zx = 34, iflx_bu_zy = 35, iflx_bu_zz = 36 !  

    integer :: nnonlin = 12
    integer :: inonlin_uu_x = 1 , inonlin_uu_y = 2 , inonlin_uu_z = 3 
    integer :: inonlin_bb_x = 4 , inonlin_bb_y = 5 , inonlin_bb_z = 6 
    integer :: inonlin_ub_x = 7 , inonlin_ub_y = 8 , inonlin_ub_z = 9 
    integer :: inonlin_bu_x = 10, inonlin_bu_y = 11, inonlin_bu_z = 12


    integer  :: i, j, k, ii, jj
    real(r8) :: filter

    real   (r8), dimension(:,:), allocatable :: trans_uu, trans_bb, trans_ub, trans_bu
    complex(r8), allocatable, dimension(:,:,:,:) :: w, wtmp 
    complex(r8), allocatable, dimension(:,:,:,:) :: w_filtd
    complex(r8), allocatable, dimension(:,:,:,:) :: flx
    complex(r8), allocatable, dimension(:,:,:,:) :: nonlin

    real   (r8), allocatable, dimension(:,:,:,:) :: w_r
    real   (r8), allocatable, dimension(:,:,:,:) :: w_filtd_r
    real   (r8), allocatable, dimension(:,:,:,:) :: flx_r

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_nltrans)

    allocate(trans_uu (nkpolar_log, nkpolar_log), source=0.d0)
    allocate(trans_bb (nkpolar_log, nkpolar_log), source=0.d0)
    allocate(trans_ub (nkpolar_log, nkpolar_log), source=0.d0)
    allocate(trans_bu (nkpolar_log, nkpolar_log), source=0.d0)

    if(.not. allocated(w   )) allocate(w   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0,0.d0))
    if(.not. allocated(wtmp)) allocate(wtmp(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0,0.d0))
    if(.not. allocated(w_r )) allocate(w_r (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nfields), source=0.d0)

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          w(i,k,j,iux) = ux(i,k,j)
          w(i,k,j,iuy) = uy(i,k,j)
          w(i,k,j,iuz) = uz(i,k,j)
          w(i,k,j,ibx) = bx(i,k,j)
          w(i,k,j,iby) = by(i,k,j)
          w(i,k,j,ibz) = bz(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! get unfiltered fields in real space
    ! for some reason c2r_many doesn't work at Fugaku
    !call p3dfft_btran_c2r_many(w, nk_local_tot, w_r, nl_local_tot, nfields, 'tff')
    wtmp = w
    do i = 1, nfields
      call p3dfft_btran_c2r(w(:,:,:,i), w_r(:,:,:,i), 'tff')
    enddo
    w = wtmp
    if(allocated(wtmp)) deallocate(wtmp)


    ! get nonlinear transfer for each kprp_log(ii)
    do ii = 1, nkpolar_log

      if(.not. allocated(w_filtd)) allocate(w_filtd(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0,0.d0))
      ! filter out where kprp != kprp_log(ii)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            if(ii == nkpolar_log) then
              if(k2t(i, k, j) >= (kpbin_log(ii))**2) then
                filter = 1.d0
              else
                filter = 0.d0
              endif
            else
              if(k2t(i, k, j) >= (kpbin_log(ii))**2 .and. k2t(i, k, j) < (kpbin_log(ii + 1))**2) then
                filter = 1.d0
              else
                filter = 0.d0
              endif
            endif

            w_filtd(i, k, j, :) = filter*w(i, k, j, :)

          end do
        end do
      end do

      if(.not. allocated(w_filtd_r)) allocate(w_filtd_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nfields), source=0.d0)
      ! get filtered fields in real space
      ! call p3dfft_btran_c2r_many(w_filtd, nk_local_tot, w_filtd_r, nl_local_tot, nfields, 'tff')
      ! for some reason c2r_many doesn't work at Fugaku
      !call p3dfft_btran_c2r_many(w_filtd, nk_local_tot, w_filtd_r, nl_local_tot, nfields, 'tff')
      do i = 1, nfields
        call p3dfft_btran_c2r(w_filtd(:,:,:,i), w_filtd_r(:,:,:,i), 'tff')
      enddo
      if(allocated(w_filtd)) deallocate(w_filtd)

      if(.not. allocated(flx_r)) allocate(flx_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nftran ), source=0.d0)
      !$omp parallel do private(j, k) schedule(static)
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            ! uu^Q (u^Q is filtered on |k| = Q) 
            flx_r(j,k,i,iflx_uu_xx) = w_r(j,k,i,iux)*w_filtd_r(j,k,i,iux)
            flx_r(j,k,i,iflx_uu_xy) = w_r(j,k,i,iux)*w_filtd_r(j,k,i,iuy)
            flx_r(j,k,i,iflx_uu_xz) = w_r(j,k,i,iux)*w_filtd_r(j,k,i,iuz)

            flx_r(j,k,i,iflx_uu_yx) = w_r(j,k,i,iuy)*w_filtd_r(j,k,i,iux)
            flx_r(j,k,i,iflx_uu_yy) = w_r(j,k,i,iuy)*w_filtd_r(j,k,i,iuy)
            flx_r(j,k,i,iflx_uu_yz) = w_r(j,k,i,iuy)*w_filtd_r(j,k,i,iuz)

            flx_r(j,k,i,iflx_uu_zx) = w_r(j,k,i,iuz)*w_filtd_r(j,k,i,iux)
            flx_r(j,k,i,iflx_uu_zy) = w_r(j,k,i,iuz)*w_filtd_r(j,k,i,iuy)
            flx_r(j,k,i,iflx_uu_zz) = w_r(j,k,i,iuz)*w_filtd_r(j,k,i,iuz)

            ! bb^Q (b^Q is filtered on |k| = Q) 
            flx_r(j,k,i,iflx_bb_xx) = w_r(j,k,i,ibx)*w_filtd_r(j,k,i,ibx)
            flx_r(j,k,i,iflx_bb_xy) = w_r(j,k,i,ibx)*w_filtd_r(j,k,i,iby)
            flx_r(j,k,i,iflx_bb_xz) = w_r(j,k,i,ibx)*w_filtd_r(j,k,i,ibz)

            flx_r(j,k,i,iflx_bb_yx) = w_r(j,k,i,iby)*w_filtd_r(j,k,i,ibx)
            flx_r(j,k,i,iflx_bb_yy) = w_r(j,k,i,iby)*w_filtd_r(j,k,i,iby)
            flx_r(j,k,i,iflx_bb_yz) = w_r(j,k,i,iby)*w_filtd_r(j,k,i,ibz)

            flx_r(j,k,i,iflx_bb_zx) = w_r(j,k,i,ibz)*w_filtd_r(j,k,i,ibx)
            flx_r(j,k,i,iflx_bb_zy) = w_r(j,k,i,ibz)*w_filtd_r(j,k,i,iby)
            flx_r(j,k,i,iflx_bb_zz) = w_r(j,k,i,ibz)*w_filtd_r(j,k,i,ibz)

            ! ub^Q (b^Q is filtered on |k| = Q) 
            flx_r(j,k,i,iflx_ub_xx) = w_r(j,k,i,iux)*w_filtd_r(j,k,i,ibx)
            flx_r(j,k,i,iflx_ub_xy) = w_r(j,k,i,iux)*w_filtd_r(j,k,i,iby)
            flx_r(j,k,i,iflx_ub_xz) = w_r(j,k,i,iux)*w_filtd_r(j,k,i,ibz)

            flx_r(j,k,i,iflx_ub_yx) = w_r(j,k,i,iuy)*w_filtd_r(j,k,i,ibx)
            flx_r(j,k,i,iflx_ub_yy) = w_r(j,k,i,iuy)*w_filtd_r(j,k,i,iby)
            flx_r(j,k,i,iflx_ub_yz) = w_r(j,k,i,iuy)*w_filtd_r(j,k,i,ibz)

            flx_r(j,k,i,iflx_ub_zx) = w_r(j,k,i,iuz)*w_filtd_r(j,k,i,ibx)
            flx_r(j,k,i,iflx_ub_zy) = w_r(j,k,i,iuz)*w_filtd_r(j,k,i,iby)
            flx_r(j,k,i,iflx_ub_zz) = w_r(j,k,i,iuz)*w_filtd_r(j,k,i,ibz)

            ! bu^Q (u^Q is filtered on |k| = Q) 
            flx_r(j,k,i,iflx_bu_xx) = w_r(j,k,i,ibx)*w_filtd_r(j,k,i,iux)
            flx_r(j,k,i,iflx_bu_xy) = w_r(j,k,i,ibx)*w_filtd_r(j,k,i,iuy)
            flx_r(j,k,i,iflx_bu_xz) = w_r(j,k,i,ibx)*w_filtd_r(j,k,i,iuz)

            flx_r(j,k,i,iflx_bu_yx) = w_r(j,k,i,iby)*w_filtd_r(j,k,i,iux)
            flx_r(j,k,i,iflx_bu_yy) = w_r(j,k,i,iby)*w_filtd_r(j,k,i,iuy)
            flx_r(j,k,i,iflx_bu_yz) = w_r(j,k,i,iby)*w_filtd_r(j,k,i,iuz)

            flx_r(j,k,i,iflx_bu_zx) = w_r(j,k,i,ibz)*w_filtd_r(j,k,i,iux)
            flx_r(j,k,i,iflx_bu_zy) = w_r(j,k,i,ibz)*w_filtd_r(j,k,i,iuy)
            flx_r(j,k,i,iflx_bu_zz) = w_r(j,k,i,ibz)*w_filtd_r(j,k,i,iuz)

          enddo
        enddo
      enddo
      !$omp end parallel do
      if(allocated(w_filtd_r)) deallocate(w_filtd_r)

      if(.not. allocated(flx)) allocate(flx(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran ), source=(0.d0,0.d0))
      ! for some reason r2c_many doesn't work at Fugaku
      !call p3dfft_ftran_r2c_many(flx_r, nl_local_tot, flx, nk_local_tot, nftran, 'fft')
      do i = 1, nftran
        call p3dfft_ftran_r2c(flx_r(:,:,:,i), flx(:,:,:,i), 'fft')
      enddo
      if(allocated(flx_r)) deallocate(flx_r)

      !$omp workshare
      flx = flx/nlx/nly/nlz
      !$omp end workshare

      if(.not. allocated(nonlin)) allocate(nonlin(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nnonlin), source=(0.d0,0.d0))
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            nonlin(i,k,j,inonlin_uu_x) = zi*(  kxt(i,j)*flx(i,k,j,iflx_uu_xx) &
                                             + ky (j)  *flx(i,k,j,iflx_uu_yx) &
                                             + kz (k)  *flx(i,k,j,iflx_uu_zx) ) 
            nonlin(i,k,j,inonlin_uu_y) = zi*(  kxt(i,j)*flx(i,k,j,iflx_uu_xy) &
                                             + ky (j)  *flx(i,k,j,iflx_uu_yy) &
                                             + kz (k)  *flx(i,k,j,iflx_uu_zy) ) 
            nonlin(i,k,j,inonlin_uu_z) = zi*(  kxt(i,j)*flx(i,k,j,iflx_uu_xz) &
                                             + ky (j)  *flx(i,k,j,iflx_uu_yz) &
                                             + kz (k)  *flx(i,k,j,iflx_uu_zz) ) 

            nonlin(i,k,j,inonlin_bb_x) = zi*(  kxt(i,j)*flx(i,k,j,iflx_bb_xx) &
                                             + ky (j)  *flx(i,k,j,iflx_bb_yx) &
                                             + kz (k)  *flx(i,k,j,iflx_bb_zx) ) 
            nonlin(i,k,j,inonlin_bb_y) = zi*(  kxt(i,j)*flx(i,k,j,iflx_bb_xy) &
                                             + ky (j)  *flx(i,k,j,iflx_bb_yy) &
                                             + kz (k)  *flx(i,k,j,iflx_bb_zy) ) 
            nonlin(i,k,j,inonlin_bb_z) = zi*(  kxt(i,j)*flx(i,k,j,iflx_bb_xz) &
                                             + ky (j)  *flx(i,k,j,iflx_bb_yz) &
                                             + kz (k)  *flx(i,k,j,iflx_bb_zz) ) 

            nonlin(i,k,j,inonlin_ub_x) = zi*(  kxt(i,j)*flx(i,k,j,iflx_ub_xx) &
                                             + ky (j)  *flx(i,k,j,iflx_ub_yx) &
                                             + kz (k)  *flx(i,k,j,iflx_ub_zx) ) 
            nonlin(i,k,j,inonlin_ub_y) = zi*(  kxt(i,j)*flx(i,k,j,iflx_ub_xy) &
                                             + ky (j)  *flx(i,k,j,iflx_ub_yy) &
                                             + kz (k)  *flx(i,k,j,iflx_ub_zy) ) 
            nonlin(i,k,j,inonlin_ub_z) = zi*(  kxt(i,j)*flx(i,k,j,iflx_ub_xz) &
                                             + ky (j)  *flx(i,k,j,iflx_ub_yz) &
                                             + kz (k)  *flx(i,k,j,iflx_ub_zz) ) 

            nonlin(i,k,j,inonlin_bu_x) = zi*(  kxt(i,j)*flx(i,k,j,iflx_bu_xx) &
                                             + ky (j)  *flx(i,k,j,iflx_bu_yx) &
                                             + kz (k)  *flx(i,k,j,iflx_bu_zx) ) 
            nonlin(i,k,j,inonlin_bu_y) = zi*(  kxt(i,j)*flx(i,k,j,iflx_bu_xy) &
                                             + ky (j)  *flx(i,k,j,iflx_bu_yy) &
                                             + kz (k)  *flx(i,k,j,iflx_bu_zy) ) 
            nonlin(i,k,j,inonlin_bu_z) = zi*(  kxt(i,j)*flx(i,k,j,iflx_bu_xz) &
                                             + ky (j)  *flx(i,k,j,iflx_bu_yz) &
                                             + kz (k)  *flx(i,k,j,iflx_bu_zz) ) 

          enddo
        enddo
      enddo
      !$omp end parallel do
      if(allocated(flx)) deallocate(flx)

      ! get nonlinear transfer for each kprp_log(jj)
      do jj = 1, nkpolar_log
        do j = iky_st, iky_en
          do k = ikz_st, ikz_en
            do i = ikx_st, ikx_en

              if(jj == nkpolar_log) then
                if(k2t(i, k, j) >= (kpbin_log(jj))**2) then
                  filter = 1.d0
                else
                  filter = 0.d0
                endif
              else
                if(k2t(i, k, j) >= (kpbin_log(jj))**2 .and. k2t(i, k, j) < (kpbin_log(jj + 1))**2) then
                  filter = 1.d0
                else
                  filter = 0.d0
                endif
              endif

              ! The reason for the following treatment for kx == 0 mode is the following. Compile it with LaTeX.
              !-----------------------------------------------------------------------------------------------------------------------------------
              ! The volume integral of a quadratic function is
              ! \[\int \mathrm{d}^3\mathbf{r}\, f(x,y)^2 = \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2}\sum_{k_y = -n_{k_y}/2}^{n_{k_y}/2} 
              ! |f_{k_x, k_y}|^2 = \left( \sum_{k_y = -n_{k_y}/2}^{-1}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 0}
              ! \sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} + \sum_{k_y = 1}^{n_{k_y}/2}\sum_{k_x = -n_{k_x}/2}^{n_{k_x}/2} \right) |f_{k_x, k_y}|^2\]
              ! Since FFTW only computes the second and third terms, we need to compensate the first term, which is equivalent to the third term.
              !-----------------------------------------------------------------------------------------------------------------------------------
              if (j /= 1) filter = filter * 2

              trans_uu(jj, ii) = trans_uu(jj, ii) - filter*dble( &
                                    w(i,k,j,iux)*conjg(nonlin(i,k,j,inonlin_uu_x)) &
                                  + w(i,k,j,iuy)*conjg(nonlin(i,k,j,inonlin_uu_y)) &
                                  + w(i,k,j,iuz)*conjg(nonlin(i,k,j,inonlin_uu_z)) &
                                 )                                                                                                 
                                                                                                                                   
              trans_bb(jj, ii) = trans_bb(jj, ii) - filter*dble( &                                                                 
                                    w(i,k,j,ibx)*conjg(nonlin(i,k,j,inonlin_ub_x)) &
                                  + w(i,k,j,iby)*conjg(nonlin(i,k,j,inonlin_ub_y)) &
                                  + w(i,k,j,ibz)*conjg(nonlin(i,k,j,inonlin_ub_z)) &
                                 )                                                                                                 
                                                                                                                                   
              trans_ub(jj, ii) = trans_ub(jj, ii) + filter*dble( &                                                                 
                                    w(i,k,j,ibx)*conjg(nonlin(i,k,j,inonlin_bu_x)) &
                                  + w(i,k,j,iby)*conjg(nonlin(i,k,j,inonlin_bu_y)) &
                                  + w(i,k,j,ibz)*conjg(nonlin(i,k,j,inonlin_bu_z)) &
                                 )                                                                                                 
                                                                                                                                   
              trans_bu(jj, ii) = trans_bu(jj, ii) + filter*dble( &                                                                 
                                    w(i,k,j,iux)*conjg(nonlin(i,k,j,inonlin_bb_x)) &
                                  + w(i,k,j,iuy)*conjg(nonlin(i,k,j,inonlin_bb_y)) &
                                  + w(i,k,j,iuz)*conjg(nonlin(i,k,j,inonlin_bb_z)) &
                                 )

            enddo
          enddo
        enddo
      enddo
      if(allocated(nonlin)) deallocate(nonlin)

    enddo

    call sum_reduce(trans_uu, 0)
    call sum_reduce(trans_bb, 0)
    call sum_reduce(trans_ub, 0)
    call sum_reduce(trans_bu, 0)

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    if (proc0) call put_time_stamp(timer_diagnostics_nltrans)

    call loop_io_nltrans(nkpolar_log, trans_uu, trans_bb, trans_ub, trans_bu)

    deallocate(trans_uu )
    deallocate(trans_bb )
    deallocate(trans_ub )
    deallocate(trans_bu )
    if(allocated(w))   deallocate(w  )
    if(allocated(w_r)) deallocate(w_r)

  end subroutine loop_diagnostics_nltrans

end module diagnostics




