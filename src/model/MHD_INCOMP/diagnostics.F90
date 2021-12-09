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
  public :: loop_diagnostics, loop_diagnostics_fields_secion, loop_diagnostics_SF2

  private

  integer  :: nl, nsamp_r0, nsamp_ang
  real(r8), allocatable :: ll(:), lpar(:), lper(:)
contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of diagnostics
!-----------------------------------------------!
  subroutine init_diagnostics
    use diagnostics_common, only: init_polar_spectrum_2d, init_polar_spectrum_3d
    use diagnostics_common, only: init_series_modes
    use io, only: init_io 
    use grid, only: nlz
    implicit none

    if(nlz == 2) then
      call init_polar_spectrum_2d
    else
      call init_polar_spectrum_3d
    endif
    call init_SF2
    call init_io(nkpolar, kpbin, nl, lpar, lper)
    call init_series_modes
  end subroutine init_diagnostics


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
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use fields, only: ux_old1, uy_old1, uz_old1
    use fields, only: bx_old1, by_old1, bz_old1
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics
    use params, only: zi, nu, nu_exp, eta, eta_exp, shear, q
    use force, only: fux, fuy, fuz, fux_old1, fuy_old1, fuz_old1
    use shearing_box, only: shear_flg, tsc, nremap, k2t
    implicit none
    integer :: i, j, k

    real(r8)   , allocatable, dimension(:,:,:) :: u2, ux2, uy2, uz2
    real(r8)   , allocatable, dimension(:,:,:) :: b2, bx2, by2, bz2
    real(r8)   , allocatable, dimension(:,:,:) :: u2old1, b2old1
    real(r8)   , allocatable, dimension(:,:,:) :: u2dissip, b2dissip
    real(r8)   , allocatable, dimension(:,:,:) :: p_u
    real(r8)   , allocatable, dimension(:,:,:) :: zp2, zm2

    real(r8) :: u2_sum, b2_sum
    real(r8) :: u2dot_sum, b2dot_sum
    real(r8) :: u2dissip_sum, b2dissip_sum
    real(r8) :: p_u_sum
    real(r8) :: zp2_sum, zm2_sum
    real(r8) :: bx0, by0, bz0 ! mean magnetic field

    real(r8), dimension(:), allocatable :: u2_bin, ux2_bin, uy2_bin, uz2_bin
    real(r8), dimension(:), allocatable :: b2_bin, bx2_bin, by2_bin, bz2_bin
    real(r8), dimension(:), allocatable :: zp2_bin, zm2_bin
    complex(r8) ::  ux_mid,  uy_mid,  uz_mid
    complex(r8) ::  bx_mid,  by_mid,  bz_mid
    complex(r8) :: fux_mid, fuy_mid, fuz_mid

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics)

    allocate(u2      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); u2       = 0.d0
    allocate(ux2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); ux2      = 0.d0
    allocate(uy2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); uy2      = 0.d0
    allocate(uz2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); uz2      = 0.d0
    allocate(b2      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); b2       = 0.d0
    allocate(bx2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bx2      = 0.d0
    allocate(by2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); by2      = 0.d0
    allocate(bz2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bz2      = 0.d0
    allocate(u2old1  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); u2old1   = 0.d0
    allocate(b2old1  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); b2old1   = 0.d0
    allocate(u2dissip(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); u2dissip = 0.d0
    allocate(b2dissip(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); b2dissip = 0.d0
    allocate(p_u     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); p_u      = 0.d0
    allocate(zp2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); zp2      = 0.d0
    allocate(zm2     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); zm2      = 0.d0

    allocate( u2_bin(1:nkpolar));  u2_bin = 0.d0
    allocate(ux2_bin(1:nkpolar)); ux2_bin = 0.d0
    allocate(uy2_bin(1:nkpolar)); uy2_bin = 0.d0
    allocate(uz2_bin(1:nkpolar)); uz2_bin = 0.d0
    allocate( b2_bin(1:nkpolar));  b2_bin = 0.d0
    allocate(bx2_bin(1:nkpolar)); bx2_bin = 0.d0
    allocate(by2_bin(1:nkpolar)); by2_bin = 0.d0
    allocate(bz2_bin(1:nkpolar)); bz2_bin = 0.d0
    allocate(zp2_bin(1:nkpolar)); zp2_bin = 0.d0
    allocate(zm2_bin(1:nkpolar)); zm2_bin = 0.d0

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
            bx0 = abs(bx(i,k,j))
            by0 = abs(by(i,k,j))
            bz0 = abs(bz(i,k,j))
          endif
          u2    (i, k, j) = 0.5d0*(abs(ux     (i, k, j))**2 + abs(uy     (i, k, j))**2 + abs(uz     (i, k, j))**2)
          ux2   (i, k, j) = 0.5d0*(abs(ux     (i, k, j))**2)
          uy2   (i, k, j) = 0.5d0*(abs(uy     (i, k, j))**2)
          uz2   (i, k, j) = 0.5d0*(abs(uz     (i, k, j))**2)
          b2    (i, k, j) = 0.5d0*(abs(bx     (i, k, j))**2 + abs(by     (i, k, j))**2 + abs(bz     (i, k, j))**2)
          bx2   (i, k, j) = 0.5d0*(abs(bx     (i, k, j))**2)
          by2   (i, k, j) = 0.5d0*(abs(by     (i, k, j))**2)
          bz2   (i, k, j) = 0.5d0*(abs(bz     (i, k, j))**2)
          u2old1(i, k, j) = 0.5d0*(abs(ux_old1(i, k, j))**2 + abs(uy_old1(i, k, j))**2 + abs(uz_old1(i, k, j))**2)
          b2old1(i, k, j) = 0.5d0*(abs(bx_old1(i, k, j))**2 + abs(by_old1(i, k, j))**2 + abs(bz_old1(i, k, j))**2)
          zp2   (i, k, j) = abs(ux(i,k,j) + bx(i,k,j))**2 + abs(uy(i,k,j) + by(i,k,j))**2 + abs(uz(i,k,j) + bz(i,k,j))**2
          zm2   (i, k, j) = abs(ux(i,k,j) - bx(i,k,j))**2 + abs(uy(i,k,j) - by(i,k,j))**2 + abs(uz(i,k,j) - bz(i,k,j))**2

           ux_mid  = 0.5d0*( ux(i, k, j) +  ux_old1(i, k, j))
           uy_mid  = 0.5d0*( uy(i, k, j) +  uy_old1(i, k, j))
           uz_mid  = 0.5d0*( uz(i, k, j) +  uz_old1(i, k, j))
           bx_mid  = 0.5d0*( bx(i, k, j) +  bx_old1(i, k, j))
           by_mid  = 0.5d0*( by(i, k, j) +  by_old1(i, k, j))
           bz_mid  = 0.5d0*( bz(i, k, j) +  bz_old1(i, k, j))
          fux_mid  = 0.5d0*(fux(i, k, j) + fux_old1(i, k, j))
          fuy_mid  = 0.5d0*(fuy(i, k, j) + fuy_old1(i, k, j))
          fuz_mid  = 0.5d0*(fuz(i, k, j) + fuz_old1(i, k, j))

          u2dissip(i, k, j) = nu *(k2t(i, k, j)/k2_max)**nu_exp *(abs(ux_mid)**2 + abs(uy_mid)**2 + abs(uz_mid)**2)
          b2dissip(i, k, j) = eta*(k2t(i, k, j)/k2_max)**eta_exp*(abs(bx_mid)**2 + abs(by_mid)**2 + abs(bz_mid)**2)
          p_u     (i, k, j) = 0.5d0*( &
                                  (fux_mid*conjg(ux_mid) + conjg(fux_mid)*ux_mid) &
                                + (fuy_mid*conjg(uy_mid) + conjg(fuy_mid)*uy_mid) &
                                + (fuz_mid*conjg(uz_mid) + conjg(fuz_mid)*uz_mid) &
                                + q*shear_flg*(   ux_mid*conjg(uy_mid) + conjg(ux_mid)*uy_mid &
                                                - bx_mid*conjg(by_mid) - conjg(bx_mid)*by_mid &
                                              ) &
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
            u2      (i, k, j) = 2.0d0*u2      (i, k, j)
            ux2     (i, k, j) = 2.0d0*ux2     (i, k, j)
            uy2     (i, k, j) = 2.0d0*uy2     (i, k, j)
            uz2     (i, k, j) = 2.0d0*uz2     (i, k, j)
            b2      (i, k, j) = 2.0d0*b2      (i, k, j)
            bx2     (i, k, j) = 2.0d0*bx2     (i, k, j)
            by2     (i, k, j) = 2.0d0*by2     (i, k, j)
            bz2     (i, k, j) = 2.0d0*bz2     (i, k, j)
            u2old1  (i, k, j) = 2.0d0*u2old1  (i, k, j)
            b2old1  (i, k, j) = 2.0d0*b2old1  (i, k, j)
            u2dissip(i, k, j) = 2.0d0*u2dissip(i, k, j)
            b2dissip(i, k, j) = 2.0d0*b2dissip(i, k, j)
            p_u     (i, k, j) = 2.0d0*p_u     (i, k, j)
            zp2     (i, k, j) = 2.0d0*zp2     (i, k, j)
            zm2     (i, k, j) = 2.0d0*zm2     (i, k, j)
          endif

        end do
      end do
    end do

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    u2_sum = sum(u2); call sum_reduce(u2_sum, 0)
    b2_sum = sum(b2); call sum_reduce(b2_sum, 0)

    u2dot_sum = sum((u2 - u2old1)/dt); call sum_reduce(u2dot_sum, 0)
    b2dot_sum = sum((b2 - b2old1)/dt); call sum_reduce(b2dot_sum, 0)

    u2dissip_sum = sum(u2dissip); call sum_reduce(u2dissip_sum, 0)
    b2dissip_sum = sum(b2dissip); call sum_reduce(b2dissip_sum, 0)

    p_u_sum      = sum(p_u); call sum_reduce(p_u_sum, 0)

    zp2_sum = sum(zp2); call sum_reduce(zp2_sum, 0)
    zm2_sum = sum(zm2); call sum_reduce(zm2_sum, 0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv       bin over (kx, ky, kz)       vvvvvvvvvvvvvvvvvv!
    call get_polar_spectrum_3d( u2,  u2_bin)
    call get_polar_spectrum_3d(ux2, ux2_bin)
    call get_polar_spectrum_3d(uy2, uy2_bin)
    call get_polar_spectrum_3d(uz2, uz2_bin)
    call get_polar_spectrum_3d( b2,  b2_bin)
    call get_polar_spectrum_3d(bx2, bx2_bin)
    call get_polar_spectrum_3d(by2, by2_bin)
    call get_polar_spectrum_3d(bz2, bz2_bin)
    call get_polar_spectrum_3d(zp2, zp2_bin)
    call get_polar_spectrum_3d(zm2, zm2_bin)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    call loop_io( &
                  u2_sum, b2_sum, &
                  u2dot_sum, b2dot_sum, &
                  u2dissip_sum, b2dissip_sum, &
                  p_u_sum, &
                  zp2_sum, zm2_sum, &
                  bx0, by0, bz0, &
                  !
                  nkpolar, &
                  u2_bin, ux2_bin, uy2_bin, uz2_bin, &
                  b2_bin, bx2_bin, by2_bin, bz2_bin, &
                  zp2_bin, zm2_bin &
                )

    deallocate(u2)
    deallocate(ux2)
    deallocate(uy2)
    deallocate(uz2)
    deallocate(b2)
    deallocate(bx2)
    deallocate(by2)
    deallocate(bz2)
    deallocate(u2old1)
    deallocate(b2old1)
    deallocate(u2dissip)
    deallocate(b2dissip)
    deallocate(p_u)
    deallocate(zp2)
    deallocate(zm2)

    deallocate( u2_bin)
    deallocate(ux2_bin)
    deallocate(uy2_bin)
    deallocate(uz2_bin)
    deallocate( b2_bin)
    deallocate(bx2_bin)
    deallocate(by2_bin)
    deallocate(bz2_bin)
    deallocate(zp2_bin)
    deallocate(zm2_bin)

    if (proc0) call put_time_stamp(timer_diagnostics)
  end subroutine loop_diagnostics


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Diagnostics for cross section of fileds
!-----------------------------------------------!
  subroutine loop_diagnostics_fields_secion
    use io, only: loop_io_fields_section
    use utils, only: curl
    use mp, only: proc0
    use grid, only: nlx, nly, nlz, nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics
    use params, only: shear
    use shearing_box, only: to_non_shearing_coordinate, tsc, nremap
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

    integer :: i, j, k

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics)

    allocate(f (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f   = 0.d0
    allocate(fr(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); fr  = 0.d0

    allocate(ux_r_z0(nlx, nly)); ux_r_z0 = 0.d0
    allocate(ux_r_x0(nly, nlz)); ux_r_x0 = 0.d0
    allocate(ux_r_y0(nlx, nlz)); ux_r_y0 = 0.d0

    allocate(uy_r_z0(nlx, nly)); uy_r_z0 = 0.d0
    allocate(uy_r_x0(nly, nlz)); uy_r_x0 = 0.d0
    allocate(uy_r_y0(nlx, nlz)); uy_r_y0 = 0.d0

    allocate(uz_r_z0(nlx, nly)); uz_r_z0 = 0.d0
    allocate(uz_r_x0(nly, nlz)); uz_r_x0 = 0.d0
    allocate(uz_r_y0(nlx, nlz)); uz_r_y0 = 0.d0

    allocate(wx_r_z0(nlx, nly)); wx_r_z0 = 0.d0
    allocate(wx_r_x0(nly, nlz)); wx_r_x0 = 0.d0
    allocate(wx_r_y0(nlx, nlz)); wx_r_y0 = 0.d0

    allocate(wy_r_z0(nlx, nly)); wy_r_z0 = 0.d0
    allocate(wy_r_x0(nly, nlz)); wy_r_x0 = 0.d0
    allocate(wy_r_y0(nlx, nlz)); wy_r_y0 = 0.d0

    allocate(wz_r_z0(nlx, nly)); wz_r_z0 = 0.d0
    allocate(wz_r_x0(nly, nlz)); wz_r_x0 = 0.d0
    allocate(wz_r_y0(nlx, nlz)); wz_r_y0 = 0.d0

    allocate(bx_r_z0(nlx, nly)); bx_r_z0 = 0.d0
    allocate(bx_r_x0(nly, nlz)); bx_r_x0 = 0.d0
    allocate(bx_r_y0(nlx, nlz)); bx_r_y0 = 0.d0

    allocate(by_r_z0(nlx, nly)); by_r_z0 = 0.d0
    allocate(by_r_x0(nly, nlz)); by_r_x0 = 0.d0
    allocate(by_r_y0(nlx, nlz)); by_r_y0 = 0.d0

    allocate(bz_r_z0(nlx, nly)); bz_r_z0 = 0.d0
    allocate(bz_r_x0(nly, nlz)); bz_r_x0 = 0.d0
    allocate(bz_r_y0(nlx, nlz)); bz_r_y0 = 0.d0

    allocate(jx_r_z0(nlx, nly)); jx_r_z0 = 0.d0
    allocate(jx_r_x0(nly, nlz)); jx_r_x0 = 0.d0
    allocate(jx_r_y0(nlx, nlz)); jx_r_y0 = 0.d0

    allocate(jy_r_z0(nlx, nly)); jy_r_z0 = 0.d0
    allocate(jy_r_x0(nly, nlz)); jy_r_x0 = 0.d0
    allocate(jy_r_y0(nlx, nlz)); jy_r_y0 = 0.d0

    allocate(jz_r_z0(nlx, nly)); jz_r_z0 = 0.d0
    allocate(jz_r_x0(nly, nlz)); jz_r_x0 = 0.d0
    allocate(jz_r_y0(nlx, nlz)); jz_r_y0 = 0.d0

    allocate(wx (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); wx  = 0.d0
    allocate(wy (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); wy  = 0.d0
    allocate(wz (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); wz  = 0.d0
    allocate(jx (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); jx  = 0.d0
    allocate(jy (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); jy  = 0.d0
    allocate(jz (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); jz  = 0.d0

    allocate(u2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); u2  = 0.d0
    allocate(b2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); b2  = 0.d0

    allocate(u2_kxy(nkx, nky)); u2_kxy = 0.d0
    allocate(u2_kyz(nky, nkz)); u2_kyz = 0.d0
    allocate(u2_kxz(nkx, nkz)); u2_kxz = 0.d0
    allocate(b2_kxy(nkx, nky)); b2_kxy = 0.d0
    allocate(b2_kyz(nky, nkz)); b2_kyz = 0.d0
    allocate(b2_kxz(nkx, nkz)); b2_kxz = 0.d0

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
          u2(i, k, j) = 0.5d0*(abs(ux(i, k, j))**2 + abs(uy(i, k, j))**2 + abs(uz(i, k, j))**2)
          b2(i, k, j) = 0.5d0*(abs(bx(i, k, j))**2 + abs(by(i, k, j))**2 + abs(bz(i, k, j))**2)
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

    call loop_io_fields_section( &
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

    if (proc0) call put_time_stamp(timer_diagnostics)
  end subroutine loop_diagnostics_fields_secion


!-----------------------------------------------!
!> @author  YK
!! @date    4 Jul 2021
!! @brief   Initialize order structure function
!-----------------------------------------------!
  subroutine init_SF2
    use grid, only: lx, ly, lz, dlx, dly, dlz
    use params, only: SF2_nsample
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

    allocate(lpar(nl), lper(nl), source=ll)

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
    use time_stamp, only: put_time_stamp, timer_diagnostics, timer_diagnostics_SF2
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

    if (proc0) call put_time_stamp(timer_diagnostics)
    if (proc0) call put_time_stamp(timer_diagnostics_SF2)

    allocate( f   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f    = 0.d0
    allocate( bx_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); bx_r = 0.d0
    allocate( by_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); by_r = 0.d0
    allocate( bz_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); bz_r = 0.d0
    allocate( ux_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); ux_r = 0.d0
    allocate( uy_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); uy_r = 0.d0
    allocate( uz_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); uz_r = 0.d0

    f = bx ; call p3dfft_btran_c2r(f, bx_r , 'tff')
    f = by ; call p3dfft_btran_c2r(f, by_r , 'tff')
    f = bz ; call p3dfft_btran_c2r(f, bz_r , 'tff')
    f = ux ; call p3dfft_btran_c2r(f, ux_r , 'tff')
    f = uy ; call p3dfft_btran_c2r(f, uy_r , 'tff')
    f = uz ; call p3dfft_btran_c2r(f, uz_r , 'tff')

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

    call loop_io_SF2(nl, sf2b, sf2u)

    deallocate(f)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)

    if (proc0) call put_time_stamp(timer_diagnostics)
    if (proc0) call put_time_stamp(timer_diagnostics_SF2)
  end subroutine loop_diagnostics_SF2

end module diagnostics




