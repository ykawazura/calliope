!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../diagnostics_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Diagnostics for MHD_COMP_ISOTH
!-----------------------------------------------!
module diagnostics
  use diagnostics_common
  use p3dfft
  implicit none

  public :: init_diagnostics, finish_diagnostics
  public :: loop_diagnostics, loop_diagnostics_2D, loop_diagnostics_kpar, loop_diagnostics_SF2

  private
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
    call init_io(nkpolar, kpbin)
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
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use fields, only: rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: rho_old
    use fields, only: mx_old, my_old, mz_old
    use fields, only: bx_old, by_old, bz_old
    use fields, only: m_to_u
    use mp, only: sum_reduce
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: zi, nu, nu_exp, eta, eta_exp, lmd, lmd_exp, cs2va2, shear, q
    use force, only: fmx, fmy, fmz, fmx_old, fmy_old, fmz_old
    use shearing_box, only: shear_flg, tsc, nremap, k2t
    implicit none
    integer :: i, j, k


    real(r8)   , allocatable, dimension(:,:,:) :: wkin, wmag
    real(r8)   , allocatable, dimension(:,:,:) :: wkin_old, wmag_old
    real(r8)   , allocatable, dimension(:,:,:) :: wkin_dissip, wmag_dissip, wrho_dissip
    real(r8)   , allocatable, dimension(:,:,:) :: wrho
    real(r8)   , allocatable, dimension(:,:,:) :: wrho_old
    real(r8)   , allocatable, dimension(:,:,:) :: u2, ux2, uy2, uz2
    real(r8)   , allocatable, dimension(:,:,:) :: b2, bx2, by2, bz2
    real(r8)   , allocatable, dimension(:,:,:) :: p_ext, p_re, p_ma 
    real(r8)   , allocatable, dimension(:,:,:) :: zp2, zm2
    real(r8)   , allocatable, dimension(:,:,:) :: src_r

    complex(r8), allocatable, dimension(:,:,:) :: ux, uy, uz
    complex(r8), allocatable, dimension(:,:,:) :: ux_old, uy_old, uz_old
    complex(r8), allocatable, dimension(:,:,:) :: lnrho, lnrho_old
    complex(r8), allocatable, dimension(:,:,:) :: u2half, u2half_old
    complex(r8), allocatable, dimension(:,:,:) :: f
    complex(r8), allocatable, dimension(:,:,:) :: src_c

    real(r8)   , allocatable, dimension(:,:,:) :: rho_r, ux_r, uy_r, uz_r, bx_r, by_r, bz_r

    real(r8) :: wkin_sum, wmag_sum
    real(r8) :: wkin_dot_sum, wmag_dot_sum
    real(r8) :: wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum
    real(r8) :: wrho_sum
    real(r8) :: wrho_dot_sum
    real(r8) :: p_ext_sum, p_re_sum, p_ma_sum
    real(r8) :: zp2_sum, zm2_sum

    real(r8) :: rho_rms, u_rms, b_rms
    real(r8) :: smach_rms, amach_rms, beta_rms
    real(r8) :: bx0, by0, bz0 ! mean magnetic field

    real(r8), dimension(:), allocatable :: u2_bin, ux2_bin, uy2_bin, uz2_bin
    real(r8), dimension(:), allocatable :: b2_bin, bx2_bin, by2_bin, bz2_bin
    real(r8), dimension(:), allocatable :: rho_bin
    real(r8), dimension(:), allocatable :: zp2_bin, zm2_bin
    complex(r8) ::  mx_mid,    my_mid,  mz_mid
    complex(r8) ::  ux_mid,    uy_mid,  uz_mid
    complex(r8) ::  bx_mid,    by_mid,  bz_mid
    complex(r8) :: rho_mid, lnrho_mid,  u2half_mid
    complex(r8) :: fmx_mid, fmy_mid, fmz_mid

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    allocate(ux         , source=src_c)
    allocate(uy         , source=src_c)
    allocate(uz         , source=src_c)
    allocate(ux_old     , source=src_c)
    allocate(uy_old     , source=src_c)
    allocate(uz_old     , source=src_c)
    allocate(lnrho      , source=src_c)
    allocate(lnrho_old  , source=src_c)
    allocate(u2half     , source=src_c)
    allocate(u2half_old , source=src_c)
    allocate(f          , source=src_c)
    deallocate(src_c)

    allocate(src_r(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(wkin       , source=src_r)
    allocate(wkin_old   , source=src_r)
    allocate(wkin_dissip, source=src_r)
    allocate(wmag       , source=src_r)
    allocate(wmag_old   , source=src_r)
    allocate(wmag_dissip, source=src_r)
    allocate(wrho       , source=src_r)
    allocate(wrho_old   , source=src_r)
    allocate(wrho_dissip, source=src_r)
    allocate(u2         , source=src_r)
    allocate(ux2        , source=src_r)
    allocate(uy2        , source=src_r)
    allocate(uz2        , source=src_r)
    allocate(b2         , source=src_r)
    allocate(bx2        , source=src_r)
    allocate(by2        , source=src_r)
    allocate(bz2        , source=src_r)
    allocate(p_ext      , source=src_r)
    allocate(p_re       , source=src_r)
    allocate(p_ma       , source=src_r)
    allocate(zp2        , source=src_r)
    allocate(zm2        , source=src_r)
    deallocate(src_r)


    allocate( u2_bin(1:nkpolar));  u2_bin = 0.d0
    allocate(ux2_bin(1:nkpolar)); ux2_bin = 0.d0
    allocate(uy2_bin(1:nkpolar)); uy2_bin = 0.d0
    allocate(uz2_bin(1:nkpolar)); uz2_bin = 0.d0
    allocate( b2_bin(1:nkpolar));  b2_bin = 0.d0
    allocate(bx2_bin(1:nkpolar)); bx2_bin = 0.d0
    allocate(by2_bin(1:nkpolar)); by2_bin = 0.d0
    allocate(bz2_bin(1:nkpolar)); bz2_bin = 0.d0
    allocate(rho_bin(1:nkpolar)); rho_bin = 0.d0
    allocate(zp2_bin(1:nkpolar)); zp2_bin = 0.d0
    allocate(zm2_bin(1:nkpolar)); zm2_bin = 0.d0

    allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(rho_r, source=src_r)
    allocate( ux_r, source=src_r)
    allocate( uy_r, source=src_r)
    allocate( uz_r, source=src_r)
    allocate( bx_r, source=src_r)
    allocate( by_r, source=src_r)
    allocate( bz_r, source=src_r)
    deallocate(src_r)

    call rho_to_lnrho(rho    , lnrho    )
    call rho_to_lnrho(rho_old, lnrho_old)
    call m_to_u(rho    , mx    , my    , mz    , ux    , uy    , uz    )
    call m_to_u(rho_old, mx_old, my_old, mz_old, ux_old, uy_old, uz_old)
    call get_u2half(ux    , uy    , uz    , u2half    )
    call get_u2half(ux_old, uy_old, uz_old, u2half_old)

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
            bx0 = abs(bx(i,k,j))
            by0 = abs(by(i,k,j))
            bz0 = abs(bz(i,k,j))
          endif
          wkin    (i, k, j) = 0.25d0*( & 
                                        mx    (i, k, j)*conjg(ux    (i, k, j)) + ux    (i, k, j)*conjg(mx    (i, k, j)) &
                                      + my    (i, k, j)*conjg(uy    (i, k, j)) + uy    (i, k, j)*conjg(my    (i, k, j)) &
                                      + mz    (i, k, j)*conjg(uz    (i, k, j)) + uz    (i, k, j)*conjg(mz    (i, k, j)) &
                                     )
          wkin_old(i, k, j) = 0.25d0*( & 
                                        mx_old(i, k, j)*conjg(ux_old(i, k, j)) + ux_old(i, k, j)*conjg(mx_old(i, k, j)) &
                                      + my_old(i, k, j)*conjg(uy_old(i, k, j)) + uy_old(i, k, j)*conjg(my_old(i, k, j)) &
                                      + mz_old(i, k, j)*conjg(uz_old(i, k, j)) + uz_old(i, k, j)*conjg(mz_old(i, k, j)) &
                                     )
          wmag    (i, k, j) = 0.5d0*( & 
                                        bx    (i, k, j)*conjg(bx    (i, k, j)) &
                                      + by    (i, k, j)*conjg(by    (i, k, j)) &
                                      + bz    (i, k, j)*conjg(bz    (i, k, j)) &
                                     )
          wmag_old(i, k, j) = 0.5d0*( & 
                                        bx_old(i, k, j)*conjg(bx_old(i, k, j)) &
                                      + by_old(i, k, j)*conjg(by_old(i, k, j)) &
                                      + bz_old(i, k, j)*conjg(bz_old(i, k, j)) &
                                     )
          u2      (i, k, j) = 0.5d0*(abs(ux(i, k, j))**2 + abs(uy(i, k, j))**2 + abs(uz(i, k, j))**2)
          ux2     (i, k, j) = 0.5d0*(abs(ux(i, k, j))**2)
          uy2     (i, k, j) = 0.5d0*(abs(uy(i, k, j))**2)
          uz2     (i, k, j) = 0.5d0*(abs(uz(i, k, j))**2)
          b2      (i, k, j) = 0.5d0*(abs(bx(i, k, j))**2 + abs(by(i, k, j))**2 + abs(bz(i, k, j))**2)
          bx2     (i, k, j) = 0.5d0*(abs(bx(i, k, j))**2)
          by2     (i, k, j) = 0.5d0*(abs(by(i, k, j))**2)
          bz2     (i, k, j) = 0.5d0*(abs(bz(i, k, j))**2)

          wrho    (i, k, j) = 0.5d0*cs2va2*( & 
                                   rho    (i, k, j)*conjg(lnrho    (i, k, j)) + lnrho    (i, k, j)*conjg(rho    (i, k, j)) &
                                 )
          wrho_old(i, k, j) = 0.5d0*cs2va2*( & 
                                   rho_old(i, k, j)*conjg(lnrho_old(i, k, j)) + lnrho_old(i, k, j)*conjg(rho_old(i, k, j)) &
                                 )
          zp2     (i, k, j) = abs(ux(i,k,j) + bx(i,k,j))**2 + abs(uy(i,k,j) + by(i,k,j))**2 + abs(uz(i,k,j) + bz(i,k,j))**2
          zm2     (i, k, j) = abs(ux(i,k,j) - bx(i,k,j))**2 + abs(uy(i,k,j) - by(i,k,j))**2 + abs(uz(i,k,j) - bz(i,k,j))**2

           mx_mid    = 0.5d0*( mx   (i, k, j) +  mx_old   (i, k, j))
           my_mid    = 0.5d0*( my   (i, k, j) +  my_old   (i, k, j))
           mz_mid    = 0.5d0*( mz   (i, k, j) +  mz_old   (i, k, j))
           ux_mid    = 0.5d0*( ux   (i, k, j) +  ux_old   (i, k, j))
           uy_mid    = 0.5d0*( uy   (i, k, j) +  uy_old   (i, k, j))
           uz_mid    = 0.5d0*( uz   (i, k, j) +  uz_old   (i, k, j))
           bx_mid    = 0.5d0*( bx   (i, k, j) +  bx_old   (i, k, j))
           by_mid    = 0.5d0*( by   (i, k, j) +  by_old   (i, k, j))
           bz_mid    = 0.5d0*( bz   (i, k, j) +  bz_old   (i, k, j))
          fmx_mid    = 0.5d0*(fmx   (i, k, j) + fmx_old   (i, k, j))
          fmy_mid    = 0.5d0*(fmy   (i, k, j) + fmy_old   (i, k, j))
          fmz_mid    = 0.5d0*(fmz   (i, k, j) + fmz_old   (i, k, j))
          rho_mid    = 0.5d0*(rho   (i, k, j) + rho_old   (i, k, j))
          lnrho_mid  = 0.5d0*(lnrho (i, k, j) + lnrho_old (i, k, j))
          u2half_mid = 0.5d0*(u2half(i, k, j) + u2half_old(i, k, j))

          wkin_dissip(i, k, j) = nu*(k2t(i, k, j)/k2_max)**nu_exp*0.5d0*( &
                                                                      mx_mid*conjg(ux_mid) + ux_mid*conjg(mx_mid) &
                                                                    + my_mid*conjg(uy_mid) + uy_mid*conjg(my_mid) &
                                                                    + mz_mid*conjg(uz_mid) + uz_mid*conjg(mz_mid) &
                                                                   )
          wmag_dissip(i, k, j) = eta*(k2t(i, k, j)/k2_max)**eta_exp*(abs(bx_mid)**2 + abs(by_mid)**2 + abs(bz_mid)**2)
          wrho_dissip(i, k, j) = lmd*(k2t(i, k, j)/k2_max)**lmd_exp*0.5d0*( &
                                                                        (-u2half_mid + cs2va2*lnrho_mid)*conjg(rho_mid) &
                                                                      + rho_mid*conjg(-u2half_mid + cs2va2*lnrho_mid) &
                                                                    )
          p_ext   (i, k, j) = 0.5d0*( &
                                        (fmx_mid*conjg(ux_mid) + conjg(fmx_mid)*ux_mid) &
                                      + (fmy_mid*conjg(uy_mid) + conjg(fmy_mid)*uy_mid) &
                                      + (fmz_mid*conjg(uz_mid) + conjg(fmz_mid)*uz_mid) &
                                    )

          p_re    (i, k, j) = + 0.5d0*q*shear_flg*(ux_mid*conjg(my_mid) + conjg(ux_mid)*my_mid)
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
            wkin       (i, k, j) = 2.0d0*wkin       (i, k, j)
            wkin_old   (i, k, j) = 2.0d0*wkin_old   (i, k, j)
            wkin_dissip(i, k, j) = 2.0d0*wkin_dissip(i, k, j)
            wmag       (i, k, j) = 2.0d0*wmag       (i, k, j)
            wmag_old   (i, k, j) = 2.0d0*wmag_old   (i, k, j)
            wmag_dissip(i, k, j) = 2.0d0*wmag_dissip(i, k, j)
            wrho       (i, k, j) = 2.0d0*wrho       (i, k, j)
            wrho_old   (i, k, j) = 2.0d0*wrho_old   (i, k, j)
            wrho_dissip(i, k, j) = 2.0d0*wrho_dissip(i, k, j)
            u2         (i, k, j) = 2.0d0*u2         (i, k, j)
            ux2        (i, k, j) = 2.0d0*ux2        (i, k, j)
            uy2        (i, k, j) = 2.0d0*uy2        (i, k, j)
            uz2        (i, k, j) = 2.0d0*uz2        (i, k, j)
            b2         (i, k, j) = 2.0d0*b2         (i, k, j)
            bx2        (i, k, j) = 2.0d0*bx2        (i, k, j)
            by2        (i, k, j) = 2.0d0*by2        (i, k, j)
            bz2        (i, k, j) = 2.0d0*bz2        (i, k, j)
            p_ext      (i, k, j) = 2.0d0*p_ext      (i, k, j)
            p_re       (i, k, j) = 2.0d0*p_re       (i, k, j)
            p_ma       (i, k, j) = 2.0d0*p_ma       (i, k, j)
            zp2        (i, k, j) = 2.0d0*zp2        (i, k, j)
            zm2        (i, k, j) = 2.0d0*zm2        (i, k, j)
          endif

        end do
      end do
    end do

    !vvvvvvvvvvvvvvvvvv     integrate over kx, ky, kz     vvvvvvvvvvvvvvvvvv!
    wkin_sum        = sum(wkin); call sum_reduce(wkin_sum, 0)
    wmag_sum        = sum(wmag); call sum_reduce(wmag_sum, 0)
    wrho_sum        = sum(wrho); call sum_reduce(wrho_sum, 0)

    wkin_dot_sum    = sum((wkin - wkin_old)/dt); call sum_reduce(wkin_dot_sum, 0)
    wmag_dot_sum    = sum((wmag - wmag_old)/dt); call sum_reduce(wmag_dot_sum, 0)
    wrho_dot_sum    = sum((wrho - wrho_old)/dt); call sum_reduce(wrho_dot_sum, 0)

    wkin_dissip_sum = sum(wkin_dissip); call sum_reduce(wkin_dissip_sum, 0)
    wmag_dissip_sum = sum(wmag_dissip); call sum_reduce(wmag_dissip_sum, 0)
    wrho_dissip_sum = sum(wrho_dissip); call sum_reduce(wrho_dissip_sum, 0)

    p_ext_sum       = sum(p_ext); call sum_reduce(p_ext_sum, 0)
    p_re_sum        = sum(p_re ); call sum_reduce(p_re_sum , 0)
    p_ma_sum        = sum(p_ma ); call sum_reduce(p_ma_sum , 0)

    zp2_sum         = sum(zp2); call sum_reduce(zp2_sum, 0)
    zm2_sum         = sum(zm2); call sum_reduce(zm2_sum, 0)
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
    call get_polar_spectrum_3d(abs(rho), rho_bin)
    call get_polar_spectrum_3d(zp2, zp2_bin)
    call get_polar_spectrum_3d(zm2, zm2_bin)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvv         maximum mach values       vvvvvvvvvvvvvvvvvv!
    f = rho; call p3dfft_btran_c2r(f, rho_r, 'tff')
    f = ux ; call p3dfft_btran_c2r(f, ux_r , 'tff')
    f = uy ; call p3dfft_btran_c2r(f, uy_r , 'tff')
    f = uz ; call p3dfft_btran_c2r(f, uz_r , 'tff')
    f = bx ; call p3dfft_btran_c2r(f, bx_r , 'tff')
    f = by ; call p3dfft_btran_c2r(f, by_r , 'tff')
    f = bz ; call p3dfft_btran_c2r(f, bz_r , 'tff')

    rho_rms = sum(rho_r**2)
    call sum_reduce(rho_rms, 0)
    rho_rms = sqrt(rho_rms/nlx/nly/nlz)

    u_rms = sum(ux_r**2 + uy_r**2 + uz_r**2)
    call sum_reduce(u_rms, 0)
    u_rms = sqrt(u_rms/nlx/nly/nlz)

    b_rms = sum(bx_r**2 + by_r**2 + bz_r**2)
    call sum_reduce(b_rms, 0)
    b_rms = sqrt(b_rms/nlx/nly/nlz)

    smach_rms = u_rms/sqrt(cs2va2)
    amach_rms = sqrt(rho_rms)*u_rms/b_rms
    beta_rms  = 2.d0*cs2va2*rho_rms/b_rms**2
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  
    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io( &
                  wkin_sum, wmag_sum, wrho_sum, &
                  wkin_dot_sum, wmag_dot_sum, wrho_dot_sum, &
                  wkin_dissip_sum, wmag_dissip_sum, wrho_dissip_sum, &
                  p_ext_sum, p_re_sum, p_ma_sum, &
                  smach_rms, amach_rms, beta_rms, &
                  zp2_sum, zm2_sum, &
                  bx0, by0, bz0, &
                  !
                  nkpolar, &
                  rho_bin, &
                  u2_bin, ux2_bin, uy2_bin, uz2_bin, &
                  b2_bin, bx2_bin, by2_bin, bz2_bin, &
                  zp2_bin, zm2_bin &
                  !
                )

    deallocate(ux         )
    deallocate(uy         )
    deallocate(uz         )
    deallocate(ux_old     )
    deallocate(uy_old     )
    deallocate(uz_old     )
    deallocate(lnrho      )
    deallocate(lnrho_old  )
    deallocate(u2half     )
    deallocate(u2half_old )

    deallocate(wkin       )
    deallocate(wkin_old   )
    deallocate(wkin_dissip)
    deallocate(wmag       )
    deallocate(wmag_old   )
    deallocate(wmag_dissip)
    deallocate(wrho       )
    deallocate(wrho_old   )
    deallocate(wrho_dissip)
    deallocate(u2         )
    deallocate(ux2        )
    deallocate(uy2        )
    deallocate(uz2        )
    deallocate(b2         )
    deallocate(bx2        )
    deallocate(by2        )
    deallocate(bz2        )
    deallocate(p_ext      )
    deallocate(p_re       )
    deallocate(p_ma       )
    deallocate(zp2        )
    deallocate(zm2        )


    deallocate( u2_bin)
    deallocate(ux2_bin)
    deallocate(uy2_bin)
    deallocate(uz2_bin)
    deallocate( b2_bin)
    deallocate(bx2_bin)
    deallocate(by2_bin)
    deallocate(bz2_bin)
    deallocate(rho_bin)
    deallocate(zp2_bin)
    deallocate(zm2_bin)

    deallocate(f )

    deallocate(rho_r)
    deallocate( ux_r)
    deallocate( uy_r)
    deallocate( uz_r)
    deallocate( bx_r)
    deallocate( by_r)
    deallocate( bz_r)
  end subroutine loop_diagnostics


!---------------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Diagnostics for cross section of fileds
!---------------------------------------------------!
  subroutine loop_diagnostics_2D
    use io, only: loop_io_2D
    use utils, only: curl
    use mp, only: proc0
    use grid, only: nlx, nly, nlz, nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use fields, only: rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: m_to_u
    use time, only: dt
    use time_stamp, only: put_time_stamp, timer_diagnostics_total
    use params, only: shear
    use shearing_box, only: to_non_shearing_coordinate, tsc, nremap
    implicit none

    complex(r8), allocatable, dimension(:,:,:) :: f
    real(r8)   , allocatable, dimension(:,:,:) :: fr
    complex(r8), allocatable, dimension(:,:,:) :: ux, uy, uz
    complex(r8), allocatable, dimension(:,:,:) :: wx, wy, wz, jx, jy, jz
    real(r8)   , allocatable, dimension(:,:)   :: rho_r_z0, rho_r_x0, rho_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: mx_r_z0, mx_r_x0, mx_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: my_r_z0, my_r_x0, my_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: mz_r_z0, mz_r_x0, mz_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: wx_r_z0, wx_r_x0, wx_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: wy_r_z0, wy_r_x0, wy_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: wz_r_z0, wz_r_x0, wz_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: bx_r_z0, bx_r_x0, bx_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: by_r_z0, by_r_x0, by_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: bz_r_z0, bz_r_x0, bz_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: jx_r_z0, jx_r_x0, jx_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: jy_r_z0, jy_r_x0, jy_r_y0
    real(r8)   , allocatable, dimension(:,:)   :: jz_r_z0, jz_r_x0, jz_r_y0

    real(r8)   , allocatable, dimension(:,:,:) :: u2, b2
    real(r8)   , allocatable, dimension(:,:)   :: rho_kxy, rho_kyz, rho_kxz
    real(r8)   , allocatable, dimension(:,:)   ::  u2_kxy,  u2_kyz,  u2_kxz
    real(r8)   , allocatable, dimension(:,:)   ::  b2_kxy,  b2_kyz,  b2_kxz

    real(r8)   , allocatable, dimension(:,:)   :: src1, src2, src3
    complex(r8), allocatable, dimension(:,:,:) :: src4

    integer :: i, j, k

    if(nremap > 0 .and. tsc <= 5.*dt) then
      return !skip 5 loops after remapping
    endif
    if (proc0) call put_time_stamp(timer_diagnostics_total)

    allocate(f (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f   = 0.d0
    allocate(fr(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); fr  = 0.d0

    allocate(src1(nlx, nly), source=0.d0); allocate(src2(nly, nlz), source=0.d0); allocate(src3(nlx, nlz), source=0.d0)
    allocate(rho_r_z0, source=src1)
    allocate(rho_r_x0, source=src2)
    allocate(rho_r_y0, source=src3)

    allocate(mx_r_z0, source=src1)
    allocate(mx_r_x0, source=src2)
    allocate(mx_r_y0, source=src3)
                                  
    allocate(my_r_z0, source=src1)
    allocate(my_r_x0, source=src2)
    allocate(my_r_y0, source=src3)
                                  
    allocate(mz_r_z0, source=src1)
    allocate(mz_r_x0, source=src2)
    allocate(mz_r_y0, source=src3)
                                  
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
    allocate(ux, source=src4)
    allocate(uy, source=src4)
    allocate(uz, source=src4)
    allocate(wx, source=src4)
    allocate(wy, source=src4)
    allocate(wz, source=src4)
    allocate(jx, source=src4)
    allocate(jy, source=src4)
    allocate(jz, source=src4)
    deallocate(src4)

    allocate(u2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)
    allocate(b2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)

    allocate(src1(nkx, nky), source=0.d0); allocate(src2(nky, nkz), source=0.d0); allocate(src3(nkx, nkz), source=0.d0)
    allocate(rho_kxy, source=src1)
    allocate(rho_kyz, source=src2)
    allocate(rho_kxz, source=src3)
    allocate( u2_kxy, source=src1)
    allocate( u2_kyz, source=src2)
    allocate( u2_kxz, source=src3)
    allocate( b2_kxy, source=src1)
    allocate( b2_kyz, source=src2)
    allocate( b2_kxz, source=src3)
    deallocate(src1, src2, src3)

    !vvvvvvvvvvvvvvvvvv         2D cut of fields          vvvvvvvvvvvvvvvvvv!
    call m_to_u(rho, mx, my, mz, ux, uy, uz)
    call curl(ux, uy, uz, wx, wy, wz)
    call curl(bx, by, bz, jx, jy, jz)

    f = rho; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, rho_r_z0, rho_r_x0, rho_r_y0)
    f = mx ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, mx_r_z0, mx_r_x0, mx_r_y0)
    f = my ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, my_r_z0, my_r_x0, my_r_y0)
    f = mz ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, mz_r_z0, mz_r_x0, mz_r_y0)
    f = bx ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, bx_r_z0, bx_r_x0, bx_r_y0)
    f = by ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, by_r_z0, by_r_x0, by_r_y0)
    f = bz ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, bz_r_z0, bz_r_x0, bz_r_y0)
    f = wx ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, wx_r_z0, wx_r_x0, wx_r_y0)
    f = wy ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, wy_r_z0, wy_r_x0, wy_r_y0)
    f = wz ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, wz_r_z0, wz_r_x0, wz_r_y0)
    f = jx ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, jx_r_z0, jx_r_x0, jx_r_y0)
    f = jy ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
        call cut_2d_r(fr, jy_r_z0, jy_r_x0, jy_r_y0)
    f = jz ; call p3dfft_btran_c2r(f, fr, 'tff'); if(shear) call to_non_shearing_coordinate(fr); 
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
    call sum_2d_k(abs(rho), rho_kxy, rho_kyz, rho_kxz)
    call sum_2d_k(     u2 ,  u2_kxy,  u2_kyz,  u2_kxz)
    call sum_2d_k(     b2 ,  b2_kxy,  b2_kyz,  b2_kxz)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if (proc0) call put_time_stamp(timer_diagnostics_total)
    call loop_io_2D( &
                  rho_r_z0, rho_r_x0, rho_r_y0, &
                  
                  mx_r_z0, mx_r_x0, mx_r_y0, &
                  my_r_z0, my_r_x0, my_r_y0, &
                  mz_r_z0, mz_r_x0, mz_r_y0, &
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
                  rho_kxy, rho_kyz, rho_kxz, &
                   u2_kxy,  u2_kyz,  u2_kxz, &
                   b2_kxy,  b2_kyz,  b2_kxz  &
                )

    deallocate(f )
    deallocate(fr)

    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(wx)
    deallocate(wy)
    deallocate(wz)
    deallocate(jx)
    deallocate(jy)
    deallocate(jz)

    deallocate(rho_r_z0)
    deallocate(rho_r_x0)
    deallocate(rho_r_y0)

    deallocate(mx_r_z0)
    deallocate(mx_r_x0)
    deallocate(mx_r_y0)
                                           
    deallocate(my_r_z0)
    deallocate(my_r_x0)
    deallocate(my_r_y0)
                                           
    deallocate(mz_r_z0)
    deallocate(mz_r_x0)
    deallocate(mz_r_y0)
                                           
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


    deallocate(u2)
    deallocate(b2)

    deallocate(rho_kxy)
    deallocate(rho_kyz)
    deallocate(rho_kxz)
    deallocate( u2_kxy)
    deallocate( u2_kyz)
    deallocate( u2_kxz)
    deallocate( b2_kxy)
    deallocate( b2_kyz)
    deallocate( b2_kxz)
  end subroutine loop_diagnostics_2D

!-----------------------------------------------!
!> @author  YK
!! @date    15 Dec 2020
!! @brief   rho to ln(rho)
!-----------------------------------------------!
  subroutine rho_to_lnrho(rho, lnrho)
    use grid, only: nlx, nly, nlz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: rho
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(out) :: lnrho
    real   (r8), allocatable, dimension(:,:,:) :: rho_r, lnrho_r
    complex(r8), allocatable, dimension(:,:,:) :: rho_

    allocate(rho_   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(rho_r  (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(lnrho_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    !$omp workshare
    rho_ = rho
    !$omp end workshare

    call p3dfft_btran_c2r(rho_, rho_r, 'tff')

    !$omp workshare
    lnrho_r = dlog(rho_r)
    !$omp end workshare

    call p3dfft_ftran_r2c(lnrho_r, lnrho, 'fft'); lnrho = lnrho/nlx/nly/nlz 

    deallocate(rho_)
    deallocate(rho_r)
    deallocate(lnrho_r)
  end subroutine rho_to_lnrho

!-----------------------------------------------!
!> @author  YK
!! @date    24 Mar 2021
!! @brief   get u2/2
!-----------------------------------------------!
  subroutine get_u2half(ux, uy, uz, u2half)
    use grid, only: nlx, nly, nlz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: ux, uy, uz
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(out) :: u2half
    real   (r8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r, u2half_r
    complex(r8), allocatable, dimension(:,:,:) :: ux_ , uy_ , uz_
    real   (r8), allocatable, dimension(:,:,:) :: src_r
    complex(r8), allocatable, dimension(:,:,:) :: src_c

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(ux_, source=src_c)
    allocate(uy_, source=src_c)
    allocate(uz_, source=src_c)
    deallocate(src_c)
    
    allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(ux_r    , source=src_r)
    allocate(uy_r    , source=src_r)
    allocate(uz_r    , source=src_r)
    allocate(u2half_r, source=src_r)
    deallocate(src_r)

    !$omp workshare
    ux_ = ux
    uy_ = uy
    uz_ = uz
    !$omp end workshare

    call p3dfft_btran_c2r(ux_, ux_r, 'tff')
    call p3dfft_btran_c2r(uy_, uy_r, 'tff')
    call p3dfft_btran_c2r(uz_, uz_r, 'tff')

    !$omp workshare
    u2half_r = 0.5d0*(ux_r**2 + uy_r**2 + uz_r**2)
    !$omp end workshare

    call p3dfft_ftran_r2c(u2half_r, u2half, 'fft'); u2half = u2half/nlx/nly/nlz 

    deallocate(ux_)
    deallocate(uy_)
    deallocate(uz_)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
    deallocate(u2half_r)
  end subroutine get_u2half


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

end module diagnostics


