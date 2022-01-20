!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Diagnostics for RMHD
!-----------------------------------------------!
module diagnostics_common
  use p3dfft
  implicit none

  public :: init_polar_spectrum_2d, init_polar_spectrum_3d, init_series_modes, finish_diagnostics
  public :: get_polar_spectrum_2d, get_polar_spectrum_3d, cut_2d_r, sum_2d_k, zz_sum
  public :: nkpolar, kpbin
  public :: series_modes_unit

  private

  integer :: nkpolar !    = min(nkx, nky) for 2D bin
                     ! or = min(nkx, nky, nkz) for 3D bin
  real(r8), dimension(:), allocatable :: kpbin
  integer :: series_modes_unit

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of polar spectrum
!!          variables, based on AGK subroutine
!-----------------------------------------------!
  ! bin over (x, y)
  subroutine init_polar_spectrum_2d
    use grid, only: kx, ky
    use params, only: pi
    implicit none
    real(r8) :: dkp
    integer :: i

    dkp = max(kx(2), ky(2))
    nkpolar = int(min(maxval(kx), maxval(ky))/dkp)

    allocate (kpbin(1:nkpolar))

    do i = 1, nkpolar
      kpbin(i) = (i - 1)*dkp
    enddo
  end subroutine init_polar_spectrum_2d

  ! bin over (x, y, z)
  subroutine init_polar_spectrum_3d
    use grid, only: kx, ky, kz
    use params, only: pi
    implicit none
    real(r8) :: dkp
    integer :: i

    dkp = max(kx(2), ky(2), kz(2))
    nkpolar = int(min(maxval(kx), maxval(ky), maxval(kz))/dkp)

    allocate (kpbin(1:nkpolar))

    do i = 1, nkpolar
      kpbin(i) = (i - 1)*dkp
    enddo
  end subroutine init_polar_spectrum_3d


!-----------------------------------------------!
!> @author  YK
!! @date    13 Jul 2021
!! @brief   Initialization series output of modes
!-----------------------------------------------!
  subroutine init_series_modes
    use params, only: series_output
    use file, only: open_output_file
    use params, only: runname
    implicit none

    if(series_output) call open_output_file (series_modes_unit, trim(runname)//'.modes.out')

  end subroutine init_series_modes


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Return kpolar binned and log averaged
!           spectrum, based on AGK subroutine
!-----------------------------------------------!
  ! bin over (x, y)
  subroutine get_polar_spectrum_2d(ee, ebin)
    use params, only: pi
    use mp, only: sum_reduce
    use grid, only: kx, ky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    real(r8), dimension (ikx_st:ikx_en, &
                         ikz_st:ikz_en, &
                         iky_st:iky_en), intent (in)  :: ee   !variable ee by (kx,ky) mode
    real(r8), dimension (1:nkpolar, nkz)  , intent (out) :: ebin 
    integer :: idx, i, j, k

    !Initialize ebin
    ebin = 0.d0

    !Loop through all modes and sum number at each kprp
    do idx = 1, nkpolar - 1
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if (kx(i)**2 + ky(j)**2 >= kpbin(idx)**2 .and. kx(i)**2 + ky(j)**2 < kpbin(idx+1)**2) then
              ebin(idx, k) = ebin(idx, k) + ee(i, k, j)
            endif
          enddo
        enddo
      enddo
    enddo
    call sum_reduce(ebin, 0)

  end subroutine get_polar_spectrum_2d

  ! bin over (x, y, z)
  subroutine get_polar_spectrum_3d(ee, ebin)
    use params, only: pi
    use mp, only: sum_reduce
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: shear, q
    use shearing_box, only: shear_flg, tsc
    implicit none
    real(r8), dimension (ikx_st:ikx_en, &
                         ikz_st:ikz_en, &
                         iky_st:iky_en), intent (in)  :: ee   !variable ee by (kx,ky) mode
    real(r8), dimension (1:nkpolar)  , intent (out) :: ebin 
    real(r8) :: kxt
    integer :: idx, i, j, k

    !Initialize ebin
    ebin = 0.d0

    !Loop through all modes and sum number at each kprp
    do idx = 1, nkpolar - 1
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if(shear) then
              kxt = kx(i) + q*shear_flg*tsc*ky(j)
            else
              kxt = kx(i)
            endif

            if (      kxt**2 + ky(j)**2 + kz(k)**2 >= kpbin(idx)**2 &
                .and. kxt**2 + ky(j)**2 + kz(k)**2 <  kpbin(idx+1)**2) then
              ebin(idx) = ebin(idx) + ee(i, k, j)
            endif
          enddo
        enddo
      enddo
    enddo
    call sum_reduce(ebin, 0)

  end subroutine get_polar_spectrum_3d


!-----------------------------------------------!
!> @author  YK
!! @date    25 Mar 2019
!! @brief   2D cut of real field 
!!                           along z=0, x=0, y=0
!-----------------------------------------------!
  subroutine cut_2d_r(f, fr_z0, fr_x0, fr_y0)
    use mp, only: sum_reduce
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    real(r8), dimension (ily_st:ily_en, &
                         ilz_st:ilz_en, &
                         ilx_st:ilx_en), intent(in) :: f
    real(r8), dimension (nlx, nly), intent(out) :: fr_z0
    real(r8), dimension (nly, nlz), intent(out) :: fr_x0
    real(r8), dimension (nlx, nlz), intent(out) :: fr_y0

    integer :: i, j, k

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          if(k == 1) fr_z0(i, j) = f(j, k, i)
          if(i == 1) fr_x0(j, k) = f(j, k, i)
          if(j == 1) fr_y0(i, k) = f(j, k, i)
        end do
      end do
    end do

    call sum_reduce(fr_z0, 0)
    call sum_reduce(fr_x0, 0)
    call sum_reduce(fr_y0, 0)

  end subroutine cut_2d_r


!-----------------------------------------------!
!> @author  YK
!! @date    8 Sep 2021
!! @brief   sum of spectrum over kz, kx, and ky
!-----------------------------------------------!
  subroutine sum_2d_k(fk, fk_kxy, fk_kyz, fk_kxz)
    use mp, only: sum_reduce
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    real(r8), dimension (ikx_st:ikx_en, &
                         ikz_st:ikz_en, &
                         iky_st:iky_en), intent(in) :: fk
    real(r8), dimension (nkx, nky), intent(out) :: fk_kxy
    real(r8), dimension (nky, nkz), intent(out) :: fk_kyz
    real(r8), dimension (nkx, nkz), intent(out) :: fk_kxz

    integer :: i, j, k

    fk_kxy(:, :) = 0.d0
    fk_kyz(:, :) = 0.d0
    fk_kxz(:, :) = 0.d0

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          fk_kxy(i, j) = fk_kxy(i, j) + fk(i, k, j)
          fk_kyz(j, k) = fk_kyz(j, k) + fk(i, k, j)
          fk_kxz(i, k) = fk_kxz(i, k) + fk(i, k, j)
        end do
      end do
    end do

    call sum_reduce(fk_kxy, 0)
    call sum_reduce(fk_kyz, 0)
    call sum_reduce(fk_kxz, 0)

  end subroutine sum_2d_k


!-----------------------------------------------!
!> @author  YK
!! @date    25 Mar 2019
!! @brief   real field integrated over z
!-----------------------------------------------!
  subroutine zz_sum(f, fr_zsum)
    use mp, only: sum_reduce
    use grid, only: nlx, nly
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: f
    real   (r8), dimension (nlx, nly), intent(out) :: fr_zsum

    complex(r8), allocatable, dimension(:,:,:) :: f_
    real(r8)   , allocatable, dimension(:,:,:) :: fr

    integer :: i, j

    allocate(f_(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); f_ = f
    allocate(fr(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); fr = 0.d0

    ! compute r2c transform
    call p3dfft_ftran_r2c(fr, f_, 'fft')

    do i = ilx_st, ilx_en
      do j = ily_st, ily_en
        fr_zsum(i, j) = sum(fr(j, :, i))
      end do
    end do

    call sum_reduce(fr_zsum, 0)

    deallocate(f_, fr)

  end subroutine zz_sum


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Finalization of diagnostics
!-----------------------------------------------!
  subroutine finish_diagnostics
    use io, only: finish_io
    implicit none

    deallocate(kpbin)

    call finish_io
  end subroutine finish_diagnostics

end module diagnostics_common

