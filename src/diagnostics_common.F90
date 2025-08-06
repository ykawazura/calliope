!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Diagnostics for RMHD
!-----------------------------------------------!
module diagnostics_common
  use p3dfft
  implicit none

  public :: write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, write_intvl_SF2, SF2_nsample
  public :: write_intvl_nltrans
  public :: series_output, n_series_modes, series_modes
  public :: init_polar_spectrum_2d, init_polar_spectrum_3d, init_series_modes, finish_diagnostics
  public :: get_polar_spectrum_2d, get_polar_spectrum_3d, write_polar_spectrum_2d_in_3d, cut_2d_r, sum_2d_k
  public :: nkpolar, kpbin
  public :: nkpolar_log, kpolar_log_sep, kpbin_log
  public :: series_modes_unit
  public :: read_parameters

  private

  real(r8)             :: write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, write_intvl_SF2, SF2_nsample
  real(r8)             :: write_intvl_nltrans
  integer              :: n_series_modes
  integer, allocatable :: series_modes(:, :)
  logical              :: series_output = .false.

  integer :: nkpolar !    = min(nkx, nky) for 2D bin
                     ! or = min(nkx, nky, nkz) for 3D bin
  real(r8), dimension(:), allocatable :: kpbin

  ! kpbin_log = 0, k0, A*k0, A^2*k0, ..., A^nkpolar_log*k0 ~= k_max
  ! This A is set in inputfile
  integer  :: nkpolar_log
  real(r8) :: kpolar_log_sep
  real(r8), dimension(:), allocatable :: kpbin_log

  integer :: series_modes_unit

contains


!-----------------------------------------------!
!> @author  YK
!! @date    19 Sep 2019
!! @brief   Read inputfile for various parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    integer :: unit
    
    character(len=100), intent(in) :: filename
    integer :: output_modes(100)
    integer  :: i, ierr

    namelist /diagnostics_parameters/ write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, &
                                      write_intvl_SF2, SF2_nsample, write_intvl_nltrans, output_modes, kpolar_log_sep

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    write_intvl         = 0.d0
    write_intvl_2D      = 0.d0
    write_intvl_3D      = 0.d0
    write_intvl_kpar    = 0.d0
    write_intvl_SF2     = 0.d0
    SF2_nsample         = 100000
    write_intvl_nltrans = 0.d0
    output_modes(:)  = 0
    kpolar_log_sep   = 2.d0**(1.d0/3.d0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=diagnostics_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading diagnostics_parameters failed"
    close(unit)

    if (write_intvl         == 0.d0) write_intvl         = -dble(huge(1))
    if (write_intvl_2D      == 0.d0) write_intvl_2D      = -dble(huge(1))
    if (write_intvl_3D      == 0.d0) write_intvl_3D      = -dble(huge(1))
    if (write_intvl_kpar    == 0.d0) write_intvl_kpar    = -dble(huge(1))
    if (write_intvl_SF2     == 0.d0) write_intvl_SF2     = -dble(huge(1))
    if (write_intvl_nltrans == 0.d0) write_intvl_nltrans = -dble(huge(1))

    ! Time series output of modes
    if(count(output_modes /= 0) > 0) then
      if(mod(count(output_modes /= 0), 3) /= 0) then
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print *, '!              Error!              !'
        print *, '!  output modes must be multiple   !'
        print *, '!  of 3; set of (ikx, iky, ikz)    !'
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        stop
      endif

      n_series_modes = count(output_modes /= 0) / 3
      allocate(series_modes(n_series_modes, 3))

      do i = 1, n_series_modes
        series_modes(i, 1) = output_modes((i-1)*3 + 1)
        series_modes(i, 2) = output_modes((i-1)*3 + 2)
        series_modes(i, 3) = output_modes((i-1)*3 + 3)
      enddo

      series_output = .true.
    endif

  end subroutine read_parameters


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
    nkpolar = int(max(maxval(kx), maxval(ky))/dkp)

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
    real(r8) :: dkp, k0, kmax
    integer :: i

    ! set kpbin
    dkp = max(kx(2), ky(2), kz(2))
    nkpolar = int(max(maxval(kx), maxval(ky), maxval(kz))/dkp)

    allocate (kpbin(1:nkpolar))

    do i = 1, nkpolar
      kpbin(i) = (i - 1)*dkp
    enddo

    ! set kpbin_log
    k0   = min(kx(2), ky(2), kz(2))
    kmax = max(maxval(kx), maxval(ky), maxval(kz))
    nkpolar_log = floor( dlog(kmax/k0)/dlog(kpolar_log_sep) ) + 2

    allocate (kpbin_log(1:nkpolar_log))

    kpbin_log(1) = 0.d0
    do i = 1, nkpolar_log - 1
      kpbin_log(i+1) = k0*kpolar_log_sep**(i-1)
    enddo

  end subroutine init_polar_spectrum_3d


!-----------------------------------------------!
!> @author  YK
!! @date    13 Jul 2021
!! @brief   Initialization series output of modes
!-----------------------------------------------!
  subroutine init_series_modes
    use file, only: open_output_file
    use params, only: runname
    use mp, only: proc0
    implicit none

    if(series_output .and. proc0) call open_output_file (series_modes_unit, trim(runname)//'.modes.out')

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

  ! bin over (x, y) leaving z
  !       or (x, z) leaving y
  ! no x leaving because kx will be kxt when shear is on.
  subroutine write_polar_spectrum_2d_in_3d(ee, direction, unit)
    use mp, only: proc0, sum_reduce
    use grid, only: nky, nkz
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: shear, q
    use shearing_box, only: shear_flg, tsc
    implicit none
    real(r8), dimension (ikx_st:ikx_en, &
                         ikz_st:ikz_en, &
                         iky_st:iky_en), intent (in)  :: ee   !variable ee by (kx,ky) mode
    character(len=1), intent(in) :: direction
    integer, intent(in) :: unit

    real(r8), dimension (1:nkpolar) :: ebin 
    real(r8) :: kxt
    integer :: idx, i, j, k

    if(trim(direction) == 'z') then
      do k = 1, nkz
        ebin = 0.d0

        do idx = 1, nkpolar - 1
          if(k >= ikz_st .and. k <= ikz_en) then
            
            do j = iky_st, iky_en
              do i = ikx_st, ikx_en
                if(shear) then
                  kxt = kx(i) + q*shear_flg*tsc*ky(j)
                else
                  kxt = kx(i)
                endif

                if (      kxt**2 + ky(j)**2 >= kpbin(idx)**2 &
                    .and. kxt**2 + ky(j)**2 <  kpbin(idx+1)**2) then
                  ebin(idx) = ebin(idx) + ee(i, k, j)
                endif

              enddo
            enddo

          endif
        enddo
        call sum_reduce(ebin, 0)

        if(proc0) write (unit=unit) ebin

      enddo
    elseif(direction == 'y') then
      do j = 1, nky
        ebin = 0.d0

        do idx = 1, nkpolar - 1
          if(j >= iky_st .and. j <= iky_en) then
            
            do k = ikz_st, ikz_en
              do i = ikx_st, ikx_en
                if(shear) then
                  kxt = kx(i) + q*shear_flg*tsc*ky(j)
                else
                  kxt = kx(i)
                endif

                if (      kxt**2 + kz(k)**2 >= kpbin(idx)**2 &
                    .and. kxt**2 + kz(k)**2 <  kpbin(idx+1)**2) then
                  ebin(idx) = ebin(idx) + ee(i, k, j)
                endif

              enddo
            enddo

          endif
        enddo
        call sum_reduce(ebin, 0)

        if(proc0) write (unit=unit) ebin

      enddo
    endif

    if(proc0) call flush(unit) 

  end subroutine write_polar_spectrum_2d_in_3d


!-----------------------------------------------!
!> @author  YK
!! @date    25 Mar 2019
!! @brief   2D cut of real field 
!!                           along z=0, x=0, y=0
!-----------------------------------------------!
  subroutine cut_2d_r(f, fr_z0, fr_x0, fr_y0)
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    real(r8), dimension (ily_st:ily_en, &
                         ilz_st:ilz_en, &
                         ilx_st:ilx_en), intent(in) :: f
    real(r8), dimension (ilx_st:ilx_en, ily_st:ily_en), intent(out) :: fr_z0
    real(r8), dimension (ily_st:ily_en, ilz_st:ilz_en), intent(out) :: fr_x0
    real(r8), dimension (ilx_st:ilx_en, ilz_st:ilz_en), intent(out) :: fr_y0

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

  end subroutine cut_2d_r


!-----------------------------------------------!
!> @author  YK
!! @date    8 Sep 2021
!! @brief   sum of spectrum over kz, kx, and ky
!-----------------------------------------------!
  subroutine sum_2d_k(fk, fk_kxy, fk_kyz, fk_kxz)
    use grid, only: nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: idx2procid
    use mp, only: sum_reduce
    implicit none
    real(r8), dimension (ikx_st:ikx_en, &
                         ikz_st:ikz_en, &
                         iky_st:iky_en), intent(in) :: fk
    real(r8), dimension (ikx_st:ikx_en, iky_st:iky_en), intent(out) :: fk_kxy
    real(r8), dimension (iky_st:iky_en, ikz_st:ikz_en), intent(out) :: fk_kyz
    real(r8), dimension (ikx_st:ikx_en, ikz_st:ikz_en), intent(out) :: fk_kxz
    real(r8), allocatable :: message(:)

    integer :: i, j, k
    integer :: proc_j1, proc_k1

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

    allocate(message(ikx_st:ikx_en), source=0.d0)
    do j = 1, nky
      proc_k1 = idx2procid(1, j, 'c')
      if(j >= iky_st .and. j <= iky_en) then
        message(:) = fk_kxy(:, j)
      else
        message(:) = 0.d0
      endif

      call sum_reduce(message, proc_k1)
      if(j >= iky_st .and. j <= iky_en) then
        fk_kxy(:, j) = message(:)
      endif
    enddo

    do k = 1, nkz
      proc_j1 = idx2procid(k, 1, 'c')
      if(k >= ikz_st .and. k <= ikz_en) then
        message(:) = fk_kxz(:, k)
      else
        message(:) = 0.d0
      endif

      call sum_reduce(message, proc_j1)
      if(k >= ikz_st .and. k <= ikz_en) then
        fk_kxz(:, k) = message(:)
      endif
    enddo
    deallocate(message)

  end subroutine sum_2d_k


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

