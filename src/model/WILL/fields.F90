!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Field setting for RRMHD
!-----------------------------------------------!
module fields
  use p3dfft
  implicit none

  public :: init_fields, finish_fields
  public :: phi, omg, psi, aaa
  public :: phi_old, omg_old, psi_old, aaa_old  
  public :: nfields
  public :: iomg, ipsi, iaaa

  private

  complex(r8), allocatable, dimension(:,:,:) :: phi, omg, psi, aaa
  complex(r8), allocatable, dimension(:,:,:) :: phi_old, omg_old, psi_old, aaa_old
  character(100) :: init_type

  ! Field index
  integer, parameter :: nfields = 4
  integer, parameter :: iomg = 1, ipsi = 2, iaaa = 3

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialization of fields
!-----------------------------------------------!
  subroutine init_fields
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: inputfile
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(phi    , source=src)
    allocate(omg    , source=src)
    allocate(psi    , source=src)
    allocate(aaa    , source=src)
    allocate(phi_old, source=src)
    allocate(omg_old, source=src)
    allocate(psi_old, source=src)
    allocate(aaa_old, source=src)
    deallocate(src)

    call read_parameters(inputfile)

    if(init_type == 'zero') then
      call init_zero
    endif
    if(init_type == 'single_mode') then
      call init_single_mode
    endif
    if(init_type == 'OT2') then
      call init_OT2
    endif
    if(init_type == 'OT3') then
      call init_OT3
    endif
    if(init_type == 'random') then
      call init_random
    endif
    if(init_type == 'restart') then
      call restart
    endif
    ! if(init_type == 'test') then
      ! call init_test
    ! endif

  end subroutine init_fields


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Read inputfile for initial condition
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /initial_condition/ init_type

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    init_type = 'zero'
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=initial_condition,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Zero initialization
!-----------------------------------------------!
  subroutine init_zero
    use mp, only: proc0
    implicit none

    if(proc0) then
      print *, 'Zero initialization'
    endif

    phi     = 0.d0
    psi     = 0.d0
    aaa     = 0.d0
    phi_old = 0.d0
    psi_old = 0.d0
    aaa_old = 0.d0

  end subroutine init_zero


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Single mode initialization
!-----------------------------------------------!
  subroutine init_single_mode
    use p3dfft
    use mp, only: proc0
    use grid, only: kprp2 
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    integer :: i, j, k

    if(proc0) then
      print '("Single mode initialization")'
      print *
    endif

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          if (i == 1 .and. j == 1 .and. k == 2) then
            phi(i, k, j) = 1.d0
            psi(i, k, j) = 1.d0
            aaa(i, k, j) = 1.d0
          endif
        end do
      end do
    end do

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          omg(i, k, j) = phi(i, k, j)*(-kprp2(i, k, j))
        enddo
      enddo
    enddo

    phi_old = phi
    omg_old = omg
    psi_old = psi
    aaa_old = aaa

  end subroutine init_single_mode


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Random initialization
!-----------------------------------------------!
  subroutine init_random
    use p3dfft
    use grid, only: kx, ky, kz, nlx, nly, nlz, kprp2 
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, proc_id
    use time, only: microsleep
    implicit none
    real   (r8), allocatable, dimension(:,:,:) :: phi_r, psi_r, aaa_r
    complex(r8), allocatable, dimension(:,:,:) :: phi_
    integer :: seedsize
    integer, allocatable :: seed(:)
    integer :: i, j, k

    if(proc0) then
      print '("Random initialization")'
      print *
    endif

    allocate(phi_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(psi_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(aaa_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(phi_ (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v             create random number             v
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call random_seed(size=seedsize)
    allocate(seed(seedsize)) 

    do i = 1, seedsize
      call system_clock(count=seed(i))
      call microsleep(1000)
      call system_clock(count=seed(i))
    end do
    call random_seed(put=proc_id*seed(:)) 

    call random_number(phi_r)
    call random_number(psi_r)
    call random_number(aaa_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    ! compute r2c transform
    call p3dfft_ftran_r2c(phi_r*1d-1, phi, 'fft'); phi = phi/nlx/nly/nlz 
    call p3dfft_ftran_r2c(psi_r*1d-1, psi, 'fft'); psi = psi/nlx/nly/nlz 
    call p3dfft_ftran_r2c(aaa_r*1d-1, aaa, 'fft'); aaa = aaa/nlx/nly/nlz 

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
            phi(i, k, j) = phi(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            psi(i, k, j) = psi(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            aaa(i, k, j) = aaa(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
          endif
          if (kx(i)**2 + ky(j)**2 + kz(k)**2 >= (0.8*maxval(kx))**2) then
            phi(i, k, j) = 0.d0
            psi(i, k, j) = 0.d0
            aaa(i, k, j) = 0.d0
          endif
        end do
      end do
    end do

    phi_ = phi

    ! compute c2r transform
    call p3dfft_btran_c2r(phi_, phi_r, 'tff')
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          omg(i, k, j) = phi(i, k, j)*(-kprp2(i, k, j))
        enddo
      enddo
    enddo

    do k = ikz_st, ikz_en
      if(kz(k) == 0.d0) then
        phi(:,k,:) = 0.d0
        psi(:,k,:) = 0.d0
        omg(:,k,:) = 0.d0
        aaa(:,k,:) = 0.d0
      endif
    enddo

    phi_old = phi
    omg_old = omg
    psi_old = psi
    aaa_old = aaa

    deallocate(phi_r)
    deallocate(psi_r)
    deallocate(aaa_r)
    deallocate(phi_ )
  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @date    18 Feb 2021
!! @brief   2D Orszag Tang problem initialization
!-----------------------------------------------!
  subroutine init_OT2
    use params, only: pi, zi
    use grid, only: lx, ly, xx, yy, kprp2
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    implicit none
    real(r8), allocatable, dimension(:,:,:) :: phi_r, psi_r, aaa_r
    real(r8) :: x0, y0
    integer :: i, j, k

    if(proc0) then
      print *, '2D Orszag Tang problem initialization'
    endif

    x0 = lx/(2.d0*pi)
    y0 = ly/(2.d0*pi)

    allocate(phi_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(psi_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(aaa_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          phi_r(j, k, i) = -cos(xx(i)/x0) - cos(yy(j)/y0)
          psi_r(j, k, i) = -0.5d0*cos(2.d0*(xx(i)/x0)) - cos(yy(j)/y0)
        end do
      end do
    end do

    ! compute r2c transform
    call p3dfft_ftran_r2c(phi_r, phi, 'fft'); phi = phi/nlx/nly/nlz
    call p3dfft_ftran_r2c(psi_r, psi, 'fft'); psi = psi/nlx/nly/nlz
    call p3dfft_ftran_r2c(aaa_r, aaa, 'fft'); aaa = aaa/nlx/nly/nlz
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          omg(i, k, j) = phi(i, k, j)*(-kprp2(i, k, j))
        enddo
      enddo
    enddo

    phi_old = phi
    omg_old = omg
    psi_old = psi
    aaa_old = aaa

    deallocate(phi_r)
    deallocate(psi_r)
  end subroutine init_OT2


!-----------------------------------------------!
!> @author  YK
!! @date    18 Feb 2021
!! @brief   3D Orszag Tang problem initialization
!-----------------------------------------------!
  subroutine init_OT3
    use params, only: pi, zi
    use grid, only: lx, ly, lz, xx, yy, zz, kprp2
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    implicit none
    real(r8), allocatable, dimension(:,:,:) :: phi_r, psi_r, aaa_r
    real(r8) :: x0, y0, z0
    integer :: i, j, k

    if(proc0) then
      print *, '3D Orszag Tang problem initialization'
    endif

    x0 = lx/(2.d0*pi)
    y0 = ly/(2.d0*pi)
    z0 = lz/(2.d0*pi)

    allocate(phi_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(psi_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(aaa_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    phi_r = 0.d0
    psi_r = 0.d0
    aaa_r = 0.d0

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          phi_r(j, k, i) = -(cos(xx(i)/x0 + 1.4d0 + zz(k)/z0) + cos(yy(j)/y0 + 0.5d0 + zz(k)/z0))
          psi_r(j, k, i) = -(0.5d0*cos(2.d0*(xx(i)/x0) + 2.3d0 + zz(k)/z0) + cos(yy(j)/x0 + 4.1d0 + zz(k)/z0))
        end do
      end do
    end do

    ! compute r2c transform
    call p3dfft_ftran_r2c(phi_r, phi, 'fft'); phi = phi/nlx/nly/nlz
    call p3dfft_ftran_r2c(psi_r, psi, 'fft'); psi = psi/nlx/nly/nlz
    call p3dfft_ftran_r2c(aaa_r, psi, 'fft'); aaa = aaa/nlx/nly/nlz
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          omg(i, k, j) = phi(i, k, j)*(-kprp2(i, k, j))
        enddo
      enddo
    enddo

    phi_old = phi
    omg_old = omg
    psi_old = psi
    aaa_old = aaa

    deallocate(phi_r)
    deallocate(psi_r)
    deallocate(aaa_r)
  end subroutine init_OT3


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Restart
!-----------------------------------------------!
  subroutine restart
    use p3dfft
    use mp, only: proc0
    use time, only: tt, dt
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: restart_dir
    use file, only: open_input_file, close_file
    use mpiio, only: mpiio_read_one
    use MPI
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts

    if(proc0) then
      print *, 'Restart'
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                   Read time                 v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    call open_input_file (time_unit, trim(restart_dir)//'time.dat')
    read (unit=time_unit, fmt=*)
    read (unit=time_unit, fmt="(100es30.21)") tt
    call close_file (time_unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v               Read Binary file              v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    sizes(1) = nkx
    sizes(2) = nkz
    sizes(3) = nky
    subsizes(1) = ikx_en - ikx_st + 1
    subsizes(2) = ikz_en - ikz_st + 1
    subsizes(3) = iky_en - iky_st + 1
    starts(1) = ikx_st - 1
    starts(2) = ikz_st - 1
    starts(3) = iky_st - 1

    call mpiio_read_one(phi, sizes, subsizes, starts, trim(restart_dir)//'phi.dat')
    call mpiio_read_one(omg, sizes, subsizes, starts, trim(restart_dir)//'omg.dat')
    call mpiio_read_one(psi, sizes, subsizes, starts, trim(restart_dir)//'psi.dat')
    call mpiio_read_one(aaa, sizes, subsizes, starts, trim(restart_dir)//'aaa.dat')

    phi_old = phi
    omg_old = omg
    psi_old = psi
    aaa_old = aaa
  end subroutine restart


! !-----------------------------------------------!
! !> @author  YK
! !! @date    29 Dec 2018
! !! @brief   Test
! !-----------------------------------------------!
  ! subroutine init_test
    ! use decomp_2d_fft
    ! use params, only: pi, zi
    ! use grid, only: xx, yy
    ! use grid, only: nlx, nly, nlz
    ! use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    ! use mpienv, only: ph, rank0
    ! implicit none
    ! real(mytype), allocatable, dimension(:,:,:) :: phi_r
    ! integer :: i, j, k

    ! if(rank0) then
      ! print *, '2D Orszag Tang problem initialization'
    ! endif

    ! call alloc_x(phi_r, ph, opt_global=.true.)
    ! phi_r = 0._mytype

    ! do k = ilz_st, ilz_en
      ! do j = ily_st, ily_en
        ! do i = ilx_st, ilx_en
          ! ! phi_r(i, j, k) = cos(xx(i) + yy(j))
          ! ! phi_r(i, j, k) = sin(xx(i) + yy(j))
          ! ! phi_r(i, j, k) = -(cos(xx(i) + 1.4_mytype) + cos(yy(j) + 0.5_mytype))
          ! phi_r(i, j, k) = -(cos(10*xx(i)) + cos(yy(j)))
          ! ! phi_r(i, j, k) = -(cos(xx(i)) + cos(yy(j) + 0.5_mytype))
          ! ! phi_r(i, j, k) = -(cos(xx(i) + 1.4_mytype) + cos(yy(j)))
        ! end do
      ! end do
    ! end do

    ! ! compute r2c transform
    ! call decomp_2d_fft_3d(phi_r, phi); phi = phi/nlx/nly/nlz
    ! ! do k = sp%zst(3), sp%zen(3)
      ! ! k = 1
      ! ! do j = sp%zst(2), sp%zen(2)
        ! ! do i = sp%zst(1), sp%zen(1)
          ! ! omg(i, j, k) = phi(i, j, k)*(-kprp2(i, j, k))
          ! ! if(abs(phi(i, j, k)) > 1e-10_mytype) then
            ! ! print *, i, j, kx(i), ky(j), phi(i, j, k)
          ! ! else
            ! ! print *, i, j, kx(i), ky(j), (0.0, 0.0)
          ! ! endif
        ! ! enddo
      ! ! enddo
      ! ! stop
    ! ! enddo

    ! deallocate(phi_r)

  ! end subroutine init_test


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    implicit none

    deallocate(phi)
    deallocate(omg)
    deallocate(psi)
    deallocate(aaa)
    deallocate(phi_old)
    deallocate(omg_old)
    deallocate(psi_old)
    deallocate(aaa_old)

  end subroutine finish_fields

end module fields

