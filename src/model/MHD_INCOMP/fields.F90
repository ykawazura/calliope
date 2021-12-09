!-----------------------------------------------!
!> @author  YK
!! @date    25 Feb 2021
!! @brief   Field setting for MHD_INCOMP
!-----------------------------------------------!
module fields
  use p3dfft
  implicit none

  public :: init_fields, finish_fields
  public :: ux, uy, uz, p
  public :: bx, by, bz
  public :: ux_old1, uy_old1, uz_old1, p_old1
  public :: bx_old1, by_old1, bz_old1

  private

  complex(r8), allocatable, dimension(:,:,:) :: ux, uy, uz, p
  complex(r8), allocatable, dimension(:,:,:) :: bx, by, bz
  complex(r8), allocatable, dimension(:,:,:) :: ux_old1, uy_old1, uz_old1, p_old1
  complex(r8), allocatable, dimension(:,:,:) :: bx_old1, by_old1, bz_old1
  character(100) :: init_type
  real   (r8) :: b0(3)

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

    allocate(ux     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(uy     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(uz     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(bx     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(by     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(bz     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(p      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(ux_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(uy_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(uz_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(bx_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(by_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(bz_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate( p_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))

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

    ux = 0.d0
    uy = 0.d0
    uz = 0.d0
    bx = 0.d0
    by = 0.d0
    bz = 0.d0
     p = 0.d0

  end subroutine init_zero


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Single mode initialization
!-----------------------------------------------!
  subroutine init_single_mode
    use p3dfft
    use mp, only: nproc, proc0
    use grid, only: nlx, nly, nlz, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: zi, q, inputfile
    use file, only: get_unused_unit
    implicit none
    character(20) :: mode_type
    integer :: mode(3)
    real(r8) :: u1(3), b1(3)
    real(r8), allocatable, dimension(:,:,:) :: bz_r
    real(r8) :: kpara, bz0
    complex(r8) :: omega
    integer :: i, j, k

    integer  :: unit, ierr
    namelist /initial_condition_params/ mode_type, mode, b0, b1, u1

    if(proc0) then
      print '("Single mode initialization")'
      print *
    endif

    if(nproc /= 1 .and. proc0) then
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      print *, '!              Error!              !'
      print *, '!  nproc must one for single_mode  !'
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      stop
    endif

    allocate(bz_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); bz_r = 0.d0


    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                read inputfile               v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    mode_type = 'arbitrary'
    mode = (/1, 1, 1/)
    b0  = (/0.d0, 0.d0, 0.d0/)
    b1   = 0.d0
    u1   = 0.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')

    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if(trim(mode_type) /= 'MRI' .and. trim(mode_type) /= 'arbitrary' .and. proc0) then
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      print *, '!              Error!              !'
      print *, '!      mode_type must be either    !'
      print *, '!      "MRI" or "arbitrary".       !'
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      stop
    endif

    ! MRI-eigen mode
    if(trim(mode_type) == 'MRI') then
      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0
      bz = 0.d0

      i = 1; j = 1; k = 2

      bz0  = b0(3)
      bz_r = bz0
      call p3dfft_ftran_r2c(bz_r, bz, 'fft'); bz = bz/nlx/nly/nlz 

      kpara = kz(k)
      omega = dsqrt(4.d0*(kpara*bz0)**2 + (2.d0 - q)**2) - (kpara*bz0)**2 - (2.d0 - q)

      if(real(omega) > 0.d0) then
        omega = zi*dsqrt(real(omega))
        write (*, "('MRI growth rate = ', f5.4)") dimag(omega)
      else
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print *, '!  MRI is stable for these params. !'
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        stop
      endif

      bx(i, k, j) = 1.d-1
      by(i, k, j) = 2.d0*zi*omega/((kpara*bz0)**2 - omega**2)
      ux(i, k, j) = -1.d0/kpara*omega*bx(i, k, j)
      uy(i, k, j) = -1.d0/kpara*(omega*by(i, k, j) + zi*q*bx(i, k, j))

    endif

    ! Arbitrary mode
    if(trim(mode_type) == 'arbitrary') then
      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0
      bz = 0.d0
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if (i == 1 .and. j == 2 .and. k == 2) then
              ux(i, k, j) = 1.d0
              uy(i, k, j) = 1.d0
              uz(i, k, j) = 1.d0
            endif
            if (i == 2 .and. j == 1 .and. k == 2) then
              bx(i, k, j) = 1.d0
              by(i, k, j) = 1.d0
              bz(i, k, j) = 1.d0
            endif
          end do
        end do
      end do
    endif

    ux_old1 = ux
    uy_old1 = uy
    uz_old1 = uz
    bx_old1 = bx
    by_old1 = by
    bz_old1 = bz

    deallocate(bz_r)

  end subroutine init_single_mode


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Random initialization
!-----------------------------------------------!
  subroutine init_random
    use p3dfft
    use grid, only: kx, ky, kz, nlx, nly, nlz, nkx, nky, nkz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, proc_id
    use params, only: zi
    use mp, only: sum_allreduce, sum_reduce
    use params, only: inputfile
    use file, only: get_unused_unit
    use time, only: microsleep
    implicit none
    real(r8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r
    real(r8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r
    integer :: seedsize
    integer, allocatable :: seed(:)
    real(r8) :: rms, kmin(3), kmax(3)
    real(r8) :: u1(3), b1(3)
    integer :: i, j, k

    integer  :: unit, ierr
    namelist /initial_condition_params/ kmin, kmax, b0, b1, u1

    if(proc0) then
      print *, 'Random initialization'
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                read inputfile               v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    kmin = (/0.d0, 0.d0, 0.d0/)
    kmax = (/maxval(kx), maxval(ky), maxval(kz)/)
    b0  = (/0.d0, 0.d0, 0.d0/)
    b1   = 0.d0
    u1   = 0.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')

    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    allocate(ux_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); ux_r = 0.d0
    allocate(uy_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); uy_r = 0.d0
    allocate(uz_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); uz_r = 0.d0
    allocate(bx_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); bx_r = 0.d0
    allocate(by_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); by_r = 0.d0
    allocate(bz_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); bz_r = 0.d0

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
    call random_seed(put=(proc_id+1)*seed(:)) 

    call random_number(ux_r)
    call random_number(uy_r)
    call random_number(uz_r)
    call random_number(bx_r)
    call random_number(by_r)
    call random_number(bz_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    ! compute r2c transform
    call p3dfft_ftran_r2c(ux_r, ux, 'fft'); ux = ux/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uy_r, uy, 'fft'); uy = uy/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uz_r, uz, 'fft'); uz = uz/nlx/nly/nlz 
    call p3dfft_ftran_r2c(bx_r, bx, 'fft'); bx = bx/nlx/nly/nlz 
    call p3dfft_ftran_r2c(by_r, by, 'fft'); by = by/nlx/nly/nlz 
    call p3dfft_ftran_r2c(bz_r, bz, 'fft'); bz = bz/nlx/nly/nlz 

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          ! remove mean fields
          if (kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
            ux(i, k, j) = 0.d0
            uy(i, k, j) = 0.d0
            uz(i, k, j) = 0.d0
            bx(i, k, j) = 0.d0
            by(i, k, j) = 0.d0
            bz(i, k, j) = 0.d0
          endif

          if(nkx > 1)then
            if (kx(i)**2 <= kmin(1)**2) then
              ux(i, k, j) = 0.d0
              uy(i, k, j) = 0.d0
              uz(i, k, j) = 0.d0
              bx(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
              bz(i, k, j) = 0.d0
            endif
            if (kx(i)**2 >= kmax(2)**2) then
              ux(i, k, j) = 0.d0
              uy(i, k, j) = 0.d0
              uz(i, k, j) = 0.d0
              bx(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
              bz(i, k, j) = 0.d0
            endif
          endif

          if(nky > 1)then
            if (ky(j)**2 <= kmin(1)**2) then
              ux(i, k, j) = 0.d0
              uy(i, k, j) = 0.d0
              uz(i, k, j) = 0.d0
              bx(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
              bz(i, k, j) = 0.d0
            endif
            if (ky(j)**2 >= kmax(3)**2) then
              ux(i, k, j) = 0.d0
              uy(i, k, j) = 0.d0
              uz(i, k, j) = 0.d0
              bx(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
              bz(i, k, j) = 0.d0
            endif
          endif

          if(nkz > 1)then
            if (kz(k)**2 <= kmin(2)**2) then
              ux(i, k, j) = 0.d0
              uy(i, k, j) = 0.d0
              uz(i, k, j) = 0.d0
              bx(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
              bz(i, k, j) = 0.d0
            endif

            if (kz(k)**2 >= kmax(3)**2) then
              ux(i, k, j) = 0.d0
              uy(i, k, j) = 0.d0
              uz(i, k, j) = 0.d0
              bx(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
              bz(i, k, j) = 0.d0
            endif
          endif


          if(nky > 1)then
            if (ky(j) == 0.d0) then
              uy(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
            else
              uy(i, k, j) = -(kx(i)*ux(i, k, j) + kz(k)*uz(i, k, j))/ky(j)
              by(i, k, j) = -(kx(i)*bx(i, k, j) + kz(k)*bz(i, k, j))/ky(j)
            endif
          elseif(nkx > 1)then
            if (kx(i) == 0.d0) then
              uy(i, k, j) = 0.d0
              by(i, k, j) = 0.d0
            else
              ux(i, k, j) = -(ky(j)*uy(i, k, j) + kz(k)*uz(i, k, j))/kx(i)
              bx(i, k, j) = -(ky(j)*by(i, k, j) + kz(k)*bz(i, k, j))/kx(i)
            endif
          endif
          ! endif
        end do
      end do
    end do

    ! normalize ux, uy, uz by rms of |u|
    call p3dfft_btran_c2r(ux, ux_r, 'tff')
    call p3dfft_btran_c2r(uy, uy_r, 'tff')
    call p3dfft_btran_c2r(uz, uz_r, 'tff')
    rms = sum(ux_r**2 + uy_r**2 + uz_r**2)
    call sum_allreduce(rms)
    rms = sqrt(rms/nlx/nly/nlz)

    ux_r = u1(1)*ux_r/rms
    uy_r = u1(2)*uy_r/rms
    uz_r = u1(3)*uz_r/rms
    call p3dfft_ftran_r2c(ux_r, ux, 'fft'); ux = ux/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uy_r, uy, 'fft'); uy = uy/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uz_r, uz, 'fft'); uz = uz/nlx/nly/nlz 

    ! normalize bx, by, bz by rms of |b|
    call p3dfft_btran_c2r(bx, bx_r, 'tff')
    call p3dfft_btran_c2r(by, by_r, 'tff')
    call p3dfft_btran_c2r(bz, bz_r, 'tff')
    rms = sum(bx_r**2 + by_r**2 + bz_r**2)
    call sum_allreduce(rms)
    rms = sqrt(rms/nlx/nly/nlz)

    bx_r = b1(1)*bx_r/rms + b0(1)
    by_r = b1(2)*by_r/rms + b0(2)
    bz_r = b1(3)*bz_r/rms + b0(3)

    call p3dfft_ftran_r2c(bx_r, bx, 'fft'); bx = bx/nlx/nly/nlz 
    call p3dfft_ftran_r2c(by_r, by, 'fft'); by = by/nlx/nly/nlz 
    call p3dfft_ftran_r2c(bz_r, bz, 'fft'); bz = bz/nlx/nly/nlz 

    ! check div free of u and b
    call is_div_free('u', ux, uy, uz)
    call is_div_free('b', bx, by, bz)

    deallocate(ux_r )
    deallocate(uy_r )
    deallocate(uz_r )
    deallocate(bx_r )
    deallocate(by_r )
    deallocate(bz_r )
  end subroutine init_random


!-----------------------------------------------!
!> @author  YK
!! @date    18 Feb 2021
!! @brief   2D Orszag Tang problem initialization
!-----------------------------------------------!
  subroutine init_OT2
    use grid, only: lx, ly, xx, yy, nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use mp, only: proc0
    use params, only: zi, pi
    implicit none
    real(r8), allocatable, dimension(:,:,:) :: ux_r, uy_r, bx_r, by_r
    real(r8) :: x0, y0
    integer :: i, j, k

    if(proc0) then
      print *, 'OT2 initialization'
    endif

    x0 = lx/(2.d0*pi)
    y0 = ly/(2.d0*pi)

    allocate(ux_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); ux_r = 0.d0
    allocate(uy_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); uy_r = 0.d0
    allocate(bx_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); bx_r = 0.d0
    allocate(by_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); by_r = 0.d0

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          ux_r(j, k, i) = -sin(yy(j))
          uy_r(j, k, i) =  sin(xx(i))
          bx_r(j, k, i) = -sin(yy(j))
          by_r(j, k, i) =  sin(2.d0*xx(i))
        end do
      end do
    end do

    ! compute r2c transform
    call p3dfft_ftran_r2c(ux_r, ux, 'fft'); ux = ux/nlx/nly/nlz
    call p3dfft_ftran_r2c(uy_r, uy, 'fft'); uy = uy/nlx/nly/nlz
    call p3dfft_ftran_r2c(bx_r, bx, 'fft'); bx = bx/nlx/nly/nlz
    call p3dfft_ftran_r2c(by_r, by, 'fft'); by = by/nlx/nly/nlz

    ux_old1 = ux
    uy_old1 = uy
    uz_old1 = uz
    bx_old1 = bx
    by_old1 = by
    bz_old1 = bz

    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(bx_r)
    deallocate(by_r)

    ! check div free of u and b
    call is_div_free('u', ux, uy, uz)
    call is_div_free('b', bx, by, bz)

  end subroutine init_OT2


!-----------------------------------------------!
!> @author  YK
!! @date    18 Feb 2021
!! @brief   3D Orszag Tang problem initialization
!-----------------------------------------------!
  subroutine init_OT3
    use grid, only: lx, ly, xx, yy, zz, nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use mp, only: proc0
    use params, only: zi, pi
    implicit none
    real(r8), allocatable, dimension(:,:,:) :: ux_r, uy_r, bx_r, by_r
    real(r8) :: x0, y0
    integer :: i, j, k

    if(proc0) then
      print *, 'OT3 initialization'
    endif

    x0 = lx/(2.d0*pi)
    y0 = ly/(2.d0*pi)

    allocate(ux_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); ux_r = 0.d0
    allocate(uy_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); uy_r = 0.d0
    allocate(bx_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); bx_r = 0.d0
    allocate(by_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); by_r = 0.d0

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          ux_r(j, k, i) = -sin(yy(j) + zz(k))
          uy_r(j, k, i) =  sin(xx(i) + zz(k))
          bx_r(j, k, i) = -sin(yy(j) + zz(k))
          by_r(j, k, i) =  sin(2.d0*xx(i) + zz(k))
        end do
      end do
    end do

    ! compute r2c transform
    call p3dfft_ftran_r2c(ux_r, ux, 'fft'); ux = ux/nlx/nly/nlz
    call p3dfft_ftran_r2c(uy_r, uy, 'fft'); uy = uy/nlx/nly/nlz
    call p3dfft_ftran_r2c(bx_r, bx, 'fft'); bx = bx/nlx/nly/nlz
    call p3dfft_ftran_r2c(by_r, by, 'fft'); by = by/nlx/nly/nlz

    ux_old1 = ux
    uy_old1 = uy
    uz_old1 = uz
    bx_old1 = bx
    by_old1 = by
    bz_old1 = bz

    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(bx_r)
    deallocate(by_r)

    ! check div free of u and b
    call is_div_free('u', ux, uy, uz)
    call is_div_free('b', bx, by, bz)

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
    use shearing_box, only: tsc
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
    read (unit=time_unit, fmt="(100es30.21)") tt, tsc
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

    call mpiio_read_one(ux, sizes, subsizes, starts, trim(restart_dir)//'ux.dat')
    call mpiio_read_one(uy, sizes, subsizes, starts, trim(restart_dir)//'uy.dat')
    call mpiio_read_one(uz, sizes, subsizes, starts, trim(restart_dir)//'uz.dat')
    call mpiio_read_one(bx, sizes, subsizes, starts, trim(restart_dir)//'bx.dat')
    call mpiio_read_one(by, sizes, subsizes, starts, trim(restart_dir)//'by.dat')
    call mpiio_read_one(bz, sizes, subsizes, starts, trim(restart_dir)//'bz.dat')
    call mpiio_read_one(p , sizes, subsizes, starts, trim(restart_dir)//'p.dat' )

    ux_old1 = ux
    uy_old1 = uy
    uz_old1 = uz
    bx_old1 = bx
    by_old1 = by
    bz_old1 = bz
     p_old1 =  p
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @date    3 Oct 2020
!! @brief   Check if divergence free is satisfied
!-----------------------------------------------!
  subroutine is_div_free(name, fx, fy, fz)
    use grid, only: kx, ky, kz, k2
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use params, only: zi
    use mp, only: sum_reduce
    implicit none
    character(*) :: name
    integer :: i, j, k
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: fx, fy, fz
    complex(r8), allocatable, dimension(:,:,:) :: abs_div_f
    real(r8) :: abs_div_f_sum, abs_f_sum, abs_k_sum

    allocate(abs_div_f(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); abs_div_f = 0.d0

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          abs_div_f(i,k,j) = abs(kx(i)*fx(i,k,j) + ky(j)*fy(i,k,j) + kz(k)*fz(i,k,j))
        end do
      end do
    end do

    abs_div_f_sum = sum(abs(abs_div_f)); call sum_reduce(abs_div_f_sum, 0)
    abs_f_sum = sum(sqrt(fx**2 + fy**2 + fz**2)); call sum_reduce(abs_f_sum, 0)
    abs_k_sum = sum(sqrt(k2)); call sum_reduce(abs_k_sum, 0)

    if(proc0) write(*, "('<|k.',A2,'|>/(<|k|><|',A2,'|>) = ', es10.3)") trim(name), trim(name), abs_div_f_sum/(abs_f_sum*abs_k_sum)

    deallocate(abs_div_f)

  end subroutine is_div_free


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Finalization of Fields
!-----------------------------------------------!
  subroutine finish_fields
    implicit none

    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
    deallocate( p)
    deallocate(ux_old1)
    deallocate(uy_old1)
    deallocate(uz_old1)
    deallocate(bx_old1)
    deallocate(by_old1)
    deallocate(bz_old1)
    deallocate( p_old1)

  end subroutine finish_fields

end module fields

