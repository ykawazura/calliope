!-----------------------------------------------!
!> @author  YK
!! @date    25 Feb 2021
!! @brief   Field setting for MHD_COMP_ISOTH
!-----------------------------------------------!
module fields
  use p3dfft
  implicit none

  public :: init_fields, finish_fields, m_to_u, u_to_m
  public :: lnrho_to_rho, rho_to_lnrho
  public :: lnrho
  public :: mx, my, mz
  public :: ux, uy, uz
  public :: bx, by, bz
  public :: lnrho_old
  public :: mx_old, my_old, mz_old
  public :: bx_old, by_old, bz_old
  public :: ux_old, uy_old, uz_old
  public :: nfields
  public :: ilnrho
  public :: imx, imy, imz
  public :: ibx, iby, ibz

  private

  complex(r8), allocatable, dimension(:,:,:) :: lnrho ! state variable s = ln(rho); rho = exp(s) > 0 by construction
  complex(r8), allocatable, dimension(:,:,:) :: mx, my, mz
  complex(r8), allocatable, dimension(:,:,:) :: ux, uy, uz
  complex(r8), allocatable, dimension(:,:,:) :: bx, by, bz
  complex(r8), allocatable, dimension(:,:,:) :: lnrho_old
  complex(r8), allocatable, dimension(:,:,:) :: mx_old, my_old, mz_old
  complex(r8), allocatable, dimension(:,:,:) :: ux_old, uy_old, uz_old
  complex(r8), allocatable, dimension(:,:,:) :: bx_old, by_old, bz_old
  character(100) :: init_type
  real   (r8) :: b0(3)

  ! Field index
  integer, parameter :: nfields = 7
  integer, parameter :: ilnrho = 1
  integer, parameter :: imx  = 2, imy = 3, imz = 4
  integer, parameter :: ibx  = 5, iby = 6, ibz = 7

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
    allocate(lnrho  , source=src)
    allocate(mx     , source=src)
    allocate(my     , source=src)
    allocate(mz     , source=src)
    allocate(ux     , source=src)
    allocate(uy     , source=src)
    allocate(uz     , source=src)
    allocate(bx     , source=src)
    allocate(by     , source=src)
    allocate(bz     , source=src)
    allocate(lnrho_old, source=src)
    allocate(mx_old , source=src)
    allocate(my_old , source=src)
    allocate(mz_old , source=src)
    allocate(ux_old , source=src)
    allocate(uy_old , source=src)
    allocate(uz_old , source=src)
    allocate(bx_old , source=src)
    allocate(by_old , source=src)
    allocate(bz_old , source=src)
    deallocate(src)

    call read_parameters(inputfile)

    if(init_type == 'zero') then
      call init_zero
    endif
    if(init_type == 'single_mode') then
      call init_single_mode
    endif
    if(init_type == 'random') then
      call init_random
    endif
    if(init_type == 'sound_wave' .or. init_type == 'simple_wave') then
      call init_wave1d
    endif
    if(init_type == 'riemann') then
      call init_riemann
    endif
    if(init_type == 'profile1d') then
      call init_profile1d
    endif
    if(init_type == 'restart') then
      call restart
    endif

    if(       init_type /= 'zero'       .and. init_type /= 'single_mode' .and. init_type /= 'random' &
        .and. init_type /= 'sound_wave' .and. init_type /= 'simple_wave' .and. init_type /= 'riemann' &
        .and. init_type /= 'profile1d'  .and. init_type /= 'restart') then
      print *, 'Unknown init_type: ', trim(init_type)
      stop
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

    lnrho = 0.d0 ! rho = 1

    mx = 0.d0
    my = 0.d0
    mz = 0.d0
    bx = 0.d0
    by = 0.d0
    bz = 0.d0

    lnrho_old = lnrho
    mx_old    = mx
    my_old    = my
    mz_old    = mz
    bx_old    = bx
    by_old    = by
    bz_old    = bz

  end subroutine init_zero


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Single mode initialization
!-----------------------------------------------!
  subroutine init_single_mode
    use p3dfft
    use mp, only: nproc, proc0
    use grid, only: nlx, nly, nlz, kx, kz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: inputfile, cs2va2
    use file, only: get_unused_unit
    implicit none
    character(20) :: mode_type
    integer :: mode(3)
    real(r8) :: u1(3), b1(3), rho0, rho1
    real(r8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r, rho_r
    real(r8), allocatable, dimension(:,:,:) :: src
    complex(r8), allocatable, dimension(:,:,:) :: rho
    real(r8) :: kprp, kpar, omega, alpha
    integer :: i, j, k

    integer  :: unit, ierr
    namelist /initial_condition_params/ mode_type, mode, b0, b1, u1, rho0, rho1

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


    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                read inputfile               v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    mode_type = 'arbitrary'
    mode = (/1, 1, 1/)
    b0  = (/0.d0, 0.d0, 0.d0/)
    b1   = 0.d0
    u1   = 0.d0
    rho0 = 1.d0
    rho1 = 0.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')

    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    if(trim(mode_type) /= 'Alfven'    .and. &
       trim(mode_type) /= 'slow'      .and. &
       trim(mode_type) /= 'fast'      .and. &
       trim(mode_type) /= 'arbitrary' .and. proc0) then
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      print *, '!              Error!              !'
      print *, '!  mode_type must be either        !'
      print *, '!  "Alfven", "slow", "fast".       !'
      print *, '!  or "arbitrary".                 !'
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      stop
    endif

    allocate(src(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(rho_r, source=src)
    allocate(bx_r , source=src)
    allocate(by_r , source=src)
    allocate(bz_r , source=src)
    deallocate(src)
    allocate(rho(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))

    ! Alfven eigen mode
    if(trim(mode_type) == 'Alfven') then
      rho_r = 1.d0
      bz_r  = 1.d0
      call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 
      call p3dfft_ftran_r2c(bz_r , bz , 'fft'); bz  = bz /nlx/nly/nlz 

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0

      i = 2; j = 1; k = 2

      uy(i, k, j) = 1.d-1
      by(i, k, j) = -uy(i, k, j)

      ! uy(i, nkz, j) = dconjg(uy(i, k, j)) ! reality condition
      ! by(i, nkz, j) = dconjg(by(i, k, j)) ! reality condition
    endif

    ! slow eigen mode
    if(trim(mode_type) == 'slow') then
      kprp  = kx(2)
      kpar  = kz(2)
      omega = dsqrt( &
                      0.5d0*(kprp**2 + kpar**2)* &
                      (1.d0 + cs2va2 - dsqrt((1.d0 + cs2va2)**2 - 4.d0*cs2va2*kpar**2/(kprp**2 + kpar**2))) &
                   )
      alpha = 1.d0 - (kprp**2 + kpar**2)/omega**2

      rho_r = 1.d0
      bz_r  = 1.d0
      call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 
      call p3dfft_ftran_r2c(bz_r , bz , 'fft'); bz  = bz /nlx/nly/nlz 

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0

      i = 2; j = 1; k = 2

      ux (i, k, j) = 1.d-3
      uz (i, k, j) = alpha*kpar/kprp*ux(i, k, j)
      bx (i, k, j) = -kpar/omega*ux(i, k, j)
      bz (i, k, j) =  kprp/omega*ux(i, k, j)
      rho(i, k, j) = alpha*omega/(kprp*cs2va2)*ux(i, k, j)

    endif

    ! fast eigen mode
    if(trim(mode_type) == 'fast') then
      kprp  = kx(2)
      kpar  = kz(2)
      omega = dsqrt( &
                      0.5d0*(kprp**2 + kpar**2)* &
                      (1.d0 + cs2va2 + dsqrt((1.d0 + cs2va2)**2 - 4.d0*cs2va2*kpar**2/(kprp**2 + kpar**2))) &
                   )
      alpha = 1.d0 - (kprp**2 + kpar**2)/omega**2

      rho_r = 1.d0
      bz_r  = 1.d0
      call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 
      call p3dfft_ftran_r2c(bz_r , bz , 'fft'); bz  = bz /nlx/nly/nlz 

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      bx = 0.d0
      by = 0.d0

      i = 2; j = 1; k = 2

      ux (i, k, j) = 1.d-3
      uz (i, k, j) = alpha*kpar/kprp*ux(i, k, j)
      bx (i, k, j) = -kpar/omega*ux(i, k, j)
      bz (i, k, j) =  kprp/omega*ux(i, k, j)
      rho(i, k, j) = alpha*omega/(kprp*dsqrt(cs2va2))*ux(i, k, j)

    endif

    ! Arbitrary mode
    if(trim(mode_type) == 'arbitrary') then
      rho_r = rho0
      bx_r  = b0(1)
      by_r  = b0(2)
      bz_r  = b0(3)
      call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 
      call p3dfft_ftran_r2c(bx_r , bx , 'fft'); bx  = bx /nlx/nly/nlz 
      call p3dfft_ftran_r2c(by_r , by , 'fft'); by  = by /nlx/nly/nlz 
      call p3dfft_ftran_r2c(bz_r , bz , 'fft'); bz  = bz /nlx/nly/nlz 

      ux = 0.d0
      uy = 0.d0
      uz = 0.d0
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if (i == mode(1) .and. j == mode(2) .and. k == mode(3)) then
              rho(i, k, j) = rho1
              ux (i, k, j) = u1(1)
              uy (i, k, j) = u1(2)
              uz (i, k, j) = u1(3)
              bx (i, k, j) = b1(1)
              by (i, k, j) = b1(2)
              bz (i, k, j) = b1(3)
            endif
          end do
        end do
      end do
    endif

    call u_to_m(rho, ux, uy, uz, mx, my, mz)
    call rho_to_lnrho(rho, lnrho)

    lnrho_old = lnrho
    mx_old  = mx
    my_old  = my
    mz_old  = mz
    bx_old  = bx
    by_old  = by
    bz_old  = bz

    deallocate(rho_r)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)
    deallocate(rho)

  end subroutine init_single_mode


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Random initialization
!-----------------------------------------------!
  subroutine init_random
    use p3dfft
    use grid, only: kx, ky, kz, nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, proc_id
    use params, only: zi, cs2va2
    use mp, only: sum_allreduce
    use params, only: inputfile
    use file, only: get_unused_unit
    use time, only: microsleep
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: phi, psi, ux, uy, uz, rho
    real   (r8), allocatable, dimension(:,:,:) :: phi_r, psi_r, rho_r
    real   (r8), allocatable, dimension(:,:,:) :: mx_r, my_r, mz_r
    real   (r8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r
    real   (r8), allocatable, dimension(:,:,:) :: bx_r, by_r, bz_r
    real   (r8), allocatable, dimension(:,:,:) :: src_r
    complex(r8), allocatable, dimension(:,:,:) :: src_c
    integer :: seedsize
    integer, allocatable :: seed(:)
    real(r8) :: b_rms, u_rms, rho_rms, kmin(3), kmax(3)
    real(r8) :: u1(3), b1(3), rho0, rho1
    real(r8) :: smach_rms, amach_rms, beta_rms
    integer :: i, j, k

    integer  :: unit, ierr
    namelist /initial_condition_params/ kmin, kmax, b0, b1, u1, rho0, rho1

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
    rho0 = 1.d0
    rho1 = 0.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')

    read(unit,nml=initial_condition_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading initial_condition failed"
    close(unit)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(phi_r, source=src_r)
    allocate(psi_r, source=src_r)
    allocate(rho_r, source=src_r)
    allocate(mx_r , source=src_r)
    allocate(my_r , source=src_r)
    allocate(mz_r , source=src_r)
    allocate(ux_r , source=src_r)
    allocate(uy_r , source=src_r)
    allocate(uz_r , source=src_r)
    allocate(bx_r , source=src_r)
    allocate(by_r , source=src_r)
    allocate(bz_r , source=src_r)
    deallocate(src_r)

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(phi, source=src_c)
    allocate(psi, source=src_c)
    allocate(ux , source=src_c)
    allocate(uy , source=src_c)
    allocate(uz , source=src_c)
    allocate(rho, source=src_c)
    deallocate(src_c)

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

    call random_number(phi_r)
    call random_number(psi_r)
    call random_number(rho_r)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    ! compute r2c transform
    call p3dfft_ftran_r2c(phi_r, phi, 'fft'); phi = phi/nlx/nly/nlz 
    call p3dfft_ftran_r2c(psi_r, psi, 'fft'); psi = psi/nlx/nly/nlz 
    call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 

    if(nlz <= 2) then
      ! 2D case
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
              phi(i, k, j) = phi(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              psi(i, k, j) = psi(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              rho(i, k, j) = rho(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            endif

            if (kx(i)**2 < kmin(1)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
            if (ky(j)**2 < kmin(2)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif

            if (kx(i)**2 > kmax(1)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
            if (ky(j)**2 > kmax(2)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif

            if (kz(k) /= 0.d0) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
          end do
        end do
      end do
    else
      ! 3D case
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if (kx(i)**2 + ky(j)**2 + kz(k)**2 /= 0.d0) then
              phi(i, k, j) = phi(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              psi(i, k, j) = psi(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
              rho(i, k, j) = rho(i, k, j)/dsqrt(kx(i)**2 + ky(j)**2 + kz(k)**2)
            endif

            if (kx(i)**2 < kmin(1)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
            if (ky(j)**2 < kmin(2)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
            if (kz(k)**2 < kmin(3)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif

            if (kx(i)**2 > kmax(1)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
            if (ky(j)**2 > kmax(2)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
            if (kz(k)**2 > kmax(3)**2) then
              phi(i, k, j) = 0.d0
              psi(i, k, j) = 0.d0
              rho(i, k, j) = 0.d0
            endif
          end do
        end do
      end do
    endif

    ! u = z x grad(phi)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          ux(i, k, j) = -zi*ky(j)*phi(i, k, j)
          uy(i, k, j) =  zi*kx(i)*phi(i, k, j)
          uz(i, k, j) =  0.d0
          bx(i, k, j) = -zi*ky(j)*psi(i, k, j)
          by(i, k, j) =  zi*kx(i)*psi(i, k, j)
          bz(i, k, j) =  0.d0
        enddo
      enddo
    enddo

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                   set rho                   v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! normalize rho by rms of |rho|
    call p3dfft_btran_c2r(rho, rho_r, 'tff')
    rho_rms = sum(rho_r**2)
    call sum_allreduce(rho_rms)
    rho_rms = sqrt(rho_rms/nlx/nly/nlz)

    rho_r = rho1*rho_r + rho0

    ! recompute rho_rms after adding mean field
    rho_rms = sum(rho_r**2)
    call sum_allreduce(rho_rms)
    rho_rms = sqrt(rho_rms/nlx/nly/nlz)

    call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                    set u                    v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! normalize ux, uy, uz by rms of |u|
    call p3dfft_btran_c2r(ux, ux_r, 'tff')
    call p3dfft_btran_c2r(uy, uy_r, 'tff')
    call p3dfft_btran_c2r(uz, uz_r, 'tff')
    u_rms = sum(ux_r**2 + uy_r**2 + uz_r**2)
    call sum_allreduce(u_rms)
    u_rms = sqrt(u_rms/nlx/nly/nlz)

    ux_r = u1(1)*ux_r/u_rms
    uy_r = u1(2)*uy_r/u_rms
    uz_r = u1(3)*uz_r/u_rms
    mx_r = rho_r*ux_r
    my_r = rho_r*uy_r
    mz_r = rho_r*uz_r

    ! recompute u_rms
    u_rms = sum(ux_r**2 + uy_r**2 + uz_r**2)
    call sum_allreduce(u_rms)
    u_rms = sqrt(u_rms/nlx/nly/nlz)

    call p3dfft_ftran_r2c(ux_r, ux, 'fft'); ux = ux/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uy_r, uy, 'fft'); uy = uy/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uz_r, uz, 'fft'); uz = uz/nlx/nly/nlz 
    call p3dfft_ftran_r2c(mx_r, mx, 'fft'); mx = mx/nlx/nly/nlz 
    call p3dfft_ftran_r2c(my_r, my, 'fft'); my = my/nlx/nly/nlz 
    call p3dfft_ftran_r2c(mz_r, mz, 'fft'); mz = mz/nlx/nly/nlz 
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                    set b                    v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    ! normalize bx, by, bz by rms of |b|
    call p3dfft_btran_c2r(bx, bx_r, 'tff')
    call p3dfft_btran_c2r(by, by_r, 'tff')
    call p3dfft_btran_c2r(bz, bz_r, 'tff')
    b_rms = sum(bx_r**2 + by_r**2 + bz_r**2)
    call sum_allreduce(b_rms)
    b_rms = sqrt(b_rms/nlx/nly/nlz)

    bx_r = b1(1)*bx_r/b_rms + b0(1)
    by_r = b1(2)*by_r/b_rms + b0(2)
    bz_r = b1(3)*bz_r/b_rms + b0(3)

    ! recompute b_rms after adding mean field
    b_rms = sum(bx_r**2 + by_r**2 + bz_r**2)
    call sum_allreduce(b_rms)
    b_rms = sqrt(b_rms/nlx/nly/nlz)

    call p3dfft_ftran_r2c(bx_r, bx, 'fft'); bx = bx/nlx/nly/nlz 
    call p3dfft_ftran_r2c(by_r, by, 'fft'); by = by/nlx/nly/nlz 
    call p3dfft_ftran_r2c(bz_r, bz, 'fft'); bz = bz/nlx/nly/nlz 
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    ! check div free of u and b
    call is_div_free('b', bx, by, bz)

    smach_rms  = u_rms/sqrt(cs2va2)
    amach_rms  = sqrt(rho_rms)*u_rms/b_rms
    beta_rms   = 2.d0*cs2va2*rho_rms/b_rms**2
    if(proc0) write(*, "('rms of (1) sound mach num = ', es10.3, ',  (2) Alfven mach num = ', es10.3, ',  (3) beta = ', es10.3)") &
                & smach_rms, amach_rms, beta_rms

    call rho_to_lnrho(rho, lnrho)

    lnrho_old = lnrho
    mx_old  = mx
    my_old  = my
    mz_old  = mz
    bx_old  = bx
    by_old  = by
    bz_old  = bz

    deallocate(phi_r)
    deallocate(phi  )
    deallocate(psi_r)
    deallocate(psi  )
    deallocate(rho_r)
    deallocate(mx_r)
    deallocate(my_r)
    deallocate(mz_r)
    deallocate(ux  )
    deallocate(uy  )
    deallocate(uz  )
    deallocate(rho )
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)
  end subroutine init_random


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

    call mpiio_read_one(lnrho, sizes, subsizes, starts, trim(restart_dir)//'lnrho.dat')
    call mpiio_read_one(mx , sizes, subsizes, starts, trim(restart_dir)//'mx.dat')
    call mpiio_read_one(my , sizes, subsizes, starts, trim(restart_dir)//'my.dat')
    call mpiio_read_one(mz , sizes, subsizes, starts, trim(restart_dir)//'mz.dat')
    call mpiio_read_one(bx , sizes, subsizes, starts, trim(restart_dir)//'bx.dat')
    call mpiio_read_one(by , sizes, subsizes, starts, trim(restart_dir)//'by.dat')
    call mpiio_read_one(bz , sizes, subsizes, starts, trim(restart_dir)//'bz.dat')

    lnrho_old = lnrho
    mx_old  = mx
    my_old  = my
    mz_old  = mz
    bx_old  = bx
    by_old  = by
    bz_old  = bz
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @date    15 Dec 2020
!! @brief   m to u (u = m/exp(lnrho))
!-----------------------------------------------!
  subroutine m_to_u(lnrho, mx, my, mz, ux, uy, uz)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: &
                            lnrho, mx, my, mz
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(out) :: &
                            ux, uy, uz
    real   (r8), allocatable, dimension(:,:,:) :: rho_r
    real   (r8), allocatable, dimension(:,:,:) :: mx_r, my_r, mz_r
    real   (r8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r
    complex(r8), allocatable, dimension(:,:,:) :: rho_
    complex(r8), allocatable, dimension(:,:,:) :: mx_, my_, mz_
    complex(r8), allocatable, dimension(:,:,:) :: ux_, uy_, uz_
    real   (r8), allocatable, dimension(:,:,:) :: src_r
    complex(r8), allocatable, dimension(:,:,:) :: src_c

    allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(rho_r, source=src_r)
    allocate(mx_r , source=src_r)
    allocate(my_r , source=src_r)
    allocate(mz_r , source=src_r)
    allocate(ux_r , source=src_r)
    allocate(uy_r , source=src_r)
    allocate(uz_r , source=src_r)
    deallocate(src_r)

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(rho_, source=src_c)
    allocate(mx_ , source=src_c)
    allocate(my_ , source=src_c)
    allocate(mz_ , source=src_c)
    allocate(ux_ , source=src_c)
    allocate(uy_ , source=src_c)
    allocate(uz_ , source=src_c)
    deallocate(src_c)

    rho_ = lnrho
    mx_  = mx
    my_  = my
    mz_  = mz

    call p3dfft_btran_c2r(rho_, rho_r, 'tff')
    call p3dfft_btran_c2r(mx_ , mx_r , 'tff')
    call p3dfft_btran_c2r(my_ , my_r , 'tff')
    call p3dfft_btran_c2r(mz_ , mz_r , 'tff')

    rho_r = exp(rho_r)
    ux_r = mx_r/rho_r
    uy_r = my_r/rho_r
    uz_r = mz_r/rho_r

    call p3dfft_ftran_r2c(ux_r, ux, 'fft'); ux = ux/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uy_r, uy, 'fft'); uy = uy/nlx/nly/nlz 
    call p3dfft_ftran_r2c(uz_r, uz, 'fft'); uz = uz/nlx/nly/nlz 

    deallocate(rho_r)
    deallocate(mx_r)
    deallocate(my_r)
    deallocate(mz_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
    deallocate(rho_)
    deallocate(mx_)
    deallocate(my_)
    deallocate(mz_)
    deallocate(ux_)
    deallocate(uy_)
    deallocate(uz_)
  end subroutine m_to_u


!-----------------------------------------------!
!> @author  YK
!! @date    15 Dec 2020
!! @brief   u to m (rho is the density, not ln(rho); used only for initial conditions)
!-----------------------------------------------!
  subroutine u_to_m(rho, ux, uy, uz, mx, my, mz)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: &
                            rho, ux, uy, uz
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(out) :: &
                            mx, my, mz
    real   (r8), allocatable, dimension(:,:,:) :: rho_r
    real   (r8), allocatable, dimension(:,:,:) :: ux_r, uy_r, uz_r
    real   (r8), allocatable, dimension(:,:,:) :: mx_r, my_r, mz_r
    complex(r8), allocatable, dimension(:,:,:) :: rho_
    complex(r8), allocatable, dimension(:,:,:) :: ux_, uy_, uz_
    complex(r8), allocatable, dimension(:,:,:) :: mx_, my_, mz_
    real   (r8), allocatable, dimension(:,:,:) :: src_r
    complex(r8), allocatable, dimension(:,:,:) :: src_c

    allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(rho_r, source=src_r)
    allocate(mx_r , source=src_r)
    allocate(my_r , source=src_r)
    allocate(mz_r , source=src_r)
    allocate(ux_r , source=src_r)
    allocate(uy_r , source=src_r)
    allocate(uz_r , source=src_r)
    deallocate(src_r)

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(rho_, source=src_c)
    allocate(mx_ , source=src_c)
    allocate(my_ , source=src_c)
    allocate(mz_ , source=src_c)
    allocate(ux_ , source=src_c)
    allocate(uy_ , source=src_c)
    allocate(uz_ , source=src_c)
    deallocate(src_c)

    rho_ = rho
    ux_  = ux
    uy_  = uy
    uz_  = uz

    call p3dfft_btran_c2r(rho_, rho_r, 'tff')
    call p3dfft_btran_c2r(ux_ , ux_r , 'tff')
    call p3dfft_btran_c2r(uy_ , uy_r , 'tff')
    call p3dfft_btran_c2r(uz_ , uz_r , 'tff')

    mx_r = ux_r*rho_r
    my_r = uy_r*rho_r
    mz_r = uz_r*rho_r

    call p3dfft_ftran_r2c(mx_r, mx, 'fft'); mx = mx/nlx/nly/nlz 
    call p3dfft_ftran_r2c(my_r, my, 'fft'); my = my/nlx/nly/nlz 
    call p3dfft_ftran_r2c(mz_r, mz, 'fft'); mz = mz/nlx/nly/nlz 

    deallocate(rho_r)
    deallocate(mx_r)
    deallocate(my_r)
    deallocate(mz_r)
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
    deallocate(rho_)
    deallocate(mx_)
    deallocate(my_)
    deallocate(mz_)
    deallocate(ux_)
    deallocate(uy_)
    deallocate(uz_)
  end subroutine u_to_m


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   rho to ln(rho) (Fourier to Fourier).
!!          Stops if rho <= 0 at any grid point; no floor is applied.
!-----------------------------------------------!
  subroutine rho_to_lnrho(rho, lnrho)
    use grid, only: nlx, nly, nlz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use mp, only: proc0, min_allreduce
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: rho
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(out) :: lnrho
    real   (r8), allocatable, dimension(:,:,:) :: rho_r
    complex(r8), allocatable, dimension(:,:,:) :: rho_
    real   (r8) :: rho_min_r

    allocate(rho_ (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(rho_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    rho_ = rho
    call p3dfft_btran_c2r(rho_, rho_r, 'tff')

    rho_min_r = minval(rho_r); call min_allreduce(rho_min_r)
    if(rho_min_r <= 0.d0) then
      if(proc0) then
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print *, '!              Error!              !'
        print *, '!  initial rho <= 0 somewhere      !'
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print '("  min(rho) = ", es12.4e3)', rho_min_r
      endif
      stop
    endif

    rho_r = log(rho_r)
    call p3dfft_ftran_r2c(rho_r, lnrho, 'fft'); lnrho = lnrho/nlx/nly/nlz 

    deallocate(rho_)
    deallocate(rho_r)
  end subroutine rho_to_lnrho


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   ln(rho) to rho (Fourier to Fourier)
!-----------------------------------------------!
  subroutine lnrho_to_rho(lnrho, rho)
    use grid, only: nlx, nly, nlz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: lnrho
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(out) :: rho
    real   (r8), allocatable, dimension(:,:,:) :: rho_r
    complex(r8), allocatable, dimension(:,:,:) :: lnrho_

    allocate(lnrho_(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(rho_r (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    lnrho_ = lnrho
    call p3dfft_btran_c2r(lnrho_, rho_r, 'tff')
    rho_r = exp(rho_r)
    call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 

    deallocate(lnrho_)
    deallocate(rho_r)
  end subroutine lnrho_to_rho


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   1D wave along dir ('x', 'y', or 'z') on a uniform background (rho0, b0)
!!          'sound_wave' : standing wave,  ln(rho) = ln(rho0) + amp*cos(k xi), u = 0
!!          'simple_wave': right-going simple wave, ln(rho) = ln(rho0) + amp*sin(k xi),
!!                         u_dir = c_s ln(rho/rho0) (Riemann invariant u - c_s ln(rho) is uniform).
!!                         Exact solution of isothermal Euler (B parallel to dir) until breaking.
!-----------------------------------------------!
  subroutine init_wave1d
    use mp, only: proc0
    use grid, only: nlx, nly, nlz, lx, ly, lz, xx, yy, zz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: inputfile, cs2va2, pi
    use file, only: get_unused_unit
    implicit none
    character(1) :: dir
    integer  :: mode
    real(r8) :: amp, rho0
    real(r8) :: kk, xi, cs
    real(r8), allocatable, dimension(:,:,:) :: s_r, u_r, m_r, b_r
    integer :: i, j, k
    integer :: unit, ierr
    namelist /wave1d_params/ dir, mode, amp, rho0, b0

    dir  = 'x'
    mode = 1
    amp  = 1.d-3
    rho0 = 1.d0
    b0   = (/0.d0, 0.d0, 0.d0/)

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')
    read(unit,nml=wave1d_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading wave1d_params failed"
    close(unit)

    if(proc0) then
      print '("1D wave initialization: ", A, ", dir = ", A, ", mode = ", I4, ", amp = ", es12.4e3)', &
            trim(init_type), dir, mode, amp
    endif
    if(rho0 <= 0.d0 .or. (dir /= 'x' .and. dir /= 'y' .and. dir /= 'z')) then
      if(proc0) print *, 'init_wave1d: rho0 must be positive and dir must be x, y, or z'
      stop
    endif

    allocate(s_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(u_r, source=s_r)
    allocate(m_r, source=s_r)
    allocate(b_r, source=s_r)

    cs = dsqrt(cs2va2)
    if(dir == 'x') kk = 2.d0*pi*mode/lx
    if(dir == 'y') kk = 2.d0*pi*mode/ly
    if(dir == 'z') kk = 2.d0*pi*mode/lz

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          if(dir == 'x') xi = xx(i)
          if(dir == 'y') xi = yy(j)
          if(dir == 'z') xi = zz(k)
          if(init_type == 'sound_wave') then
            s_r(j, k, i) = log(rho0) + amp*cos(kk*xi)
            u_r(j, k, i) = 0.d0
          else
            s_r(j, k, i) = log(rho0) + amp*sin(kk*xi)
            u_r(j, k, i) = cs*amp*sin(kk*xi)
          endif
        end do
      end do
    end do
    m_r = exp(s_r)*u_r

    call p3dfft_ftran_r2c(s_r, lnrho, 'fft'); lnrho = lnrho/nlx/nly/nlz 
    mx = 0.d0; my = 0.d0; mz = 0.d0
    if(dir == 'x') then
      call p3dfft_ftran_r2c(m_r, mx, 'fft'); mx = mx/nlx/nly/nlz 
    endif
    if(dir == 'y') then
      call p3dfft_ftran_r2c(m_r, my, 'fft'); my = my/nlx/nly/nlz 
    endif
    if(dir == 'z') then
      call p3dfft_ftran_r2c(m_r, mz, 'fft'); mz = mz/nlx/nly/nlz 
    endif
    b_r = b0(1); call p3dfft_ftran_r2c(b_r, bx, 'fft'); bx = bx/nlx/nly/nlz 
    b_r = b0(2); call p3dfft_ftran_r2c(b_r, by, 'fft'); by = by/nlx/nly/nlz 
    b_r = b0(3); call p3dfft_ftran_r2c(b_r, bz, 'fft'); bz = bz/nlx/nly/nlz 

    lnrho_old = lnrho
    mx_old    = mx
    my_old    = my
    mz_old    = mz
    bx_old    = bx
    by_old    = by
    bz_old    = bz

    deallocate(s_r)
    deallocate(u_r)
    deallocate(m_r)
    deallocate(b_r)
  end subroutine init_wave1d


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Periodic Riemann problem along dir ('x', 'y', or 'z').
!!          The R state occupies L/4 < xi < 3L/4 and the L state the rest:
!!            H(xi) = [tanh((xi - L/4)/w) - tanh((xi - 3L/4)/w)]/2,  w = width*dl (width <= 0: step, H = 1/2 at the interfaces),
!!            (rho, u, B) = (1 - H)*(L state) + H*(R state).
!!          The interface at L/4 is the Riemann problem (L | R); the one at 3L/4 is its mirror (R | L).
!!          B_dir must be the same on both sides (div B = 0).
!-----------------------------------------------!
  subroutine init_riemann
    use mp, only: proc0
    use grid, only: nlx, nly, nlz, lx, ly, lz, xx, yy, zz, dlx, dly, dlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: inputfile
    use file, only: get_unused_unit
    implicit none
    character(1) :: dir
    real(r8) :: rho_L, rho_R, u_L(3), u_R(3), b_L(3), b_R(3), width
    real(r8) :: ll, ww, xi, hh
    real(r8), allocatable, dimension(:,:,:) :: s_r, den_r, hh_r, f_r ! density is den_r: rho_r is the same name as rho_R
    integer :: i, j, k, idir
    integer :: unit, ierr
    namelist /riemann_params/ dir, rho_L, rho_R, u_L, u_R, b_L, b_R, width

    dir   = 'x'
    rho_L = 1.d0
    rho_R = 1.d0
    u_L   = 0.d0
    u_R   = 0.d0
    b_L   = 0.d0
    b_R   = 0.d0
    width = 1.d0

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')
    read(unit,nml=riemann_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading riemann_params failed"
    close(unit)

    if(dir == 'x') then
      idir = 1; ll = lx; ww = width*dlx
    elseif(dir == 'y') then
      idir = 2; ll = ly; ww = width*dly
    elseif(dir == 'z') then
      idir = 3; ll = lz; ww = width*dlz
    else
      if(proc0) print *, 'init_riemann: dir must be x, y, or z'
      stop
    endif
    if(rho_L <= 0.d0 .or. rho_R <= 0.d0) then
      if(proc0) print *, 'init_riemann: rho_L and rho_R must be positive'
      stop
    endif
    if(b_L(idir) /= b_R(idir)) then
      if(proc0) print *, 'init_riemann: b_L(dir) must be equal to b_R(dir) (div B = 0)'
      stop
    endif

    if(proc0) then
      print '("Riemann initialization: dir = ", A, ", width = ", f6.2, " cells")', dir, width
      print '("  L: rho = ", es12.4e3, ", u = ", 3es12.4e3, ", b = ", 3es12.4e3)', rho_L, u_L, b_L
      print '("  R: rho = ", es12.4e3, ", u = ", 3es12.4e3, ", b = ", 3es12.4e3)', rho_R, u_R, b_R
    endif

    allocate(s_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(den_r, source=s_r)
    allocate(hh_r , source=s_r)
    allocate(f_r  , source=s_r)

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          if(dir == 'x') xi = xx(i)
          if(dir == 'y') xi = yy(j)
          if(dir == 'z') xi = zz(k)
          if(ww > 0.d0) then
            hh = 0.5d0*(tanh((xi - 0.25d0*ll)/ww) - tanh((xi - 0.75d0*ll)/ww))
          else
            hh = 0.d0
            if(xi > 0.25d0*ll .and. xi < 0.75d0*ll) hh = 1.d0
            if(xi == 0.25d0*ll .or. xi == 0.75d0*ll) hh = 0.5d0
          endif
          hh_r(j, k, i) = hh
        end do
      end do
    end do

    den_r = (1.d0 - hh_r)*rho_L + hh_r*rho_R
    s_r   = log(den_r)
    call p3dfft_ftran_r2c(s_r, lnrho, 'fft'); lnrho = lnrho/nlx/nly/nlz 

    f_r = den_r*((1.d0 - hh_r)*u_L(1) + hh_r*u_R(1)); call p3dfft_ftran_r2c(f_r, mx, 'fft'); mx = mx/nlx/nly/nlz 
    f_r = den_r*((1.d0 - hh_r)*u_L(2) + hh_r*u_R(2)); call p3dfft_ftran_r2c(f_r, my, 'fft'); my = my/nlx/nly/nlz 
    f_r = den_r*((1.d0 - hh_r)*u_L(3) + hh_r*u_R(3)); call p3dfft_ftran_r2c(f_r, mz, 'fft'); mz = mz/nlx/nly/nlz 
    f_r =        (1.d0 - hh_r)*b_L(1) + hh_r*b_R(1) ; call p3dfft_ftran_r2c(f_r, bx, 'fft'); bx = bx/nlx/nly/nlz 
    f_r =        (1.d0 - hh_r)*b_L(2) + hh_r*b_R(2) ; call p3dfft_ftran_r2c(f_r, by, 'fft'); by = by/nlx/nly/nlz 
    f_r =        (1.d0 - hh_r)*b_L(3) + hh_r*b_R(3) ; call p3dfft_ftran_r2c(f_r, bz, 'fft'); bz = bz/nlx/nly/nlz 

    lnrho_old = lnrho
    mx_old    = mx
    my_old    = my
    mz_old    = mz
    bx_old    = bx
    by_old    = by
    bz_old    = bz

    deallocate(s_r)
    deallocate(den_r)
    deallocate(hh_r)
    deallocate(f_r)
  end subroutine init_riemann


!-----------------------------------------------!
!> @author  YK
!! @date    2 Nov 2026
!! @brief   1D profile along dir ('x', 'y', or 'z') read from the ASCII table 'table', one row per grid point:
!!            xi, rho, u_x, u_y, u_z, b_x, b_y, b_z
!!          Blank lines and lines starting with '#' are skipped. Every rank reads the whole table and takes its
!!          own range. The rows must be the grid points xi = (l - 1)*L/n in order, rho > 0, and b_dir must be
!!          uniform (div B = 0).
!-----------------------------------------------!
  subroutine init_profile1d
    use mp, only: proc0
    use grid, only: nlx, nly, nlz, lx, ly, lz, xx, yy, zz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: inputfile
    use file, only: get_unused_unit
    implicit none
    character(1)    :: dir
    character(256)  :: table
    character(1024) :: line
    real(r8) :: ll
    real(r8), allocatable :: tab(:,:), xi(:)
    real(r8), allocatable, dimension(:,:,:) :: s_r, den_r, f_r ! density is den_r, as in init_riemann
    integer :: i, j, k, l, n, nrow, idir
    integer :: unit, ierr
    namelist /profile1d_params/ dir, table

    dir   = 'x'
    table = 'profile1d.dat'

    call get_unused_unit (unit)
    open(unit=unit,file=inputfile,status='old')
    read(unit,nml=profile1d_params,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading profile1d_params failed"
    close(unit)

    if(dir == 'x') then
      idir = 1; ll = lx; n = nlx
    elseif(dir == 'y') then
      idir = 2; ll = ly; n = nly
    elseif(dir == 'z') then
      idir = 3; ll = lz; n = nlz
    else
      if(proc0) print *, 'init_profile1d: dir must be x, y, or z'
      stop
    endif
    allocate(xi(n))
    if(idir == 1) xi = xx
    if(idir == 2) xi = yy
    if(idir == 3) xi = zz

    allocate(tab(8, n), source=0.d0)
    call get_unused_unit (unit)
    open(unit=unit,file=trim(table),status='old',action='read',iostat=ierr)
    if(ierr /= 0) then
      if(proc0) print *, 'init_profile1d: cannot open ', trim(table)
      stop
    endif
    nrow = 0
    do
      read(unit,'(A)',iostat=ierr) line
      if(ierr /= 0) exit
      line = adjustl(line)
      if(len_trim(line) == 0 .or. line(1:1) == '#') cycle
      nrow = nrow + 1
      if(nrow > n) exit
      read(line,*,iostat=ierr) tab(:, nrow)
      if(ierr /= 0) then
        if(proc0) print '("init_profile1d: cannot read 8 numbers from data row ", I8, " of ", A)', nrow, trim(table)
        stop
      endif
    end do
    close(unit)

    if(nrow /= n) then
      if(proc0) print '("init_profile1d: ", A, " must have ", I8, " data rows (one per grid point along ", A, ")")', &
                      trim(table), n, dir
      stop
    endif
    if(maxval(abs(tab(1, :) - xi)) > 1.d-10*ll) then
      if(proc0) print '("init_profile1d: the first column of ", A, " is not the grid (max error ", es10.2e3, ")")', &
                      trim(table), maxval(abs(tab(1, :) - xi))
      stop
    endif
    if(minval(tab(2, :)) <= 0.d0) then
      if(proc0) print *, 'init_profile1d: rho must be positive'
      stop
    endif
    if(maxval(tab(5 + idir, :)) /= minval(tab(5 + idir, :))) then
      if(proc0) print *, 'init_profile1d: b_dir must be uniform (div B = 0)'
      stop
    endif

    if(proc0) then
      print '("1D profile initialization: dir = ", A, ", table = ", A)', dir, trim(table)
      print '("  rho in [", es12.4e3, ", ", es12.4e3, "], b_dir = ", es12.4e3)', &
            minval(tab(2, :)), maxval(tab(2, :)), tab(5 + idir, 1)
    endif

    allocate(s_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(den_r, source=s_r)
    allocate(f_r  , source=s_r)

    call fill(2); den_r = f_r
    s_r = log(den_r)
    call p3dfft_ftran_r2c(s_r, lnrho, 'fft'); lnrho = lnrho/nlx/nly/nlz

    call fill(3); f_r = den_r*f_r; call p3dfft_ftran_r2c(f_r, mx, 'fft'); mx = mx/nlx/nly/nlz
    call fill(4); f_r = den_r*f_r; call p3dfft_ftran_r2c(f_r, my, 'fft'); my = my/nlx/nly/nlz
    call fill(5); f_r = den_r*f_r; call p3dfft_ftran_r2c(f_r, mz, 'fft'); mz = mz/nlx/nly/nlz
    call fill(6);                  call p3dfft_ftran_r2c(f_r, bx, 'fft'); bx = bx/nlx/nly/nlz
    call fill(7);                  call p3dfft_ftran_r2c(f_r, by, 'fft'); by = by/nlx/nly/nlz
    call fill(8);                  call p3dfft_ftran_r2c(f_r, bz, 'fft'); bz = bz/nlx/nly/nlz

    lnrho_old = lnrho
    mx_old    = mx
    my_old    = my
    mz_old    = mz
    bx_old    = bx
    by_old    = by
    bz_old    = bz

    deallocate(s_r)
    deallocate(den_r)
    deallocate(f_r)
    deallocate(tab)
    deallocate(xi)

  contains

    !> f_r = column col of the table at the local grid points
    subroutine fill(col)
      integer, intent(in) :: col
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            if(idir == 1) l = i
            if(idir == 2) l = j
            if(idir == 3) l = k
            f_r(j, k, i) = tab(col, l)
          end do
        end do
      end do
    end subroutine fill

  end subroutine init_profile1d


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
    use shock_viscosity, only: report_timers
    implicit none

    call report_timers

    deallocate(lnrho)
    deallocate(mx)
    deallocate(my)
    deallocate(mz)
    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
    deallocate(lnrho_old)
    deallocate(mx_old)
    deallocate(my_old)
    deallocate(mz_old)
    deallocate(ux_old)
    deallocate(uy_old)
    deallocate(uz_old)
    deallocate(bx_old)
    deallocate(by_old)
    deallocate(bz_old)

  end subroutine finish_fields

end module fields


