!-----------------------------------------------!
!> @author  YK
!! @date    25 Feb 2021
!! @brief   Field setting for MHD_INCOMP
!-----------------------------------------------!
module fields
  use p3dfft
  implicit none

  public :: init_fields, finish_fields, m_to_u, u_to_m
  public :: rho
  public :: mx, my, mz
  public :: ux, uy, uz
  public :: bx, by, bz
  public :: rho_old1
  public :: mx_old1, my_old1, mz_old1
  public :: bx_old1, by_old1, bz_old1
  public :: ux_old1, uy_old1, uz_old1

  private

  complex(r8), allocatable, dimension(:,:,:) :: rho
  complex(r8), allocatable, dimension(:,:,:) :: mx, my, mz
  complex(r8), allocatable, dimension(:,:,:) :: ux, uy, uz
  complex(r8), allocatable, dimension(:,:,:) :: bx, by, bz
  complex(r8), allocatable, dimension(:,:,:) :: rho_old1
  complex(r8), allocatable, dimension(:,:,:) :: mx_old1, my_old1, mz_old1
  complex(r8), allocatable, dimension(:,:,:) :: ux_old1, uy_old1, uz_old1
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
    complex(r8), allocatable, dimension(:,:,:) :: src

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(rho     , source=src)
    allocate(mx      , source=src)
    allocate(my      , source=src)
    allocate(mz      , source=src)
    allocate(ux      , source=src)
    allocate(uy      , source=src)
    allocate(uz      , source=src)
    allocate(bx      , source=src)
    allocate(by      , source=src)
    allocate(bz      , source=src)
    allocate(rho_old1, source=src)
    allocate(mx_old1 , source=src)
    allocate(my_old1 , source=src)
    allocate(mz_old1 , source=src)
    allocate(ux_old1 , source=src)
    allocate(uy_old1 , source=src)
    allocate(uz_old1 , source=src)
    allocate(bx_old1 , source=src)
    allocate(by_old1 , source=src)
    allocate(bz_old1 , source=src)
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
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    real(r8), allocatable, dimension(:,:,:) :: rho_r

    if(proc0) then
      print *, 'Zero initialization'
    endif

    allocate(rho_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en)); rho_r = 0.d0
    rho_r = 1.d0
    call p3dfft_ftran_r2c(rho_r, rho, 'fft'); rho = rho/nlx/nly/nlz 

    mx = 0.d0
    my = 0.d0
    mz = 0.d0
    bx = 0.d0
    by = 0.d0
    bz = 0.d0

    deallocate(rho_r)

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

    rho_old1 = rho
    mx_old1  = mx
    my_old1  = my
    mz_old1  = mz
    bx_old1  = bx
    by_old1  = by
    bz_old1  = bz

    deallocate(rho_r)
    deallocate(bx_r)
    deallocate(by_r)
    deallocate(bz_r)

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
    complex(r8), allocatable, dimension(:,:,:) :: phi, psi, ux, uy, uz
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

    rho_old1 = rho
    mx_old1  = mx
    my_old1  = my
    mz_old1  = mz
    bx_old1  = bx
    by_old1  = by
    bz_old1  = bz

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

    call mpiio_read_one(rho, sizes, subsizes, starts, trim(restart_dir)//'rho.dat' )
    call mpiio_read_one(mx , sizes, subsizes, starts, trim(restart_dir)//'mx.dat')
    call mpiio_read_one(my , sizes, subsizes, starts, trim(restart_dir)//'my.dat')
    call mpiio_read_one(mz , sizes, subsizes, starts, trim(restart_dir)//'mz.dat')
    call mpiio_read_one(bx , sizes, subsizes, starts, trim(restart_dir)//'bx.dat')
    call mpiio_read_one(by , sizes, subsizes, starts, trim(restart_dir)//'by.dat')
    call mpiio_read_one(bz , sizes, subsizes, starts, trim(restart_dir)//'bz.dat')

    rho_old1 = rho
    mx_old1  = mx
    my_old1  = my
    mz_old1  = mz
    bx_old1  = bx
    by_old1  = by
    bz_old1  = bz
  end subroutine restart


!-----------------------------------------------!
!> @author  YK
!! @date    15 Dec 2020
!! @brief   m to u
!-----------------------------------------------!
  subroutine m_to_u(rho, mx, my, mz, ux, uy, uz)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: &
                            rho, mx, my, mz
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

    rho_ = rho
    mx_  = mx
    my_  = my
    mz_  = mz

    call p3dfft_btran_c2r(rho_, rho_r, 'tff')
    call p3dfft_btran_c2r(mx_ , mx_r , 'tff')
    call p3dfft_btran_c2r(my_ , my_r , 'tff')
    call p3dfft_btran_c2r(mz_ , mz_r , 'tff')

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
!! @brief   m to u
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

    deallocate(rho)
    deallocate(mx)
    deallocate(my)
    deallocate(mz)
    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
    deallocate(rho_old1)
    deallocate(mx_old1)
    deallocate(my_old1)
    deallocate(mz_old1)
    deallocate(ux_old1)
    deallocate(uy_old1)
    deallocate(uz_old1)
    deallocate(bx_old1)
    deallocate(by_old1)
    deallocate(bz_old1)

  end subroutine finish_fields

end module fields


