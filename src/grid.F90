!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Grid setting
!-----------------------------------------------!
module grid
  use p3dfft
  implicit none

  public  init_grid
  public  lx, ly, lz
  public  nlx, nly, nlz, nlxc, nlyc, nlzc
  public  nkx, nky, nkz
  public  nk_local_tot, nl_local_tot
  public  xx, yy, zz
  public  kx, ky, kz, kperp2, kperp2inv, kz2, kperp2_max, kz2_max
  public  ikx, iky, ikz
  public  dlx, dly, dlz, dkx, dky, dkz
  public  k2, k2inv, k2_max
  public  ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
  public  ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
  private read_parameters

  real(r8) :: lx, ly, lz
  integer  :: nlx, nly, nlz, nlxc, nlyc, nlzc
  integer  :: nkx, nky, nkz
  integer  :: nk_local_tot, nl_local_tot
  real(r8), allocatable :: xx(:), yy(:), zz(:)
  real(r8), allocatable :: kx(:), ky(:), kz(:), kz2(:)
  real(r8), allocatable :: k2(:, :, :), k2inv(:, :, :), kperp2(:, :, :), kperp2inv(:, :, :)
  integer , allocatable :: ikx(:), iky(:), ikz(:)
  real(r8) :: dlx, dly, dlz, dkx, dky, dkz
  integer  :: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
  integer  :: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
  real(r8) :: kperp2_max, kz2_max, k2_max

contains


!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Initialization of grid
!-----------------------------------------------!
  subroutine init_grid
    use mp, only: dims, proc0
#ifdef DEBG
    use mp, only: proc_id, nproc
#endif
    use params, only: pi, inputfile
    implicit none
    include 'mpif.h'
    integer istart(3), iend(3), isize(3)
    integer fstart(3), fend(3), fsize(3)
    integer :: i, j, k, ierr

    call read_parameters(inputfile)

    nkx = nlxc
    nky = nlyc/2 + 1
    nkz = nlzc

    allocate(xx(nlx))
    allocate(yy(nly))
    allocate(zz(nlz))
    allocate(kx(nkx))
    allocate(ky(nky))
    allocate(kz(nkz))
    allocate(ikx(nkx))
    allocate(iky(nky))
    allocate(ikz(nkz))
    allocate(kz2(nkz))

    do i = 1, nlx
      xx(i) = lx*dble(1.d0/nlx*(i-1))
    enddo
    do j = 1, nly
      yy(j) = ly*dble(1.d0/nly*(j-1))
    enddo
    do k = 1, nlz
      zz(k) = lz*dble(1.d0/nlz*(k-1))
    enddo

    do j = 1, nky
      iky(j) = j - 1
       ky(j) = 2.d0*pi*dble(iky(j))/ly
    enddo
    do i = 1, nkx
      if (i <= nkx/2 + 1) then
        ikx(i) = i - 1
      else
        ikx(i) = i - nkx - 1
      endif

      kx(i) = 2.d0*pi*dble(ikx(i))/lx
    enddo
    do k = 1, nkz
      if (k <= nkz/2 + 1) then
        ikz(k) = k - 1
      else
        ikz(k) = k - nkz - 1
      endif

      kz(k) = 2.d0*pi*dble(ikz(k))/lz
      kz2(k) = kz(k)**2
    enddo

    dlx = abs(xx(2) - xx(1))
    dly = abs(yy(2) - yy(1))
    dlz = abs(zz(2) - zz(1))
    dkx = abs(kx(2) - kx(1))
    dky = abs(ky(2) - ky(1))
    dkz = abs(kz(2) - kz(1))

    kperp2_max = maxval(abs(kx))**2 + maxval(abs(ky))**2
    kz2_max    = maxval(kz2)
    k2_max     = maxval(abs(kx))**2 + maxval(abs(ky))**2 + maxval(abs(kz))**2

    ! P3DFFT initialization
    call p3dfft_setup (dims, nly, nlz, nlx, MPI_COMM_WORLD, nlyc, nlzc, nlxc, .true.)
    call p3dfft_get_dims(istart, iend, isize, 1)
    call p3dfft_get_dims(fstart, fend, fsize, 2)

    ily_st = istart(1)
    ily_en = iend  (1)
    ilz_st = istart(2)
    ilz_en = iend  (2)
    ilx_st = istart(3)
    ilx_en = iend  (3)

    ikx_st = fstart(1)
    ikx_en = fend  (1)
    ikz_st = fstart(2)
    ikz_en = fend  (2)
    iky_st = fstart(3)
    iky_en = fend  (3)

    nl_local_tot = isize(1)*isize(2)*isize(3)
    nk_local_tot = fsize(1)*fsize(2)*fsize(3)

    if(ikx_en /= nkx) then
      print *, 'Hmm... This is strange because when STRIDE1 is enabled, &
              & ikx_en must be equal to nkx for all processes.'
      stop
    endif

    if(proc0) then
      print *
      print '("P3DFFT initialization is done.")'
      print '("Keep in mind that the array layouts are")'
      print '("  u (*ily, ilz, ilx) for the real space,")'
      print '("  uk(*ikx, ikz, iky) for the Fourier space.")'
      print '("* means that the slot is local in each MPI proc.")'
      print *
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#ifdef DEBG
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(proc0) print *
    if(proc0) print '("Grid decomposition between procs are")'
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call sleep(1)

    do i = 0, nproc
      if(i == proc_id) then
        print '("proc_id = ", i4, ", ily_st:ily_en = ", i4, ":", i4, &
                                  ", ilz_st:ilz_en = ", i4, ":", i4, &
                                  ", ilx_st:ilx_en = ", i4, ":", i4)', &
                proc_id, ily_st, ily_en, ilz_st, ilz_en, ilx_st, ilx_en
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(proc0) print *
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    do i = 0, nproc
      if(i == proc_id) then
        print '("proc_id = ", i4, ", ikx_st:ikx_en = ", i4, ":", i4, &
                                  ", ikz_st:ikz_en = ", i4, ":", i4, &
                                  ", iky_st:iky_en = ", i4, ":", i4)', &
                proc_id, ikx_st, ikx_en, ikz_st, ikz_en, iky_st, iky_en
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if(proc0) print *
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call sleep(1)
#endif

    allocate(k2       (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(k2inv    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(kperp2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(kperp2inv(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          kperp2(i, k, j) = kx(i)**2 + ky(j)**2
          if(kperp2(i, k, j) == 0.d0) then
            kperp2inv(i, k, j) = 0.d0
          else
            kperp2inv(i, k, j) = 1.0d0/kperp2(i, k, j)
          endif

          k2(i, k, j) = kx(i)**2 + ky(j)**2 + kz(k)**2
          if(k2(i, k, j) == 0.d0) then
            k2inv(i, k, j) = 0.d0
          else
            k2inv(i, k, j) = 1.0d0/k2(i, k, j)
          endif
        enddo
      enddo
    enddo

  end subroutine init_grid


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Read inputfile for box parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use mp, only: nproc, proc0, iproc, jproc
    use file, only: get_unused_unit
    use params, only: pi, dealias
    implicit none
    
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /box_parameters/ lx, ly, lz, nlx, nly, nlz

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    lx = 1.d0
    ly = 1.d0
    lz = 1.d0

    nlx = 32
    nly = 32
    nlz = 32
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=box_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading box_parameters failed"
    close(unit)

    lx = 2.d0*pi*lx
    ly = 2.d0*pi*ly
    lz = 2.d0*pi*lz

    if(trim(dealias) == '2/3') then
      nlxc = int(nlx*2.d0/3.d0) ! This is for dealiasing with 2/3 pruned FFT
      nlyc = int(nly*2.d0/3.d0) ! This is for dealiasing with 2/3 pruned FFT
      nlzc = int(nlz*2.d0/3.d0) ! This is for dealiasing with 2/3 pruned FFT
    else
      nlxc = nlx
      nlyc = nly
      nlzc = nlz
    endif

    if((nlx <= 2 .or. nly <= 2 .or. nlz <=2) .and. nproc > 1) then
      if(proc0) then
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, '!                 Only a serial run is allowed for 2D.               !'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *
        stop
      endif
    endif

    if(proc0) then
      print '("nproc  = ", i6, ",  iproc  = ", i6, ",  jproc  = ", i6)', nproc, iproc, jproc
      print '("  nlx  = ", i6, ",    nly  = ", i6, ",    nlz  = ", i6)', nlx , nly , nlz
      print '("  nlxc = ", i6, ",    nlyc = ", i6, ",    nlzc = ", i6)', nlxc, nlyc, nlzc
    endif

  end subroutine read_parameters

end module grid

