program fft3d

  use p3dfft
  implicit none
  include 'mpif.h'

  integer, parameter :: nx = 256, ny = 256, nz = 256
  integer, parameter :: ndim = 2
  integer, parameter :: nloop = 10
  integer :: ierr, i, j, k

  real   (p3dfft_type), dimension(:,:,:),  allocatable :: u, u2, u_orig
  complex(p3dfft_type), dimension(:,:,:),  allocatable :: uk
  real   (p3dfft_type) :: twopi
  real   (p3dfft_type) :: time

  integer(i8) Ntot
  real(p3dfft_type) factor
  real(r8) Nglob
  integer dims(2), nproc, proc_id
  integer istart(3), iend(3), isize(3)
  integer fstart(3), fend(3), fsize(3)
  integer iproc, jproc, nxc, nyc, nzc

  call MPI_INIT (ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

  twopi=atan(1.0d0)*8.0d0

  if (proc_id == 0) then
    print *
    print *,'P3DFFT test, random input'
    print '(" procs = ", i7, ",  nx = ", i5, ",  ny = ", i5, ",  nz = ", i5, ",  &
    ndim = ", i5, ",  repeat = ", i5)', nproc, nx, ny, nz, ndim, nloop

    if(p3dfft_type == 4) then
      print *,'Single precision version'
    else if(p3dfft_type == 8) then
      print *,'Double precision version'
    endif
  endif

  ! nproc is devided into a iproc x jproc stencle
  if(ndim == 1) then
    dims(1) = 1
    dims(2) = nproc
  else if(ndim == 2) then
    if (proc_id==0) print *, 'Creating proc. grid with mpi_dims_create'
    dims(1) = 0
    dims(2) = 0
    call MPI_Dims_create(nproc,2,dims,ierr)
    if(dims(1) > dims(2)) then
      dims(1) = dims(2)
      dims(2) = nproc / dims(1)
    endif
  endif

  iproc = dims(1)
  jproc = dims(2)

  if(proc_id == 0) then
    print '(" Using processor grid: ", i3, " x ", i3)', iproc, jproc
  endif

  nxc = nx
  nyc = ny
  nzc = nz

  ! Set up work structures for P3DFFT
  call p3dfft_setup (dims, nx, ny, nz, MPI_COMM_WORLD, nxc, nyc, nzc, .true.)
  if(proc_id == 0) print *

  ! Get dimensions for the original array of real numbers, X-pencils
  call p3dfft_get_dims(istart, iend, isize, 1)

  ! Get dimensions for the R2C-forward-transformed array of complex numbers
  !   Z-pencils (depending on how the library was compiled, the first
  !   dimension could be either X or Z)
  !
  call p3dfft_get_dims(fstart, fend, fsize, 2)

  allocate (u     (istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), stat=ierr)
  if(ierr /= 0)  print *,'Error ',ierr,' allocating array u'
  allocate (u_orig(istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), stat=ierr)
  if(ierr /= 0)  print *,'Error ',ierr,' allocating array u_orig'
  allocate (uk    (fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)), stat=ierr)
  if(ierr /= 0)  print *,'Error ',ierr,' allocating array uk'
  allocate (u2    (istart(1):iend(1), istart(2):iend(2), istart(3):iend(3)), stat=ierr)
  if(ierr /= 0)  print *,'Error ',ierr,' allocating array u2'

  !
  ! Initialize the array to be transformed
  !
  call random_number(u)
  u_orig = u

  ! A couple of Warm-up calls (for more accurate timing reporting)
  call p3dfft_ftran_r2c (u, uk, 'fft')
  call p3dfft_ftran_r2c (u, uk, 'fft')

  Ntot = fsize(1)*fsize(2)*fsize(3)
  Nglob = nx * ny
  Nglob = Nglob * nz
  factor = 1.0d0/Nglob

  if(proc_id == 0) time = mpi_wtime()

  ! do i = 0, nproc - 1
    ! if(nrank == i) then
      ! print *, 'nrank = ', nrank
      ! print *, xstart(1), ':', xend(1), ', ', xstart(2), ':', xend(2), ', ', xstart(3), ':', xend(3)
      ! print *, fstart(1), ':', fend(1), ', ', fstart(2), ':', fend(2), ', ', fstart(3), ':', fend(3)
    ! endif
    ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ! enddo

  ! Repeat n times
  do i = 1, nloop
    ! Forward transform
    call p3dfft_ftran_r2c (u , uk, 'fft')
    ! Backward transform
    call p3dfft_btran_c2r (uk, u2, 'tff')
  end do

  if(proc_id == 0) then
    print '(" time per step  = ", es10.3, " seconds")', (mpi_wtime() - time)/nloop
    print *
  endif

  ! Check results
  call p3dfft_ftran_r2c (u , uk, 'fft')
  uk = uk*factor
  call p3dfft_btran_c2r (uk, u2, 'tff')
  call check_res

  ! Free work space
  call p3dfft_clean

  call MPI_FINALIZE (ierr)

  contains

  !=========================================================
  subroutine check_res
  !=========================================================
    implicit none

    real(p3dfft_type) cdiff, ccdiff, prec

    cdiff = 0.0d0
    ccdiff = 0.0d0
    do 20 k = istart(3), iend(3)
      do 20 j = istart(2), iend(2)
        do 20 i = istart(1), iend(1)
          if(cdiff < abs(u_orig(i,j,k) - u2(i,j,k))) then
            cdiff = abs(u_orig(i,j,k) - u2(i,j,k))
          endif
  20   continue
    call MPI_Reduce(cdiff, ccdiff, 1, p3dfft_mpireal, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if(proc_id == 0) then
      if(p3dfft_type == 8) then
        prec = 1e-14
      else
        prec = 1e-5
      endif
      if(ccdiff > prec * Nglob*0.25) then
        print *, 'Results are incorrect'
      else
        print *, 'Results are correct'
      endif
      print '(" max diff  = ", es10.3)', ccdiff
    endif

    return
  end subroutine

end

