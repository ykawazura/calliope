module mp
!
! Easier Fortran90 interface to the MPI Message Passing Library.
!
!     (c) Copyright 1991 to 1998 by Michael A. Beer, William D. Dorland, 
!     P. B. Snyder, Q. P. Liu, and Gregory W. Hammett. ALL RIGHTS RESERVED.
!     
! Note: mp_mpi_r8.f90 is a version of mp_mpi.f90 to use when compiling 
! with -r8 (where the default real type is taken to be 8 bytes).  Just 
! replaced all occurances of MPI_REAL with MPI_DOUBLE_PRECISION and 
! MPI_COMPLEX with MPI_DOUBLE_COMPLEX.
!
  use p3dfft
  implicit none
  private
  include 'mpif.h'

  public :: init_mp, finish_mp
  public :: proc_id, nproc, dims, iproc, jproc, proc0
  public :: broadcast, sum_reduce, sum_allreduce
  public :: max_reduce, max_allreduce
  public :: min_reduce, min_allreduce
  public :: send, receive, isend, ireceive
  public :: barrier

  integer :: dim_decomp, dim1
  integer :: proc_id, nproc, dims(2), iproc, jproc
  logical :: proc0

  interface broadcast
     module procedure broadcast_integer 
     module procedure broadcast_integer_array 

     module procedure broadcast_real    
     module procedure broadcast_real_array    
     module procedure broadcast_real_2d_array

     module procedure broadcast_complex 
     module procedure broadcast_complex_array
     module procedure broadcast_complex_2d_array

     module procedure broadcast_logical 
     module procedure broadcast_logical_array 

     module procedure bcastfrom_integer 
     module procedure bcastfrom_integer_array 

     module procedure bcastfrom_real    
     module procedure bcastfrom_real_array    

     module procedure bcastfrom_complex 
     module procedure bcastfrom_complex_array 

     module procedure bcastfrom_logical 
     module procedure bcastfrom_logical_array 

     module procedure broadcast_character
     module procedure broadcast_character_array
     module procedure bcastfrom_character
  end interface

  interface sum_reduce
     module procedure sum_reduce_integer
     module procedure sum_reduce_integer_array
     module procedure sum_reduce_integer_2array

     module procedure sum_reduce_real
     module procedure sum_reduce_real_array
     module procedure sum_reduce_real_2array
     module procedure sum_reduce_real_3array

     module procedure sum_reduce_complex
     module procedure sum_reduce_complex_array
  end interface

  interface sum_allreduce
     module procedure sum_allreduce_integer
     module procedure sum_allreduce_integer_array
     module procedure sum_allreduce_integer_2array

     module procedure sum_allreduce_real
     module procedure sum_allreduce_real_array
     module procedure sum_allreduce_real_2array

     module procedure sum_allreduce_complex
     module procedure sum_allreduce_complex_array
  end interface

  interface max_reduce
     module procedure max_reduce_integer
     module procedure max_reduce_integer_array

     module procedure max_reduce_real
     module procedure max_reduce_real_array
  end interface

  interface max_allreduce
     module procedure max_allreduce_integer
     module procedure max_allreduce_integer_array

     module procedure max_allreduce_real
     module procedure max_allreduce_real_array
  end interface

  interface min_reduce
     module procedure min_reduce_integer
     module procedure min_reduce_integer_array

     module procedure min_reduce_real
     module procedure min_reduce_real_array
  end interface

  interface min_allreduce
     module procedure min_allreduce_integer
     module procedure min_allreduce_integer_array

     module procedure min_allreduce_real
     module procedure min_allreduce_real_array
  end interface

  interface send
     module procedure send_integer
     module procedure send_integer_array

     module procedure send_real
     module procedure send_real_array

     module procedure send_complex
     module procedure send_complex_array

     module procedure send_logical
     module procedure send_logical_array
  end interface

  interface receive
     module procedure receive_integer
     module procedure receive_integer_array

     module procedure receive_real
     module procedure receive_real_array

     module procedure receive_complex
     module procedure receive_complex_array

     module procedure receive_logical
     module procedure receive_logical_array
  end interface

contains

!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Initialization of MPI
!-----------------------------------------------!
  subroutine init_mp
    use params, only: inputfile
    implicit none
    integer :: ierr

    call MPI_INIT (ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

    proc0 = proc_id == 0 

    call read_parameters(inputfile)

    ! nproc is devided into a iproc x jproc stencle
    if(dim_decomp == 1) then
      dims(1) = 1
      dims(2) = nproc
    else if(dim_decomp == 2) then
      if(dim1 == 0) then
        dims(1) = 0
        dims(2) = 0
        call MPI_Dims_create(nproc,2,dims,ierr)
        if(dims(1) > dims(2)) then
          dims(1) = dims(2)
          dims(2) = nproc / dims(1)
        endif
      else
        dims(1) = dim1
        dims(2) = nproc / dims(1)
      endif
    endif

    iproc = dims(1)
    jproc = dims(2)

#ifndef STRIDE1
    if(proc0) then
      print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *, '!                                                                         !'
      print *, '!  This code is compatible only with P3DFFT with STRIDE1 option enabled.  !'
      print *, '!                The execution is now terminated. Sorry.                  !'
      print *, '!                                                                         !'
      print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *
    endif
    call sleep(1)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    stop
#endif

  end subroutine init_mp


!-----------------------------------------------!
!> @author  YK
!! @date    18 Feb 2021
!! @brief   Read inputfile for MPI settings
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /mpi_settings/ dim_decomp, dim1

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    dim_decomp = 2
    dim1 = 0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=mpi_settings,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading mpi_settings failed"
    close(unit)

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Finalization of MPI & P3DFFT
!-----------------------------------------------!
  subroutine finish_mp
    implicit none
    integer :: ierr

    call p3dfft_clean
    call MPI_FINALIZE (ierr)
  end subroutine finish_mp

! ************** broadcasts *****************************

  subroutine broadcast_character (char)
    implicit none
    character(*), intent (in out) :: char
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_character

  subroutine broadcast_character_array(arr)
    use mpi
    implicit none
    character(len=*), dimension(:), intent(in out) :: arr
    integer :: ierror
    integer :: total_len
    total_len = size(arr) * len(arr(1))
    call mpi_bcast(arr, total_len, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_character_array

  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer

  subroutine broadcast_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer_array

  subroutine broadcast_real (x)
    implicit none
    real(r8), intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real

  subroutine broadcast_real_array (x)
    implicit none
    real(r8), dimension (:), intent (in out) :: x
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_array


!***********************************
  subroutine broadcast_real_2d_array (x)
    implicit none
    real(r8),dimension(:,:), intent (in out) :: x
    integer :: ierror
    call mpi_bcast(x, size(x), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_real_2d_array

!*********************************

  subroutine broadcast_complex (z)
    implicit none
    complex(r8), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex

  subroutine broadcast_complex_array (z)
    implicit none
    complex(r8), dimension (:), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_array

  subroutine broadcast_complex_2d_array (z)
    implicit none
    complex(r8), dimension (:,:), intent (in out) :: z
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_complex_2d_array

  subroutine broadcast_logical (f)
    implicit none
    logical, intent (in out) :: f
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical

  subroutine broadcast_logical_array (f)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_logical_array

  subroutine bcastfrom_logical (f, src)
    implicit none
    logical, intent (in out) :: f
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (f, 1, MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical

  subroutine bcastfrom_logical_array (f, src)
    implicit none
    logical, dimension (:), intent (in out) :: f
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (f, size(f), MPI_LOGICAL, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_logical_array

  subroutine bcastfrom_character (c, src)
    implicit none
    character(*), intent (in out) :: c
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (c, len(c), MPI_CHARACTER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_character

  subroutine bcastfrom_integer (i, src)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer

  subroutine bcastfrom_integer_array (i, src)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (i, size(i), MPI_INTEGER, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_integer_array

  subroutine bcastfrom_real (x, src)
    implicit none
    real(r8), intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, 1, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real

  subroutine bcastfrom_real_array (x, src)
    implicit none
    real(r8), dimension (:), intent (in out) :: x
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (x, size(x), MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_real_array

  subroutine bcastfrom_complex (z, src)
    implicit none
    complex(r8), intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, 1, MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex

  subroutine bcastfrom_complex_array (z, src)
    implicit none
    complex(r8), dimension (:), intent (in out) :: z
    integer, intent (in) :: src
    integer :: ierror
    call mpi_bcast (z, size(z), MPI_DOUBLE_COMPLEX, src, MPI_COMM_WORLD, ierror)
  end subroutine bcastfrom_complex_array

! ************** reductions ***********************

  subroutine sum_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer

  subroutine sum_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_integer_array

  subroutine sum_reduce_integer_2array (i, dest)
    use p3dfft
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    integer, intent (in) :: dest
    integer :: ierror
    if(proc_id.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
    else
       call mpi_reduce &
         (i, i, size(i), MPI_INTEGER, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
    endif
  end subroutine sum_reduce_integer_2array

  subroutine sum_reduce_real (a, dest)
    implicit none
    real(r8), intent (in out) :: a
    integer, intent (in) :: dest
    real(r8) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real

  subroutine sum_reduce_real_array (a, dest)
    implicit none
    real(r8), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(r8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_real_array

  subroutine sum_reduce_real_2array (a, dest)
    use p3dfft
    implicit none
    real(r8), dimension (:,:), intent (in out) :: a
    integer, intent (in) :: dest
    integer :: ierror
    if(proc_id.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
    endif
  end subroutine sum_reduce_real_2array

  subroutine sum_reduce_real_3array (a, dest)
    use p3dfft
    implicit none
    real(r8), dimension (:,:,:), intent (in out) :: a
    integer, intent (in) :: dest
    integer :: ierror
    if(proc_id.eq.dest)then
       call mpi_reduce &
         (MPI_IN_PLACE, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
    else
       call mpi_reduce &
         (a, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
    endif
  end subroutine sum_reduce_real_3array

  subroutine sum_reduce_complex (z, dest)
    implicit none
    complex(r8), intent (in out) :: z
    integer, intent (in) :: dest
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex

  subroutine sum_reduce_complex_array (z, dest)
    implicit none
    complex(r8), dimension (:), intent (in out) :: z
    integer, intent (in) :: dest
    complex(r8), dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_reduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, dest, MPI_COMM_WORLD, ierror)
  end subroutine sum_reduce_complex_array

  subroutine sum_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer

  subroutine sum_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_array

  subroutine sum_allreduce_integer_2array (i)
    implicit none
    integer, dimension (:,:), intent (in out) :: i
    integer, dimension (size(i,1), size(i,2)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_integer_2array

  subroutine sum_allreduce_real (a)
    implicit none
    real(r8), intent (in out) :: a
    real(r8) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real

  subroutine sum_allreduce_real_array (a)
    implicit none
    real(r8), dimension (:), intent (in out) :: a
    real(r8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_array

  subroutine sum_allreduce_real_2array (a)
    implicit none
    real(r8), dimension (:,:), intent (in out) :: a
    real(r8), dimension (size(a,1), size(a,2)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_real_2array

  subroutine sum_allreduce_complex (z)
    implicit none
    complex(r8), intent (in out) :: z
    complex :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex

  subroutine sum_allreduce_complex_array (z)
    implicit none
    complex(r8), dimension (:), intent (in out) :: z
    complex(r8), dimension (size(z)) :: z1
    integer :: ierror
    z1 = z
    call mpi_allreduce &
         (z1, z, size(z), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierror)
  end subroutine sum_allreduce_complex_array

  subroutine max_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer

  subroutine max_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_integer_array

  subroutine max_reduce_real (a, dest)
    implicit none
    real(r8), intent (in out) :: a
    integer, intent (in) :: dest
    real(r8) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real

  subroutine max_reduce_real_array (a, dest)
    implicit none
    real(r8), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(r8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, dest, MPI_COMM_WORLD, ierror)
  end subroutine max_reduce_real_array

  subroutine max_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer

  subroutine max_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_integer_array

  subroutine max_allreduce_real (a)
    implicit none
    real(r8), intent (in out) :: a
    real(r8) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real

  subroutine max_allreduce_real_array (a)
    implicit none
    real(r8), dimension (:), intent (in out) :: a
    real(r8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine max_allreduce_real_array

  subroutine min_reduce_integer (i, dest)
    implicit none
    integer, intent (in out) :: i
    integer, intent (in) :: dest
    integer :: i1, ierror
    i1 = i
    call mpi_reduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer

  subroutine min_reduce_integer_array (i, dest)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, intent (in) :: dest
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_reduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_integer_array

  subroutine min_reduce_real (a, dest)
    implicit none
    real(r8), intent (in out) :: a
    integer, intent (in) :: dest
    real(r8) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real

  subroutine min_reduce_real_array (a, dest)
    implicit none
    real(r8), dimension (:), intent (in out) :: a
    integer, intent (in) :: dest
    real(r8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_reduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, dest, MPI_COMM_WORLD, ierror)
  end subroutine min_reduce_real_array

  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer

  subroutine min_allreduce_integer_array (i)
    implicit none
    integer, dimension (:), intent (in out) :: i
    integer, dimension (size(i)) :: i1
    integer :: ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, size(i), MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer_array

  subroutine min_allreduce_real (a)
    implicit none
    real(r8), intent (in out) :: a
    real(r8) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real

  subroutine min_allreduce_real_array (a)
    implicit none
    real(r8), dimension (:), intent (in out) :: a
    real(r8), dimension (size(a)) :: a1
    integer :: ierror
    a1 = a
    call mpi_allreduce &
         (a1, a, size(a), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_real_array

! ********************* barrier **********************

  subroutine barrier
    implicit none
    integer :: ierror
    call mpi_barrier (MPI_COMM_WORLD, ierror)
  end subroutine barrier

! ********************* sends **********************

  subroutine send_integer (i, dest, tag)
    implicit none
    integer, intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, 1, MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer

  subroutine send_integer_array (i, dest, tag)
    implicit none
    integer, dimension (:), intent (in) :: i
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (i, size(i), MPI_INTEGER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_integer_array

  subroutine send_real (a, dest, tag)
    implicit none
    real(r8), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, 1, MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real

  subroutine send_real_array (a, dest, tag)
    implicit none
    real(r8), dimension (:), intent (in) :: a
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (a, size(a), MPI_DOUBLE_PRECISION, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_real_array

  subroutine send_complex (z, dest, tag)
    implicit none
    complex(r8), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex

 subroutine isend (z, dest, tag)
    implicit none
    complex(r8), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_isend (z, 1, MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine isend

  subroutine send_complex_array (z, dest, tag)
    implicit none
    complex(r8), dimension (:), intent (in) :: z
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (z, size(z), MPI_DOUBLE_COMPLEX, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_complex_array

  subroutine send_logical (f, dest, tag)
    implicit none
    logical, intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, 1, MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical

  subroutine send_logical_array (f, dest, tag)
    implicit none
    logical, dimension (:), intent (in) :: f
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send (f, size(f), MPI_LOGICAL, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_logical_array

  subroutine send_character (s, dest, tag)
    implicit none
    character(*), intent (in) :: s
    integer, intent (in) :: dest
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_send &
         (s, len(s), MPI_CHARACTER, dest, tagp, MPI_COMM_WORLD, ierror)
  end subroutine send_character

! ********************* receives  **********************

  subroutine receive_integer (i, src, tag)
    implicit none
    integer, intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, 1, MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer

  subroutine receive_integer_array (i, src, tag)
    implicit none
    integer, dimension (:), intent (out) :: i
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (i, size(i), MPI_INTEGER, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_integer_array

  subroutine receive_real (a, src, tag)
    implicit none
    real(r8), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, 1, MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real

  subroutine receive_real_array (a, src, tag)
    implicit none
    real(r8), dimension (:), intent (out) :: a
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (a, size(a), MPI_DOUBLE_PRECISION, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_real_array

  subroutine receive_complex (z, src, tag)
    implicit none
    complex(r8), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, 1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex

 subroutine ireceive(z, src, tag)
    implicit none
    complex(r8), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_irecv (z, 1, MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine ireceive

  subroutine receive_complex_array (z, src, tag)
    implicit none
    complex(r8), dimension (:), intent (out) :: z
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (z, size(z), MPI_DOUBLE_COMPLEX, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_complex_array

  subroutine receive_logical (f, src, tag)
    implicit none
    logical, intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, 1, MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical

  subroutine receive_logical_array (f, src, tag)
    implicit none
    logical, dimension (:), intent (out) :: f
    integer, intent (in) :: src
    integer, intent (in), optional :: tag
    integer :: ierror
    integer :: tagp
    integer, dimension (MPI_STATUS_SIZE) :: status
    tagp = 0
    if (present(tag)) tagp = tag
    call mpi_recv (f, size(f), MPI_LOGICAL, src, tagp, MPI_COMM_WORLD, &
        status, ierror)
  end subroutine receive_logical_array

end module mp
