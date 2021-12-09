!=======================================================================
! This is based on the 2DECOMP&FFT library
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!=======================================================================

module mpiio
  use p3dfft
  use MPI
  implicit none

  private        ! Make everything private unless declared public

  public :: mpiio_write_one, mpiio_read_one, mpiio_write_var

  interface mpiio_write_one
    module procedure mpiio_write_one_complex
  end interface mpiio_write_one

  interface mpiio_read_one
     module procedure mpiio_read_one_complex
  end interface mpiio_read_one

  interface mpiio_write_var
     module procedure mpiio_write_var_complex
  end interface mpiio_write_var

  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_one_complex(var, sizes, subsizes, starts, filename)
    
    implicit none
    
    complex(r8), dimension(:,:,:), intent(IN) :: var
    integer, dimension(3), intent(IN) :: sizes, subsizes, starts
    character(len=*), intent(IN) :: filename

    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer :: ierror, newtype, fh, data_type

    data_type = complex_type

#include "mpiio_write_one.F90"
    
    return
  end subroutine mpiio_write_one_complex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to read from a file a single 3D array
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_read_one_complex(var, sizes, subsizes, starts, filename)
    
    implicit none
    
    complex(r8), dimension(:,:,:), intent(INOUT) :: var
    integer, dimension(3), intent(IN) :: sizes, subsizes, starts
    character(len=*), intent(IN) :: filename

    integer(kind=MPI_OFFSET_KIND) :: disp
    integer :: ierror, newtype, fh, data_type
    
    data_type = complex_type

#include "mpiio_read_one.F90"
    
    return
  end subroutine mpiio_read_one_complex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the writing
  !  operation to prepare the writing of next chunk of data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_var_complex(fh, disp, sizes, subsizes, starts, var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, dimension(3), intent(IN) :: sizes, subsizes, starts
    complex(r8), dimension(:,:,:), intent(IN) :: var

    integer :: ierror, newtype, data_type, bytes

    data_type = complex_type
    call MPI_TYPE_SIZE(data_type,bytes,ierror)

#include "mpiio_write_var.F90"

    return
  end subroutine mpiio_write_var_complex

end module mpiio
