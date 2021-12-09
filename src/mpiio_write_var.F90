!=======================================================================
! This is based on the 2DECOMP&FFT library
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!=======================================================================

! This file contain common code to be included by subroutines 
! 'write_var_...' in io.f90

  ! Using MPI-IO to write a distributed 3D variable to a file. File 
  ! operations (open/close) need to be done in calling application. This
  ! allows multiple variables to be written to a single file. Together 
  ! with the corresponding read operation, this is the perfect solution
  ! for applications to perform restart/checkpointing.
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for the next write operation
    disp = disp + sizes(1)*sizes(2)*sizes(3)*bytes
    if (data_type == complex_type) then
       disp = disp + sizes(1)*sizes(2)*sizes(3)*bytes
    end if
