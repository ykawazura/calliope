!=======================================================================
! This is based on the 2DECOMP&FFT library
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!=======================================================================

! This file contain common code to be included by subroutines 
! 'mpiio_read_one_...' in io.f90

    ! Using MPI-IO to write a distributed 3D array into a file
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,data_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
