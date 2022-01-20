!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the 'generic' implementation of the FFT library

module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  use glassman

  !YK: 11 Aug 2020
  use time_stamp, only: put_time_stamp, timer_fft, timer_transpose
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  complex(mytype), allocatable, dimension(:) :: buf, scratch

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    integer :: cbuf_size

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the generic FFT engine *****'
       write(*,*) ' '
    end if

    cbuf_size = max(ph%xsz(1), ph%ysz(2))
    cbuf_size = max(cbuf_size, ph%zsz(3))
    allocate(buf(cbuf_size))
    allocate(scratch(cbuf_size))

    return
  end subroutine init_fft_engine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    deallocate(buf,scratch)

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k
    
    do k=1,decomp%xsz(3)
       do j=1,decomp%xsz(2)
          do i=1,decomp%xsz(1)
             buf(i) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%xsz(1),isign,scratch)
          do i=1,decomp%xsz(1)
             inout(i,j,k) = buf(i)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_x

  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    do k=1,decomp%ysz(3)
       do i=1,decomp%ysz(1)
          do j=1,decomp%ysz(2)
             buf(j) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%ysz(2),isign,scratch)
          do j=1,decomp%ysz(2)
             inout(i,j,k) = buf(j)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_y

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    do j=1,decomp%zsz(2)
       do i=1,decomp%zsz(1)
          do k=1,decomp%zsz(3)
             buf(k) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%zsz(3),isign,scratch)
          do k=1,decomp%zsz(3)
             inout(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_z

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d1

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d1 = size(output,1)

    do k=1,s3
       do j=1,s2
          ! Glassman's FFT is c2c only, 
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do i=1,s1
             buf(i) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,s1,-1,scratch)
          ! note d1 ~ s1/2+1
          ! simply drop the redundant part of the complex output
          do i=1,d1
             output(i,j,k) = buf(i)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_x

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d3

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d3 = size(output,3)

    do j=1,s2
       do i=1,s1
          ! Glassman's FFT is c2c only, 
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do k=1,s3
             buf(k) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,s3,-1,scratch)
          ! note d3 ~ s3/2+1
          ! simply drop the redundant part of the complex output
          do k=1,d3
             output(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_z

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    do k=1,d3
       do j=1,d2
          ! Glassman's FFT is c2c only, 
          ! needing some pre- and post-processing for c2r
          do i=1,d1/2+1
             buf(i) = input(i,j,k)
          end do
          ! expanding to a full-size complex array
          ! For odd N, the storage is:
          !  1, 2, ...... N/2+1   integer division rounded down
          !     N, ...... N/2+2   => a(i) is conjugate of a(N+2-i)
          ! For even N, the storage is:
          !  1, 2, ...... N/2  , N/2+1
          !     N, ...... N/2+2  again a(i) conjugate of a(N+2-i)
          do i=d1/2+2,d1
             buf(i) =  conjg(buf(d1+2-i))
          end do
          call spcfft(buf,d1,1,scratch)
          do i=1,d1
             ! simply drop imaginary part
             output(i,j,k) = real(buf(i), kind=mytype)
          end do
       end do
    end do

    return

  end subroutine c2r_1m_x

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    do j=1,d2
       do i=1,d1
          do k=1,d3/2+1
             buf(k) = input(i,j,k)
          end do
          do k=d3/2+2,d3
             buf(k) =  conjg(buf(d3+2-k))
          end do
          call spcfft(buf,d3,1,scratch)
          do k=1,d3
             output(i,j,k) = real(buf(k), kind=mytype)
          end do
       end do
    end do

    return

  end subroutine c2r_1m_z


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign

#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then
       
       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in,isign,ph)
#else
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       wk1 = in
       call c2c_1m_x(wk1,isign,ph)
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====

       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in,wk2_c2c,ph)
#else
          call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,ph)
       else
#ifdef OVERWRITE
          call c2c_1m_y(in,isign,ph)
#else
          call c2c_1m_y(wk1,isign,ph)
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_c2c,out,ph)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in,out,ph)
#else
          call transpose_y_to_z(wk1,out,ph)
#endif
       end if
       call c2c_1m_z(out,isign,ph)

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in,isign,ph)
#else
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       wk1 = in
       call c2c_1m_z(wk1,isign,ph)
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_z_to_y(in,wk2_c2c,ph)
#else
          call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,ph)
       else  ! out==wk2_c2c if 1D decomposition
#ifdef OVERWRITE
          call transpose_z_to_y(in,out,ph)
#else
          call transpose_z_to_y(wk1,out,ph)
#endif
          call c2c_1m_y(out,isign,ph)
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_c2c,out,ph)
       end if
       call c2c_1m_x(out,isign,ph)
       
    end if

    return
  end subroutine fft_3d_c2c

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D forward FFT - real to complex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_r2c(in_r, out_c)
    
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====

       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call r2c_1m_x(in_r,wk13)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_x_to_y(wk13,wk2_r2c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK

          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2c_1m_y(wk2_r2c,-1,sp)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       else
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2c_1m_y(wk13,-1,sp)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_y_to_z(wk2_r2c,out_c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
       else
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_y_to_z(wk13,out_c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
       end if

       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2c_1m_z(out_c,-1,sp)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
                
    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call r2c_1m_z(in_r,wk13)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_z_to_y(wk13,wk2_r2c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK

          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2c_1m_y(wk2_r2c,-1,sp)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       else  ! out_c==wk2_r2c if 1D decomposition
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_z_to_y(wk13,out_c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK

          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2c_1m_y(out_c,-1,sp)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_y_to_x(wk2_r2c,out_c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
       end if

       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2c_1m_x(out_c,-1,sp)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK

    end if
    
    return
  end subroutine fft_3d_r2c
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D inverse FFT - complex to real
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2r(in_c, out_r)
    
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r

#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2c_1m_z(in_c,1,sp)       
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
#else
       allocate(wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       wk1 = in_c

       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2c_1m_z(wk1,1,sp)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
       call transpose_z_to_y(in_c,wk2_r2c,sp)
       if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
#else
       if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
       call transpose_z_to_y(wk1,wk2_r2c,sp)
       if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
#endif
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2c_1m_y(wk2_r2c,1,sp)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_y_to_x(wk2_r2c,wk13,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK

          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2r_1m_x(wk13,out_r)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       else
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2r_1m_x(wk2_r2c,out_r)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       end if

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2c_1m_x(in_c,1,sp)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
#else
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       wk1 = in_c

       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2c_1m_x(wk1,1,sp)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_x_to_y(in_c,wk2_r2c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
#else
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_x_to_y(wk1,wk2_r2c,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
#endif
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2c_1m_y(wk2_r2c,1,sp)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       else  ! in_c==wk2_r2c if 1D decomposition
#ifdef OVERWRITE
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2c_1m_y(in_c,1,sp)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
#else
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
          call c2c_1m_y(wk1,1,sp)
          if (nrank == 0) call put_time_stamp(timer_fft) ! YK
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_y_to_z(wk2_r2c,wk13,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
       else
#ifdef OVERWRITE
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_y_to_z(in_c,wk13,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
#else
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
          call transpose_y_to_z(wk1,wk13,sp)
          if (nrank == 0) call put_time_stamp(timer_transpose) ! YK
#endif
       end if
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK
       call c2r_1m_z(wk13,out_r)
       if (nrank == 0) call put_time_stamp(timer_fft) ! YK

    end if

    return
  end subroutine fft_3d_c2r

  
end module decomp_2d_fft
