!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Time stamp
!-----------------------------------------------!
module time_stamp
  use p3dfft
  implicit none

  public  put_time_stamp
  public  timer_total, timer_init, timer_diagnostics, timer_diagnostics_SF2
  public  timer_advance, timer_nonlinear_terms 
  public  timer_fft
  public  microsleep

  real(r8) :: timer_total          (2) = 0.d0
  real(r8) :: timer_init           (2) = 0.d0
  real(r8) :: timer_diagnostics    (2) = 0.d0
  real(r8) :: timer_diagnostics_SF2(2) = 0.d0
  real(r8) :: timer_advance        (2) = 0.d0
  real(r8) :: timer_nonlinear_terms(2) = 0.d0
  real(r8) :: timer_fft            (2) = 0.d0


  interface
    subroutine usleep(useconds) bind(c)
      use, intrinsic :: iso_c_binding, only: c_int32_t
      implicit none
      integer(c_int32_t), value :: useconds
    end subroutine
  end interface

contains

!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Counts elapse time between two calls
!-----------------------------------------------!
  subroutine put_time_stamp(targ)
    use MPI
    implicit none
    real(8), intent(in out) :: targ(2) ! tsum and told
    real(8) :: tnew
    real(8), parameter :: small_number=1.d-10

    tnew = mpi_wtime()

    if (targ(2) == 0.d0) then
       if (tnew == 0.d0) tnew=small_number
       targ(2) = tnew
    else
       targ(1) = targ(1) + tnew - targ(2)
       targ(2) = 0.d0
    end if

  end subroutine put_time_stamp


  subroutine microsleep(useconds)
    !! Wrapper for `usleep()` that converts integer to c_int32_t.
    use, intrinsic :: iso_c_binding, only: c_int32_t
    implicit none
    integer :: useconds
    call usleep(int(useconds, kind=c_int32_t))
  end subroutine microsleep

end module time_stamp

