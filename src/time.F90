!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Time setting
!-----------------------------------------------!
module time
  use p3dfft
  implicit none

  public  init_time
  public  tt, dt, cfl, reset_method, increase_dt
  public  track_safety, track_raise, track_digits, t_end
  public  nstep
  public  microsleep
  public  reset_dt
  public  t_end_reached
  private read_parameters

  real(r8) :: tt, dt, cfl, increase_dt
  ! reset_method = 'track': dt = track_safety*dt_cfl (rounded down to track_digits significant digits) when
  ! dt_cfl < dt or dt < track_raise*dt_cfl
  real(r8) :: track_safety, track_raise
  integer  :: track_digits
  ! the main loop stops when tt reaches t_end (within dt/2)
  real(r8) :: t_end
  character(20) :: reset_method
  integer :: nstep = 0


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
!! @date    15 Feb 2021
!! @brief   Initialization of time parameters,
!!          followed by input file reading
!-----------------------------------------------!
  subroutine init_time
    use params, only: inputfile
    use grid, only: lx, ly
    implicit none

    tt  = 0.d0

    call read_parameters(inputfile)

  end subroutine init_time


!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Read inputfile for time parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    use mp, only: proc0
    implicit none
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /time_parameters/ dt, nstep, cfl, reset_method, increase_dt, track_safety, track_raise, &
                               track_digits, t_end

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    dt    = 1.d0
    nstep = 100
    cfl   = 0.5d0
    reset_method = 'multiply'
    increase_dt = huge(1.d0)
    ! track_safety: the one-iteration lambda_svv of MHD_COMP_ISOTH was at most 1.193 times below the converged
    ! value (turb1 on 128^3), and 0.8*1.193 < 1. track_raise: on dt_cfl(t) of that run 0.7 takes 12,251 steps
    ! to 3 tau0 and 0.5 13,388 (tests/MHD_COMP_ISOTH/REPORT_stage6.md)
    ! track_digits: dt_cfl carries a round-off that depends on the MPI decomposition (lambda_svv of
    ! MHD_COMP_ISOTH by 1e-10, relative); 3 digits keep dt the same on any decomposition unless that moves
    ! track_safety*dt_cfl across a rounding point (chance ~1e-7 per change of dt), and lose < 1% of it
    track_safety = 0.8d0
    track_raise  = 0.7d0
    track_digits = 3
    t_end = huge(1.d0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=time_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading time_parameters failed"
    close(unit)

    if(reset_method == 'track' .and. .not. (0.d0 < track_raise .and. track_raise < track_safety &
                                            .and. track_safety <= 1.d0 .and. track_digits >= 1)) then
      if(proc0) print '(" reset_method = ''track'' needs 0 < track_raise < track_safety <= 1 and track_digits >= 1", &
                      &" (track_raise = ", es12.4, ", track_safety = ", es12.4, ", track_digits = ", i0, ")")', &
                      track_raise, track_safety, track_digits
      stop
    endif

  end subroutine read_parameters


!-----------------------------------------------!
!> @brief   True when tt has reached t_end (within dt/2)
!-----------------------------------------------!
  logical function t_end_reached()
    implicit none

    t_end_reached = tt >= t_end - 0.5d0*dt

  end function t_end_reached


!-----------------------------------------------!
!> @author  YK
!! @date    20 Aug 2025
!! @brief   Reset dt
!-----------------------------------------------!
  subroutine reset_dt(dt_cfl, counter, can_increase)
    use mp, only: proc0
    implicit none
    real(8), intent(in) :: dt_cfl
    integer, intent(inout) :: counter
    logical, intent(in), optional :: can_increase
    real   (r8) :: dt_digit, x, q
    character(9) :: change
    logical :: raise

    ! 'track': dt follows dt_cfl, down to track_safety*dt_cfl when dt_cfl < dt, and up to it when
    ! dt < track_raise*dt_cfl and can_increase (default true), rounded down to track_digits significant
    ! digits; increase_dt caps dt
    if (reset_method == 'track') then
      if(.not. (dt_cfl > 0.d0 .and. dt_cfl <= huge(dt_cfl))) return
      raise = .true.
      if(present(can_increase)) raise = can_increase
      if(dt_cfl < dt) then
        change = 'decreased'
      elseif(raise .and. dt < track_raise*dt_cfl .and. dt < increase_dt) then
        change = 'increased'
      else
        return
      endif

      if(proc0) then
        print *
        write (*, '("dt is ", a, " from ", es12.4e3)', advance='no') change, dt
      endif

      x  = track_safety*dt_cfl
      q  = 10.d0**(floor(log10(x)) - track_digits + 1)
      dt = min(aint(x/q)*q, increase_dt)

      counter = 1

      if(proc0) then
        print '("  to ", es12.4e3)', dt
        print *
      endif
      return
    endif

    if(dt_cfl < dt) then
      if(proc0) then
        print *
        write (*, '("dt is decreased from ", es12.4e3)', advance='no') dt
      endif

      dt_digit = (log10(dt_cfl)/abs(log10(dt_cfl)))*ceiling(abs(log10(dt_cfl)))
      ! dt = floor(dt_cfl*10.d0**(-dt_digit))*10.d0**dt_digit

      if (reset_method == 'multiply') then
        dt = 0.5d0*dt
      elseif (reset_method == 'decrement') then
        dt_digit = (log10(dt)/abs(log10(dt)))*ceiling(abs(log10(dt)))

        ! when dt = 0.0**01***
        if (dt*10.d0**(-dt_digit) - 1.0d0 < 1.0d0) then
          dt = 0.9d0*10.d0**dt_digit
        else
          dt = (dt*10.d0**(-dt_digit) - 1.0d0)*10.d0**dt_digit
        endif
      elseif (reset_method == 'stop') then
        if(proc0) then
          print *
          print '("  This run is stopped...")'
          print '("  Bye bye!")'
        endif
        stop
      endif

      counter = 1

      if(proc0) then
        print '("  to ", es12.4e3)', dt
        print *
      endif
    endif
    if(dt_cfl > 0.d0 .and. dt_cfl > increase_dt .and. dt < increase_dt) then
      if(proc0) then
        print *
        write (*, '("dt is increased from ", es12.4e3)', advance='no') dt
      endif

      dt = increase_dt

      counter = 1

      if(proc0) then
        print '("  to ", es12.4e3)', dt
        print *
      endif
    endif

  end subroutine reset_dt


  subroutine microsleep(useconds)
    !! Wrapper for `usleep()` that converts integer to c_int32_t.
    use, intrinsic :: iso_c_binding, only: c_int32_t
    implicit none
    integer :: useconds
    call usleep(int(useconds, kind=c_int32_t))
  end subroutine microsleep

end module time

