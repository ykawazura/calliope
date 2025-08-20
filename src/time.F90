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
  public  nstep
  public  microsleep
  public  reset_dt
  private read_parameters

  real(r8) :: tt, dt, cfl, increase_dt
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
    implicit none
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /time_parameters/ dt, nstep, cfl, reset_method, increase_dt

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

    dt    = 1.d0
    nstep = 100
    cfl   = 0.5d0
    reset_method = 'multiply'
    increase_dt = huge(1.d0)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=time_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading time_parameters failed"
    close(unit)

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @date    20 Aug 2025
!! @brief   Reset dt
!-----------------------------------------------!
  subroutine reset_dt(dt_cfl, counter)
    use mp, only: proc0
    implicit none
    real(8), intent(in) :: dt_cfl
    integer, intent(inout) :: counter
    real   (r8) :: dt_digit

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

