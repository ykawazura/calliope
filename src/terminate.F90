!-----------------------------------------------!
!> @author  YK
!! @date    3 Jun 2021
!! @brief   Terminate the simulation when something happened.
!-----------------------------------------------!
module terminate
  use p3dfft
  implicit none

  public  monitor_terminate
  public  max_wall_time
  public  terminate_file_name
  public  terminated

  logical :: monitor_terminate
  real(r8) :: max_wall_time
  character(len=100) :: terminate_file_name
  logical :: terminated = .false.

contains


!-----------------------------------------------!
!> @author  YK
!! @date    3 Jun 2021
!! @brief   Initialization of terminate
!-----------------------------------------------!
  subroutine init_terminate
    use params, only: inputfile
    implicit none

    call read_parameters(inputfile)

  end subroutine init_terminate


!-----------------------------------------------!
!> @author  YK
!! @date    3 Jun 2021
!! @brief   Read inputfile for various parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    use mp, only: proc0
    implicit none
    integer :: unit
    
    character(len=100), intent(in) :: filename
    integer  :: ierr

    namelist /terminate/ max_wall_time, terminate_file_name

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    max_wall_time = 0.d0
    terminate_file_name = ''
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=terminate,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading terminate failed"
    close(unit)

    monitor_terminate = .false.

    if(max_wall_time == 0.d0) then
      max_wall_time = huge(1)
    else
      monitor_terminate = .true.

      if(proc0) then
        write(*, '(" This run will be terminated when the wall time exceeds ", f5.2, " hours.")') max_wall_time
      endif
    endif

    if(trim(terminate_file_name) /= '') then
      monitor_terminate = .true.

      if(proc0) then
        write(*, '(" This run will be terminated when a file named {", 100A)', advance='no') trim(terminate_file_name)
        write(*, '("} is created in the root directory.")')
      endif
    endif

    if(monitor_terminate .and. proc0) write(*, *)

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @date    3 Jun 2021
!! @brief   Initialization of terminate
!-----------------------------------------------!
  subroutine check_terminate
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_total
    use mp, only: broadcast
    use file, only: get_unused_unit
    implicit none
    real(8) :: wall_time ! in hours
    integer :: unit
    logical :: file_exists


    if (proc0) then
      ! First check if the wall time exceeded the set value
      call put_time_stamp(timer_total)
      wall_time = timer_total(1)/60/60.
      call put_time_stamp(timer_total)

      if(wall_time >= max_wall_time) then
        terminated = .true.
        print *
        write(*, '("!-----------------------------------------------------------------------!")')
        write(*, '("!  Wall time has exceeded the set value. Terminating the simulation...  !")')
        write(*, '("!-----------------------------------------------------------------------!")')
      endif

      ! Second check if there is the termination file is created. 
      call get_unused_unit (unit)
      inquire(file=trim(terminate_file_name), exist=file_exists)

      if(file_exists) then
        ! When the termination file is detected, delete it.
        open(unit=unit,file=trim(terminate_file_name),status='old')
        close(unit, status='delete')

        terminated = .true.
        print *
        write(*, '("!-------------------------------------------------------------------------------------------------------!")')
        write(*, '("!  Termination file has been detected. Terminating the simulation and deleting the termination file...  !")')
        write(*, '("!-------------------------------------------------------------------------------------------------------!")')
      endif

    endif
    
    call broadcast(terminated)


  end subroutine check_terminate

end module terminate

