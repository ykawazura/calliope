!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Read and set various run parameters
!-----------------------------------------------!
module params_common
  use p3dfft
  implicit none

  public  init_params_common
  public  zi, pi
  public  runname, inputfile
  public  time_step_scheme
  public  save_restart_intvl
  public  restart_dir
  private read_parameters

  !> square root of -1.
  complex(r8), parameter :: zi = ( 0.d0 , 1.d0 )

  !> pi to quad precision
  real(r8), parameter :: pi = &
       3.14159265358979323846264338327950288419716939938d0

  character(len=100) :: runname = 'calliope'
  character(len=100) :: inputfile

  character(len=100) :: time_step_scheme

  real(r8)           :: save_restart_intvl
  character(len=100) :: restart_dir

contains


!-----------------------------------------------!
!> @author  YK
!! @date    6 Dec 2018
!! @brief   Initialization of run parameters,
!!          followed by input file reading
!-----------------------------------------------!
  subroutine init_params_common
    implicit none

    character(len=100) :: arg

    ! set runname
    call getarg(1,arg)
    if (len(trim(arg)) /= 0) runname = arg

    ! read inputfile = $(runname).in
    inputfile = trim(runname)//".in"

    call read_parameters(inputfile)

  end subroutine init_params_common


!-----------------------------------------------!
!> @author  YK
!! @date    19 Sep 2019
!! @brief   Read inputfile for various parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    integer :: unit
    
    character(len=100), intent(in) :: filename
    integer :: output_modes(100)
    integer :: ierr

    namelist /scheme_parameters/ time_step_scheme
    namelist /restart_parameters/ save_restart_intvl, restart_dir

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    time_step_scheme = 'eSSPIFRK3'
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=scheme_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading scheme_parameters failed"
    close(unit)

    if(time_step_scheme /= 'eSSPIFRK3' .and. time_step_scheme /= 'gear3') then
      print *, 'time_step_scheme must be eSSPIFRK3 or gear3'
      stop
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    save_restart_intvl = 0.d0
    restart_dir        = './restart/'
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=restart_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading restart_parameters failed"
    close(unit)

    if (save_restart_intvl == 0.d0) save_restart_intvl = dble(huge(1))

  end subroutine read_parameters

end module params_common

