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
  public  dealias
  public  write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, write_intvl_SF2, SF2_nsample
  public  save_restart_intvl
  public  restart_dir
  public  series_output, n_series_modes, series_modes
  private read_parameters

  !> square root of -1.
  complex(r8), parameter :: zi = ( 0.d0 , 1.d0 )

  !> pi to quad precision
  real(r8), parameter :: pi = &
       3.14159265358979323846264338327950288419716939938d0

  character(len=100) :: runname = 'calliope'
  character(len=100) :: inputfile

  character(len=100) :: time_step_scheme
  character(len=100) :: dealias

  real(r8)           :: write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, write_intvl_SF2, SF2_nsample
  real(r8)           :: save_restart_intvl
  character(len=100) :: restart_dir
  integer            :: n_series_modes
  integer, allocatable :: series_modes(:, :)
  logical            :: series_output = .false.

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
    integer  :: i, ierr

    namelist /scheme_parameters/ time_step_scheme
    namelist /dealias_parameters/ dealias
    namelist /diagnostics_parameters/ write_intvl, write_intvl_2D, write_intvl_3D, write_intvl_kpar, &
                                      write_intvl_SF2, SF2_nsample, save_restart_intvl, restart_dir, output_modes

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
    dealias = '2/3'
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=dealias_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading dealias_parameters failed"
    close(unit)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    write_intvl        = 0.d0
    write_intvl_2D     = 0.d0
    write_intvl_3D     = 0.d0
    write_intvl_kpar   = 0.d0
    write_intvl_SF2    = 0.d0
    SF2_nsample        = 100000
    save_restart_intvl = 0.d0
    restart_dir        = './restart/'
    output_modes(:)    = 0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=diagnostics_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading diagnostics_parameters failed"
    close(unit)

    if (write_intvl        == 0.d0) write_intvl        = dble(huge(1))
    if (write_intvl_2D     == 0.d0) write_intvl_2D     = dble(huge(1))
    if (write_intvl_3D     == 0.d0) write_intvl_3D     = dble(huge(1))
    if (write_intvl_kpar   == 0.d0) write_intvl_kpar   = dble(huge(1))
    if (write_intvl_SF2    == 0.d0) write_intvl_SF2    = dble(huge(1))
    if (save_restart_intvl == 0.d0) save_restart_intvl = dble(huge(1))

    ! Time series output of modes
    if(count(output_modes /= 0) > 0) then
      if(mod(count(output_modes /= 0), 3) /= 0) then
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        print *, '!              Error!              !'
        print *, '!  output modes must be multiple   !'
        print *, '!  of 3; set of (ikx, iky, ikz)    !'
        print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
        stop
      endif

      n_series_modes = count(output_modes /= 0) / 3
      allocate(series_modes(n_series_modes, 3))

      do i = 1, n_series_modes
        series_modes(i, 1) = output_modes((i-1)*3 + 1)
        series_modes(i, 2) = output_modes((i-1)*3 + 2)
        series_modes(i, 3) = output_modes((i-1)*3 + 3)
      enddo

      series_output = .true.
    endif

  end subroutine read_parameters

end module params_common

