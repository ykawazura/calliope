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
  public  dealias
  public  nwrite, nwrite_fld_section, nwrite_fld_3D, nwrite_kpar, nwrite_SF2, SF2_nsample
  public  nsave_restart
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

  character(len=100) :: dealias

  integer            :: nwrite, nwrite_fld_section, nwrite_fld_3D, nwrite_kpar, nwrite_SF2, SF2_nsample
  integer            :: nsave_restart
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

    namelist /dealias_parameters/ dealias
    namelist /diagnostics_parameters/ nwrite, nwrite_fld_section, nwrite_fld_3D, nwrite_kpar, &
                                      nwrite_SF2, SF2_nsample, nsave_restart, restart_dir, output_modes

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
    nwrite             = 0
    nwrite_fld_section = 0
    nwrite_fld_3D      = 0
    nwrite_kpar        = 0
    nwrite_SF2         = 0
    SF2_nsample        = 100000
    nsave_restart      = 0
    restart_dir = './restart/'
    output_modes(:)    = 0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=diagnostics_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading diagnostics_parameters failed"
    close(unit)

    if (nwrite             == 0) nwrite             = huge(1)
    if (nwrite_fld_section == 0) nwrite_fld_section = huge(1)
    if (nwrite_fld_3D      == 0) nwrite_fld_3D      = huge(1)
    if (nwrite_kpar        == 0) nwrite_kpar        = huge(1)
    if (nwrite_SF2         == 0) nwrite_SF2         = huge(1)
    if (nsave_restart      == 0) nsave_restart      = huge(1)

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

