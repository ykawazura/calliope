!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Parameter setting specific to RRMHD
!!          'params_common' is inherited
!-----------------------------------------------!
module params
  use p3dfft
  use params_common
  implicit none

  public  init_params
  public  nonlinear
  public  nupe , nupe_exp , nuz , nuz_exp
  public  etape, etape_exp, etaz, etaz_exp
  public  chipe, chipe_exp, chiz, chiz_exp
  public  shear, q
  private read_parameters

  logical  :: nonlinear
  real(r8) :: nupe , nuz
  real(r8) :: etape, etaz
  real(r8) :: chipe, chiz
  integer  :: nupe_exp , nuz_exp
  integer  :: etape_exp, etaz_exp
  integer  :: chipe_exp, chiz_exp
  logical, parameter  :: shear = .false.
  real(r8) :: q = 0.d0

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialization of run parameters,
!!          followed by input file reading
!-----------------------------------------------!
  subroutine init_params
    implicit none

    call init_params_common
    call read_parameters(inputfile)

  end subroutine init_params


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Read inputfile for various parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use file, only: get_unused_unit
    implicit none
    
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /operation_parameters/ nonlinear
    namelist /physical_parameters/ nupe , nupe_exp , nuz , nuz_exp , &
                                   etape, etape_exp, etaz, etaz_exp, &
                                   chipe, chipe_exp, chiz, chiz_exp

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    nonlinear = .false.
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=operation_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading operation_parameters failed"
    close(unit)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    nupe  = 0.d0
    etape = 0.d0
    chipe = 0.d0
    nuz   = 0.d0
    etaz  = 0.d0
    chiz  = 0.d0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=physical_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading physical_parameters failed"
    close(unit)

  end subroutine read_parameters

end module params
