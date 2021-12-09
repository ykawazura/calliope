!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    25 Feb 2021
!! @brief   Parameter setting specific to MHD_INCOMP
!!          'params_common' is inherited
!-----------------------------------------------!
module params
  use p3dfft
  use params_common
  implicit none

  public  init_params
  public  nonlinear
  public  nu , eta, lmd
  public  nu_exp, eta_exp, lmd_exp
  public  beta0, cs2va2
  public  rho_min
  public  shear, q
  private read_parameters

  logical  :: nonlinear
  real(r8) :: nu, eta, lmd
  integer  :: nu_exp, eta_exp, lmd_exp
  real(r8) :: beta0, cs2va2
  real(r8) :: rho_min
  logical  :: shear
  real(r8) :: q

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
    namelist /physical_parameters/ nu, nu_exp, eta, eta_exp, lmd, lmd_exp, beta0, rho_min, shear, q

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
    nu      = 0.d0
    eta     = 0.d0
    lmd     = 0.d0
    beta0   = 1.d0
    rho_min = 1.d-2
    shear   = .false.
    q       = 0.d0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=physical_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading physical_parameters failed"
    close(unit)

    cs2va2 = beta0/2.d0

  end subroutine read_parameters

end module params
