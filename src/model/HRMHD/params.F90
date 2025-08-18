!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   Parameter setting specific to HRMHD
!!          'params_common' is inherited
!-----------------------------------------------!
module params
  use p3dfft
  use params_common
  implicit none

  public  init_params
  public  nonlinear
  public  rho, sgm
  public  nupe_x , nupe_x_exp , nupe_z , nupe_z_exp
  public  nupa_x , nupa_x_exp , nupa_z , nupa_z_exp
  public  etape_x, etape_x_exp, etape_z, etape_z_exp
  public  etapa_x, etapa_x_exp, etapa_z, etapa_z_exp
  public  shear, q
  private read_parameters

  logical  :: nonlinear
  real(r8) :: rho, sgm
  real(r8) :: nupe_x , nupe_z
  real(r8) :: etape_x, etape_z
  real(r8) :: nupa_x , nupa_z
  real(r8) :: etapa_x, etapa_z
  integer  :: nupe_x_exp , nupe_z_exp
  integer  :: etape_x_exp, etape_z_exp
  integer  :: nupa_x_exp , nupa_z_exp
  integer  :: etapa_x_exp, etapa_z_exp
  logical, parameter  :: shear = .false.
  real(r8), parameter :: q = 0.d0

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
    namelist /physical_parameters/ nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                                   nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                                   etape_x, etape_x_exp, etape_z, etape_z_exp, &
                                   etapa_x, etapa_x_exp, etapa_z, etapa_z_exp, rho, sgm

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
    nupe_x  = 0.d0
    nupe_z  = 0.d0
    etape_x = 0.d0
    etape_z = 0.d0
    nupa_x  = 0.d0
    nupa_z  = 0.d0
    etapa_x = 0.d0
    etapa_z = 0.d0
    rho     = 0.d0
    sgm     = 1.d0
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=physical_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading physical_parameters failed"
    close(unit)

  end subroutine read_parameters

end module params
