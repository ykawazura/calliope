!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../params_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    25 Feb 2021
!! @brief   Parameter setting specific to MHD_COMP_ISOTH
!!          'params_common' is inherited
!-----------------------------------------------!
module params
  use p3dfft
  use params_common
  implicit none

  public  init_params
  public  nonlinear
  public  nu   , eta  , lmd  
  public  nu_h , eta_h, lmd_h
  public  nu_h_exp, eta_h_exp, lmd_h_exp
  public  beta0, cs2va2
  public  rho_min
  public  shear, q
  public  formulation, visc_form
  public  filter_lnrho, filter_mom, filter_mag
  public  cfl_visc, mass_correction, s_abort, debug_nan_step
  public  C_shock, shock_nmax, shock_nsmooth, shock_dl, shock_selftest, C_eta
  public  svv_eps0, svv_mN, svv_above_kc, svv_form, cfl_svv, svv_pi_init, svv_pi_check
  public  write_intvl_1D, write_intvl_cons, budget
  private read_parameters

  logical  :: nonlinear
  real(r8) :: nu  , eta  , lmd  
  real(r8) :: nu_h, eta_h, lmd_h
  integer  :: nu_h_exp, eta_h_exp, lmd_h_exp
  real(r8) :: beta0, cs2va2
  real(r8) :: rho_min
  logical  :: shear
  real(r8) :: q

  ! numerics (ln(rho) formulation)
  character(20) :: formulation   ! 'C': state = (ln rho, m, B)
  character(20) :: visc_form     ! 'flux': q = rho zeta div(u) on the grid (only 'flux' is implemented)
  logical  :: filter_lnrho, filter_mom, filter_mag ! apply dealias filter to each state after each RK stage
  real(r8) :: cfl_visc           ! dt_visc = cfl_visc/(largest damping rate of zeta and the SVV)
  logical  :: mass_correction    ! uniform correction ln rho <- ln rho + ln(M0/M) after each step
  real(r8) :: s_abort            ! stop if |ln rho| exceeds this value
  integer  :: debug_nan_step     ! inject NaN at this step to test the abort check (< 0: off)

  ! shock viscosity: zeta = C_shock*Delta^2*S[max(0, -div u)]
  real(r8) :: C_shock
  integer  :: shock_nmax, shock_nsmooth
  real(r8) :: shock_dl           ! Delta; <= 0 means min of the grid spacing in the active directions
  logical  :: shock_selftest     ! write zeta of an analytic D (shock_selftest.dat) for the decomposition test
  ! shock resistivity: eta_shock = C_eta*Delta^2*S[max(0, -div u)] (the detector of zeta), dB/dt = -curl_N(eta_shock J)
  real(r8) :: C_eta              ! <= 0: off

  ! spectral vanishing viscosity (Tadmor 1989) with the kernel Q(|k_d|) and eps_d = svv_eps0*dl_d*max(|u_d| + cf_d)
  ! (see shock_viscosity); B: rate lambda(k) = sum_d eps_d k_d^2 Q(|k_d|)
  real(r8) :: svv_eps0           ! <= 0: off
  real(r8) :: svv_mN             ! Q = 0 for |k_d| <= svv_mN*k_c,d (k_c,d = 2/3 of the Nyquist wavenumber)
  logical  :: svv_above_kc       ! Q = 1 (T) or 0 (F) for |k_d| >= k_c,d (only ln rho has such modes, when filter_lnrho = F)
  character(20) :: svv_form      ! rho and m: 'entropy' (flux form in the entropy variables; dE/dt <= 0 for any state)
                                 !            or 'linear' (-lambda(k) applied to m and, mass-conserving, to rho)
  real(r8) :: cfl_svv            ! 'entropy': dt <= cfl_svv/lambda_svv, lambda_svv by the power method (<= 2.51 of eSSPIFRK3)
  integer  :: svv_pi_init        ! 'entropy': power iterations at the first step (one per step afterwards)
  integer  :: svv_pi_check       ! 'entropy': every svv_pi_check steps the power method is run until converged and
                                 !            the run stops if dt*lambda_svv > 2.51 (<= 0: off)

  ! diagnostics specific to this model
  real(r8) :: write_intvl_1D     ! interval of line output (out1d/); <= 0: off
  real(r8) :: write_intvl_cons   ! interval of conservation.dat; 0: every step, < 0: off
  logical  :: budget             ! output energy/mass budget at every RK stage (budget.dat)

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
    namelist /physical_parameters/ nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, lmd, lmd_h, lmd_h_exp, &
                                   beta0, rho_min, shear, q
    namelist /numerics_parameters/ formulation, visc_form, filter_lnrho, filter_mom, filter_mag, &
                                   cfl_visc, mass_correction, s_abort, debug_nan_step
    namelist /shock_parameters/ C_shock, shock_nmax, shock_nsmooth, shock_dl, shock_selftest, C_eta
    namelist /svv_parameters/ svv_eps0, svv_mN, svv_above_kc, svv_form, cfl_svv, svv_pi_init, svv_pi_check
    namelist /isoth_diagnostics/ write_intvl_1D, write_intvl_cons, budget

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
    nu        = 0.d0
    nu_h      = 0.d0
    eta       = 0.d0
    eta_h     = 0.d0
    lmd       = 0.d0
    lmd_h     = 0.d0
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

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    formulation     = 'C'
    visc_form       = 'flux'
    filter_lnrho    = .false. ! the projection of ln rho changes sum exp(s) (mass)
    filter_mom      = .true.
    filter_mag      = .true.
    cfl_visc        = 1.d0
    mass_correction = .false.
    s_abort         = 1.d2
    debug_nan_step  = -1
    C_shock         = 2.d0
    shock_nmax      = 1
    shock_nsmooth   = 1
    shock_dl        = 0.d0
    shock_selftest  = .false.
    C_eta           = 0.d0
    svv_eps0        = 0.5d0
    svv_mN          = 0.5d0
    svv_above_kc    = .true.
    svv_form        = 'entropy'
    cfl_svv         = 2.d0
    svv_pi_init     = 50
    svv_pi_check    = 100
    write_intvl_1D  = 0.d0
    write_intvl_cons = 0.d0
    budget          = .false.
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')
    read(unit,nml=numerics_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading numerics_parameters failed (defaults are used)"
    rewind(unit)
    read(unit,nml=shock_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading shock_parameters failed (defaults are used)"
    rewind(unit)
    read(unit,nml=svv_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading svv_parameters failed (defaults are used)"
    rewind(unit)
    read(unit,nml=isoth_diagnostics,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading isoth_diagnostics failed (defaults are used)"
    close(unit)

  end subroutine read_parameters

end module params
