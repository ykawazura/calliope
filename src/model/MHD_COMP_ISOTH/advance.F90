!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Time stepping for MHD_COMP_ISOTH
!!
!!          State (Fourier): s = ln(rho), m = rho*u, B
!!            ds/dt = -(div_N m)/rho
!!            dm/dt = -div_N [m u + (cs^2 rho + B^2/2 - q) I - B B]
!!            dB/dt = -curl_N (B x u)
!!          q = rho*zeta*D, D = div_N u (see shock_viscosity).
!!          div_N = F^-1 i k_d M F with k_d = 0 at Nyquist and M the
!!          filter of the differentiated field.
!!
!!          rho = exp(s) is positive by construction; no floor is
!!          applied. A non-finite value or |s| > s_abort stops the run.
!!          Only eSSPIFRK3 is implemented. Its effective weights are
!!          b = (1/4, 3/16, 9/16) at c = (0, 2/3, 2/3), which are used
!!          for the energy budget.
!-----------------------------------------------!
module advance
  use p3dfft
  use fields, only: nfields
  use fields, only: ilnrho
  use fields, only: imx, imy, imz
  use fields, only: ibx, iby, ibz
  implicit none

  public solve

  logical :: is_allocated = .false.

  integer :: counter = 0
  integer :: istep_adv = 0 ! number of steps completed in this run
  complex(r8), allocatable, dimension(:,:,:)   :: lnrho_new
  complex(r8), allocatable, dimension(:,:,:)   ::    mx_new
  complex(r8), allocatable, dimension(:,:,:)   ::    my_new
  complex(r8), allocatable, dimension(:,:,:)   ::    mz_new
  complex(r8), allocatable, dimension(:,:,:)   ::    bx_new
  complex(r8), allocatable, dimension(:,:,:)   ::    by_new
  complex(r8), allocatable, dimension(:,:,:)   ::    bz_new
  complex(r8), allocatable, dimension(:,:,:,:) :: w
  real   (r8), allocatable, dimension(:,:,:,:) :: w_r, flx_r
  real   (r8), allocatable, dimension(:,:,:)   :: rho_r, ux_r, uy_r, uz_r
  real   (r8), allocatable, dimension(:,:,:)   :: div_r, zeta_r, q_r ! D, zeta, q = rho*zeta*D
  real   (r8), allocatable, dimension(:,:,:)   :: svv_r ! SVV term of the s equation, (d rho/dt)_svv/rho
  complex(r8), allocatable, dimension(:,:,:,:) :: svv_mh ! SVV term of m ('entropy'), filtered and normalized as the state
  complex(r8), allocatable, dimension(:,:,:,:) :: svv_uh ! F u/(nlx nly nlz) of the current stage ('entropy')
  complex(r8), allocatable, dimension(:,:,:,:) :: flx, exp_terms
  real   (r8), allocatable, dimension(:,:)     :: kxt ! sheared kx
  real   (r8), allocatable, dimension(:,:)     :: kxd ! kx for derivatives: kx_d (shear = F) or kxt (shear = T)

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: lnrho_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::    mx_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::    my_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::    mz_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::    bx_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::    by_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::    bz_tmp
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms0
  real(r8), parameter :: b_rk(3) = (/0.25d0, 0.1875d0, 0.5625d0/) ! effective weights of eSSPIFRK3

  ! per-field settings: filter after each stage, linear dissipation (lmd for s, nu for m, eta for B)
  logical  :: filt_on  (nfields)
  real(r8) :: dcoef    (nfields), dcoef_h(nfields)
  integer  :: dcoef_exp(nfields)

  integer :: cfl_unit, cons_unit, budget_unit
  integer :: nline1d = 0

  ! Forward FFT variables
  integer, parameter :: nftran = 10
  integer, parameter :: iflx_s   = 1                               ! -(div_N m)/rho
  integer, parameter :: iflx_mxx = 2, iflx_mxy = 3, iflx_mxz = 4   ! m u - B B + (cs^2 rho + B^2/2 - q) I (tensor)
  integer, parameter ::               iflx_myy = 5, iflx_myz = 6   ! m u - B B + (cs^2 rho + B^2/2 - q) I (tensor)
  integer, parameter ::                             iflx_mzz = 7   ! m u - B B + (cs^2 rho + B^2/2 - q) I (tensor)
  integer, parameter :: iflx_bx  = 8, iflx_by  = 9, iflx_bz  = 10  ! B x u

  ! Volume integrals of the state (conservation.dat)
  integer, parameter :: nmon = 7
  integer, parameter :: imon_M  = 1                         ! mass
  integer, parameter :: imon_px = 2, imon_py = 3, imon_pz = 4 ! momentum
  integer, parameter :: imon_ek = 5, imon_em = 6, imon_ei = 7 ! kinetic, magnetic, internal energy
  real(r8) :: mon_n(nmon)          ! at the beginning of the current step (after the mass correction)
  real(r8) :: mass0 = 0.d0         ! reference mass = mass at the first step of this run
  real(r8) :: dmass_corr     = 0.d0 ! ln(M0/M) applied at the current step
  real(r8) :: dmass_corr_tot = 0.d0

  ! Of the latest RHS evaluation (volume integrals)
  real(r8) :: dMdt_raw  ! -sum div_N m (zero in exact arithmetic)
  real(r8) :: dMdt_abs  !  sum |div_N m| (scale of dMdt_raw)
  real(r8) :: visc_diss !  sum rho zeta D^2 (>= 0)

  ! SVV coefficients eps_d of the current step (set at stage 1 and fixed over the step)
  real(r8) :: eps_svv(3) = 0.d0
  ! 'entropy': D_svv = sum [a t^2 + b |g_m|^2 + kap rho g_h^2] dV (>= 0) of the latest RHS evaluation
  real(r8) :: svv_diss = 0.d0
  ! dt before reset_dt of the current step (for the check of dt*lambda_svv)
  real(r8) :: dt_prev = 0.d0

  ! Energy and mass budget accumulated over the RK stages with weights b_rk (budget.dat)
  integer, parameter :: nbgt = 16
  integer, parameter :: ib_dMdt_raw = 1, ib_dMdt_proj = 2, ib_dMdt_abs = 3
  integer, parameter :: ib_eps_N    = 4, ib_eps_shock = 5, ib_visc     = 6, ib_R_zeta = 7
  integer, parameter :: ib_eps_force = 8, ib_eps_res  = 9, ib_eps_nu   = 10, ib_eps_lmd = 11
  integer, parameter :: ib_eps_svv  = 12, ib_eps_svv_b = 13 ! SVV of (rho, m) and of B
  integer, parameter :: ib_D_svv    = 14, ib_R_svv     = 15 ! D_svv and R_svv = eps_svv + D_svv ('entropy')
  integer, parameter :: ib_A_svv    = 16 ! sum |dE/dU.U_svv| dV over (rho, m): rounding scale of eps_svv
  real(r8) :: bgt(nbgt) = 0.d0
  real(r8) :: R_zeta_ratio_max = 0.d0
  real(r8) :: R_svv_ratio_max  = 0.d0
  real(r8) :: R_svv_scale_max  = 0.d0
  real(r8) :: t_bgt, dt_bgt
  logical  :: bgt_pending = .false.

contains


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Solve the time evolution (eSSPIFRK3)
!-----------------------------------------------!
  subroutine solve
    use fields, only: lnrho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: lnrho_old
    use fields, only: mx_old, my_old, mz_old
    use fields, only: bx_old, by_old, bz_old
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, tt, nstep
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: shear, budget
    use shearing_box, only: k2t, k2t_inv, tsc, tremap
    use force, only: fmx, fmz, fmy, fmx_old, fmy_old, fmz_old, driven, update_force, get_force
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use diagnostics_common, only: series_output
    implicit none
    real(r8) :: imp_terms_tintg0(nfields), imp_terms_tintg1(nfields)
    real(r8) :: imp_terms_tintg2(nfields), imp_terms_tintg3(nfields)
    integer :: i, j, k, n

    if (proc0) call put_time_stamp(timer_advance)

    ! initialize tmp fields
    if(counter == 0) then
      call init_work_fields
      counter = 1
    endif

    ! Calcualte force terms
    if (driven) then
      ! at n
      call get_force('mx', fmx)
      call get_force('my', fmy)
      call get_force('mz', fmz)
    endif

    !---------------  RK 1st step  ---------------
    ! Calcualte nonlinear terms (also CFL, monitor, and mass correction of the real-space s)
    call get_nonlinear_terms(lnrho, mx, my, mz, bx, by, bz, 1)
    if(dmass_corr /= 0.d0) call shift_mean(lnrho, dmass_corr)

    call get_explicit_terms(mx, my, mz, bx, by, bz)
    if(budget) call stage_budget(1, lnrho, mx, my, mz, bx, by, bz)

    !$omp parallel do private(j, k, i, n, imp_terms_tintg0, imp_terms_tintg1)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          ! Calculate time integral of implicit terms
          do n = 1, nfields
            call get_imp_terms_tintg(imp_terms_tintg0(n), tsc               , kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg1(n), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
          enddo

          ! update s
          call eSSPIFRK1(lnrho_tmp(i,k,j), lnrho(i,k,j), &
             exp_terms(i,k,j,ilnrho), &
             imp_terms_tintg1(ilnrho), imp_terms_tintg0(ilnrho) &
          )

          ! update m
          call eSSPIFRK1(mx_tmp(i,k,j), mx(i,k,j), &
             exp_terms(i,k,j,imx), &
             imp_terms_tintg1(imx), imp_terms_tintg0(imx) &
          )
          call eSSPIFRK1(my_tmp(i,k,j), my(i,k,j), &
             exp_terms(i,k,j,imy), &
             imp_terms_tintg1(imy), imp_terms_tintg0(imy) &
          )
          call eSSPIFRK1(mz_tmp(i,k,j), mz(i,k,j), &
             exp_terms(i,k,j,imz), &
             imp_terms_tintg1(imz), imp_terms_tintg0(imz) &
          )

          ! update b
          call eSSPIFRK1(bx_tmp(i,k,j), bx(i,k,j), &
             exp_terms(i,k,j,ibx), &
             imp_terms_tintg1(ibx), imp_terms_tintg0(ibx) &
          )
          call eSSPIFRK1(by_tmp(i,k,j), by(i,k,j), &
             exp_terms(i,k,j,iby), &
             imp_terms_tintg1(iby), imp_terms_tintg0(iby) &
          )
          call eSSPIFRK1(bz_tmp(i,k,j), bz(i,k,j), &
             exp_terms(i,k,j,ibz), &
             imp_terms_tintg1(ibz), imp_terms_tintg0(ibz) &
          )
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! save explicit terms at the previous step
    !$omp workshare
    exp_terms0 = exp_terms
    !$omp end workshare

    ! Dealiasing
    call apply_filter(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)

    !---------------  RK 2nd step  ---------------
    ! Calcualte force terms
    if (driven) then
      ! at n + 2/3
      call update_force(2.d0/3.d0*dt)
      call get_force('mx', fmx)
      call get_force('my', fmy)
      call get_force('mz', fmz)

      ! go to n + 1
      call update_force(1.d0/3.d0*dt)
    endif

    ! Calcualte kxt at n + 2/3
    if(shear) call set_kx(tsc + 2.d0/3.d0*dt)

    ! Calcualte nonlinear terms
    call get_nonlinear_terms(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, 2)

    call get_explicit_terms(mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
    if(budget) call stage_budget(2, lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)

    !$omp parallel do private(j, k, i, n, imp_terms_tintg0, imp_terms_tintg2)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          ! Calculate time integral of implicit terms
          do n = 1, nfields
            call get_imp_terms_tintg(imp_terms_tintg0(n), tsc               , kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg2(n), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
          enddo

          ! update s
          call eSSPIFRK2(lnrho_tmp(i,k,j), lnrho_tmp(i,k,j), lnrho(i,k,j), &
             exp_terms(i,k,j,ilnrho), &
             imp_terms_tintg2(ilnrho), imp_terms_tintg0(ilnrho) &
          )

          ! update m
          call eSSPIFRK2(mx_tmp(i,k,j), mx_tmp(i,k,j), mx(i,k,j), &
             exp_terms(i,k,j,imx), &
             imp_terms_tintg2(imx), imp_terms_tintg0(imx) &
          )
          call eSSPIFRK2(my_tmp(i,k,j), my_tmp(i,k,j), my(i,k,j), &
             exp_terms(i,k,j,imy), &
             imp_terms_tintg2(imy), imp_terms_tintg0(imy) &
          )
          call eSSPIFRK2(mz_tmp(i,k,j), mz_tmp(i,k,j), mz(i,k,j), &
             exp_terms(i,k,j,imz), &
             imp_terms_tintg2(imz), imp_terms_tintg0(imz) &
          )

          ! update b
          call eSSPIFRK2(bx_tmp(i,k,j), bx_tmp(i,k,j), bx(i,k,j), &
             exp_terms(i,k,j,ibx), &
             imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
          )
          call eSSPIFRK2(by_tmp(i,k,j), by_tmp(i,k,j), by(i,k,j), &
             exp_terms(i,k,j,iby), &
             imp_terms_tintg2(iby), imp_terms_tintg0(iby) &
          )
          call eSSPIFRK2(bz_tmp(i,k,j), bz_tmp(i,k,j), bz(i,k,j), &
             exp_terms(i,k,j,ibz), &
             imp_terms_tintg2(ibz), imp_terms_tintg0(ibz) &
          )
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! Dealiasing
    call apply_filter(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)

    !---------------  RK 3rd step  ---------------
    ! Calcualte nonlinear terms
    call get_nonlinear_terms(lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, 3)

    call get_explicit_terms(mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)
    if(budget) call stage_budget(3, lnrho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp)

    !$omp parallel do private(j, k, i, n, imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          ! Calculate time integral of implicit terms
          do n = 1, nfields
            call get_imp_terms_tintg(imp_terms_tintg0(n), tsc               , kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg2(n), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
            call get_imp_terms_tintg(imp_terms_tintg3(n), tsc +           dt, kx(i), ky(j), kz(k), &
                                     dcoef(n), dcoef_h(n), dcoef_exp(n))
          enddo

          ! update s
          call eSSPIFRK3(lnrho_new(i,k,j), lnrho_tmp(i,k,j), lnrho(i,k,j), &
             exp_terms(i,k,j,ilnrho), exp_terms0(i,k,j,ilnrho), &
             imp_terms_tintg3(ilnrho), imp_terms_tintg2(ilnrho), imp_terms_tintg0(ilnrho) &
          )

          ! update m
          call eSSPIFRK3(mx_new(i,k,j), mx_tmp(i,k,j), mx(i,k,j), &
             exp_terms(i,k,j,imx), exp_terms0(i,k,j,imx), &
             imp_terms_tintg3(imx), imp_terms_tintg2(imx), imp_terms_tintg0(imx) &
          )
          call eSSPIFRK3(my_new(i,k,j), my_tmp(i,k,j), my(i,k,j), &
             exp_terms(i,k,j,imy), exp_terms0(i,k,j,imy), &
             imp_terms_tintg3(imy), imp_terms_tintg2(imy), imp_terms_tintg0(imy) &
          )
          call eSSPIFRK3(mz_new(i,k,j), mz_tmp(i,k,j), mz(i,k,j), &
             exp_terms(i,k,j,imz), exp_terms0(i,k,j,imz), &
             imp_terms_tintg3(imz), imp_terms_tintg2(imz), imp_terms_tintg0(imz) &
          )

          ! update b
          call eSSPIFRK3(bx_new(i,k,j), bx_tmp(i,k,j), bx(i,k,j), &
             exp_terms(i,k,j,ibx), exp_terms0(i,k,j,ibx), &
             imp_terms_tintg3(ibx), imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
          )
          call eSSPIFRK3(by_new(i,k,j), by_tmp(i,k,j), by(i,k,j), &
             exp_terms(i,k,j,iby), exp_terms0(i,k,j,iby), &
             imp_terms_tintg3(iby), imp_terms_tintg2(iby), imp_terms_tintg0(iby) &
          )
          call eSSPIFRK3(bz_new(i,k,j), bz_tmp(i,k,j), bz(i,k,j), &
             exp_terms(i,k,j,ibz), exp_terms0(i,k,j,ibz), &
             imp_terms_tintg3(ibz), imp_terms_tintg2(ibz), imp_terms_tintg0(ibz) &
          )
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! save fields at the previous step
    !$omp workshare
    lnrho_old = lnrho
       mx_old = mx
       my_old = my
       mz_old = mz
       bx_old = bx
       by_old = by
       bz_old = bz
    !$omp end workshare

    if (driven) then
      !$omp workshare
      fmx_old = fmx
      fmy_old = fmy
      fmz_old = fmz
      !$omp end workshare
    endif

    ! Dealiasing
    call apply_filter(lnrho_new, mx_new, my_new, mz_new, bx_new, by_new, bz_new)

    !$omp workshare
    lnrho = lnrho_new
       mx =    mx_new
       my =    my_new
       mz =    mz_new
       bx =    bx_new
       by =    by_new
       bz =    bz_new
    !$omp end workshare

    tt  = tt  + dt
    tsc = tsc + dt

    if(shear .and. tsc > tremap) call remap

    if(shear) then
      call set_kx(tsc)
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do i = ikx_st, ikx_en
          do k = ikz_st, ikz_en
            k2t(i,k,j) = kxt(i,j)**2 + ky(j)**2 + kz(k)**2

            if(k2t(i,k,j) == 0.d0) then
              k2t_inv(i,k,j) = 0.d0
            else
              k2t_inv(i,k,j) = 1.d0/k2t(i,k,j)
            endif
          enddo
        enddo
      enddo
      !$omp end parallel do
    endif

    ! Div b cleaing
    call div_free(bx, by, bz)

    if(series_output) call output_series_modes

    istep_adv = istep_adv + 1

    ! Monitor and line output of the final state (the RHS is not evaluated)
    if(istep_adv == nstep) call get_nonlinear_terms(lnrho, mx, my, mz, bx, by, bz, 0)

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Check the run parameters, allocate
!!          fields used only here, and open the
!!          output files of this module
!-----------------------------------------------!
  subroutine init_work_fields
    use fields, only: lnrho, mx, my, mz, bx, by, bz
    use grid, only: kx, ky, kz, pruned
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: nonlinear, shear, time_step_scheme
    use params, only: nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, lmd, lmd_h, lmd_h_exp
    use params, only: formulation, visc_form, filter_lnrho, filter_mom, filter_mag
    use params, only: C_shock, s_abort, budget, write_intvl_cons, svv_eps0, svv_form, cfl_svv
    use shearing_box, only: tsc, k2t, k2t_inv
    use shock_viscosity, only: init_shock_viscosity, abort_run
    use file, only: open_output_file
    use mp, only: proc0
    implicit none
    integer :: i, j, k

    if(time_step_scheme /= 'eSSPIFRK3') &
      call abort_run('MHD_COMP_ISOTH: time_step_scheme = '//trim(time_step_scheme)//' is not implemented (use eSSPIFRK3)')
    if(trim(formulation) /= 'C') &
      call abort_run('MHD_COMP_ISOTH: formulation = '//trim(formulation)//' is not implemented (use C)')
    if(trim(visc_form) /= 'flux') &
      call abort_run('MHD_COMP_ISOTH: visc_form = '//trim(visc_form)//' is not implemented yet (use flux)')
    if(.not. nonlinear) &
      call abort_run('MHD_COMP_ISOTH: nonlinear = F is not supported (mass flux and pressure are nonlinear terms)')
    if(shear .and. C_shock > 0.d0) &
      call abort_run('MHD_COMP_ISOTH: shear = T with C_shock > 0 has not been verified')
    if(shear .and. svv_eps0 > 0.d0) &
      call abort_run('MHD_COMP_ISOTH: shear = T with svv_eps0 > 0 is not implemented (the SVV kernel uses the unsheared kx)')
    if(trim(svv_form) /= 'entropy' .and. trim(svv_form) /= 'linear') &
      call abort_run('MHD_COMP_ISOTH: svv_form = '//trim(svv_form)//' is not implemented (use entropy or linear)')
    if(trim(svv_form) == 'entropy' .and. cfl_svv <= 0.d0) &
      call abort_run('MHD_COMP_ISOTH: cfl_svv must be > 0')
    if(pruned .and. .not. (filter_lnrho .and. filter_mom .and. filter_mag)) &
      call abort_run('MHD_COMP_ISOTH: pruned = T truncates every field; set filter_lnrho = filter_mom = filter_mag = T')
    if(s_abort > 7.d2) &
      call abort_run('MHD_COMP_ISOTH: s_abort must be <= 700 so that exp(s) does not overflow')
    if(proc0 .and. (lmd /= 0.d0 .or. lmd_h /= 0.d0)) then
      print *
      print *, 'WARNING: lmd /= 0 diffuses ln(rho); this term does not conserve mass'
      print *
    endif
    if(proc0 .and. svv_eps0 > 0.d0 .and. trim(svv_form) == 'entropy' .and. filter_lnrho) then
      print *
      print *, 'WARNING: filter_lnrho = T filters the SVV term of ln(rho); dE/dt <= 0 of the entropy-stable SVV'
      print *, '         is then not guaranteed (see R_svv in budget.dat)'
      print *
    endif

    filt_on(ilnrho) = filter_lnrho
    filt_on(imx:imz) = filter_mom
    filt_on(ibx:ibz) = filter_mag
    dcoef(ilnrho) = lmd; dcoef_h(ilnrho) = lmd_h; dcoef_exp(ilnrho) = lmd_h_exp
    dcoef(imx:imz) = nu ; dcoef_h(imx:imz) = nu_h ; dcoef_exp(imx:imz) = nu_h_exp
    dcoef(ibx:ibz) = eta; dcoef_h(ibx:ibz) = eta_h; dcoef_exp(ibx:ibz) = eta_h_exp

    call init_shock_viscosity

    allocate(kxt(ikx_st:ikx_en, iky_st:iky_en))
    allocate(kxd(ikx_st:ikx_en, iky_st:iky_en))
    call set_kx(tsc)

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do i = ikx_st, ikx_en
        do k = ikz_st, ikz_en
          k2t(i, k, j) = kxt(i, j)**2 + ky(j)**2 + kz(k)**2
          if(k2t(i, k, j) == 0.d0) then
            k2t_inv(i, k, j) = 0.d0
          else
            k2t_inv(i, k, j) = 1.0d0/k2t(i, k, j)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

    call allocate_advance

    ! The RK stages keep the filtered fields; start from the filtered initial condition
    if(proc0) then
      print *
      print '(" MHD_COMP_ISOTH: the initial condition is filtered (filter_lnrho, filter_mom, filter_mag = ", 3L2, ")")', &
            filter_lnrho, filter_mom, filter_mag
      print *
    endif
    call apply_filter(lnrho, mx, my, mz, bx, by, bz)

    if(proc0) then
      call open_output_file (cfl_unit, 'cfl.dat')
      write(cfl_unit, '(a)') '# 1:t 2:dt 3:dt_adv 4:dt_visc 5:max(|ux|+cf_x) 6:max(|uy|+cf_y) 7:max(|uz|+cf_z) ' &
                           //'8:zeta_max 9:max(rho*zeta)/(rho_min*zeta_max) 10:eps_svv_x 11:eps_svv_y 12:eps_svv_z ' &
                           //'13:lam_svv 14:F_svv(=13/sum_d eps_d max k_d^2 Q) 15:dt*lam_svv 16:n_iter'
      write(cfl_unit, '(a)') '# 13-16: power method of the entropy-stable SVV (0 for svv_form = linear)'

      if(write_intvl_cons >= 0.d0) then
        call open_output_file (cons_unit, 'conservation.dat')
        write(cons_unit, '(a)') '# 1:t 2:dt 3:M 4:(M-M0)/M0 5:Px 6:Py 7:Pz 8:E_kin 9:E_mag 10:E_int 11:E ' &
                              //'12:rho_min 13:rho_max 14:zeta_max 15:D_min 16:sum(rho*zeta*D^2)dV ' &
                              //'17:dmass_corr 18:dmass_corr_tot'
        write(cons_unit, '(a)') '# M and E are before the mass correction of this step; M0 = M at the first step of this run'
      endif

      if(budget) then
        call open_output_file (budget_unit, 'budget.dat')
        write(budget_unit, '(a)') '# one line per step [t, t+dt]; rates are sum_i b_i X_i over the RK stages, b = (1/4,3/16,9/16)'
        write(budget_unit, '(a)') '# 1:t 2:dt 3:dMdt_actual 4:dMdt_raw 5:dMdt_proj 6:dMdt_abs 7:dEdt_actual 8:dEdt_sum ' &
                                //'9:err_t(=7-8) 10:eps_ideal(=eps_N-11-15-21) 11:eps_shock 12:-visc_diss 13:R_zeta ' &
                                //'14:max|R_zeta|/visc_diss 15:eps_force 16:eps_res 17:eps_nu 18:eps_lmd 19:E_n 20:E_n+1 ' &
                                //'21:eps_svv(all) 22:eps_svv_B 23:-D_svv 24:R_svv(=21-22-23) 25:max|R_svv|/D_svv ' &
                                //'26:max|R_svv|/sum|dE/dU.U_svv|dV'
        write(budget_unit, '(a)') '# 23-26: entropy-stable SVV of (rho, m) (0 for svv_form = linear); 26 is relative to the ' &
                                //'rounding scale of 21-22 (25 is not when D_svv is itself at the rounding level)'
      endif
    endif

  end subroutine init_work_fields


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   kxt = kx + q*shear_flg*t*ky and the
!!          x-wavenumber for derivatives
!-----------------------------------------------!
  subroutine set_kx(t)
    use grid, only: kx, ky
    use grid, only: ikx_st, iky_st, ikx_en, iky_en
    use params, only: q, shear
    use shearing_box, only: shear_flg
    use shock_viscosity, only: kx_d
    implicit none
    real(r8), intent(in) :: t
    integer :: i, j

    do j = iky_st, iky_en
      do i = ikx_st, ikx_en
        kxt(i,j) = kx(i) + q*shear_flg*t*ky(j)
        if(shear) then
          kxd(i,j) = kxt(i,j)
        else
          kxd(i,j) = kx_d(i)
        endif
      enddo
    enddo

  end subroutine set_kx


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Calculate nonlinear terms via
!!          1. Inverse FFT of (s, m, B)
!!          2. Validity check, rho = exp(s),
!!             monitor and mass correction
!!             (istage <= 1), u = m/rho
!!          3. D, zeta, q; SVV coefficients and
!!             CFL (istage = 1)
!!          4. Fluxes in real space, SVV of ln rho
!!          5. Forward FFT
!!          istage = 1, 2, 3: RK stage
!!          istage = 0      : steps 1-3 only, for
!!                            the final state
!-----------------------------------------------!
  subroutine get_nonlinear_terms(lnrho, mx, my, mz, bx, by, bz, istage)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: zi, cs2va2, C_shock, cfl_visc, filter_mom, mass_correction, debug_nan_step
    use params, only: write_intvl_cons, write_intvl_1D, cfl_svv, svv_pi_check
    use dealias, only: filter
    use mp, only: proc0, max_allreduce, sum_allreduce
    use time, only: cfl, dt, tt, reset_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    use advance_common, only: dt_adjust_while_running
    use shock_viscosity, only: get_div_u, get_zeta, ky_d, kz_d, kd2_max
    use shock_viscosity, only: active_x, active_y, active_z, is_nonfinite, abort_run
    use shock_viscosity, only: svv_on, kq2x, kq2y, kq2z, svv_lam_m, svv_lam_s, svv_lam_b, get_svv_eps
    use shock_viscosity, only: svv_es, get_svv_es, svv_power, svv_lam
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: lnrho, mx, my, mz, bx, by, bz
    integer, intent(in) :: istage

    real(r8) :: mon_pre(nmon), sums(3), maxs(8)
    real(r8) :: vad2_r(3), va2_r, cf2_r, b2h_r, ptot_r
    real(r8) :: dV, mfilt, dt_adv, dt_visc, dt_cfl, rate_visc, F_svv
    integer  :: i, j, k, n, d, pi_nit
    logical  :: first, pi_check, pi_conv
    real(r8), parameter :: rk3_lim = 2.51d0 ! real-axis stability limit of SSPRK3 (2.5127)
    character(len=256) :: msg

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    dV    = dlx*dly*dlz
    first = istage == 1 .and. istep_adv == 0

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          w(i,k,j,ilnrho) = lnrho(i,k,j)
          w(i,k,j,imx   ) = mx   (i,k,j)
          w(i,k,j,imy   ) = my   (i,k,j)
          w(i,k,j,imz   ) = mz   (i,k,j)
          w(i,k,j,ibx   ) = bx   (i,k,j)
          w(i,k,j,iby   ) = by   (i,k,j)
          w(i,k,j,ibz   ) = bz   (i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 1. Inverse FFT
    if (proc0) call put_time_stamp(timer_fft)
    ! for some reason c2r_many doesn't work at Fugaku
    do n = 1, nfields
      call p3dfft_btran_c2r(w(:,:,:,n), w_r(:,:,:,n), 'tff')
    enddo
    if (proc0) call put_time_stamp(timer_fft)

    ! 2. Validity check (no floor), rho, monitor, mass correction
    if(istage == 1 .and. istep_adv == debug_nan_step .and. proc0) then
      print *, 'debug_nan_step: NaN is injected into mx at the first grid point of proc 0'
      w_r(ily_st, ilz_st, ilx_st, imx) = ieee_value(1.d0, ieee_quiet_nan)
    endif
    call check_state(istage)

    rho_r = exp(w_r(:,:,:,ilnrho))

    if(istage <= 1) then
      call get_monitor_sums(mon_pre)
      if(first) mass0 = mon_pre(imon_M)

      ! close the budget of the previous step with the state before the mass correction
      if(bgt_pending) call write_budget(mon_pre)

      dmass_corr = 0.d0
      if(istage == 1 .and. mass_correction .and. .not. first) then
        ! uniform shift s <- s + ln(M0/M); this does not restore the local conservation
        dmass_corr     = log(mass0/mon_pre(imon_M))
        dmass_corr_tot = dmass_corr_tot + dmass_corr
        w_r(:,:,:,ilnrho) = w_r(:,:,:,ilnrho) + dmass_corr
        rho_r = exp(w_r(:,:,:,ilnrho))
        call get_monitor_sums(mon_n)
      else
        mon_n = mon_pre
      endif
    endif

    !$omp workshare
    ux_r = w_r(:,:,:,imx)/rho_r
    uy_r = w_r(:,:,:,imy)/rho_r
    uz_r = w_r(:,:,:,imz)/rho_r
    !$omp end workshare

    ! 3. Shock viscosity: D = div_N u, zeta, q = rho*zeta*D
    if(svv_es) then
      call get_div_u(ux_r, uy_r, uz_r, kxd, div_r, svv_uh)
    elseif(C_shock > 0.d0 .or. istage <= 1) then
      call get_div_u(ux_r, uy_r, uz_r, kxd, div_r)
    endif
    if(C_shock > 0.d0) then
      call get_zeta(div_r, zeta_r)
      q_r = rho_r*zeta_r*div_r
    else
      zeta_r = 0.d0
      q_r    = 0.d0
    endif

    ! div_N m (w is used as a work array; its content is destroyed by btran anyway)
    !$omp parallel do private(i, k, mfilt) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          mfilt = merge(filter(i,k,j), 1.d0, filter_mom)
          w(i,k,j,1) = zi*(kxd(i,j)*mx(i,k,j) + ky_d(j)*my(i,k,j) + kz_d(k)*mz(i,k,j))*mfilt
        enddo
      enddo
    enddo
    !$omp end parallel do
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_btran_c2r(w(:,:,:,1), flx_r(:,:,:,iflx_s), 'tff')
    if (proc0) call put_time_stamp(timer_fft)

    sums(1) = -sum(flx_r(:,:,:,iflx_s))
    sums(2) =  sum(abs(flx_r(:,:,:,iflx_s)))
    sums(3) =  sum(q_r*div_r)
    call sum_allreduce(sums)
    dMdt_raw  = sums(1)*dV
    dMdt_abs  = sums(2)*dV
    visc_diss = sums(3)*dV

    if(istage <= 1) then
      ! maxs: max rho, -min rho, max zeta, -min D, max rho*zeta, max(|u_d| + cf_d) for d = x, y, z
      maxs(1) =  maxval(rho_r)
      maxs(2) = -minval(rho_r)
      maxs(3) =  maxval(zeta_r)
      maxs(4) = -minval(div_r)
      maxs(5) =  maxval(rho_r*zeta_r)
      maxs(6:8) = 0.d0
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            vad2_r(1) = w_r(j,k,i,ibx)**2/rho_r(j,k,i)
            vad2_r(2) = w_r(j,k,i,iby)**2/rho_r(j,k,i)
            vad2_r(3) = w_r(j,k,i,ibz)**2/rho_r(j,k,i)
            va2_r     = sum(vad2_r)
            do d = 1, 3
              cf2_r = 0.5d0*(cs2va2 + va2_r + sqrt(max(0.d0, (cs2va2 + va2_r)**2 - 4.d0*cs2va2*vad2_r(d))))
              select case(d)
                case(1); maxs(6) = max(maxs(6), abs(ux_r(j,k,i)) + sqrt(cf2_r))
                case(2); maxs(7) = max(maxs(7), abs(uy_r(j,k,i)) + sqrt(cf2_r))
                case(3); maxs(8) = max(maxs(8), abs(uz_r(j,k,i)) + sqrt(cf2_r))
              end select
            enddo
          enddo
        enddo
      enddo
      call max_allreduce(maxs)
    endif

    ! SVV: eps_d (fixed over the step) and, for svv_form = 'entropy', the terms of rho and m
    ! and the power method for lambda_svv (one iteration per step; until converged at the check steps)
    if(istage == 1) call get_svv_eps(maxs(6:8), eps_svv)
    pi_nit   = 0
    pi_check = .false.
    if(svv_es .and. istage >= 1) then
      call get_svv_es(w_r(:,:,:,ilnrho), rho_r, ux_r, uy_r, uz_r, eps_svv, svv_r, svv_mh, svv_diss, svv_uh)
      if(istage == 1) then
        pi_check = svv_pi_check > 0 .and. istep_adv > 0 .and. mod(istep_adv, svv_pi_check) == 0
        if(pi_check) then
          call svv_power(w_r(:,:,:,ilnrho), rho_r, w_r(:,:,:,imx), w_r(:,:,:,imy), w_r(:,:,:,imz), eps_svv, &
                         svv_r, svv_mh, 200, 1.d-3, pi_nit, pi_conv)
          if(.not. pi_conv .and. proc0) &
            print '(" WARNING: the power method for lambda_svv did not converge in ", i4, " iterations")', pi_nit
        else
          call svv_power(w_r(:,:,:,ilnrho), rho_r, w_r(:,:,:,imx), w_r(:,:,:,imy), w_r(:,:,:,imz), eps_svv, &
                         svv_r, svv_mh, 1, 0.d0, pi_nit, pi_conv)
        endif
      endif
    endif

    ! CFL: advection by fast waves over the active directions, and the explicit viscous terms
    if(istage == 1) then
      dt_prev = dt

      ! When the file 'dt_adjust' including a float number is created,
      ! dt will be manually adjusted to that value while running.
      call dt_adjust_while_running()

      dt_adv = 0.d0
      if(active_x) dt_adv = dt_adv + maxs(6)/dlx
      if(active_y) dt_adv = dt_adv + maxs(7)/dly
      if(active_z) dt_adv = dt_adv + maxs(8)/dlz
      dt_adv = cfl/dt_adv

      if(svv_es) then
        ! 1/dt_visc: shock viscosity/cfl_visc + lambda_svv/cfl_svv in (rho, m), SVV in B/cfl_visc
        rate_visc = max(maxs(3)*kd2_max/cfl_visc + svv_lam/cfl_svv, sum(eps_svv*svv_lam_b)/cfl_visc)
        if(rate_visc > 0.d0) then
          dt_visc = 1.d0/rate_visc
        else
          dt_visc = huge(1.d0)
        endif
      else
        ! largest damping rate: shock viscosity + SVV in m, SVV in ln rho and in B
        rate_visc = max(maxs(3)*kd2_max + sum(eps_svv*svv_lam_m), sum(eps_svv*svv_lam_s), sum(eps_svv*svv_lam_b))
        if(rate_visc > 0.d0) then
          dt_visc = cfl_visc/rate_visc
        else
          dt_visc = huge(1.d0)
        endif
      endif

      dt_cfl = min(dt_adv, dt_visc)
      if(is_nonfinite(dt_cfl) .or. dt_cfl <= 0.d0) call abort_run('get_nonlinear_terms: invalid dt_cfl')

      call reset_dt(dt_cfl, counter)
      n = 0
      do while(dt_cfl < dt .and. n < 100)
        call reset_dt(dt_cfl, counter)
        n = n + 1
      enddo
      if(dt_cfl < dt) call abort_run('get_nonlinear_terms: dt could not be reduced below dt_cfl (check reset_method)')

      ! check of the explicit SVV step with the converged lambda_svv: dt of the previous and of this step
      if(pi_check .and. max(dt_prev, dt)*svv_lam > rk3_lim) then
        write(msg, '("get_nonlinear_terms: dt*lambda_svv = ", es12.4, " > ", f5.2, " (dt_prev = ", es12.4, &
                    &", dt = ", es12.4, ", lambda_svv = ", es12.4, ", step = ", i10, ")")') &
              max(dt_prev, dt)*svv_lam, rk3_lim, dt_prev, dt, svv_lam, istep_adv
        call abort_run(trim(msg))
      endif

      F_svv = 0.d0
      if(svv_es .and. sum(eps_svv*svv_lam_s) > 0.d0) F_svv = svv_lam/sum(eps_svv*svv_lam_s)
      if(proc0) then
        write (unit=cfl_unit, fmt="(100es26.17e3)") tt, dt, dt_adv, dt_visc, maxs(6:8), maxs(3), &
                                                  maxs(5)/max(-maxs(2)*maxs(3), tiny(1.d0)), eps_svv, &
                                                  merge(svv_lam, 0.d0, svv_es), F_svv, &
                                                  merge(dt*svv_lam, 0.d0, svv_es), dble(pi_nit)
        call flush(cfl_unit)
      endif
    endif

    ! Output of the state at tt
    if(istage <= 1) then
      if(write_intvl_cons == 0.d0 .or. (write_intvl_cons > 0.d0 .and. &
         (first .or. istage == 0 .or. write_due(write_intvl_cons)))) then
        if(proc0) then
          write (unit=cons_unit, fmt="(100es26.17e3)") tt, dt, mon_pre(imon_M), (mon_pre(imon_M) - mass0)/mass0, &
                mon_pre(imon_px:imon_pz), mon_pre(imon_ek:imon_ei), sum(mon_pre(imon_ek:imon_ei)), &
                -maxs(2), maxs(1), maxs(3), -maxs(4), visc_diss, dmass_corr, dmass_corr_tot
          call flush(cons_unit)
        endif
      endif

      if(write_intvl_1D > 0.d0 .and. (first .or. istage == 0 .or. write_due(write_intvl_1D))) call write_line
    endif

    if(istage == 0) then
      if (proc0) call put_time_stamp(timer_nonlinear_terms)
      return
    endif

    ! 4. Calculate nonlinear terms in real space
    flx_r(:,:,:,iflx_s) = -flx_r(:,:,:,iflx_s)/rho_r

    ! SVV of the mass equation: ds/dt += (d rho/dt)_svv/rho, so that sum rho ds/dt = 0 is kept
    if(svv_es) then
      ! 'entropy': (d rho/dt)_svv from get_svv_es
      svv_r = svv_r/rho_r
      flx_r(:,:,:,iflx_s) = flx_r(:,:,:,iflx_s) + svv_r
    elseif(svv_on) then
      ! 'linear': d rho/dt = F^-1(-lambda F rho) (lambda(0) = 0)
      svv_r = rho_r
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_ftran_r2c(svv_r, w(:,:,:,1), 'fft')
      if (proc0) call put_time_stamp(timer_fft)
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            w(i,k,j,1) = -(eps_svv(1)*kq2x(i) + eps_svv(2)*kq2y(j) + eps_svv(3)*kq2z(k))*w(i,k,j,1)/nlx/nly/nlz
          enddo
        enddo
      enddo
      !$omp end parallel do
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_btran_c2r(w(:,:,:,1), svv_r, 'tff')
      if (proc0) call put_time_stamp(timer_fft)
      svv_r = svv_r/rho_r
      flx_r(:,:,:,iflx_s) = flx_r(:,:,:,iflx_s) + svv_r
    endif

    !$omp parallel do private(j, k, b2h_r, ptot_r) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          b2h_r  = 0.5d0*(w_r(j,k,i,ibx)**2 + w_r(j,k,i,iby)**2 + w_r(j,k,i,ibz)**2)
          ptot_r = cs2va2*rho_r(j,k,i) + b2h_r - q_r(j,k,i)

          flx_r(j,k,i,iflx_mxx) = rho_r(j,k,i)*ux_r(j,k,i)*ux_r(j,k,i) - w_r(j,k,i,ibx)*w_r(j,k,i,ibx) + ptot_r
          flx_r(j,k,i,iflx_mxy) = rho_r(j,k,i)*ux_r(j,k,i)*uy_r(j,k,i) - w_r(j,k,i,ibx)*w_r(j,k,i,iby)
          flx_r(j,k,i,iflx_mxz) = rho_r(j,k,i)*ux_r(j,k,i)*uz_r(j,k,i) - w_r(j,k,i,ibx)*w_r(j,k,i,ibz)
          flx_r(j,k,i,iflx_myy) = rho_r(j,k,i)*uy_r(j,k,i)*uy_r(j,k,i) - w_r(j,k,i,iby)*w_r(j,k,i,iby) + ptot_r
          flx_r(j,k,i,iflx_myz) = rho_r(j,k,i)*uy_r(j,k,i)*uz_r(j,k,i) - w_r(j,k,i,iby)*w_r(j,k,i,ibz)
          flx_r(j,k,i,iflx_mzz) = rho_r(j,k,i)*uz_r(j,k,i)*uz_r(j,k,i) - w_r(j,k,i,ibz)*w_r(j,k,i,ibz) + ptot_r

          flx_r(j,k,i,iflx_bx ) = w_r(j,k,i,iby)*uz_r(j,k,i) - w_r(j,k,i,ibz)*uy_r(j,k,i)
          flx_r(j,k,i,iflx_by ) = w_r(j,k,i,ibz)*ux_r(j,k,i) - w_r(j,k,i,ibx)*uz_r(j,k,i)
          flx_r(j,k,i,iflx_bz ) = w_r(j,k,i,ibx)*uy_r(j,k,i) - w_r(j,k,i,iby)*ux_r(j,k,i)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 5. Forward FFT
    if (proc0) call put_time_stamp(timer_fft)
    ! for some reason r2c_many doesn't work at Fugaku
    do i = 1, nftran
      call p3dfft_ftran_r2c(flx_r(:,:,:,i), flx(:,:,:,i), 'fft')
    enddo
    if (proc0) call put_time_stamp(timer_fft)
    flx = flx/nlx/nly/nlz

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Stop if a field is non-finite or
!!          |ln rho| > s_abort (no correction)
!-----------------------------------------------!
  subroutine check_state(istage)
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: s_abort
    use mp, only: proc_id
    use time, only: tt
    use shock_viscosity, only: is_nonfinite, abort_run
    implicit none
    integer, intent(in) :: istage
    character(len=5), parameter :: fname(nfields) = (/'lnrho', 'mx   ', 'my   ', 'mz   ', 'bx   ', 'by   ', 'bz   '/)
    integer :: i, j, k, n

    if(.not. (any(is_nonfinite(w_r)) .or. any(abs(w_r(:,:,:,ilnrho)) > s_abort))) return

    do n = 1, nfields
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            if(is_nonfinite(w_r(j,k,i,n)) .or. (n == ilnrho .and. abs(w_r(j,k,i,n)) > s_abort)) then
              write(*, '("[proc ", i6, "] invalid state at t = ", es24.16, ", step ", i10, ", stage ", i2, &
                       &": ", a, "(iy, iz, ix) = (", i8, ",", i8, ",", i8, ") = ", es24.16)') &
                    proc_id, tt, istep_adv + 1, istage, trim(fname(n)), j, k, i, w_r(j,k,i,n)
              call abort_run('non-finite value or |ln rho| > s_abort')
            endif
          enddo
        enddo
      enddo
    enddo

  end subroutine check_state


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Volume integrals of the state in
!!          w_r and rho_r (grid sums times dV)
!-----------------------------------------------!
  subroutine get_monitor_sums(mon)
    use grid, only: dlx, dly, dlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: cs2va2
    use mp, only: sum_allreduce
    implicit none
    real(r8), intent(out) :: mon(nmon)
    integer :: i, j, k

    mon = 0.d0
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          mon(imon_M ) = mon(imon_M ) + rho_r(j,k,i)
          mon(imon_px) = mon(imon_px) + w_r(j,k,i,imx)
          mon(imon_py) = mon(imon_py) + w_r(j,k,i,imy)
          mon(imon_pz) = mon(imon_pz) + w_r(j,k,i,imz)
          mon(imon_ek) = mon(imon_ek) + 0.5d0*(w_r(j,k,i,imx)**2 + w_r(j,k,i,imy)**2 + w_r(j,k,i,imz)**2)/rho_r(j,k,i)
          mon(imon_em) = mon(imon_em) + 0.5d0*(w_r(j,k,i,ibx)**2 + w_r(j,k,i,iby)**2 + w_r(j,k,i,ibz)**2)
          mon(imon_ei) = mon(imon_ei) + cs2va2*rho_r(j,k,i)*(w_r(j,k,i,ilnrho) - 1.d0)
        enddo
      enddo
    enddo
    call sum_allreduce(mon)
    mon = mon*dlx*dly*dlz

  end subroutine get_monitor_sums


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Same criterion as check_write_now
!!          in calliope.F90
!-----------------------------------------------!
  function write_due(intvl)
    use time, only: tt, dt
    implicit none
    real(r8), intent(in) :: intvl
    logical  :: write_due
    integer, parameter :: i16 = selected_int_kind(16)

    write_due = floor(tt/intvl, i16) /= floor((tt+dt)/intvl, i16)

  end function write_due


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Add a to the k = 0 mode of f
!!          (uniform shift in real space)
!-----------------------------------------------!
  subroutine shift_mean(f, a)
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: f
    real(r8), intent(in) :: a
    integer :: i, j, k

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) f(i,k,j) = f(i,k,j) + a
        enddo
      enddo
    enddo

  end subroutine shift_mean


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Multiply each field by the dealiasing
!!          filter if filter_* of the field is T
!-----------------------------------------------!
  subroutine apply_filter(lnrho, mx, my, mz, bx, by, bz)
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use dealias, only: filter
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: lnrho, mx, my, mz, bx, by, bz
    integer :: i, j, k

    !$omp parallel do private(j, i) schedule(static)
    do k = ikz_st, ikz_en
      do j = iky_st, iky_en
        do i = ikx_st, ikx_en
          lnrho(i,k,j) = lnrho(i,k,j)*merge(filter(i,k,j), 1.d0, filt_on(ilnrho))
          mx   (i,k,j) = mx   (i,k,j)*merge(filter(i,k,j), 1.d0, filt_on(imx   ))
          my   (i,k,j) = my   (i,k,j)*merge(filter(i,k,j), 1.d0, filt_on(imy   ))
          mz   (i,k,j) = mz   (i,k,j)*merge(filter(i,k,j), 1.d0, filt_on(imz   ))
          bx   (i,k,j) = bx   (i,k,j)*merge(filter(i,k,j), 1.d0, filt_on(ibx   ))
          by   (i,k,j) = by   (i,k,j)*merge(filter(i,k,j), 1.d0, filt_on(iby   ))
          bz   (i,k,j) = bz   (i,k,j)*merge(filter(i,k,j), 1.d0, filt_on(ibz   ))
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine apply_filter


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Explicit terms of all modes from flx
!!          (filled by get_nonlinear_terms) and
!!          the SVV of m and B: -lambda(k) B^, and
!!          svv_mh ('entropy') or -lambda(k) m^
!!          ('linear')
!-----------------------------------------------!
  subroutine get_explicit_terms(mx, my, mz, bx, by, bz)
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use force, only: fmx, fmy, fmz
    use shock_viscosity, only: ky_d, kz_d
    use shock_viscosity, only: svv_on, svv_es, kq2x, kq2y, kq2z
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: mx, my, mz, bx, by, bz
    real(r8) :: lam
    integer :: i, j, k

    !$omp parallel do private(j, k, i, lam)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          call get_ext_terms(exp_terms(i,k,j,:), &
                             mx(i,k,j), my(i,k,j), mz(i,k,j), bx(i,k,j), by(i,k,j), bz(i,k,j), &
                             flx(i,k,j,:), &
                             fmx(i,k,j), fmy(i,k,j), fmz(i,k,j), &
                             kxd(i,j), ky_d(j), kz_d(k))
          if(svv_on) then
            lam = eps_svv(1)*kq2x(i) + eps_svv(2)*kq2y(j) + eps_svv(3)*kq2z(k)
            if(svv_es) then
              exp_terms(i,k,j,imx) = exp_terms(i,k,j,imx) + svv_mh(i,k,j,1)
              exp_terms(i,k,j,imy) = exp_terms(i,k,j,imy) + svv_mh(i,k,j,2)
              exp_terms(i,k,j,imz) = exp_terms(i,k,j,imz) + svv_mh(i,k,j,3)
            else
              exp_terms(i,k,j,imx) = exp_terms(i,k,j,imx) - lam*mx(i,k,j)
              exp_terms(i,k,j,imy) = exp_terms(i,k,j,imy) - lam*my(i,k,j)
              exp_terms(i,k,j,imz) = exp_terms(i,k,j,imz) - lam*mz(i,k,j)
            endif
            exp_terms(i,k,j,ibx) = exp_terms(i,k,j,ibx) - lam*bx(i,k,j)
            exp_terms(i,k,j,iby) = exp_terms(i,k,j,iby) - lam*by(i,k,j)
            exp_terms(i,k,j,ibz) = exp_terms(i,k,j,ibz) - lam*bz(i,k,j)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine get_explicit_terms


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms, &
                           mx, my, mz, bx, by, bz, &
                           flx, &
                           fmx, fmy, fmz, &
                           kxd, kyd, kzd)
    use params, only: zi, q
    use shearing_box, only: shear_flg
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: mx, my, mz, bx, by, bz
    complex(r8), intent(in ) :: flx(nftran)
    complex(r8), intent(in ) :: fmx, fmy, fmz
    real(r8)   , intent(in)  :: kxd, kyd, kzd
    complex(r8)              :: nl(nfields)

    ! -(div_N m)/rho (computed in real space)
    nl(ilnrho) = flx(iflx_s)

    ! div {m u - B B + [cs^2 rho + B^2/2 - q]*I}
    nl(imx ) = -zi*( kxd*flx(iflx_mxx) + kyd*flx(iflx_mxy) + kzd*flx(iflx_mxz) )
    nl(imy ) = -zi*( kxd*flx(iflx_mxy) + kyd*flx(iflx_myy) + kzd*flx(iflx_myz) )
    nl(imz ) = -zi*( kxd*flx(iflx_mxz) + kyd*flx(iflx_myz) + kzd*flx(iflx_mzz) )

    ! curl (b x u)
    nl(ibx ) = -zi*( kyd*flx(iflx_bz) - kzd*flx(iflx_by) )
    nl(iby ) = -zi*( kzd*flx(iflx_bx) - kxd*flx(iflx_bz) )
    nl(ibz ) = -zi*( kxd*flx(iflx_by) - kyd*flx(iflx_bx) )

    exp_terms(ilnrho) = nl(ilnrho)

    exp_terms(imx ) = nl(imx) + fmx + 2.d0*shear_flg*my
    exp_terms(imy ) = nl(imy) + fmy - (2.d0 - q)*shear_flg*mx
    exp_terms(imz ) = nl(imz) + fmz

    exp_terms(ibx ) = nl(ibx)
    exp_terms(iby ) = nl(iby) - q*shear_flg*bx
    exp_terms(ibz ) = nl(ibz)

  end subroutine get_ext_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Time integral of hyperdissipation
!-----------------------------------------------!
  subroutine get_imp_terms_tintg(imp_terms_tintg, t, kx, ky, kz, coeff, coeff_h, nexp)
    use grid, only: k2_max
    use params, only: shear
    use shearing_box, only: get_imp_terms_tintg_with_shear
    implicit none
    real(r8), intent(out) :: imp_terms_tintg
    real(r8), intent(in) :: t, kx, ky, kz, coeff, coeff_h
    integer, intent(in) :: nexp

    if(shear) then
      call get_imp_terms_tintg_with_shear(imp_terms_tintg, t, kx, ky, kz, coeff, coeff_h, nexp )
    else
      imp_terms_tintg = -(coeff*((kx**2 + ky**2 + kz**2)/k2_max) + coeff_h*((kx**2 + ky**2 + kz**2)/k2_max)**nexp)*t
    endif
  end subroutine get_imp_terms_tintg


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Energy and mass rates of one RK stage
!!          dE/dt = sum [dE/ds Pds + u.Pdm + B.PdB] dV
!!          with E = m^2/(2 rho) + B^2/2
!!                 + cs^2 rho (ln rho - 1),
!!          dE/ds = rho (cs^2 s - u^2/2), and P the
!!          filter applied after the stage.
!!          eps_shock = sum u.F^-1(P i k_d q^) dV;
!!          R_zeta = eps_shock + sum rho zeta D^2 dV
!!          vanishes in exact arithmetic.
!!          eps_svv, eps_svv_b: energy rates of the SVV
!!          of (rho, m) and of B (included in eps_N).
!!          eps_svv_b <= 0. 'linear': eps_svv is not
!!          sign-definite. 'entropy': R_svv = eps_svv
!!          + D_svv vanishes in exact arithmetic and
!!          D_svv >= 0 (a sum of squares). A_svv =
!!          sum |dE/dU.U_svv| dV bounds the rounding
!!          error of eps_svv.
!!          The linear terms (integrating factor) are
!!          evaluated with the stage state.
!-----------------------------------------------!
  subroutine stage_budget(istage, lnrho, mx, my, mz, bx, by, bz)
    use grid, only: dlx, dly, dlz, nlx, nly, nlz, k2_max
    use grid, only: ky, kz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: zi, cs2va2, C_shock, filter_mom
    use dealias, only: filter
    use force, only: driven, fmx, fmy, fmz
    use mp, only: sum_allreduce
    use time, only: tt, dt
    use time_stamp, only: put_time_stamp, timer_fft
    use mp, only: proc0
    use shock_viscosity, only: ky_d, kz_d
    use shock_viscosity, only: svv_on, svv_es, kq2x, kq2y, kq2z
    implicit none
    integer, intent(in) :: istage
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: lnrho, mx, my, mz, bx, by, bz
    complex(r8), allocatable, dimension(:,:,:) :: wk, qh
    real   (r8), allocatable, dimension(:,:,:) :: wk_r, dEds_r
    real(r8) :: loc(nbgt), kd, mfilt, k2r, lrate
    integer  :: i, j, k, n, d

    allocate(wk  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(wk_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))
    allocate(dEds_r, mold=wk_r)

    dEds_r = rho_r*(cs2va2*w_r(:,:,:,ilnrho) - 0.5d0*(ux_r**2 + uy_r**2 + uz_r**2))

    loc = 0.d0

    ! Rates of the explicit terms after the filter of each field
    do n = 1, nfields
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            wk(i,k,j) = exp_terms(i,k,j,n)*merge(filter(i,k,j), 1.d0, filt_on(n))
          enddo
        enddo
      enddo
      call btran_wk
      select case(n)
        case(ilnrho)
          loc(ib_eps_N    ) = loc(ib_eps_N) + sum(dEds_r*wk_r)
          loc(ib_dMdt_proj) = sum(rho_r*wk_r)
        case(imx); loc(ib_eps_N) = loc(ib_eps_N) + sum(ux_r*wk_r)
        case(imy); loc(ib_eps_N) = loc(ib_eps_N) + sum(uy_r*wk_r)
        case(imz); loc(ib_eps_N) = loc(ib_eps_N) + sum(uz_r*wk_r)
        case(ibx); loc(ib_eps_N) = loc(ib_eps_N) + sum(w_r(:,:,:,ibx)*wk_r)
        case(iby); loc(ib_eps_N) = loc(ib_eps_N) + sum(w_r(:,:,:,iby)*wk_r)
        case(ibz); loc(ib_eps_N) = loc(ib_eps_N) + sum(w_r(:,:,:,ibz)*wk_r)
      end select
    enddo

    ! Shock viscosity: u.F^-1(P i k_d q^)
    if(C_shock > 0.d0) then
      allocate(qh, mold=wk)
      wk_r = q_r
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_ftran_r2c(wk_r, qh, 'fft')
      if (proc0) call put_time_stamp(timer_fft)
      qh = qh/nlx/nly/nlz
      do d = 1, 3
        do j = iky_st, iky_en
          do k = ikz_st, ikz_en
            do i = ikx_st, ikx_en
              select case(d)
                case(1); kd = kxd(i,j)
                case(2); kd = ky_d(j)
                case(3); kd = kz_d(k)
              end select
              mfilt = merge(filter(i,k,j), 1.d0, filter_mom)
              wk(i,k,j) = zi*kd*qh(i,k,j)*mfilt
            enddo
          enddo
        enddo
        call btran_wk
        select case(d)
          case(1); loc(ib_eps_shock) = loc(ib_eps_shock) + sum(ux_r*wk_r)
          case(2); loc(ib_eps_shock) = loc(ib_eps_shock) + sum(uy_r*wk_r)
          case(3); loc(ib_eps_shock) = loc(ib_eps_shock) + sum(uz_r*wk_r)
        end select
      enddo
      deallocate(qh)
    endif

    ! SVV: dE/ds.F^-1(P F svv_r) + u.F^-1(P svv_m^) (eps_svv) + B.F^-1(P(-lambda B^)) (eps_svv_b),
    ! svv_m^ = svv_mh ('entropy') or -lambda m^ ('linear')
    if(svv_on) then
      wk_r = svv_r
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_ftran_r2c(wk_r, wk, 'fft')
      if (proc0) call put_time_stamp(timer_fft)
      wk = wk/nlx/nly/nlz*merge(filter, 1.d0, filt_on(ilnrho))
      call btran_wk
      loc(ib_eps_svv) = loc(ib_eps_svv) + sum(dEds_r*wk_r)
      loc(ib_A_svv  ) = loc(ib_A_svv  ) + sum(abs(dEds_r*wk_r))
      do n = imx, ibz
        do j = iky_st, iky_en
          do k = ikz_st, ikz_en
            do i = ikx_st, ikx_en
              lrate = (eps_svv(1)*kq2x(i) + eps_svv(2)*kq2y(j) + eps_svv(3)*kq2z(k)) &
                     *merge(filter(i,k,j), 1.d0, filt_on(n))
              select case(n)
                case(imx:imz)
                  if(svv_es) then
                    wk(i,k,j) = svv_mh(i,k,j,n-imx+1)*merge(filter(i,k,j), 1.d0, filt_on(n))
                  else
                    select case(n)
                      case(imx); wk(i,k,j) = -lrate*mx(i,k,j)
                      case(imy); wk(i,k,j) = -lrate*my(i,k,j)
                      case(imz); wk(i,k,j) = -lrate*mz(i,k,j)
                    end select
                  endif
                case(ibx); wk(i,k,j) = -lrate*bx(i,k,j)
                case(iby); wk(i,k,j) = -lrate*by(i,k,j)
                case(ibz); wk(i,k,j) = -lrate*bz(i,k,j)
              end select
            enddo
          enddo
        enddo
        call btran_wk
        select case(n)
          case(imx); loc(ib_eps_svv  ) = loc(ib_eps_svv  ) + sum(ux_r*wk_r); loc(ib_A_svv) = loc(ib_A_svv) + sum(abs(ux_r*wk_r))
          case(imy); loc(ib_eps_svv  ) = loc(ib_eps_svv  ) + sum(uy_r*wk_r); loc(ib_A_svv) = loc(ib_A_svv) + sum(abs(uy_r*wk_r))
          case(imz); loc(ib_eps_svv  ) = loc(ib_eps_svv  ) + sum(uz_r*wk_r); loc(ib_A_svv) = loc(ib_A_svv) + sum(abs(uz_r*wk_r))
          case(ibx); loc(ib_eps_svv_b) = loc(ib_eps_svv_b) + sum(w_r(:,:,:,ibx)*wk_r)
          case(iby); loc(ib_eps_svv_b) = loc(ib_eps_svv_b) + sum(w_r(:,:,:,iby)*wk_r)
          case(ibz); loc(ib_eps_svv_b) = loc(ib_eps_svv_b) + sum(w_r(:,:,:,ibz)*wk_r)
        end select
      enddo
      ! svv_diss is a local sum (get_svv_es); it is reduced with loc
      if(svv_es) loc(ib_D_svv) = svv_diss
    endif

    ! Forcing: u.F^-1(P f^)
    if(driven) then
      wk = fmx*merge(filter, 1.d0, filter_mom); call btran_wk; loc(ib_eps_force) = loc(ib_eps_force) + sum(ux_r*wk_r)
      wk = fmy*merge(filter, 1.d0, filter_mom); call btran_wk; loc(ib_eps_force) = loc(ib_eps_force) + sum(uy_r*wk_r)
      wk = fmz*merge(filter, 1.d0, filter_mom); call btran_wk; loc(ib_eps_force) = loc(ib_eps_force) + sum(uz_r*wk_r)
    endif

    ! Linear dissipation: -L(k) w^ with L = coeff*k^2/k2_max + coeff_h*(k^2/k2_max)^nexp
    do n = 1, nfields
      if(dcoef(n) == 0.d0 .and. dcoef_h(n) == 0.d0) cycle
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            k2r   = (kxt(i,j)**2 + ky(j)**2 + kz(k)**2)/k2_max
            lrate = dcoef(n)*k2r + dcoef_h(n)*k2r**dcoef_exp(n)
            select case(n)
              case(ilnrho); wk(i,k,j) = -lrate*lnrho(i,k,j)
              case(imx   ); wk(i,k,j) = -lrate*mx   (i,k,j)
              case(imy   ); wk(i,k,j) = -lrate*my   (i,k,j)
              case(imz   ); wk(i,k,j) = -lrate*mz   (i,k,j)
              case(ibx   ); wk(i,k,j) = -lrate*bx   (i,k,j)
              case(iby   ); wk(i,k,j) = -lrate*by   (i,k,j)
              case(ibz   ); wk(i,k,j) = -lrate*bz   (i,k,j)
            end select
          enddo
        enddo
      enddo
      call btran_wk
      select case(n)
        case(ilnrho); loc(ib_eps_lmd) = loc(ib_eps_lmd) + sum(dEds_r*wk_r)
        case(imx   ); loc(ib_eps_nu ) = loc(ib_eps_nu ) + sum(ux_r*wk_r)
        case(imy   ); loc(ib_eps_nu ) = loc(ib_eps_nu ) + sum(uy_r*wk_r)
        case(imz   ); loc(ib_eps_nu ) = loc(ib_eps_nu ) + sum(uz_r*wk_r)
        case(ibx   ); loc(ib_eps_res) = loc(ib_eps_res) + sum(w_r(:,:,:,ibx)*wk_r)
        case(iby   ); loc(ib_eps_res) = loc(ib_eps_res) + sum(w_r(:,:,:,iby)*wk_r)
        case(ibz   ); loc(ib_eps_res) = loc(ib_eps_res) + sum(w_r(:,:,:,ibz)*wk_r)
      end select
    enddo

    call sum_allreduce(loc)
    loc = loc*dlx*dly*dlz

    ! already reduced in get_nonlinear_terms
    loc(ib_dMdt_raw) = dMdt_raw
    loc(ib_dMdt_abs) = dMdt_abs
    loc(ib_visc    ) = visc_diss
    loc(ib_R_zeta  ) = loc(ib_eps_shock) + visc_diss
    if(svv_es) loc(ib_R_svv) = loc(ib_eps_svv) + loc(ib_D_svv)

    if(istage == 1) then
      t_bgt  = tt
      dt_bgt = dt
    endif
    bgt = bgt + b_rk(istage)*loc
    if(visc_diss > 0.d0) R_zeta_ratio_max = max(R_zeta_ratio_max, abs(loc(ib_R_zeta))/visc_diss)
    if(loc(ib_D_svv) > 0.d0) R_svv_ratio_max = max(R_svv_ratio_max, abs(loc(ib_R_svv))/loc(ib_D_svv))
    if(svv_es .and. loc(ib_A_svv) > 0.d0) R_svv_scale_max = max(R_svv_scale_max, abs(loc(ib_R_svv))/loc(ib_A_svv))
    if(istage == 3) bgt_pending = .true.

    deallocate(wk, wk_r, dEds_r)

  contains

    subroutine btran_wk
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_btran_c2r(wk, wk_r, 'tff')
      if (proc0) call put_time_stamp(timer_fft)
    end subroutine btran_wk

  end subroutine stage_budget


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Write the budget of the step that has
!!          just finished. mon_np1: state at the
!!          end of the step (before the mass
!!          correction). E_n = mon_n.
!-----------------------------------------------!
  subroutine write_budget(mon_np1)
    use mp, only: proc0
    implicit none
    real(r8), intent(in) :: mon_np1(nmon)
    real(r8) :: E_n, E_np1, dMdt_act, dEdt_act, dEdt_sum, eps_ideal

    E_n   = sum(mon_n  (imon_ek:imon_ei))
    E_np1 = sum(mon_np1(imon_ek:imon_ei))
    dMdt_act  = (mon_np1(imon_M) - mon_n(imon_M))/dt_bgt
    dEdt_act  = (E_np1 - E_n)/dt_bgt
    dEdt_sum  = bgt(ib_eps_N) + bgt(ib_eps_res) + bgt(ib_eps_nu) + bgt(ib_eps_lmd)
    eps_ideal = bgt(ib_eps_N) - bgt(ib_eps_shock) - bgt(ib_eps_force) - bgt(ib_eps_svv) - bgt(ib_eps_svv_b)

    if(proc0) then
      write (unit=budget_unit, fmt="(100es26.17e3)") t_bgt, dt_bgt, &
            dMdt_act, bgt(ib_dMdt_raw), bgt(ib_dMdt_proj), bgt(ib_dMdt_abs), &
            dEdt_act, dEdt_sum, dEdt_act - dEdt_sum, &
            eps_ideal, bgt(ib_eps_shock), -bgt(ib_visc), bgt(ib_R_zeta), R_zeta_ratio_max, &
            bgt(ib_eps_force), bgt(ib_eps_res), bgt(ib_eps_nu), bgt(ib_eps_lmd), E_n, E_np1, &
            bgt(ib_eps_svv) + bgt(ib_eps_svv_b), bgt(ib_eps_svv_b), -bgt(ib_D_svv), bgt(ib_R_svv), R_svv_ratio_max, &
            R_svv_scale_max
      call flush(budget_unit)
    endif

    bgt = 0.d0
    R_zeta_ratio_max = 0.d0
    R_svv_ratio_max  = 0.d0
    R_svv_scale_max  = 0.d0
    bgt_pending = .false.

  end subroutine write_budget


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   ASCII output of rho, u, B, zeta, D on
!!          the line along the first active
!!          direction through grid index 1 of the
!!          other directions (out1d/line_NNNNN.dat)
!-----------------------------------------------!
  subroutine write_line
    use grid, only: nlx, nly, nlz, xx, yy, zz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use mp, only: proc0, sum_reduce
    use time, only: tt
    use file, only: get_unused_unit
    use shock_viscosity, only: active_x, active_y, active_z, abort_run
    implicit none
    real(r8), allocatable :: buf(:,:)
    character(len=1)   :: dir
    character(len=100) :: fn
    integer :: i, j, k, nl, unit, ierr

    if(active_x) then
      dir = 'x'; nl = nlx
    elseif(active_y) then
      dir = 'y'; nl = nly
    else
      dir = 'z'; nl = nlz
    endif
    allocate(buf(nl, 10), source=0.d0)

    select case(dir)
      case('x')
        if(ily_st <= 1 .and. 1 <= ily_en .and. ilz_st <= 1 .and. 1 <= ilz_en) then
          do i = ilx_st, ilx_en
            buf(i,:) = (/xx(i), rho_r(1,1,i), ux_r(1,1,i), uy_r(1,1,i), uz_r(1,1,i), &
                         w_r(1,1,i,ibx), w_r(1,1,i,iby), w_r(1,1,i,ibz), zeta_r(1,1,i), div_r(1,1,i)/)
          enddo
        endif
      case('y')
        if(ilz_st <= 1 .and. 1 <= ilz_en .and. ilx_st <= 1 .and. 1 <= ilx_en) then
          do j = ily_st, ily_en
            buf(j,:) = (/yy(j), rho_r(j,1,1), ux_r(j,1,1), uy_r(j,1,1), uz_r(j,1,1), &
                         w_r(j,1,1,ibx), w_r(j,1,1,iby), w_r(j,1,1,ibz), zeta_r(j,1,1), div_r(j,1,1)/)
          enddo
        endif
      case('z')
        if(ily_st <= 1 .and. 1 <= ily_en .and. ilx_st <= 1 .and. 1 <= ilx_en) then
          do k = ilz_st, ilz_en
            buf(k,:) = (/zz(k), rho_r(1,k,1), ux_r(1,k,1), uy_r(1,k,1), uz_r(1,k,1), &
                         w_r(1,k,1,ibx), w_r(1,k,1,iby), w_r(1,k,1,ibz), zeta_r(1,k,1), div_r(1,k,1)/)
          enddo
        endif
    end select

    ! each point is owned by one rank; the other ranks add 0, which is exact
    call sum_reduce(buf, 0)

    if(proc0) then
      write(fn, '("out1d/line_", i5.5, ".dat")') nline1d
      call get_unused_unit(unit)
      open(unit=unit, file=trim(fn), status='replace', iostat=ierr)
      if(ierr /= 0) call abort_run('write_line: cannot open '//trim(fn)//' (create the directory out1d)')
      write(unit, '("# t = ", es26.17e3, ", dir = ", a1, ", step = ", i10)') tt, dir, istep_adv
      write(unit, '(a)') '# 1:'//dir//' 2:rho 3:ux 4:uy 5:uz 6:bx 7:by 8:bz 9:zeta 10:D'
      do i = 1, nl
        write(unit, '(10es26.17e3)') buf(i,:)
      enddo
      close(unit)
    endif
    nline1d = nline1d + 1

    deallocate(buf)

  end subroutine write_line


!-----------------------------------------------!
!> @author  YK
!! @date    15 Jul 2021
!! @brief   Output series modes
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: lnrho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, sum_reduce
    use time, only: tt
    use diagnostics_common, only: n_series_modes, series_modes
    use diagnostics_common, only: series_modes_unit
    implicit none
    complex(r8), dimension(n_series_modes) :: lnrho_modes
    complex(r8), dimension(n_series_modes) ::    mx_modes, my_modes, mz_modes
    complex(r8), dimension(n_series_modes) ::    bx_modes, by_modes, bz_modes
    integer :: n, i, j, k

    lnrho_modes(:) = 0.d0
       mx_modes(:) = 0.d0
       my_modes(:) = 0.d0
       mz_modes(:) = 0.d0
       bx_modes(:) = 0.d0
       by_modes(:) = 0.d0
       bz_modes(:) = 0.d0

    do n = 1, n_series_modes
      i = series_modes(n, 1)
      j = series_modes(n, 2)
      k = series_modes(n, 3)

      if(       (i >= ikx_st .and. i <= ikx_en) &
          .and. (j >= iky_st .and. j <= iky_en) &
          .and. (k >= ikz_st .and. k <= ikz_en) &
        ) then

        lnrho_modes(n) = lnrho(i, k, j)
           mx_modes(n) =    mx(i, k, j)
           my_modes(n) =    my(i, k, j)
           mz_modes(n) =    mz(i, k, j)

           bx_modes(n) =    bx(i, k, j)
           by_modes(n) =    by(i, k, j)
           bz_modes(n) =    bz(i, k, j)

      endif
    enddo

    call sum_reduce(lnrho_modes, 0)
    call sum_reduce(   mx_modes, 0)
    call sum_reduce(   my_modes, 0)
    call sum_reduce(   mz_modes, 0)
    call sum_reduce(   bx_modes, 0)
    call sum_reduce(   by_modes, 0)
    call sum_reduce(   bz_modes, 0)

    do n = 1, n_series_modes
      if(proc0) then
        i = series_modes(n, 1)
        j = series_modes(n, 2)
        k = series_modes(n, 3)
999 format(es30.21, A6, 5es30.21e3)
        write (unit=series_modes_unit, fmt=999) tt, 'lnrho', kx(i), ky(j), kz(k), lnrho_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,    'mx', kx(i), ky(j), kz(k),    mx_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,    'my', kx(i), ky(j), kz(k),    my_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,    'mz', kx(i), ky(j), kz(k),    mz_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,    'bx', kx(i), ky(j), kz(k),    bx_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,    'by', kx(i), ky(j), kz(k),    by_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,    'bz', kx(i), ky(j), kz(k),    bz_modes(n)
        call flush(series_modes_unit)

      endif
    enddo

  end subroutine output_series_modes


!-----------------------------------------------!
!> @author  YK
!! @date    3 Mar 2021
!! @brief   Remap
!-----------------------------------------------!
  subroutine remap
    use fields, only: lnrho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp
    use shearing_box, only: tsc, nremap, k2t, k2t_inv
    use shearing_box, only: timer_remap, remap_fld
    implicit none
    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_remap)

    if(proc0) then
      print *
      print *, 'remapping...'
      print *
    endif

    call remap_fld(mx)
    call remap_fld(my)
    call remap_fld(mz)
    call remap_fld(bx)
    call remap_fld(by)
    call remap_fld(bz)
    call remap_fld(lnrho)

    tsc = 0.d0
    nremap = nremap + 1
    counter = 1

    ! Dealiasing
    call apply_filter(lnrho, mx, my, mz, bx, by, bz)

    call set_kx(tsc)

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          k2t(i, k, j) = kx(i)**2 + ky(j)**2 + kz(k)**2

          if(k2t(i, k, j) == 0.d0) then
            k2t_inv(i, k, j) = 0.d0
          else
            k2t_inv(i, k, j) = 1.d0/k2t(i, k, j)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

    if (proc0) call put_time_stamp(timer_remap)
  end subroutine remap


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Enforce div free for (wx, wy, wz)
!-----------------------------------------------!
  subroutine div_free(wx, wy, wz)
    use grid, only: ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: zi
    use shearing_box, only: k2t_inv
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: wx, wy, wz
    complex(r8), allocatable, dimension(:,:,:)   :: nbl2inv_div_w ! nabla^-2 (div w)
    integer :: i, j, k

    allocate(nbl2inv_div_w(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          nbl2inv_div_w(i,k,j) = -zi*(   kxt(i,j)*wx(i,k,j) &
                                       + ky(j)   *wy(i,k,j) &
                                       + kz(k)   *wz(i,k,j) )*k2t_inv(i,k,j)

          wx(i,k,j) = wx(i,k,j) - zi*kxt(i,j)*nbl2inv_div_w(i,k,j)
          wy(i,k,j) = wy(i,k,j) - zi*ky(j)   *nbl2inv_div_w(i,k,j)
          wz(i,k,j) = wz(i,k,j) - zi*kz(k)   *nbl2inv_div_w(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    deallocate(nbl2inv_div_w)

  end subroutine div_free


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Allocates arrays for this module
!-----------------------------------------------!
  subroutine allocate_advance
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use shock_viscosity, only: svv_es
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src
    real   (r8), allocatable, dimension(:,:,:) :: src_r

    if (is_allocated) then
      return
    else
      is_allocated = .true.
    endif

    allocate(src  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    allocate(src_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    allocate(lnrho_new, source=src)
    allocate(   mx_new, source=src)
    allocate(   my_new, source=src)
    allocate(   mz_new, source=src)
    allocate(   bx_new, source=src)
    allocate(   by_new, source=src)
    allocate(   bz_new, source=src)

    allocate(w        (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))
    allocate(flx      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran ), source=(0.d0, 0.d0))
    allocate(w_r      (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nfields), source=0.d0)
    allocate(flx_r    (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nftran ), source=0.d0)
    allocate(exp_terms(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))

    allocate(rho_r , source=src_r)
    allocate(ux_r  , source=src_r)
    allocate(uy_r  , source=src_r)
    allocate(uz_r  , source=src_r)
    allocate(div_r , source=src_r)
    allocate(zeta_r, source=src_r)
    allocate(q_r   , source=src_r)
    allocate(svv_r , source=src_r)
    if(svv_es) then
      allocate(svv_mh(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, 3), source=(0.d0, 0.d0))
      allocate(svv_uh(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, 3), source=(0.d0, 0.d0))
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    allocate(lnrho_tmp, source=src)
    allocate(   mx_tmp, source=src)
    allocate(   my_tmp, source=src)
    allocate(   mz_tmp, source=src)
    allocate(   bx_tmp, source=src)
    allocate(   by_tmp, source=src)
    allocate(   bz_tmp, source=src)

    allocate(exp_terms0(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms0  = 0.d0

    deallocate(src, src_r)

  end subroutine allocate_advance


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Deallocates arrays for this module
!-----------------------------------------------!
  subroutine deallocate_advance
    implicit none

    if (is_allocated) then
      is_allocated = .false.
    else
      return
    endif

    deallocate(lnrho_new)
    deallocate(   mx_new)
    deallocate(   my_new)
    deallocate(   mz_new)
    deallocate(   bx_new)
    deallocate(   by_new)
    deallocate(   bz_new)

    deallocate(w        )
    deallocate(flx      )
    deallocate(w_r      )
    deallocate(flx_r    )
    deallocate(exp_terms)

    deallocate(rho_r )
    deallocate(ux_r  )
    deallocate(uy_r  )
    deallocate(uz_r  )
    deallocate(div_r )
    deallocate(zeta_r)
    deallocate(q_r   )
    deallocate(svv_r )
    if(allocated(svv_mh)) deallocate(svv_mh)
    if(allocated(svv_uh)) deallocate(svv_uh)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    deallocate(lnrho_tmp)
    deallocate(   mx_tmp)
    deallocate(   my_tmp)
    deallocate(   mz_tmp)
    deallocate(   bx_tmp)
    deallocate(   by_tmp)
    deallocate(   bz_tmp)

    deallocate(exp_terms0)

  end subroutine deallocate_advance

end module advance

