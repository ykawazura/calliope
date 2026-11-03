!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Local shock viscosity for MHD_COMP_ISOTH
!!
!!          zeta = C_shock*Delta^2*S[max(0, -D)],  D = div_N u
!!
!!          S = (3-point local max)^shock_nmax followed by
!!              ([1,2,1]/4 smoothing)^shock_nsmooth.
!!          Each pass is applied direction by direction in the fixed
!!          order y, z, x, only in the active directions (nl > 2).
!!          Both operations have non-negative weights, so zeta >= 0,
!!          and zeta = 0 farther than shock_nmax + shock_nsmooth cells
!!          from compression.
!!
!!          The halo (width 1) is exchanged with MPI_Sendrecv.
!!          At every grid point each pass uses the same operands in the
!!          same order, irrespective of the MPI decomposition, so zeta
!!          for a given D does not depend on the decomposition (bitwise).
!!
!!          Derivative: div_N = F^-1 i k_d M F, where k_d is set to 0 at
!!          the Nyquist mode (skew-adjointness) and M is the dealiasing
!!          filter applied to the momentum (identity if filter_mom = F).
!!
!!          Spectral vanishing viscosity (Tadmor 1989; Maday, Ould Kaber
!!          & Tadmor 1993) is also set up here. Its kernel and coefficients are
!!            Q = 0                                for |k_d| <= m_N,
!!                exp(-((|k_d| - k_c)/(|k_d| - m_N))^2) for m_N < |k_d| < k_c,
!!                1 (svv_above_kc = T) or 0 (F)    for |k_d| >= k_c,
!!            eps_d = svv_eps0*dl_d*max(|u_d| + cf_d)   (active d only),
!!          with k_c,d = 2/3 of the Nyquist wavenumber (the edge of the
!!          dealiasing filter) and m_N = svv_mN*k_c,d. Here k_d is the
!!          wavenumber itself (the Nyquist mode is damped).
!!          B: dB^/dt = -lambda(k) B^, lambda(k) = sum_d eps_d k_d^2 Q(|k_d|).
!!          rho and m, svv_form = 'linear': the same -lambda(k) on m^ and,
!!            in the mass-conserving form, on rho^ (the energy rate is not
!!            sign-definite).
!!          rho and m, svv_form = 'entropy': flux form in the entropy
!!            variables w = (cs^2 s - u^2/2, u) (get_svv_es),
!!              A_L,d = F^-1 sig_d M F,  A_H,d = F^-1 sig_d (1 - M) F,
!!              sig_d = i k_d sqrt(Q) (|k_d| sqrt(Q) at Nyquist, so that A
!!                      maps a real field to a real field),
!!              g_r = A_L,d w_r, g_m = A_L,d u, g_h = A_H,d w_r,
!!              f_r = a t, t = g_r + u.g_m, f_m = u f_r + b g_m, f_h = kap rho g_h,
!!              a = eps_d rho/cs^2, b = eps_d rho, kap = eps_d/(cs^2 + max|u|^2),
!!              d rho/dt = -sum_d (A_L,d^* f_r + A_H,d^* f_h),
!!              dm/dt    = -sum_d  A_L,d^* f_m,
!!            where M is the filter of m (identity if filter_mom = F; then
!!            A_H = 0; M is 0 or 1, so the filter of m after each stage
!!            leaves dm/dt unchanged). Then sum d rho/dt = sum dm/dt = 0 (sig_d = 0 at
!!            k_d = 0) and dE/dt = -D_svv with
!!              D_svv = sum_d sum [a t^2 + b |g_m|^2 + kap rho g_h^2] dV >= 0
!!            for any state, if the s tendency is not filtered
!!            (filter_lnrho = F). The explicit step is limited by the
!!            largest |eigenvalue| of the Jacobian of this term, estimated by
!!            the power method (svv_power).
!-----------------------------------------------!
module shock_viscosity
  use p3dfft
  implicit none
  private
  include 'mpif.h'

  public init_shock_viscosity
  public get_div_u, get_zeta
  public kx_d, ky_d, kz_d, kd2_max
  public delta_shock, active_x, active_y, active_z
  public is_nonfinite, abort_run
  public svv_on, kq2x, kq2y, kq2z, svv_lam_m, svv_lam_s, svv_lam_b, get_svv_eps
  public svv_es, get_svv_es, svv_power, svv_lam
  public tmr_on, tmr_off, report_timers
  public itmr_step, itmr_zeta, itmr_svv_es, itmr_svv_pi, itmr_budget, itmr_planes

  logical  :: initialized = .false.
  real(r8), allocatable :: kx_d(:), ky_d(:), kz_d(:) ! wavenumbers for derivatives (0 at Nyquist)
  real(r8) :: kd2_max     ! sum_d (max |k_d| over the modes kept by M)^2; used for the viscous dt
  real(r8) :: delta_shock ! Delta in zeta = C_shock*Delta^2*S[...]
  logical  :: active_x, active_y, active_z
  logical  :: split_x, split_z
  integer  :: proc_xlo, proc_xhi, proc_zlo, proc_zhi

  ! spectral vanishing viscosity
  logical  :: svv_on = .false.
  real(r8), allocatable :: kq2x(:), kq2y(:), kq2z(:) ! k_d^2 Q(|k_d|); 0 in inactive directions
  real(r8) :: svv_lam_m(3), svv_lam_s(3), svv_lam_b(3) ! max_k k_d^2 Q over the modes kept in m, ln rho and B
  logical  :: svv_es = .false.                        ! svv_form = 'entropy'
  complex(r8), allocatable :: sgx(:), sgy(:), sgz(:)  ! sig_d = i k_d sqrt(Q) (|k_d| sqrt(Q) at Nyquist); 0 if inactive
  real(r8) :: svv_lam = 0.d0                          ! power-method estimate of the largest |eigenvalue| (svv_power)
  real(r8), allocatable :: pi_v(:,:,:,:)              ! power-method vector (s, mx, my, mz) in real space

  ! wall-clock time [s] of parts of the step, summed over the run on each rank (see report_timers)
  integer, parameter :: ntmr = 6
  integer, parameter :: itmr_step = 1, itmr_zeta = 2, itmr_svv_es = 3, itmr_svv_pi = 4, itmr_budget = 5, itmr_planes = 6
  character(len=28), parameter :: tmr_name(ntmr) = [character(len=28) :: &
    'step (solve)', 'div u, zeta, eta (halo)', 'get_svv_es', 'svv_power', 'stage_budget', 'write_planes']
  real(r8) :: tmr(ntmr) = 0.d0, tmr_t0(ntmr) = 0.d0

contains


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Initialization: derivative wavenumbers,
!!          neighbor ranks, Delta. Must be called
!!          after init_filter.
!-----------------------------------------------!
  subroutine init_shock_viscosity
    use grid, only: nlx, nly, nlz, dlx, dly, dlz
    use grid, only: kx, ky, kz, ikx, iky, ikz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: idx2procid
    use dealias, only: filter
    use mp, only: proc0, max_allreduce
    use grid, only: kx_max_noprune, ky_max_noprune, kz_max_noprune
    use params, only: C_shock, shock_dl, shock_nmax, shock_nsmooth, filter_mom, shock_selftest, C_eta
    use params, only: svv_eps0, svv_mN, svv_above_kc, filter_lnrho, filter_mag, svv_form, cfl_svv
    implicit none
    real(r8) :: kmx2, kmy2, kmz2, mfilt, sfilt, bfilt
    integer :: i, j, k

    if(initialized) return
    initialized = .true.

    ! wavenumbers for derivatives: Nyquist mode set to 0
    allocate(kx_d, source=kx)
    allocate(ky_d, source=ky)
    allocate(kz_d, source=kz)
    do i = 1, size(kx)
      if(abs(ikx(i)) == nlx/2) kx_d(i) = 0.d0
    enddo
    do j = 1, size(ky)
      if(abs(iky(j)) == nly/2) ky_d(j) = 0.d0
    enddo
    do k = 1, size(kz)
      if(abs(ikz(k)) == nlz/2) kz_d(k) = 0.d0
    enddo

    ! largest |k_d| among the modes kept in the momentum
    kmx2 = 0.d0; kmy2 = 0.d0; kmz2 = 0.d0
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          mfilt = merge(filter(i,k,j), 1.d0, filter_mom)
          if(mfilt > 0.d0) then
            kmx2 = max(kmx2, kx_d(i)**2)
            kmy2 = max(kmy2, ky_d(j)**2)
            kmz2 = max(kmz2, kz_d(k)**2)
          endif
        enddo
      enddo
    enddo
    call max_allreduce(kmx2)
    call max_allreduce(kmy2)
    call max_allreduce(kmz2)
    kd2_max = kmx2 + kmy2 + kmz2

    ! active directions
    active_x = nlx > 2
    active_y = nly > 2
    active_z = nlz > 2

    ! Delta
    if(shock_dl > 0.d0) then
      delta_shock = shock_dl
    else
      delta_shock = huge(1.d0)
      if(active_x) delta_shock = min(delta_shock, dlx)
      if(active_y) delta_shock = min(delta_shock, dly)
      if(active_z) delta_shock = min(delta_shock, dlz)
      if(.not. (active_x .or. active_y .or. active_z)) delta_shock = 0.d0
    endif

    ! spectral vanishing viscosity: k_d^2 Q per direction and its maxima over the kept modes
    svv_on = svv_eps0 > 0.d0
    allocate(kq2x, mold=kx)
    allocate(kq2y, mold=ky)
    allocate(kq2z, mold=kz)
    call svv_kernel(kx, kx_max_noprune, active_x, kq2x)
    call svv_kernel(ky, ky_max_noprune, active_y, kq2y)
    call svv_kernel(kz, kz_max_noprune, active_z, kq2z)
    svv_es = svv_on .and. trim(svv_form) == 'entropy'
    allocate(sgx(size(kx)), sgy(size(ky)), sgz(size(kz)))
    call svv_symbol(kx, ikx, nlx, kq2x, sgx)
    call svv_symbol(ky, iky, nly, kq2y, sgy)
    call svv_symbol(kz, ikz, nlz, kq2z, sgz)
    svv_lam_m = 0.d0
    svv_lam_s = 0.d0
    svv_lam_b = 0.d0
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          mfilt = merge(filter(i,k,j), 1.d0, filter_mom)
          sfilt = merge(filter(i,k,j), 1.d0, filter_lnrho)
          bfilt = merge(filter(i,k,j), 1.d0, filter_mag)
          if(mfilt > 0.d0) svv_lam_m = max(svv_lam_m, (/kq2x(i), kq2y(j), kq2z(k)/))
          if(sfilt > 0.d0) svv_lam_s = max(svv_lam_s, (/kq2x(i), kq2y(j), kq2z(k)/))
          if(bfilt > 0.d0) svv_lam_b = max(svv_lam_b, (/kq2x(i), kq2y(j), kq2z(k)/))
        enddo
      enddo
    enddo
    call max_allreduce(svv_lam_m)
    call max_allreduce(svv_lam_s)
    call max_allreduce(svv_lam_b)

    ! y must be local in real space
    if(ily_st /= 1 .or. ily_en /= nly) call abort_run('shock_viscosity: y is not local in real space')

    ! neighbor ranks (real space is a tensor product of z and x slabs)
    split_z = (ilz_en - ilz_st + 1) /= nlz
    split_x = (ilx_en - ilx_st + 1) /= nlx
    proc_zlo = idx2procid(mod(ilz_st - 2 + nlz, nlz) + 1, ilx_st, 'r')
    proc_zhi = idx2procid(mod(ilz_en         , nlz) + 1, ilx_st, 'r')
    proc_xlo = idx2procid(ilz_st, mod(ilx_st - 2 + nlx, nlx) + 1, 'r')
    proc_xhi = idx2procid(ilz_st, mod(ilx_en         , nlx) + 1, 'r')

    if(proc0) then
      print *
      print '("  shock viscosity: C_shock = ", es12.4, ", Delta = ", es12.4, &
             &", n_max = ", i3, ", n_smooth = ", i3)', C_shock, delta_shock, shock_nmax, shock_nsmooth
      print '("                   active (x, y, z) = ", 3l2, ", sum_d k_d,max^2 = ", es12.4)', &
             active_x, active_y, active_z, kd2_max
      if(C_eta > 0.d0) print '("  shock resistivity: C_eta = ", es12.4, " (eta = C_eta*Delta^2*S[max(0, -D)])")', C_eta
      if(svv_on) then
        print '("  spectral vanishing viscosity: eps0 = ", es12.4, ", m_N/k_c = ", es12.4, &
               &", Q(|k| >= k_c) = ", i2)', svv_eps0, svv_mN, merge(1, 0, svv_above_kc)
        print '("                   max k_d^2 Q (x, y, z) in m: ", 3es12.4, ", in ln rho: ", 3es12.4)', &
               svv_lam_m, svv_lam_s
        if(svv_es) then
          print '("                   rho and m: entropy-stable flux form, cfl_svv = ", es12.4)', cfl_svv
        else
          print '("                   rho and m: linear form, -lambda(k)")'
        endif
      else
        print '("  spectral vanishing viscosity: off")'
      endif
      print *
    endif

    if(shock_selftest) call selftest

  end subroutine init_shock_viscosity


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   kq2 = k^2 Q(|k|) of the SVV kernel for
!!          one direction (0 if the direction is
!!          inactive); kmax_noprune = Nyquist k
!-----------------------------------------------!
  subroutine svv_kernel(k, kmax_noprune, active, kq2)
    use params, only: svv_mN, svv_above_kc
    implicit none
    real(r8), intent(in ) :: k(:), kmax_noprune
    logical , intent(in ) :: active
    real(r8), intent(out) :: kq2(:)
    real(r8) :: kc, mN, a, r, Q
    integer  :: i

    kc = kmax_noprune*2.d0/3.d0 ! same threshold as init_filter
    mN = svv_mN*kc
    do i = 1, size(k)
      a = abs(k(i))
      if(a >= kc) then
        Q = merge(1.d0, 0.d0, svv_above_kc)
      elseif(a > mN) then
        r = (kc - a)/(a - mN)
        Q = merge(exp(-r**2), 0.d0, r < 30.d0) ! no underflow next to m_N
      else
        Q = 0.d0
      endif
      kq2(i) = merge(k(i)**2*Q, 0.d0, active)
    enddo

  end subroutine svv_kernel


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   sig = i k sqrt(Q) of the entropy-stable
!!          SVV for one direction, from kq2 = k^2 Q
!!          (real |k| sqrt(Q) at Nyquist; 0 if the
!!          direction is inactive, where kq2 = 0)
!-----------------------------------------------!
  subroutine svv_symbol(k, ik, nl, kq2, sg)
    use params, only: zi
    implicit none
    real(r8)   , intent(in ) :: k(:), kq2(:)
    integer    , intent(in ) :: ik(:), nl
    complex(r8), intent(out) :: sg(:)
    integer :: i

    do i = 1, size(k)
      if(abs(ik(i)) == nl/2) then
        sg(i) = sqrt(kq2(i))
      else
        sg(i) = zi*sign(sqrt(kq2(i)), k(i))
      endif
    enddo

  end subroutine svv_symbol


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   SVV coefficients eps_d = svv_eps0*dl_d
!!          *vmax_d (active d), vmax_d = max(|u_d| + cf_d)
!-----------------------------------------------!
  subroutine get_svv_eps(vmax, eps)
    use grid, only: dlx, dly, dlz
    use params, only: svv_eps0
    implicit none
    real(r8), intent(in ) :: vmax(3)
    real(r8), intent(out) :: eps(3)

    eps = 0.d0
    if(.not. svv_on) return
    if(active_x) eps(1) = svv_eps0*dlx*vmax(1)
    if(active_y) eps(2) = svv_eps0*dly*vmax(2)
    if(active_z) eps(3) = svv_eps0*dlz*vmax(3)

  end subroutine get_svv_eps


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   D = div_N u = F^-1 i k_d M F u
!!          kxd(ikx, iky) is the x-wavenumber
!!          (sheared when shear = T).
!!          uh (optional): F u/(nlx nly nlz)
!-----------------------------------------------!
  subroutine get_div_u(ux_r, uy_r, uz_r, kxd, div_r, uh)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use dealias, only: filter
    use params, only: zi, filter_mom
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_fft
    implicit none
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(in ) :: ux_r, uy_r, uz_r
    real(r8), dimension(ikx_st:, iky_st:)         , intent(in ) :: kxd
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(out) :: div_r
    complex(r8), dimension(ikx_st:, ikz_st:, iky_st:, :), intent(out), optional :: uh
    real(r8)   , allocatable, dimension(:,:,:) :: f_r
    complex(r8), allocatable, dimension(:,:,:) :: ux, uy, uz, divu
    real(r8) :: mfilt
    integer :: i, j, k

    allocate(f_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))
    allocate(ux  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(uy  , mold=ux)
    allocate(uz  , mold=ux)
    allocate(divu, mold=ux)

    if (proc0) call put_time_stamp(timer_fft)
    f_r = ux_r; call p3dfft_ftran_r2c(f_r, ux, 'fft')
    f_r = uy_r; call p3dfft_ftran_r2c(f_r, uy, 'fft')
    f_r = uz_r; call p3dfft_ftran_r2c(f_r, uz, 'fft')
    if (proc0) call put_time_stamp(timer_fft)

    if(present(uh)) then
      uh(:,:,:,1) = ux/nlx/nly/nlz
      uh(:,:,:,2) = uy/nlx/nly/nlz
      uh(:,:,:,3) = uz/nlx/nly/nlz
    endif

    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          mfilt = merge(filter(i,k,j), 1.d0, filter_mom)
          divu(i,k,j) = zi*(kxd(i,j)*ux(i,k,j) + ky_d(j)*uy(i,k,j) + kz_d(k)*uz(i,k,j)) &
                        *mfilt/nlx/nly/nlz
        enddo
      enddo
    enddo

    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_btran_c2r(divu, div_r, 'tff')
    if (proc0) call put_time_stamp(timer_fft)

    deallocate(f_r, ux, uy, uz, divu)

  end subroutine get_div_u


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Entropy-stable SVV of rho and m for the
!!          coefficients eps (module header):
!!          rdot_r = d rho/dt (real space),
!!          mdh    = dm/dt (Fourier, normalized as
!!                   the state; filtered by M),
!!          dsum   = local sum of a t^2 + b |g_m|^2
!!                   + kap rho g_h^2 (>= 0) over the
!!                   grid points and the directions.
!!          uh (optional): F u/(nlx nly nlz); it is
!!          computed here if absent.
!-----------------------------------------------!
  subroutine get_svv_es(s_r, rho_r, ux_r, uy_r, uz_r, eps, rdot_r, mdh, dsum, uh)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_fft
    implicit none
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(in ) :: s_r, rho_r, ux_r, uy_r, uz_r
    real(r8), intent(in ) :: eps(3)
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(out) :: rdot_r
    complex(r8), dimension(ikx_st:, ikz_st:, iky_st:, :), intent(out) :: mdh
    real(r8), intent(out) :: dsum
    complex(r8), dimension(ikx_st:, ikz_st:, iky_st:, :), intent(in), optional :: uh
    real(r8)   , allocatable, dimension(:,:,:)   :: f_r
    complex(r8), allocatable, dimension(:,:,:,:) :: uh_loc

    if(present(uh)) then
      call svv_es_core(s_r, rho_r, ux_r, uy_r, uz_r, uh, eps, rdot_r, mdh, dsum)
    else
      allocate(f_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))
      allocate(uh_loc(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, 3))
      if (proc0) call put_time_stamp(timer_fft)
      f_r = ux_r; call p3dfft_ftran_r2c(f_r, uh_loc(:,:,:,1), 'fft')
      f_r = uy_r; call p3dfft_ftran_r2c(f_r, uh_loc(:,:,:,2), 'fft')
      f_r = uz_r; call p3dfft_ftran_r2c(f_r, uh_loc(:,:,:,3), 'fft')
      if (proc0) call put_time_stamp(timer_fft)
      uh_loc = uh_loc/nlx/nly/nlz
      call svv_es_core(s_r, rho_r, ux_r, uy_r, uz_r, uh_loc, eps, rdot_r, mdh, dsum)
      deallocate(f_r, uh_loc)
    endif

  end subroutine get_svv_es


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Body of get_svv_es. Per active
!!          direction: 4 (5 with A_H) inverse and
!!          4 (5) forward FFTs; in addition 1
!!          forward (w_r) and 1 inverse (d rho/dt).
!-----------------------------------------------!
  subroutine svv_es_core(s_r, rho_r, ux_r, uy_r, uz_r, uh, eps, rdot_r, mdh, dsum)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use dealias, only: filter
    use params, only: cs2va2, filter_mom
    use mp, only: proc0, max_allreduce
    use time_stamp, only: put_time_stamp, timer_fft
    implicit none
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(in ) :: s_r, rho_r, ux_r, uy_r, uz_r
    complex(r8), dimension(ikx_st:, ikz_st:, iky_st:, :), intent(in ) :: uh
    real(r8), intent(in ) :: eps(3)
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(out) :: rdot_r
    complex(r8), dimension(ikx_st:, ikz_st:, iky_st:, :), intent(out) :: mdh
    real(r8), intent(out) :: dsum
    real(r8)   , allocatable, dimension(:,:,:)   :: gr, gh
    real(r8)   , allocatable, dimension(:,:,:,:) :: gm
    complex(r8), allocatable, dimension(:,:,:)   :: wh, fh, rdh
    real(r8) :: umax2, a, b, kap, t
    integer  :: i, j, k, c, d

    allocate(gr(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))
    allocate(gh, mold=gr)
    allocate(gm(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, 3))
    allocate(wh(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))
    allocate(fh , mold=wh)
    allocate(rdh, mold=wh)

    umax2 = maxval(ux_r**2 + uy_r**2 + uz_r**2)
    call max_allreduce(umax2)

    ! w_r = cs^2 s - u^2/2
    gr = cs2va2*s_r - 0.5d0*(ux_r**2 + uy_r**2 + uz_r**2)
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_ftran_r2c(gr, wh, 'fft')
    if (proc0) call put_time_stamp(timer_fft)
    wh = wh/nlx/nly/nlz

    rdh  = 0.d0
    mdh  = 0.d0
    dsum = 0.d0
    do d = 1, 3
      if(eps(d) == 0.d0) cycle ! inactive direction or svv off

      ! gradients g_r = A_L w_r, g_m = A_L u, g_h = A_H w_r
      call grad(wh, .true., gr)
      do c = 1, 3
        call grad(uh(:,:,:,c), .true., gm(:,:,:,c))
      enddo
      if(filter_mom) call grad(wh, .false., gh)

      ! fluxes (overwrite the gradients) and the dissipation density
      kap = eps(d)/(cs2va2 + umax2)
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            a = eps(d)*rho_r(j,k,i)/cs2va2
            b = eps(d)*rho_r(j,k,i)
            t = gr(j,k,i) + ux_r(j,k,i)*gm(j,k,i,1) + uy_r(j,k,i)*gm(j,k,i,2) + uz_r(j,k,i)*gm(j,k,i,3)
            dsum = dsum + a*t**2 + b*(gm(j,k,i,1)**2 + gm(j,k,i,2)**2 + gm(j,k,i,3)**2)
            gr(j,k,i) = a*t
            gm(j,k,i,1) = ux_r(j,k,i)*gr(j,k,i) + b*gm(j,k,i,1)
            gm(j,k,i,2) = uy_r(j,k,i)*gr(j,k,i) + b*gm(j,k,i,2)
            gm(j,k,i,3) = uz_r(j,k,i)*gr(j,k,i) + b*gm(j,k,i,3)
          enddo
        enddo
      enddo
      if(filter_mom) then
        dsum = dsum + sum(kap*rho_r*gh**2)
        gh = kap*rho_r*gh
      endif

      ! adjoints: d rho/dt -= A_L^* f_r + A_H^* f_h, dm/dt -= A_L^* f_m
      call grad_adj(gr, .true.)
      rdh = rdh - fh
      if(filter_mom) then
        call grad_adj(gh, .false.)
        rdh = rdh - fh
      endif
      do c = 1, 3
        call grad_adj(gm(:,:,:,c), .true.)
        mdh(:,:,:,c) = mdh(:,:,:,c) - fh
      enddo
    enddo

    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_btran_c2r(rdh, rdot_r, 'tff')
    if (proc0) call put_time_stamp(timer_fft)

    deallocate(gr, gh, gm, wh, fh, rdh)

  contains

    ! y = F^-1 sig_d mf x, mf = M (low = T) or 1 - M (low = F); fh is the work array
    subroutine grad(x, low, y)
      complex(r8), dimension(ikx_st:, ikz_st:, iky_st:), intent(in ) :: x
      logical, intent(in) :: low
      real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(out) :: y
      fh = x
      call mult_sig(low, .false.)
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_btran_c2r(fh, y, 'tff')
      if (proc0) call put_time_stamp(timer_fft)
    end subroutine grad

    ! fh = conj(sig_d) mf F y/(nlx nly nlz) (y is destroyed)
    subroutine grad_adj(y, low)
      real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(inout) :: y
      logical, intent(in) :: low
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_ftran_r2c(y, fh, 'fft')
      if (proc0) call put_time_stamp(timer_fft)
      fh = fh/nlx/nly/nlz
      call mult_sig(low, .true.)
    end subroutine grad_adj

    ! fh <- sig_d mf fh (conj(sig_d) for adj = T)
    subroutine mult_sig(low, adj)
      logical, intent(in) :: low, adj
      complex(r8) :: sg
      real(r8) :: mf
      integer  :: i, j, k

      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            select case(d)
              case(1); sg = sgx(i)
              case(2); sg = sgy(j)
              case(3); sg = sgz(k)
            end select
            if(adj) sg = conjg(sg)
            mf = merge(filter(i,k,j), 1.d0, filter_mom)
            if(.not. low) mf = 1.d0 - mf
            fh(i,k,j) = sg*mf*fh(i,k,j)
          enddo
        enddo
      enddo
    end subroutine mult_sig

  end subroutine svv_es_core


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   SVV tendency of the state (s, m) with
!!          eps fixed: Gs = (d rho/dt)/rho,
!!          Gm = dm/dt (real space)
!-----------------------------------------------!
  subroutine svv_G(s_r, mx_r, my_r, mz_r, eps, Gs, Gm)
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_fft
    implicit none
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(in ) :: s_r, mx_r, my_r, mz_r
    real(r8), intent(in ) :: eps(3)
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:)   , intent(out) :: Gs
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:, :), intent(out) :: Gm
    real(r8)   , allocatable, dimension(:,:,:)   :: rho_r, ux_r, uy_r, uz_r
    complex(r8), allocatable, dimension(:,:,:,:) :: mdh
    real(r8) :: dsum
    integer  :: c

    allocate(rho_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))
    allocate(ux_r, mold=rho_r)
    allocate(uy_r, mold=rho_r)
    allocate(uz_r, mold=rho_r)
    allocate(mdh(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, 3))

    rho_r = exp(s_r)
    ux_r  = mx_r/rho_r
    uy_r  = my_r/rho_r
    uz_r  = mz_r/rho_r
    call get_svv_es(s_r, rho_r, ux_r, uy_r, uz_r, eps, Gs, mdh, dsum)
    Gs = Gs/rho_r
    if (proc0) call put_time_stamp(timer_fft)
    do c = 1, 3
      call p3dfft_btran_c2r(mdh(:,:,:,c), Gm(:,:,:,c), 'tff')
    enddo
    if (proc0) call put_time_stamp(timer_fft)

    deallocate(rho_r, ux_r, uy_r, uz_r, mdh)

  end subroutine svv_G


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Power method for svv_lam, the largest
!!          |eigenvalue| of the Jacobian J of
!!          G(s, m) = ((d rho/dt)/rho, dm/dt) of the
!!          entropy-stable SVV with eps fixed:
!!            J v = (G(U + h v) - G(U))/h, h = 1e-7,
!!            v <- J v/|v|, svv_lam = |J v|/|v|,
!!          |v|^2 = sum (v_s^2 + |v_m|^2)/(nlx nly nlz).
!!          The vector is kept between calls (not
!!          across a restart). Its initial value is a
!!          fixed function of the global grid index,
!!          so that it does not depend on the
!!          decomposition. niter iterations at most
!!          (max(niter, svv_pi_init) when the vector is
!!          initialized); with tol > 0 it stops when
!!          two successive estimates of this call agree
!!          within tol (relative). nit: iterations
!!          done; conv: the tol test was met.
!!          G(U) is given by rdot0_r (d rho/dt) and
!!          mdh0 (dm/dt, Fourier) of get_svv_es.
!-----------------------------------------------!
  subroutine svv_power(s_r, rho_r, mx_r, my_r, mz_r, eps, rdot0_r, mdh0, niter, tol, nit, conv)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: pi, svv_pi_init
    use mp, only: proc0, sum_allreduce
    use time_stamp, only: put_time_stamp, timer_fft
    implicit none
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(in) :: s_r, rho_r, mx_r, my_r, mz_r, rdot0_r
    complex(r8), dimension(ikx_st:, ikz_st:, iky_st:, :), intent(in) :: mdh0
    real(r8), intent(in ) :: eps(3), tol
    integer , intent(in ) :: niter
    integer , intent(out) :: nit
    logical , intent(out) :: conv
    real(r8), parameter :: h = 1.d-7
    real(r8)   , allocatable, dimension(:,:,:)   :: G0s, Gs, ps
    real(r8)   , allocatable, dimension(:,:,:,:) :: G0m, Gm, pm
    complex(r8), allocatable, dimension(:,:,:)   :: wk
    real(r8) :: amp, nv, lam, lam_old
    integer  :: i, j, k, c, it, nmax, nsg

    nit  = 0
    conv = .false.
    if(all(eps == 0.d0)) then
      svv_lam = 0.d0
      return
    endif

    allocate(G0s(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))
    allocate(Gs, mold=G0s)
    allocate(ps, mold=G0s)
    allocate(G0m(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, 3))
    allocate(Gm , mold=G0m)
    allocate(pm , mold=G0m)
    allocate(wk(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))

    nmax = niter
    if(.not. allocated(pi_v)) then
      nmax = max(niter, svv_pi_init)
      ! (-1)^(sum_d (i_d - 1)) (1 + 0.3 sum_d cos(2 pi (i_d - 1)/nl_d + 1)) in s, 0 in m, sums over the
      ! active directions d: with nl_d = 2 the modes (-1)^(i_d - 1) of an inactive direction are not coupled
      ! by J to the fields uniform in it, and their stiffness is lower (M = 0 there, so only A_H acts)
      allocate(pi_v(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, 4))
      pi_v = 0.d0
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            amp = 1.d0
            nsg = 0
            if(active_x) then
              amp = amp + 0.3d0*cos(2.d0*pi*dble(i - 1)/dble(nlx) + 1.d0)
              nsg = nsg + i - 1
            endif
            if(active_y) then
              amp = amp + 0.3d0*cos(2.d0*pi*dble(j - 1)/dble(nly) + 1.d0)
              nsg = nsg + j - 1
            endif
            if(active_z) then
              amp = amp + 0.3d0*cos(2.d0*pi*dble(k - 1)/dble(nlz) + 1.d0)
              nsg = nsg + k - 1
            endif
            pi_v(j,k,i,1) = merge(1.d0, -1.d0, mod(nsg, 2) == 0)*amp
          enddo
        enddo
      enddo
    endif

    ! G(U)
    G0s = rdot0_r/rho_r
    do c = 1, 3
      wk = mdh0(:,:,:,c)
      if (proc0) call put_time_stamp(timer_fft)
      call p3dfft_btran_c2r(wk, G0m(:,:,:,c), 'tff')
      if (proc0) call put_time_stamp(timer_fft)
    enddo

    lam_old = 0.d0
    lam     = 0.d0
    do it = 1, nmax
      nv = norm()
      if(nv == 0.d0) then
        ! J v = 0 (degenerate); restart from the initial vector at the next call
        deallocate(pi_v)
        lam = 0.d0
        exit
      endif
      pi_v = pi_v/nv
      ps = s_r + h*pi_v(:,:,:,1)
      pm(:,:,:,1) = mx_r + h*pi_v(:,:,:,2)
      pm(:,:,:,2) = my_r + h*pi_v(:,:,:,3)
      pm(:,:,:,3) = mz_r + h*pi_v(:,:,:,4)
      call svv_G(ps, pm(:,:,:,1), pm(:,:,:,2), pm(:,:,:,3), eps, Gs, Gm)
      pi_v(:,:,:,1) = (Gs - G0s)/h
      do c = 1, 3
        pi_v(:,:,:,c+1) = (Gm(:,:,:,c) - G0m(:,:,:,c))/h
      enddo
      lam  = norm()
      nit  = it
      conv = tol > 0.d0 .and. it > 1 .and. abs(lam - lam_old) <= tol*lam
      lam_old = lam
      if(conv) exit
    enddo
    svv_lam = lam

    deallocate(G0s, Gs, ps, G0m, Gm, pm, wk)

  contains

    real(r8) function norm()
      real(r8) :: ss
      ss = sum(pi_v**2)
      call sum_allreduce(ss)
      norm = sqrt(ss/(dble(nlx)*dble(nly)*dble(nlz)))
    end function norm

  end subroutine svv_power


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   zeta = C_shock*Delta^2*S[max(0, -D)]
!!          and (optional) the shock resistivity
!!          eta = C_eta*Delta^2*S[max(0, -D)]
!-----------------------------------------------!
  subroutine get_zeta(div_r, zeta_r, eta_r)
    use grid, only: ily_st, ilz_st, ilx_st
    use params, only: C_shock, shock_nmax, shock_nsmooth, C_eta
    implicit none
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(in ) :: div_r
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(out) :: zeta_r
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(out), optional :: eta_r
    integer :: n

    zeta_r = max(0.d0, -div_r)

    do n = 1, shock_nmax
      if(active_y) call pass_1d(zeta_r, 'y', 'm')
      if(active_z) call pass_1d(zeta_r, 'z', 'm')
      if(active_x) call pass_1d(zeta_r, 'x', 'm')
    enddo
    do n = 1, shock_nsmooth
      if(active_y) call pass_1d(zeta_r, 'y', 's')
      if(active_z) call pass_1d(zeta_r, 'z', 's')
      if(active_x) call pass_1d(zeta_r, 'x', 's')
    enddo

    if(present(eta_r)) eta_r = C_eta*delta_shock**2*zeta_r
    zeta_r = C_shock*delta_shock**2*zeta_r

  end subroutine get_zeta


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   One 3-point pass in direction dir
!!          op = 'm': f_j <- max(f_{j-1}, f_j, f_{j+1})
!!          op = 's': f_j <- (f_{j-1} + 2 f_j + f_{j+1})/4,
!!                    evaluated as 0.5*(f_j + 0.5*(f_{j-1} + f_{j+1}))
!!                    so that a uniform field is kept exactly
!-----------------------------------------------!
  subroutine pass_1d(f, dir, op)
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    implicit none
    real(r8), dimension(ily_st:, ilz_st:, ilx_st:), intent(inout) :: f
    character(len=1), intent(in) :: dir, op
    real(r8), allocatable, dimension(:,:,:) :: g
    integer :: i, j, k, dj, dk, di

    dj = 0; dk = 0; di = 0
    select case(dir)
    case('y')
      dj = 1
    case('z')
      dk = 1
    case('x')
      di = 1
    end select

    ! g = f extended by one halo layer in direction dir
    allocate(g(ily_st - dj:ily_en + dj, ilz_st - dk:ilz_en + dk, ilx_st - di:ilx_en + di))
    g(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en) = f

    select case(dir)
    case('y') ! always local
      g(ily_st - 1, :, :) = f(ily_en, :, :)
      g(ily_en + 1, :, :) = f(ily_st, :, :)
    case('z')
      if(split_z) then
        call exchange_halo(f(:, ilz_st, :), f(:, ilz_en, :), g(:, ilz_st - 1, :), g(:, ilz_en + 1, :), &
                           proc_zlo, proc_zhi)
      else
        g(:, ilz_st - 1, :) = f(:, ilz_en, :)
        g(:, ilz_en + 1, :) = f(:, ilz_st, :)
      endif
    case('x')
      if(split_x) then
        call exchange_halo(f(:, :, ilx_st), f(:, :, ilx_en), g(:, :, ilx_st - 1), g(:, :, ilx_en + 1), &
                           proc_xlo, proc_xhi)
      else
        g(:, :, ilx_st - 1) = f(:, :, ilx_en)
        g(:, :, ilx_en + 1) = f(:, :, ilx_st)
      endif
    end select

    if(op == 'm') then
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            f(j,k,i) = max(g(j-dj,k-dk,i-di), g(j,k,i), g(j+dj,k+dk,i+di))
          enddo
        enddo
      enddo
    else
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            f(j,k,i) = 0.5d0*(g(j,k,i) + 0.5d0*(g(j-dj,k-dk,i-di) + g(j+dj,k+dk,i+di)))
          enddo
        enddo
      enddo
    endif

    deallocate(g)

  end subroutine pass_1d


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Exchange one boundary plane with the
!!          neighbor ranks (periodic)
!!          face_lo/hi: own lowest/highest plane
!!          halo_lo/hi: plane below/above own range
!-----------------------------------------------!
  subroutine exchange_halo(face_lo, face_hi, halo_lo, halo_hi, proc_lo, proc_hi)
    implicit none
    real(r8), dimension(:,:), intent(in ) :: face_lo, face_hi
    real(r8), dimension(:,:), intent(out) :: halo_lo, halo_hi
    integer, intent(in) :: proc_lo, proc_hi
    real(r8), allocatable, dimension(:,:) :: sbuf, rbuf
    integer :: n, ierr
    integer :: status(MPI_STATUS_SIZE)

    n = size(face_lo)
    allocate(sbuf, mold=face_lo)
    allocate(rbuf, mold=face_lo)

    ! own highest plane -> lower halo of proc_hi
    sbuf = face_hi
    call MPI_Sendrecv(sbuf, n, MPI_DOUBLE_PRECISION, proc_hi, 1, &
                      rbuf, n, MPI_DOUBLE_PRECISION, proc_lo, 1, MPI_COMM_WORLD, status, ierr)
    halo_lo = rbuf

    ! own lowest plane -> upper halo of proc_lo
    sbuf = face_lo
    call MPI_Sendrecv(sbuf, n, MPI_DOUBLE_PRECISION, proc_lo, 2, &
                      rbuf, n, MPI_DOUBLE_PRECISION, proc_hi, 2, MPI_COMM_WORLD, status, ierr)
    halo_hi = rbuf

    deallocate(sbuf, rbuf)

  end subroutine exchange_halo


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Self test (shock_selftest = T):
!!          zeta of an analytic D on the grid is
!!          gathered to proc0 and written to
!!          shock_selftest.dat (es25.17). Files from
!!          different decompositions must be identical.
!-----------------------------------------------!
  subroutine selftest
    use grid, only: nlx, nly, nlz, xx, yy, zz, lx, ly, lz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: pi
    use mp, only: proc0, sum_reduce
    use file, only: open_output_file
    implicit none
    real(r8), allocatable, dimension(:,:,:) :: d_r, z_r, z_all
    integer :: i, j, k, unit

    allocate(d_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))
    allocate(z_r, mold=d_r)
    allocate(z_all(nly, nlz, nlx), source=0.d0)

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          d_r(j,k,i) = sin(2.d0*pi*xx(i)/lx + 0.3d0)*cos(4.d0*pi*yy(j)/ly + 1.1d0) &
                     + 0.7d0*sin(6.d0*pi*zz(k)/lz + 0.5d0) - 0.2d0
        enddo
      enddo
    enddo

    call get_zeta(d_r, z_r)

    z_all(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en) = z_r
    call sum_reduce(z_all, 0) ! exact: other ranks contribute 0

    if(proc0) then
      call open_output_file(unit, 'shock_selftest.dat')
      do i = 1, nlx
        do k = 1, nlz
          do j = 1, nly
            write(unit, '(3i6, es25.17)') i, j, k, z_all(j,k,i)
          enddo
        enddo
      enddo
      close(unit)
      print '("  shock_selftest.dat written: min/max zeta = ", 2es12.4)', minval(z_all), maxval(z_all)
    endif

    deallocate(d_r, z_r, z_all)

  end subroutine selftest


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   True for Inf/NaN. Tests the exponent
!!          bits so that it is not removed by
!!          value-unsafe optimizations (-O3, fast fp)
!-----------------------------------------------!
  elemental logical function is_nonfinite(x)
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    real(r8), intent(in) :: x
    integer(int64) :: bits

    bits = transfer(x, 0_int64)
    is_nonfinite = iand(ishft(bits, -52), 2047_int64) == 2047_int64

  end function is_nonfinite


!-----------------------------------------------!
!> @author  YK
!! @date    31 Oct 2026
!! @brief   Print a message with the rank and
!!          abort all ranks (non-zero exit code).
!!          Does not require the other ranks to
!!          reach the same point.
!-----------------------------------------------!
  subroutine abort_run(msg)
    use mp, only: proc_id
    implicit none
    character(len=*), intent(in) :: msg
    integer :: ierr

    write(*, '("[proc ", i6, "] ABORT: ", a)') proc_id, trim(msg)
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  end subroutine abort_run


!-----------------------------------------------!
!> @author  YK
!! @date    3 Nov 2026
!! @brief   Start and stop the timer i (MPI_Wtime
!!          on each rank; no synchronization)
!-----------------------------------------------!
  subroutine tmr_on(i)
    implicit none
    integer, intent(in) :: i

    tmr_t0(i) = MPI_Wtime()

  end subroutine tmr_on

  subroutine tmr_off(i)
    implicit none
    integer, intent(in) :: i

    tmr(i) = tmr(i) + (MPI_Wtime() - tmr_t0(i))

  end subroutine tmr_off


!-----------------------------------------------!
!> @author  YK
!! @date    3 Nov 2026
!! @brief   Print the min and max over the ranks
!!          of the timers (collective)
!-----------------------------------------------!
  subroutine report_timers
    use mp, only: proc0, max_reduce, min_reduce
    implicit none
    real(r8) :: tmin(ntmr), tmax(ntmr)
    integer :: i

    tmin = tmr
    tmax = tmr
    call min_reduce(tmin, 0)
    call max_reduce(tmax, 0)
    if(proc0) then
      print '(/, " MHD_COMP_ISOTH timers [s] (min, max over ranks; % of max step)")'
      do i = 1, ntmr
        print '(3x, a28, 2es12.4, f8.2, " %")', tmr_name(i), tmin(i), tmax(i), &
              100.d0*tmax(i)/max(tmax(itmr_step), tiny(1.d0))
      enddo
    endif

  end subroutine report_timers

end module shock_viscosity
