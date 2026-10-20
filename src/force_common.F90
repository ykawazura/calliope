!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Forcing; based on AstroGK
!-----------------------------------------------!
module force_common
  use p3dfft
  implicit none

  public update_force
  public driven
  public nk_stir, a_force, b_force
  public f_name, kx_stir, ky_stir, kz_stir
  public is_forced, zeroing_at_forcing_wavenumber

  logical :: driven
  integer :: nk_stir
  logical :: elsasser, fix_power
  real(r8):: ene_inj, xhl_inj, res_inj
  real(r8):: kmin(3), kmax(3)
  integer :: nfields
  !> Base seed for the stochastic forcing RNG. stir_seed < 0 (default) draws a
  !> wall-clock seed => nondeterministic run-to-run. stir_seed >= 0 gives a
  !> fixed, reproducible forcing realization, independent of the number of MPI ranks.
  integer :: stir_seed
  character(len=100), allocatable :: field_names(:)
  real   (r8), allocatable :: amplitudes(:)
  complex(r8), allocatable :: frequencies(:)

  complex(r8), dimension(:,:), allocatable :: a_force, b_force
  integer :: nseed                   !Length of random number  seed integer vector
  integer, dimension(:), allocatable :: init_seed, fin_seed   !Initial and final seeds

  integer, dimension(:,:), allocatable :: kx_stir, ky_stir, kz_stir
  character(len=100), dimension(:,:), allocatable :: f_name
  real   (r8)       , dimension(:,:), allocatable :: f_amplitude
  complex(r8)       , dimension(:,:), allocatable :: f_frequency

  integer, parameter :: kind_id = selected_int_kind (15)

  integer :: force_unit
contains


!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Initialization of forcing
!-----------------------------------------------!
  subroutine init_force
    use params, only: inputfile
    implicit none

    !Get random number seed length and allocate init_seed and fin_seed
    nseed = get_rnd_seed_length() 
    allocate(init_seed(1:nseed), source=0)
    allocate(fin_seed (1:nseed), source=0)

    call read_parameters(inputfile)

  end subroutine init_force


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Read inputfile for forcing parameters
!-----------------------------------------------!
  subroutine read_parameters(filename)
    use MPI
    use params, only: pi
    use mp, only: proc0
    use mp, only: broadcast
    use grid, only: nkx, nky, nkz
    use grid, only: ikx, iky, ikz
    use file, only: get_unused_unit, get_indexed_namelist_unit, open_output_file
    implicit none
    integer :: i, j, k, count, ifield, unit, in_file
    
    character(len=100), intent(in) :: filename
    integer  :: ierr

    namelist /force/ driven, elsasser, fix_power, ene_inj, xhl_inj, res_inj, kmin, kmax, nfields, stir_seed
    namelist /forced_fields/ field_names, amplitudes, frequencies

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v    used only when the corresponding value   v!
    !v    does not exist in the input file         v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    elsasser = .false.
    driven = .false.
    fix_power = .false.
    ene_inj = 0.d0
    xhl_inj = 0.d0
    res_inj = 0.d0
    kmin = (/0.d0, 0.d0, 0.d0/)
    kmax = (/0.d0, 0.d0, 0.d0/)
    nfields = 0
    stir_seed = -1   ! < 0 => wall-clock (random) seed
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=force,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading force_parameters failed"
    close(unit)

    allocate(field_names(nfields))
    allocate(amplitudes (nfields))
    allocate(frequencies(nfields))

    if (proc0) then
      call get_unused_unit (unit)
      open(unit=unit,file=filename,status='old')

      read(unit,nml=forced_fields,iostat=ierr)
          if (ierr/=0) write(*,*) "Reading force_parameters failed"
      close(unit)
    endif

    call broadcast(field_names)
    call broadcast(amplitudes)
    call broadcast(frequencies)

    ! The list of forced modes is built from the global wavenumbers, so every
    ! rank holds the same list and draws the same random numbers. On the ky = 0
    ! plane, only half of the modes are listed; get_force sets the other half
    ! by the reality condition.
    count = 0
    do j = 1, nky
      do k = 1, nkz
        do i = 1, nkx
          if(is_stirred(i, j, k)) count = count + 1
        end do
      end do
    end do

    nk_stir = count

    allocate (kx_stir    (nfields, nk_stir))
    allocate (ky_stir    (nfields, nk_stir))
    allocate (kz_stir    (nfields, nk_stir))
    allocate (f_name     (nfields, nk_stir))
    allocate (f_amplitude(nfields, nk_stir))
    allocate (f_frequency(nfields, nk_stir))
    allocate (a_force    (nfields, nk_stir))
    allocate (b_force    (nfields, nk_stir))

    do ifield = 1, nfields
      count = 1
      do j = 1, nky
        do k = 1, nkz
          do i = 1, nkx
            if(is_stirred(i, j, k)) then

              kx_stir    (ifield, count) = ikx(i)
              ky_stir    (ifield, count) = iky(j)
              kz_stir    (ifield, count) = ikz(k)
              f_name     (ifield, count) = field_names(ifield)
              f_amplitude(ifield, count) = amplitudes (ifield)
              f_frequency(ifield, count) = frequencies(ifield)
              a_force    (ifield, count) = amplitudes (ifield)*cmplx(1.d0,1.d0)/2.d0 
              b_force    (ifield, count) = amplitudes (ifield)*cmplx(1.d0,1.d0)/2.d0 

              count = count + 1
            endif
          end do
        end do
      end do
    end do

    ! Seed the stochastic forcing. Default (stir_seed < 0) draws a wall-clock
    ! seed => nondeterministic between runs; stir_seed >= 0 gives a fixed,
    ! reproducible forcing realization. proc0 builds the seed and broadcasts it
    ! so every rank shares the same forcing sequence.
    if(proc0) then
      if (stir_seed < 0) then
        call init_ranf(.true., init_seed)
      else
        call seed_from_int(stir_seed, init_seed)
      endif
    endif
    call broadcast (init_seed)
    fin_seed = init_seed

    if(proc0) print '("  forcing: nk_stir = ", i0, " modes per field, stir_seed = ", i0)', nk_stir, stir_seed

    ! if(proc0) then
    !   call open_output_file (force_unit, 'force.dat')
    ! endif

  end subroutine read_parameters


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2019
!! @brief   Update randomized forcing amplitude;
!           must be called every time step
!           when forcing
!-----------------------------------------------!
  subroutine update_force(dt)
    use params, only: pi, zi
    use time, only: tt
    use grid, only: lz
    use mp, only: proc0
    use mp, only: broadcast
    use file, only: open_output_file
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    real(r8), intent(in) :: dt
    complex(r8), dimension(:,:), allocatable :: w_stir
    complex(r8) :: fa, fb
    real(r8) :: sigma
    integer :: ifield, istir

    if (proc0) call put_time_stamp(timer_force)

    allocate (w_stir(nfields, nk_stir))

    ! Continue the forcing's own random sequence, so that other calls of
    ! random_seed/random_number (e.g., random initial conditions) do not alter it
    call random_seed(put=fin_seed)

    do istir = 1, nk_stir
      do ifield = 1, nfields

        w_stir(ifield, istir) = 2.d0*pi/lz*abs(kz_stir(ifield, istir))*f_frequency(ifield, istir)
        
        sigma = sqrt(12.d0*abs(aimag(w_stir(ifield, istir)))/dt)*f_amplitude(ifield, istir)
        fa = cmplx(ranf() - 0.5d0, ranf() - 0.5d0) * sigma
        fb = cmplx(ranf() - 0.5d0, ranf() - 0.5d0) * sigma

        a_force(ifield, istir) = a_force(ifield, istir)*exp(-zi*w_stir(ifield, istir)*dt) + fa*dt
        b_force(ifield, istir) = 0.d0
      end do
    end do

    call get_rnd_seed(fin_seed)

    ! if(proc0) then
    !   write (unit=force_unit, fmt="(100es30.21)", advance='no') tt + dt
    !   do istir = 1, nk_stir
    !     do ifield = 1, nfields
    !       write (unit=force_unit, fmt="(100es30.21)", advance='no') &
    !                                 real(a_force(ifield, istir)), imag(a_force(ifield, istir)), &
    !                                 real(b_force(ifield, istir)), imag(b_force(ifield, istir))
    !     end do
    !   end do
    !   write (unit=force_unit, fmt="(100es30.21)")
    !   flush (force_unit)
    ! endif

    deallocate (w_stir)

    if (proc0) call put_time_stamp(timer_force)
  end subroutine update_force


!-----------------------------------------------!
!> @author  YK
!! @date    5 Oct 2019
!! @brief   Return forcing field when 'name'
!!          matches with 'field_name' of the
!!          inputfile
!-----------------------------------------------!
  subroutine get_force (name, u)
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    character(*) :: name
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: u
    integer :: i, ifield, istir
    complex(r8) :: f

    if (proc0) call put_time_stamp(timer_force)

    ifield = 0
    do i = 1, nfields
      if(trim(name) == field_names(i)) ifield = i
    enddo

    if (ifield == 0) return
    
    do istir = 1, nk_stir
      f = (a_force(ifield, istir) + b_force(ifield, istir))/sqrt(2.d0)
      call set_mode(u, kx_stir(ifield, istir), ky_stir(ifield, istir), kz_stir(ifield, istir), f)

      ! reality condition on the ky = 0 plane: f(-k) = conjg(f(k))
      if(ky_stir(ifield, istir) == 0) &
        call set_mode(u, -kx_stir(ifield, istir), 0, -kz_stir(ifield, istir), conjg(f))
    end do

    if (proc0) call put_time_stamp(timer_force)

  end subroutine get_force



!-----------------------------------------------!
!> @author  YK
!! @date    19 Sep 2025
!! @brief   Return whether the field `name` 
!!          is forced or not
!-----------------------------------------------!
  function is_forced (name)
    implicit none
    logical :: is_forced
    character(*) :: name
    integer :: i

    is_forced = .false.

    do i = 1, nfields
      if(trim(name) == field_names(i)) is_forced = .true.
    enddo

  end function is_forced



!-----------------------------------------------!
!> @author  YK
!! @date    19 Sep 2025
!! @brief   Enforce zero the field at the forcing
!!          wavenumber (this field should be  
!!          unforced one )
!-----------------------------------------------!
  subroutine zeroing_at_forcing_wavenumber (u)
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: u
    integer :: ifield, istir

    if (proc0) call put_time_stamp(timer_force)

    do ifield = 1, nfields
      do istir = 1, nk_stir
        call set_mode(u, kx_stir(ifield, istir), ky_stir(ifield, istir), kz_stir(ifield, istir), (0.0_r8, 0.0_r8))

        if(ky_stir(ifield, istir) == 0) &
          call set_mode(u, -kx_stir(ifield, istir), 0, -kz_stir(ifield, istir), (0.0_r8, 0.0_r8))
      end do
    enddo

    if (proc0) call put_time_stamp(timer_force)
  end subroutine zeroing_at_forcing_wavenumber


!-----------------------------------------------!
!> @brief   Whether the mode (kx(i), ky(j), kz(k))
!!          is listed as a forced mode
!-----------------------------------------------!
  function is_stirred (i, j, k)
    use grid, only: kx, ky, kz, ikx, iky, ikz
    use grid, only: nlxc, nlyc, nlzc
    implicit none
    logical :: is_stirred
    integer, intent(in) :: i, j, k

    ! In the band kmin <= |k| <= kmax, without the Nyquist mode and, for an even
    ! pruned size, the unpaired edge mode (on the ky = 0 plane, neither has a
    ! distinct -k partner on the grid; both are left out on every plane).
    ! On the ky = 0 plane, only the half kz > 0 .or. (kz == 0 .and. kx > 0) is
    ! listed (k = 0 excluded); get_force sets the other half.
    is_stirred = in_band(kx(i), 1) .and. in_band(ky(j), 2) .and. in_band(kz(k), 3) &
      .and. 2*abs(ikx(i)) < nlxc .and. 2*iky(j) < nlyc .and. 2*abs(ikz(k)) < nlzc &
      .and. (iky(j) > 0 .or. ikz(k) > 0 .or. (ikz(k) == 0 .and. ikx(i) > 0))

  contains

    ! kmin(d) <= |k| <= kmax(d), tolerating round-off in k
    ! (e.g., kz = 0.9999999999999999 with Fujitsu -Kfast)
    logical function in_band (kk, d)
      real(r8), intent(in) :: kk
      integer , intent(in) :: d
      real(r8), parameter :: tol = 1.d-8

      in_band = kk**2 >= kmin(d)**2*(1.d0 - tol) .and. kk**2 <= kmax(d)**2*(1.d0 + tol)
    end function in_band

  end function is_stirred


!-----------------------------------------------!
!> @brief   u(k) = val if the mode with integer
!!          wavenumbers (jkx, jky, jkz), jky >= 0,
!!          is on this rank
!-----------------------------------------------!
  subroutine set_mode (u, jkx, jky, jkz, val)
    use grid, only: nkx, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: u
    integer, intent(in) :: jkx, jky, jkz
    complex(r8), intent(in) :: val
    integer :: i, j, k

    i = modulo(jkx, nkx) + 1
    j = jky + 1
    k = modulo(jkz, nkz) + 1

    if(      (i >= ikx_st .and. i <= ikx_en) &
       .and. (j >= iky_st .and. j <= iky_en) &
       .and. (k >= ikz_st .and. k <= ikz_en) ) u(i, k, j) = val

  end subroutine set_mode


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v    Below is based on 'file_utils.fpp' of AstroGK    v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  function get_rnd_seed_length () result (l)
    !  get_rnd_seed_length gets the length of the integer vector for
    !      the random number generator seed
    integer :: l
    call random_seed(size=l)
    
  end function get_rnd_seed_length

  subroutine get_rnd_seed(seed)
    !  get_rnd_seed  random number seed integer vector
    integer, dimension(:), intent(out) :: seed
    call random_seed(get=seed)
  end subroutine get_rnd_seed

  subroutine init_ranf(randomize,init_seed)
    !  init_ranf seeds the choosen random number generator.
    !  if randomize=T, a random seed based on date and time is used.
    !  (this randomizing procedure is proposed in GNU gfortran manual.)
    !  Otherwise, it sets the seed using init_seed 
    !  In either case, it outputs the initial seed in init_seed
    implicit none
    logical, intent(in) :: randomize
    integer, intent(inout), dimension(:) :: init_seed
    integer :: i, nseed, dt(8)
    integer (kind=kind_id) :: t
    
    call system_clock(t)
    if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_kind_id * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_kind_id * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_kind_id * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
    end if
       
    if (randomize) then
      call random_seed(size=nseed)
       
      do i = 1, nseed
        init_seed(i) = lcg(t)
      end do
       
      call random_seed(put=init_seed)
    else
      call random_seed(put=init_seed)
    endif

  end subroutine init_ranf

  subroutine seed_from_int(base, seed_out)
    !  Fill the RNG seed vector deterministically from an integer base,
    !  reusing the lcg mixer of the randomize path in init_ranf. Independent of
    !  the wall clock, so the forcing realization is bitwise-reproducible.
    implicit none
    integer, intent(in) :: base
    integer, intent(inout), dimension(:) :: seed_out
    integer :: i, n
    integer (kind=kind_id) :: t

    call random_seed(size=n)
    t = int(base, kind_id) + 1_kind_id   ! +1 avoids the t==0 degenerate branch in lcg
    do i = 1, n
      seed_out(i) = lcg(t)
    end do

  end subroutine seed_from_int

  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer (kind=kind_id) :: s

    if (kind_id > 0) then
      if (s == 0) then
        s = 104729
      else
        s = mod(s, 4294967296_kind_id)
      end if
      s = mod(s * 279470273_kind_id, 4294967291_kind_id)
      lcg = int(mod(s, int(huge(0), kind_id)), kind(0))
    else
      lcg = mod(16807*s, 2147483647)
    end if
    
  end function lcg

  function ranf (seed)
    !  returns a uniform deviate in [0., 1.)
    !  The generator is initialized with the given seed if exists,
    !  otherwise uses the default seed.
    integer, intent(in), optional :: seed
    real :: ranf
    integer :: l
    integer, allocatable :: seed_in(:)
    if (present(seed)) then
       call random_seed(size=l)
       allocate(seed_in(l))
       seed_in(:)=seed
       call random_seed(put=seed_in)
    endif
    call random_number(ranf)

  end function ranf

end module force_common

