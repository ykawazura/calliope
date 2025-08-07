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

  logical :: driven
  integer :: nk_stir
  logical :: elsasser, fix_power
  real(r8):: ene_inj, xhl_inj, res_inj
  real(r8):: kmin(3), kmax(3)
  integer :: nfields
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
    use grid, only: kx, ky, kz
    use grid, only: ikx, iky, ikz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use file, only: get_unused_unit, get_indexed_namelist_unit, open_output_file
    implicit none
    integer :: i, j, k, count, ifield, unit, in_file
    
    character(len=100), intent(in) :: filename
    integer  :: ierr

    namelist /force/ driven, elsasser, fix_power, ene_inj, xhl_inj, res_inj, kmin, kmax, nfields
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

    count = 0
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          if(       kx(i)**2 >= kmin(1)**2 .and. kx(i)**2 <= kmax(1)**2 &
              .and. ky(j)**2 >= kmin(2)**2 .and. ky(j)**2 <= kmax(2)**2 &
              .and. kz(k)**2 >= kmin(3)**2 .and. kz(k)**2 <= kmax(3)**2 &
            ) then
            count = count + 1
          endif
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
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if(     kx(i)**2 >= kmin(1)**2 .and. kx(i)**2 <= kmax(1)**2 &
              .and. ky(j)**2 >= kmin(2)**2 .and. ky(j)**2 <= kmax(2)**2 &
              .and. kz(k)**2 >= kmin(3)**2 .and. kz(k)**2 <= kmax(3)**2 &
              ) then

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

    if(proc0) call init_ranf(.true.,init_seed)
    call broadcast (init_seed)

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

    do istir = 1, nk_stir
      do ifield = 1, nfields

        w_stir(ifield, istir) = lz/(2.d0*pi)*abs(kz_stir(ifield, istir))*f_frequency(ifield, istir)
        
        sigma = sqrt(12.d0*abs(aimag(w_stir(ifield, istir)))/dt)*f_amplitude(ifield, istir)
        fa = cmplx(ranf() - 0.5d0, ranf() - 0.5d0) * sigma
        fb = cmplx(ranf() - 0.5d0, ranf() - 0.5d0) * sigma

        a_force(ifield, istir) = a_force(ifield, istir)*exp(-zi*w_stir(ifield, istir)*dt) + fa*dt
        b_force(ifield, istir) = 0.d0

        call get_rnd_seed(fin_seed)
      end do
    end do

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
    use grid, only: nkx, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_force
    implicit none
    character(*) :: name
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: u
    integer :: i, j, k, ifield, istir

    if (proc0) call put_time_stamp(timer_force)

    ifield = 0
    do i = 1, nfields
      if(trim(name) == field_names(i)) ifield = i
    enddo

    if (ifield == 0) return
    
    do istir = 1, nk_stir
      j = ky_stir(ifield, istir) + 1

      if(kx_stir(ifield, istir) >= 0) then
        i = kx_stir(ifield, istir) + 1
      else
        i = nkx + kx_stir(ifield, istir) + 1
      endif

      if(kz_stir(ifield, istir) >= 0) then
        k = kz_stir(ifield, istir) + 1
      else
        k = nkz + kz_stir(ifield, istir) + 1
      endif

      if(      (i >= ikx_st .and. i <= ikx_en) &
         .and. (j >= iky_st .and. j <= iky_en) &
         .and. (k >= ikz_st .and. k <= ikz_en) ) then

         u(i, k, j) = (a_force(ifield, istir) + b_force(ifield, istir))/sqrt(2.d0)
      endif
    end do

    if (proc0) call put_time_stamp(timer_force)

  end subroutine get_force


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

