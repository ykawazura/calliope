!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    25 Sep 2019
!! @brief   Time stepping for MHD_INCOMP
!-----------------------------------------------!
module advance
  use p3dfft
  use fields, only: nfields
  use fields, only: ibx, iby, ibz
  implicit none

  public solve

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: bx_new
  complex(r8), allocatable, dimension(:,:,:)   :: by_new
  complex(r8), allocatable, dimension(:,:,:)   :: bz_new
  complex(r8), allocatable, dimension(:,:,:,:) :: flx, exp_terms
  complex(r8), allocatable, dimension(:,:,:)   :: nbl2inv_div_b ! nabla^-2 (div b)
  real   (r8), allocatable, dimension(:,:,:)   :: de2k2_plus_1

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: bx_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: by_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: bz_tmp
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: bx_old2
  complex(r8), allocatable, dimension(:,:,:)   :: by_old2
  complex(r8), allocatable, dimension(:,:,:)   :: bz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fbx_old2, fby_old2, fbz_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2

  real   (r8) :: cflx, cfly, cflz
  integer :: max_vel_unit

  ! Backward FFT variables
  integer, parameter :: nbtran = 6
  integer, parameter :: iux = 4, iuy = 5, iuz = 6
  ! Forward FFT variables
  integer, parameter :: nftran = 3
  integer, parameter :: iflx_bx  = 1, iflx_by  = 2, iflx_bz  = 3 ! b^* x u

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: bx, by, bz
    use fields, only: bx_old, by_old, bz_old
    use grid, only: k2_max
    use grid, only: kx, ky, kz, k2, k2inv
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, cfl, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: eta, eta_exp, zi, nonlinear
    use params, only: dealias_scheme => dealias
    use dealias, only: filter
    use force, only: fbx, fbz, fby, fbx_old, fby_old, fbz_old, driven, update_force, get_force
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use advance_common, only: gear1, gear2, gear3
    use params, only: series_output
    use params, only: time_step_scheme
    implicit none
    real(r8) :: imp_terms_tintg0(nfields), imp_terms_tintg1(nfields)
    real(r8) :: imp_terms_tintg2(nfields), imp_terms_tintg3(nfields)  
    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_advance)

    ! initialize tmp fields
    if(counter == 0) then
      call init_multistep_fields

      cflx = maxval(abs(kx))/cfl
      cfly = maxval(abs(ky))/cfl
      cflz = maxval(abs(kz))/cfl

      counter = 1
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      ! Calcualte force terms
      if (driven) then
        ! at n
        call get_force('bx', fbx)
        call get_force('by', fby)
        call get_force('bz', fbz)
      endif

      !---------------  RK 1st step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(bx, by, bz, .true.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               bx(i,k,j), by(i,k,j), bz(i,k,j), &
                               flx(i,k,j,:), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kx(i), ky(j), kz(k), de2k2_plus_1(i, k, j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(ibx),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)

            ! update b
            call eSSPIFRK1(bx_tmp(i,k,j), bx(i,k,j), &
               exp_terms(i,k,j,ibx), &
               imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
            )
            call eSSPIFRK1(by_tmp(i,k,j), by(i,k,j), &
               exp_terms(i,k,j,iby), &
               imp_terms_tintg2(iby), imp_terms_tintg0(iby) &
            )
            call eSSPIFRK1(bz_tmp(i,k,j), bz(i,k,j), &
               exp_terms(i,k,j,ibz), &
               imp_terms_tintg2(ibz), imp_terms_tintg0(ibz) &
            )

          enddo
        enddo
      enddo
      !$omp end parallel do

      ! save explicit terms at the previous step
      !$omp workshare
      exp_terms0 = exp_terms
      !$omp end workshare

      !---------------  RK 2nd step  ---------------
      ! Calcualte force terms
      if (driven) then
        ! at n + 2/3
        call update_force(2.d0/3.d0*dt)
        call get_force('bx', fbx)
        call get_force('by', fby)
        call get_force('bz', fbz)

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(bx_tmp, by_tmp, bz_tmp, .false.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               bx_tmp(i,k,j), by_tmp(i,k,j), bz_tmp(i,k,j), &
                               flx(i,k,j,:), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kx(i), ky(j), kz(k), de2k2_plus_1(i, k, j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(ibx),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)

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

      !---------------  RK 3rd step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(bx_tmp, by_tmp, bz_tmp, .false.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               bx_tmp(i,k,j), by_tmp(i,k,j), bz_tmp(i,k,j), &
                               flx(i,k,j,:), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kx(i), ky(j), kz(k), de2k2_plus_1(i, k, j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(ibx),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibx),           dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(iby),           dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz),         0.d0, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz), 2.d0/3.d0*dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibz),           dt, k2(i,k,j), de2k2_plus_1(i,k,j), eta, eta_exp)

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
      bx_old = bx
      by_old = by
      bz_old = bz
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fbx_old = fbx
        fby_old = fby
        fbz_old = fbz
        !$omp end workshare
      endif
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      ! Calcualte force terms
      if (driven) then
        call update_force(dt)
        call get_force('bx', fbx)
        call get_force('by', fby)
        call get_force('bz', fbz)
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(bx, by, bz, .true.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               bx(i,k,j), by(i,k,j), bz(i,k,j), &
                               flx(i,k,j,:), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kx(i), ky(j), kz(k), de2k2_plus_1(i, k, j))

            ! 1st order 
            if(counter == 1) then
              ! update b
              call gear1(bx_new(i,k,j), bx(i,k,j), &
                 exp_terms(i,k,j,ibx), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )
              call gear1(by_new(i,k,j), by(i,k,j), &
                 exp_terms(i,k,j,iby), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )
              call gear1(bz_new(i,k,j), bz(i,k,j), &
                 exp_terms(i,k,j,ibz), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )

            ! 2nd order 
            elseif(counter == 2) then
              ! update b
              call gear2(bx_new(i,k,j), bx(i,k,j), bx_old(i,k,j), &
                 exp_terms    (i,k,j,ibx), &
                 exp_terms_old(i,k,j,ibx), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )
              call gear2(by_new(i,k,j), by(i,k,j), by_old(i,k,j), &
                 exp_terms    (i,k,j,iby), &
                 exp_terms_old(i,k,j,iby), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )
              call gear2(bz_new(i,k,j), bz(i,k,j), bz_old(i,k,j), &
                 exp_terms    (i,k,j,ibz), &
                 exp_terms_old(i,k,j,ibz), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )

            ! 3rd order 
            else
              ! update b
              call gear3(bx_new(i,k,j), bx(i,k,j), bx_old(i,k,j), bx_old2(i,k,j), &
                 exp_terms     (i,k,j,ibx), &
                 exp_terms_old (i,k,j,ibx), &
                 exp_terms_old2(i,k,j,ibx), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )
              call gear3(by_new(i,k,j), by(i,k,j), by_old(i,k,j), by_old2(i,k,j), &
                 exp_terms     (i,k,j,iby), &
                 exp_terms_old (i,k,j,iby), &
                 exp_terms_old2(i,k,j,iby), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )
              call gear3(bz_new(i,k,j), bz(i,k,j), bz_old(i,k,j), bz_old2(i,k,j), &
                 exp_terms     (i,k,j,ibz), &
                 exp_terms_old (i,k,j,ibz), &
                 exp_terms_old2(i,k,j,ibz), &
                 eta*(k2(i,k,j)/k2_max)**eta_exp/de2k2_plus_1(i, k, j) &
              )
            endif
          enddo
        enddo
      enddo
      !$omp end parallel do

      if(counter <= 2) counter = counter + 1
      
      ! save fields at the previous steps
      !$omp workshare
      bx_old2 = bx_old 
      bx_old  = bx

      by_old2 = by_old 
      by_old  = by

      bz_old2 = bz_old 
      bz_old  = bz

      exp_terms_old2 = exp_terms_old
      exp_terms_old  = exp_terms
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fbx_old2 = fbx_old 
        fbx_old  = fbx

        fby_old2 = fby_old 
        fby_old  = fby

        fbz_old2 = fbz_old 
        fbz_old  = fbz
        !$omp end workshare
      endif
    endif

    !$omp workshare
    bx = bx_new
    by = by_new
    bz = bz_new
    !$omp end workshare

    ! Dealiasing
    if(trim(dealias_scheme) /= '2/3') then
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            bx(i,k,j) = bx_new(i,k,j)*filter(i,k,j)
            by(i,k,j) = by_new(i,k,j)*filter(i,k,j)
            bz(i,k,j) = bz_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
    endif

    tt  = tt  + dt

    ! Div b cleaing
    allocate(nbl2inv_div_b(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          nbl2inv_div_b(i,k,j) = -zi*(   kx(i)*bx(i,k,j) &
                                       + ky(j)*by(i,k,j) &
                                       + kz(k)*bz(i,k,j) )*k2inv(i,k,j)

          bx(i,k,j) = bx(i,k,j) - zi*kx(i)*nbl2inv_div_b(i,k,j)
          by(i,k,j) = by(i,k,j) - zi*ky(j)*nbl2inv_div_b(i,k,j)
          bz(i,k,j) = bz(i,k,j) - zi*kz(k)*nbl2inv_div_b(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do
    deallocate(nbl2inv_div_b)

    if(series_output) call output_series_modes

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Allocate tmp fields for a multi 
!!          timestep method
!-----------------------------------------------!
  subroutine init_multistep_fields
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en, k2
    use file, only: open_output_file
    use params, only: time_step_scheme
    use params, only: de2
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))

    allocate(bx_new, source=src)
    allocate(by_new, source=src)
    allocate(bz_new, source=src)

    allocate(flx      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran )); flx       = 0.d0
    allocate(exp_terms(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms = 0.d0

    allocate(de2k2_plus_1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); de2k2_plus_1 = 0.d0

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do i = ikx_st, ikx_en
        do k = ikz_st, ikz_en
          de2k2_plus_1(i, k, j) = 1.d0 + de2*k2(i, k, j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(bx_tmp, source=src)
      allocate(by_tmp, source=src)
      allocate(bz_tmp, source=src)

      allocate(exp_terms0(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields )); exp_terms0  = 0.d0
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      allocate( bx_old2, source=src)
      allocate( by_old2, source=src)
      allocate( bz_old2, source=src)

      allocate(fbx_old2, source=src)
      allocate(fby_old2, source=src)
      allocate(fbz_old2, source=src)

      allocate(exp_terms_old (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old2 = 0.d0
    endif

    deallocate(src)

    call open_output_file (max_vel_unit, 'max_vel.dat')

  end subroutine init_multistep_fields


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Calculate nonlinear terms via
!!          1. Calculate grad in Fourier space
!!          2. Inverse FFT
!!          3. Calculate nonlinear terms 
!!             in real space
!!          4. Forward FFT
!-----------------------------------------------!
  subroutine get_nonlinear_terms(bx, by, bz, dt_reset)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: kx, ky, kz
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: bx, by, bz

    complex(r8), allocatable, dimension(:,:,:,:) :: wbk
    real   (r8), allocatable, dimension(:,:,:,:) :: wb , wf 
    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (r8) :: max_vel, dt_cfl, dt_digit

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    allocate(wbk(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nbtran), source=(0.d0, 0.d0))
    allocate(wb (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nbtran), source=0.d0)
    allocate(wf (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nftran), source=0.d0)

    ! 1. Calculate grad in Fourier space
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          ! b^* = (1 + de^2 nabla*2)b
          wbk(i,k,j,ibx) = bx(i,k,j)*de2k2_plus_1(i, k, j)
          wbk(i,k,j,iby) = by(i,k,j)*de2k2_plus_1(i, k, j)
          wbk(i,k,j,ibz) = bz(i,k,j)*de2k2_plus_1(i, k, j)

          ! u = -curl b : electron velocity
          wbk(i,k,j,iux) = -zi*(  ky(j)*bz(i,k,j) &
                                - kz(k)*by(i,k,j) )
                                                      
          wbk(i,k,j,iuy) = -zi*(  kz(k)*bx(i,k,j) &
                                - kx(i)*bz(i,k,j) )
                                                      
          wbk(i,k,j,iuz) = -zi*(  kx(i)*by(i,k,j) &
                                - ky(j)*bx(i,k,j) )
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 2. Inverse FFT
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_btran_c2r_many(wbk, nk_local_tot, wb, nl_local_tot, nbtran, 'tff')
    if (proc0) call put_time_stamp(timer_fft)

    ! (get max_vel for dt reset)
    if(dt_reset) then
      max_vel = max( &
                maxval(abs(wb(:,:,:,iux)))*cflx, maxval(abs(wb(:,:,:,iuy)))*cfly, maxval(abs(wb(:,:,:,iuz)))*cflz &
              )
      call max_allreduce(max_vel)
      dt_cfl = 1.d0/max_vel
      if(proc0) then
        write (unit=max_vel_unit, fmt="(100es30.21)") tt, max_vel
        call flush(max_vel_unit) 
      endif

      if(dt_cfl < dt) then
        if(proc0) then
          print *
          write (*, '("dt is decreased from ", es12.4e3)', advance='no') dt
        endif

        dt_digit = (log10(dt_cfl)/abs(log10(dt_cfl)))*ceiling(abs(log10(dt_cfl)))
        ! dt = floor(dt_cfl*10.d0**(-dt_digit))*10.d0**dt_digit

        if (reset_method == 'multiply') then
          dt = 0.5d0*dt
        elseif (reset_method == 'decrement') then
          dt_digit = (log10(dt)/abs(log10(dt)))*ceiling(abs(log10(dt)))

          ! when dt = 0.0**01***
          if (dt*10.d0**(-dt_digit) - 1.0d0 < 1.0d0) then
            dt = 0.9d0*10.d0**dt_digit
          else
            dt = (dt*10.d0**(-dt_digit) - 1.0d0)*10.d0**dt_digit
          endif
        endif

        counter = 1

        if(proc0) then
          print '("  to ", es12.4e3)', dt
          print *
        endif
      endif
      if(dt_cfl > 0.d0 .and. dt_cfl > increase_dt .and. dt < increase_dt) then
        if(proc0) then
          print *
          write (*, '("dt is increased from ", es12.4e3)', advance='no') dt
        endif

        dt = increase_dt

        counter = 1

        if(proc0) then
          print '("  to ", es12.4e3)', dt
          print *
        endif
      endif
    endif

    ! 3. Calculate nonlinear terms in real space
    !$omp parallel do private(j, k) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          wf(j,k,i,iflx_bx) = wb(j,k,i,iby)*wb(j,k,i,iuz) - wb(j,k,i,ibz)*wb(j,k,i,iuy)
          wf(j,k,i,iflx_by) = wb(j,k,i,ibz)*wb(j,k,i,iux) - wb(j,k,i,ibx)*wb(j,k,i,iuz)
          wf(j,k,i,iflx_bz) = wb(j,k,i,ibx)*wb(j,k,i,iuy) - wb(j,k,i,iby)*wb(j,k,i,iux)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 4. Forward FFT
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_ftran_r2c_many(wf, nl_local_tot, flx, nk_local_tot, nftran, 'fft')
    if (proc0) call put_time_stamp(timer_fft)
    !$omp workshare
    flx = flx/nlx/nly/nlz
    !$omp end workshare

    deallocate(wbk)
    deallocate(wb )
    deallocate(wf )

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms, &
                           bx, by, bz, &
                           flx, &
                           fbx, fby, fbz, &
                           kx, ky, kz, de2k2_plus_1)
    use params, only: zi
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: bx, by, bz
    complex(r8), intent(in ) :: flx(nftran)
    complex(r8), intent(in ) :: fbx, fby, fbz
    real(r8)   , intent(in)  :: kx, ky, kz, de2k2_plus_1
    complex(r8)              :: nl(nfields)

    ! curl (b^* x u)
    nl(ibx) = -zi*(  ky*flx(iflx_bz ) - kz*flx(iflx_by ) )
    nl(iby) = -zi*(  kz*flx(iflx_bx ) - kx*flx(iflx_bz ) )
    nl(ibz) = -zi*(  kx*flx(iflx_by ) - ky*flx(iflx_bx ) )

    exp_terms(ibx) = (nl(ibx) + fbx)/de2k2_plus_1
    exp_terms(iby) = (nl(iby) + fby)/de2k2_plus_1
    exp_terms(ibz) = (nl(ibz) + fbz)/de2k2_plus_1

  end subroutine get_ext_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Time integral of hyperdissipation
!-----------------------------------------------!
  subroutine get_imp_terms_tintg(imp_terms_tintg, t, k2, de2k2_plus_1, coeff, nexp)
    use grid, only: k2_max
    implicit none
    real(r8), intent(out) :: imp_terms_tintg
    real(r8), intent(in) :: t, k2, de2k2_plus_1, coeff
    integer, intent(in) :: nexp

    imp_terms_tintg = -coeff*(k2/k2_max)**nexp*t/de2k2_plus_1
  end subroutine get_imp_terms_tintg


!-----------------------------------------------!
!> @author  YK
!! @date    15 Jul 2021
!! @brief   Output series modes
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: bx, by, bz
    use params, only: n_series_modes, series_modes
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, sum_reduce
    use time, only: tt
    use diagnostics_common, only: series_modes_unit
    implicit none
    complex(r8), dimension(n_series_modes) :: bx_modes, by_modes, bz_modes
    integer :: n, i, j, k

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

        bx_modes(n) = bx(i, k, j)
        by_modes(n) = by(i, k, j)
        bz_modes(n) = bz(i, k, j)

      endif
    enddo

    call sum_reduce(bx_modes, 0)
    call sum_reduce(by_modes, 0)
    call sum_reduce(bz_modes, 0)

    do n = 1, n_series_modes
      if(proc0) then
        i = series_modes(n, 1)
        j = series_modes(n, 2)
        k = series_modes(n, 3)
999 format(es30.21, A6, 5es30.21e3)
        write (unit=series_modes_unit, fmt=999) tt, 'bx', kx(i), ky(j), kz(k), bx_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'by', kx(i), ky(j), kz(k), by_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bz', kx(i), ky(j), kz(k), bz_modes(n)
        call flush(series_modes_unit) 

      endif
    enddo

  end subroutine output_series_modes

end module advance

