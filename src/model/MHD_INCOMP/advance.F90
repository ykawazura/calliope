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
  use fields, only: iux, iuy, iuz
  use fields, only: ibx, iby, ibz
  use mp, only: proc0, max_allreduce, sum_allreduce
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: ux_new
  complex(r8), allocatable, dimension(:,:,:)   :: uy_new
  complex(r8), allocatable, dimension(:,:,:)   :: uz_new
  complex(r8), allocatable, dimension(:,:,:)   :: bx_new
  complex(r8), allocatable, dimension(:,:,:)   :: by_new
  complex(r8), allocatable, dimension(:,:,:)   :: bz_new
  complex(r8), allocatable, dimension(:,:,:,:) :: w
  real   (r8), allocatable, dimension(:,:,:,:) :: w_r, flx_r 
  complex(r8), allocatable, dimension(:,:,:,:) :: flx, exp_terms
  real   (r8), allocatable, dimension(:,:)     :: kxt

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: ux_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: uy_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: uz_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: bx_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: by_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: bz_tmp
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: ux_old2
  complex(r8), allocatable, dimension(:,:,:)   :: uy_old2
  complex(r8), allocatable, dimension(:,:,:)   :: uz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: bx_old2
  complex(r8), allocatable, dimension(:,:,:)   :: by_old2
  complex(r8), allocatable, dimension(:,:,:)   :: bz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fux_old2, fuy_old2, fuz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fbx_old2, fby_old2, fbz_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2
  real   (r8), allocatable, dimension(:,:)     :: kxt_old, kxt_old2

  integer :: cfl_unit, rms_unit

  ! Forward FFT variables
  integer, parameter :: nftran = 9
  integer, parameter :: iflx_uxx = 1, iflx_uxy = 2, iflx_uxz = 3 ! uu - bb 
  integer, parameter ::               iflx_uyy = 4, iflx_uyz = 5 ! uu - bb 
  integer, parameter ::                             iflx_uzz = 6 ! uu - bb 
  integer, parameter :: iflx_bx  = 7, iflx_by  = 8, iflx_bz  = 9 ! b x u

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use fields, only: ux_old, uy_old, uz_old
    use fields, only: bx_old, by_old, bz_old
    use grid, only: k2_max
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nu, nu_h, eta_h, nu_h_exp, eta, eta_h_exp, zi, nonlinear, shear, q
    use dealias, only: filter
    use shearing_box, only: shear_flg, k2t, k2t_inv, tsc, tremap
    use force, only: elsasser, driven, update_force, get_force, normalize_force, normalize_force_els
    use force, only: fux, fuy, fuz, fux_old, fuy_old, fuz_old
    use force, only: fbx, fby, fbz, fbx_old, fby_old, fbz_old
    use force, only: fzpx, fzpy, fzpz, fzmx, fzmy, fzmz
    use advance_common, only: gear1, gear2, gear3
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use diagnostics_common, only: series_output
    use params, only: time_step_scheme
    implicit none
    real(r8) :: imp_terms_tintg0(nfields), imp_terms_tintg1(nfields)
    real(r8) :: imp_terms_tintg2(nfields), imp_terms_tintg3(nfields)  
    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_advance)

    ! initialize tmp fields
    if(counter == 0) then
      call init_work_fields
      counter = 1
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      ! Calcualte force terms
      if (driven) then
        ! at n
        if(elsasser) then
          call get_force('zpx', fzpx)
          call get_force('zpy', fzpy)
          call get_force('zpz', fzpz)
          call get_force('zmx', fzmx)
          call get_force('zmy', fzmy)
          call get_force('zmz', fzmz)
          call div_free_force(fzpx, fzpy, fzpz)
          call div_free_force(fzmx, fzmy, fzmz)
          call normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          fux = fzpx + fzmx
          fuy = fzpy + fzmy
          fuz = fzpz + fzmz
          fbx = fzpx - fzmx
          fby = fzpy - fzmy
          fbz = fzpz - fzmz
        else
          call get_force('ux', fux)
          call get_force('uy', fuy)
          call get_force('uz', fuz)
          call get_force('bx', fbx)
          call get_force('by', fby)
          call get_force('bz', fbz)
          call div_free_force(fux, fuy, fuz)
          call div_free_force(fbx, fby, fbz)
          call normalize_force(fux, fuy, fuz, fbx, fby, fbz)
        endif
      endif

      !---------------  RK 1st step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux, uy, uz, bx, by, bz, .true.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg1)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               ux(i,k,j), uy(i,k,j), uz(i,k,j), bx(i,k,j), by(i,k,j), bz(i,k,j), &
                               flx(i,k,j,:), &
                               fux(i,k,j), fuy(i,k,j), fuz(i,k,j), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iux), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iux), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuy), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iuy), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuz), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iuz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
                                                                                                          
            call get_imp_terms_tintg(imp_terms_tintg0(ibx), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ibx), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(iby), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ibz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update u
            call eSSPIFRK1(ux_tmp(i,k,j), ux(i,k,j), &
               exp_terms(i,k,j,iux), &
               imp_terms_tintg1(iux), imp_terms_tintg0(iux) &
            )
            call eSSPIFRK1(uy_tmp(i,k,j), uy(i,k,j), &
               exp_terms(i,k,j,iuy), &
               imp_terms_tintg1(iuy), imp_terms_tintg0(iuy) &
            )
            call eSSPIFRK1(uz_tmp(i,k,j), uz(i,k,j), &
               exp_terms(i,k,j,iuz), &
               imp_terms_tintg1(iuz), imp_terms_tintg0(iuz) &
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
      !$omp parallel do private(j, k, i)
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            ux_tmp(i,k,j) = ux_tmp(i,k,j)*filter(i,k,j)
            uy_tmp(i,k,j) = uy_tmp(i,k,j)*filter(i,k,j)
            uz_tmp(i,k,j) = uz_tmp(i,k,j)*filter(i,k,j)
            bx_tmp(i,k,j) = bx_tmp(i,k,j)*filter(i,k,j)
            by_tmp(i,k,j) = by_tmp(i,k,j)*filter(i,k,j)
            bz_tmp(i,k,j) = bz_tmp(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
      !$omp end parallel do

      !---------------  RK 2nd step  ---------------
      ! Calcualte force terms
      if (driven) then
        ! at n + 2/3
        call update_force(2.d0/3.d0*dt)

        if(elsasser) then
          call get_force('zpx', fzpx)
          call get_force('zpy', fzpy)
          call get_force('zpz', fzpz)
          call get_force('zmx', fzmx)
          call get_force('zmy', fzmy)
          call get_force('zmz', fzmz)
          call div_free_force(fzpx, fzpy, fzpz)
          call div_free_force(fzmx, fzmy, fzmz)
          call normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          fux = fzpx + fzmx
          fuy = fzpy + fzmy
          fuz = fzpz + fzmz
          fbx = fzpx - fzmx
          fby = fzpy - fzmy
          fbz = fzpz - fzmz
        else
          call get_force('ux', fux)
          call get_force('uy', fuy)
          call get_force('uz', fuz)
          call get_force('bx', fbx)
          call get_force('by', fby)
          call get_force('bz', fbz)
          call div_free_force(fux, fuy, fuz)
          call div_free_force(fbx, fby, fbz)
          call normalize_force(fux, fuy, fuz, fbx, fby, fbz)
        endif

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      ! Calcualte kxt at n + 2/3
      if(shear) then
        !$omp parallel do private(i, k) schedule(static)
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            kxt(i,j) = kx(i) + q*shear_flg*(tsc + 2.d0/3.d0*dt)*ky(j)
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

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               ux_tmp(i,k,j), uy_tmp(i,k,j), uz_tmp(i,k,j), bx_tmp(i,k,j), by_tmp(i,k,j), bz_tmp(i,k,j), &
                               flx(i,k,j,:), &
                               fux(i,k,j), fuy(i,k,j), fuz(i,k,j), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iux), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iux), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuy), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuy), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuz), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
                                                                                                          
            call get_imp_terms_tintg(imp_terms_tintg0(ibx), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update u
            call eSSPIFRK2(ux_tmp(i,k,j), ux_tmp(i,k,j), ux(i,k,j), &
               exp_terms(i,k,j,iux), &
               imp_terms_tintg2(iux), imp_terms_tintg0(iux) &
            )
            call eSSPIFRK2(uy_tmp(i,k,j), uy_tmp(i,k,j), uy(i,k,j), &
               exp_terms(i,k,j,iuy), &
               imp_terms_tintg2(iuy), imp_terms_tintg0(iuy) &
            )
            call eSSPIFRK2(uz_tmp(i,k,j), uz_tmp(i,k,j), uz(i,k,j), &
               exp_terms(i,k,j,iuz), &
               imp_terms_tintg2(iuz), imp_terms_tintg0(iuz) &
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
      !$omp parallel do private(j, k, i)
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            ux_tmp(i,k,j) = ux_tmp(i,k,j)*filter(i,k,j)
            uy_tmp(i,k,j) = uy_tmp(i,k,j)*filter(i,k,j)
            uz_tmp(i,k,j) = uz_tmp(i,k,j)*filter(i,k,j)
            bx_tmp(i,k,j) = bx_tmp(i,k,j)*filter(i,k,j)
            by_tmp(i,k,j) = by_tmp(i,k,j)*filter(i,k,j)
            bz_tmp(i,k,j) = bz_tmp(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
      !$omp end parallel do

      !---------------  RK 3rd step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux_tmp, uy_tmp, uz_tmp, bx_tmp, by_tmp, bz_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               ux_tmp(i,k,j), uy_tmp(i,k,j), uz_tmp(i,k,j), bx_tmp(i,k,j), by_tmp(i,k,j), bz_tmp(i,k,j), &
                               flx(i,k,j,:), &
                               fux(i,k,j), fuy(i,k,j), fuz(i,k,j), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iux), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iux), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iux), tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuy), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuy), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iuy), tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(iuz), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iuz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iuz), tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
                                                                                                          
            call get_imp_terms_tintg(imp_terms_tintg0(ibx), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibx), tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(iby), tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibz), tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update u
            call eSSPIFRK3(ux_new(i,k,j), ux_tmp(i,k,j), ux(i,k,j), &
               exp_terms(i,k,j,iux), exp_terms0(i,k,j,iux), &
               imp_terms_tintg3(iux), imp_terms_tintg2(iux), imp_terms_tintg0(iux) &
            )
            call eSSPIFRK3(uy_new(i,k,j), uy_tmp(i,k,j), uy(i,k,j), &
               exp_terms(i,k,j,iuy), exp_terms0(i,k,j,iuy), &
               imp_terms_tintg3(iuy), imp_terms_tintg2(iuy), imp_terms_tintg0(iuy) &
            )
            call eSSPIFRK3(uz_new(i,k,j), uz_tmp(i,k,j), uz(i,k,j), &
               exp_terms(i,k,j,iuz), exp_terms0(i,k,j,iuz), &
               imp_terms_tintg3(iuz), imp_terms_tintg2(iuz), imp_terms_tintg0(iuz) &
            )

            ! update b
            call eSSPIFRK3(bx_new(i,k,j), bx_tmp(i,k,j), bx(i,k,j), &
               exp_terms(i,k,j,ibx), exp_terms0(i,k,j,ibx), &
               imp_terms_tintg3(ibx), imp_terms_tintg2(ibx), imp_terms_tintg0(ibx) &
            )
            call eSSPIFRK3(by_new(i,k,j), by_tmp(i,k,j), by(i,k,j), &
               exp_terms (i,k,j,iby), exp_terms0(i,k,j,iby), &
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
      ux_old = ux
      uy_old = uy
      uz_old = uz
      bx_old = bx
      by_old = by
      bz_old = bz
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fux_old = fux
        fuy_old = fuy
        fuz_old = fuz
        fbx_old = fbx
        fby_old = fby
        fbz_old = fbz
        !$omp end workshare
      endif

      ! Dealiasing
      !$omp parallel do private(j, k, i)
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            ux_new(i,k,j) = ux_new(i,k,j)*filter(i,k,j)
            uy_new(i,k,j) = uy_new(i,k,j)*filter(i,k,j)
            uz_new(i,k,j) = uz_new(i,k,j)*filter(i,k,j)
            bx_new(i,k,j) = bx_new(i,k,j)*filter(i,k,j)
            by_new(i,k,j) = by_new(i,k,j)*filter(i,k,j)
            bz_new(i,k,j) = bz_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
      !$omp end parallel do
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      ! Calcualte force terms
      if (driven) then
        call update_force(dt)

        if(elsasser) then
          call get_force('zpx', fzpx)
          call get_force('zpy', fzpy)
          call get_force('zpz', fzpz)
          call get_force('zmx', fzmx)
          call get_force('zmy', fzmy)
          call get_force('zmz', fzmz)
          call div_free_force(fzpx, fzpy, fzpz)
          call div_free_force(fzmx, fzmy, fzmz)
          call normalize_force_els(fzpx, fzpy, fzpz, fzmx, fzmy, fzmz)
          fux = fzpx + fzmx
          fuy = fzpy + fzmy
          fuz = fzpz + fzmz
          fbx = fzpx - fzmx
          fby = fzpy - fzmy
          fbz = fzpz - fzmz
        else
          call get_force('ux', fux)
          call get_force('uy', fuy)
          call get_force('uz', fuz)
          call get_force('bx', fbx)
          call get_force('by', fby)
          call get_force('bz', fbz)
          call div_free(fux, fuy, fuz)
          call div_free(fbx, fby, fbz)
          call normalize_force(fux, fuy, fuz, fbx, fby, fbz)
        endif
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(ux, uy, uz, bx, by, bz, .true.)

      !$omp parallel do private(j, k, i)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               ux(i,k,j), uy(i,k,j), uz(i,k,j), bx(i,k,j), by(i,k,j), bz(i,k,j), &
                               flx(i,k,j,:), &
                               fux(i,k,j), fuy(i,k,j), fuz(i,k,j), &
                               fbx(i,k,j), fby(i,k,j), fbz(i,k,j), &
                               kxt(i,j), ky(j), kz(k), k2t_inv(i,k,j))

            ! 1st order 
            if(counter == 1) then
              ! update u
              call gear1(ux_new(i,k,j), ux(i,k,j), &
                 exp_terms(i,k,j,iux), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear1(uy_new(i,k,j), uy(i,k,j), &
                 exp_terms(i,k,j,iuy), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear1(uz_new(i,k,j), uz(i,k,j), &
                 exp_terms(i,k,j,iuz), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )

              ! update b
              call gear1(bx_new(i,k,j), bx(i,k,j), &
                 exp_terms(i,k,j,ibx), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )
              call gear1(by_new(i,k,j), by(i,k,j), &
                 exp_terms(i,k,j,iby), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )
              call gear1(bz_new(i,k,j), bz(i,k,j), &
                 exp_terms(i,k,j,ibz), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )

            ! 2nd order 
            elseif(counter == 2) then
              ! update u
              call gear2(ux_new(i,k,j), ux(i,k,j), ux_old(i,k,j), &
                 exp_terms    (i,k,j,iux), &
                 exp_terms_old(i,k,j,iux), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear2(uy_new(i,k,j), uy(i,k,j), uy_old(i,k,j), &
                 exp_terms    (i,k,j,iuy), &
                 exp_terms_old(i,k,j,iuy), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear2(uz_new(i,k,j), uz(i,k,j), uz_old(i,k,j), &
                 exp_terms    (i,k,j,iuz), &
                 exp_terms_old(i,k,j,iuz), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )

              ! update b
              call gear2(bx_new(i,k,j), bx(i,k,j), bx_old(i,k,j), &
                 exp_terms    (i,k,j,ibx), &
                 exp_terms_old(i,k,j,ibx), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )
              call gear2(by_new(i,k,j), by(i,k,j), by_old(i,k,j), &
                 exp_terms    (i,k,j,iby), &
                 exp_terms_old(i,k,j,iby), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )
              call gear2(bz_new(i,k,j), bz(i,k,j), bz_old(i,k,j), &
                 exp_terms    (i,k,j,ibz), &
                 exp_terms_old(i,k,j,ibz), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )

            ! 3rd order 
            else
              ! update u
              call gear3(ux_new(i,k,j), ux(i,k,j), ux_old(i,k,j), ux_old2(i,k,j), &
                 exp_terms     (i,k,j,iux), &
                 exp_terms_old (i,k,j,iux), &
                 exp_terms_old2(i,k,j,iux), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear3(uy_new(i,k,j), uy(i,k,j), uy_old(i,k,j), uy_old2(i,k,j), &
                 exp_terms     (i,k,j,iuy), &
                 exp_terms_old (i,k,j,iuy), &
                 exp_terms_old2(i,k,j,iuy), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear3(uz_new(i,k,j), uz(i,k,j), uz_old(i,k,j), uz_old2(i,k,j), &
                 exp_terms     (i,k,j,iuz), &
                 exp_terms_old (i,k,j,iuz), &
                 exp_terms_old2(i,k,j,iuz), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )

              ! update b
              call gear3(bx_new(i,k,j), bx(i,k,j), bx_old(i,k,j), bx_old2(i,k,j), &
                 exp_terms     (i,k,j,ibx), &
                 exp_terms_old (i,k,j,ibx), &
                 exp_terms_old2(i,k,j,ibx), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )
              call gear3(by_new(i,k,j), by(i,k,j), by_old(i,k,j), by_old2(i,k,j), &
                 exp_terms     (i,k,j,iby), &
                 exp_terms_old (i,k,j,iby), &
                 exp_terms_old2(i,k,j,iby), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )
              call gear3(bz_new(i,k,j), bz(i,k,j), bz_old(i,k,j), bz_old2(i,k,j), &
                 exp_terms     (i,k,j,ibz), &
                 exp_terms_old (i,k,j,ibz), &
                 exp_terms_old2(i,k,j,ibz), &
                 eta*(k2t(i,k,j)/k2_max) + eta_h*(k2t(i,k,j)/k2_max)**eta_h_exp &
              )
            endif
          enddo
        enddo
      enddo
      !$omp end parallel do

      if(counter <= 2) counter = counter + 1

      ! save fields at the previous steps
      !$omp workshare
      ux_old2 = ux_old 
      ux_old  = ux

      uy_old2 = uy_old 
      uy_old  = uy

      uz_old2 = uz_old 
      uz_old  = uz

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
        fux_old2 = fux_old 
        fux_old  = fux

        fuy_old2 = fuy_old 
        fuy_old  = fuy

        fuz_old2 = fuz_old 
        fuz_old  = fuz

        fbx_old2 = fbx_old 
        fbx_old  = fbx

        fby_old2 = fby_old 
        fby_old  = fby

        fbz_old2 = fbz_old 
        fbz_old  = fbz
        !$omp end workshare
      endif

      if(shear) then
        !$omp workshare
        kxt_old2 = kxt_old 
        kxt_old  = kxt
        !$omp end workshare
      endif

      ! Dealiasing
      !$omp parallel do private(j, k, i)
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            ux_new(i,k,j) = ux_new(i,k,j)*filter(i,k,j)
            uy_new(i,k,j) = uy_new(i,k,j)*filter(i,k,j)
            uz_new(i,k,j) = uz_new(i,k,j)*filter(i,k,j)
            bx_new(i,k,j) = bx_new(i,k,j)*filter(i,k,j)
            by_new(i,k,j) = by_new(i,k,j)*filter(i,k,j)
            bz_new(i,k,j) = bz_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
      !$omp end parallel do
    endif

    !$omp workshare
    ux = ux_new
    uy = uy_new
    uz = uz_new
    bx = bx_new
    by = by_new
    bz = bz_new
    !$omp end workshare

    tt  = tt  + dt
    tsc = tsc + dt

    if(shear .and. tsc > tremap) call remap

    if(shear) then
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do i = ikx_st, ikx_en
          kxt(i,j) = kx(i) + q*shear_flg*tsc*ky(j)
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

    ! Div u & b cleaing
    call div_free(ux, uy, uz)
    call div_free(bx, by, bz)

    if(series_output) call output_series_modes

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Allocate fields used only here
!-----------------------------------------------!
  subroutine init_work_fields
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: q
    use shearing_box, only: shear_flg, tsc, k2t, k2t_inv
    use file, only: open_output_file
    use params, only: time_step_scheme
    use mp, only: proc0
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src
    integer :: i, j, k

    allocate(kxt(ikx_st:ikx_en, iky_st:iky_en))

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do i = ikx_st, ikx_en
        kxt(i, j) = kx(i) + q*shear_flg*tsc*ky(j)
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

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))

      allocate( ux_old2, source=src)
      allocate( uy_old2, source=src)
      allocate( uz_old2, source=src)
      allocate( bx_old2, source=src)
      allocate( by_old2, source=src)
      allocate( bz_old2, source=src)

      allocate(fux_old2, source=src)
      allocate(fuy_old2, source=src)
      allocate(fuz_old2, source=src)
      allocate(fbx_old2, source=src)
      allocate(fby_old2, source=src)
      allocate(fbz_old2, source=src)

      allocate(exp_terms_old (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old2 = 0.d0

      allocate(kxt_old (ikx_st:ikx_en, iky_st:iky_en)); kxt_old  = kxt
      allocate(kxt_old2(ikx_st:ikx_en, iky_st:iky_en)); kxt_old2 = kxt

      deallocate(src)
    endif


    if(proc0) call open_output_file (cfl_unit, 'cfl.dat')
    if(proc0) call open_output_file (rms_unit, 'rms.dat')

  end subroutine init_work_fields


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
  subroutine get_nonlinear_terms(ux, uy, uz, bx, by, bz, dt_reset)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi
    use mp, only: proc0, max_allreduce, sum_allreduce
    use time, only: cfl, dt, tt, reset_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    use advance_common, only: dt_adjust_while_running 
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: ux, uy, uz, bx, by, bz

    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (r8) :: max_vel_x, max_vel_y, max_vel_z , dt_cfl
    real   (r8) :: ux_rms, uy_rms, uz_rms, bx_rms, by_rms, bz_rms

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          w(i,k,j,iux) = ux(i,k,j)
          w(i,k,j,iuy) = uy(i,k,j)
          w(i,k,j,iuz) = uz(i,k,j)
          w(i,k,j,ibx) = bx(i,k,j)
          w(i,k,j,iby) = by(i,k,j)
          w(i,k,j,ibz) = bz(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 1. Inverse FFT
    if (proc0) call put_time_stamp(timer_fft)
    ! for some reason c2r_many doesn't work at Fugaku
    !call p3dfft_btran_c2r_many(w, nk_local_tot, w_r, nl_local_tot, nfields, 'tff')
    do i = 1, nfields
      call p3dfft_btran_c2r(w(:,:,:,i), w_r(:,:,:,i), 'tff')
    enddo
    if (proc0) call put_time_stamp(timer_fft)


    ! write rms
    ux_rms = sum(w_r(:,:,:,iux)**2); uy_rms = sum(w_r(:,:,:,iuy)**2); uz_rms = sum(w_r(:,:,:,iuz)**2)
    bx_rms = sum(w_r(:,:,:,ibx)**2); by_rms = sum(w_r(:,:,:,iby)**2); bz_rms = sum(w_r(:,:,:,ibz)**2)

    call sum_allreduce(ux_rms); call sum_allreduce(uy_rms); call sum_allreduce(uz_rms)
    call sum_allreduce(bx_rms); call sum_allreduce(by_rms); call sum_allreduce(bz_rms)

    ux_rms = sqrt(ux_rms/nlx/nly/nlz); uy_rms = sqrt(uy_rms/nlx/nly/nlz); uz_rms = sqrt(uz_rms/nlx/nly/nlz)
    bx_rms = sqrt(bx_rms/nlx/nly/nlz); by_rms = sqrt(by_rms/nlx/nly/nlz); bz_rms = sqrt(bz_rms/nlx/nly/nlz)

    if(proc0) then
      write (unit=rms_unit, fmt="(100es30.21)") tt, ux_rms, uy_rms, uz_rms, bx_rms, by_rms, bz_rms
      call flush(rms_unit) 
    endif


    ! (get max_vel for dt reset)
    if(dt_reset) then
      max_vel_x = max( &
                    maxval(abs(w_r(:,:,:,iux) + w_r(:,:,:,ibx))), &
                    maxval(abs(w_r(:,:,:,iux) - w_r(:,:,:,ibx)))  &
                  )
      max_vel_y = max( &
                    maxval(abs(w_r(:,:,:,iuy) + w_r(:,:,:,iby))), &
                    maxval(abs(w_r(:,:,:,iuy) - w_r(:,:,:,iby)))  &
                  )
      max_vel_z = max( &
                    maxval(abs(w_r(:,:,:,iuz) + w_r(:,:,:,ibz))), &
                    maxval(abs(w_r(:,:,:,iuz) - w_r(:,:,:,ibz)))  &
                  )
      call max_allreduce(max_vel_x)
      call max_allreduce(max_vel_y)
      call max_allreduce(max_vel_z)
      dt_cfl = cfl*min(dlx/max_vel_x, dly/max_vel_y, dlz/max_vel_z)

      if(proc0) then
        write (unit=cfl_unit, fmt="(100es30.21)") tt, dt_cfl, max_vel_x, max_vel_y, max_vel_z
        call flush(cfl_unit) 
      endif

      call reset_dt(dt_cfl, counter)

    endif

    ! When the file 'dt_adjust' including a float number is created,
    ! dt will be manually adjusted to that value while running.
    call dt_adjust_while_running() 

    ! 2. Calculate nonlinear terms in real space
    !$omp parallel do private(j, k) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          flx_r(j,k,i,iflx_uxx) = w_r(j,k,i,iux)*w_r(j,k,i,iux) - w_r(j,k,i,ibx)*w_r(j,k,i,ibx)
          flx_r(j,k,i,iflx_uxy) = w_r(j,k,i,iux)*w_r(j,k,i,iuy) - w_r(j,k,i,ibx)*w_r(j,k,i,iby)
          flx_r(j,k,i,iflx_uxz) = w_r(j,k,i,iux)*w_r(j,k,i,iuz) - w_r(j,k,i,ibx)*w_r(j,k,i,ibz)
          flx_r(j,k,i,iflx_uyy) = w_r(j,k,i,iuy)*w_r(j,k,i,iuy) - w_r(j,k,i,iby)*w_r(j,k,i,iby)
          flx_r(j,k,i,iflx_uyz) = w_r(j,k,i,iuy)*w_r(j,k,i,iuz) - w_r(j,k,i,iby)*w_r(j,k,i,ibz)
          flx_r(j,k,i,iflx_uzz) = w_r(j,k,i,iuz)*w_r(j,k,i,iuz) - w_r(j,k,i,ibz)*w_r(j,k,i,ibz)

          flx_r(j,k,i,iflx_bx ) = w_r(j,k,i,iby)*w_r(j,k,i,iuz) - w_r(j,k,i,ibz)*w_r(j,k,i,iuy)
          flx_r(j,k,i,iflx_by ) = w_r(j,k,i,ibz)*w_r(j,k,i,iux) - w_r(j,k,i,ibx)*w_r(j,k,i,iuz)
          flx_r(j,k,i,iflx_bz ) = w_r(j,k,i,ibx)*w_r(j,k,i,iuy) - w_r(j,k,i,iby)*w_r(j,k,i,iux)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 3. Forward FFT
    if (proc0) call put_time_stamp(timer_fft)
    ! for some reason r2c_many doesn't work at Fugaku
    !call p3dfft_ftran_r2c_many(flx_r, nl_local_tot, flx, nk_local_tot, nftran, 'fft')
    do i = 1, nftran
      call p3dfft_ftran_r2c(flx_r(:,:,:,i), flx(:,:,:,i), 'fft')
    enddo
    if (proc0) call put_time_stamp(timer_fft)
    !$omp workshare
    flx = flx/nlx/nly/nlz
    !$omp end workshare

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms, &
                           ux, uy, uz, bx, by, bz, &
                           flx, &
                           fux, fuy, fuz, &
                           fbx, fby, fbz, &
                           kxt, ky, kz, k2t_inv)
    use params, only: zi, q
    use shearing_box, only: shear_flg
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: ux, uy, uz, bx, by, bz
    complex(r8), intent(in ) :: flx(nftran)
    complex(r8), intent(in ) :: fux, fuy, fuz
    complex(r8), intent(in ) :: fbx, fby, fbz
    real(r8)   , intent(in)  :: kxt, ky, kz, k2t_inv
    complex(r8)              :: nl(nfields)
    complex(r8)              :: p

    ! div (uu - bb)
    nl(iux) = -zi*( kxt*flx(iflx_uxx) + ky*flx(iflx_uxy) + kz*flx(iflx_uxz) )
    nl(iuy) = -zi*( kxt*flx(iflx_uxy) + ky*flx(iflx_uyy) + kz*flx(iflx_uyz) )
    nl(iuz) = -zi*( kxt*flx(iflx_uxz) + ky*flx(iflx_uyz) + kz*flx(iflx_uzz) )

    ! curl (b x u)
    nl(ibx) = -zi*( ky *flx(iflx_bz) - kz *flx(iflx_by) )
    nl(iby) = -zi*( kz *flx(iflx_bx) - kxt*flx(iflx_bz) )
    nl(ibz) = -zi*( kxt*flx(iflx_by) - ky *flx(iflx_bx) )

    ! get pressure
    p = -zi*( kxt*nl(iux) + ky*nl(iuy) + kz*nl(iuz) + 2.d0*zi*shear_flg*(kxt*uy - ky*ux) )*k2t_inv

    exp_terms(iux) = nl(iux) + fux - zi*kxt*p + 2.d0*shear_flg*uy
    exp_terms(iuy) = nl(iuy) + fuy - zi*ky *p - (2.d0 - q)*shear_flg*ux
    exp_terms(iuz) = nl(iuz) + fuz - zi*kz *p
    
    exp_terms(ibx) = nl(ibx) + fbx
    exp_terms(iby) = nl(iby) + fby - q*shear_flg*bx
    exp_terms(ibz) = nl(ibz) + fbz

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
!! @date    15 Jul 2021
!! @brief   Output series modes
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, sum_reduce
    use time, only: tt
    use diagnostics_common, only: n_series_modes, series_modes
    use diagnostics_common, only: series_modes_unit
    implicit none
    complex(r8), dimension(n_series_modes) :: ux_modes, uy_modes, uz_modes
    complex(r8), dimension(n_series_modes) :: bx_modes, by_modes, bz_modes
    integer :: n, i, j, k

    ux_modes(:) = 0.d0
    uy_modes(:) = 0.d0
    uz_modes(:) = 0.d0
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

        ux_modes(n) = ux(i, k, j)
        uy_modes(n) = uy(i, k, j)
        uz_modes(n) = uz(i, k, j)

        bx_modes(n) = bx(i, k, j)
        by_modes(n) = by(i, k, j)
        bz_modes(n) = bz(i, k, j)

      endif
    enddo

    call sum_reduce(ux_modes, 0)
    call sum_reduce(uy_modes, 0)
    call sum_reduce(uz_modes, 0)
    call sum_reduce(bx_modes, 0)
    call sum_reduce(by_modes, 0)
    call sum_reduce(bz_modes, 0)

    do n = 1, n_series_modes
      if(proc0) then
        i = series_modes(n, 1)
        j = series_modes(n, 2)
        k = series_modes(n, 3)
999 format(es30.21, A6, 5es30.21e3)
        write (unit=series_modes_unit, fmt=999) tt, 'ux', kx(i), ky(j), kz(k), ux_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'uy', kx(i), ky(j), kz(k), uy_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'uz', kx(i), ky(j), kz(k), uz_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bx', kx(i), ky(j), kz(k), bx_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'by', kx(i), ky(j), kz(k), by_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bz', kx(i), ky(j), kz(k), bz_modes(n)
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
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp
    use shearing_box, only: tsc, nremap, k2t, k2t_inv
    use shearing_box, only: timer_remap, remap_fld
    use dealias, only: filter
    implicit none
    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_remap)

    if(proc0) then
      print *
      print *, 'remapping...'
      print *
    endif

    call remap_fld(ux)
    call remap_fld(uy)
    call remap_fld(uz)
    call remap_fld(bx)
    call remap_fld(by)
    call remap_fld(bz)

    tsc = 0.d0
    nremap = nremap + 1
    counter = 1

    do k = ikz_st, ikz_en
      do j = iky_st, iky_en
        do i = ikx_st, ikx_en
          ux(i,k,j) = ux(i,k,j)*filter(i,k,j)
          uy(i,k,j) = uy(i,k,j)*filter(i,k,j)
          uz(i,k,j) = uz(i,k,j)*filter(i,k,j)
          bx(i,k,j) = bx(i,k,j)*filter(i,k,j)
          by(i,k,j) = by(i,k,j)*filter(i,k,j)
          bz(i,k,j) = bz(i,k,j)*filter(i,k,j)
        enddo
      enddo
    enddo

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          kxt(i,j) = kx(i)
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
!! @brief   Enforce div free for (fwx, fwy, fwz)
!-----------------------------------------------!
  subroutine div_free_force(fwx, fwy, fwz)
    use grid, only: ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: zi
    use shearing_box, only: kxt
    use mp, only: max_allreduce
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: fwx, fwy, fwz
    complex(r8), allocatable, dimension(:,:,:)   :: nbl2inv_div_fw ! nabla^-2 (div fw)
    integer :: i, j, k
    real(r8) :: fwx_max, fwy_max, fwz_max
    real(r8) :: eps, sx, sy, sz, k2t, k2t_inv

    allocate(nbl2inv_div_fw(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    fwx_max = maxval(abs(fwx)); call max_allreduce(fwx_max)
    fwy_max = maxval(abs(fwy)); call max_allreduce(fwy_max)
    fwz_max = maxval(abs(fwz)); call max_allreduce(fwz_max)

    eps = 1d-10
    if(fwx_max < eps) then
      sx = 0.d0
    else
      sx = 1.d0
    endif

    if(fwy_max < eps) then
      sy = 0.d0
    else
      sy = 1.d0
    endif

    if(fwz_max < eps) then
      sz = 0.d0
    else
      sz = 1.d0
    endif

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          k2t = sx*kxt(i, j)**2 + sy*ky(j)**2 + sz*kz(k)**2
          if(k2t == 0.d0) then
            k2t_inv = 0.d0
          else
            k2t_inv = 1.d0/k2t
          endif
          nbl2inv_div_fw(i,k,j) = -zi*(  kxt(i,j)*fwx(i,k,j) &
                                       + ky(j)   *fwy(i,k,j) &
                                       + kz(k)   *fwz(i,k,j) )*k2t_inv

          fwx(i,k,j) = fwx(i,k,j) - sx*zi*kxt(i,j)*nbl2inv_div_fw(i,k,j)
          fwy(i,k,j) = fwy(i,k,j) - sy*zi*ky(j)   *nbl2inv_div_fw(i,k,j)
          fwz(i,k,j) = fwz(i,k,j) - sz*zi*kz(k)   *nbl2inv_div_fw(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    deallocate(nbl2inv_div_fw)

  end subroutine div_free_force


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Allocates arrays for this module
!-----------------------------------------------!
  subroutine allocate_advance
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: time_step_scheme
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src

    if (is_allocated) then
      return
    else
      is_allocated = .true.
    endif

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))

    allocate(ux_new, source=src)
    allocate(uy_new, source=src)
    allocate(uz_new, source=src)
    allocate(bx_new, source=src)
    allocate(by_new, source=src)
    allocate(bz_new, source=src)

    allocate(w        (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))
    allocate(flx      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran ), source=(0.d0, 0.d0))
    allocate(w_r      (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nfields), source=0.d0)
    allocate(flx_r    (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nftran ), source=0.d0)
    allocate(exp_terms(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(ux_tmp, source=src)
      allocate(uy_tmp, source=src)
      allocate(uz_tmp, source=src)
      allocate(bx_tmp, source=src)
      allocate(by_tmp, source=src)
      allocate(bz_tmp, source=src)

      allocate(exp_terms0(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms0  = 0.d0
    endif

    deallocate(src)

  end subroutine allocate_advance


!-----------------------------------------------!
!> @author  YK
!! @date    18 May 2022
!! @brief   Deallocates arrays for this module
!-----------------------------------------------!
  subroutine deallocate_advance
    use params, only: time_step_scheme
    implicit none

    if (is_allocated) then
      is_allocated = .false.
    else
      return
    endif

    deallocate(ux_new)
    deallocate(uy_new)
    deallocate(uz_new)
    deallocate(bx_new)
    deallocate(by_new)
    deallocate(bz_new)

    deallocate(w        )
    deallocate(flx      )
    deallocate(w_r      )
    deallocate(flx_r    )
    deallocate(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      deallocate(ux_tmp)
      deallocate(uy_tmp)
      deallocate(uz_tmp)
      deallocate(bx_tmp)
      deallocate(by_tmp)
      deallocate(bz_tmp)

      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance


end module advance
