!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    25 Sep 2019
!! @brief   Time stepping for MHD_COMP_ISOTH
!-----------------------------------------------!
module advance
  use p3dfft
  use fields, only: nfields
  use fields, only: irho
  use fields, only: imx, imy, imz
  use fields, only: ibx, iby, ibz
  implicit none

  public solve

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: rho_new
  complex(r8), allocatable, dimension(:,:,:)   ::  mx_new
  complex(r8), allocatable, dimension(:,:,:)   ::  my_new
  complex(r8), allocatable, dimension(:,:,:)   ::  mz_new
  complex(r8), allocatable, dimension(:,:,:)   ::  bx_new
  complex(r8), allocatable, dimension(:,:,:)   ::  by_new
  complex(r8), allocatable, dimension(:,:,:)   ::  bz_new
  complex(r8), allocatable, dimension(:,:,:,:) :: w
  real   (r8), allocatable, dimension(:,:,:,:) :: w_r, flx_r 
  real   (r8), allocatable, dimension(:,:,:)   :: ux_r, uy_r, uz_r
  real   (r8), allocatable, dimension(:,:,:)   :: cfx2_r, cfy2_r, cfz2_r ! Fast mode phase speed
  complex(r8), allocatable, dimension(:,:,:,:) :: flx, exp_terms
  real   (r8), allocatable, dimension(:,:)     :: kxt

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK4                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: rho_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::  mx_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::  my_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::  mz_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::  bx_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::  by_tmp
  complex(r8), allocatable, dimension(:,:,:)   ::  bz_tmp
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: rho_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  mx_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  my_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  mz_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  bx_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  by_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  bz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fmx_old2, fmy_old2, fmz_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2
  real   (r8), allocatable, dimension(:,:)     :: kxt_old, kxt_old2

  integer :: cfl_unit

  ! Forward FFT variables
  integer, parameter :: nftran = 9
  integer, parameter :: iflx_mxx = 1, iflx_mxy = 2, iflx_mxz = 3 ! rho*uu - bb + [(cs^2/va^2)*rho + b^2/2]*I (tensor) 
  integer, parameter ::               iflx_myy = 4, iflx_myz = 5 ! rho*uu - bb + [(cs^2/va^2)*rho + b^2/2]*I (tensor) 
  integer, parameter ::                             iflx_mzz = 6 ! rho*uu - bb + [(cs^2/va^2)*rho + b^2/2]*I (tensor) 
  integer, parameter :: iflx_bx  = 7, iflx_by  = 8, iflx_bz  = 9 ! b x u

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: rho_old
    use fields, only: mx_old, my_old, mz_old
    use fields, only: bx_old, by_old, bz_old
    use grid, only: k2_max
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nu, nu_h, nu_h_exp, eta, eta_h, eta_h_exp, lmd, lmd_h, lmd_h_exp, zi, nonlinear, shear, q
    use params, only: rho_min
    use dealias, only: filter
    use shearing_box, only: shear_flg, k2t, k2t_inv, tsc, tremap
    use force, only: fmx, fmz, fmy, fmx_old, fmy_old, fmz_old, driven, update_force, get_force
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use advance_common, only: gear1, gear2, gear3
    use utils, only: check_floor
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
        call get_force('mx', fmx)
        call get_force('my', fmy)
        call get_force('mz', fmz)
      endif

      !---------------  RK 1st step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(rho, mx, my, mz, bx, by, bz, .true.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg1)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               mx(i,k,j), my(i,k,j), mz(i,k,j), bx(i,k,j), by(i,k,j), bz(i,k,j), &
                               flx(i,k,j,:), &
                               fmx(i,k,j), fmy(i,k,j), fmz(i,k,j), &
                               kxt(i,j), ky(j), kz(k))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(irho), tsc               , kx(i), ky(j), kz(k), lmd, lmd_h, lmd_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(irho), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), lmd, lmd_h, lmd_h_exp)

            call get_imp_terms_tintg(imp_terms_tintg0(imx ), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(imx ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(imy ), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(imy ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(imz ), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(imz ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
                                                          
            call get_imp_terms_tintg(imp_terms_tintg0(ibx ), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ibx ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby ), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(iby ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz ), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ibz ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update rho
            call eSSPIFRK1(rho_tmp(i,k,j), rho(i,k,j), &
               exp_terms(i,k,j,irho), &
               imp_terms_tintg1(irho), imp_terms_tintg0(irho) &
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
      !$omp parallel do private(j, k, i)
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            rho_tmp(i,k,j) = rho_tmp(i,k,j)*filter(i,k,j)
            mx_tmp (i,k,j) = mx_tmp (i,k,j)*filter(i,k,j)
            my_tmp (i,k,j) = my_tmp (i,k,j)*filter(i,k,j)
            mz_tmp (i,k,j) = mz_tmp (i,k,j)*filter(i,k,j)
            bx_tmp (i,k,j) = bx_tmp (i,k,j)*filter(i,k,j)
            by_tmp (i,k,j) = by_tmp (i,k,j)*filter(i,k,j)
            bz_tmp (i,k,j) = bz_tmp (i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
      !$omp end parallel do

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
      if(shear) then
        !$omp parallel do private(i, k) schedule(static)
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            kxt(i,j) = kx(i) + q*shear_flg*(tsc + 2.d0/3.d0*dt)*ky(j)
          enddo
        enddo
        !$omp end parallel do
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(rho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               mx_tmp(i,k,j), my_tmp(i,k,j), mz_tmp(i,k,j), bx_tmp(i,k,j), by_tmp(i,k,j), bz_tmp(i,k,j), &
                               flx(i,k,j,:), &
                               fmx(i,k,j), fmy(i,k,j), fmz(i,k,j), &
                               kxt(i,j), ky(j), kz(k))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(irho), tsc               , kx(i), ky(j), kz(k), lmd, lmd_h, lmd_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(irho), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), lmd, lmd_h, lmd_h_exp)
                                                                                                         
            call get_imp_terms_tintg(imp_terms_tintg0(imx ), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(imx ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(imy ), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(imy ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(imz ), tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(imz ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
                                                                                                         
            call get_imp_terms_tintg(imp_terms_tintg0(ibx ), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby ), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz ), tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz ), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)


            ! update rho
            call eSSPIFRK2(rho_tmp(i,k,j), rho_tmp(i,k,j), rho(i,k,j), &
               exp_terms(i,k,j,irho), &
               imp_terms_tintg2(irho), imp_terms_tintg0(irho) &
            )

            ! update u
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
      !$omp parallel do private(j, k, i)
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            rho_tmp(i,k,j) = rho_tmp(i,k,j)*filter(i,k,j)
            mx_tmp (i,k,j) = mx_tmp (i,k,j)*filter(i,k,j)
            my_tmp (i,k,j) = my_tmp (i,k,j)*filter(i,k,j)
            mz_tmp (i,k,j) = mz_tmp (i,k,j)*filter(i,k,j)
            bx_tmp (i,k,j) = bx_tmp (i,k,j)*filter(i,k,j)
            by_tmp (i,k,j) = by_tmp (i,k,j)*filter(i,k,j)
            bz_tmp (i,k,j) = bz_tmp (i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
      !$omp end parallel do

      !---------------  RK 3rd step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(rho_tmp, mx_tmp, my_tmp, mz_tmp, bx_tmp, by_tmp, bz_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               mx_tmp(i,k,j), my_tmp(i,k,j), mz_tmp(i,k,j), bx_tmp(i,k,j), by_tmp(i,k,j), bz_tmp(i,k,j), &
                               flx(i,k,j,:), &
                               fmx(i,k,j), fmy(i,k,j), fmz(i,k,j), &
                               kxt(i,j), ky(j), kz(k))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(irho), tsc               , kx(i), ky(j), kz(k), lmd, lmd_h, lmd_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(irho), tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), lmd, lmd_h, lmd_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(irho), tsc +           dt, kx(i), ky(j), kz(k), lmd, lmd_h, lmd_h_exp)

            call get_imp_terms_tintg(imp_terms_tintg0(imx) , tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(imx) , tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(imx) , tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(imy) , tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(imy) , tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(imy) , tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(imz) , tsc               , kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(imz) , tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(imz) , tsc +           dt, kx(i), ky(j), kz(k), nu , nu_h , nu_h_exp )
                                                           
            call get_imp_terms_tintg(imp_terms_tintg0(ibx) , tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibx) , tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibx) , tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(iby) , tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(iby) , tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(iby) , tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg0(ibz) , tsc               , kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibz) , tsc + 2.d0/3.d0*dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibz) , tsc +           dt, kx(i), ky(j), kz(k), eta, eta_h, eta_h_exp)

            ! update rho
            call eSSPIFRK3(rho_new(i,k,j), rho_tmp(i,k,j), rho(i,k,j), &
               exp_terms(i,k,j,irho), exp_terms0(i,k,j,irho), &
               imp_terms_tintg3(irho), imp_terms_tintg2(irho), imp_terms_tintg0(irho) &
            )

            ! update u
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
      rho_old = rho
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
      !$omp parallel do private(j, k, i)
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            rho_new(i,k,j) = rho_new(i,k,j)*filter(i,k,j)
            mx_new (i,k,j) = mx_new (i,k,j)*filter(i,k,j)
            my_new (i,k,j) = my_new (i,k,j)*filter(i,k,j)
            mz_new (i,k,j) = mz_new (i,k,j)*filter(i,k,j)
            bx_new (i,k,j) = bx_new (i,k,j)*filter(i,k,j)
            by_new (i,k,j) = by_new (i,k,j)*filter(i,k,j)
            bz_new (i,k,j) = bz_new (i,k,j)*filter(i,k,j)
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
        call get_force('mx', fmx)
        call get_force('my', fmy)
        call get_force('mz', fmz)
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(rho, mx, my, mz, bx, by, bz, .true.)

      !$omp parallel do private(j, k, i)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               mx(i,k,j), my(i,k,j), mz(i,k,j), bx(i,k,j), by(i,k,j), bz(i,k,j), &
                               flx(i,k,j,:), &
                               fmx(i,k,j), fmy(i,k,j), fmz(i,k,j), &
                               kxt(i,j), ky(j), kz(k))

            ! 1st order 
            if(counter == 1) then
              ! update rho
              call gear1(rho_new(i,k,j), rho(i,k,j), &
                 exp_terms(i,k,j,irho), &
                 lmd*(k2t(i,k,j)/k2_max) + lmd_h*(k2t(i,k,j)/k2_max)**lmd_h_exp &
              )

              ! update m
              call gear1(mx_new(i,k,j), mx(i,k,j), &
                 exp_terms(i,k,j,imx), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear1(my_new(i,k,j), my(i,k,j), &
                 exp_terms(i,k,j,imy), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear1(mz_new(i,k,j), mz(i,k,j), &
                 exp_terms(i,k,j,imz), &
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
              ! update rho
              call gear2(rho_new(i,k,j), rho(i,k,j), rho_old(i,k,j), &
                 exp_terms    (i,k,j,irho), &
                 exp_terms_old(i,k,j,irho), &
                 lmd*(k2t(i,k,j)/k2_max) + lmd_h*(k2t(i,k,j)/k2_max)**lmd_h_exp &
              )

              ! update m
              call gear2(mx_new(i,k,j), mx(i,k,j), mx_old(i,k,j), &
                 exp_terms    (i,k,j,imx), &
                 exp_terms_old(i,k,j,imx), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear2(my_new(i,k,j), my(i,k,j), my_old(i,k,j), &
                 exp_terms    (i,k,j,imy), &
                 exp_terms_old(i,k,j,imy), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear2(mz_new(i,k,j), mz(i,k,j), mz_old(i,k,j), &
                 exp_terms    (i,k,j,imz), &
                 exp_terms_old(i,k,j,imz), &
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
              ! update rho
              call gear3(rho_new(i,k,j), rho(i,k,j), rho_old(i,k,j), rho_old2(i,k,j), &
                 exp_terms     (i,k,j,irho), &
                 exp_terms_old (i,k,j,irho), &
                 exp_terms_old2(i,k,j,irho), &
                 lmd*(k2t(i,k,j)/k2_max) + lmd_h*(k2t(i,k,j)/k2_max)**lmd_h_exp &
              )

              ! update m
              call gear3(mx_new(i,k,j), mx(i,k,j), mx_old(i,k,j), mx_old2(i,k,j), &
                 exp_terms     (i,k,j,imx), &
                 exp_terms_old (i,k,j,imx), &
                 exp_terms_old2(i,k,j,imx), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear3(my_new(i,k,j), my(i,k,j), my_old(i,k,j), my_old2(i,k,j), &
                 exp_terms     (i,k,j,imy), &
                 exp_terms_old (i,k,j,imy), &
                 exp_terms_old2(i,k,j,imy), &
                 nu*(k2t(i,k,j)/k2_max) + nu_h*(k2t(i,k,j)/k2_max)**nu_h_exp &
              )
              call gear3(mz_new(i,k,j), mz(i,k,j), mz_old(i,k,j), mz_old2(i,k,j), &
                 exp_terms     (i,k,j,imz), &
                 exp_terms_old (i,k,j,imz), &
                 exp_terms_old2(i,k,j,imz), &
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
      rho_old2 = rho_old 
      rho_old  = rho

      mx_old2 = mx_old 
      mx_old  = mx

      my_old2 = my_old 
      my_old  = my

      mz_old2 = mz_old 
      mz_old  = mz

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
        fmx_old2 = fmx_old 
        fmx_old  = fmx

        fmy_old2 = fmy_old 
        fmy_old  = fmy

        fmz_old2 = fmz_old 
        fmz_old  = fmz
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
            rho_new(i,k,j) = rho_new(i,k,j)*filter(i,k,j)
            mx_new (i,k,j) = mx_new (i,k,j)*filter(i,k,j)
            my_new (i,k,j) = my_new (i,k,j)*filter(i,k,j)
            mz_new (i,k,j) = mz_new (i,k,j)*filter(i,k,j)
            bx_new (i,k,j) = bx_new (i,k,j)*filter(i,k,j)
            by_new (i,k,j) = by_new (i,k,j)*filter(i,k,j)
            bz_new (i,k,j) = bz_new (i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
      !$omp end parallel do
    endif

    !$omp workshare
    rho = rho_new
     mx =  mx_new
     my =  my_new
     mz =  mz_new
     bx =  bx_new
     by =  by_new
     bz =  bz_new
    !$omp end workshare

    call check_floor(rho, rho_min)

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

    ! Div b cleaing
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
      allocate(rho_old2, source=src)
      allocate( mx_old2, source=src)
      allocate( my_old2, source=src)
      allocate( mz_old2, source=src)
      allocate( bx_old2, source=src)
      allocate( by_old2, source=src)
      allocate( bz_old2, source=src)

      allocate(fmx_old2, source=src)
      allocate(fmy_old2, source=src)
      allocate(fmz_old2, source=src)

      allocate(exp_terms_old (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old2 = 0.d0

      allocate(kxt_old (ikx_st:ikx_en, iky_st:iky_en)); kxt_old  = kxt
      allocate(kxt_old2(ikx_st:ikx_en, iky_st:iky_en)); kxt_old2 = kxt
    endif

    if(proc0) call open_output_file (cfl_unit, 'cfl.dat')

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
  subroutine get_nonlinear_terms(rho, mx, my, mz, bx, by, bz, dt_reset)
    use grid, only: dlx, dly, dlz
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi, cs2va2
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    use advance_common, only: dt_adjust_while_running 
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: rho, mx, my, mz, bx, by, bz

    real   (r8) :: vax2_r, vay2_r, vaz2_r ! Alfven speed
    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (r8) :: max_vel_x, max_vel_y, max_vel_z , dt_cfl, dt_digit

    if (proc0) call put_time_stamp(timer_nonlinear_terms)


    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          w(i,k,j,irho) = rho(i,k,j)
          w(i,k,j,imx ) = mx (i,k,j)
          w(i,k,j,imy ) = my (i,k,j)
          w(i,k,j,imz ) = mz (i,k,j)
          w(i,k,j,ibx ) = bx (i,k,j)
          w(i,k,j,iby ) = by (i,k,j)
          w(i,k,j,ibz ) = bz (i,k,j)
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

    !$omp workshare
    ux_r = w_r(:,:,:,imx)/w_r(:,:,:,irho)
    uy_r = w_r(:,:,:,imy)/w_r(:,:,:,irho)
    uz_r = w_r(:,:,:,imz)/w_r(:,:,:,irho)
    !$omp end workshare

    ! (get max_vel for dt reset)
    if(dt_reset) then
      !$omp parallel do private(j, k) schedule(static)
      do i = ilx_st, ilx_en
        do k = ilz_st, ilz_en
          do j = ily_st, ily_en
            vax2_r = w_r(j,k,i,ibx)**2/w_r(j,k,i,irho)
            vay2_r = w_r(j,k,i,iby)**2/w_r(j,k,i,irho)
            vaz2_r = w_r(j,k,i,ibz)**2/w_r(j,k,i,irho)

            cfx2_r(j,k,i) = 0.5d0*(cs2va2 + vax2_r +  vay2_r + vaz2_r  &
                            + sqrt( (cs2va2 + vax2_r +  vay2_r + vaz2_r)**2 - 4.d0*cs2va2*vax2_r))
            cfy2_r(j,k,i) = 0.5d0*(cs2va2 + vax2_r +  vay2_r + vaz2_r  &
                            + sqrt( (cs2va2 + vax2_r +  vay2_r + vaz2_r)**2 - 4.d0*cs2va2*vay2_r))
            cfz2_r(j,k,i) = 0.5d0*(cs2va2 + vax2_r +  vay2_r + vaz2_r  &
                            + sqrt( (cs2va2 + vax2_r +  vay2_r + vaz2_r)**2 - 4.d0*cs2va2*vaz2_r))
          enddo
        enddo
      enddo
      !$omp end parallel do
      max_vel_x = max( &
                    maxval(abs(ux_r + sqrt(cfx2_r))), &
                    maxval(abs(ux_r - sqrt(cfx2_r)))  &
                  )
      max_vel_y = max( &
                    maxval(abs(uy_r + sqrt(cfy2_r))), &
                    maxval(abs(uy_r - sqrt(cfy2_r)))  &
                  )
      max_vel_z = max( &
                    maxval(abs(uz_r + sqrt(cfz2_r))), &
                    maxval(abs(uz_r - sqrt(cfz2_r)))  &
                  )
      call max_allreduce(max_vel_x)
      call max_allreduce(max_vel_y)
      call max_allreduce(max_vel_z)
      dt_cfl = cfl*min(dlx/max_vel_x, dly/max_vel_y, dlz/max_vel_z)

      if(proc0) then
        write (unit=cfl_unit, fmt="(100es30.21)") tt, dt_cfl, max_vel_x, max_vel_y, max_vel_z
        call flush(cfl_unit) 
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

    ! When the file 'dt_adjust' including a float number is created,
    ! dt will be manually adjusted to that value while running.
    call dt_adjust_while_running() 

    ! 2. Calculate nonlinear terms in real space
    !$omp parallel do private(j, k) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          flx_r(j,k,i,iflx_mxx) = w_r(j,k,i,irho)*ux_r(j,k,i)*ux_r(j,k,i) - w_r(j,k,i,ibx)*w_r(j,k,i,ibx) &
                                 + 0.5d0*(w_r(j,k,i,ibx)**2 + w_r(j,k,i,iby)**2 + w_r(j,k,i,ibz)**2) &
                                 + cs2va2*w_r(j,k,i,irho)
          flx_r(j,k,i,iflx_mxy) = w_r(j,k,i,irho)*ux_r(j,k,i)*uy_r(j,k,i) - w_r(j,k,i,ibx)*w_r(j,k,i,iby)
          flx_r(j,k,i,iflx_mxz) = w_r(j,k,i,irho)*ux_r(j,k,i)*uz_r(j,k,i) - w_r(j,k,i,ibx)*w_r(j,k,i,ibz)
          flx_r(j,k,i,iflx_myy) = w_r(j,k,i,irho)*uy_r(j,k,i)*uy_r(j,k,i) - w_r(j,k,i,iby)*w_r(j,k,i,iby) &
                                 + 0.5d0*(w_r(j,k,i,ibx)**2 + w_r(j,k,i,iby)**2 + w_r(j,k,i,ibz)**2) &
                                 + cs2va2*w_r(j,k,i,irho)
          flx_r(j,k,i,iflx_myz) = w_r(j,k,i,irho)*uy_r(j,k,i)*uz_r(j,k,i) - w_r(j,k,i,iby)*w_r(j,k,i,ibz)
          flx_r(j,k,i,iflx_mzz) = w_r(j,k,i,irho)*uz_r(j,k,i)*uz_r(j,k,i) - w_r(j,k,i,ibz)*w_r(j,k,i,ibz) &
                                 + 0.5d0*(w_r(j,k,i,ibx)**2 + w_r(j,k,i,iby)**2 + w_r(j,k,i,ibz)**2) &
                                 + cs2va2*w_r(j,k,i,irho)

          flx_r(j,k,i,iflx_bx ) = w_r(j,k,i,iby)*uz_r(j,k,i) - w_r(j,k,i,ibz)*uy_r(j,k,i)
          flx_r(j,k,i,iflx_by ) = w_r(j,k,i,ibz)*ux_r(j,k,i) - w_r(j,k,i,ibx)*uz_r(j,k,i)
          flx_r(j,k,i,iflx_bz ) = w_r(j,k,i,ibx)*uy_r(j,k,i) - w_r(j,k,i,iby)*ux_r(j,k,i)
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
    flx = flx/nlx/nly/nlz

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms, &
                           mx, my, mz, bx, by, bz, &
                           flx, &
                           fmx, fmy, fmz, &
                           kxt, ky, kz)
    use params, only: zi, q
    use shearing_box, only: shear_flg
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: mx, my, mz, bx, by, bz
    complex(r8), intent(in ) :: flx(nftran)
    complex(r8), intent(in ) :: fmx, fmy, fmz
    real(r8)   , intent(in)  :: kxt, ky, kz
    complex(r8)              :: nl(nfields)

    ! div m
    nl(irho) = -zi*( kxt*mx + ky*my + kz*mz )

    ! div {rho*uu - bb + [(cs^2/va^2)*rho + b^2/2]*I}
    nl(imx ) = -zi*( kxt*flx(iflx_mxx) + ky*flx(iflx_mxy) + kz*flx(iflx_mxz) )
    nl(imy ) = -zi*( kxt*flx(iflx_mxy) + ky*flx(iflx_myy) + kz*flx(iflx_myz) )
    nl(imz ) = -zi*( kxt*flx(iflx_mxz) + ky*flx(iflx_myz) + kz*flx(iflx_mzz) )

    ! curl (b x u)
    nl(ibx ) = -zi*( ky *flx(iflx_bz) - kz *flx(iflx_by) )
    nl(iby ) = -zi*( kz *flx(iflx_bx) - kxt*flx(iflx_bz) )
    nl(ibz ) = -zi*( kxt*flx(iflx_by) - ky *flx(iflx_bx) )

    exp_terms(irho) = nl(irho)

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
!! @date    15 Jul 2021
!! @brief   Output series modes
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, sum_reduce
    use time, only: tt
    use diagnostics_common, only: n_series_modes, series_modes
    use diagnostics_common, only: series_modes_unit
    implicit none
    complex(r8), dimension(n_series_modes) :: rho_modes
    complex(r8), dimension(n_series_modes) ::  mx_modes, my_modes, mz_modes
    complex(r8), dimension(n_series_modes) ::  bx_modes, by_modes, bz_modes
    integer :: n, i, j, k

    rho_modes(:) = 0.d0
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

        rho_modes(n) = rho(i, k, j)
         mx_modes(n) =  mx(i, k, j)
         my_modes(n) =  my(i, k, j)
         mz_modes(n) =  mz(i, k, j)

         bx_modes(n) =  bx(i, k, j)
         by_modes(n) =  by(i, k, j)
         bz_modes(n) =  bz(i, k, j)

      endif
    enddo

    call sum_reduce(rho_modes, 0)
    call sum_reduce( mx_modes, 0)
    call sum_reduce( my_modes, 0)
    call sum_reduce( mz_modes, 0)
    call sum_reduce( bx_modes, 0)
    call sum_reduce( by_modes, 0)
    call sum_reduce( bz_modes, 0)

    do n = 1, n_series_modes
      if(proc0) then
        i = series_modes(n, 1)
        j = series_modes(n, 2)
        k = series_modes(n, 3)
999 format(es30.21, A6, 5es30.21e3)
        write (unit=series_modes_unit, fmt=999) tt, 'rho', kx(i), ky(j), kz(k), rho_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,  'mx', kx(i), ky(j), kz(k),  mx_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,  'my', kx(i), ky(j), kz(k),  my_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,  'mz', kx(i), ky(j), kz(k),  mz_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,  'bx', kx(i), ky(j), kz(k),  bx_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,  'by', kx(i), ky(j), kz(k),  by_modes(n)
        write (unit=series_modes_unit, fmt=999) tt,  'bz', kx(i), ky(j), kz(k),  bz_modes(n)
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
    use fields, only: rho
    use fields, only: mx, my, mz
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

    call remap_fld(mx)
    call remap_fld(my)
    call remap_fld(mz)
    call remap_fld(bx)
    call remap_fld(by)
    call remap_fld(bz)
    call remap_fld(rho)

    tsc = 0.d0
    nremap = nremap + 1
    counter = 1

    ! Dealiasing
    do k = ikz_st, ikz_en
      do j = iky_st, iky_en
        do i = ikx_st, ikx_en
          rho(i,k,j) = rho(i,k,j)*filter(i,k,j)
          mx (i,k,j) = mx (i,k,j)*filter(i,k,j)
          my (i,k,j) = my (i,k,j)*filter(i,k,j)
          mz (i,k,j) = mz (i,k,j)*filter(i,k,j)
          bx (i,k,j) = bx (i,k,j)*filter(i,k,j)
          by (i,k,j) = by (i,k,j)*filter(i,k,j)
          bz (i,k,j) = bz (i,k,j)*filter(i,k,j)
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

    allocate(rho_new, source=src)
    allocate(mx_new , source=src)
    allocate(my_new , source=src)
    allocate(mz_new , source=src)
    allocate(bx_new , source=src)
    allocate(by_new , source=src)
    allocate(bz_new , source=src)

    allocate(w        (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))
    allocate(flx      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran ), source=(0.d0, 0.d0))
    allocate(w_r      (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nfields), source=0.d0)
    allocate(flx_r    (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nftran ), source=0.d0)
    allocate(exp_terms(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))

    allocate(ux_r  (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(uy_r  (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(uz_r  (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(cfx2_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(cfy2_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)
    allocate(cfz2_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(rho_tmp, source=src)
      allocate( mx_tmp, source=src)
      allocate( my_tmp, source=src)
      allocate( mz_tmp, source=src)
      allocate( bx_tmp, source=src)
      allocate( by_tmp, source=src)
      allocate( bz_tmp, source=src)

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

    deallocate(rho_new)
    deallocate( mx_new)
    deallocate( my_new)
    deallocate( mz_new)
    deallocate( bx_new)
    deallocate( by_new)
    deallocate( bz_new)

    deallocate(w        )
    deallocate(flx      )
    deallocate(w_r      )
    deallocate(flx_r    )
    deallocate(exp_terms)

    deallocate(ux_r  )
    deallocate(uy_r  )
    deallocate(uz_r  )
    deallocate(cfx2_r)
    deallocate(cfy2_r)
    deallocate(cfz2_r)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      deallocate(rho_tmp)
      deallocate( mx_tmp)
      deallocate( my_tmp)
      deallocate( mz_tmp)
      deallocate( bx_tmp)
      deallocate( by_tmp)
      deallocate( bz_tmp)

      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance

end module advance



