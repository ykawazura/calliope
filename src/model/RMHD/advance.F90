!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Time stepping for RMHD
!-----------------------------------------------!
module advance
  use p3dfft
  use fields, only: nfields
  use fields, only: iomg, ipsi
  implicit none

  public solve

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: phi_new
  complex(r8), allocatable, dimension(:,:,:)   :: omg_new
  complex(r8), allocatable, dimension(:,:,:)   :: psi_new
  complex(r8), allocatable, dimension(:,:,:,:) :: nonlin, exp_terms

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: phi_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: omg_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: psi_tmp
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: phi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: omg_old2
  complex(r8), allocatable, dimension(:,:,:)   :: psi_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fomg_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fpsi_old2

  real   (r8) :: cflx, cfly
  integer :: max_vel_unit

  ! Backward FFT variables
  integer, parameter :: nbtran = 8
  integer, parameter :: idphi_dx = 1, idphi_dy = 2
  integer, parameter :: idomg_dx = 3, idomg_dy = 4
  integer, parameter :: idpsi_dx = 5, idpsi_dy = 6
  integer, parameter :: idjpa_dx = 7, idjpa_dy = 8

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: phi, omg, psi
    use fields, only: phi_old, omg_old, psi_old
    use grid, only: kprp2, kprp2inv, kz2, kprp2_max, kz2_max
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, cfl, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: dealias_scheme => dealias
    use dealias, only: filter
    use params, only: nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      zi, nonlinear
    use force, only: fomg, fpsi, fomg_old, fpsi_old, driven, update_force, get_force
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use advance_common, only: gear1, gear2, gear3
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
      counter = 1
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      ! Calcualte force terms
      if (driven) then
        ! at n
        call get_force('omg', fomg)
        call get_force('psi', fpsi)
      endif

      !---------------  RK 1st step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, .true.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi(i,k,j), psi(i,k,j), &
                               nonlin(i,k,j,:), &
                               fomg(i,k,j), fpsi(i,k,j), &
                               kz(k), kprp2(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iomg), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ipsi), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)

            call eSSPIFRK1(omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), &
               imp_terms_tintg1(iomg), imp_terms_tintg0(iomg) &
            )

            call eSSPIFRK1(psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), &
               imp_terms_tintg1(ipsi), imp_terms_tintg0(ipsi) &
            )

            phi_tmp(i,k,j) = omg_tmp(i,k,j)*(-kprp2inv(i,k,j))
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
        call get_force('omg', fomg)
        call get_force('psi', fpsi)

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, .false.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi_tmp(i,k,j), psi_tmp(i,k,j), &
                               nonlin(i,k,j,:), &
                               fomg(i,k,j), fpsi(i,k,j), &
                               kz(k), kprp2(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iomg), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ipsi), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)

            call eSSPIFRK2(omg_tmp(i,k,j), omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), &
               imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )

            call eSSPIFRK2(psi_tmp(i,k,j), psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), &
               imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )

            phi_tmp(i,k,j) = omg_tmp(i,k,j)*(-kprp2inv(i,k,j))
          enddo
        enddo
      enddo
      !$omp end parallel do

      !---------------  RK 3rd step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, .false.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi_tmp(i,k,j), psi_tmp(i,k,j), &
                               nonlin(i,k,j,:), &
                               fomg(i,k,j), fpsi(i,k,j), &
                               kz(k), kprp2(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iomg), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iomg),           dt, kprp2(i,k,j), kz2(k), &
                                     nupe_x , nupe_z , nupe_x_exp , nupe_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ipsi), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ipsi),           dt, kprp2(i,k,j), kz2(k), &
                                     etape_x, etape_z, etape_x_exp, etape_z_exp)

            call eSSPIFRK3(omg_new(i,k,j), omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms (i,k,j,iomg), exp_terms0(i,k,j,iomg), &
               imp_terms_tintg3(iomg), imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )
            call eSSPIFRK3(psi_new(i,k,j), psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms (i,k,j,ipsi), exp_terms0(i,k,j,ipsi), &
               imp_terms_tintg3(ipsi), imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )

            phi_new(i,k,j) = omg_new(i,k,j)*(-kprp2inv(i,k,j))
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! save fields at the previous step
      !$omp workshare
      phi_old = phi
      omg_old = omg
      psi_old = psi
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fomg_old = fomg
        fpsi_old = fpsi
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
        call get_force('omg', fomg)
        call get_force('psi', fpsi)
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, .true.)

      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi(i,k,j), psi(i,k,j), &
                               nonlin(i,k,j,:), &
                               fomg(i,k,j), fpsi(i,k,j), &
                               kz(k), kprp2(i,k,j))

            ! 1st order 
            if(counter == 1) then
              call gear1(omg_new(i,k,j), omg(i,k,j), &
                 exp_terms(i,k,j,iomg), &
                 nupe_x*(kprp2(i,k,j)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
              )
              call gear1(psi_new(i,k,j), psi(i,k,j), &
                 exp_terms(i,k,j,ipsi), &
                 etape_x*(kprp2(i,k,j)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
              )

            ! 2nd order 
            elseif(counter == 2) then
              call gear2(omg_new(i,k,j), omg(i,k,j), omg_old(i,k,j), &
                 exp_terms    (i,k,j,iomg), &
                 exp_terms_old(i,k,j,iomg), &
                 nupe_x*(kprp2(i,k,j)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
              )
              call gear2(psi_new(i,k,j), psi(i,k,j), psi_old(i,k,j), &
                 exp_terms    (i,k,j,ipsi), &
                 exp_terms_old(i,k,j,ipsi), &
                 etape_x*(kprp2(i,k,j)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
              )

            ! 3rd order 
            else
              call gear3(omg_new(i,k,j), omg(i,k,j), omg_old(i,k,j), omg_old2(i,k,j), &
                 exp_terms     (i,k,j,iomg), &
                 exp_terms_old (i,k,j,iomg), &
                 exp_terms_old2(i,k,j,iomg), &
                 nupe_x*(kprp2(i,k,j)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
              )
              call gear3(psi_new(i,k,j), psi(i,k,j), psi_old(i,k,j), psi_old2(i,k,j), &
                 exp_terms     (i,k,j,ipsi), &
                 exp_terms_old (i,k,j,ipsi), &
                 exp_terms_old2(i,k,j,ipsi), &
                 etape_x*(kprp2(i,k,j)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
              )

            endif
            phi_new(i,k,j) = omg_new(i,k,j)*(-kprp2inv(i,k,j))
          enddo
        enddo
      enddo
      !$omp end parallel do

      if(counter <= 2) counter = counter + 1

      ! save fields at the previous steps
      !$omp workshare
      phi_old2 = phi_old
      phi_old  = phi

      omg_old2 = omg_old 
      omg_old  = omg

      psi_old2 = psi_old 
      psi_old  = psi

      exp_terms_old2 = exp_terms_old
      exp_terms_old  = exp_terms
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fomg_old2 = fomg_old 
        fomg_old  = fomg

        fpsi_old2 = fpsi_old 
        fpsi_old  = fpsi
        !$omp end workshare
      endif
    endif

    !$omp workshare
    phi = phi_new
    omg = omg_new
    psi = psi_new
    !$omp end workshare

    ! Dealiasing
    if(trim(dealias_scheme) /= '2/3') then
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi(i,k,j) = phi_new(i,k,j)*filter(i,k,j)
            psi(i,k,j) = psi_new(i,k,j)*filter(i,k,j)
            omg(i,k,j) = omg_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
    endif

    tt = tt  + dt

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Allocate tmp fields for a multi 
!!          timestep method
!-----------------------------------------------!
  subroutine init_multistep_fields
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use file, only: open_output_file
    use params, only: time_step_scheme
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))

    allocate(phi_new  , source=src)
    allocate(omg_new  , source=src)
    allocate(psi_new  , source=src)

    allocate(nonlin   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); nonlin    = 0.d0
    allocate(exp_terms(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms = 0.d0

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(phi_tmp, source=src)
      allocate(omg_tmp, source=src)
      allocate(psi_tmp, source=src)

      allocate(exp_terms0(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms0 = 0.d0
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      allocate(phi_old2 , source=src)
      allocate(omg_old2 , source=src)
      allocate(psi_old2 , source=src)

      allocate(fomg_old2, source=src)
      allocate(fpsi_old2, source=src)

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
  subroutine get_nonlinear_terms(phi, psi, dt_reset)
    use grid, only: kprp2, nlx, nly, nlz, kx, ky
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: phi, psi

    complex(r8), allocatable, dimension(:,:,:,:) :: wbk
    real   (r8), allocatable, dimension(:,:,:,:) :: wb , wf 
    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (r8) :: max_vel, dt_cfl, dt_digit

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    allocate(wbk(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nbtran ), source=(0.d0, 0.d0))
    allocate(wb (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nbtran ), source=0.d0)
    allocate(wf (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nfields), source=0.d0)

    ! 1. Calculate grad in Fourier space
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          wbk(i,k,j,idphi_dx) = zi*kx(i)                  *phi(i, k, j)
          wbk(i,k,j,idphi_dy) = zi*ky(j)                  *phi(i, k, j)
          wbk(i,k,j,idomg_dx) = zi*kx(i)*(-kprp2(i, k, j))*phi(i, k, j)
          wbk(i,k,j,idomg_dy) = zi*ky(j)*(-kprp2(i, k, j))*phi(i, k, j)
                                 
          wbk(i,k,j,idpsi_dx) = zi*kx(i)                  *psi(i, k, j)
          wbk(i,k,j,idpsi_dy) = zi*ky(j)                  *psi(i, k, j)
          wbk(i,k,j,idjpa_dx) = zi*kx(i)*(-kprp2(i, k, j))*psi(i, k, j)
          wbk(i,k,j,idjpa_dy) = zi*ky(j)*(-kprp2(i, k, j))*psi(i, k, j)
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
                maxval(abs(wb(:,:,:,idphi_dx)))*cfly, maxval(abs(wb(:,:,:,idphi_dy))*cflx) &
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
          wf(j,k,i,iomg) = - wb(j,k,i,idphi_dx)*wb(j,k,i,idomg_dy) &
                           + wb(j,k,i,idphi_dy)*wb(j,k,i,idomg_dx) &
                           + wb(j,k,i,idpsi_dx)*wb(j,k,i,idjpa_dy) &
                           - wb(j,k,i,idpsi_dy)*wb(j,k,i,idjpa_dx)
          wf(j,k,i,ipsi) = - wb(j,k,i,idphi_dx)*wb(j,k,i,idpsi_dy) &
                           + wb(j,k,i,idphi_dy)*wb(j,k,i,idpsi_dx)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 4. Forward FFT
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_ftran_r2c_many(wf, nl_local_tot, nonlin, nk_local_tot, nfields, 'fft')
    if (proc0) call put_time_stamp(timer_fft)
    !$omp workshare
    nonlin = nonlin/nlx/nly/nlz
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
                           phi, psi, &
                           nonlin, &
                           fomg, fpsi, &
                           kz, kprp2)
    use params, only: zi
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: phi, psi
    complex(r8), intent(in ) :: nonlin(nfields)
    complex(r8), intent(in ) :: fomg, fpsi
    real(r8)   , intent(in)  :: kz, kprp2

    exp_terms(iomg) = nonlin(iomg) + fomg - zi*kz*kprp2*psi
    exp_terms(ipsi) = nonlin(ipsi) + fpsi + zi*kz*phi

  end subroutine get_ext_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Time integral of hyperdissipation
!-----------------------------------------------!
  subroutine get_imp_terms_tintg(imp_terms_tintg, t, kprp2, kz2, coeff_x, coeff_z, nexp_x, nexp_z)
    use grid, only: kprp2_max, kz2_max
    implicit none
    real(r8), intent(out) :: imp_terms_tintg
    real(r8), intent(in) :: t, kprp2, kz2, coeff_x, coeff_z
    integer, intent(in) :: nexp_x, nexp_z

    imp_terms_tintg = -(coeff_x*(kprp2/kprp2_max)**nexp_x + coeff_z*(kz2/kz2_max)**nexp_z)*t
  end subroutine get_imp_terms_tintg

end module advance

