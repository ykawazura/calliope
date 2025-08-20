!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Time stepping for RRMHD
!-----------------------------------------------!
module advance
  use p3dfft
  use fields, only: nfields
  use fields, only: iomg, ipsi, iaaa
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: phi_new
  complex(r8), allocatable, dimension(:,:,:)   :: omg_new
  complex(r8), allocatable, dimension(:,:,:)   :: psi_new
  complex(r8), allocatable, dimension(:,:,:)   :: aaa_new
  complex(r8), allocatable, dimension(:,:,:,:) :: w
  real   (r8), allocatable, dimension(:,:,:,:) :: w_r, nonlin_r 
  complex(r8), allocatable, dimension(:,:,:,:) :: nonlin
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                For eSSPIFRK3                v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: phi_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: omg_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: psi_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: aaa_tmp
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: phi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: omg_old2
  complex(r8), allocatable, dimension(:,:,:)   :: psi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: aaa_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2

  integer :: cfl_unit

  ! Backward FFT variables
  integer, parameter :: nbtran = 12
  integer, parameter :: idphi_dx = 1 , idphi_dy = 2
  integer, parameter :: idomg_dx = 3 , idomg_dy = 4
  integer, parameter :: idpsi_dx = 5 , idpsi_dy = 6
  integer, parameter :: idjpa_dx = 7 , idjpa_dy = 8
  integer, parameter :: idaaa_dx = 9 , idaaa_dy = 10

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: phi, omg, psi, aaa
    use fields, only: phi_old, omg_old, psi_old, aaa_old
    use grid, only: kprp2, kprp2inv, kz2, kprp2_max, kz2_max
    use grid, only: ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use dealias, only: filter
    use params, only: zi, nonlinear,  nupe , nupe_exp , nuz , nuz_exp , &
                                      etape, etape_exp, etaz, etaz_exp, &
                                      chipe, chipe_exp, chiz, chiz_exp   
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
      call init_work_fields
      counter = 1
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      !---------------  RK 1st step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, aaa, .true.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg1)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi(i,k,j), psi(i,k,j), aaa(i,k,j), &
                               nonlin(i,k,j,:), &
                               ky(j), kz(k), kprp2(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupe , nuz , nupe_exp , nuz_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iomg), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupe , nuz , nupe_exp , nuz_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etape, etaz, etape_exp, etaz_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ipsi), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etape, etaz, etape_exp, etaz_exp)

            call get_imp_terms_tintg(imp_terms_tintg0(iaaa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     chipe , chiz , chipe_exp , chiz_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iaaa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     chipe , chiz , chipe_exp , chiz_exp )

            call eSSPIFRK1(omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), &
               imp_terms_tintg1(iomg), imp_terms_tintg0(iomg) &
            )
            call eSSPIFRK1(psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), &
               imp_terms_tintg1(ipsi), imp_terms_tintg0(ipsi) &
            )
            call eSSPIFRK1(aaa_tmp(i,k,j), aaa(i,k,j), &
               exp_terms(i,k,j,iaaa), &
               imp_terms_tintg1(iaaa), imp_terms_tintg0(iaaa) &
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

      ! Dealiasing
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi_tmp(i,k,j) = phi_tmp(i,k,j)*filter(i,k,j)
            psi_tmp(i,k,j) = psi_tmp(i,k,j)*filter(i,k,j)
            omg_tmp(i,k,j) = omg_tmp(i,k,j)*filter(i,k,j)
            aaa_tmp(i,k,j) = aaa_tmp(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo

      !---------------  RK 2nd step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, aaa_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi_tmp(i,k,j), psi_tmp(i,k,j), aaa_tmp(i,k,j), &
                               nonlin(i,k,j,:), &
                               ky(j), kz(k), kprp2(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupe , nuz , nupe_exp , nuz_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iomg), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupe , nuz , nupe_exp , nuz_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etape, etaz, etape_exp, etaz_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ipsi), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etape, etaz, etape_exp, etaz_exp)

            call get_imp_terms_tintg(imp_terms_tintg0(iaaa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     chipe , chiz , chipe_exp , chiz_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iaaa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     chipe , chiz , chipe_exp , chiz_exp )

            call eSSPIFRK2(omg_tmp(i,k,j), omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), &
               imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )

            call eSSPIFRK2(psi_tmp(i,k,j), psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), &
               imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )

            call eSSPIFRK2(aaa_tmp(i,k,j), aaa_tmp(i,k,j), aaa(i,k,j), &
               exp_terms(i,k,j,iaaa), &
               imp_terms_tintg2(iaaa), imp_terms_tintg0(iaaa) &
            )

            phi_tmp(i,k,j) = omg_tmp(i,k,j)*(-kprp2inv(i,k,j))
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! Dealiasing
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi_tmp(i,k,j) = phi_tmp(i,k,j)*filter(i,k,j)
            psi_tmp(i,k,j) = psi_tmp(i,k,j)*filter(i,k,j)
            omg_tmp(i,k,j) = omg_tmp(i,k,j)*filter(i,k,j)
            aaa_tmp(i,k,j) = aaa_tmp(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo

      !---------------  RK 3rd step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, aaa_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi_tmp(i,k,j), psi_tmp(i,k,j), aaa_tmp(i,k,j), &
                               nonlin(i,k,j,:), &
                               ky(j), kz(k), kprp2(i,k,j))

            ! Calculate time integral of explicit terms
            call get_imp_terms_tintg(imp_terms_tintg0(iomg),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupe , nuz , nupe_exp , nuz_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iomg), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupe , nuz , nupe_exp , nuz_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iomg),           dt, kprp2(i,k,j), kz2(k), &
                                     nupe , nuz , nupe_exp , nuz_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ipsi),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etape, etaz, etape_exp, etaz_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ipsi), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etape, etaz, etape_exp, etaz_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ipsi),           dt, kprp2(i,k,j), kz2(k), &
                                     etape, etaz, etape_exp, etaz_exp)

            call get_imp_terms_tintg(imp_terms_tintg0(iaaa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     chipe , chiz , chipe_exp , chiz_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iaaa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     chipe , chiz , chipe_exp , chiz_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iaaa),           dt, kprp2(i,k,j), kz2(k), &
                                     chipe , chiz , chipe_exp , chiz_exp )


            call eSSPIFRK3(omg_new(i,k,j), omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), exp_terms0(i,k,j,iomg), &
               imp_terms_tintg3(iomg), imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )
            call eSSPIFRK3(psi_new(i,k,j), psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), exp_terms0(i,k,j,ipsi), &
               imp_terms_tintg3(ipsi), imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )
            call eSSPIFRK3(aaa_new(i,k,j), aaa_tmp(i,k,j), aaa(i,k,j), &
               exp_terms(i,k,j,iaaa), exp_terms0(i,k,j,iaaa), &
               imp_terms_tintg3(iaaa), imp_terms_tintg2(iaaa), imp_terms_tintg0(iaaa) &
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
      aaa_old = aaa
      !$omp end workshare

      ! Dealiasing
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi_new(i,k,j) = phi_new(i,k,j)*filter(i,k,j)
            psi_new(i,k,j) = psi_new(i,k,j)*filter(i,k,j)
            omg_new(i,k,j) = omg_new(i,k,j)*filter(i,k,j)
            aaa_new(i,k,j) = aaa_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, aaa, .true.)

      !$omp parallel do private(j, k, i)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            exp_terms(i,k,j,iomg) = nonlin(i,k,j,iomg) - zi*kz(k)*kprp2(i,k,j)*psi(i,k,j) &
                                         - zi*ky(j)*aaa(i,k,j)
            exp_terms(i,k,j,ipsi) = nonlin(i,k,j,ipsi) + zi*kz(k)*phi(i,k,j)
            exp_terms(i,k,j,iaaa) = nonlin(i,k,j,iaaa) - zi*ky(j)*phi(i,k,j)
            ! 1st order 
            if(counter == 1) then
              call gear1(omg_new(i,k,j), omg(i,k,j), &
                 exp_terms(i,k,j,iomg), &
                 nupe*(kprp2(i,k,j)/kprp2_max)**nupe_exp + nuz*(kz2(k)/kz2_max)**nuz_exp &
              )
              call gear1(psi_new(i,k,j), psi(i,k,j), &
                 exp_terms(i,k,j,ipsi), &
                 etape*(kprp2(i,k,j)/kprp2_max)**etape_exp + etaz*(kz2(k)/kz2_max)**etaz_exp &
              )
              call gear1(aaa_new(i,k,j), aaa(i,k,j), &
                 exp_terms(i,k,j,iaaa), &
                 chipe*(kprp2(i,k,j)/kprp2_max)**chipe_exp + chiz*(kz2(k)/kz2_max)**chiz_exp &
              )

            ! 2nd order 
            elseif(counter == 2) then
              call gear2(omg_new(i,k,j), omg(i,k,j), omg_old(i,k,j), &
                 exp_terms    (i,k,j,iomg), &
                 exp_terms_old(i,k,j,iomg), &
                 nupe*(kprp2(i,k,j)/kprp2_max)**nupe_exp + nuz*(kz2(k)/kz2_max)**nuz_exp &
              )
              call gear2(psi_new(i,k,j), psi(i,k,j), psi_old(i,k,j), &
                 exp_terms    (i,k,j,ipsi), &
                 exp_terms_old(i,k,j,ipsi), &
                 etape*(kprp2(i,k,j)/kprp2_max)**etape_exp + etaz*(kz2(k)/kz2_max)**etaz_exp &
              )
              call gear2(aaa_new(i,k,j), aaa(i,k,j), aaa_old(i,k,j), &
                 exp_terms    (i,k,j,iaaa), &
                 exp_terms_old(i,k,j,iaaa), &
                 chipe*(kprp2(i,k,j)/kprp2_max)**chipe_exp + chiz*(kz2(k)/kz2_max)**chiz_exp &
              )

            ! 3rd order 
            else
              call gear3(omg_new(i,k,j), omg(i,k,j), omg_old(i,k,j), omg_old2(i,k,j), &
                 exp_terms     (i,k,j,iomg), &
                 exp_terms_old (i,k,j,iomg), &
                 exp_terms_old2(i,k,j,iomg), &
                 nupe*(kprp2(i,k,j)/kprp2_max)**nupe_exp + nuz*(kz2(k)/kz2_max)**nuz_exp &
              )
              call gear3(psi_new(i,k,j), psi(i,k,j), psi_old(i,k,j), psi_old2(i,k,j), &
                 exp_terms     (i,k,j,ipsi), &
                 exp_terms_old (i,k,j,ipsi), &
                 exp_terms_old2(i,k,j,ipsi), &
                 etape*(kprp2(i,k,j)/kprp2_max)**etape_exp + etaz*(kz2(k)/kz2_max)**etaz_exp &
              )
              call gear3(aaa_new(i,k,j), aaa(i,k,j), aaa_old(i,k,j), aaa_old2(i,k,j), &
                 exp_terms     (i,k,j,iaaa), &
                 exp_terms_old (i,k,j,iaaa), &
                 exp_terms_old2(i,k,j,iaaa), &
                 chipe*(kprp2(i,k,j)/kprp2_max)**chipe_exp + chiz*(kz2(k)/kz2_max)**chiz_exp &
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

      aaa_old2 = aaa_old 
      aaa_old  = aaa

      exp_terms_old2 = exp_terms_old
      exp_terms_old  = exp_terms
      !$omp end workshare

      ! Dealiasing
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi_new(i,k,j) = phi_new(i,k,j)*filter(i,k,j)
            psi_new(i,k,j) = psi_new(i,k,j)*filter(i,k,j)
            omg_new(i,k,j) = omg_new(i,k,j)*filter(i,k,j)
            aaa_new(i,k,j) = aaa_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
    endif

    do k = ikz_st, ikz_en
      if(kz(k) == 0.d0) then
        phi_new(:,k,:) = 0.d0
        psi_new(:,k,:) = 0.d0
        omg_new(:,k,:) = 0.d0
        aaa_new(:,k,:) = 0.d0
      endif
    enddo

    !$omp workshare
    phi = phi_new
    omg = omg_new
    psi = psi_new
    aaa = aaa_new
    !$omp end workshare

    tt = tt  + dt

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Allocate fields used only here
!-----------------------------------------------!
  subroutine init_work_fields
    use file, only: open_output_file
    use mp, only: proc0
    use params, only: time_step_scheme
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src

    call allocate_advance

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))

      allocate(phi_old2, source=src)
      allocate(omg_old2, source=src)
      allocate(psi_old2, source=src)
      allocate(aaa_old2, source=src)

      allocate(exp_terms_old (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old2 = 0.d0

      deallocate(src)
    endif

    if(proc0) call open_output_file (cfl_unit, 'cfl.dat')

  end subroutine init_work_fields


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Calculate nonlinear terms via
!!          1. Calculate grad in Fourier space
!!          2. Inverse FFT
!!          3. Calculate nonlineaer terms 
!!             in real space
!!          4. Forward FFT
!-----------------------------------------------!
  subroutine get_nonlinear_terms(phi, psi, aaa, dt_reset)
    use grid, only: dlx, dly
    use grid, only: kprp2, nlx, nly, nlz, kx, ky
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    use advance_common, only: dt_adjust_while_running 
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: phi, psi, aaa

    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (r8) :: max_vel_x, max_vel_y, dt_cfl

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    ! 1. Calculate grad in Fourier space
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          w(i,k,j,idphi_dx) = zi*kx(i)                  *phi(i, k, j)
          w(i,k,j,idphi_dy) = zi*ky(j)                  *phi(i, k, j)
          w(i,k,j,idomg_dx) = zi*kx(i)*(-kprp2(i, k, j))*phi(i, k, j)
          w(i,k,j,idomg_dy) = zi*ky(j)*(-kprp2(i, k, j))*phi(i, k, j)
                               
          w(i,k,j,idpsi_dx) = zi*kx(i)                  *psi(i, k, j)
          w(i,k,j,idpsi_dy) = zi*ky(j)                  *psi(i, k, j)
          w(i,k,j,idjpa_dx) = zi*kx(i)*(-kprp2(i, k, j))*psi(i, k, j)
          w(i,k,j,idjpa_dy) = zi*ky(j)*(-kprp2(i, k, j))*psi(i, k, j)

          w(i,k,j,idaaa_dx) = zi*kx(i)                  *aaa(i, k, j)
          w(i,k,j,idaaa_dy) = zi*ky(j)                  *aaa(i, k, j)
        enddo
      enddo
    enddo
    !$omp end parallel do


    ! 2. Inverse FFT
    if (proc0) call put_time_stamp(timer_fft)
    ! for some reason c2r_many doesn't work at Fugaku
    !call p3dfft_btran_c2r_many(w, nk_local_tot, w_r, nl_local_tot, nbtran, 'tff')
    do i = 1, nbtran
      call p3dfft_btran_c2r(w(:,:,:,i), w_r(:,:,:,i), 'tff')
    enddo
    if (proc0) call put_time_stamp(timer_fft)

    ! (get max_vel for dt reset)
    if(dt_reset) then
      max_vel_x = maxval(abs(w_r(:,:,:,idphi_dy)))
      max_vel_y = maxval(abs(w_r(:,:,:,idphi_dx)))
      call max_allreduce(max_vel_x)
      call max_allreduce(max_vel_y)
      dt_cfl = cfl*min(dlx/max_vel_x, dly/max_vel_y)

      if(proc0) then
        write (unit=cfl_unit, fmt="(100es30.21)") tt, dt_cfl, max_vel_x, max_vel_y
        call flush(cfl_unit) 
      endif

      call reset_dt(dt_cfl, counter)

    endif

    ! When the file 'dt_adjust' including a float number is created,
    ! dt will be manually adjusted to that value while running.
    call dt_adjust_while_running() 

    ! 3. Calculate nonlineaer terms in real space
    !$omp parallel do private(j, k) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          nonlin_r(j,k,i,iomg) = - w_r(j,k,i,idphi_dx)*w_r(j,k,i,idomg_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idomg_dx) &
                                 + w_r(j,k,i,idpsi_dx)*w_r(j,k,i,idjpa_dy) &
                                 - w_r(j,k,i,idpsi_dy)*w_r(j,k,i,idjpa_dx)
          nonlin_r(j,k,i,ipsi) = - w_r(j,k,i,idphi_dx)*w_r(j,k,i,idpsi_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idpsi_dx)
          nonlin_r(j,k,i,iaaa) = - w_r(j,k,i,idphi_dx)*w_r(j,k,i,idaaa_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idaaa_dx)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 4. Forward FFT
    if (proc0) call put_time_stamp(timer_fft)
    ! for some reason r2c_many doesn't work at Fugaku
    !call p3dfft_ftran_r2c_many(nonlin_r, nl_local_tot, nonlin, nk_local_tot, nfields, 'fft')
    do i = 1, nfields
      call p3dfft_ftran_r2c(nonlin_r(:,:,:,i), nonlin(:,:,:,i), 'fft')
    enddo
    if (proc0) call put_time_stamp(timer_fft)
    !$omp workshare
    nonlin = nonlin/nlx/nly/nlz
    !$omp end workshare

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms, &
                           phi, psi, aaa, &
                           nonlin, &
                           ky, kz, kprp2)
    use params, only: zi
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: phi, psi, aaa
    complex(r8), intent(in ) :: nonlin(nfields)
    real(r8)   , intent(in)  :: ky, kz, kprp2

    exp_terms(iomg) = nonlin(iomg) - zi*kz*kprp2*psi - zi*ky*aaa
    exp_terms(ipsi) = nonlin(ipsi) + zi*kz*phi
    exp_terms(iaaa) = nonlin(iaaa) - zi*ky*phi

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

    allocate(phi_new, source=src)
    allocate(omg_new, source=src)
    allocate(psi_new, source=src)
    allocate(aaa_new, source=src)

    allocate(w        (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nbtran ), source=(0.d0, 0.d0))
    allocate(nonlin   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))
    allocate(w_r      (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nbtran ), source=0.d0)
    allocate(nonlin_r (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nfields), source=0.d0)
    allocate(exp_terms(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields), source=(0.d0, 0.d0))

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      allocate(phi_tmp, source=src)
      allocate(omg_tmp, source=src)
      allocate(psi_tmp, source=src)
      allocate(aaa_tmp, source=src)

      allocate(exp_terms0(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms0 = 0.d0
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

    deallocate(phi_new)
    deallocate(omg_new)
    deallocate(psi_new)
    deallocate(aaa_new)

    deallocate(w        )
    deallocate(nonlin   )
    deallocate(w_r      )
    deallocate(nonlin_r )
    deallocate(exp_terms)

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                For eSSPIFRK3                v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'eSSPIFRK3') then
      deallocate(phi_tmp)
      deallocate(omg_tmp)
      deallocate(psi_tmp)
      deallocate(aaa_tmp)

      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance

end module advance

