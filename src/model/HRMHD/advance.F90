!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!
include "../../advance_common.F90"
!*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*!

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Time stepping for HRMHD
!-----------------------------------------------!
module advance
  use p3dfft
  use fields, only: nfields
  use fields, only: iomg, ipsi, iupa, ibpa
  implicit none

  public solve, is_allocated, allocate_advance, deallocate_advance

  logical :: is_allocated = .false.

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: phi_new
  complex(r8), allocatable, dimension(:,:,:)   :: omg_new
  complex(r8), allocatable, dimension(:,:,:)   :: psi_new
  complex(r8), allocatable, dimension(:,:,:)   :: upa_new
  complex(r8), allocatable, dimension(:,:,:)   :: bpa_new
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
  complex(r8), allocatable, dimension(:,:,:)   :: upa_tmp
  complex(r8), allocatable, dimension(:,:,:)   :: bpa_tmp
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms0

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !v                  For Gear3                  v!
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  complex(r8), allocatable, dimension(:,:,:)   :: phi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: omg_old2
  complex(r8), allocatable, dimension(:,:,:)   :: psi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: upa_old2
  complex(r8), allocatable, dimension(:,:,:)   :: bpa_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: exp_terms_old, exp_terms_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fphi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fpsi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fupa_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fbpa_old2

  integer :: cfl_unit, rms_unit

  ! Backward FFT variables
  integer, parameter :: nbtran = 12
  integer, parameter :: idphi_dx = 1 , idphi_dy = 2
  integer, parameter :: idomg_dx = 3 , idomg_dy = 4
  integer, parameter :: idpsi_dx = 5 , idpsi_dy = 6
  integer, parameter :: idjpa_dx = 7 , idjpa_dy = 8
  integer, parameter :: idupa_dx = 9 , idupa_dy = 10
  integer, parameter :: idbpa_dx = 11, idbpa_dy = 12

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Solve the time evolution
!-----------------------------------------------!
  subroutine solve
    use fields, only: phi, omg, psi, upa, bpa
    use fields, only: phi_old, omg_old, psi_old, upa_old, bpa_old
    use grid, only: kprp2, kprp2inv, kz2, kprp2_max, kz2_max
    use grid, only: kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use dealias, only: filter
    use params, only: nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                      nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                      etape_x, etape_x_exp, etape_z, etape_z_exp, &
                      etapa_x, etapa_x_exp, etapa_z, etapa_z_exp, &
                      zi, nonlinear
    use force, only: fphi, fpsi, fupa, fbpa, fphi_old, fpsi_old, fupa_old, fbpa_old
    use force, only: driven, elsasser, update_force, get_force, normalize_force, get_force, normalize_force_els
    use force, only: fzppe, fzmpe, fzppa, fzmpa
    use advance_common, only: eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
    use advance_common, only: gear1, gear2, gear3
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
          call get_force('zppe', fzppe)
          call get_force('zmpe', fzmpe)
          call normalize_force_els(fzppe, fzmpe)
          fphi = fzppe + fzmpe
          fpsi = fzppe - fzmpe
        else
          call get_force('phi', fphi)
          call get_force('psi', fpsi)
          call get_force('upa', fupa)
          call get_force('bpa', fbpa)
          call normalize_force(fphi, fpsi)
        endif
      endif

      !---------------  RK 1st step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, upa, bpa, .true.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg1)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi(i,k,j), psi(i,k,j), upa(i,k,j), bpa(i,k,j), &
                               nonlin(i,k,j,:), &
                               fphi(i,k,j), fpsi(i,k,j), fupa(i,k,j), fbpa(i,k,j), &
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

            call get_imp_terms_tintg(imp_terms_tintg0(iupa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupa_x , nupa_z , nupa_x_exp , nupa_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg1(iupa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupa_x , nupa_z , nupa_x_exp , nupa_z_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ibpa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etapa_x, etapa_z, etapa_x_exp, etapa_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg1(ibpa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etapa_x, etapa_z, etapa_x_exp, etapa_z_exp)

            call eSSPIFRK1(omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), &
               imp_terms_tintg1(iomg), imp_terms_tintg0(iomg) &
            )
            call eSSPIFRK1(psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), &
               imp_terms_tintg1(ipsi), imp_terms_tintg0(ipsi) &
            )
            call eSSPIFRK1(upa_tmp(i,k,j), upa(i,k,j), &
               exp_terms(i,k,j,iupa), &
               imp_terms_tintg1(iupa), imp_terms_tintg0(iupa) &
            )
            call eSSPIFRK1(bpa_tmp(i,k,j), bpa(i,k,j), &
               exp_terms(i,k,j,ibpa), &
               imp_terms_tintg1(ibpa), imp_terms_tintg0(ibpa) &
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
            upa_tmp(i,k,j) = upa_tmp(i,k,j)*filter(i,k,j)
            bpa_tmp(i,k,j) = bpa_tmp(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo

      !---------------  RK 2nd step  ---------------
      ! Calcualte force terms
      if (driven) then
        ! at n + 2/3
        call update_force(2.d0/3.d0*dt)

        if(elsasser) then
          call get_force('zppe', fzppe)
          call get_force('zmpe', fzmpe)
          call normalize_force_els(fzppe, fzmpe)
          fphi = fzppe + fzmpe
          fpsi = fzppe - fzmpe
        else
          call get_force('phi', fphi)
          call get_force('psi', fpsi)
          call get_force('upa', fupa)
          call get_force('bpa', fbpa)
          call normalize_force(fphi, fpsi)
        endif

        ! go to n + 1
        call update_force(1.d0/3.d0*dt)
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, upa_tmp, bpa_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi_tmp(i,k,j), psi_tmp(i,k,j), upa_tmp(i,k,j), bpa_tmp(i,k,j), &
                               nonlin(i,k,j,:), &
                               fphi(i,k,j), fpsi(i,k,j), fupa(i,k,j), fbpa(i,k,j), &
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

            call get_imp_terms_tintg(imp_terms_tintg0(iupa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupa_x , nupa_z , nupa_x_exp , nupa_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iupa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupa_x , nupa_z , nupa_x_exp , nupa_z_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ibpa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etapa_x, etapa_z, etapa_x_exp, etapa_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibpa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etapa_x, etapa_z, etapa_x_exp, etapa_z_exp)

            call eSSPIFRK2(omg_tmp(i,k,j), omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), &
               imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )

            call eSSPIFRK2(psi_tmp(i,k,j), psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), &
               imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )
            call eSSPIFRK2(upa_tmp(i,k,j), upa_tmp(i,k,j), upa(i,k,j), &
               exp_terms(i,k,j,iupa), &
               imp_terms_tintg2(iupa), imp_terms_tintg0(iupa) &
            )
            call eSSPIFRK2(bpa_tmp(i,k,j), bpa_tmp(i,k,j), bpa(i,k,j), &
               exp_terms(i,k,j,ibpa), &
               imp_terms_tintg2(ibpa), imp_terms_tintg0(ibpa) &
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
            upa_tmp(i,k,j) = upa_tmp(i,k,j)*filter(i,k,j)
            bpa_tmp(i,k,j) = bpa_tmp(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo

      !---------------  RK 3rd step  ---------------
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi_tmp, psi_tmp, upa_tmp, bpa_tmp, .false.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg2, imp_terms_tintg3)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi_tmp(i,k,j), psi_tmp(i,k,j), upa_tmp(i,k,j), bpa_tmp(i,k,j), &
                               nonlin(i,k,j,:), &
                               fphi(i,k,j), fpsi(i,k,j), fupa(i,k,j), fbpa(i,k,j), &
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

            call get_imp_terms_tintg(imp_terms_tintg0(iupa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     nupa_x , nupa_z , nupa_x_exp , nupa_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg2(iupa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     nupa_x , nupa_z , nupa_x_exp , nupa_z_exp )
            call get_imp_terms_tintg(imp_terms_tintg3(iupa),           dt, kprp2(i,k,j), kz2(k), &
                                     nupa_x , nupa_z , nupa_x_exp , nupa_z_exp )

            call get_imp_terms_tintg(imp_terms_tintg0(ibpa),         0.d0, kprp2(i,k,j), kz2(k), &
                                     etapa_x, etapa_z, etapa_x_exp, etapa_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg2(ibpa), 2.d0/3.d0*dt, kprp2(i,k,j), kz2(k), &
                                     etapa_x, etapa_z, etapa_x_exp, etapa_z_exp)
            call get_imp_terms_tintg(imp_terms_tintg3(ibpa),           dt, kprp2(i,k,j), kz2(k), &
                                     etapa_x, etapa_z, etapa_x_exp, etapa_z_exp)

            call eSSPIFRK3(omg_new(i,k,j), omg_tmp(i,k,j), omg(i,k,j), &
               exp_terms(i,k,j,iomg), exp_terms0(i,k,j,iomg), &
               imp_terms_tintg3(iomg), imp_terms_tintg2(iomg), imp_terms_tintg0(iomg) &
            )
            call eSSPIFRK3(psi_new(i,k,j), psi_tmp(i,k,j), psi(i,k,j), &
               exp_terms(i,k,j,ipsi), exp_terms0(i,k,j,ipsi), &
               imp_terms_tintg3(ipsi), imp_terms_tintg2(ipsi), imp_terms_tintg0(ipsi) &
            )
            call eSSPIFRK3(upa_new(i,k,j), upa_tmp(i,k,j), upa(i,k,j), &
               exp_terms(i,k,j,iupa), exp_terms0(i,k,j,iupa), &
               imp_terms_tintg3(iupa), imp_terms_tintg2(iupa), imp_terms_tintg0(iupa) &
            )
            call eSSPIFRK3(bpa_new(i,k,j), bpa_tmp(i,k,j), bpa(i,k,j), &
               exp_terms(i,k,j,ibpa), exp_terms0(i,k,j,ibpa), &
               imp_terms_tintg3(ibpa), imp_terms_tintg2(ibpa), imp_terms_tintg0(ibpa) &
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
      upa_old = upa
      bpa_old = bpa
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fphi_old = fphi
        fpsi_old = fpsi
        fupa_old = fupa
        fbpa_old = fbpa
        !$omp end workshare
      endif

      ! Dealiasing
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi_new(i,k,j) = phi_new(i,k,j)*filter(i,k,j)
            psi_new(i,k,j) = psi_new(i,k,j)*filter(i,k,j)
            omg_new(i,k,j) = omg_new(i,k,j)*filter(i,k,j)
            upa_new(i,k,j) = upa_new(i,k,j)*filter(i,k,j)
            bpa_new(i,k,j) = bpa_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
    endif

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !v                  For Gear3                  v!
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    if(time_step_scheme == 'gear3') then
      ! Calcualte force terms
      if (driven) then
        call update_force(dt)

        if(elsasser) then
          call get_force('zppe', fzppe)
          call get_force('zmpe', fzmpe)
          call normalize_force_els(fzppe, fzmpe)
          fphi = fzppe + fzmpe
          fpsi = fzppe - fzmpe
        else
          call get_force('phi', fphi)
          call get_force('psi', fpsi)
          call get_force('upa', fupa)
          call get_force('bpa', fbpa)
          call normalize_force(fphi, fpsi)
        endif
      endif

      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, upa, bpa, .true.)

      !$omp parallel do private(j, k, i)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi(i,k,j), psi(i,k,j), upa(i,k,j), bpa(i,k,j), &
                               nonlin(i,k,j,:), &
                               fphi(i,k,j), fpsi(i,k,j), fupa(i,k,j), fbpa(i,k,j), &
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
              call gear1(upa_new(i,k,j), upa(i,k,j), &
                 exp_terms(i,k,j,iupa), &
                 nupa_x*(kprp2(i,k,j)/kprp2_max)**nupa_x_exp + nupa_z*(kz2(k)/kz2_max)**nupa_z_exp &
              )
              call gear1(bpa_new(i,k,j), bpa(i,k,j), &
                 exp_terms(i,k,j,ibpa), &
                 etapa_x*(kprp2(i,k,j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp &
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
              call gear2(upa_new(i,k,j), upa(i,k,j), upa_old(i,k,j), &
                 exp_terms    (i,k,j,iupa), &
                 exp_terms_old(i,k,j,iupa), &
                 nupa_x*(kprp2(i,k,j)/kprp2_max)**nupa_x_exp + nupa_z*(kz2(k)/kz2_max)**nupa_z_exp &
              )
              call gear2(bpa_new(i,k,j), bpa(i,k,j), bpa_old(i,k,j), &
                 exp_terms    (i,k,j,ibpa), &
                 exp_terms_old(i,k,j,ibpa), &
                 etapa_x*(kprp2(i,k,j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp &
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
              call gear3(upa_new(i,k,j), upa(i,k,j), upa_old(i,k,j), upa_old2(i,k,j), &
                 exp_terms     (i,k,j,iupa), &
                 exp_terms_old (i,k,j,iupa), &
                 exp_terms_old2(i,k,j,iupa), &
                 nupa_x*(kprp2(i,k,j)/kprp2_max)**nupa_x_exp + nupa_z*(kz2(k)/kz2_max)**nupa_z_exp &
              )
              call gear3(bpa_new(i,k,j), bpa(i,k,j), bpa_old(i,k,j), bpa_old2(i,k,j), &
                 exp_terms     (i,k,j,ibpa), &
                 exp_terms_old (i,k,j,ibpa), &
                 exp_terms_old2(i,k,j,ibpa), &
                 etapa_x*(kprp2(i,k,j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp &
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

      upa_old2 = upa_old 
      upa_old  = upa

      bpa_old2 = bpa_old 
      bpa_old  = bpa

      exp_terms_old2 = exp_terms_old
      exp_terms_old  = exp_terms
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fphi_old2 = fphi_old 
        fphi_old  = fphi

        fpsi_old2 = fpsi_old 
        fpsi_old  = fpsi

        fupa_old2 = fupa_old 
        fupa_old  = fupa

        fbpa_old2 = fbpa_old 
        fbpa_old  = fbpa
        !$omp end workshare
      endif

      ! Dealiasing
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi_new(i,k,j) = phi_new(i,k,j)*filter(i,k,j)
            psi_new(i,k,j) = psi_new(i,k,j)*filter(i,k,j)
            omg_new(i,k,j) = omg_new(i,k,j)*filter(i,k,j)
            upa_new(i,k,j) = upa_new(i,k,j)*filter(i,k,j)
            bpa_new(i,k,j) = bpa_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
    endif

    !$omp workshare
    phi = phi_new
    omg = omg_new
    psi = psi_new
    upa = upa_new
    bpa = bpa_new
    !$omp end workshare

    tt = tt  + dt

    if(series_output) call output_series_modes

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Allocate fields used only here
!-----------------------------------------------!
  subroutine init_work_fields
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: time_step_scheme
    use file, only: open_output_file
    use mp, only: proc0
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
      allocate(upa_old2, source=src)
      allocate(bpa_old2, source=src)

      allocate(fphi_old2, source=src)
      allocate(fpsi_old2, source=src)

      allocate(exp_terms_old (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old  = 0.d0
      allocate(exp_terms_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nfields)); exp_terms_old2 = 0.d0

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
!!          3. Calculate nonlineaer terms 
!!             in real space
!!          4. Forward FFT
!-----------------------------------------------!
  subroutine get_nonlinear_terms(phi, psi, upa, bpa, dt_reset)
    use grid, only: dlx, dly, dlz
    use grid, only: kprp2, nlx, nly, nlz, kx, ky
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi, sgm
    use mp, only: proc0, max_allreduce, sum_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    use advance_common, only: dt_adjust_while_running 
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: phi, psi, upa, bpa

    logical, intent(in) :: dt_reset

    integer :: i, j, k
    real   (r8), allocatable :: v_KAW(:,:,:)
    real   (r8) :: max_vel_x, max_vel_y, max_vel_KAW, dt_cfl, dt_digit
    real   (r8) :: ux_rms, uy_rms, bx_rms, by_rms

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    allocate(v_KAW(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=0.d0)

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

          w(i,k,j,idupa_dx) = zi*kx(i)                  *upa(i, k, j)
          w(i,k,j,idupa_dy) = zi*ky(j)                  *upa(i, k, j)
          w(i,k,j,idbpa_dx) = zi*kx(i)                  *bpa(i, k, j)
          w(i,k,j,idbpa_dy) = zi*ky(j)                  *bpa(i, k, j)

          v_KAW(i,k,j) = sqrt(0.5d0*( &
                                    1.d0 + sgm**2 + kprp2(i,k,j) &
                                  + sqrt( (kprp2(i,k,j) + (1.d0 + sgm)**2)*(kprp2(i,k,j) + (1.d0 - sgm)**2) ) &
                         ))
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


    ! write rms
    ux_rms = sum(w_r(:,:,:,idphi_dy)**2); uy_rms = sum(w_r(:,:,:,idphi_dx)**2)
    bx_rms = sum(w_r(:,:,:,idpsi_dy)**2); by_rms = sum(w_r(:,:,:,idpsi_dx)**2)

    call sum_allreduce(ux_rms); call sum_allreduce(uy_rms)
    call sum_allreduce(bx_rms); call sum_allreduce(by_rms)

    ux_rms = sqrt(ux_rms/nlx/nly/nlz); uy_rms = sqrt(uy_rms/nlx/nly/nlz)
    bx_rms = sqrt(bx_rms/nlx/nly/nlz); by_rms = sqrt(by_rms/nlx/nly/nlz)

    if(proc0) then
      write (unit=rms_unit, fmt="(100es30.21)") tt, ux_rms, uy_rms, bx_rms, by_rms
      call flush(rms_unit) 
    endif


    ! (get max_vel for dt reset)
    if(dt_reset) then
      max_vel_x   = maxval(abs(w_r(:,:,:,idphi_dy)))
      max_vel_y   = maxval(abs(w_r(:,:,:,idphi_dx)))
      max_vel_KAW = maxval(abs(v_KAW(:,:,:)      ))
      call max_allreduce(max_vel_x  )
      call max_allreduce(max_vel_y  )
      call max_allreduce(max_vel_KAW)
      dt_cfl = cfl*min(dlx/max_vel_x, dly/max_vel_y, dlz/max_vel_KAW)

      if(proc0) then
        write (unit=cfl_unit, fmt="(100es30.21)") tt, dt_cfl, max_vel_x, max_vel_y
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

    ! 3. Calculate nonlineaer terms in real space
    !$omp parallel do private(j, k) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          nonlin_r(j,k,i,iomg) = - w_r(j,k,i,idphi_dx)*w_r(j,k,i,idomg_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idomg_dx) &
                                 + w_r(j,k,i,idpsi_dx)*w_r(j,k,i,idjpa_dy) &
                                 - w_r(j,k,i,idpsi_dy)*w_r(j,k,i,idjpa_dx)
          nonlin_r(j,k,i,ipsi) =   w_r(j,k,i,idpsi_dx)*(w_r(j,k,i,idphi_dy) + w_r(j,k,i,idbpa_dy))&
                                 - w_r(j,k,i,idpsi_dy)*(w_r(j,k,i,idphi_dx) + w_r(j,k,i,idbpa_dx))
          nonlin_r(j,k,i,iupa) = - w_r(j,k,i,idphi_dx)*w_r(j,k,i,idupa_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idupa_dx) &
                                 + w_r(j,k,i,idpsi_dx)*sgm*w_r(j,k,i,idbpa_dy) &
                                 - w_r(j,k,i,idpsi_dy)*sgm*w_r(j,k,i,idbpa_dx)
          nonlin_r(j,k,i,ibpa) = - w_r(j,k,i,idphi_dx)*w_r(j,k,i,idbpa_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idbpa_dx) &
                                 + w_r(j,k,i,idpsi_dx)*(sgm*w_r(j,k,i,idupa_dy) - w_r(j,k,i,idjpa_dy))&
                                 - w_r(j,k,i,idpsi_dy)*(sgm*w_r(j,k,i,idupa_dx) - w_r(j,k,i,idjpa_dx))
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
    !
    deallocate(v_KAW)

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Calculate explicit terms
!-----------------------------------------------!
  subroutine get_ext_terms(exp_terms, &
                           phi, psi, upa, bpa, &
                           nonlin, &
                           fphi, fpsi, fupa, fbpa, &
                           kz, kprp2)
    use params, only: zi, sgm
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: phi, psi, upa, bpa
    complex(r8), intent(in ) :: nonlin(nfields)
    complex(r8), intent(in ) :: fphi, fpsi, fupa, fbpa
    real(r8)   , intent(in)  :: kz, kprp2

    exp_terms(iomg) = nonlin(iomg) - kprp2*fphi - zi*kz*kprp2*psi
    exp_terms(ipsi) = nonlin(ipsi) + fpsi + zi*kz*(phi + bpa)
    exp_terms(iupa) = nonlin(iupa) + fupa + zi*kz*sgm*bpa
    exp_terms(ibpa) = nonlin(ibpa) + fbpa + zi*kz*(sgm*upa + kprp2*psi)

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
!! @date    15 Jul 2021
!! @brief   Output series modes
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: phi, psi, upa, bpa
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, sum_reduce
    use time, only: tt
    use diagnostics_common, only: n_series_modes, series_modes
    use diagnostics_common, only: series_modes_unit
    use diagnostics, only: split_KAW_ICW
    implicit none
    complex(r8), dimension(n_series_modes) :: phi_modes, psi_modes, upa_modes, bpa_modes
    complex(r8), dimension(n_series_modes) :: phi_KAW_modes, psi_KAW_modes, upa_KAW_modes, bpa_KAW_modes
    complex(r8), dimension(n_series_modes) :: phi_ICW_modes, psi_ICW_modes, upa_ICW_modes, bpa_ICW_modes
    complex(r8), allocatable, dimension(:,:,:) :: src_c
    complex(r8), allocatable, dimension(:,:,:) :: phi_KAW, psi_KAW, upa_KAW, bpa_KAW
    complex(r8), allocatable, dimension(:,:,:) :: phi_ICW, psi_ICW, upa_ICW, bpa_ICW
    integer :: n, i, j, k

    allocate(src_c(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))
    allocate(phi_KAW, source=src_c)
    allocate(psi_KAW, source=src_c)
    allocate(upa_KAW, source=src_c)
    allocate(bpa_KAW, source=src_c)
    allocate(phi_ICW, source=src_c)
    allocate(psi_ICW, source=src_c)
    allocate(upa_ICW, source=src_c)
    allocate(bpa_ICW, source=src_c)
    deallocate(src_c)

    phi_modes(:) = 0.d0
    psi_modes(:) = 0.d0
    upa_modes(:) = 0.d0
    bpa_modes(:) = 0.d0

    call split_KAW_ICW( &
                        phi    , psi    , upa    , bpa    , &
                        phi_KAW, psi_KAW, upa_KAW, bpa_KAW, &
                        phi_ICW, psi_ICW, upa_ICW, bpa_ICW  &
                      )

    do n = 1, n_series_modes
      i = series_modes(n, 1)
      j = series_modes(n, 2)
      k = series_modes(n, 3)

      if(       (i >= ikx_st .and. i <= ikx_en) &
          .and. (j >= iky_st .and. j <= iky_en) &
          .and. (k >= ikz_st .and. k <= ikz_en) &
        ) then

        phi_modes    (n) = phi    (i, k, j)
        psi_modes    (n) = psi    (i, k, j)
        upa_modes    (n) = upa    (i, k, j)
        bpa_modes    (n) = bpa    (i, k, j)

        phi_KAW_modes(n) = phi_KAW(i, k, j)
        psi_KAW_modes(n) = psi_KAW(i, k, j)
        upa_KAW_modes(n) = upa_KAW(i, k, j)
        bpa_KAW_modes(n) = bpa_KAW(i, k, j)

        phi_ICW_modes(n) = phi_ICW(i, k, j)
        psi_ICW_modes(n) = psi_ICW(i, k, j)
        upa_ICW_modes(n) = upa_ICW(i, k, j)
        bpa_ICW_modes(n) = bpa_ICW(i, k, j)

      endif
    enddo

    call sum_reduce(phi_modes    , 0)
    call sum_reduce(psi_modes    , 0)
    call sum_reduce(upa_modes    , 0)
    call sum_reduce(bpa_modes    , 0)

    call sum_reduce(phi_KAW_modes, 0)
    call sum_reduce(psi_KAW_modes, 0)
    call sum_reduce(upa_KAW_modes, 0)
    call sum_reduce(bpa_KAW_modes, 0)

    call sum_reduce(phi_ICW_modes, 0)
    call sum_reduce(psi_ICW_modes, 0)
    call sum_reduce(upa_ICW_modes, 0)
    call sum_reduce(bpa_ICW_modes, 0)

    do n = 1, n_series_modes
      if(proc0) then
        i = series_modes(n, 1)
        j = series_modes(n, 2)
        k = series_modes(n, 3)
999 format(es30.21, A10, 5es30.21e3)
        write (unit=series_modes_unit, fmt=999) tt, 'phi'    , kx(i), ky(j), kz(k), phi_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'psi'    , kx(i), ky(j), kz(k), psi_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'upa'    , kx(i), ky(j), kz(k), upa_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bpa'    , kx(i), ky(j), kz(k), bpa_modes(n)

        write (unit=series_modes_unit, fmt=999) tt, 'phi_KAW', kx(i), ky(j), kz(k), phi_KAW_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'psi_KAW', kx(i), ky(j), kz(k), psi_KAW_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'upa_KAW', kx(i), ky(j), kz(k), upa_KAW_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bpa_KAW', kx(i), ky(j), kz(k), bpa_KAW_modes(n)

        write (unit=series_modes_unit, fmt=999) tt, 'phi_ICW', kx(i), ky(j), kz(k), phi_ICW_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'psi_ICW', kx(i), ky(j), kz(k), psi_ICW_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'upa_ICW', kx(i), ky(j), kz(k), upa_ICW_modes(n)
        write (unit=series_modes_unit, fmt=999) tt, 'bpa_ICW', kx(i), ky(j), kz(k), bpa_ICW_modes(n)

        call flush(series_modes_unit) 

      endif
    enddo

    deallocate(phi_KAW)
    deallocate(psi_KAW)
    deallocate(upa_KAW)
    deallocate(bpa_KAW)
    deallocate(phi_ICW)
    deallocate(psi_ICW)
    deallocate(upa_ICW)
    deallocate(bpa_ICW)

  end subroutine output_series_modes


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
    allocate(upa_new, source=src)
    allocate(bpa_new, source=src)

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
      allocate(upa_tmp, source=src)
      allocate(bpa_tmp, source=src)

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
    deallocate(upa_new)
    deallocate(bpa_new)

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
      deallocate(upa_tmp)
      deallocate(bpa_tmp)

      deallocate(exp_terms0)
    endif

  end subroutine deallocate_advance

end module advance

