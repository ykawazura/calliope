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

  integer :: cfl_unit

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
    use grid, only: ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use dealias, only: filter
    use params, only: va2cs2_plus_1, nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                                     nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                                     etape_x, etape_x_exp, etape_z, etape_z_exp, &
                                     etapa_x, etapa_x_exp, etapa_z, etapa_z_exp, &
                      zi, nonlinear, q
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
      if(nonlinear) call get_nonlinear_terms(phi, psi, upa, bpa, .true.)

      !$omp parallel do private(j, k, i, imp_terms_tintg0, imp_terms_tintg1)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            ! Calculate explicit terms
            call get_ext_terms(exp_terms(i,k,j,:), &
                               phi(i,k,j), psi(i,k,j), upa(i,k,j), bpa(i,k,j), &
                               nonlin(i,k,j,:), &
                               ky(j), kz(k), kprp2(i,k,j))

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
            imp_terms_tintg0(ibpa) = imp_terms_tintg0(ibpa)/va2cs2_plus_1
            imp_terms_tintg1(ibpa) = imp_terms_tintg1(ibpa)/va2cs2_plus_1

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
                               ky(j), kz(k), kprp2(i,k,j))

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
            imp_terms_tintg0(ibpa) = imp_terms_tintg0(ibpa)/va2cs2_plus_1
            imp_terms_tintg2(ibpa) = imp_terms_tintg2(ibpa)/va2cs2_plus_1

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
                               ky(j), kz(k), kprp2(i,k,j))

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
            imp_terms_tintg0(ibpa) = imp_terms_tintg0(ibpa)/va2cs2_plus_1
            imp_terms_tintg2(ibpa) = imp_terms_tintg2(ibpa)/va2cs2_plus_1
            imp_terms_tintg3(ibpa) = imp_terms_tintg3(ibpa)/va2cs2_plus_1

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
      ! Calcualte nonlinear terms
      if(nonlinear) call get_nonlinear_terms(phi, psi, upa, bpa, .true.)

      !$omp parallel do private(j, k, i)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en

            exp_terms(i,k,j,iomg) = nonlin(i,k,j,iomg) - zi*kz(k)*kprp2(i,k,j)*psi(i,k,j) &
                                         - 2.d0*zi*ky(j)*upa(i,k,j)
            exp_terms(i,k,j,ipsi) = nonlin(i,k,j,ipsi) + zi*kz(k)*phi(i,k,j)
            exp_terms(i,k,j,iupa) = nonlin(i,k,j,iupa) + zi*kz(k)*bpa(i,k,j) + (2.d0 - q)*zi*ky(j)*phi(i,k,j)
            exp_terms(i,k,j,ibpa) = (nonlin(i,k,j,ibpa) + zi*kz(k)*upa(i,k,j) + q*zi*ky(j)*psi(i,k,j) ) &
                                         /va2cs2_plus_1
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
                 (etapa_x*(kprp2(i,k,j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp)/va2cs2_plus_1 &
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
                 (etapa_x*(kprp2(i,k,j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp)/va2cs2_plus_1 &
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
                 (etapa_x*(kprp2(i,k,j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp)/va2cs2_plus_1 &
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

      allocate(phi_old2 , source=src)
      allocate(omg_old2 , source=src)
      allocate(psi_old2 , source=src)
      allocate(upa_old2, source=src)
      allocate(bpa_old2, source=src)

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
  subroutine get_nonlinear_terms(phi, psi, upa, bpa, dt_reset)
    use grid, only: dlx, dly
    use grid, only: kprp2, nlx, nly, nlz, kx, ky
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi, va2cs2_plus_1
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    use advance_common, only: dt_adjust_while_running 
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in) :: phi, psi, upa, bpa

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

          w(i,k,j,idupa_dx) = zi*kx(i)                  *upa(i, k, j)
          w(i,k,j,idupa_dy) = zi*ky(j)                  *upa(i, k, j)
          w(i,k,j,idbpa_dx) = zi*kx(i)                  *bpa(i, k, j)
          w(i,k,j,idbpa_dy) = zi*ky(j)                  *bpa(i, k, j)
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
          nonlin_r(j,k,i,iupa) = - w_r(j,k,i,idphi_dx)*w_r(j,k,i,idupa_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idupa_dx) &
                                 + w_r(j,k,i,idpsi_dx)*w_r(j,k,i,idbpa_dy) &
                                 - w_r(j,k,i,idpsi_dy)*w_r(j,k,i,idbpa_dx)
          nonlin_r(j,k,i,ibpa) =(- w_r(j,k,i,idphi_dx)*w_r(j,k,i,idbpa_dy) &
                                 + w_r(j,k,i,idphi_dy)*w_r(j,k,i,idbpa_dx))*va2cs2_plus_1 &
                                 + w_r(j,k,i,idpsi_dx)*w_r(j,k,i,idupa_dy) &
                                 - w_r(j,k,i,idpsi_dy)*w_r(j,k,i,idupa_dx)
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
                           phi, psi, upa, bpa, &
                           nonlin, &
                           ky, kz, kprp2)
    use params, only: zi, va2cs2_plus_1, q
    implicit none
    complex(r8), intent(out) :: exp_terms(nfields)
    complex(r8), intent(in ) :: phi, psi, upa, bpa
    complex(r8), intent(in ) :: nonlin(nfields)
    real(r8)   , intent(in)  :: ky, kz, kprp2

    exp_terms(iomg) = nonlin(iomg) - zi*kz*kprp2*psi - 2.d0*zi*ky*upa
    exp_terms(ipsi) = nonlin(ipsi) + zi*kz*phi
    exp_terms(iupa) = nonlin(iupa) + zi*kz*bpa + (2.d0 - q)*zi*ky*phi
    exp_terms(ibpa) = ( nonlin(ibpa) + zi*kz*upa + q*zi*ky*psi )/va2cs2_plus_1

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

    allocate(phi_new  , source=src)
    allocate(omg_new  , source=src)
    allocate(psi_new  , source=src)
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

