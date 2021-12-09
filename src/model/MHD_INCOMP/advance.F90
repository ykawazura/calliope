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
  implicit none

  public solve

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: ux_new, ux_old2
  complex(r8), allocatable, dimension(:,:,:)   :: uy_new, uy_old2
  complex(r8), allocatable, dimension(:,:,:)   :: uz_new, uz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: bx_new, bx_old2
  complex(r8), allocatable, dimension(:,:,:)   :: by_new, by_old2
  complex(r8), allocatable, dimension(:,:,:)   :: bz_new, bz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: p_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: flx
  complex(r8), allocatable, dimension(:,:,:)   :: nl_ux, nl_ux_old1, nl_ux_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_uy, nl_uy_old1, nl_uy_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_uz, nl_uz_old1, nl_uz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_bx, nl_bx_old1, nl_bx_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_by, nl_by_old1, nl_by_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_bz, nl_bz_old1, nl_bz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fux_old2, fuy_old2, fuz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nbl2inv_div_u ! nabla^-2 (div u)
  complex(r8), allocatable, dimension(:,:,:)   :: nbl2inv_div_b ! nabla^-2 (div b)
  real   (r8), allocatable, dimension(:,:)     :: kxt, kxt_old1, kxt_old2
  real   (r8) :: cflx, cfly, cflz
  integer :: max_vel_unit

  ! Backward FFT variables
  integer, parameter :: nbtran = 6
  integer, parameter :: iux = 1, iuy = 2, iuz = 3
  integer, parameter :: ibx = 4, iby = 5, ibz = 6
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
!! @brief   Solve the equation of motion
!!          3rd Order Gear's method
!!               linear terms: explicit
!!            nonlinear terms: explicit
!!          dissipation terms: implicit
!!          [Karniadakis and Israeli, JCP 1991]
!-----------------------------------------------!
  subroutine solve
    use fields, only: ux, uy, uz, p
    use fields, only: bx, by, bz
    use fields, only: ux_old1, uy_old1, uz_old1, p_old1
    use fields, only: bx_old1, by_old1, bz_old1
    use grid, only: k2_max
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, cfl, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nu, eta, nu_exp, eta_exp, zi, nonlinear, shear, q
    use params, only: dealias_scheme => dealias
    use dealias, only: filter
    use shearing_box, only: shear_flg, k2t, k2t_inv, tsc, tremap
    use force, only: fux, fuz, fuy, fux_old1, fuy_old1, fuz_old1, driven, update_force, get_force
    use advance_common, only: gear1, gear2, gear3
    use params, only: series_output
    implicit none
    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_advance)

    ! initialize tmp fields
    if(counter == 0) then
      call init_multistep_fields

      cflx = maxval(abs(kx))/cfl
      cfly = maxval(abs(ky))/cfl
      cflz = maxval(abs(kz))/cfl

      call init_pressure
      counter = 1
    endif

    ! Calcualte nonlinear terms
    if(nonlinear) call get_nonlinear_terms

    ! Calcualte force terms
    if (driven) then
      call update_force
      call get_force('ux', fux)
      call get_force('uy', fuy)
      call get_force('uz', fuz)
    endif

    ! 1st order 
    if(counter == 1) then
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            ! get p
            p(i,k,j) = -zi*( kxt(i,j)*nl_ux(i,k,j) &
                           + ky(j)   *nl_uy(i,k,j) &
                           + kz(k)   *nl_uz(i,k,j) )*k2t_inv(i,k,j)

            ! div (uu - bb)
            nl_ux(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxx) &
                                + ky(j)   *flx(i,k,j,iflx_uxy) &
                                + kz(k)   *flx(i,k,j,iflx_uxz) )
                                                          
            nl_uy(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxy) &
                                + ky(j)   *flx(i,k,j,iflx_uyy) &
                                + kz(k)   *flx(i,k,j,iflx_uyz) )
                                                          
            nl_uz(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxz) &
                                + ky(j)   *flx(i,k,j,iflx_uyz) &
                                + kz(k)   *flx(i,k,j,iflx_uzz) )

            ! curl (b x u)
            nl_bx(i,k,j) = -zi*(  ky(j)   *flx(i,k,j,iflx_bz ) &
                                - kz(k)   *flx(i,k,j,iflx_by ) )
                                                          
            nl_by(i,k,j) = -zi*(  kz(k)   *flx(i,k,j,iflx_bx ) &
                                - kxt(i,j)*flx(i,k,j,iflx_bz ) )
                                                          
            nl_bz(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_by ) &
                                - ky(j)   *flx(i,k,j,iflx_bx ) )

            ! update u
            call gear1(dt, ux_new(i,k,j), ux(i,k,j), &
               nl_ux(i,k,j) - zi*kxt(i,j)*p(i,k,j) + fux(i,k,j)  &
                 + 2.d0*shear_flg*uy(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )
            call gear1(dt, uy_new(i,k,j), uy(i,k,j), &
               nl_uy(i,k,j) - zi*ky(j)   *p(i,k,j) + fuy(i,k,j)  &
                 - (2.d0 - q)*shear_flg*ux(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )
            call gear1(dt, uz_new(i,k,j), uz(i,k,j), &
               nl_uz(i,k,j) - zi*kz(k)   *p(i,k,j) + fuz(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )

            ! update b
            call gear1(dt, bx_new(i,k,j), bx(i,k,j), &
               nl_bx(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
            call gear1(dt, by_new(i,k,j), by(i,k,j), &
               nl_by(i,k,j)  &
                 - q*shear_flg*bx(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
            call gear1(dt, bz_new(i,k,j), bz(i,k,j), &
               nl_bz(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      ux_old1 = ux
      nl_ux_old1 = nl_ux

      uy_old1 = uy
      nl_uy_old1 = nl_uy

      uz_old1 = uz
      nl_uz_old1 = nl_uz

      bx_old1 = bx
      nl_bx_old1 = nl_bx

      by_old1 = by
      nl_by_old1 = nl_by

      bz_old1 = bz
      nl_bz_old1 = nl_bz

      p_old1 =  p
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fux_old1 = fux
        fuy_old1 = fuy
        fuz_old1 = fuz
        !$omp end workshare
      endif

      if(shear) then
        !$omp workshare
        kxt_old1 = kxt
        !$omp end workshare
      endif

      counter = counter + 1
    ! 2nd order 
    elseif(counter == 2) then
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            ! get p
            p(i,k,j) = -zi*( kxt(i,j)*nl_ux(i,k,j) &
                           + ky(j)   *nl_uy(i,k,j) &
                           + kz(k)   *nl_uz(i,k,j) )*k2t_inv(i,k,j)

            ! div (uu - bb)
            nl_ux(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxx) &
                                + ky(j)   *flx(i,k,j,iflx_uxy) &
                                + kz(k)   *flx(i,k,j,iflx_uxz) )
                                                          
            nl_uy(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxy) &
                                + ky(j)   *flx(i,k,j,iflx_uyy) &
                                + kz(k)   *flx(i,k,j,iflx_uyz) )
                                                          
            nl_uz(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxz) &
                                + ky(j)   *flx(i,k,j,iflx_uyz) &
                                + kz(k)   *flx(i,k,j,iflx_uzz) )

            ! curl (b x u)
            nl_bx(i,k,j) = -zi*(  ky(j)   *flx(i,k,j,iflx_bz ) &
                                - kz(k)   *flx(i,k,j,iflx_by ) )
                                                          
            nl_by(i,k,j) = -zi*(  kz(k)   *flx(i,k,j,iflx_bx ) &
                                - kxt(i,j)*flx(i,k,j,iflx_bz ) )
                                                          
            nl_bz(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_by ) &
                                - ky(j)   *flx(i,k,j,iflx_bx ) )

            ! update u
            call gear2(dt, ux_new(i,k,j), ux(i,k,j), ux_old1(i,k,j), &
               nl_ux     (i,k,j) - zi*kxt     (i,j)*p     (i,k,j) + fux     (i,k,j)  &
                 + 2.d0*shear_flg*uy     (i,k,j), &
               nl_ux_old1(i,k,j) - zi*kxt_old1(i,j)*p_old1(i,k,j) + fux_old1(i,k,j)  &
                 + 2.d0*shear_flg*uy_old1(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )
            call gear2(dt, uy_new(i,k,j), uy(i,k,j), uy_old1(i,k,j), &
               nl_uy     (i,k,j) - zi*ky(j)*p     (i,k,j) + fuy     (i,k,j)  &
                 - (2.d0 - q)*shear_flg*ux     (i,k,j), &
               nl_uy_old1(i,k,j) - zi*ky(j)*p_old1(i,k,j) + fuy_old1(i,k,j)  &
                 - (2.d0 - q)*shear_flg*ux_old1(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )
            call gear2(dt, uz_new(i,k,j), uz(i,k,j), uz_old1(i,k,j), &
               nl_uz     (i,k,j) - zi*kz(k)*p     (i,k,j) + fuz     (i,k,j), &
               nl_uz_old1(i,k,j) - zi*kz(k)*p_old1(i,k,j) + fuz_old1(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )

            ! update b
            call gear2(dt, bx_new(i,k,j), bx(i,k,j), bx_old1(i,k,j), &
               nl_bx     (i,k,j), &
               nl_bx_old1(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
            call gear2(dt, by_new(i,k,j), by(i,k,j), by_old1(i,k,j), &
               nl_by     (i,k,j)  &
                 - q*shear_flg*bx     (i,k,j), &
               nl_by_old1(i,k,j)  &
                 - q*shear_flg*bx_old1(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
            call gear2(dt, bz_new(i,k,j), bz(i,k,j), bz_old1(i,k,j), &
               nl_bz     (i,k,j), &
               nl_bz_old1(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      ux_old2 = ux_old1
      ux_old1 = ux
      nl_ux_old2 = nl_ux_old1
      nl_ux_old1 = nl_ux

      uy_old2 = uy_old1
      uy_old1 = uy
      nl_uy_old2 = nl_uy_old1
      nl_uy_old1 = nl_uy

      uz_old2 = uz_old1
      uz_old1 = uz
      nl_uz_old2 = nl_uz_old1
      nl_uz_old1 = nl_uz

      bx_old2 = bx_old1
      bx_old1 = bx
      nl_bx_old2 = nl_bx_old1
      nl_bx_old1 = nl_bx

      by_old2 = by_old1
      by_old1 = by
      nl_by_old2 = nl_by_old1
      nl_by_old1 = nl_by

      bz_old2 = bz_old1
      bz_old1 = bz
      nl_bz_old2 = nl_bz_old1
      nl_bz_old1 = nl_bz

      p_old2 =  p_old1
      p_old1 =  p
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fux_old2 = fux_old1
        fux_old1 = fux

        fuy_old2 = fuy_old1
        fuy_old1 = fuy

        fuz_old2 = fuz_old1
        fuz_old1 = fuz
        !$omp end workshare
      endif

      if(shear) then
        !$omp workshare
        kxt_old2 = kxt_old1
        kxt_old1 = kxt
        !$omp end workshare
      endif

      counter = counter + 1
    ! 3rd order 
    else
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            ! get p
            p(i,k,j) = -zi*( kxt(i,j)*nl_ux(i,k,j) &
                           + ky(j)   *nl_uy(i,k,j) &
                           + kz(k)   *nl_uz(i,k,j) )*k2t_inv(i,k,j)

            ! div (uu - bb)
            nl_ux(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxx) &
                                + ky(j)   *flx(i,k,j,iflx_uxy) &
                                + kz(k)   *flx(i,k,j,iflx_uxz) )
                                                          
            nl_uy(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxy) &
                                + ky(j)   *flx(i,k,j,iflx_uyy) &
                                + kz(k)   *flx(i,k,j,iflx_uyz) )
                                                          
            nl_uz(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_uxz) &
                                + ky(j)   *flx(i,k,j,iflx_uyz) &
                                + kz(k)   *flx(i,k,j,iflx_uzz) )

            ! curl (b x u)
            nl_bx(i,k,j) = -zi*(  ky(j)   *flx(i,k,j,iflx_bz ) &
                                - kz(k)   *flx(i,k,j,iflx_by ) )
                                                          
            nl_by(i,k,j) = -zi*(  kz(k)   *flx(i,k,j,iflx_bx ) &
                                - kxt(i,j)*flx(i,k,j,iflx_bz ) )
                                                          
            nl_bz(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_by ) &
                                - ky(j)   *flx(i,k,j,iflx_bx ) )

            ! update u
            call gear3(dt, ux_new(i,k,j), ux(i,k,j), ux_old1(i,k,j), ux_old2(i,k,j), &
               nl_ux     (i,k,j) - zi*kxt     (i,j)*p     (i,k,j) + fux     (i,k,j)  &
                 + 2.d0*shear_flg*uy     (i,k,j), &
               nl_ux_old1(i,k,j) - zi*kxt_old1(i,j)*p_old1(i,k,j) + fux_old1(i,k,j)  &
                 + 2.d0*shear_flg*uy_old1(i,k,j), &
               nl_ux_old2(i,k,j) - zi*kxt_old2(i,j)*p_old2(i,k,j) + fux_old2(i,k,j)  &
                 + 2.d0*shear_flg*uy_old2(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )
            call gear3(dt, uy_new(i,k,j), uy(i,k,j), uy_old1(i,k,j), uy_old2(i,k,j), &
               nl_uy     (i,k,j) - zi*ky(j)*p     (i,k,j) + fuy     (i,k,j)  &
                 - (2.d0 - q)*shear_flg*ux     (i,k,j), &
               nl_uy_old1(i,k,j) - zi*ky(j)*p_old1(i,k,j) + fuy_old1(i,k,j)  &
                 - (2.d0 - q)*shear_flg*ux_old1(i,k,j), &
               nl_uy_old2(i,k,j) - zi*ky(j)*p_old2(i,k,j) + fuy_old2(i,k,j)  &
                 - (2.d0 - q)*shear_flg*ux_old2(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )
            call gear3(dt, uz_new(i,k,j), uz(i,k,j), uz_old1(i,k,j), uz_old2(i,k,j), &
               nl_uz     (i,k,j) - zi*kz(k)*p     (i,k,j) + fuz     (i,k,j), &
               nl_uz_old1(i,k,j) - zi*kz(k)*p_old1(i,k,j) + fuz_old1(i,k,j), &
               nl_uz_old2(i,k,j) - zi*kz(k)*p_old2(i,k,j) + fuz_old2(i,k,j), &
               nu*(k2t(i,k,j)/k2_max)**nu_exp &
            )

            ! update b
            call gear3(dt, bx_new(i,k,j), bx(i,k,j), bx_old1(i,k,j), bx_old2(i,k,j), &
               nl_bx     (i,k,j), &
               nl_bx_old1(i,k,j), &
               nl_bx_old2(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
            call gear3(dt, by_new(i,k,j), by(i,k,j), by_old1(i,k,j), by_old2(i,k,j), &
               nl_by     (i,k,j)  &
                 - q*shear_flg*bx     (i,k,j), &
               nl_by_old1(i,k,j)  &
                 - q*shear_flg*bx_old1(i,k,j), &
               nl_by_old2(i,k,j)  &
                 - q*shear_flg*bx_old2(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )
            call gear3(dt, bz_new(i,k,j), bz(i,k,j), bz_old1(i,k,j), bz_old2(i,k,j), &
               nl_bz     (i,k,j), &
               nl_bz_old1(i,k,j), &
               nl_bz_old2(i,k,j), &
               eta*(k2t(i,k,j)/k2_max)**eta_exp &
            )

          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      ux_old2 = ux_old1
      ux_old1 = ux
      nl_ux_old2 = nl_ux_old1
      nl_ux_old1 = nl_ux

      uy_old2 = uy_old1
      uy_old1 = uy
      nl_uy_old2 = nl_uy_old1
      nl_uy_old1 = nl_uy

      uz_old2 = uz_old1
      uz_old1 = uz
      nl_uz_old2 = nl_uz_old1
      nl_uz_old1 = nl_uz

      bx_old2 = bx_old1
      bx_old1 = bx
      nl_bx_old2 = nl_bx_old1
      nl_bx_old1 = nl_bx

      by_old2 = by_old1
      by_old1 = by
      nl_by_old2 = nl_by_old1
      nl_by_old1 = nl_by

      bz_old2 = bz_old1
      bz_old1 = bz
      nl_bz_old2 = nl_bz_old1
      nl_bz_old1 = nl_bz

      p_old2 =  p_old1
      p_old1 =  p
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fux_old2 = fux_old1
        fux_old1 = fux

        fuy_old2 = fuy_old1
        fuy_old1 = fuy

        fuz_old2 = fuz_old1
        fuz_old1 = fuz
        !$omp end workshare
      endif

      if(shear) then
        !$omp workshare
        kxt_old2 = kxt_old1
        kxt_old1 = kxt
        !$omp end workshare
      endif
    endif

    !$omp workshare
    ux = ux_new
    uy = uy_new
    uz = uz_new
    bx = bx_new
    by = by_new
    bz = bz_new
    !$omp end workshare

    ! Dealiasing
    if(trim(dealias_scheme) /= '2/3') then
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            ux(i,k,j) = ux_new(i,k,j)*filter(i,k,j)
            uy(i,k,j) = uy_new(i,k,j)*filter(i,k,j)
            uz(i,k,j) = uz_new(i,k,j)*filter(i,k,j)
            bx(i,k,j) = bx_new(i,k,j)*filter(i,k,j)
            by(i,k,j) = by_new(i,k,j)*filter(i,k,j)
            bz(i,k,j) = bz_new(i,k,j)*filter(i,k,j)
          enddo
        enddo
      enddo
    endif

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
    allocate(nbl2inv_div_u(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nbl2inv_div_u = 0.d0
    allocate(nbl2inv_div_b(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nbl2inv_div_b = 0.d0
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          nbl2inv_div_u(i,k,j) = -zi*(   kxt(i,j)*ux(i,k,j) &
                                       + ky(j)   *uy(i,k,j) &
                                       + kz(k)   *uz(i,k,j) )*k2t_inv(i,k,j)
          nbl2inv_div_b(i,k,j) = -zi*(   kxt(i,j)*bx(i,k,j) &
                                       + ky(j)   *by(i,k,j) &
                                       + kz(k)   *bz(i,k,j) )*k2t_inv(i,k,j)

          ux(i,k,j) = ux(i,k,j) - zi*kxt(i,j)*nbl2inv_div_u(i,k,j)
          uy(i,k,j) = uy(i,k,j) - zi*ky(j)   *nbl2inv_div_u(i,k,j)
          uz(i,k,j) = uz(i,k,j) - zi*kz(k)   *nbl2inv_div_u(i,k,j)
          bx(i,k,j) = bx(i,k,j) - zi*kxt(i,j)*nbl2inv_div_b(i,k,j)
          by(i,k,j) = by(i,k,j) - zi*ky(j)   *nbl2inv_div_b(i,k,j)
          bz(i,k,j) = bz(i,k,j) - zi*kz(k)   *nbl2inv_div_b(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do
    deallocate(nbl2inv_div_u)
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
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: q
    use shearing_box, only: shear_flg, tsc, k2t, k2t_inv
    use file, only: open_output_file
    implicit none
    integer :: i, j, k

    allocate(ux_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); ux_new     = 0.d0
    allocate(ux_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); ux_old2    = 0.d0

    allocate(uy_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); uy_new     = 0.d0
    allocate(uy_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); uy_old2    = 0.d0

    allocate(uz_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); uz_new     = 0.d0
    allocate(uz_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); uz_old2    = 0.d0

    allocate(bx_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bx_new     = 0.d0
    allocate(bx_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bx_old2    = 0.d0

    allocate(by_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); by_new     = 0.d0
    allocate(by_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); by_old2    = 0.d0

    allocate(bz_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bz_new     = 0.d0
    allocate(bz_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bz_old2    = 0.d0

    allocate( p_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en));  p_old2    = 0.d0

    allocate(nl_ux     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_ux      = 0.d0
    allocate(nl_ux_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_ux_old1 = 0.d0
    allocate(nl_ux_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_ux_old2 = 0.d0

    allocate(nl_uy     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_uy      = 0.d0
    allocate(nl_uy_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_uy_old1 = 0.d0
    allocate(nl_uy_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_uy_old2 = 0.d0

    allocate(nl_uz     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_uz      = 0.d0
    allocate(nl_uz_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_uz_old1 = 0.d0
    allocate(nl_uz_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_uz_old2 = 0.d0

    allocate(nl_bx     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bx      = 0.d0
    allocate(nl_bx_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bx_old1 = 0.d0
    allocate(nl_bx_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bx_old2 = 0.d0

    allocate(nl_by     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_by      = 0.d0
    allocate(nl_by_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_by_old1 = 0.d0
    allocate(nl_by_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_by_old2 = 0.d0

    allocate(nl_bz     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bz      = 0.d0
    allocate(nl_bz_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bz_old1 = 0.d0
    allocate(nl_bz_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bz_old2 = 0.d0

    allocate(flx       (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran)); flx = 0.d0

    allocate(fux_old2  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fux_old2   = 0.d0
    allocate(fuy_old2  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fuy_old2   = 0.d0
    allocate(fuz_old2  (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fuz_old2   = 0.d0

    allocate(kxt     (ikx_st:ikx_en, iky_st:iky_en))
    allocate(kxt_old1(ikx_st:ikx_en, iky_st:iky_en))
    allocate(kxt_old2(ikx_st:ikx_en, iky_st:iky_en))

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

    kxt_old1 = kxt
    kxt_old2 = kxt

    call open_output_file (max_vel_unit, 'max_vel.dat')

  end subroutine init_multistep_fields


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Compute pressure from initial u
!-----------------------------------------------!
  subroutine init_pressure
    use fields, only: p
    use params, only: zi, nonlinear
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use shearing_box, only: k2t_inv
    implicit none
    integer :: i, j, k

    if(nonlinear) call get_nonlinear_terms
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          ! get p
          p(i,k,j) = -zi*( kx(i)*nl_ux(i,k,j) &
                         + ky(j)*nl_uy(i,k,j) &
                         + kz(k)*nl_uz(i,k,j) )*k2t_inv(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do
  end subroutine init_pressure


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Calculate nonlinear terms via
!!          1. Calculate grad in Fourier space
!!          2. Inverse FFT
!!          3. Calculate poisson brackets 
!!             in real space
!!          4. Forward FFT
!-----------------------------------------------!
  subroutine get_nonlinear_terms
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    implicit none

    complex(r8), allocatable, dimension(:,:,:,:) :: wbk
    real   (r8), allocatable, dimension(:,:,:,:) :: wb , wf 

    integer :: i, j, k
    real   (r8) :: max_vel, dt_cfl, dt_digit

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    allocate(wbk(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nbtran)); wbk = 0.d0
    allocate(wb (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nbtran)); wb  = 0.d0
    allocate(wf (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nftran)); wf  = 0.d0

    ! 1. Calculate grad in Fourier space
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          wbk(i,k,j,iux) = ux(i,k,j)
          wbk(i,k,j,iuy) = uy(i,k,j)
          wbk(i,k,j,iuz) = uz(i,k,j)
          wbk(i,k,j,ibx) = bx(i,k,j)
          wbk(i,k,j,iby) = by(i,k,j)
          wbk(i,k,j,ibz) = bz(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 2. Inverse FFT
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_btran_c2r_many(wbk, nk_local_tot, wb, nl_local_tot, nbtran, 'tff')
    if (proc0) call put_time_stamp(timer_fft)

    ! (get max_vel for dt reset)
    max_vel = max( &
              maxval(abs(wb(:,:,:,iux)))*cflx, maxval(abs(wb(:,:,:,iuy))*cfly), maxval(abs(wb(:,:,:,iuz))*cflz) &
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
      dt = floor(dt_cfl*10.d0**(-dt_digit))*10.d0**dt_digit

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
    if(dt_cfl > increase_dt .and. dt < increase_dt) then
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

    ! 3. Calculate poisson brackets in real space
    !$omp parallel do private(j, k) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          wf(j,k,i,iflx_uxx) = wb(j,k,i,iux)*wb(j,k,i,iux) - wb(j,k,i,ibx)*wb(j,k,i,ibx)
          wf(j,k,i,iflx_uxy) = wb(j,k,i,iux)*wb(j,k,i,iuy) - wb(j,k,i,ibx)*wb(j,k,i,iby)
          wf(j,k,i,iflx_uxz) = wb(j,k,i,iux)*wb(j,k,i,iuz) - wb(j,k,i,ibx)*wb(j,k,i,ibz)
          wf(j,k,i,iflx_uyy) = wb(j,k,i,iuy)*wb(j,k,i,iuy) - wb(j,k,i,iby)*wb(j,k,i,iby)
          wf(j,k,i,iflx_uyz) = wb(j,k,i,iuy)*wb(j,k,i,iuz) - wb(j,k,i,iby)*wb(j,k,i,ibz)
          wf(j,k,i,iflx_uzz) = wb(j,k,i,iuz)*wb(j,k,i,iuz) - wb(j,k,i,ibz)*wb(j,k,i,ibz)

          wf(j,k,i,iflx_bx ) = wb(j,k,i,iby)*wb(j,k,i,iuz) - wb(j,k,i,ibz)*wb(j,k,i,iuy)
          wf(j,k,i,iflx_by ) = wb(j,k,i,ibz)*wb(j,k,i,iux) - wb(j,k,i,ibx)*wb(j,k,i,iuz)
          wf(j,k,i,iflx_bz ) = wb(j,k,i,ibx)*wb(j,k,i,iuy) - wb(j,k,i,iby)*wb(j,k,i,iux)
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
!! @date    15 Jul 2021
!! @brief   Output series modes
!-----------------------------------------------!
  subroutine output_series_modes
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use params, only: n_series_modes, series_modes
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0, sum_reduce
    use time, only: tt
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
    use params, only: dealias_scheme => dealias
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

    if(trim(dealias_scheme) /= '2/3') then
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
    endif

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

end module advance

