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
  complex(r8), allocatable, dimension(:,:,:)   :: rho_new, rho_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  mx_new,  mx_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  my_new,  my_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  mz_new,  mz_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  bx_new,  bx_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  by_new,  by_old2
  complex(r8), allocatable, dimension(:,:,:)   ::  bz_new,  bz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: sgm_new, sgm_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: flx
  complex(r8), allocatable, dimension(:,:,:)   :: nl_rho, nl_rho_old1, nl_rho_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_mx , nl_mx_old1 , nl_mx_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_my , nl_my_old1 , nl_my_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_mz , nl_mz_old1 , nl_mz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_bx , nl_bx_old1 , nl_bx_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_by , nl_by_old1 , nl_by_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_bz , nl_bz_old1 , nl_bz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nl_sgm, nl_sgm_old1, nl_sgm_old2
  complex(r8), allocatable, dimension(:,:,:)   :: fmx_old2, fmy_old2, fmz_old2
  complex(r8), allocatable, dimension(:,:,:)   :: nbl2inv_div_b ! nabla^-2 (div b)
  real   (r8), allocatable, dimension(:,:)     :: kxt, kxt_old1, kxt_old2
  real   (r8) :: cflx, cfly, cflz
  integer :: max_vel_unit

  ! Backward FFT variables
  integer, parameter :: nbtran = 8
  integer, parameter :: irho = 1
  integer, parameter :: imx  = 2, imy = 3, imz = 4
  integer, parameter :: ibx  = 5, iby = 6, ibz = 7
  integer, parameter :: isgm = 8
  ! Forward FFT variables
  integer, parameter :: nftran = 12
  integer, parameter :: iflx_mxx  = 1 , iflx_mxy  = 2 , iflx_mxz  = 3  ! rho*uu - bb + [(cs^2/va^2)*p + b^2/2]*I (tensor) 
  integer, parameter ::                 iflx_myy  = 4 , iflx_myz  = 5  ! rho*uu - bb + [(cs^2/va^2)*p + b^2/2]*I (tensor) 
  integer, parameter ::                                 iflx_mzz  = 6  ! rho*uu - bb + [(cs^2/va^2)*p + b^2/2]*I (tensor) 
  integer, parameter :: iflx_bx   = 7 , iflx_by   = 8 , iflx_bz   = 9  ! b x u
  integer, parameter :: iflx_sgmx = 10, iflx_sgmy = 11, iflx_sgmz = 12 ! sgm*u

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
    use fields, only: rho, sgm
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use fields, only: rho_old1, sgm_old1
    use fields, only: mx_old1, my_old1, mz_old1
    use fields, only: bx_old1, by_old1, bz_old1
    use grid, only: k2, k2inv, k2_max
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, cfl, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: nu, eta, nu_exp, eta_exp, zi, nonlinear, shear, q
    use shearing_box, only: shear_flg, k2_old1, tsc, tremap
    use force, only: fmx, fmz, fmy, fmx_old1, fmy_old1, fmz_old1, driven, update_force, get_force
    use advance_common, only: gear1, gear2, gear3
    implicit none
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

    ! Calcualte nonlinear terms
    if(nonlinear) call get_nonlinear_terms

    ! Calcualte force terms
    if (driven) then
      call update_force
      call get_force('mx', fmx)
      call get_force('my', fmy)
      call get_force('mz', fmz)
    endif

    ! 1st order 
    if(counter == 1) then
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            nl_rho(i,k,j) = -zi*( kxt(i,j)*mx(i,k,j) &
                                + ky(j)   *my(i,k,j) &
                                + kz(k)   *mz(i,k,j) )

            ! div {rho*uu - bb + [(cs^2/va^2)*rho + b^2/2]*I}
            nl_mx(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxx) &
                                 + ky(j)   *flx(i,k,j,iflx_mxy) &
                                 + kz(k)   *flx(i,k,j,iflx_mxz) )
                                                          
            nl_my(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxy) &
                                 + ky(j)   *flx(i,k,j,iflx_myy) &
                                 + kz(k)   *flx(i,k,j,iflx_myz) )
                                                          
            nl_mz(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxz) &
                                 + ky(j)   *flx(i,k,j,iflx_myz) &
                                 + kz(k)   *flx(i,k,j,iflx_mzz) )

            ! curl (b x u)
            nl_bx(i,k,j)  = -zi*(  ky(j)   *flx(i,k,j,iflx_bz ) &
                                 - kz(k)   *flx(i,k,j,iflx_by ) )
                                                          
            nl_by(i,k,j)  = -zi*(  kz(k)   *flx(i,k,j,iflx_bx ) &
                                 - kxt(i,j)*flx(i,k,j,iflx_bz ) )
                                                          
            nl_bz(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_by ) &
                                 - ky(j)   *flx(i,k,j,iflx_bx ) )
            ! div (sigma*u)
            nl_sgm(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_sgmx) &
                                 + ky(j)   *flx(i,k,j,iflx_sgmy) &
                                 + kz(k)   *flx(i,k,j,iflx_sgmz) )


            ! update rho
            call gear1(dt, rho_new(i,k,j), rho(i,k,j), &
               nl_rho(i,k,j), &
               0*nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

            ! update m
            call gear1(dt, mx_new(i,k,j), mx(i,k,j), &
               nl_mx(i,k,j) + fmx(i,k,j) + 2.d0*shear_flg*my(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )
            call gear1(dt, my_new(i,k,j), my(i,k,j), &
               nl_my(i,k,j) + fmy(i,k,j) - (2.d0 - q)*shear_flg*mx(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )
            call gear1(dt, mz_new(i,k,j), mz(i,k,j), &
               nl_mz(i,k,j) + fmz(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

            ! update b
            call gear1(dt, bx_new(i,k,j), bx(i,k,j), &
               nl_bx(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )
            call gear1(dt, by_new(i,k,j), by(i,k,j), &
               nl_by(i,k,j) - q*shear_flg*bx(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )
            call gear1(dt, bz_new(i,k,j), bz(i,k,j), &
               nl_bz(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )

            ! update sgm
            call gear1(dt, sgm_new(i,k,j), sgm(i,k,j), &
               nl_sgm(i,k,j), &
               0*nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      rho_old1 = rho
      nl_rho_old1 = nl_rho

      mx_old1 = mx
      nl_mx_old1 = nl_mx

      my_old1 = my
      nl_my_old1 = nl_my

      mz_old1 = mz
      nl_mz_old1 = nl_mz

      bx_old1 = bx
      nl_bx_old1 = nl_bx

      by_old1 = by
      nl_by_old1 = nl_by

      bz_old1 = bz
      nl_bz_old1 = nl_bz

      sgm_old1 = sgm
      nl_sgm_old1 = nl_sgm
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fmx_old1 = fmx
        fmy_old1 = fmy
        fmz_old1 = fmz
        !$omp end workshare
      endif

      if(shear) then
        !$omp workshare
        k2_old1 = k2
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
            nl_rho(i,k,j) = -zi*( kxt(i,j)*mx(i,k,j) &
                                + ky(j)   *my(i,k,j) &
                                + kz(k)   *mz(i,k,j) )

            ! div {rho*uu - bb + [(cs^2/va^2)*rho + b^2/2]*I}
            nl_mx(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxx) &
                                 + ky(j)   *flx(i,k,j,iflx_mxy) &
                                 + kz(k)   *flx(i,k,j,iflx_mxz) )
                                                          
            nl_my(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxy) &
                                 + ky(j)   *flx(i,k,j,iflx_myy) &
                                 + kz(k)   *flx(i,k,j,iflx_myz) )
                                                          
            nl_mz(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxz) &
                                 + ky(j)   *flx(i,k,j,iflx_myz) &
                                 + kz(k)   *flx(i,k,j,iflx_mzz) )

            ! curl (b x u)
            nl_bx(i,k,j)  = -zi*(  ky(j)   *flx(i,k,j,iflx_bz ) &
                                 - kz(k)   *flx(i,k,j,iflx_by ) )
                                                          
            nl_by(i,k,j)  = -zi*(  kz(k)   *flx(i,k,j,iflx_bx ) &
                                 - kxt(i,j)*flx(i,k,j,iflx_bz ) )
                                                          
            nl_bz(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_by ) &
                                 - ky(j)   *flx(i,k,j,iflx_bx ) )
            ! div (sigma*u)
            nl_sgm(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_sgmx) &
                                 + ky(j)   *flx(i,k,j,iflx_sgmy) &
                                 + kz(k)   *flx(i,k,j,iflx_sgmz) )


            ! update rho
            call gear2(dt, rho_new(i,k,j), rho(i,k,j), rho_old1(i,k,j), &
               nl_rho     (i,k,j), &
               nl_rho_old1(i,k,j), &
               0*nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

            ! update m
            call gear2(dt, mx_new(i,k,j), mx(i,k,j), mx_old1(i,k,j), &
               nl_mx     (i,k,j) + fmx     (i,k,j) + 2.d0*shear_flg*my     (i,k,j), &
               nl_mx_old1(i,k,j) + fmx_old1(i,k,j) + 2.d0*shear_flg*my_old1(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )
            call gear2(dt, my_new(i,k,j), my(i,k,j), my_old1(i,k,j), &
               nl_my     (i,k,j) + fmy     (i,k,j) - (2.d0 - q)*shear_flg*mx     (i,k,j), &
               nl_my_old1(i,k,j) + fmy_old1(i,k,j) - (2.d0 - q)*shear_flg*mx_old1(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )
            call gear2(dt, mz_new(i,k,j), mz(i,k,j), mz_old1(i,k,j), &
               nl_mz     (i,k,j) + fmz     (i,k,j), &
               nl_mz_old1(i,k,j) + fmz_old1(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

            ! update b
            call gear2(dt, bx_new(i,k,j), bx(i,k,j), bx_old1(i,k,j), &
               nl_bx     (i,k,j), &
               nl_bx_old1(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )
            call gear2(dt, by_new(i,k,j), by(i,k,j), by_old1(i,k,j), &
               nl_by     (i,k,j) - q*shear_flg*bx     (i,k,j), &
               nl_by_old1(i,k,j) - q*shear_flg*bx_old1(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )
            call gear2(dt, bz_new(i,k,j), bz(i,k,j), bz_old1(i,k,j), &
               nl_bz     (i,k,j), &
               nl_bz_old1(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )

            ! update sgm
            call gear2(dt, sgm_new(i,k,j), sgm(i,k,j), sgm_old1(i,k,j), &
               nl_sgm     (i,k,j), &
               nl_sgm_old1(i,k,j), &
               0*nu*(k2(i,k,j)/k2_max)**nu_exp &
            )
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      rho_old2 = rho_old1
      rho_old1 = rho
      nl_rho_old2 = nl_rho_old1
      nl_rho_old1 = nl_rho

      mx_old2 = mx_old1
      mx_old1 = mx
      nl_mx_old2 = nl_mx_old1
      nl_mx_old1 = nl_mx

      my_old2 = my_old1
      my_old1 = my
      nl_my_old2 = nl_my_old1
      nl_my_old1 = nl_my

      mz_old2 = mz_old1
      mz_old1 = mz
      nl_mz_old2 = nl_mz_old1
      nl_mz_old1 = nl_mz

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

      sgm_old2 = sgm_old1
      sgm_old1 = sgm
      nl_sgm_old2 = nl_sgm_old1
      nl_sgm_old1 = nl_sgm
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fmx_old2 = fmx_old1
        fmx_old1 = fmx

        fmy_old2 = fmy_old1
        fmy_old1 = fmy

        fmz_old2 = fmz_old1
        fmz_old1 = fmz
        !$omp end workshare
      endif

      if(shear) then
        !$omp workshare
        k2_old1 = k2

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
            nl_rho(i,k,j) = -zi*( kxt(i,j)*mx(i,k,j) &
                                + ky(j)   *my(i,k,j) &
                                + kz(k)   *mz(i,k,j) )

            ! div {rho*uu - bb + [(cs^2/va^2)*rho + b^2/2]*I}
            nl_mx(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxx) &
                                 + ky(j)   *flx(i,k,j,iflx_mxy) &
                                 + kz(k)   *flx(i,k,j,iflx_mxz) )
                                                          
            nl_my(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxy) &
                                 + ky(j)   *flx(i,k,j,iflx_myy) &
                                 + kz(k)   *flx(i,k,j,iflx_myz) )
                                                          
            nl_mz(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_mxz) &
                                 + ky(j)   *flx(i,k,j,iflx_myz) &
                                 + kz(k)   *flx(i,k,j,iflx_mzz) )

            ! curl (b x u)
            nl_bx(i,k,j)  = -zi*(  ky(j)   *flx(i,k,j,iflx_bz ) &
                                 - kz(k)   *flx(i,k,j,iflx_by ) )
                                                          
            nl_by(i,k,j)  = -zi*(  kz(k)   *flx(i,k,j,iflx_bx ) &
                                 - kxt(i,j)*flx(i,k,j,iflx_bz ) )
                                                          
            nl_bz(i,k,j)  = -zi*(  kxt(i,j)*flx(i,k,j,iflx_by ) &
                                 - ky(j)   *flx(i,k,j,iflx_bx ) )
            ! div (sigma*u)
            nl_sgm(i,k,j) = -zi*(  kxt(i,j)*flx(i,k,j,iflx_sgmx) &
                                 + ky(j)   *flx(i,k,j,iflx_sgmy) &
                                 + kz(k)   *flx(i,k,j,iflx_sgmz) )


            ! update rho
            call gear3(dt, rho_new(i,k,j), rho(i,k,j), rho_old1(i,k,j), rho_old2(i,k,j), &
               nl_rho     (i,k,j), &
               nl_rho_old1(i,k,j), &
               nl_rho_old2(i,k,j), &
               0*nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

            ! update m
            call gear3(dt, mx_new(i,k,j), mx(i,k,j), mx_old1(i,k,j), mx_old2(i,k,j), &
               nl_mx     (i,k,j) + fmx     (i,k,j) + 2.d0*shear_flg*my     (i,k,j), &
               nl_mx_old1(i,k,j) + fmx_old1(i,k,j) + 2.d0*shear_flg*my_old1(i,k,j), &
               nl_mx_old2(i,k,j) + fmx_old2(i,k,j) + 2.d0*shear_flg*my_old2(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )
            call gear3(dt, my_new(i,k,j), my(i,k,j), my_old1(i,k,j), my_old2(i,k,j), &
               nl_my     (i,k,j) + fmy     (i,k,j) - (2.d0 - q)*shear_flg*mx     (i,k,j), &
               nl_my_old1(i,k,j) + fmy_old1(i,k,j) - (2.d0 - q)*shear_flg*mx_old1(i,k,j), &
               nl_my_old2(i,k,j) + fmy_old2(i,k,j) - (2.d0 - q)*shear_flg*mx_old2(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )
            call gear3(dt, mz_new(i,k,j), mz(i,k,j), mz_old1(i,k,j), mz_old2(i,k,j), &
               nl_mz     (i,k,j) + fmz     (i,k,j), &
               nl_mz_old1(i,k,j) + fmz_old1(i,k,j), &
               nl_mz_old2(i,k,j) + fmz_old2(i,k,j), &
               nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

            ! update b
            call gear3(dt, bx_new(i,k,j), bx(i,k,j), bx_old1(i,k,j), bx_old2(i,k,j), &
               nl_bx     (i,k,j), &
               nl_bx_old1(i,k,j), &
               nl_bx_old2(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )
            call gear3(dt, by_new(i,k,j), by(i,k,j), by_old1(i,k,j), by_old2(i,k,j), &
               nl_by     (i,k,j) - q*shear_flg*bx     (i,k,j), &
               nl_by_old1(i,k,j) - q*shear_flg*bx_old1(i,k,j), &
               nl_by_old2(i,k,j) - q*shear_flg*bx_old2(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )
            call gear3(dt, bz_new(i,k,j), bz(i,k,j), bz_old1(i,k,j), bz_old2(i,k,j), &
               nl_bz     (i,k,j), &
               nl_bz_old1(i,k,j), &
               nl_bz_old2(i,k,j), &
               eta*(k2(i,k,j)/k2_max)**eta_exp &
            )

            ! update sgm
            call gear3(dt, sgm_new(i,k,j), sgm(i,k,j), sgm_old1(i,k,j), sgm_old2(i,k,j), &
               nl_sgm     (i,k,j), &
               nl_sgm_old1(i,k,j), &
               nl_sgm_old2(i,k,j), &
               0*nu*(k2(i,k,j)/k2_max)**nu_exp &
            )

          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      rho_old2 = rho_old1
      rho_old1 = rho
      nl_rho_old2 = nl_rho_old1
      nl_rho_old1 = nl_rho

      mx_old2 = mx_old1
      mx_old1 = mx
      nl_mx_old2 = nl_mx_old1
      nl_mx_old1 = nl_mx

      my_old2 = my_old1
      my_old1 = my
      nl_my_old2 = nl_my_old1
      nl_my_old1 = nl_my

      mz_old2 = mz_old1
      mz_old1 = mz
      nl_mz_old2 = nl_mz_old1
      nl_mz_old1 = nl_mz

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

      sgm_old2 = sgm_old1
      sgm_old1 = sgm
      nl_sgm_old2 = nl_sgm_old1
      nl_sgm_old1 = nl_sgm
      !$omp end workshare

      if (driven) then
        !$omp workshare
        fmx_old2 = fmx_old1
        fmx_old1 = fmx

        fmy_old2 = fmy_old1
        fmy_old1 = fmy

        fmz_old2 = fmz_old1
        fmz_old1 = fmz
        !$omp end workshare
      endif

      if(shear) then
        !$omp workshare
        k2_old1 = k2

        kxt_old2 = kxt_old1
        kxt_old1 = kxt
        !$omp end workshare
      endif
    endif

    !$omp workshare
    rho = rho_new
     mx =  mx_new
     my =  my_new
     mz =  mz_new
     bx =  bx_new
     by =  by_new
     bz =  bz_new
    sgm = sgm_new
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
            k2(i,k,j) = kxt(i,j)**2 + ky(j)**2 + kz(k)**2

            if(k2(i,k,j) == 0.d0) then
              k2inv(i,k,j) = 0.d0
            else
              k2inv(i,k,j) = 1.d0/k2(i,k,j)
            endif
          enddo
        enddo
      enddo
      !$omp end parallel do
    endif

    ! Div u & b cleaing
    allocate(nbl2inv_div_b(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nbl2inv_div_b = 0.d0
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          nbl2inv_div_b(i,k,j) = -zi*(   kxt(i,j)*bx(i,k,j) &
                                       + ky(j)   *by(i,k,j) &
                                       + kz(k)   *bz(i,k,j) )*k2inv(i,k,j)

          bx(i,k,j) = bx(i,k,j) - zi*kxt(i,j)*nbl2inv_div_b(i,k,j)
          by(i,k,j) = by(i,k,j) - zi*ky(j)   *nbl2inv_div_b(i,k,j)
          bz(i,k,j) = bz(i,k,j) - zi*kz(k)   *nbl2inv_div_b(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do
    deallocate(nbl2inv_div_b)

    if (proc0) call put_time_stamp(timer_advance)
  end subroutine solve


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Allocate tmp fields for a multi 
!!          timestep method
!-----------------------------------------------!
  subroutine init_multistep_fields
    use grid, only: kx, ky, kz, k2, k2inv
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: q
    use shearing_box, only: shear_flg, tsc
    use file, only: open_output_file
    implicit none
    integer :: i, j, k

    allocate(rho_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); rho_new     = 0.d0
    allocate(rho_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); rho_old2    = 0.d0

    allocate(mx_new     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); mx_new      = 0.d0
    allocate(mx_old2    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); mx_old2     = 0.d0
                                                                                   
    allocate(my_new     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); my_new      = 0.d0
    allocate(my_old2    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); my_old2     = 0.d0
                                                                                   
    allocate(mz_new     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); mz_new      = 0.d0
    allocate(mz_old2    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); mz_old2     = 0.d0
                                                                                   
    allocate(bx_new     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bx_new      = 0.d0
    allocate(bx_old2    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bx_old2     = 0.d0
                                                                                   
    allocate(by_new     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); by_new      = 0.d0
    allocate(by_old2    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); by_old2     = 0.d0
                                                                                   
    allocate(bz_new     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bz_new      = 0.d0
    allocate(bz_old2    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); bz_old2     = 0.d0

    allocate(sgm_new    (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); sgm_new     = 0.d0
    allocate(sgm_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); sgm_old2    = 0.d0
                                                                                   
    allocate(nl_rho     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_rho      = 0.d0
    allocate(nl_rho_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_rho_old1 = 0.d0
    allocate(nl_rho_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_rho_old2 = 0.d0
                                                                                   
    allocate(nl_mx      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_mx       = 0.d0
    allocate(nl_mx_old1 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_mx_old1  = 0.d0
    allocate(nl_mx_old2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_mx_old2  = 0.d0
                                                                                   
    allocate(nl_my      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_my       = 0.d0
    allocate(nl_my_old1 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_my_old1  = 0.d0
    allocate(nl_my_old2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_my_old2  = 0.d0
                                                                                   
    allocate(nl_mz      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_mz       = 0.d0
    allocate(nl_mz_old1 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_mz_old1  = 0.d0
    allocate(nl_mz_old2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_mz_old2  = 0.d0
                                                                                   
    allocate(nl_bx      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bx       = 0.d0
    allocate(nl_bx_old1 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bx_old1  = 0.d0
    allocate(nl_bx_old2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bx_old2  = 0.d0
                                                                                   
    allocate(nl_by      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_by       = 0.d0
    allocate(nl_by_old1 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_by_old1  = 0.d0
    allocate(nl_by_old2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_by_old2  = 0.d0
                                                                                   
    allocate(nl_bz      (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bz       = 0.d0
    allocate(nl_bz_old1 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bz_old1  = 0.d0
    allocate(nl_bz_old2 (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_bz_old2  = 0.d0
                                                                                   
    allocate(nl_sgm     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_sgm      = 0.d0
    allocate(nl_sgm_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_sgm_old1 = 0.d0
    allocate(nl_sgm_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); nl_sgm_old2 = 0.d0

    allocate(flx        (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran)); flx = 0.d0

    allocate(fmx_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmx_old2    = 0.d0
    allocate(fmy_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmy_old2    = 0.d0
    allocate(fmz_old2   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en)); fmz_old2    = 0.d0

    allocate(kxt     (ikx_st:ikx_en, iky_st:iky_en))
    allocate(kxt_old1(ikx_st:ikx_en, iky_st:iky_en))
    allocate(kxt_old2(ikx_st:ikx_en, iky_st:iky_en))

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do i = ikx_st, ikx_en
        kxt(i, j) = kx(i) + q*shear_flg*tsc*ky(j)
        do k = ikz_st, ikz_en
          k2(i, k, j) = kxt(i, j)**2 + ky(j)**2 + kz(k)**2
          if(k2(i, k, j) == 0.d0) then
            k2inv(i, k, j) = 0.d0
          else
            k2inv(i, k, j) = 1.0d0/k2(i, k, j)
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
!! @brief   Calculate nonlinear terms via
!!          1. Calculate grad in Fourier space
!!          2. Inverse FFT
!!          3. Calculate poisson brackets 
!!             in real space
!!          4. Forward FFT
!-----------------------------------------------!
  subroutine get_nonlinear_terms
    use fields, only: rho, sgm
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi, gamma, cs2va2
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    implicit none

    complex(r8), allocatable, dimension(:,:,:,:) :: wbk
    real   (r8), allocatable, dimension(:,:,:,:) :: wb , wf 
    real   (r8), allocatable, dimension(:,:,:)   :: ux_r, uy_r, uz_r, p_r
    real   (r8), allocatable, dimension(:,:,:) :: cfx2_r, cfy2_r, cfz2_r ! Fast mode phase speed
    real   (r8) :: vax2_r, vay2_r, vaz2_r ! Alfven speed

    integer :: i, j, k
    real   (r8) :: max_vel, dt_cfl, dt_digit

    if (proc0) call put_time_stamp(timer_nonlinear_terms)

    allocate(wbk   (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nbtran)); wbk    = 0.d0
    allocate(wb    (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nbtran)); wb     = 0.d0
    allocate(wf    (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en, nftran)); wf     = 0.d0
    allocate(ux_r  (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))        ; ux_r   = 0.d0
    allocate(uy_r  (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))        ; uy_r   = 0.d0
    allocate(uz_r  (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))        ; uz_r   = 0.d0
    allocate(p_r   (ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))        ; p_r    = 0.d0
    allocate(cfx2_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))        ; cfx2_r = 0.d0
    allocate(cfy2_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))        ; cfy2_r = 0.d0
    allocate(cfz2_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en))        ; cfx2_r = 0.d0

    ! 1. Calculate grad in Fourier space
    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          wbk(i,k,j,irho) = rho(i,k,j)
          wbk(i,k,j,imx ) = mx (i,k,j)
          wbk(i,k,j,imy ) = my (i,k,j)
          wbk(i,k,j,imz ) = mz (i,k,j)
          wbk(i,k,j,ibx ) = bx (i,k,j)
          wbk(i,k,j,iby ) = by (i,k,j)
          wbk(i,k,j,ibz ) = bz (i,k,j)
          wbk(i,k,j,isgm) = sgm(i,k,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 2. Inverse FFT
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_btran_c2r_many(wbk, nk_local_tot, wb, nl_local_tot, nbtran, 'tff')
    if (proc0) call put_time_stamp(timer_fft)

    !$omp workshare
    ux_r = wb(:,:,:,imx)/wb(:,:,:,irho)
    uy_r = wb(:,:,:,imy)/wb(:,:,:,irho)
    uz_r = wb(:,:,:,imz)/wb(:,:,:,irho)
    p_r  = wb(:,:,:,isgm)/wb(:,:,:,irho)**(-gamma + 1.d0)
    !$omp end workshare

    ! (get max_vel for dt reset)
    !$omp parallel do private(j, k) schedule(static)
    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          vax2_r = wb(j,k,i,ibx)**2/wb(j,k,i,irho)
          vay2_r = wb(j,k,i,iby)**2/wb(j,k,i,irho)
          vaz2_r = wb(j,k,i,ibz)**2/wb(j,k,i,irho)

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
    max_vel = max( &
              maxval(abs(ux_r + sqrt(cfx2_r)))*cflx, &
              maxval(abs(ux_r - sqrt(cfx2_r)))*cflx, &
              maxval(abs(uy_r + sqrt(cfy2_r)))*cfly, &
              maxval(abs(uy_r - sqrt(cfy2_r)))*cfly, &
              maxval(abs(uz_r + sqrt(cfz2_r)))*cflz, &
              maxval(abs(uz_r - sqrt(cfz2_r)))*cflz  &
            )
    call max_allreduce(max_vel)
    dt_cfl = 1.d0/max_vel
    if(proc0) write (unit=max_vel_unit, fmt="(100es30.21)") tt, max_vel

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
    if(dt_cfl > 10.0d0*dt .and. increase_dt) then
      if(proc0) then
        print *
        write (*, '("dt is increased from ", es12.4e3)', advance='no') dt
      endif

      dt = 5.0d0*dt

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
          wf(j,k,i,iflx_mxx) = wb(j,k,i,irho)*ux_r(j,k,i)*ux_r(j,k,i) - wb(j,k,i,ibx)*wb(j,k,i,ibx) &
                               + 0.5d0*(wb(j,k,i,ibx)**2 + wb(j,k,i,iby)**2 + wb(j,k,i,ibz)**2) &
                               + cs2va2/gamma*p_r(j,k,i)
          wf(j,k,i,iflx_mxy) = wb(j,k,i,irho)*ux_r(j,k,i)*uy_r(j,k,i) - wb(j,k,i,ibx)*wb(j,k,i,iby)
          wf(j,k,i,iflx_mxz) = wb(j,k,i,irho)*ux_r(j,k,i)*uz_r(j,k,i) - wb(j,k,i,ibx)*wb(j,k,i,ibz)
          wf(j,k,i,iflx_myy) = wb(j,k,i,irho)*uy_r(j,k,i)*uy_r(j,k,i) - wb(j,k,i,iby)*wb(j,k,i,iby) &
                               + 0.5d0*(wb(j,k,i,ibx)**2 + wb(j,k,i,iby)**2 + wb(j,k,i,ibz)**2) &
                               + cs2va2/gamma*p_r(j,k,i)
          wf(j,k,i,iflx_myz) = wb(j,k,i,irho)*uy_r(j,k,i)*uz_r(j,k,i) - wb(j,k,i,iby)*wb(j,k,i,ibz)
          wf(j,k,i,iflx_mzz) = wb(j,k,i,irho)*uz_r(j,k,i)*uz_r(j,k,i) - wb(j,k,i,ibz)*wb(j,k,i,ibz) &
                               + 0.5d0*(wb(j,k,i,ibx)**2 + wb(j,k,i,iby)**2 + wb(j,k,i,ibz)**2) &
                               + cs2va2/gamma*p_r(j,k,i)

          wf(j,k,i,iflx_bx ) = wb(j,k,i,iby)*uz_r(j,k,i) - wb(j,k,i,ibz)*uy_r(j,k,i)
          wf(j,k,i,iflx_by ) = wb(j,k,i,ibz)*ux_r(j,k,i) - wb(j,k,i,ibx)*uz_r(j,k,i)
          wf(j,k,i,iflx_bz ) = wb(j,k,i,ibx)*uy_r(j,k,i) - wb(j,k,i,iby)*ux_r(j,k,i)

          wf(j,k,i,iflx_sgmx) = wb(j,k,i,isgm)*ux_r(j,k,i)
          wf(j,k,i,iflx_sgmy) = wb(j,k,i,isgm)*uy_r(j,k,i)
          wf(j,k,i,iflx_sgmz) = wb(j,k,i,isgm)*uz_r(j,k,i)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 4. Forward FFT
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_ftran_r2c_many(wf, nl_local_tot, flx, nk_local_tot, nftran, 'fft')
    if (proc0) call put_time_stamp(timer_fft)
    flx = flx/nlx/nly/nlz

    deallocate(wbk)
    deallocate(wb )
    deallocate(wf )
    deallocate(ux_r)
    deallocate(uy_r)
    deallocate(uz_r)
    deallocate(p_r)
    deallocate(cfx2_r)
    deallocate(cfy2_r)
    deallocate(cfz2_r)

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms


!-----------------------------------------------!
!> @author  YK
!! @date    3 Mar 2021
!! @brief   Remap
!-----------------------------------------------!
  subroutine remap
    use fields, only: rho
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use grid, only: kx, ky, kz, k2, k2inv
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use time_stamp, only: put_time_stamp
    use shearing_box, only: tsc, nremap
    use shearing_box, only: timer_remap, remap_fld
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

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          kxt(i,j) = kx(i)
          k2(i, k, j) = kx(i)**2 + ky(j)**2 + kz(k)**2

          if(k2(i, k, j) == 0.d0) then
            k2inv(i, k, j) = 0.d0
          else
            k2inv(i, k, j) = 1.d0/k2(i, k, j)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

    if (proc0) call put_time_stamp(timer_remap)
  end subroutine remap

end module advance


