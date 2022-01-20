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
  implicit none

  public solve

  integer :: counter = 0
  complex(r8), allocatable, dimension(:,:,:)   :: phi_new, phi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: omg_new, omg_old2
  complex(r8), allocatable, dimension(:,:,:)   :: psi_new, psi_old2
  complex(r8), allocatable, dimension(:,:,:)   :: upa_new, upa_old2
  complex(r8), allocatable, dimension(:,:,:)   :: bpa_new, bpa_old2
  complex(r8), allocatable, dimension(:,:,:,:) :: nonlin, nonlin_old1, nonlin_old2
  real   (r8) :: cflx, cfly
  integer :: max_vel_unit

  ! Backward FFT variables
  integer, parameter :: nbtran = 12
  integer, parameter :: idphi_dx = 1 , idphi_dy = 2
  integer, parameter :: idomg_dx = 3 , idomg_dy = 4
  integer, parameter :: idpsi_dx = 5 , idpsi_dy = 6
  integer, parameter :: idjpa_dx = 7 , idjpa_dy = 8
  integer, parameter :: idupa_dx = 9 , idupa_dy = 10
  integer, parameter :: idbpa_dx = 11, idbpa_dy = 12
  ! Forward FFT variables
  integer, parameter :: nftran = 4
  integer, parameter :: inonlin_omg = 1, inonlin_psi = 2, inonlin_upa = 3, inonlin_bpa = 4

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
    use fields, only: phi, omg, psi, upa, bpa
    use fields, only: phi_old1, omg_old1, psi_old1, upa_old1, bpa_old1
    use grid, only: kprp2, kprp2inv, kz2, kprp2_max, kz2_max
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: dt, cfl, tt
    use time_stamp, only: put_time_stamp, timer_advance
    use mp, only: proc0
    use params, only: dealias_scheme => dealias
    use dealias, only: filter
    use params, only: va2cs2_plus_1, nupe_x , nupe_x_exp , nupe_z , nupe_z_exp , &
                                     nupa_x , nupa_x_exp , nupa_z , nupa_z_exp , &
                                     etape_x, etape_x_exp, etape_z, etape_z_exp, &
                                     etapa_x, etapa_x_exp, etapa_z, etapa_z_exp, &
                      zi, nonlinear, q
    use advance_common, only: gear1, gear2, gear3
    implicit none
    integer :: i, j, k

    if (proc0) call put_time_stamp(timer_advance)

    ! initialize tmp fields
    if(counter == 0) then
      call init_multistep_fields

      cflx = maxval(abs(kx))/cfl
      cfly = maxval(abs(ky))/cfl
      counter = 1
    endif

    ! Calcualte nonlinear terms
    if(nonlinear) call get_nonlinear_terms

    ! 1st order 
    if(counter == 1) then
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            call gear1(dt, omg_new(i, k, j), omg(i, k, j), &
               nonlin(i, k, j, inonlin_omg) - zi*kz(k)*kprp2(i, k, j)*psi(i, k, j) &
                 - 2.d0*zi*ky(j)*upa(i, k, j), &
               nupe_x*(kprp2(i, k, j)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
            )
            call gear1(dt, psi_new(i, k, j), psi(i, k, j), &
               nonlin(i, k, j, inonlin_psi) + zi*kz(k)*phi(i, k, j), &
               etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
            )
            call gear1(dt, upa_new(i, k, j), upa(i, k, j), &
               nonlin(i, k, j, inonlin_upa) + zi*kz(k)*bpa(i, k, j) &
                 + (2.d0 - q)*zi*ky(j)*phi(i, k, j), &
               nupa_x*(kprp2(i, k, j)/kprp2_max)**nupa_x_exp + nupa_z*(kz2(k)/kz2_max)**nupa_z_exp &
            )
            call gear1(dt, bpa_new(i, k, j), bpa(i, k, j), &
               (nonlin(i, k, j, inonlin_bpa) + zi*kz(k)*upa(i, k, j) &
                 + q*zi*ky(j)*psi(i, k, j) )/va2cs2_plus_1, &
               (etapa_x*(kprp2(i, k, j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp)/va2cs2_plus_1 &
            )

            phi_new(i, k, j) = omg_new(i, k, j)*(-kprp2inv(i, k, j))
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      phi_old1 = phi
      omg_old1 = omg
      psi_old1 = psi
      upa_old1 = upa
      bpa_old1 = bpa

      nonlin_old1 = nonlin
      !$omp end workshare

      counter = counter + 1
    ! 2nd order 
    elseif(counter == 2) then
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            call gear2(dt, omg_new(i, k, j), omg(i, k, j), omg_old1(i, k, j), &
               nonlin     (i, k, j, inonlin_omg) - zi*kz(k)*kprp2(i, k, j)*psi     (i, k, j) &
                 - 2.d0*zi*ky(j)*upa(i, k, j), &
               nonlin_old1(i, k, j, inonlin_omg) - zi*kz(k)*kprp2(i, k, j)*psi_old1(i, k, j) &
                 - 2.d0*zi*ky(j)*upa_old1(i, k, j), &
               nupe_x*(kprp2(i, k, j)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
            )
            call gear2(dt, psi_new(i, k, j), psi(i, k, j), psi_old1(i, k, j), &
               nonlin     (i, k, j, inonlin_psi) + zi*kz(k)*phi     (i, k, j), &
               nonlin_old1(i, k, j, inonlin_psi) + zi*kz(k)*phi_old1(i, k, j), &
               etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
            )
            call gear2(dt, upa_new(i, k, j), upa(i, k, j), upa_old1(i, k, j), &
               nonlin     (i, k, j, inonlin_upa) + zi*kz(k)*bpa     (i, k, j) &
                 + (2.d0 - q)*zi*ky(j)*phi     (i, k, j), &
               nonlin_old1(i, k, j, inonlin_upa) + zi*kz(k)*bpa_old1(i, k, j) &
                 + (2.d0 - q)*zi*ky(j)*phi_old1(i, k, j), &
               nupa_x*(kprp2(i, k, j)/kprp2_max)**nupa_x_exp + nupa_z*(kz2(k)/kz2_max)**nupa_z_exp &
            )
            call gear2(dt, bpa_new(i, k, j), bpa(i, k, j), bpa_old1(i, k, j), &
               (nonlin     (i, k, j, inonlin_bpa) + zi*kz(k)*upa     (i, k, j) &
                 + q*zi*ky(j)*psi     (i, k, j) )/va2cs2_plus_1, &
               (nonlin_old1(i, k, j, inonlin_bpa) + zi*kz(k)*upa_old1(i, k, j) &
                 + q*zi*ky(j)*psi_old1(i, k, j) )/va2cs2_plus_1, &
               (etapa_x*(kprp2(i, k, j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp)/va2cs2_plus_1 &
            )

            phi_new(i, k, j) = omg_new(i, k, j)*(-kprp2inv(i, k, j))
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      phi_old2 = phi_old1
      phi_old1 = phi

      omg_old2 = omg_old1
      omg_old1 = omg

      psi_old2 = psi_old1
      psi_old1 = psi

      upa_old2 = upa_old1
      upa_old1 = upa

      bpa_old2 = bpa_old1
      bpa_old1 = bpa

      nonlin_old2 = nonlin_old1
      nonlin_old1 = nonlin
      !$omp end workshare

      counter = counter + 1
    ! 3rd order 
    else
      !$omp parallel do private(i, k) schedule(static)
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            call gear3(dt, omg_new(i, k, j), omg(i, k, j), omg_old1(i, k, j), omg_old2(i, k, j), &
               nonlin     (i, k, j, inonlin_omg) - zi*kz(k)*kprp2(i, k, j)*psi     (i, k, j) &
                 - 2.d0*zi*ky(j)*upa(i, k, j), &
               nonlin_old1(i, k, j, inonlin_omg) - zi*kz(k)*kprp2(i, k, j)*psi_old1(i, k, j) &
                 - 2.d0*zi*ky(j)*upa_old1(i, k, j), &
               nonlin_old2(i, k, j, inonlin_omg) - zi*kz(k)*kprp2(i, k, j)*psi_old2(i, k, j) &
                 - 2.d0*zi*ky(j)*upa_old2(i, k, j), &
               nupe_x*(kprp2(i, k, j)/kprp2_max)**nupe_x_exp + nupe_z*(kz2(k)/kz2_max)**nupe_z_exp &
            )
            call gear3(dt, psi_new(i, k, j), psi(i, k, j), psi_old1(i, k, j), psi_old2(i, k, j), &
               nonlin     (i, k, j, inonlin_psi) + zi*kz(k)*phi     (i, k, j), &
               nonlin_old1(i, k, j, inonlin_psi) + zi*kz(k)*phi_old1(i, k, j), &
               nonlin_old2(i, k, j, inonlin_psi) + zi*kz(k)*phi_old2(i, k, j), &
               etape_x*(kprp2(i, k, j)/kprp2_max)**etape_x_exp + etape_z*(kz2(k)/kz2_max)**etape_z_exp &
            )
            call gear3(dt, upa_new(i, k, j), upa(i, k, j), upa_old1(i, k, j), upa_old2(i, k, j), &
               nonlin     (i, k, j, inonlin_upa) + zi*kz(k)*bpa     (i, k, j) &
                 + (2.d0 - q)*zi*ky(j)*phi     (i, k, j), &
               nonlin_old1(i, k, j, inonlin_upa) + zi*kz(k)*bpa_old1(i, k, j) &
                 + (2.d0 - q)*zi*ky(j)*phi_old1(i, k, j), &
               nonlin_old2(i, k, j, inonlin_upa) + zi*kz(k)*bpa_old2(i, k, j) &
                 + (2.d0 - q)*zi*ky(j)*phi_old2(i, k, j), &
               nupa_x*(kprp2(i, k, j)/kprp2_max)**nupa_x_exp + nupa_z*(kz2(k)/kz2_max)**nupa_z_exp &
            )
            call gear3(dt, bpa_new(i, k, j), bpa(i, k, j), bpa_old1(i, k, j), bpa_old2(i, k, j), &
               (nonlin     (i, k, j, inonlin_bpa) + zi*kz(k)*upa     (i, k, j) &
                 + q*zi*ky(j)*psi     (i, k, j) )/va2cs2_plus_1, &
               (nonlin_old1(i, k, j, inonlin_bpa) + zi*kz(k)*upa_old1(i, k, j) &
                 + q*zi*ky(j)*psi_old1(i, k, j) )/va2cs2_plus_1, &
               (nonlin_old2(i, k, j, inonlin_bpa) + zi*kz(k)*upa_old2(i, k, j) &
                 + q*zi*ky(j)*psi_old2(i, k, j) )/va2cs2_plus_1, &
               (etapa_x*(kprp2(i, k, j)/kprp2_max)**etapa_x_exp + etapa_z*(kz2(k)/kz2_max)**etapa_z_exp)/va2cs2_plus_1 &
            )

            phi_new(i, k, j) = omg_new(i, k, j)*(-kprp2inv(i, k, j))
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! values at the previous steps
      !$omp workshare
      phi_old2 = phi_old1
      phi_old1 = phi

      omg_old2 = omg_old1
      omg_old1 = omg

      psi_old2 = psi_old1
      psi_old1 = psi

      upa_old2 = upa_old1
      upa_old1 = upa

      bpa_old2 = bpa_old1
      bpa_old1 = bpa

      nonlin_old2 = nonlin_old1
      nonlin_old1 = nonlin
      !$omp end workshare
    endif

    !$omp workshare
    phi = phi_new
    omg = omg_new
    psi = psi_new
    upa = upa_new
    bpa = bpa_new
    !$omp end workshare

    ! Dealiasing
    if(trim(dealias_scheme) /= '2/3') then
      do k = ikz_st, ikz_en
        do j = iky_st, iky_en
          do i = ikx_st, ikx_en
            phi(i,k,j) = phi_new(i,k,j)*filter(i,k,j)
            psi(i,k,j) = psi_new(i,k,j)*filter(i,k,j)
            omg(i,k,j) = omg_new(i,k,j)*filter(i,k,j)
            upa(i,k,j) = upa_new(i,k,j)*filter(i,k,j)
            bpa(i,k,j) = bpa_new(i,k,j)*filter(i,k,j)
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
    implicit none
    complex(r8), allocatable, dimension(:,:,:) :: src

    allocate(src(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0,0.d0))
    allocate(phi_new  , source=src)
    allocate(phi_old2 , source=src)

    allocate(omg_new  , source=src)
    allocate(omg_old2 , source=src)

    allocate(psi_new  , source=src)
    allocate(psi_old2 , source=src)

    allocate(upa_new  , source=src)
    allocate(upa_old2 , source=src)

    allocate(bpa_new  , source=src)
    allocate(bpa_old2 , source=src)
    deallocate(src)

    allocate(nonlin     (ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran)); nonlin      = 0.d0
    allocate(nonlin_old1(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran)); nonlin_old1 = 0.d0
    allocate(nonlin_old2(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en, nftran)); nonlin_old2 = 0.d0

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
    use fields, only: phi, psi, upa, bpa
    use grid, only: kprp2, nlx, nly, nlz, kx, ky
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: nl_local_tot, nk_local_tot
    use params, only: zi, va2cs2_plus_1
    use mp, only: proc0, max_allreduce
    use time, only: cfl, dt, tt, reset_method, increase_dt
    use time_stamp, only: put_time_stamp, timer_nonlinear_terms, timer_fft
    implicit none

    complex(r8), allocatable, dimension(:,:,:,:) :: wbk
    real   (r8), allocatable, dimension(:,:,:,:) :: wb , wf 

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
          wbk(i,k,j,idphi_dx) = zi*kx(i)                  *phi(i, k, j)
          wbk(i,k,j,idphi_dy) = zi*ky(j)                  *phi(i, k, j)
          wbk(i,k,j,idomg_dx) = zi*kx(i)*(-kprp2(i, k, j))*phi(i, k, j)
          wbk(i,k,j,idomg_dy) = zi*ky(j)*(-kprp2(i, k, j))*phi(i, k, j)
                                 
          wbk(i,k,j,idpsi_dx) = zi*kx(i)                  *psi(i, k, j)
          wbk(i,k,j,idpsi_dy) = zi*ky(j)                  *psi(i, k, j)
          wbk(i,k,j,idjpa_dx) = zi*kx(i)*(-kprp2(i, k, j))*psi(i, k, j)
          wbk(i,k,j,idjpa_dy) = zi*ky(j)*(-kprp2(i, k, j))*psi(i, k, j)

          wbk(i,k,j,idupa_dx) = zi*kx(i)                  *upa(i, k, j)
          wbk(i,k,j,idupa_dy) = zi*ky(j)                  *upa(i, k, j)
          wbk(i,k,j,idbpa_dx) = zi*kx(i)                  *bpa(i, k, j)
          wbk(i,k,j,idbpa_dy) = zi*ky(j)                  *bpa(i, k, j)
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
          wf(j,k,i,inonlin_omg) = - wb(j,k,i,idphi_dx)*wb(j,k,i,idomg_dy) &
                                  + wb(j,k,i,idphi_dy)*wb(j,k,i,idomg_dx) &
                                  + wb(j,k,i,idpsi_dx)*wb(j,k,i,idjpa_dy) &
                                  - wb(j,k,i,idpsi_dy)*wb(j,k,i,idjpa_dx)
          wf(j,k,i,inonlin_psi) = - wb(j,k,i,idphi_dx)*wb(j,k,i,idpsi_dy) &
                                  + wb(j,k,i,idphi_dy)*wb(j,k,i,idpsi_dx)
          wf(j,k,i,inonlin_upa) = - wb(j,k,i,idphi_dx)*wb(j,k,i,idupa_dy) &
                                  + wb(j,k,i,idphi_dy)*wb(j,k,i,idupa_dx) &
                                  + wb(j,k,i,idpsi_dx)*wb(j,k,i,idbpa_dy) &
                                  - wb(j,k,i,idpsi_dy)*wb(j,k,i,idbpa_dx)
          wf(j,k,i,inonlin_bpa) =(- wb(j,k,i,idphi_dx)*wb(j,k,i,idbpa_dy) &
                                  + wb(j,k,i,idphi_dy)*wb(j,k,i,idbpa_dx))*va2cs2_plus_1 &
                                  + wb(j,k,i,idpsi_dx)*wb(j,k,i,idupa_dy) &
                                  - wb(j,k,i,idpsi_dy)*wb(j,k,i,idupa_dx)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! 4. Forward FFT
    if (proc0) call put_time_stamp(timer_fft)
    call p3dfft_ftran_r2c_many(wf, nl_local_tot, nonlin, nk_local_tot, nftran, 'fft')
    if (proc0) call put_time_stamp(timer_fft)
    !$omp workshare
    nonlin = nonlin/nlx/nly/nlz
    !$omp end workshare

    deallocate(wbk)
    deallocate(wb )
    deallocate(wf )

    if (proc0) call put_time_stamp(timer_nonlinear_terms)
  end subroutine get_nonlinear_terms

end module advance

