!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Gear's method
!-----------------------------------------------!
module advance_common
  use p3dfft
  use time, only: dt
  implicit none
  public gear1, gear2, gear3
  public eSSPIFRK1, eSSPIFRK2, eSSPIFRK3
contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Mar 2022
!! @brief   Solve the equation of motion
!!          Strong Stability Preserving 
!!          Integrating Factor Runge–Kutta Method
!!          [Isherwood, Grant, & Gottlieb
!!           Siam J. Numer. Anal, 2018]
!-----------------------------------------------!
  subroutine eSSPIFRK1(u1, u0, exp_terms0, imp_tintg1, imp_tintg0)
    implicit none
    complex(r8), intent (out) :: u1
    complex(r8), intent (in ) :: u0, exp_terms0
    real   (r8), intent (in ) :: imp_tintg1, imp_tintg0

    u1 = dexp(imp_tintg1 - imp_tintg0)*(u0 + 2.d0/3.d0*dt*exp_terms0)
  end subroutine eSSPIFRK1

  subroutine eSSPIFRK2(u2, u1, u0, exp_terms1, imp_tintg2, imp_tintg0)
    implicit none
    complex(r8), intent (out) :: u2
    complex(r8), intent (in ) :: u1, u0, exp_terms1
    real   (r8), intent (in ) :: imp_tintg2, imp_tintg0

    u2 = 2.d0/3.d0*dexp(imp_tintg2 - imp_tintg0)*u0 + 1.d0/3.d0*u1 + 4.d0/9.d0*dt*exp_terms1
  end subroutine eSSPIFRK2

  subroutine eSSPIFRK3(unew, u2, u0, exp_terms2, exp_terms0, imp_tintg3, imp_tintg2, imp_tintg0)
    implicit none
    complex(r8), intent (out) :: unew
    complex(r8), intent (in ) :: u2, u0, exp_terms2, exp_terms0
    real   (r8), intent (in ) :: imp_tintg3, imp_tintg2, imp_tintg0

    unew = dexp(imp_tintg3 - imp_tintg0)*(37.d0/64.d0*u0 + 5.d0/32.d0*dt*exp_terms0) &
           + dexp(imp_tintg3 - imp_tintg2)*(27.d0/64.d0*u2 + 9.d0/16.d0*dt*exp_terms2)
  end subroutine eSSPIFRK3


!-----------------------------------------------!
!> @author  YK
!! @date    17 Feb 2021
!! @brief   Solve the equation of motion
!!          3rd Order Gear's method
!!               linear terms: explicit
!!            nonlinear terms: explicit
!!          dissipation terms: implicit
!!          [Karniadakis and Israeli, JCP 1991]
!-----------------------------------------------!
  subroutine gear1(unew, u, exp_terms, imp_terms)
    implicit none
    complex(r8), intent (out) :: unew
    complex(r8), intent (in ) :: u, exp_terms
    real   (r8), intent (in ) :: imp_terms

    unew = (u &
              + dt*exp_terms & 
            )/(1.d0 + dt*imp_terms)
  end subroutine gear1

  subroutine gear2(unew, u, uold1, exp_terms, exp_terms_old1, imp_terms)
    implicit none
    complex(r8), intent (out) :: unew
    complex(r8), intent (in ) :: u, uold1, exp_terms, exp_terms_old1
    real   (r8), intent (in ) :: imp_terms

    unew = (2.d0*u - 0.5d0*uold1 &
             + dt*(2.d0*exp_terms - exp_terms_old1) &
           )/(1.5d0 + dt*imp_terms)
  end subroutine gear2

  subroutine gear3(unew, u, uold1, uold2, exp_terms, exp_terms_old1, exp_terms_old2, imp_terms)
    implicit none
    complex(r8), intent (out) :: unew
    complex(r8), intent (in ) :: u, uold1, uold2, exp_terms, exp_terms_old1, exp_terms_old2
    real   (r8), intent (in ) :: imp_terms

    unew = (3.d0*u - 1.5d0*uold1 + 1.d0/3.d0*uold2 &
             + dt*(3.d0*exp_terms - 3.d0*exp_terms_old1 + exp_terms_old2) & 
           )/(11.d0/6.d0 + dt*imp_terms)
  end subroutine gear3


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Manually adjust dt while running
!-----------------------------------------------!
  subroutine dt_adjust_while_running
    use mp, only: proc0
    use mp, only: broadcast
    use file, only: get_unused_unit
    use time, only: dt
    implicit none
    include 'mpif.h'
    real(r8) :: dt_new
    integer :: ierror, unit
    logical :: file_exists

    call mpi_barrier (MPI_COMM_WORLD, ierror)

    if (proc0) then
      call get_unused_unit (unit)
      inquire(file='dt_adjust', exist=file_exists)

      if(file_exists) then
        print *
        write (*, '("Manually dt is changed from ", es12.4e3)', advance='no') dt

        ! When the dt adjust file is detected, delete it.
        open(unit=unit,file='dt_adjust',status='old')
        read (unit,*) dt_new
        close(unit, status='delete')

        dt = dt_new
        print '("  to ", es12.4e3)', dt
        print *

      endif

    endif
    
    call broadcast(dt)

  end subroutine dt_adjust_while_running

end module advance_common

