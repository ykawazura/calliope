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
contains


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
  subroutine gear1(dt, unew, u, exp_terms, imp_terms)
    implicit none
    complex(r8), intent (out) :: unew
    complex(r8), intent (in ) :: u, exp_terms
    real   (r8), intent (in ) :: dt, imp_terms

    unew = (u &
              + dt*exp_terms & 
            )/(1.d0 + dt*imp_terms)
  end subroutine gear1

  subroutine gear2(dt, unew, u, uold1, exp_terms, exp_terms_old1, imp_terms)
    implicit none
    complex(r8), intent (out) :: unew
    complex(r8), intent (in ) :: u, uold1, exp_terms, exp_terms_old1
    real   (r8), intent (in ) :: dt, imp_terms

    unew = (2.d0*u - 0.5d0*uold1 &
             + dt*(2.d0*exp_terms - exp_terms_old1) &
           )/(1.5d0 + dt*imp_terms)
  end subroutine gear2

  subroutine gear3(dt, unew, u, uold1, uold2, exp_terms, exp_terms_old1, exp_terms_old2, imp_terms)
    implicit none
    complex(r8), intent (out) :: unew
    complex(r8), intent (in ) :: u, uold1, uold2, exp_terms, exp_terms_old1, exp_terms_old2
    real   (r8), intent (in ) :: dt, imp_terms

    unew = (3.d0*u - 1.5d0*uold1 + 1.d0/3.d0*uold2 &
             + dt*(3.d0*exp_terms - 3.d0*exp_terms_old1 + exp_terms_old2) & 
           )/(11.d0/6.d0 + dt*imp_terms)
  end subroutine gear3

end module advance_common

