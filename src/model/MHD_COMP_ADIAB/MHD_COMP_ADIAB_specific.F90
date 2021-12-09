!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Model specific subroutines
!-----------------------------------------------!
module model_specific
  use force_common
  implicit none

  public init_model_specific

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialize model specific subroutines
!-----------------------------------------------!
  subroutine init_model_specific
    use mp, only: proc0
    use grid, only: kx, ky, kz
    use force, only :init_force, alloc_force
    use shearing_box, only: init_shearing_time
    use params, only: shear
    use fields, only: bx, by, bz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: broadcast
    implicit none
    real(r8) :: bx0, by0, bz0 ! mean magnetic field
    integer :: i, j, k

    call init_shearing_time

    ! if(proc0 .and. maxval(kx**2) /= maxval(ky**2)) then
      ! print *, maxval(kx**2) , maxval(ky**2)
      ! print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      ! print *, '!              Error!              !'
      ! print *, '!  kx_max must be equal to ky_max  !'
      ! print *, '!  i.e., lx/nlx == ly/nly          !'
      ! print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      ! stop
    ! endif

    call init_force
    call alloc_force

    ! show minimum k in unit of omega0/vA
    if(shear) then
      do j = iky_st, iky_en
        do k = ikz_st, ikz_en
          do i = ikx_st, ikx_en
            if(kx(i) == 0.d0 .and. ky(j) == 0.d0 .and. kz(k) == 0.d0) then
              bx0 = abs(bx(i,k,j))
              by0 = abs(by(i,k,j))
              bz0 = abs(bz(i,k,j))
            endif
          end do
        end do
      end do
      call broadcast(bx0)
      call broadcast(by0)
      call broadcast(bz0)

      if(proc0) then
        write(*, "('kx0*vA_x/omega0 = ', f5.2, &
              & ',  ky0*vA_y/omega0 = ', f5.2, &
              & ',  kz0*vA_z/omega0 = ', f5.2)") &
              & bx0*kx(2), by0*ky(2), bz0*kz(2)
        write(*, "(' vA is computed by the initial mean magnetic fields.')")
      endif
    endif

  end subroutine init_model_specific
end module model_specific

