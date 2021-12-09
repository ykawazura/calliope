!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Dealiasing setting
!-----------------------------------------------!
module dealias
  use p3dfft
  implicit none

  public  init_filter

  ! de-aliasing filter; 2/3 or Hou-Li depending on input
  real(r8), allocatable, dimension(:,:,:) :: filter 

contains


!-----------------------------------------------!
!> @author  YK
!! @date    16 Feb 2021
!! @brief   Initialization of filter
!-----------------------------------------------!
  subroutine init_filter
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: proc0
    use params, only: dealias
    implicit none
    integer :: i, j, k

    allocate(filter(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))

    if(proc0 .and. trim(dealias) /= '2/3' .and. trim(dealias) /= '3/2') then
      print *
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      print *, '!             Warning!             !'
      print *, '!      You are not dealising!      !'
      print *, '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!'
      print *
    endif

    ! dealiasing filter
    filter = 1.d0

    if(trim(dealias) == '3/2') then
      do i = ikx_st, ikx_en
        if(abs(kx(i)) >= maxval(abs(kx))*2.d0/3.d0) then
          filter(i, :, :) = 0.d0
        endif
      enddo
      do j = iky_st, iky_en
        if(abs(ky(j)) >= maxval(abs(ky))*2.d0/3.d0) then
          filter(:, :, j) = 0.d0
        endif
      enddo
      do k = ikz_st, ikz_en
        if(abs(kz(k)) >= maxval(abs(kz))*2.d0/3.d0) then
          filter(:, k, :) = 0.d0
        endif
      enddo
    endif

  end subroutine init_filter

end module dealias

