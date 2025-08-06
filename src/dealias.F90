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
    use grid, only: kx_max_noprune, ky_max_noprune, kz_max_noprune
    implicit none
    integer :: i, j, k

    allocate(filter(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en))

    ! dealiasing filter
    filter = 1.d0

    do i = ikx_st, ikx_en
      if(abs(kx(i)) >= kx_max_noprune*2.d0/3.d0) then
        filter(i, :, :) = 0.d0
      endif
    enddo
    do j = iky_st, iky_en
      if(abs(ky(j)) >= ky_max_noprune*2.d0/3.d0) then
        filter(:, :, j) = 0.d0
      endif
    enddo
    do k = ikz_st, ikz_en
      if(abs(kz(k)) >= kz_max_noprune*2.d0/3.d0) then
        filter(:, k, :) = 0.d0
      endif
    enddo

  end subroutine init_filter

end module dealias

