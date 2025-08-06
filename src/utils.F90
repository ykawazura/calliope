!-----------------------------------------------!
!> @author  YK
!! @date    15 Apr 2020
!! @brief   Useful functions
!-----------------------------------------------!
module utils
  use p3dfft
  implicit none

  public :: cabs2
  public :: ranf
  public :: curl
  public :: check_floor

  private

  interface cabs2
    module procedure cabs2_scaler
    module procedure cabs2_array1
    module procedure cabs2_array2
    module procedure cabs2_array3
  end interface cabs2


contains


!-----------------------------------------------!
!> @author  YK
!! @date    39 Nov 2024
!! @brief   Get square of complex
!!          At Fugaku, sometimes abs(z) gives NaN,
!!          but z*conjg(z) is OK
!-----------------------------------------------!
  function cabs2_scaler(z) result(r)
    implicit none
    complex(8), intent(in) :: z
    real(8) :: r

    r = dble(z*conjg(z))

  end function cabs2_scaler

  function cabs2_array1(z) result(r)
    complex(8), intent(in) :: z(:)
    real(8) :: r(size(z))

    r = dble(z*conjg(z))
  end function cabs2_array1

  function cabs2_array2(z) result(r)
    complex(8), intent(in) :: z(:, :)
    real(8) :: r(size(z, 1), size(z, 2))

    r = dble(z*conjg(z))
  end function cabs2_array2

  function cabs2_array3(z) result(r)
    complex(8), intent(in) :: z(:, :, :)
    real(8) :: r(size(z, 1), size(z, 2), size(z, 3))

    r = dble(z*conjg(z))
  end function cabs2_array3


!-----------------------------------------------!
!> @author  YK
!! @date    29 Jun 2021
!! @brief   Get random number
!-----------------------------------------------!
  function ranf()
    use mp, only: proc_id
    implicit none
    integer :: seedsize
    integer, allocatable :: seed(:)
    integer :: i
    real(r8) :: ranf

    call random_seed(size=seedsize)
    allocate(seed(seedsize)) 

    do i = 1, seedsize
      call system_clock(count=seed(i))
    end do
    call random_seed(put=(proc_id+1)*seed(:)) 

    call random_number(ranf)

    deallocate(seed)

  end function ranf


!-----------------------------------------------!
!> @author  YK
!! @date    15 Apr 2020
!! @brief   Get curl
!-----------------------------------------------!
  subroutine curl(fx, fy, fz, curlfx, curlfy, curlfz)
    use grid, only: kx, ky, kz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use params, only: zi
    implicit none
    integer :: i, j, k
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(in ) :: fx, fy, fz
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(out) :: curlfx, curlfy, curlfz

    !$omp parallel do private(i, k) schedule(static)
    do j = iky_st, iky_en
      do k = ikz_st, ikz_en
        do i = ikx_st, ikx_en
          curlfx(i, k, j) = zi*(ky(j)*fz(i, k, j) - kz(k)*fy(i, k, j))
          curlfy(i, k, j) = zi*(kz(k)*fx(i, k, j) - kx(i)*fz(i, k, j))
          curlfz(i, k, j) = zi*(kx(i)*fy(i, k, j) - ky(j)*fx(i, k, j))
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine curl


!-----------------------------------------------!
!> @author  YK
!! @date    25 May 2021
!! @brief   Check if the value is less than the floor value
!-----------------------------------------------!
  subroutine check_floor(u, u_min)
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mp, only: sum_allreduce
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: &
                            u
    real(r8), intent(in) :: u_min
    real(r8):: u_avg
    real(r8), allocatable, dimension(:,:,:) :: u_r
    integer :: i, j, k

    allocate(u_r(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    call p3dfft_btran_c2r(u, u_r, 'tff')

    u_avg = sum(u_r)
    call sum_allreduce(u_avg)
    u_avg = u_avg/nlx/nly/nlz

    do i = ilx_st, ilx_en
      do k = ilz_st, ilz_en
        do j = ily_st, ily_en
          if(u_r(j, k, i) <= 0.d0) print *, '!!!  CAUTION negaive value  !!!'
          if(u_r(j, k, i) < u_min*u_avg) u_r(j, k, i) = u_min*u_avg
        enddo
      enddo
    enddo

    call p3dfft_ftran_r2c(u_r, u, 'fft'); u = u/nlx/nly/nlz 

    deallocate(u_r)

  end subroutine check_floor

end module utils

