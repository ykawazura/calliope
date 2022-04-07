!-----------------------------------------------!
!> @author  YK
!! @date    14 Apr 2020
!! @brief   Shearing box setting
!-----------------------------------------------!
module shearing_box
  use p3dfft
  implicit none

  public init_shearing_time, remap_fld, get_imp_terms_tintg_with_shear
  public tsc ! time in shearing coordinate. must be 0 < tsc < tremap
  public tremap
  public nremap
  public k2t, k2t_inv, kxt, kxt_old1
  public timer_remap 
  public to_non_shearing_coordinate

  integer :: shear_flg = 0 ! 0 for wo shear, 1 for w shear
  real(r8) :: tsc = 0.d0, tremap = 0.d0
  integer :: nremap = 0

  real(r8) :: timer_remap(2) = 0.d0

  real(r8), allocatable ::  k2t(:, :, :), k2t_inv(:, :, :), kxt(:, :), kxt_old1(:, :)

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of time parameters,
!!          followed by input file reading
!-----------------------------------------------!
  subroutine init_shearing_time
    use mp, only: proc0
    use params, only: shear, q
    use grid, only: lx, ly, kx, ky, kz, k2, k2inv
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    implicit none
    integer :: i, j, k

    allocate(k2t    , source=k2)
    allocate(k2t_inv, source=k2inv)

    allocate(kxt     (ikx_st:ikx_en, iky_st:iky_en))
    allocate(kxt_old1(ikx_st:ikx_en, iky_st:iky_en))
    do j = iky_st, iky_en
      do i = ikx_st, ikx_en
        kxt(i, j) = kx(i)
      enddo
    enddo
    kxt_old1 = kxt

    if(shear) then
      shear_flg = 1
      tremap = ly/(lx*q*shear_flg)
      if(proc0) print '("shear is on; tremap = ", f5.3)', tremap

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
      kxt_old1 = kxt
    endif

  end subroutine init_shearing_time



!-----------------------------------------------!
!> @author  YK
!! @date    12 Jan 2019
!! @brief   Remap fields
!-----------------------------------------------!
  subroutine remap_fld(fld)
    use params, only: pi, q
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use grid, only: ky, dkx, ikx
    implicit none
    complex(r8), dimension (ikx_st:ikx_en, &
                            ikz_st:ikz_en, &
                            iky_st:iky_en), intent(inout) :: fld
    complex(r8), allocatable, dimension(:,:,:) :: tmp
    integer :: i, j, inew, ikxnew, dikx

    allocate(tmp(ikx_st:ikx_en, ikz_st:ikz_en, iky_st:iky_en), source=(0.d0, 0.d0))

    do j = iky_st, iky_en
      dikx = floor(q*shear_flg*tremap*ky(j)/dkx)
        do i = ikx_st, ikx_en
        ikxnew = ikx(i) + dikx
        if(ikxnew <= maxval(ikx) .and. ikxnew >= minval(ikx)) then
          inew = minloc(abs(ikx - ikxnew), 1)
          tmp(inew, :, j) = fld(i, :, j)
        endif
      enddo
    enddo

    fld = tmp

    deallocate(tmp)
  end subroutine remap_fld


!-----------------------------------------------!
!> @author  YK
!! @date    25 Aug 2020
!! @brief   Transform to non-shearing coordinate
!-----------------------------------------------!
  subroutine to_non_shearing_coordinate(fld)
    use grid, only: xx, dly, nly
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use params, only: q
    implicit none
    real(r8), dimension (ily_st:ily_en, &
                         ilz_st:ilz_en, &
                         ilx_st:ilx_en), intent(inout) :: fld
    real(r8), allocatable, dimension(:,:,:) :: tmp
    integer :: i, j, dj, j_nsc

    allocate(tmp(ily_st:ily_en, ilz_st:ilz_en, ilx_st:ilx_en), source=0.d0)

    do i = ilx_st, ilx_en
      dj = floor(q*shear_flg*tsc*xx(i)/dly)
      do j = ily_st, ily_en
        j_nsc = j + dj
        if (j_nsc < 1  ) j_nsc = nly + j_nsc
        if (j_nsc > nly) j_nsc = 1 + (j_nsc - nly)
        tmp(j, :, i) = fld(j_nsc, :, i)
      enddo
    enddo

    fld = tmp

    deallocate(tmp)
  end subroutine to_non_shearing_coordinate


!-----------------------------------------------!
!> @author  YK
!! @date    4 Apr 2022
!! @brief   Time integral of hyperdissipation
!!          -\int \mathrm{d}t\, q[k_x(t)^2 
!!                 + k_y^2 + k_z^2]^{2n}
!-----------------------------------------------!
  subroutine get_imp_terms_tintg_with_shear(imp_terms_tintg, t, kx, ky, kz, coeff, nexp)
    use params, only: q
    use grid, only: k2_max
    implicit none
    real(r8), intent(out) :: imp_terms_tintg
    real(r8), intent(in) :: t, kx, ky, kz, coeff
    integer, intent(in) :: nexp
    real(r8) :: kxt, k2yz

    if(ky == 0.d0) then
      imp_terms_tintg = (kx**2 + ky**2 + kz**2)**nexp*t
    else
      kxt  = kx + q*t*ky
      k2yz = ky**2 + kz**2
      select case(nexp)
      case (1) 
        imp_terms_tintg = (k2yz*kxt + kxt**3/3.d0)/(q*ky)
      case (2) 
        imp_terms_tintg = (k2yz**2*kxt + 2.d0/3.d0*k2yz*kxt**3 + kxt**5/5.d0)/(q*ky)
      case (3) 
        imp_terms_tintg = (k2yz**3*kxt + k2yz**2*kxt**3 + 3.d0/5.d0*k2yz*kxt**5 & 
             + kxt**7/7.d0)/(q*ky)
      case (4) 
        imp_terms_tintg = (k2yz**4*kxt + 4.d0/3.d0*k2yz**3*kxt**3 + 6.d0/5.d0*k2yz**2*kxt**5 & 
             + 4.d0/7.d0*k2yz*kxt**7 + kxt**9/9.d0)/(q*ky)
      end select
    endif
    imp_terms_tintg = -coeff*imp_terms_tintg/k2_max**nexp
  end subroutine get_imp_terms_tintg_with_shear

end module shearing_box

