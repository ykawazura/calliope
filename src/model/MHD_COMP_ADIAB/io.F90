!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   IO for RMHD
!-----------------------------------------------!
module io
  use netcdf
  use p3dfft
  use MPI
  implicit none

  public :: init_io, finish_io, loop_io, loop_io_fields, save_restart

  private

  ! MPIIO
  integer :: fh_rho, fh_mx, fh_my, fh_mz, fh_bx, fh_by, fh_bz, fh_sgm
  character(len=100) :: filename
  integer (kind=MPI_OFFSET_KIND) :: filesize_rho, disp_rho
  integer (kind=MPI_OFFSET_KIND) :: filesize_mx , disp_mx
  integer (kind=MPI_OFFSET_KIND) :: filesize_my , disp_my
  integer (kind=MPI_OFFSET_KIND) :: filesize_mz , disp_mz
  integer (kind=MPI_OFFSET_KIND) :: filesize_bx , disp_bx
  integer (kind=MPI_OFFSET_KIND) :: filesize_by , disp_by
  integer (kind=MPI_OFFSET_KIND) :: filesize_bz , disp_bz
  integer (kind=MPI_OFFSET_KIND) :: filesize_sgm, disp_sgm
  integer :: field_time_unit

  ! NETCDF
  integer :: status
  integer, parameter :: kind_nf = kind (NF90_NOERR)
  integer (kind_nf) :: ncid
  integer :: run_id
  integer (kind_nf) :: char10_dim
  ! parameter
  integer :: beta0_id, gamma_id
  integer :: shear_flg_id, q_id
  integer :: nu_id , nu_exp_id
  integer :: eta_id, eta_exp_id
  ! coordinate
  integer :: xx_id, yy_id, zz_id, kx_id, ky_id, kz_id, kpbin_id, tt_id
  ! total energy and max history
  integer :: wkin_sum_id, wmag_sum_id, wprs_sum_id
  integer :: wkin_dot_sum_id, wmag_dot_sum_id, wprs_dot_sum_id
  integer :: wkin_dissip_sum_id, wmag_dissip_sum_id
  integer :: p_u_sum_id
  integer :: smach_rms_id, amach_rms_id, beta_rms_id
  integer :: zp2_sum_id, zm2_sum_id
  ! mean magnetic field
  integer :: b0_id
  ! polar spectrum
  integer :: u2_bin_id, b2_bin_id, rho_bin_id, sgm_bin_id
  integer :: zp2_bin_id, zm2_bin_id
  ! 2d cut
  integer :: rho_r_z0_id, rho_r_x0_id, rho_r_y0_id

  integer :: mx_r_z0_id, mx_r_x0_id, mx_r_y0_id
  integer :: my_r_z0_id, my_r_x0_id, my_r_y0_id
  integer :: mz_r_z0_id, mz_r_x0_id, mz_r_y0_id

  integer :: wx_r_z0_id, wx_r_x0_id, wx_r_y0_id
  integer :: wy_r_z0_id, wy_r_x0_id, wy_r_y0_id
  integer :: wz_r_z0_id, wz_r_x0_id, wz_r_y0_id

  integer :: bx_r_z0_id, bx_r_x0_id, bx_r_y0_id
  integer :: by_r_z0_id, by_r_x0_id, by_r_y0_id
  integer :: bz_r_z0_id, bz_r_x0_id, bz_r_y0_id

  integer :: jx_r_z0_id, jx_r_x0_id, jx_r_y0_id
  integer :: jy_r_z0_id, jy_r_x0_id, jy_r_y0_id
  integer :: jz_r_z0_id, jz_r_x0_id, jz_r_y0_id

  integer :: sgm_r_z0_id, sgm_r_x0_id, sgm_r_y0_id

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim, vec_dim
  integer, dimension (2) :: mean_fld_dim
  integer, dimension (2) :: bin_dim
  integer, dimension (3) :: z0_dim, x0_dim, y0_dim

  integer :: nout

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of IO
!-----------------------------------------------!
  subroutine init_io(nkpolar, kpbin)
    implicit none
    integer, intent(in) :: nkpolar
    real(r8), intent(in) :: kpbin(1:nkpolar)

    call init_io_decomp
    call init_io_netcdf(nkpolar, kpbin)
  end subroutine init_io


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of MPIIO
!-----------------------------------------------!
  subroutine init_io_decomp
    use mp, only: proc0
    use file, only: open_output_file
    implicit none
    integer :: ierr

    ! open file for IO
    filename = 'rho.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_rho, ierr)
    filesize_rho = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_rho, filesize_rho, ierr)  ! guarantee overwriting
    disp_rho = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'mx.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_mx, ierr)
    filesize_mx = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_mx, filesize_mx, ierr)  ! guarantee overwriting
    disp_mx = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'my.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_my, ierr)
    filesize_my = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_my, filesize_my, ierr)  ! guarantee overwriting
    disp_my = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'mz.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_mz, ierr)
    filesize_mz = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_mz, filesize_mz, ierr)  ! guarantee overwriting
    disp_mz = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'bx.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_bx, ierr)
    filesize_bx = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_bx, filesize_bx, ierr)  ! guarantee overwriting
    disp_bx = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'by.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_by, ierr)
    filesize_by = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_by, filesize_by, ierr)  ! guarantee overwriting
    disp_by = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'bz.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_bz, ierr)
    filesize_bz = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_bz, filesize_bz, ierr)  ! guarantee overwriting
    disp_bz = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'sgm.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_sgm, ierr)
    filesize_sgm = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_sgm, filesize_sgm, ierr)  ! guarantee overwriting
    disp_sgm = 0_MPI_OFFSET_KIND

    if(proc0) then
      call open_output_file (field_time_unit, 'field_time.dat')
    endif
  end subroutine init_io_decomp


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of NETCDF
!-----------------------------------------------!
  subroutine init_io_netcdf(nkpolar, kpbin)
    use grid, only: nlx, nly, nlz 
    use grid, only: xx, yy, zz, kx, ky, kz
    use mp, only: proc0
    use params, only: runname, nu, nu_exp, eta, eta_exp, beta0, gamma, q
    use shearing_box, only: shear_flg
    implicit none
    integer, intent(in) :: nkpolar
    real(r8), intent(in) :: kpbin(1:nkpolar)

    if(proc0) then
      filename = trim(runname)//'.out.nc' ! File name
      status = nf90_create (filename, NF90_CLOBBER, ncid)

      status = nf90_put_att (ncid, NF90_GLOBAL, 'title', 'calliope simulation data')
      status = nf90_def_dim (ncid, 'char10', 10, char10_dim)
      status = nf90_def_var (ncid, 'run_info', NF90_CHAR, char10_dim, run_id)
      status = nf90_put_att (ncid, run_id, 'model', _MODEL_)

      status = nf90_def_dim (ncid, 'xx', size(xx), xx_dim)
      status = nf90_def_dim (ncid, 'yy', size(yy), yy_dim)
      status = nf90_def_dim (ncid, 'zz', size(zz), zz_dim)
      status = nf90_def_dim (ncid, 'kx', size(kx), kx_dim)
      status = nf90_def_dim (ncid, 'ky', size(ky), ky_dim)
      status = nf90_def_dim (ncid, 'kz', size(kz), kz_dim)
      status = nf90_def_dim (ncid, 'kpbin', size(kpbin), kpbin_dim)
      status = nf90_def_dim (ncid, 'tt', NF90_UNLIMITED, tt_dim)
      status = nf90_def_dim (ncid, '3', 3, vec_dim)

      status = nf90_def_var (ncid, 'beta0', NF90_DOUBLE, beta0_id)
      status = nf90_def_var (ncid, 'gamma', NF90_DOUBLE, gamma_id)
      status = nf90_def_var (ncid, 'shear_flg', NF90_INT, shear_flg_id)
      status = nf90_def_var (ncid, 'q', NF90_DOUBLE, q_id)
      status = nf90_def_var (ncid, 'nu' , NF90_DOUBLE, nu_id )
      status = nf90_def_var (ncid, 'eta', NF90_DOUBLE, eta_id)
      status = nf90_def_var (ncid, 'nu_exp' , NF90_DOUBLE, nu_exp_id )
      status = nf90_def_var (ncid, 'eta_exp', NF90_DOUBLE, eta_exp_id)

      status = nf90_def_var (ncid, 'xx', NF90_DOUBLE, xx_dim, xx_id)
      status = nf90_def_var (ncid, 'yy', NF90_DOUBLE, yy_dim, yy_id)
      status = nf90_def_var (ncid, 'zz', NF90_DOUBLE, zz_dim, zz_id)
      status = nf90_def_var (ncid, 'kx', NF90_DOUBLE, kx_dim, kx_id)
      status = nf90_def_var (ncid, 'ky', NF90_DOUBLE, ky_dim, ky_id)
      status = nf90_def_var (ncid, 'kz', NF90_DOUBLE, kz_dim, kz_id)
      status = nf90_def_var (ncid, 'kpbin', NF90_DOUBLE, kpbin_dim, kpbin_id)
      status = nf90_def_var (ncid, 'tt', NF90_DOUBLE, tt_dim, tt_id)
      ! energy and max values history
      status = nf90_def_var (ncid, 'wkin_sum', NF90_DOUBLE, tt_dim, wkin_sum_id)
      status = nf90_def_var (ncid, 'wmag_sum', NF90_DOUBLE, tt_dim, wmag_sum_id)
      status = nf90_def_var (ncid, 'wprs_sum', NF90_DOUBLE, tt_dim, wprs_sum_id)
      status = nf90_def_var (ncid, 'wkin_dot_sum', NF90_DOUBLE, tt_dim, wkin_dot_sum_id)
      status = nf90_def_var (ncid, 'wmag_dot_sum', NF90_DOUBLE, tt_dim, wmag_dot_sum_id)
      status = nf90_def_var (ncid, 'wprs_dot_sum', NF90_DOUBLE, tt_dim, wprs_dot_sum_id)
      status = nf90_def_var (ncid, 'wkin_dissip_sum', NF90_DOUBLE, tt_dim, wkin_dissip_sum_id)
      status = nf90_def_var (ncid, 'wmag_dissip_sum', NF90_DOUBLE, tt_dim, wmag_dissip_sum_id)
      status = nf90_def_var (ncid, 'p_u_sum'     , NF90_DOUBLE, tt_dim, p_u_sum_id   )
      status = nf90_def_var (ncid, 'smach_rms'   , NF90_DOUBLE, tt_dim, smach_rms_id )
      status = nf90_def_var (ncid, 'amach_rms'   , NF90_DOUBLE, tt_dim, amach_rms_id )
      status = nf90_def_var (ncid, 'beta_rms'    , NF90_DOUBLE, tt_dim, beta_rms_id  )
      status = nf90_def_var (ncid, 'zp2_sum', NF90_DOUBLE, tt_dim, zp2_sum_id)
      status = nf90_def_var (ncid, 'zm2_sum', NF90_DOUBLE, tt_dim, zm2_sum_id)
      ! polar spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = tt_dim
      status = nf90_def_var (ncid, 'u2_bin', NF90_DOUBLE, bin_dim, u2_bin_id)
      status = nf90_def_var (ncid, 'b2_bin', NF90_DOUBLE, bin_dim, b2_bin_id)
      status = nf90_def_var (ncid, 'rho_bin', NF90_DOUBLE, bin_dim, rho_bin_id)
      status = nf90_def_var (ncid, 'sgm_bin', NF90_DOUBLE, bin_dim, sgm_bin_id)
      status = nf90_def_var (ncid, 'zp2_bin', NF90_DOUBLE, bin_dim, zp2_bin_id)
      status = nf90_def_var (ncid, 'zm2_bin', NF90_DOUBLE, bin_dim, zm2_bin_id)
      ! mean magnetic field
      mean_fld_dim (1) = vec_dim
      mean_fld_dim (2) = tt_dim
      status = nf90_def_var (ncid, 'b0', NF90_DOUBLE, mean_fld_dim, b0_id)
      ! 2d cut
      z0_dim (1) = xx_dim
      z0_dim (2) = yy_dim
      z0_dim (3) = tt_dim

      x0_dim (1) = yy_dim
      x0_dim (2) = zz_dim
      x0_dim (3) = tt_dim

      y0_dim (1) = xx_dim
      y0_dim (2) = zz_dim
      y0_dim (3) = tt_dim
      status = nf90_def_var (ncid, 'rho_r_z0', NF90_DOUBLE, z0_dim, rho_r_z0_id)
      status = nf90_def_var (ncid, 'rho_r_x0', NF90_DOUBLE, x0_dim, rho_r_x0_id)
      status = nf90_def_var (ncid, 'rho_r_y0', NF90_DOUBLE, y0_dim, rho_r_y0_id)

      status = nf90_def_var (ncid, 'mx_r_z0', NF90_DOUBLE, z0_dim, mx_r_z0_id)
      status = nf90_def_var (ncid, 'mx_r_x0', NF90_DOUBLE, x0_dim, mx_r_x0_id)
      status = nf90_def_var (ncid, 'mx_r_y0', NF90_DOUBLE, y0_dim, mx_r_y0_id)
      status = nf90_def_var (ncid, 'my_r_z0', NF90_DOUBLE, z0_dim, my_r_z0_id)
      status = nf90_def_var (ncid, 'my_r_x0', NF90_DOUBLE, x0_dim, my_r_x0_id)
      status = nf90_def_var (ncid, 'my_r_y0', NF90_DOUBLE, y0_dim, my_r_y0_id)
      status = nf90_def_var (ncid, 'mz_r_z0', NF90_DOUBLE, z0_dim, mz_r_z0_id)
      status = nf90_def_var (ncid, 'mz_r_x0', NF90_DOUBLE, x0_dim, mz_r_x0_id)
      status = nf90_def_var (ncid, 'mz_r_y0', NF90_DOUBLE, y0_dim, mz_r_y0_id)

      status = nf90_def_var (ncid, 'wx_r_z0', NF90_DOUBLE, z0_dim, wx_r_z0_id)
      status = nf90_def_var (ncid, 'wx_r_x0', NF90_DOUBLE, x0_dim, wx_r_x0_id)
      status = nf90_def_var (ncid, 'wx_r_y0', NF90_DOUBLE, y0_dim, wx_r_y0_id)
      status = nf90_def_var (ncid, 'wy_r_z0', NF90_DOUBLE, z0_dim, wy_r_z0_id)
      status = nf90_def_var (ncid, 'wy_r_x0', NF90_DOUBLE, x0_dim, wy_r_x0_id)
      status = nf90_def_var (ncid, 'wy_r_y0', NF90_DOUBLE, y0_dim, wy_r_y0_id)
      status = nf90_def_var (ncid, 'wz_r_z0', NF90_DOUBLE, z0_dim, wz_r_z0_id)
      status = nf90_def_var (ncid, 'wz_r_x0', NF90_DOUBLE, x0_dim, wz_r_x0_id)
      status = nf90_def_var (ncid, 'wz_r_y0', NF90_DOUBLE, y0_dim, wz_r_y0_id)

      status = nf90_def_var (ncid, 'bx_r_z0', NF90_DOUBLE, z0_dim, bx_r_z0_id)
      status = nf90_def_var (ncid, 'bx_r_x0', NF90_DOUBLE, x0_dim, bx_r_x0_id)
      status = nf90_def_var (ncid, 'bx_r_y0', NF90_DOUBLE, y0_dim, bx_r_y0_id)
      status = nf90_def_var (ncid, 'by_r_z0', NF90_DOUBLE, z0_dim, by_r_z0_id)
      status = nf90_def_var (ncid, 'by_r_x0', NF90_DOUBLE, x0_dim, by_r_x0_id)
      status = nf90_def_var (ncid, 'by_r_y0', NF90_DOUBLE, y0_dim, by_r_y0_id)
      status = nf90_def_var (ncid, 'bz_r_z0', NF90_DOUBLE, z0_dim, bz_r_z0_id)
      status = nf90_def_var (ncid, 'bz_r_x0', NF90_DOUBLE, x0_dim, bz_r_x0_id)
      status = nf90_def_var (ncid, 'bz_r_y0', NF90_DOUBLE, y0_dim, bz_r_y0_id)

      status = nf90_def_var (ncid, 'jx_r_z0', NF90_DOUBLE, z0_dim, jx_r_z0_id)
      status = nf90_def_var (ncid, 'jx_r_x0', NF90_DOUBLE, x0_dim, jx_r_x0_id)
      status = nf90_def_var (ncid, 'jx_r_y0', NF90_DOUBLE, y0_dim, jx_r_y0_id)
      status = nf90_def_var (ncid, 'jy_r_z0', NF90_DOUBLE, z0_dim, jy_r_z0_id)
      status = nf90_def_var (ncid, 'jy_r_x0', NF90_DOUBLE, x0_dim, jy_r_x0_id)
      status = nf90_def_var (ncid, 'jy_r_y0', NF90_DOUBLE, y0_dim, jy_r_y0_id)
      status = nf90_def_var (ncid, 'jz_r_z0', NF90_DOUBLE, z0_dim, jz_r_z0_id)
      status = nf90_def_var (ncid, 'jz_r_x0', NF90_DOUBLE, x0_dim, jz_r_x0_id)
      status = nf90_def_var (ncid, 'jz_r_y0', NF90_DOUBLE, y0_dim, jz_r_y0_id)

      status = nf90_def_var (ncid, 'sgm_r_z0', NF90_DOUBLE, z0_dim, sgm_r_z0_id)
      status = nf90_def_var (ncid, 'sgm_r_x0', NF90_DOUBLE, x0_dim, sgm_r_x0_id)
      status = nf90_def_var (ncid, 'sgm_r_y0', NF90_DOUBLE, y0_dim, sgm_r_y0_id)


      status = nf90_enddef (ncid)  ! out of definition mode

      status = nf90_put_var (ncid, beta0_id, beta0)
      status = nf90_put_var (ncid, gamma_id, gamma)
      status = nf90_put_var (ncid, shear_flg_id, shear_flg)
      status = nf90_put_var (ncid, q_id, q)
      status = nf90_put_var (ncid, nu_id , nu )
      status = nf90_put_var (ncid, eta_id, eta)
      status = nf90_put_var (ncid, nu_exp_id , dble(nu_exp ))
      status = nf90_put_var (ncid, eta_exp_id, dble(eta_exp))

      status = nf90_put_var (ncid, xx_id, xx)
      status = nf90_put_var (ncid, yy_id, yy)
      status = nf90_put_var (ncid, zz_id, zz)
      status = nf90_put_var (ncid, kx_id, kx)
      status = nf90_put_var (ncid, ky_id, ky)
      status = nf90_put_var (ncid, kz_id, kz)
      status = nf90_put_var (ncid, kpbin_id, kpbin)

      nout = 1
    endif
  end subroutine init_io_netcdf


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Append variables to NETCDF
!-----------------------------------------------!
  subroutine loop_io( &
                      wkin_sum, wmag_sum, wprs_sum, &
                      wkin_dot_sum, wmag_dot_sum, wprs_dot_sum, &
                      wkin_dissip_sum, wmag_dissip_sum, &
                      p_u_sum, &
                      smach_rms, amach_rms, beta_rms, &
                      zp2_sum, zm2_sum, &
                      bx0, by0, bz0, &
                      !
                      nkpolar, &
                      u2_bin, b2_bin, rho_bin, sgm_bin, &
                      zp2_bin, zm2_bin, &
                      !
                      rho_r_z0, rho_r_x0, rho_r_y0, &
                      !
                      mx_r_z0, mx_r_x0, mx_r_y0, &
                      my_r_z0, my_r_x0, my_r_y0, &
                      mz_r_z0, mz_r_x0, mz_r_y0, &
                      !
                      wx_r_z0, wx_r_x0, wx_r_y0, &
                      wy_r_z0, wy_r_x0, wy_r_y0, &
                      wz_r_z0, wz_r_x0, wz_r_y0, &
                      !
                      bx_r_z0, bx_r_x0, bx_r_y0, &
                      by_r_z0, by_r_x0, by_r_y0, &
                      bz_r_z0, bz_r_x0, bz_r_y0, &
                      !
                      jx_r_z0, jx_r_x0, jx_r_y0, &
                      jy_r_z0, jy_r_x0, jy_r_y0, &
                      jz_r_z0, jz_r_x0, jz_r_y0, &
                      !
                      sgm_r_z0, sgm_r_x0, sgm_r_y0 &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz
    use mp, only: proc0
    implicit none
    real(r8), intent(in) :: wkin_sum, wmag_sum, wprs_sum
    real(r8), intent(in) :: wkin_dot_sum, wmag_dot_sum, wprs_dot_sum
    real(r8), intent(in) :: wkin_dissip_sum, wmag_dissip_sum
    real(r8), intent(in) :: p_u_sum
    real(r8), intent(in) :: zp2_sum, zm2_sum
    real(r8), intent(in) :: smach_rms, amach_rms, beta_rms
    real(r8), intent(in) :: bx0, by0, bz0

    integer , intent(in) :: nkpolar
    real(r8), intent(in) :: u2_bin(:), b2_bin(:), rho_bin(:), sgm_bin(:)
    real(r8), intent(in) :: zp2_bin(:), zm2_bin(:)

    real(r8), intent(in) :: rho_r_z0(:,:), rho_r_x0(:,:), rho_r_y0(:,:)

    real(r8), intent(in) :: mx_r_z0(:,:), mx_r_x0(:,:), mx_r_y0(:,:)
    real(r8), intent(in) :: my_r_z0(:,:), my_r_x0(:,:), my_r_y0(:,:)
    real(r8), intent(in) :: mz_r_z0(:,:), mz_r_x0(:,:), mz_r_y0(:,:)

    real(r8), intent(in) :: wx_r_z0(:,:), wx_r_x0(:,:), wx_r_y0(:,:)
    real(r8), intent(in) :: wy_r_z0(:,:), wy_r_x0(:,:), wy_r_y0(:,:)
    real(r8), intent(in) :: wz_r_z0(:,:), wz_r_x0(:,:), wz_r_y0(:,:)

    real(r8), intent(in) :: bx_r_z0(:,:), bx_r_x0(:,:), bx_r_y0(:,:)
    real(r8), intent(in) :: by_r_z0(:,:), by_r_x0(:,:), by_r_y0(:,:)
    real(r8), intent(in) :: bz_r_z0(:,:), bz_r_x0(:,:), bz_r_y0(:,:)

    real(r8), intent(in) :: jx_r_z0(:,:), jx_r_x0(:,:), jx_r_y0(:,:)
    real(r8), intent(in) :: jy_r_z0(:,:), jy_r_x0(:,:), jy_r_y0(:,:)
    real(r8), intent(in) :: jz_r_z0(:,:), jz_r_x0(:,:), jz_r_y0(:,:)
    
    real(r8), intent(in) :: sgm_r_z0(:,:), sgm_r_x0(:,:), sgm_r_y0(:,:)

    integer, dimension (2) :: start2, count2
    integer, dimension (3) :: start3, count3

    ! output via NETCDF
    if(proc0) then
      ! energy and max values history
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))
      status = nf90_put_var (ncid, wkin_sum_id, wkin_sum, start=(/nout/))
      status = nf90_put_var (ncid, wmag_sum_id, wmag_sum, start=(/nout/))
      status = nf90_put_var (ncid, wprs_sum_id, wprs_sum, start=(/nout/))
      status = nf90_put_var (ncid, wkin_dot_sum_id, wkin_dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, wmag_dot_sum_id, wmag_dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, wprs_dot_sum_id, wprs_dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, wkin_dissip_sum_id, wkin_dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, wmag_dissip_sum_id, wmag_dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_u_sum_id, p_u_sum, start=(/nout/))
      status = nf90_put_var (ncid, zp2_sum_id, zp2_sum, start=(/nout/))
      status = nf90_put_var (ncid, zm2_sum_id, zm2_sum, start=(/nout/))
      status = nf90_put_var (ncid, smach_rms_id, smach_rms, start=(/nout/))
      status = nf90_put_var (ncid, amach_rms_id, amach_rms, start=(/nout/))
      status = nf90_put_var (ncid, beta_rms_id , beta_rms , start=(/nout/))
      ! mean magnetic field
      start2(1) = 1
      start2(2) = nout

      count2(1) = 3
      count2(2) = 1
      status = nf90_put_var (ncid, b0_id, (/bx0, by0, bz0/), start=start2, count=count2)
      ! polar spectrum
      start2(1) = 1
      start2(2) = nout

      count2(1) = nkpolar
      count2(2) = 1
      status = nf90_put_var (ncid, u2_bin_id, u2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, b2_bin_id, b2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, rho_bin_id, rho_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, sgm_bin_id, sgm_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, zp2_bin_id, zp2_bin, start=start2, count=count2)
      status = nf90_put_var (ncid, zm2_bin_id, zm2_bin, start=start2, count=count2)
      ! z=0 cut
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout

      count3(1) = nlx
      count3(2) = nly
      count3(3) = 1
      status = nf90_put_var (ncid, rho_r_z0_id, rho_r_z0, start=start3, count=count3)

      status = nf90_put_var (ncid, mx_r_z0_id, mx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, my_r_z0_id, my_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, mz_r_z0_id, mz_r_z0, start=start3, count=count3)

      status = nf90_put_var (ncid, wx_r_z0_id, wx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, wy_r_z0_id, wy_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, wz_r_z0_id, wz_r_z0, start=start3, count=count3)

      status = nf90_put_var (ncid, bx_r_z0_id, bx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, by_r_z0_id, by_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, bz_r_z0_id, bz_r_z0, start=start3, count=count3)

      status = nf90_put_var (ncid, jx_r_z0_id, jx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, jy_r_z0_id, jy_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid, jz_r_z0_id, jz_r_z0, start=start3, count=count3)

      status = nf90_put_var (ncid, sgm_r_z0_id, sgm_r_z0, start=start3, count=count3)
      ! x=0 cut
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout

      count3(1) = nly
      count3(2) = nlz
      count3(3) = 1
      status = nf90_put_var (ncid, rho_r_x0_id, rho_r_x0, start=start3, count=count3)

      status = nf90_put_var (ncid, mx_r_x0_id, mx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, my_r_x0_id, my_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, mz_r_x0_id, mz_r_x0, start=start3, count=count3)

      status = nf90_put_var (ncid, wx_r_x0_id, wx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, wy_r_x0_id, wy_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, wz_r_x0_id, wz_r_x0, start=start3, count=count3)

      status = nf90_put_var (ncid, bx_r_x0_id, bx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, by_r_x0_id, by_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, bz_r_x0_id, bz_r_x0, start=start3, count=count3)

      status = nf90_put_var (ncid, jx_r_x0_id, jx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, jy_r_x0_id, jy_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid, jz_r_x0_id, jz_r_x0, start=start3, count=count3)

      status = nf90_put_var (ncid, sgm_r_x0_id, sgm_r_x0, start=start3, count=count3)
      ! y=0 cut
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout

      count3(1) = nlx
      count3(2) = nlz
      count3(3) = 1
      status = nf90_put_var (ncid, rho_r_y0_id, rho_r_y0, start=start3, count=count3)

      status = nf90_put_var (ncid, mx_r_y0_id, mx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, my_r_y0_id, my_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, mz_r_y0_id, mz_r_y0, start=start3, count=count3)

      status = nf90_put_var (ncid, wx_r_y0_id, wx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, wy_r_y0_id, wy_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, wz_r_y0_id, wz_r_y0, start=start3, count=count3)

      status = nf90_put_var (ncid, bx_r_y0_id, bx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, by_r_y0_id, by_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, bz_r_y0_id, bz_r_y0, start=start3, count=count3)

      status = nf90_put_var (ncid, jx_r_y0_id, jx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, jy_r_y0_id, jy_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid, jz_r_y0_id, jz_r_y0, start=start3, count=count3)

      status = nf90_put_var (ncid, sgm_r_y0_id, sgm_r_y0, start=start3, count=count3)

      nout = nout + 1

      status = nf90_sync (ncid)
    endif
  end subroutine loop_io


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Append field variables via MPIIO
!-----------------------------------------------!
  subroutine loop_io_fields
    use fields, only: rho, sgm
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use mp, only: proc0
    use time, only: tt
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mpiio, only: mpiio_write_var
    implicit none
    integer, dimension(3) :: sizes, subsizes, starts

    sizes(1) = nkx
    sizes(2) = nkz
    sizes(3) = nky
    subsizes(1) = ikx_en - ikx_st + 1
    subsizes(2) = ikz_en - ikz_st + 1
    subsizes(3) = iky_en - iky_st + 1
    starts(1) = ikx_st - 1
    starts(2) = ikz_st - 1
    starts(3) = iky_st - 1

    call mpiio_write_var(fh_rho, disp_rho, sizes, subsizes, starts, rho)
    call mpiio_write_var(fh_mx , disp_mx , sizes, subsizes, starts, mx )
    call mpiio_write_var(fh_my , disp_my , sizes, subsizes, starts, my )
    call mpiio_write_var(fh_mz , disp_mz , sizes, subsizes, starts, mz )
    call mpiio_write_var(fh_bx , disp_bx , sizes, subsizes, starts, bx )
    call mpiio_write_var(fh_by , disp_by , sizes, subsizes, starts, by )
    call mpiio_write_var(fh_bz , disp_bz , sizes, subsizes, starts, bz )
    call mpiio_write_var(fh_sgm, disp_sgm, sizes, subsizes, starts, sgm)

    if(proc0) then
      write (unit=field_time_unit, fmt="(100es30.21)") tt
      flush (field_time_unit)
    endif
  end subroutine loop_io_fields


!-----------------------------------------------!
!> @author  YK
!! @date    16 Jan 2019
!! @brief   Save restart file via MPIIO
!-----------------------------------------------!
  subroutine save_restart
    use fields, only: rho, sgm
    use fields, only: mx, my, mz
    use fields, only: bx, by, bz
    use mp, only: proc0
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: tt, dt
    use params, only: restart_dir
    use file, only: open_output_file, close_file
    use mpiio, only: mpiio_write_one
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts

    sizes(1) = nkx
    sizes(2) = nkz
    sizes(3) = nky
    subsizes(1) = ikx_en - ikx_st + 1
    subsizes(2) = ikz_en - ikz_st + 1
    subsizes(3) = iky_en - iky_st + 1
    starts(1) = ikx_st - 1
    starts(2) = ikz_st - 1
    starts(3) = iky_st - 1

    call mpiio_write_one(rho, sizes, subsizes, starts, trim(restart_dir)//'rho.dat')
    call mpiio_write_one(mx , sizes, subsizes, starts, trim(restart_dir)//'mx.dat' )
    call mpiio_write_one(my , sizes, subsizes, starts, trim(restart_dir)//'my.dat' )
    call mpiio_write_one(mz , sizes, subsizes, starts, trim(restart_dir)//'mz.dat' )
    call mpiio_write_one(bx , sizes, subsizes, starts, trim(restart_dir)//'bx.dat' )
    call mpiio_write_one(by , sizes, subsizes, starts, trim(restart_dir)//'by.dat' )
    call mpiio_write_one(bz , sizes, subsizes, starts, trim(restart_dir)//'bz.dat' )
    call mpiio_write_one(sgm, sizes, subsizes, starts, trim(restart_dir)//'sgm.dat')

    if(proc0) then
      call open_output_file (time_unit, trim(restart_dir)//'time.dat')
      write (unit=time_unit, fmt="(3X, 'tt')")
      write (unit=time_unit, fmt="(100es30.21)") tt
      call close_file (time_unit)
    endif
  end subroutine save_restart


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Finalization of NETCDF
!-----------------------------------------------!
  subroutine finish_io
    use mp, only: proc0
    use file, only: close_file
    implicit none
    integer :: ierr

    call MPI_FILE_CLOSE(fh_rho,ierr)
    call MPI_FILE_CLOSE(fh_mx ,ierr)
    call MPI_FILE_CLOSE(fh_my ,ierr)
    call MPI_FILE_CLOSE(fh_mz ,ierr)
    call MPI_FILE_CLOSE(fh_bx ,ierr)
    call MPI_FILE_CLOSE(fh_by ,ierr)
    call MPI_FILE_CLOSE(fh_bz ,ierr)
    call MPI_FILE_CLOSE(fh_sgm,ierr)
    if(proc0) then
      call close_file (field_time_unit)
    endif

    if(proc0) then
      status = nf90_close (ncid)
    endif

  end subroutine finish_io

end module io

