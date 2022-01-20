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

  public :: init_io, finish_io, loop_io, loop_io_fields_section, loop_io_kpar, &
            loop_io_SF2, loop_io_fields_3D, save_restart

  private

  ! MPIIO
  integer :: fh_ux, fh_uy, fh_uz, fh_bx, fh_by, fh_bz
  character(len=100) :: filename
  integer (kind=MPI_OFFSET_KIND) :: filesize_ux, disp_ux
  integer (kind=MPI_OFFSET_KIND) :: filesize_uy, disp_uy
  integer (kind=MPI_OFFSET_KIND) :: filesize_uz, disp_uz
  integer (kind=MPI_OFFSET_KIND) :: filesize_bx, disp_bx
  integer (kind=MPI_OFFSET_KIND) :: filesize_by, disp_by
  integer (kind=MPI_OFFSET_KIND) :: filesize_bz, disp_bz
  integer :: field_time_unit

  ! NETCDF for regular output file
  integer :: status
  integer, parameter :: kind_nf = kind (NF90_NOERR)
  integer (kind_nf) :: ncid
  integer :: run_id
  integer (kind_nf) :: char10_dim
  ! parameter
  integer :: shear_flg_id, q_id
  integer :: nu_id , nu_exp_id
  integer :: eta_id, eta_exp_id
  ! coordinate
  integer :: xx_id, yy_id, zz_id, kx_id, ky_id, kz_id, kpbin_id, tt_id
  ! total energy
  integer :: u2_sum_id, b2_sum_id
  integer :: u2dot_sum_id, b2dot_sum_id
  integer :: u2dissip_sum_id, b2dissip_sum_id
  integer :: p_u_sum_id
  integer :: zp2_sum_id, zm2_sum_id
  ! mean magnetic field
  integer :: b0_id
  ! polar spectrum
  integer :: u2_bin_id, ux2_bin_id, uy2_bin_id, uz2_bin_id
  integer :: b2_bin_id, bx2_bin_id, by2_bin_id, bz2_bin_id
  integer :: zp2_bin_id, zm2_bin_id

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim, vec_dim
  integer, dimension (2) :: mean_fld_dim
  integer, dimension (2) :: bin_dim

  integer :: nout

  ! NETCDF for cross section output file
  integer (kind_nf) :: ncid_fld_section
  integer :: xx_fld_section_id, yy_fld_section_id, zz_fld_section_id, tt_fld_section_id
  integer :: kx_fld_section_id, ky_fld_section_id, kz_fld_section_id

  integer :: ux_r_z0_id, ux_r_x0_id, ux_r_y0_id
  integer :: uy_r_z0_id, uy_r_x0_id, uy_r_y0_id
  integer :: uz_r_z0_id, uz_r_x0_id, uz_r_y0_id

  integer :: wx_r_z0_id, wx_r_x0_id, wx_r_y0_id
  integer :: wy_r_z0_id, wy_r_x0_id, wy_r_y0_id
  integer :: wz_r_z0_id, wz_r_x0_id, wz_r_y0_id

  integer :: bx_r_z0_id, bx_r_x0_id, bx_r_y0_id
  integer :: by_r_z0_id, by_r_x0_id, by_r_y0_id
  integer :: bz_r_z0_id, bz_r_x0_id, bz_r_y0_id

  integer :: jx_r_z0_id, jx_r_x0_id, jx_r_y0_id
  integer :: jy_r_z0_id, jy_r_x0_id, jy_r_y0_id
  integer :: jz_r_z0_id, jz_r_x0_id, jz_r_y0_id

  integer :: u2_kxy_id, u2_kyz_id, u2_kxz_id
  integer :: b2_kxy_id, b2_kyz_id, b2_kxz_id

  integer (kind_nf) :: xx_fld_section_dim, yy_fld_section_dim, zz_fld_section_dim, tt_fld_section_dim
  integer (kind_nf) :: kx_fld_section_dim, ky_fld_section_dim, kz_fld_section_dim
  integer, dimension (3) ::  z0_dim,  x0_dim,  y0_dim
  integer, dimension (3) :: kxy_dim, kyz_dim, kxz_dim

  integer :: nout_fld_section

  ! NETCDF for kpar output file
  integer (kind_nf) :: ncid_kpar
  ! coordinate
  integer :: tt_kpar_id, kpbin_kpar_id
  integer :: kpar_b_id, kpar_u_id, b1_ovr_b0_id
  integer (kind_nf) :: tt_kpar_dim, kpbin_kpar_dim
  integer, dimension (2) :: kpar_dim
  integer :: nout_kpar

  ! NETCDF for 2nd order structure function output file
  integer (kind_nf) :: ncid_SF2
  integer :: lpar_SF2_id, lper_SF2_id, tt_SF2_id

  integer :: SF2b_id, SF2u_id

  integer (kind_nf) :: nl_SF2_dim, tt_SF2_dim
  integer, dimension (3) :: SF2_dim

  integer :: nout_SF2

contains


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of IO
!-----------------------------------------------!
  subroutine init_io(nkpolar, kpbin, nl, lpar, lper)
    implicit none
    integer, intent(in) :: nkpolar, nl
    real(r8), intent(in) :: kpbin(1:nkpolar), lpar(nl), lper(nl)

    call init_io_decomp
    call init_io_netcdf(nkpolar, kpbin, nl, lpar, lper)
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
    filename = 'ux.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_ux, ierr)
    filesize_ux = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_ux, filesize_ux, ierr)  ! guarantee overwriting
    disp_ux = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'uy.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_uy, ierr)
    filesize_uy = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_uy, filesize_uy, ierr)  ! guarantee overwriting
    disp_uy = 0_MPI_OFFSET_KIND

    ! open file for IO
    filename = 'uz.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_uz, ierr)
    filesize_uz = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh_uz, filesize_uz, ierr)  ! guarantee overwriting
    disp_uz = 0_MPI_OFFSET_KIND

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

    if(proc0) then
      call open_output_file (field_time_unit, 'field_time.dat')
    endif
  end subroutine init_io_decomp


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of NETCDF
!-----------------------------------------------!
  subroutine init_io_netcdf(nkpolar, kpbin, nl, lpar, lper)
    use grid, only: nlx, nly, nlz 
    use grid, only: xx, yy, zz, kx, ky, kz
    use mp, only: proc0
    use params, only: runname, nu, nu_exp, eta, eta_exp, q
    use shearing_box, only: shear_flg
    implicit none
    integer, intent(in) :: nkpolar, nl
    real(r8), intent(in) :: kpbin(1:nkpolar), lpar(nl), lper(nl)

    if(proc0) then
      !--------------------------------------------------!
      ! Output for parameters, time history, and spectra
      !--------------------------------------------------!
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
      ! total energy
      status = nf90_def_var (ncid, 'u2_sum', NF90_DOUBLE, tt_dim, u2_sum_id)
      status = nf90_def_var (ncid, 'b2_sum', NF90_DOUBLE, tt_dim, b2_sum_id)
      status = nf90_def_var (ncid, 'u2dot_sum', NF90_DOUBLE, tt_dim, u2dot_sum_id)
      status = nf90_def_var (ncid, 'b2dot_sum', NF90_DOUBLE, tt_dim, b2dot_sum_id)
      status = nf90_def_var (ncid, 'u2dissip_sum', NF90_DOUBLE, tt_dim, u2dissip_sum_id)
      status = nf90_def_var (ncid, 'b2dissip_sum', NF90_DOUBLE, tt_dim, b2dissip_sum_id)
      status = nf90_def_var (ncid, 'p_u_sum'     , NF90_DOUBLE, tt_dim, p_u_sum_id   )
      status = nf90_def_var (ncid, 'zp2_sum', NF90_DOUBLE, tt_dim, zp2_sum_id)
      status = nf90_def_var (ncid, 'zm2_sum', NF90_DOUBLE, tt_dim, zm2_sum_id)
      ! polar spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = tt_dim
      status = nf90_def_var (ncid,  'u2_bin'  , NF90_DOUBLE, bin_dim,  u2_bin_id  )
      status = nf90_def_var (ncid, 'ux2_bin'  , NF90_DOUBLE, bin_dim, ux2_bin_id  )
      status = nf90_def_var (ncid, 'uy2_bin'  , NF90_DOUBLE, bin_dim, uy2_bin_id  )
      status = nf90_def_var (ncid, 'uz2_bin'  , NF90_DOUBLE, bin_dim, uz2_bin_id  )
      status = nf90_def_var (ncid,  'b2_bin'  , NF90_DOUBLE, bin_dim,  b2_bin_id  )
      status = nf90_def_var (ncid, 'bx2_bin'  , NF90_DOUBLE, bin_dim, bx2_bin_id  )
      status = nf90_def_var (ncid, 'by2_bin'  , NF90_DOUBLE, bin_dim, by2_bin_id  )
      status = nf90_def_var (ncid, 'bz2_bin'  , NF90_DOUBLE, bin_dim, bz2_bin_id  )
      status = nf90_def_var (ncid, 'zp2_bin'  , NF90_DOUBLE, bin_dim, zp2_bin_id  )
      status = nf90_def_var (ncid, 'zm2_bin'  , NF90_DOUBLE, bin_dim, zm2_bin_id  )
      ! mean magnetic field
      mean_fld_dim (1) = vec_dim
      mean_fld_dim (2) = tt_dim
      status = nf90_def_var (ncid, 'b0', NF90_DOUBLE, mean_fld_dim, b0_id)

      status = nf90_enddef (ncid)  ! out of definition mode

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

      !--------------------------------------------------!
      ! Output for cross sections of fields
      !--------------------------------------------------!
      filename = trim(runname)//'.out.fields_section.nc' ! File name
      status = nf90_create (filename, NF90_CLOBBER, ncid_fld_section)

      status = nf90_put_att (ncid_fld_section, NF90_GLOBAL, 'title', 'calliope simulation data')
      status = nf90_def_dim (ncid_fld_section, 'char10', 10, char10_dim)
      status = nf90_def_var (ncid_fld_section, 'run_info', NF90_CHAR, char10_dim, run_id)
      status = nf90_put_att (ncid_fld_section, run_id, 'model', _MODEL_)

      status = nf90_def_dim (ncid_fld_section, 'xx', size(xx), xx_fld_section_dim)
      status = nf90_def_dim (ncid_fld_section, 'yy', size(yy), yy_fld_section_dim)
      status = nf90_def_dim (ncid_fld_section, 'zz', size(zz), zz_fld_section_dim)
      status = nf90_def_dim (ncid_fld_section, 'kx', size(kx), kx_fld_section_dim)
      status = nf90_def_dim (ncid_fld_section, 'ky', size(ky), ky_fld_section_dim)
      status = nf90_def_dim (ncid_fld_section, 'kz', size(kz), kz_fld_section_dim)
      status = nf90_def_dim (ncid_fld_section, 'tt', NF90_UNLIMITED, tt_fld_section_dim)

      status = nf90_def_var (ncid_fld_section, 'xx', NF90_DOUBLE, xx_fld_section_dim, xx_fld_section_id)
      status = nf90_def_var (ncid_fld_section, 'yy', NF90_DOUBLE, yy_fld_section_dim, yy_fld_section_id)
      status = nf90_def_var (ncid_fld_section, 'zz', NF90_DOUBLE, zz_fld_section_dim, zz_fld_section_id)
      status = nf90_def_var (ncid_fld_section, 'kx', NF90_DOUBLE, kx_fld_section_dim, kx_fld_section_id)
      status = nf90_def_var (ncid_fld_section, 'ky', NF90_DOUBLE, ky_fld_section_dim, ky_fld_section_id)
      status = nf90_def_var (ncid_fld_section, 'kz', NF90_DOUBLE, kz_fld_section_dim, kz_fld_section_id)
      status = nf90_def_var (ncid_fld_section, 'tt', NF90_DOUBLE, tt_fld_section_dim, tt_fld_section_id)

      z0_dim (1) = xx_fld_section_dim
      z0_dim (2) = yy_fld_section_dim
      z0_dim (3) = tt_fld_section_dim

      x0_dim (1) = yy_fld_section_dim
      x0_dim (2) = zz_fld_section_dim
      x0_dim (3) = tt_fld_section_dim

      y0_dim (1) = xx_fld_section_dim
      y0_dim (2) = zz_fld_section_dim
      y0_dim (3) = tt_fld_section_dim

      kxy_dim(1) = kx_fld_section_dim
      kxy_dim(2) = ky_fld_section_dim
      kxy_dim(3) = tt_fld_section_dim

      kyz_dim(1) = ky_fld_section_dim
      kyz_dim(2) = kz_fld_section_dim
      kyz_dim(3) = tt_fld_section_dim

      kxz_dim(1) = kx_fld_section_dim
      kxz_dim(2) = kz_fld_section_dim
      kxz_dim(3) = tt_fld_section_dim
      status = nf90_def_var (ncid_fld_section, 'ux_r_z0', NF90_DOUBLE, z0_dim, ux_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'ux_r_x0', NF90_DOUBLE, x0_dim, ux_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'ux_r_y0', NF90_DOUBLE, y0_dim, ux_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'uy_r_z0', NF90_DOUBLE, z0_dim, uy_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'uy_r_x0', NF90_DOUBLE, x0_dim, uy_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'uy_r_y0', NF90_DOUBLE, y0_dim, uy_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'uz_r_z0', NF90_DOUBLE, z0_dim, uz_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'uz_r_x0', NF90_DOUBLE, x0_dim, uz_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'uz_r_y0', NF90_DOUBLE, y0_dim, uz_r_y0_id)
                                             
      status = nf90_def_var (ncid_fld_section, 'wx_r_z0', NF90_DOUBLE, z0_dim, wx_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'wx_r_x0', NF90_DOUBLE, x0_dim, wx_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'wx_r_y0', NF90_DOUBLE, y0_dim, wx_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'wy_r_z0', NF90_DOUBLE, z0_dim, wy_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'wy_r_x0', NF90_DOUBLE, x0_dim, wy_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'wy_r_y0', NF90_DOUBLE, y0_dim, wy_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'wz_r_z0', NF90_DOUBLE, z0_dim, wz_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'wz_r_x0', NF90_DOUBLE, x0_dim, wz_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'wz_r_y0', NF90_DOUBLE, y0_dim, wz_r_y0_id)
                                             
      status = nf90_def_var (ncid_fld_section, 'bx_r_z0', NF90_DOUBLE, z0_dim, bx_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'bx_r_x0', NF90_DOUBLE, x0_dim, bx_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'bx_r_y0', NF90_DOUBLE, y0_dim, bx_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'by_r_z0', NF90_DOUBLE, z0_dim, by_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'by_r_x0', NF90_DOUBLE, x0_dim, by_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'by_r_y0', NF90_DOUBLE, y0_dim, by_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'bz_r_z0', NF90_DOUBLE, z0_dim, bz_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'bz_r_x0', NF90_DOUBLE, x0_dim, bz_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'bz_r_y0', NF90_DOUBLE, y0_dim, bz_r_y0_id)
                                             
      status = nf90_def_var (ncid_fld_section, 'jx_r_z0', NF90_DOUBLE, z0_dim, jx_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'jx_r_x0', NF90_DOUBLE, x0_dim, jx_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'jx_r_y0', NF90_DOUBLE, y0_dim, jx_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'jy_r_z0', NF90_DOUBLE, z0_dim, jy_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'jy_r_x0', NF90_DOUBLE, x0_dim, jy_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'jy_r_y0', NF90_DOUBLE, y0_dim, jy_r_y0_id)
      status = nf90_def_var (ncid_fld_section, 'jz_r_z0', NF90_DOUBLE, z0_dim, jz_r_z0_id)
      status = nf90_def_var (ncid_fld_section, 'jz_r_x0', NF90_DOUBLE, x0_dim, jz_r_x0_id)
      status = nf90_def_var (ncid_fld_section, 'jz_r_y0', NF90_DOUBLE, y0_dim, jz_r_y0_id)

      status = nf90_def_var (ncid_fld_section, 'u2_kxy', NF90_DOUBLE, kxy_dim, u2_kxy_id)
      status = nf90_def_var (ncid_fld_section, 'u2_kyz', NF90_DOUBLE, kyz_dim, u2_kyz_id)
      status = nf90_def_var (ncid_fld_section, 'u2_kxz', NF90_DOUBLE, kxz_dim, u2_kxz_id)
      status = nf90_def_var (ncid_fld_section, 'b2_kxy', NF90_DOUBLE, kxy_dim, b2_kxy_id)
      status = nf90_def_var (ncid_fld_section, 'b2_kyz', NF90_DOUBLE, kyz_dim, b2_kyz_id)
      status = nf90_def_var (ncid_fld_section, 'b2_kxz', NF90_DOUBLE, kxz_dim, b2_kxz_id)

      status = nf90_enddef (ncid_fld_section)  ! out of definition mode

      status = nf90_put_var (ncid_fld_section, xx_fld_section_id, xx)
      status = nf90_put_var (ncid_fld_section, yy_fld_section_id, yy)
      status = nf90_put_var (ncid_fld_section, zz_fld_section_id, zz)
      status = nf90_put_var (ncid_fld_section, kx_fld_section_id, kx)
      status = nf90_put_var (ncid_fld_section, ky_fld_section_id, ky)
      status = nf90_put_var (ncid_fld_section, kz_fld_section_id, kz)

      nout_fld_section = 1

      !--------------------------------------------------!
      ! Output for kpar
      !--------------------------------------------------!
      filename = trim(runname)//'.out.kpar.nc' ! File name
      status = nf90_create (filename, NF90_CLOBBER, ncid_kpar)

      status = nf90_put_att (ncid_kpar, NF90_GLOBAL, 'title', 'calliope simulation data')
      status = nf90_def_dim (ncid_kpar, 'char10', 10, char10_dim)
      status = nf90_def_var (ncid_kpar, 'run_info', NF90_CHAR, char10_dim, run_id)
      status = nf90_put_att (ncid_kpar, run_id, 'model', _MODEL_)

      status = nf90_def_dim (ncid_kpar, 'tt', NF90_UNLIMITED, tt_kpar_dim)
      status = nf90_def_dim (ncid_kpar, 'kpbin', size(kpbin), kpbin_kpar_dim)

      status = nf90_def_var (ncid_kpar, 'tt'  , NF90_DOUBLE, tt_kpar_dim, tt_kpar_id)
      status = nf90_def_var (ncid_kpar, 'kpbin', NF90_DOUBLE, kpbin_kpar_dim, kpbin_kpar_id)

      kpar_dim (1) = kpbin_kpar_dim
      kpar_dim (2) = tt_kpar_dim

      status = nf90_def_var (ncid_kpar, 'kpar_b'   , NF90_DOUBLE, kpar_dim, kpar_b_id)
      status = nf90_def_var (ncid_kpar, 'kpar_u'   , NF90_DOUBLE, kpar_dim, kpar_u_id)
      status = nf90_def_var (ncid_kpar, 'b1_ovr_b0', NF90_DOUBLE, kpar_dim, b1_ovr_b0_id)

      status = nf90_enddef (ncid_kpar)  ! out of definition mode

      status = nf90_put_var (ncid_kpar, kpbin_kpar_id, kpbin)

      nout_kpar = 1

      !--------------------------------------------------!
      ! Output for 2D order structure function
      !--------------------------------------------------!
      filename = trim(runname)//'.out.SF2.nc' ! File name
      status = nf90_create (filename, NF90_CLOBBER, ncid_SF2)

      status = nf90_put_att (ncid_SF2, NF90_GLOBAL, 'title', 'calliope simulation data')
      status = nf90_def_dim (ncid_SF2, 'char10', 10, char10_dim)
      status = nf90_def_var (ncid_SF2, 'run_info', NF90_CHAR, char10_dim, run_id)
      status = nf90_put_att (ncid_SF2, run_id, 'model', _MODEL_)

      status = nf90_def_dim (ncid_SF2, 'nl', nl, nl_SF2_dim)
      status = nf90_def_dim (ncid_SF2, 'tt', NF90_UNLIMITED, tt_SF2_dim)

      status = nf90_def_var (ncid_SF2, 'lpar', NF90_DOUBLE, nl_SF2_dim, lpar_SF2_id)
      status = nf90_def_var (ncid_SF2, 'lper', NF90_DOUBLE, nl_SF2_dim, lper_SF2_id)
      status = nf90_def_var (ncid_SF2, 'tt'  , NF90_DOUBLE, tt_SF2_dim, tt_SF2_id)

      SF2_dim (1) = nl_SF2_dim
      SF2_dim (2) = nl_SF2_dim
      SF2_dim (3) = tt_SF2_dim

      status = nf90_def_var (ncid_SF2, 'SF2b', NF90_DOUBLE, SF2_dim, SF2b_id)
      status = nf90_def_var (ncid_SF2, 'SF2u', NF90_DOUBLE, SF2_dim, SF2u_id)

      status = nf90_enddef (ncid_SF2)  ! out of definition mode

      status = nf90_put_var (ncid_SF2, lpar_SF2_id, lpar)
      status = nf90_put_var (ncid_SF2, lper_SF2_id, lper)

      nout_SF2 = 1
    endif
  end subroutine init_io_netcdf


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Append variables to NETCDF
!           for params, time history & spectra
!-----------------------------------------------!
  subroutine loop_io( &
                      u2_sum, b2_sum, &
                      u2dot_sum, b2dot_sum, &
                      u2dissip_sum, b2dissip_sum, &
                      p_u_sum, &
                      zp2_sum, zm2_sum, &
                      bx0, by0, bz0, &
                      !
                      nkpolar, &
                      u2_bin, ux2_bin, uy2_bin, uz2_bin, &
                      b2_bin, bx2_bin, by2_bin, bz2_bin, &
                      zp2_bin, zm2_bin &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz
    use mp, only: proc0
    implicit none
    real(r8), intent(in) :: u2_sum, b2_sum
    real(r8), intent(in) :: u2dot_sum, b2dot_sum
    real(r8), intent(in) :: u2dissip_sum, b2dissip_sum
    real(r8), intent(in) :: p_u_sum
    real(r8), intent(in) :: zp2_sum, zm2_sum
    real(r8), intent(in) :: bx0, by0, bz0

    integer , intent(in) :: nkpolar
    real(r8), intent(in) :: u2_bin (1:nkpolar), ux2_bin(1:nkpolar), uy2_bin(1:nkpolar), uz2_bin(1:nkpolar)
    real(r8), intent(in) :: b2_bin (1:nkpolar), bx2_bin(1:nkpolar), by2_bin(1:nkpolar), bz2_bin(1:nkpolar)
    real(r8), intent(in) :: zp2_bin(1:nkpolar), zm2_bin(1:nkpolar)

    integer, dimension (2) :: start2, count2

    ! output via NETCDF
    if(proc0) then
      ! total energy
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))
      status = nf90_put_var (ncid, u2_sum_id, u2_sum, start=(/nout/))
      status = nf90_put_var (ncid, b2_sum_id, b2_sum, start=(/nout/))
      status = nf90_put_var (ncid, u2dot_sum_id, u2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, b2dot_sum_id, b2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, u2dissip_sum_id, u2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, b2dissip_sum_id, b2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_u_sum_id, p_u_sum, start=(/nout/))
      status = nf90_put_var (ncid, zp2_sum_id, zp2_sum, start=(/nout/))
      status = nf90_put_var (ncid, zm2_sum_id, zm2_sum, start=(/nout/))
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
      status = nf90_put_var (ncid,  u2_bin_id  ,  u2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, ux2_bin_id  , ux2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, uy2_bin_id  , uy2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, uz2_bin_id  , uz2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid,  b2_bin_id  ,  b2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, bx2_bin_id  , bx2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, by2_bin_id  , by2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, bz2_bin_id  , bz2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, zp2_bin_id  , zp2_bin  , start=start2, count=count2)
      status = nf90_put_var (ncid, zm2_bin_id  , zm2_bin  , start=start2, count=count2)

      status = nf90_sync (ncid)

      nout = nout + 1
    endif
  end subroutine loop_io


!-----------------------------------------------!
!> @author  YK
!! @date    28 Jun 2021
!! @brief   Append variables to NETCDF
!           for cross section of fields
!-----------------------------------------------!
  subroutine loop_io_fields_section( &
                      ux_r_z0, ux_r_x0, ux_r_y0, &
                      uy_r_z0, uy_r_x0, uy_r_y0, &
                      uz_r_z0, uz_r_x0, uz_r_y0, &
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
                      u2_kxy, u2_kyz, u2_kxz, &
                      b2_kxy, b2_kyz, b2_kxz  &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz, nkx, nky, nkz
    use mp, only: proc0
    implicit none

    real(r8), intent(in) :: ux_r_z0(nlx, nly), ux_r_x0(nly, nlz), ux_r_y0(nlx, nlz)
    real(r8), intent(in) :: uy_r_z0(nlx, nly), uy_r_x0(nly, nlz), uy_r_y0(nlx, nlz)
    real(r8), intent(in) :: uz_r_z0(nlx, nly), uz_r_x0(nly, nlz), uz_r_y0(nlx, nlz)
                                                                                  
    real(r8), intent(in) :: wx_r_z0(nlx, nly), wx_r_x0(nly, nlz), wx_r_y0(nlx, nlz)
    real(r8), intent(in) :: wy_r_z0(nlx, nly), wy_r_x0(nly, nlz), wy_r_y0(nlx, nlz)
    real(r8), intent(in) :: wz_r_z0(nlx, nly), wz_r_x0(nly, nlz), wz_r_y0(nlx, nlz)
                                                                                  
    real(r8), intent(in) :: bx_r_z0(nlx, nly), bx_r_x0(nly, nlz), bx_r_y0(nlx, nlz)
    real(r8), intent(in) :: by_r_z0(nlx, nly), by_r_x0(nly, nlz), by_r_y0(nlx, nlz)
    real(r8), intent(in) :: bz_r_z0(nlx, nly), bz_r_x0(nly, nlz), bz_r_y0(nlx, nlz)
                                                                                  
    real(r8), intent(in) :: jx_r_z0(nlx, nly), jx_r_x0(nly, nlz), jx_r_y0(nlx, nlz)
    real(r8), intent(in) :: jy_r_z0(nlx, nly), jy_r_x0(nly, nlz), jy_r_y0(nlx, nlz)
    real(r8), intent(in) :: jz_r_z0(nlx, nly), jz_r_x0(nly, nlz), jz_r_y0(nlx, nlz)

    real(r8), intent(in) :: u2_kxy(nkx, nky), u2_kyz(nky, nkz), u2_kxz(nkx, nkz)
    real(r8), intent(in) :: b2_kxy(nkx, nky), b2_kyz(nky, nkz), b2_kxz(nkx, nkz)

    integer, dimension (3) :: start3, count3

    ! output via NETCDF
    if(proc0) then
      ! z=0 cut
      status = nf90_put_var (ncid_fld_section, tt_fld_section_id, tt, start=(/nout_fld_section/))
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_fld_section

      count3(1) = nlx
      count3(2) = nly
      count3(3) = 1
      status = nf90_put_var (ncid_fld_section, ux_r_z0_id, ux_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, uy_r_z0_id, uy_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, uz_r_z0_id, uz_r_z0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, wx_r_z0_id, wx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, wy_r_z0_id, wy_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, wz_r_z0_id, wz_r_z0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, bx_r_z0_id, bx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, by_r_z0_id, by_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, bz_r_z0_id, bz_r_z0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, jx_r_z0_id, jx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, jy_r_z0_id, jy_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, jz_r_z0_id, jz_r_z0, start=start3, count=count3)
      ! x=0 cut
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_fld_section

      count3(1) = nly
      count3(2) = nlz
      count3(3) = 1
      status = nf90_put_var (ncid_fld_section, ux_r_x0_id, ux_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, uy_r_x0_id, uy_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, uz_r_x0_id, uz_r_x0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, wx_r_x0_id, wx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, wy_r_x0_id, wy_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, wz_r_x0_id, wz_r_x0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, bx_r_x0_id, bx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, by_r_x0_id, by_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, bz_r_x0_id, bz_r_x0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, jx_r_x0_id, jx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, jy_r_x0_id, jy_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, jz_r_x0_id, jz_r_x0, start=start3, count=count3)
      ! y=0 cut
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_fld_section

      count3(1) = nlx
      count3(2) = nlz
      count3(3) = 1
      status = nf90_put_var (ncid_fld_section, ux_r_y0_id, ux_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, uy_r_y0_id, uy_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, uz_r_y0_id, uz_r_y0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, wx_r_y0_id, wx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, wy_r_y0_id, wy_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, wz_r_y0_id, wz_r_y0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, bx_r_y0_id, bx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, by_r_y0_id, by_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, bz_r_y0_id, bz_r_y0, start=start3, count=count3)
                                             
      status = nf90_put_var (ncid_fld_section, jx_r_y0_id, jx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, jy_r_y0_id, jy_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, jz_r_y0_id, jz_r_y0, start=start3, count=count3)

      ! kz sum
      status = nf90_put_var (ncid_fld_section, tt_fld_section_id, tt, start=(/nout_fld_section/))
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_fld_section

      count3(1) = nkx
      count3(2) = nky
      count3(3) = 1
      status = nf90_put_var (ncid_fld_section, u2_kxy_id, u2_kxy, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, b2_kxy_id, b2_kxy, start=start3, count=count3)
      ! kx sum
      status = nf90_put_var (ncid_fld_section, tt_fld_section_id, tt, start=(/nout_fld_section/))
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_fld_section

      count3(1) = nky
      count3(2) = nkz
      count3(3) = 1
      status = nf90_put_var (ncid_fld_section, u2_kyz_id, u2_kyz, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, b2_kyz_id, b2_kyz, start=start3, count=count3)
      ! ky sum
      status = nf90_put_var (ncid_fld_section, tt_fld_section_id, tt, start=(/nout_fld_section/))
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_fld_section

      count3(1) = nkx
      count3(2) = nkz
      count3(3) = 1
      status = nf90_put_var (ncid_fld_section, u2_kxz_id, u2_kxz, start=start3, count=count3)
      status = nf90_put_var (ncid_fld_section, b2_kxz_id, b2_kxz, start=start3, count=count3)

      status = nf90_sync (ncid_fld_section)

      nout_fld_section = nout_fld_section + 1
    endif
  end subroutine loop_io_fields_section


!-----------------------------------------------!
!> @author  YK
!! @date    3 Jan 2022
!! @brief   Append variables to NETCDF
!           for kpar
!-----------------------------------------------!
  subroutine loop_io_kpar(nkpolar, kpar_b, kpar_u, b1_ovr_b0)
    use time, only: tt
    use mp, only: proc0
    implicit none
    integer , intent(in) :: nkpolar
    real(r8), intent(in) :: kpar_b(1:nkpolar), kpar_u(1:nkpolar), b1_ovr_b0(1:nkpolar)

    integer, dimension (2) :: start2, count2

    ! output via NETCDF
    if(proc0) then
      status = nf90_put_var (ncid_kpar, tt_kpar_id, tt, start=(/nout_kpar/))
      start2(1) = 1
      start2(2) = nout_kpar

      count2(1) = nkpolar
      count2(2) = 1
      status = nf90_put_var (ncid_kpar, kpar_b_id   , kpar_b   , start=start2, count=count2)
      status = nf90_put_var (ncid_kpar, kpar_u_id   , kpar_u   , start=start2, count=count2)
      status = nf90_put_var (ncid_kpar, b1_ovr_b0_id, b1_ovr_b0, start=start2, count=count2)

      status = nf90_sync (ncid_kpar)

      nout_kpar = nout_kpar + 1
    endif
  end subroutine loop_io_kpar


!-----------------------------------------------!
!> @author  YK
!! @date    5 Jul 2021
!! @brief   Append variables to NETCDF
!           for 2nd order structure function
!-----------------------------------------------!
  subroutine loop_io_SF2(nl, sf2b, sf2u)
    use time, only: tt
    use mp, only: proc0
    implicit none
    integer, intent(in) :: nl
    real(r8), intent(in) :: sf2b(nl, nl), sf2u(nl, nl)

    integer, dimension (3) :: start3, count3

    ! output via NETCDF
    if(proc0) then
      status = nf90_put_var (ncid_SF2, tt_SF2_id, tt, start=(/nout_SF2/))
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_SF2

      count3(1) = nl
      count3(2) = nl
      count3(3) = 1
      status = nf90_put_var (ncid_SF2, SF2b_id, sf2b, start=start3, count=count3)
      status = nf90_put_var (ncid_SF2, SF2u_id, sf2u, start=start3, count=count3)

      status = nf90_sync (ncid_SF2)

      nout_SF2 = nout_SF2 + 1
    endif
  end subroutine loop_io_SF2


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Append field variables via MPIIO
!-----------------------------------------------!
  subroutine loop_io_fields_3D
    use fields, only: ux, uy, uz
    use fields, only: bx, by, bz
    use mp, only: proc0
    use time, only: tt
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mpiio, only: mpiio_write_var
    use shearing_box, only: tsc
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

    call mpiio_write_var(fh_ux, disp_ux, sizes, subsizes, starts, ux)
    call mpiio_write_var(fh_uy, disp_uy, sizes, subsizes, starts, uy)
    call mpiio_write_var(fh_uz, disp_uz, sizes, subsizes, starts, uz)
    call mpiio_write_var(fh_bx, disp_bx, sizes, subsizes, starts, bx)
    call mpiio_write_var(fh_by, disp_by, sizes, subsizes, starts, by)
    call mpiio_write_var(fh_bz, disp_bz, sizes, subsizes, starts, bz)

    if(proc0) then
      write (unit=field_time_unit, fmt="(100es30.21)") tt, tsc
      flush (field_time_unit)
    endif
  end subroutine loop_io_fields_3D


!-----------------------------------------------!
!> @author  YK
!! @date    16 Jan 2019
!! @brief   Save restart file via MPIIO
!-----------------------------------------------!
  subroutine save_restart
    use fields, only: ux, uy, uz, p
    use fields, only: bx, by, bz
    use mp, only: proc0
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: tt, dt
    use params, only: restart_dir
    use file, only: open_output_file, close_file
    use mpiio, only: mpiio_write_one
    use shearing_box, only: tsc
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

    call mpiio_write_one(ux, sizes, subsizes, starts, trim(restart_dir)//'ux.dat')
    call mpiio_write_one(uy, sizes, subsizes, starts, trim(restart_dir)//'uy.dat')
    call mpiio_write_one(uz, sizes, subsizes, starts, trim(restart_dir)//'uz.dat')
    call mpiio_write_one(bx, sizes, subsizes, starts, trim(restart_dir)//'bx.dat')
    call mpiio_write_one(by, sizes, subsizes, starts, trim(restart_dir)//'by.dat')
    call mpiio_write_one(bz, sizes, subsizes, starts, trim(restart_dir)//'bz.dat')
    call mpiio_write_one( p, sizes, subsizes, starts, trim(restart_dir)//'p.dat' )

    if(proc0) then
      call open_output_file (time_unit, trim(restart_dir)//'time.dat')
      write (unit=time_unit, fmt="(3X, 'tt', 28X, 'tst')")
      write (unit=time_unit, fmt="(100es30.21)") tt, tsc
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

    call MPI_FILE_CLOSE(fh_ux,ierr)
    call MPI_FILE_CLOSE(fh_uy,ierr)
    call MPI_FILE_CLOSE(fh_uz,ierr)
    call MPI_FILE_CLOSE(fh_bx,ierr)
    call MPI_FILE_CLOSE(fh_by,ierr)
    call MPI_FILE_CLOSE(fh_bz,ierr)
    if(proc0) then
      call close_file (field_time_unit)
    endif

    if(proc0) then
      status = nf90_close (ncid)
    endif

  end subroutine finish_io

end module io


