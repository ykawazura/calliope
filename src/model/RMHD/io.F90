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

  public :: init_io, finish_io, loop_io, loop_io_2D, loop_io_3D, save_restart

  private

  ! MPIIO
  integer :: fh_phi, fh_psi
  character(len=100) :: filename
  integer (kind=MPI_OFFSET_KIND) :: disp_phi
  integer (kind=MPI_OFFSET_KIND) :: disp_psi
  integer :: field_time_unit

  ! NETCDF for regular output file
  integer :: status
  integer, parameter :: kind_nf = kind (NF90_NOERR)
  integer (kind_nf) :: ncid
  integer :: run_id
  integer (kind_nf) :: char10_dim
  ! parameter
  integer :: nupe_x_id, nupe_x_exp_id, nupe_z_id, nupe_z_exp_id
  integer :: etape_x_id, etape_x_exp_id, etape_z_id, etape_z_exp_id
  ! coordinate
  integer :: xx_id, yy_id, zz_id, kx_id, ky_id, kz_id, kpbin_id, tt_id
  ! total energy
  integer :: upe2_sum_id, bpe2_sum_id
  integer :: upe2dot_sum_id, bpe2dot_sum_id
  integer :: upe2dissip_sum_id, bpe2dissip_sum_id
  integer :: p_omg_sum_id, p_psi_sum_id
  ! polar spectrum
  integer :: upe2_bin_id, bpe2_bin_id

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim
  integer, dimension (3) :: bin_dim, z0_dim, x0_dim, y0_dim

  integer :: nout

  ! NETCDF for cross section output file
  integer (kind_nf) :: ncid_2D
  integer :: xx_2D_id, yy_2D_id, zz_2D_id, tt_2D_id

  integer :: phi_r_z0_id, phi_r_x0_id, phi_r_y0_id
  integer :: psi_r_z0_id, psi_r_x0_id, psi_r_y0_id
  integer :: omg_r_z0_id, omg_r_x0_id, omg_r_y0_id
  integer :: jpa_r_z0_id, jpa_r_x0_id, jpa_r_y0_id
  integer ::  ux_r_z0_id,  ux_r_x0_id,  ux_r_y0_id
  integer ::  uy_r_z0_id,  uy_r_x0_id,  uy_r_y0_id
  integer ::  bx_r_z0_id,  bx_r_x0_id,  bx_r_y0_id
  integer ::  by_r_z0_id,  by_r_x0_id,  by_r_y0_id

  integer (kind_nf) :: xx_2D_dim, yy_2D_dim, zz_2D_dim, tt_2D_dim

  integer :: nout_2D

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
    use params, only: restart_dir
    implicit none
    integer :: ierr

    filename = 'phi.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_phi, ierr)
    call MPI_FILE_SET_SIZE(fh_phi, 0_MPI_OFFSET_KIND, ierr)  ! guarantee overwriting
    disp_phi = 0_MPI_OFFSET_KIND

    filename = 'psi.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_psi, ierr)
    call MPI_FILE_SET_SIZE(fh_psi, 0_MPI_OFFSET_KIND, ierr)  ! guarantee overwriting
    disp_psi = 0_MPI_OFFSET_KIND

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
    use params, only: runname, &
                      nupe_x , nupe_x_exp , nupe_z , nupe_z_exp, &
                      etape_x, etape_x_exp, etape_z, etape_z_exp
    implicit none
    integer , intent(in) :: nkpolar
    real(r8), intent(in) :: kpbin(1:nkpolar)

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

      status = nf90_def_var (ncid, 'nupe_x', NF90_DOUBLE, nupe_x_id)
      status = nf90_def_var (ncid, 'nupe_x_exp', NF90_DOUBLE, nupe_x_exp_id)
      status = nf90_def_var (ncid, 'nupe_z', NF90_DOUBLE, nupe_z_id)
      status = nf90_def_var (ncid, 'nupe_z_exp', NF90_DOUBLE, nupe_z_exp_id)
      status = nf90_def_var (ncid, 'etape_x', NF90_DOUBLE, etape_x_id)
      status = nf90_def_var (ncid, 'etape_x_exp', NF90_DOUBLE, etape_x_exp_id)
      status = nf90_def_var (ncid, 'etape_z', NF90_DOUBLE, etape_z_id)
      status = nf90_def_var (ncid, 'etape_z_exp', NF90_DOUBLE, etape_z_exp_id)

      status = nf90_def_var (ncid, 'xx', NF90_DOUBLE, xx_dim, xx_id)
      status = nf90_def_var (ncid, 'yy', NF90_DOUBLE, yy_dim, yy_id)
      status = nf90_def_var (ncid, 'zz', NF90_DOUBLE, zz_dim, zz_id)
      status = nf90_def_var (ncid, 'kx', NF90_DOUBLE, kx_dim, kx_id)
      status = nf90_def_var (ncid, 'ky', NF90_DOUBLE, ky_dim, ky_id)
      status = nf90_def_var (ncid, 'kz', NF90_DOUBLE, kz_dim, kz_id)
      status = nf90_def_var (ncid, 'kpbin', NF90_DOUBLE, kpbin_dim, kpbin_id)
      status = nf90_def_var (ncid, 'tt', NF90_DOUBLE, tt_dim, tt_id)
      ! total energy
      status = nf90_def_var (ncid, 'upe2_sum', NF90_DOUBLE, tt_dim, upe2_sum_id)
      status = nf90_def_var (ncid, 'bpe2_sum', NF90_DOUBLE, tt_dim, bpe2_sum_id)
      status = nf90_def_var (ncid, 'upe2dot_sum', NF90_DOUBLE, tt_dim, upe2dot_sum_id)
      status = nf90_def_var (ncid, 'bpe2dot_sum', NF90_DOUBLE, tt_dim, bpe2dot_sum_id)
      status = nf90_def_var (ncid, 'upe2dissip_sum', NF90_DOUBLE, tt_dim, upe2dissip_sum_id)
      status = nf90_def_var (ncid, 'bpe2dissip_sum', NF90_DOUBLE, tt_dim, bpe2dissip_sum_id)
      status = nf90_def_var (ncid, 'p_omg_sum'     , NF90_DOUBLE, tt_dim, p_omg_sum_id   )
      status = nf90_def_var (ncid, 'p_psi_sum'     , NF90_DOUBLE, tt_dim, p_psi_sum_id   )
      ! polar spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = kz_dim
      bin_dim (3) = tt_dim
      status = nf90_def_var (ncid, 'upe2_bin', NF90_DOUBLE, bin_dim, upe2_bin_id)
      status = nf90_def_var (ncid, 'bpe2_bin', NF90_DOUBLE, bin_dim, bpe2_bin_id)

      status = nf90_enddef (ncid)  ! out of definition mode

      status = nf90_put_var (ncid, nupe_x_id, nupe_x)
      status = nf90_put_var (ncid, nupe_x_exp_id, dble(nupe_x_exp))
      status = nf90_put_var (ncid, nupe_z_id, nupe_z)
      status = nf90_put_var (ncid, nupe_z_exp_id, dble(nupe_z_exp))
      status = nf90_put_var (ncid, etape_x_id, etape_x)
      status = nf90_put_var (ncid, etape_x_exp_id, dble(etape_x_exp))
      status = nf90_put_var (ncid, etape_z_id, etape_z)
      status = nf90_put_var (ncid, etape_z_exp_id, dble(etape_z_exp))

      status = nf90_put_var (ncid, xx_id, xx)
      status = nf90_put_var (ncid, yy_id, yy)
      status = nf90_put_var (ncid, zz_id, zz)
      status = nf90_put_var (ncid, kx_id, kx)
      status = nf90_put_var (ncid, ky_id, ky)
      status = nf90_put_var (ncid, kz_id, kz)
      status = nf90_put_var (ncid, kpbin_id, kpbin)

      nout = 1

      !--------------------------------------------------!
      ! Output for 2D cross sections of fields
      !--------------------------------------------------!
      filename = trim(runname)//'.out.2D.nc' ! File name
      status = nf90_create (filename, NF90_CLOBBER, ncid_2D)

      status = nf90_put_att (ncid_2D, NF90_GLOBAL, 'title', 'calliope simulation data')
      status = nf90_def_dim (ncid_2D, 'char10', 10, char10_dim)
      status = nf90_def_var (ncid_2D, 'run_info', NF90_CHAR, char10_dim, run_id)
      status = nf90_put_att (ncid_2D, run_id, 'model', _MODEL_)

      status = nf90_def_dim (ncid_2D, 'xx', size(xx), xx_2D_dim)
      status = nf90_def_dim (ncid_2D, 'yy', size(yy), yy_2D_dim)
      status = nf90_def_dim (ncid_2D, 'zz', size(zz), zz_2D_dim)
      status = nf90_def_dim (ncid_2D, 'tt', NF90_UNLIMITED, tt_2D_dim)

      status = nf90_def_var (ncid_2D, 'xx', NF90_DOUBLE, xx_2D_dim, xx_2D_id)
      status = nf90_def_var (ncid_2D, 'yy', NF90_DOUBLE, yy_2D_dim, yy_2D_id)
      status = nf90_def_var (ncid_2D, 'zz', NF90_DOUBLE, zz_2D_dim, zz_2D_id)
      status = nf90_def_var (ncid_2D, 'tt', NF90_DOUBLE, tt_2D_dim, tt_2D_id)

      z0_dim (1) = xx_2D_dim
      z0_dim (2) = yy_2D_dim
      z0_dim (3) = tt_2D_dim

      x0_dim (1) = yy_2D_dim
      x0_dim (2) = zz_2D_dim
      x0_dim (3) = tt_2D_dim

      y0_dim (1) = xx_2D_dim
      y0_dim (2) = zz_2D_dim
      y0_dim (3) = tt_2D_dim
      status = nf90_def_var (ncid_2D, 'phi_r_z0', NF90_DOUBLE, z0_dim, phi_r_z0_id)
      status = nf90_def_var (ncid_2D, 'phi_r_x0', NF90_DOUBLE, x0_dim, phi_r_x0_id)
      status = nf90_def_var (ncid_2D, 'phi_r_y0', NF90_DOUBLE, y0_dim, phi_r_y0_id)
      status = nf90_def_var (ncid_2D, 'psi_r_z0', NF90_DOUBLE, z0_dim, psi_r_z0_id)
      status = nf90_def_var (ncid_2D, 'psi_r_x0', NF90_DOUBLE, x0_dim, psi_r_x0_id)
      status = nf90_def_var (ncid_2D, 'psi_r_y0', NF90_DOUBLE, y0_dim, psi_r_y0_id)
      status = nf90_def_var (ncid_2D, 'omg_r_z0', NF90_DOUBLE, z0_dim, omg_r_z0_id)
      status = nf90_def_var (ncid_2D, 'omg_r_x0', NF90_DOUBLE, x0_dim, omg_r_x0_id)
      status = nf90_def_var (ncid_2D, 'omg_r_y0', NF90_DOUBLE, y0_dim, omg_r_y0_id)
      status = nf90_def_var (ncid_2D, 'jpa_r_z0', NF90_DOUBLE, z0_dim, jpa_r_z0_id)
      status = nf90_def_var (ncid_2D, 'jpa_r_x0', NF90_DOUBLE, x0_dim, jpa_r_x0_id)
      status = nf90_def_var (ncid_2D, 'jpa_r_y0', NF90_DOUBLE, y0_dim, jpa_r_y0_id)
      status = nf90_def_var (ncid_2D,  'ux_r_z0', NF90_DOUBLE, z0_dim,  ux_r_z0_id)
      status = nf90_def_var (ncid_2D,  'ux_r_x0', NF90_DOUBLE, x0_dim,  ux_r_x0_id)
      status = nf90_def_var (ncid_2D,  'ux_r_y0', NF90_DOUBLE, y0_dim,  ux_r_y0_id)
      status = nf90_def_var (ncid_2D,  'uy_r_z0', NF90_DOUBLE, z0_dim,  uy_r_z0_id)
      status = nf90_def_var (ncid_2D,  'uy_r_x0', NF90_DOUBLE, x0_dim,  uy_r_x0_id)
      status = nf90_def_var (ncid_2D,  'uy_r_y0', NF90_DOUBLE, y0_dim,  uy_r_y0_id)
      status = nf90_def_var (ncid_2D,  'bx_r_z0', NF90_DOUBLE, z0_dim,  bx_r_z0_id)
      status = nf90_def_var (ncid_2D,  'bx_r_x0', NF90_DOUBLE, x0_dim,  bx_r_x0_id)
      status = nf90_def_var (ncid_2D,  'bx_r_y0', NF90_DOUBLE, y0_dim,  bx_r_y0_id)
      status = nf90_def_var (ncid_2D,  'by_r_z0', NF90_DOUBLE, z0_dim,  by_r_z0_id)
      status = nf90_def_var (ncid_2D,  'by_r_x0', NF90_DOUBLE, x0_dim,  by_r_x0_id)
      status = nf90_def_var (ncid_2D,  'by_r_y0', NF90_DOUBLE, y0_dim,  by_r_y0_id)

      status = nf90_enddef (ncid_2D)  ! out of definition mode

      status = nf90_put_var (ncid_2D, xx_2D_id, xx)
      status = nf90_put_var (ncid_2D, yy_2D_id, yy)
      status = nf90_put_var (ncid_2D, zz_2D_id, zz)

      nout_2D = 1
    endif
  end subroutine init_io_netcdf


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Append variables to NETCDF
!           for params, time history & spectra
!-----------------------------------------------!
  subroutine loop_io( &
                      upe2_sum, bpe2_sum, &
                      upe2dot_sum, bpe2dot_sum, &
                      upe2dissip_sum, bpe2dissip_sum, &
                      p_omg_sum, p_psi_sum, &
                      !
                      nkpolar, &
                      upe2_bin, bpe2_bin &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz, nkz
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_2D
    implicit none
    real(r8), intent(in) :: upe2_sum, bpe2_sum
    real(r8), intent(in) :: upe2dot_sum, bpe2dot_sum
    real(r8), intent(in) :: upe2dissip_sum, bpe2dissip_sum
    real(r8), intent(in) :: p_omg_sum, p_psi_sum

    integer, intent(in) :: nkpolar
    real(r8), intent(in) :: upe2_bin(1:nkpolar, nkz), bpe2_bin(1:nkpolar, nkz)

    integer, dimension (3) :: start3, count3

    if (proc0) call put_time_stamp(timer_io_total)

    ! output via NETCDF
    if(proc0) then
      ! total energy
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))
      status = nf90_put_var (ncid, upe2_sum_id, upe2_sum, start=(/nout/))
      status = nf90_put_var (ncid, bpe2_sum_id, bpe2_sum, start=(/nout/))
      status = nf90_put_var (ncid, upe2dot_sum_id, upe2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, bpe2dot_sum_id, bpe2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, upe2dissip_sum_id, upe2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, bpe2dissip_sum_id, bpe2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_omg_sum_id, p_omg_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_psi_sum_id, p_psi_sum, start=(/nout/))
      ! polar spectrum
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout

      count3(1) = nkpolar
      count3(2) = nkz
      count3(3) = 1
      status = nf90_put_var (ncid, upe2_bin_id, upe2_bin, start=start3, count=count3)
      status = nf90_put_var (ncid, bpe2_bin_id, bpe2_bin, start=start3, count=count3)

      status = nf90_sync (ncid)

      nout = nout + 1
    endif

    if (proc0) call put_time_stamp(timer_io_total)
  end subroutine loop_io


!-----------------------------------------------!
!> @author  YK
!! @date    28 Jun 2021
!! @brief   Append variables to NETCDF
!           for cross section of fields
!-----------------------------------------------!
  subroutine loop_io_2D( &
                      phi_r_z0, phi_r_x0, phi_r_y0, &
                      psi_r_z0, psi_r_x0, psi_r_y0, &
                      omg_r_z0, omg_r_x0, omg_r_y0, &
                      jpa_r_z0, jpa_r_x0, jpa_r_y0, &
                       ux_r_z0,  ux_r_x0,  ux_r_y0, &
                       uy_r_z0,  uy_r_x0,  uy_r_y0, &
                       bx_r_z0,  bx_r_x0,  bx_r_y0, &
                       by_r_z0,  by_r_x0,  by_r_y0  &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_2D
    implicit none

    real(r8), intent(in) :: phi_r_z0(nlx, nly), phi_r_x0(nly, nlz), phi_r_y0(nlx, nlz)
    real(r8), intent(in) :: psi_r_z0(nlx, nly), psi_r_x0(nly, nlz), psi_r_y0(nlx, nlz)
    real(r8), intent(in) :: omg_r_z0(nlx, nly), omg_r_x0(nly, nlz), omg_r_y0(nlx, nlz)
    real(r8), intent(in) :: jpa_r_z0(nlx, nly), jpa_r_x0(nly, nlz), jpa_r_y0(nlx, nlz)
    real(r8), intent(in) ::  ux_r_z0(nlx, nly),  ux_r_x0(nly, nlz),  ux_r_y0(nlx, nlz)
    real(r8), intent(in) ::  uy_r_z0(nlx, nly),  uy_r_x0(nly, nlz),  uy_r_y0(nlx, nlz)
    real(r8), intent(in) ::  bx_r_z0(nlx, nly),  bx_r_x0(nly, nlz),  bx_r_y0(nlx, nlz)
    real(r8), intent(in) ::  by_r_z0(nlx, nly),  by_r_x0(nly, nlz),  by_r_y0(nlx, nlz)

    integer, dimension (3) :: start3, count3

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)

    ! output via NETCDF
    if(proc0) then
      ! z=0 cut
      status = nf90_put_var (ncid_2D, tt_2D_id, tt, start=(/nout_2D/))
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_2D

      count3(1) = nlx
      count3(2) = nly
      count3(3) = 1
      status = nf90_put_var (ncid_2D, phi_r_z0_id, phi_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, psi_r_z0_id, psi_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, omg_r_z0_id, omg_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, jpa_r_z0_id, jpa_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  ux_r_z0_id,  ux_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  uy_r_z0_id,  uy_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  bx_r_z0_id,  bx_r_z0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  by_r_z0_id,  by_r_z0, start=start3, count=count3)
      ! x=0 cut
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_2D

      count3(1) = nly
      count3(2) = nlz
      count3(3) = 1
      status = nf90_put_var (ncid_2D, phi_r_x0_id, phi_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, psi_r_x0_id, psi_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, omg_r_x0_id, omg_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, jpa_r_x0_id, jpa_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  ux_r_x0_id,  ux_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  uy_r_x0_id,  uy_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  bx_r_x0_id,  bx_r_x0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  by_r_x0_id,  by_r_x0, start=start3, count=count3)
      ! y=0 cut
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout_2D

      count3(1) = nlx
      count3(2) = nlz
      count3(3) = 1
      status = nf90_put_var (ncid_2D, phi_r_y0_id, phi_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, psi_r_y0_id, psi_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, omg_r_y0_id, omg_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D, jpa_r_y0_id, jpa_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  ux_r_y0_id,  ux_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  uy_r_y0_id,  uy_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  bx_r_y0_id,  bx_r_y0, start=start3, count=count3)
      status = nf90_put_var (ncid_2D,  by_r_y0_id,  by_r_y0, start=start3, count=count3)

      status = nf90_sync (ncid_2D)

      nout_2D = nout_2D + 1
    endif

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)
  end subroutine loop_io_2D


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Append field variables via MPIIO
!-----------------------------------------------!
  subroutine loop_io_3D
    use fields, only: phi, psi
    use mp, only: proc0
    use time, only: tt
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use mpiio, only: mpiio_write_var
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_3D
    implicit none
    integer, dimension(3) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_3D)

    sizes(1) = nkx
    sizes(2) = nkz
    sizes(3) = nky
    subsizes(1) = ikx_en - ikx_st + 1
    subsizes(2) = ikz_en - ikz_st + 1
    subsizes(3) = iky_en - iky_st + 1
    starts(1) = ikx_st - 1
    starts(2) = ikz_st - 1
    starts(3) = iky_st - 1

    call mpiio_write_var(fh_phi, disp_phi, sizes, subsizes, starts, phi)
    call mpiio_write_var(fh_psi, disp_psi, sizes, subsizes, starts, psi)

    if(proc0) then
      write (unit=field_time_unit, fmt="(100es30.21)") tt
      flush (field_time_unit)
    endif

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_3D)
  end subroutine loop_io_3D


!-----------------------------------------------!
!> @author  YK
!! @date    16 Jan 2019
!! @brief   Save restart file via MPIIO
!-----------------------------------------------!
  subroutine save_restart
    use fields, only: phi, omg, psi
    use mp, only: proc0
    use grid, only: nkx, nky, nkz
    use grid, only: ikx_st, iky_st, ikz_st, ikx_en, iky_en, ikz_en
    use time, only: tt, dt
    use params, only: restart_dir
    use file, only: open_output_file, close_file
    use mpiio, only: mpiio_write_one
    use time_stamp, only: put_time_stamp, timer_save_restart
    implicit none
    integer :: time_unit
    integer, dimension(3) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_save_restart)

    sizes(1) = nkx
    sizes(2) = nkz
    sizes(3) = nky
    subsizes(1) = ikx_en - ikx_st + 1
    subsizes(2) = ikz_en - ikz_st + 1
    subsizes(3) = iky_en - iky_st + 1
    starts(1) = ikx_st - 1
    starts(2) = ikz_st - 1
    starts(3) = iky_st - 1

    call mpiio_write_one(phi, sizes, subsizes, starts, trim(restart_dir)//'phi.dat')
    call mpiio_write_one(omg, sizes, subsizes, starts, trim(restart_dir)//'omg.dat')
    call mpiio_write_one(psi, sizes, subsizes, starts, trim(restart_dir)//'psi.dat')

    if(proc0) then
      call open_output_file (time_unit, trim(restart_dir)//'time.dat')
      write (unit=time_unit, fmt="(3X, 'tt')")
      write (unit=time_unit, fmt="(100es30.21)") tt
      call close_file (time_unit)
    endif

    if (proc0) call put_time_stamp(timer_save_restart)
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

    call MPI_FILE_CLOSE(fh_phi,ierr)
    call MPI_FILE_CLOSE(fh_psi,ierr)
    if(proc0) then
      call close_file (field_time_unit)
    endif

    if(proc0) then
      status = nf90_close (ncid)
    endif

  end subroutine finish_io

end module io

