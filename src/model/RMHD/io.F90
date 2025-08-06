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
  integer :: out3d_time_unit

  ! MPIIO for 2D
  integer :: fh_phi_r_z0, fh_phi_r_x0, fh_phi_r_y0
  integer :: fh_psi_r_z0, fh_psi_r_x0, fh_psi_r_y0
  integer :: fh_omg_r_z0, fh_omg_r_x0, fh_omg_r_y0
  integer :: fh_jpa_r_z0, fh_jpa_r_x0, fh_jpa_r_y0
  integer :: fh_ux_r_z0 , fh_ux_r_x0 , fh_ux_r_y0, &
             fh_uy_r_z0 , fh_uy_r_x0 , fh_uy_r_y0
  integer :: fh_bx_r_z0 , fh_bx_r_x0 , fh_bx_r_y0, &
             fh_by_r_z0 , fh_by_r_x0 , fh_by_r_y0

  integer (kind=MPI_OFFSET_KIND) :: disp_phi_r_z0, disp_phi_r_x0, disp_phi_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_psi_r_z0, disp_psi_r_x0, disp_psi_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_omg_r_z0, disp_omg_r_x0, disp_omg_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_jpa_r_z0, disp_jpa_r_x0, disp_jpa_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_ux_r_z0, disp_ux_r_x0, disp_ux_r_y0, &
                                    disp_uy_r_z0, disp_uy_r_x0, disp_uy_r_y0
  integer (kind=MPI_OFFSET_KIND) :: disp_bx_r_z0, disp_bx_r_x0, disp_bx_r_y0, &
                                    disp_by_r_z0, disp_by_r_x0, disp_by_r_y0
  integer :: out2d_time_unit

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
  integer :: upe2_sum_id , bpe2_sum_id
  integer :: zppe2_sum_id, zmpe2_sum_id
  integer :: upe2dot_sum_id, bpe2dot_sum_id
  integer :: upe2dissip_sum_id, bpe2dissip_sum_id
  integer :: p_phi_sum_id, p_psi_sum_id, p_xhl_sum_id
  ! polar spectrum
  integer :: upe2_bin_id , bpe2_bin_id
  integer :: zppe2_bin_id, zmpe2_bin_id
  integer :: hp_bin_id, hm_bin_id

  integer (kind_nf) :: xx_dim, yy_dim, zz_dim, kx_dim, ky_dim, kz_dim, kpbin_dim, tt_dim
  integer, dimension (3) :: bin_dim

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

    call init_io_decomp_2d
    call init_io_decomp_3d
    call init_io_netcdf(nkpolar, kpbin)
  end subroutine init_io


!-----------------------------------------------!
!> @author  YK
!! @date    29 Dec 2018
!! @brief   Initialization of MPIIO for 3D
!-----------------------------------------------!
  subroutine init_io_decomp_3d
    use mp, only: proc0
    use file, only: open_output_file
    implicit none

    call set_file_handle('out3d/phi.dat', fh_phi, disp_phi)
    call set_file_handle('out3d/psi.dat', fh_psi, disp_psi)

    if(proc0) then
      call open_output_file (out3d_time_unit, 'out3d/time.dat')
    endif
  end subroutine init_io_decomp_3d


!-----------------------------------------------!
!> @author  YK
!! @date    7 May 2022
!! @brief   Initialization of MPIIO for 2D
!-----------------------------------------------!
  subroutine init_io_decomp_2d
    use mp, only: proc0
    use file, only: open_output_file
    implicit none

    call set_file_handle('out2d/phi_r_z0.dat', fh_phi_r_z0, disp_phi_r_z0)
    call set_file_handle('out2d/phi_r_x0.dat', fh_phi_r_x0, disp_phi_r_x0)
    call set_file_handle('out2d/phi_r_y0.dat', fh_phi_r_y0, disp_phi_r_y0)

    call set_file_handle('out2d/psi_r_z0.dat', fh_psi_r_z0, disp_psi_r_z0)
    call set_file_handle('out2d/psi_r_x0.dat', fh_psi_r_x0, disp_psi_r_x0)
    call set_file_handle('out2d/psi_r_y0.dat', fh_psi_r_y0, disp_psi_r_y0)

    call set_file_handle('out2d/omg_r_z0.dat', fh_omg_r_z0, disp_omg_r_z0)
    call set_file_handle('out2d/omg_r_x0.dat', fh_omg_r_x0, disp_omg_r_x0)
    call set_file_handle('out2d/omg_r_y0.dat', fh_omg_r_y0, disp_omg_r_y0)

    call set_file_handle('out2d/jpa_r_z0.dat', fh_jpa_r_z0, disp_jpa_r_z0)
    call set_file_handle('out2d/jpa_r_x0.dat', fh_jpa_r_x0, disp_jpa_r_x0)
    call set_file_handle('out2d/jpa_r_y0.dat', fh_jpa_r_y0, disp_jpa_r_y0)

    call set_file_handle('out2d/ux_r_z0.dat' , fh_ux_r_z0 , disp_ux_r_z0 )
    call set_file_handle('out2d/ux_r_x0.dat' , fh_ux_r_x0 , disp_ux_r_x0 )
    call set_file_handle('out2d/ux_r_y0.dat' , fh_ux_r_y0 , disp_ux_r_y0 )
                                                                         
    call set_file_handle('out2d/uy_r_z0.dat' , fh_uy_r_z0 , disp_uy_r_z0 )
    call set_file_handle('out2d/uy_r_x0.dat' , fh_uy_r_x0 , disp_uy_r_x0 )
    call set_file_handle('out2d/uy_r_y0.dat' , fh_uy_r_y0 , disp_uy_r_y0 )
                                                                         
    call set_file_handle('out2d/bx_r_z0.dat' , fh_bx_r_z0 , disp_bx_r_z0 )
    call set_file_handle('out2d/bx_r_x0.dat' , fh_bx_r_x0 , disp_bx_r_x0 )
    call set_file_handle('out2d/bx_r_y0.dat' , fh_bx_r_y0 , disp_bx_r_y0 )
                                                                         
    call set_file_handle('out2d/by_r_z0.dat' , fh_by_r_z0 , disp_by_r_z0 )
    call set_file_handle('out2d/by_r_x0.dat' , fh_by_r_x0 , disp_by_r_x0 )
    call set_file_handle('out2d/by_r_y0.dat' , fh_by_r_y0 , disp_by_r_y0 )


    if(proc0) then
      call open_output_file (out2d_time_unit, 'out2d/time.dat')
    endif
  end subroutine init_io_decomp_2d


!-----------------------------------------------!
!> @author  YK
!! @date    20 Jun 2022
!! @brief   Set file handle for MPIIO
!-----------------------------------------------!
  subroutine set_file_handle(fn, fh, disp)
    implicit none
    character(*) :: fn
    integer :: fh
    integer (kind=MPI_OFFSET_KIND) :: disp
    integer :: ierr

    call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(fn), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, ierr)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
  end subroutine set_file_handle


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
      status = nf90_def_var (ncid, 'upe2_sum' , NF90_DOUBLE, tt_dim, upe2_sum_id )
      status = nf90_def_var (ncid, 'bpe2_sum' , NF90_DOUBLE, tt_dim, bpe2_sum_id )
      status = nf90_def_var (ncid, 'zppe2_sum', NF90_DOUBLE, tt_dim, zppe2_sum_id)
      status = nf90_def_var (ncid, 'zmpe2_sum', NF90_DOUBLE, tt_dim, zmpe2_sum_id)
      status = nf90_def_var (ncid, 'upe2dot_sum', NF90_DOUBLE, tt_dim, upe2dot_sum_id)
      status = nf90_def_var (ncid, 'bpe2dot_sum', NF90_DOUBLE, tt_dim, bpe2dot_sum_id)
      status = nf90_def_var (ncid, 'upe2dissip_sum', NF90_DOUBLE, tt_dim, upe2dissip_sum_id)
      status = nf90_def_var (ncid, 'bpe2dissip_sum', NF90_DOUBLE, tt_dim, bpe2dissip_sum_id)
      status = nf90_def_var (ncid, 'p_phi_sum'     , NF90_DOUBLE, tt_dim, p_phi_sum_id   )
      status = nf90_def_var (ncid, 'p_psi_sum'     , NF90_DOUBLE, tt_dim, p_psi_sum_id   )
      status = nf90_def_var (ncid, 'p_xhl_sum'     , NF90_DOUBLE, tt_dim, p_xhl_sum_id   )
      ! polar spectrum
      bin_dim (1) = kpbin_dim
      bin_dim (2) = kz_dim
      bin_dim (3) = tt_dim
      status = nf90_def_var (ncid, 'upe2_bin' , NF90_DOUBLE, bin_dim, upe2_bin_id )
      status = nf90_def_var (ncid, 'bpe2_bin' , NF90_DOUBLE, bin_dim, bpe2_bin_id )
      status = nf90_def_var (ncid, 'zppe2_bin', NF90_DOUBLE, bin_dim, zppe2_bin_id)
      status = nf90_def_var (ncid, 'zmpe2_bin', NF90_DOUBLE, bin_dim, zmpe2_bin_id)
      status = nf90_def_var (ncid, 'hp_bin'   , NF90_DOUBLE, bin_dim, hp_bin_id   )
      status = nf90_def_var (ncid, 'hm_bin'   , NF90_DOUBLE, bin_dim, hm_bin_id   )

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
                      p_phi_sum, p_psi_sum, p_xhl_sum, &
                      zppe2_sum, zmpe2_sum, &
                      !
                      nkpolar, &
                      upe2_bin , bpe2_bin , &
                      zppe2_bin, zmpe2_bin, &
                      hp_bin   , hm_bin     &
                    )
    use time, only: tt
    use grid, only: nlx, nly, nlz, nkz
    use mp, only: proc0
    use time_stamp, only: put_time_stamp, timer_io_total
    implicit none
    real(r8), intent(in) :: upe2_sum , bpe2_sum
    real(r8), intent(in) :: zppe2_sum, zmpe2_sum
    real(r8), intent(in) :: upe2dot_sum, bpe2dot_sum
    real(r8), intent(in) :: upe2dissip_sum, bpe2dissip_sum
    real(r8), intent(in) :: p_phi_sum, p_psi_sum, p_xhl_sum

    integer, intent(in) :: nkpolar
    real(r8), intent(in) :: upe2_bin (1:nkpolar, nkz), bpe2_bin (1:nkpolar, nkz)
    real(r8), intent(in) :: zppe2_bin(1:nkpolar, nkz), zmpe2_bin(1:nkpolar, nkz)
    real(r8), intent(in) :: hp_bin   (1:nkpolar, nkz), hm_bin   (1:nkpolar, nkz)

    integer, dimension (3) :: start3, count3

    if (proc0) call put_time_stamp(timer_io_total)

    ! output via NETCDF
    if(proc0) then
      ! total energy
      status = nf90_put_var (ncid, tt_id, tt, start=(/nout/))
      status = nf90_put_var (ncid, upe2_sum_id , upe2_sum , start=(/nout/))
      status = nf90_put_var (ncid, bpe2_sum_id , bpe2_sum , start=(/nout/))
      status = nf90_put_var (ncid, zppe2_sum_id, zppe2_sum, start=(/nout/))
      status = nf90_put_var (ncid, zmpe2_sum_id, zmpe2_sum, start=(/nout/))
      status = nf90_put_var (ncid, upe2dot_sum_id, upe2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, bpe2dot_sum_id, bpe2dot_sum, start=(/nout/))
      status = nf90_put_var (ncid, upe2dissip_sum_id, upe2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, bpe2dissip_sum_id, bpe2dissip_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_phi_sum_id, p_phi_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_psi_sum_id, p_psi_sum, start=(/nout/))
      status = nf90_put_var (ncid, p_xhl_sum_id, p_xhl_sum, start=(/nout/))
      ! polar spectrum
      start3(1) = 1
      start3(2) = 1
      start3(3) = nout

      count3(1) = nkpolar
      count3(2) = nkz
      count3(3) = 1
      status = nf90_put_var (ncid, upe2_bin_id , upe2_bin , start=start3, count=count3)
      status = nf90_put_var (ncid, bpe2_bin_id , bpe2_bin , start=start3, count=count3)
      status = nf90_put_var (ncid, zppe2_bin_id, zppe2_bin, start=start3, count=count3)
      status = nf90_put_var (ncid, zmpe2_bin_id, zmpe2_bin, start=start3, count=count3)
      status = nf90_put_var (ncid, hp_bin_id   , hp_bin   , start=start3, count=count3)
      status = nf90_put_var (ncid, hm_bin_id   , hm_bin   , start=start3, count=count3)

      status = nf90_sync (ncid)

      nout = nout + 1
    endif

    if (proc0) call put_time_stamp(timer_io_total)
  end subroutine loop_io


!-----------------------------------------------!
!> @author  YK
!! @date    28 Jun 2021
!! @brief   Append cross section via MPIIO
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
    use grid, only: nlx, nly, nlz
    use grid, only: ilx_st, ily_st, ilz_st, ilx_en, ily_en, ilz_en
    use time, only: tt
    use mp, only: proc0
    use mpiio, only: mpiio_write_var_2d
    use time_stamp, only: put_time_stamp, timer_io_total, timer_io_2D
    implicit none

    real(r8), intent(in) :: phi_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            phi_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            phi_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)

    real(r8), intent(in) :: psi_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            psi_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            psi_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)

    real(r8), intent(in) :: omg_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            omg_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            omg_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)

    real(r8), intent(in) :: jpa_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            jpa_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            jpa_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)

    real(r8), intent(in) :: ux_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            ux_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            ux_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)
    real(r8), intent(in) :: uy_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            uy_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            uy_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)
                                                                                  
    real(r8), intent(in) :: bx_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            bx_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            bx_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)
    real(r8), intent(in) :: by_r_z0(ilx_st:ilx_en, ily_st:ily_en), &
                            by_r_x0(ily_st:ily_en, ilz_st:ilz_en), &
                            by_r_y0(ilx_st:ilx_en, ilz_st:ilz_en)

    integer, dimension(2) :: sizes, subsizes, starts

    if (proc0) call put_time_stamp(timer_io_total)
    if (proc0) call put_time_stamp(timer_io_2D)

    !--------------------------------------------------!
    !                    z = 0 cut
    !--------------------------------------------------!
    if(ilz_st == 1) then ! only the processe that has z = 0 write
      sizes(1) = nlx
      sizes(2) = nly
      subsizes(1) = ilx_en - ilx_st + 1
      subsizes(2) = ily_en - ily_st + 1
      starts(1) = ilx_st - 1
      starts(2) = ily_st - 1
    else
      sizes(1) = nlx
      sizes(2) = nly
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_phi_r_z0, disp_phi_r_z0, sizes, subsizes, starts, phi_r_z0)
    call mpiio_write_var_2d(fh_psi_r_z0, disp_psi_r_z0, sizes, subsizes, starts, psi_r_z0)
    call mpiio_write_var_2d(fh_omg_r_z0, disp_omg_r_z0, sizes, subsizes, starts, omg_r_z0)
    call mpiio_write_var_2d(fh_jpa_r_z0, disp_jpa_r_z0, sizes, subsizes, starts, jpa_r_z0)

    call mpiio_write_var_2d(fh_ux_r_z0 , disp_ux_r_z0 , sizes, subsizes, starts, ux_r_z0 )
    call mpiio_write_var_2d(fh_uy_r_z0 , disp_uy_r_z0 , sizes, subsizes, starts, uy_r_z0 )
    call mpiio_write_var_2d(fh_bx_r_z0 , disp_bx_r_z0 , sizes, subsizes, starts, bx_r_z0 )
    call mpiio_write_var_2d(fh_by_r_z0 , disp_by_r_z0 , sizes, subsizes, starts, by_r_z0 )

    !--------------------------------------------------!
    !                    x = 0 cut
    !--------------------------------------------------!
    if(ilx_st == 1) then ! only the processe that has x = 0 write
      sizes(1) = nly
      sizes(2) = nlz
      subsizes(1) = ily_en - ily_st + 1
      subsizes(2) = ilz_en - ilz_st + 1
      starts(1) = ily_st - 1
      starts(2) = ilz_st - 1
    else
      sizes(1) = nly
      sizes(2) = nlz
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_phi_r_x0, disp_phi_r_x0, sizes, subsizes, starts, phi_r_x0)
    call mpiio_write_var_2d(fh_psi_r_x0, disp_psi_r_x0, sizes, subsizes, starts, psi_r_x0)
    call mpiio_write_var_2d(fh_omg_r_x0, disp_omg_r_x0, sizes, subsizes, starts, omg_r_x0)
    call mpiio_write_var_2d(fh_jpa_r_x0, disp_jpa_r_x0, sizes, subsizes, starts, jpa_r_x0)

    call mpiio_write_var_2d(fh_ux_r_x0 , disp_ux_r_x0 , sizes, subsizes, starts, ux_r_x0 )
    call mpiio_write_var_2d(fh_uy_r_x0 , disp_uy_r_x0 , sizes, subsizes, starts, uy_r_x0 )
    call mpiio_write_var_2d(fh_bx_r_x0 , disp_bx_r_x0 , sizes, subsizes, starts, bx_r_x0 )
    call mpiio_write_var_2d(fh_by_r_x0 , disp_by_r_x0 , sizes, subsizes, starts, by_r_x0 )

    !--------------------------------------------------!
    !                    y = 0 cut
    !--------------------------------------------------!
    if(ily_st == 1) then ! only the processe that has y = 0 write
      sizes(1) = nlx
      sizes(2) = nlz
      subsizes(1) = ilx_en - ilx_st + 1
      subsizes(2) = ilz_en - ilz_st + 1
      starts(1) = ilx_st - 1
      starts(2) = ilz_st - 1
    else
      sizes(1) = nlx
      sizes(2) = nlz
      subsizes(1) = 1
      subsizes(2) = 1
      starts(1) = 0
      starts(2) = 0
    endif

    call mpiio_write_var_2d(fh_phi_r_y0, disp_phi_r_y0, sizes, subsizes, starts, phi_r_y0)
    call mpiio_write_var_2d(fh_psi_r_y0, disp_psi_r_y0, sizes, subsizes, starts, psi_r_y0)
    call mpiio_write_var_2d(fh_omg_r_y0, disp_omg_r_y0, sizes, subsizes, starts, omg_r_y0)
    call mpiio_write_var_2d(fh_jpa_r_y0, disp_jpa_r_y0, sizes, subsizes, starts, jpa_r_y0)

    call mpiio_write_var_2d(fh_ux_r_y0 , disp_ux_r_y0 , sizes, subsizes, starts, ux_r_y0 )
    call mpiio_write_var_2d(fh_uy_r_y0 , disp_uy_r_y0 , sizes, subsizes, starts, uy_r_y0 )
    call mpiio_write_var_2d(fh_bx_r_y0 , disp_bx_r_y0 , sizes, subsizes, starts, bx_r_y0 )
    call mpiio_write_var_2d(fh_by_r_y0 , disp_by_r_y0 , sizes, subsizes, starts, by_r_y0 )

    if(proc0) then
      write (unit=out2d_time_unit, fmt="(100es30.21)") tt
      flush (out2d_time_unit)
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
      write (unit=out3d_time_unit, fmt="(100es30.21)") tt
      flush (out3d_time_unit)
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

    !3D
    call MPI_FILE_CLOSE(fh_phi,ierr)
    call MPI_FILE_CLOSE(fh_psi,ierr)


    !2D
    call MPI_FILE_CLOSE(fh_phi_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_phi_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_phi_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_psi_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_psi_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_psi_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_omg_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_omg_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_omg_r_y0,ierr)

    call MPI_FILE_CLOSE(fh_jpa_r_z0,ierr)
    call MPI_FILE_CLOSE(fh_jpa_r_x0,ierr)
    call MPI_FILE_CLOSE(fh_jpa_r_y0,ierr)

    if(proc0) then
      call close_file (out2d_time_unit)
      call close_file (out3d_time_unit)
    endif

    if(proc0) then
      status = nf90_close (ncid)
    endif

  end subroutine finish_io

end module io

