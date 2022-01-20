program expand
  use MPI
  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_io
  implicit none

  ! for the list of field file names 
  integer :: nfld
  character(len=100), dimension(100) :: fldname

  ! size of before and after
  integer :: nlx_in , nly_in , nlz_in
  integer :: nkx_in , nky_in , nkz_in
  integer :: nlx_out, nly_out, nlz_out
  integer :: nkx_out, nky_out, nkz_out

  integer :: ierr, p_row, irank
  logical :: rank0
  character(len=100) :: runname = 'expand', arg, inputfile, input_dir, output_dir
  complex(mytype), allocatable, dimension(:,:,:) :: u_in  
  complex(mytype), allocatable, dimension(:,:,:) :: u_out  
  TYPE(DECOMP_INFO), save :: sp_in, sp_out  ! spectral space
  integer :: i

  call init

  if(rank0) then
    print * 
    print '("expanding: nlx : ", I5, " -> ", I5)', nlx_in, nlx_out 
    print '("           nly : ", I5, " -> ", I5)', nly_in, nly_out 
    print '("           nlz : ", I5, " -> ", I5)', nlz_in, nlz_out 
    print *, '    input_dir :  '//trim( input_dir)
    print *, '   output_dir :  '//trim(output_dir)
    print * 
  endif

  ! loop each field and expand it
  do i = 1, nfld
    call read_and_expand(trim(fldname(i)))
  end do

  if(rank0)  call system('cp '//trim(input_dir)//'/time.dat '//trim(output_dir)//'/time.dat')

  call finish

  if(rank0) then
    print '("Done!")'
    print *
  endif

contains

  subroutine init
    implicit none
    call MPI_INIT(ierr)

    ! set runname
    call getarg(1,arg)
    if (len(trim(arg)) /= 0) runname = arg

    ! read inputfile = $(runname).in
    inputfile = trim(runname)//".in"
    call read_parameters(inputfile)

    nkx_in  = nlx_in/2 + 1
    nky_in  = nly_in
    nkz_in  = nlz_in
    nkx_out = nlx_out/2 + 1
    nky_out = nly_out
    nkz_out = nlz_out

    call decomp_2d_init(nlx_in, nly_in, nlz_in, p_row, 1)
    call decomp_info_init(nkx_in , nky_in , nkz_in , sp_in )
    call decomp_info_init(nkx_out, nky_out, nkz_out, sp_out)
    rank0 = nrank == 0

    call alloc_z(u_in , sp_in , opt_global=.true.); u_in  = 0._mytype
    call alloc_z(u_out, sp_out, opt_global=.true.); u_out = 0._mytype

    ! get the list of filed file names and its size
    if(rank0) call get_restart_file_list
    call broadcast_integer(nfld)
    do i = 1, nfld
      call broadcast_character(fldname(i))
    end do
  end subroutine init


  subroutine get_restart_file_list
    implicit none
    integer :: i
    character(len=100) line

    call system('cd '//trim( input_dir)//'; ls -1 *.dat > ../tmp')

    open(10, file='tmp', status='old')

    nfld = 0
    do i = 1, 9999
      read(10, "(a)", end=999) line

      if(trim(line) /= 'time.dat') then
        nfld = nfld + 1
        fldname(nfld) = trim(line)
      endif
    end do
999 close(10)

    call system('rm tmp')

  end subroutine get_restart_file_list


  subroutine read_parameters(filename)
    implicit none
    
    character(len=100), intent(in) :: filename
    integer  :: unit, ierr

    namelist /mpi_parameters/ p_row
    namelist /parameters/ nlx_in , nly_in , nlz_in , &
                          nlx_out, nly_out, nlz_out, &
                          input_dir, output_dir
    p_row = 0

    call get_unused_unit (unit)
    open(unit=unit,file=filename,status='old')

    read(unit,nml=mpi_parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading mpi_parameters failed"
    read(unit,nml=parameters,iostat=ierr)
        if (ierr/=0) write(*,*) "Reading box_parameters failed"
    close(unit)
  end subroutine read_parameters


  subroutine open_input_file (unit, fn)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: fn
    call get_unused_unit (unit)
    open (unit=unit, file=trim(fn), action="read")
  end subroutine open_input_file


  subroutine open_output_file (unit, fn)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: fn
    call get_unused_unit (unit)
    open (unit=unit, file=trim(fn), status="replace", action="write")
  end subroutine open_output_file


  subroutine get_unused_unit (unit)
    implicit none
    integer, intent (out) :: unit
    logical :: od
    unit = 50
    do
       inquire (unit=unit, opened=od)
       if (.not.od) return
       unit = unit + 1
    end do
  end subroutine get_unused_unit


  subroutine close_file (unit)
    implicit none
    integer, intent (in) :: unit
    close (unit=unit)
  end subroutine close_file


  subroutine read_and_expand(fname)
    implicit none
    character(*), intent(in) :: fname

    ! read file !
    call decomp_2d_read_one(3, u_in, trim(input_dir)//trim(fname), sp_in)
    ! expand field !
    call expand_field(u_in, u_out)
    ! write field !
    call decomp_2d_write_one(3, u_out, trim(output_dir)//trim(fname), sp_out)

  end subroutine read_and_expand


  subroutine expand_field(fld_in, fld_out)
    implicit none
    complex(mytype), dimension (sp_in %zst(1):sp_in %zen(1), &
                                sp_in %zst(2):sp_in %zen(2), &
                                sp_in %zst(3):sp_in %zen(3)), intent(inout) :: fld_in
    complex(mytype), dimension (sp_out%zst(1):sp_out%zen(1), &
                                sp_out%zst(2):sp_out%zen(2), &
                                sp_out%zst(3):sp_out%zen(3)), intent(inout) :: fld_out
    complex(mytype), allocatable, dimension(:,:,:) :: fld_y_pencil, fld_x_pencil, tmp_x_pencil 
    integer :: i, j, k, jj, kk
    integer :: idx, xsize_in, xsize_out

    call alloc_y(fld_y_pencil, sp_out, opt_global=.true.); fld_y_pencil = 0._mytype
    call alloc_x(fld_x_pencil, sp_out, opt_global=.true.); fld_x_pencil = 0._mytype
    call alloc_x(tmp_x_pencil, sp_out, opt_global=.true.); tmp_x_pencil = 0._mytype

    do k = sp_in%zst(3), sp_in%zen(3)
      if (k > nkz_in/2 + 1) then
        kk = nkz_out/2 + k 
      else
        kk = k
      endif
      do j = sp_in%zst(2), sp_in%zen(2)
        if (j > nky_in/2 + 1) then
          jj = nky_out/2 + j 
        else
          jj = j
        endif
        idx = sp_out%zst(1)
        do i = sp_in%zst(1), sp_in%zen(1)
          fld_out(idx, jj, kk) = fld_in(i, j, k)
          idx = idx + 1
        end do
      end do
    end do

    call transpose_z_to_y(fld_out, fld_y_pencil, sp_out)
    call transpose_y_to_x(fld_y_pencil, fld_x_pencil, sp_out)

    idx = 1
    xsize_in  = sp_in %zsz(1)
    xsize_out = sp_out%zsz(1)
    call min_allreduce_integer(xsize_in )
    call min_allreduce_integer(xsize_out)
    do irank = 0, p_row - 1
      do i = 1, xsize_in
        tmp_x_pencil(idx, :, :) = fld_x_pencil(irank*xsize_out + i, :, :)
        idx = idx + 1
      enddo
    enddo

    call transpose_x_to_y(tmp_x_pencil, fld_y_pencil, sp_out)
    call transpose_y_to_z(fld_y_pencil, fld_out     , sp_out)
  end subroutine expand_field


  subroutine finish
    implicit none
    call decomp_info_finalize(sp_in)
    call decomp_info_finalize(sp_out)
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
  end subroutine finish

  subroutine broadcast_character (char)
    implicit none
    character(*), intent (in out) :: char
    integer :: ierror
    call mpi_bcast (char, len(char), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_character

  subroutine broadcast_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: ierror
    call mpi_bcast (i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  end subroutine broadcast_integer

  subroutine min_allreduce_integer (i)
    implicit none
    integer, intent (in out) :: i
    integer :: i1, ierror
    i1 = i
    call mpi_allreduce &
         (i1, i, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
  end subroutine min_allreduce_integer

end program expand

