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
  logical :: pruned
  integer :: nlx_in , nly_in , nlz_in
  integer :: nkx_in , nky_in , nkz_in
  integer :: nlx_out, nly_out, nlz_out
  integer :: nkx_out, nky_out, nkz_out

  integer :: ierr, p_row
  logical :: rank0
  character(len=100) :: runname = 'expand', arg, inputfile, input_dir, output_dir
  complex(mytype), allocatable, dimension(:,:,:) :: u_in  
  complex(mytype), allocatable, dimension(:,:,:) :: u_out  
  TYPE(DECOMP_INFO), save :: sp_in, sp_med, sp_out  ! spectral space;
                                                    !  sp_in  -> (nlx_in , nlz_in , nly_in )
                                                    !  sp_med -> (nlx_in , nlz_out, nly_out)
                                                    !  sp_out -> (nlx_out, nlz_out, nly_out)
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
    integer  :: nlxc_in , nlyc_in , nlzc_in
    integer  :: nlxc_out, nlyc_out, nlzc_out

    call MPI_INIT(ierr)

    ! set runname
    call getarg(1,arg)
    if (len(trim(arg)) /= 0) runname = arg

    ! read inputfile = $(runname).in
    inputfile = trim(runname)//".in"
    call read_parameters(inputfile)

    if(pruned) then
      nlxc_in  = int(nlx_in *2.d0/3.d0) + 2 ! This is for dealiasing with 2/3 pruned FFT
      nlyc_in  = int(nly_in *2.d0/3.d0) + 2 ! This is for dealiasing with 2/3 pruned FFT
      nlzc_in  = int(nlz_in *2.d0/3.d0) + 2 ! This is for dealiasing with 2/3 pruned FFT
      nlxc_out = int(nlx_out*2.d0/3.d0) + 2 ! This is for dealiasing with 2/3 pruned FFT
      nlyc_out = int(nly_out*2.d0/3.d0) + 2 ! This is for dealiasing with 2/3 pruned FFT
      nlzc_out = int(nlz_out*2.d0/3.d0) + 2 ! This is for dealiasing with 2/3 pruned FFT
    else
      nlxc_in  = nlx_in 
      nlyc_in  = nly_in 
      nlzc_in  = nlz_in 
      nlxc_out = nlx_out
      nlyc_out = nly_out
      nlzc_out = nlz_out
    endif

    nkx_in  = nlxc_in
    nky_in  = nlyc_in/2 + 1
    nkz_in  = nlzc_in
    nkx_out = nlxc_out
    nky_out = nlyc_out/2 + 1
    nkz_out = nlzc_out


    call decomp_2d_init(nlx_in, nlz_in, nly_in, p_row, 1)
    call decomp_info_init(nkx_in , nkz_in , nky_in , sp_in )
    call decomp_info_init(nkx_in , nkz_out, nky_out, sp_med)
    call decomp_info_init(nkx_out, nkz_out, nky_out, sp_out)
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
    namelist /parameters/ pruned, &
                          nlx_in , nly_in , nlz_in , &
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
    complex(mytype), allocatable, dimension(:,:,:) :: fld_med_z_pencil, fld_med_y_pencil, fld_med_x_pencil

    complex(mytype), allocatable, dimension(:,:,:) :: fld_y_pencil, fld_x_pencil
    integer, allocatable :: ikx_in(:), ikz_in(:), ikx_out(:), ikz_out(:)
    integer :: i, j, k, ii, kk

    call alloc_z(fld_med_z_pencil, sp_med, opt_global=.true.); fld_med_z_pencil = 0._mytype
    call alloc_y(fld_med_y_pencil, sp_med, opt_global=.true.); fld_med_y_pencil = 0._mytype
    call alloc_x(fld_med_x_pencil, sp_med, opt_global=.true.); fld_med_x_pencil = 0._mytype

    call alloc_y(fld_y_pencil, sp_out, opt_global=.true.); fld_y_pencil = 0._mytype
    call alloc_x(fld_x_pencil, sp_out, opt_global=.true.); fld_x_pencil = 0._mytype

    allocate(ikx_in (nkx_in ))
    allocate(ikz_in (nkz_in ))
    allocate(ikx_out(nkx_out))
    allocate(ikz_out(nkz_out))

    do i = 1, nkx_in
      if (i <= nkx_in/2 + 1) then
        ikx_in(i) = i - 1
      else
        ikx_in(i) = i - nkx_in - 1
      endif

    enddo
    do k = 1, nkz_in
      if (k <= nkz_in/2 + 1) then
        ikz_in(k) = k - 1
      else
        ikz_in(k) = k - nkz_in - 1
      endif
    enddo

    do i = 1, nkx_out
      if (i <= nkx_out/2 + 1) then
        ikx_out(i) = i - 1
      else
        ikx_out(i) = i - nkx_out - 1
      endif

    enddo
    do k = 1, nkz_out
      if (k <= nkz_out/2 + 1) then
        ikz_out(k) = k - 1
      else
        ikz_out(k) = k - nkz_out - 1
      endif
    enddo

    ! (nkx_in, *nkz_in, *nky_in) -> (nkx_in, *nkz_out, *nky_out) : * means the slot is localized
    do k = 1, nkz_in
      kk = minloc(abs(ikz_out - ikz_in(k)), 1)
      do j = 1, nky_in
        fld_med_z_pencil(:, kk, j) = fld_in(:, k, j)
      end do
    enddo

    ! (nkx_in, *nkz_out, *nky_out) -> (*nkx_in, nkz_out, *nky_out) : * means the slot is localized
    call transpose_z_to_y(fld_med_z_pencil, fld_med_y_pencil, sp_med)
    call transpose_y_to_x(fld_med_y_pencil, fld_med_x_pencil, sp_med)

    ! (*nkx_in, nkz_out, *nky_out) -> (*nkx_out, nkz_out, *nky_out) : * means the slot is localized
    do i = 1, nkx_in
      ii = minloc(abs(ikx_out - ikx_in(i)), 1)
      fld_x_pencil(ii, :, :) = fld_med_x_pencil(i, :, :)
    enddo

    ! (*nkx_out, nkz_out, *nky_out) -> (nkx_out, *nkz_out, *nky_out) : * means the slot is localized
    call transpose_x_to_y(fld_x_pencil, fld_y_pencil, sp_out)
    call transpose_y_to_z(fld_y_pencil, fld_out     , sp_out)

    deallocate(fld_med_z_pencil)
    deallocate(fld_med_y_pencil)
    deallocate(fld_med_x_pencil)
    deallocate(fld_y_pencil)
    deallocate(fld_x_pencil)
    deallocate(ikx_in )
    deallocate(ikz_in )
    deallocate(ikx_out)
    deallocate(ikz_out)


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

