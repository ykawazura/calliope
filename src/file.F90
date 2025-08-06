!-----------------------------------------------!
!> @author  YK
!! @date    15 Feb 2021
!! @brief   File utils. This module is mostly
!!          based on 'file_utils.fpp' of AstroGK
!-----------------------------------------------!
module file

  public open_input_file, open_output_file, open_output_file_binary
  public close_file
  public get_indexed_namelist_unit

contains

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Open a file for reading and return 
!!          its unit number in UNIT
!-----------------------------------------------!
  subroutine open_input_file (unit, fn)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: fn
    call get_unused_unit (unit)
    open (unit=unit, file=trim(fn), action="read")
  end subroutine open_input_file

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Open a file for writing and return 
!!          its unit number in UNIT
!-----------------------------------------------!
  subroutine open_output_file (unit, fn)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: fn
    call get_unused_unit (unit)
    open (unit=unit, file=trim(fn), status="replace", action="write")
  end subroutine open_output_file

!-----------------------------------------------!
!> @author  YK
!! @date    28 Jul 2022
!! @brief   Open a binary file for writing and return 
!!          its unit number in UNIT
!-----------------------------------------------!
  subroutine open_output_file_binary (unit, fn)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: fn
    call get_unused_unit (unit)
    open (unit=unit, file=trim(fn), form="unformatted", status="replace", access='stream')
  end subroutine open_output_file_binary

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Return a unit number not associated
!!          with any file
!-----------------------------------------------!
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

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Close the file associated with UNIT 
!!          from open_out(or in)put_file
!-----------------------------------------------!
  subroutine close_file (unit)
    implicit none
    integer, intent (in) :: unit
    close (unit=unit)
  end subroutine close_file

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Copy namelist, NML // '_' // INDEX, 
!!          from the input file FILENAME 
!!          to namelist, NML, in a temporary 
!!          file, UNIT 
!-----------------------------------------------!
  subroutine get_indexed_namelist_unit (unit, filename, nml, index_in)
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: nml
    character(len=100), intent(in) :: filename
    integer, intent (in) :: index_in
    character(500) :: line
    integer :: input_unit_no, iunit, iostat

    call get_unused_unit (input_unit_no)
    open(unit=input_unit_no,file=filename,status='old')

    call get_unused_unit (unit)
    open (unit=unit, file=trim(filename)//".tmp")

    write (line, *) index_in
    line = nml//"_"//trim(adjustl(line))
    iunit = input_unit(trim(line), input_unit_no)

    read (unit=iunit, fmt="(a)") line
    write (unit=unit, fmt="('&',a)") nml

    do
      read (unit=iunit, fmt="(a)", iostat=iostat) line
      if (iostat /= 0 .or. trim(adjustl(line)) == "/") exit
      write (unit=unit, fmt="(a)") trim(line)
    end do
    write (unit=unit, fmt="('/')")
    rewind (unit=unit)

    close(input_unit_no)
  end subroutine get_indexed_namelist_unit

!-----------------------------------------------!
!> @author  YK
!! @date    20 Sep 2019
!! @brief   Rewind the input file to start of
!!          namelist NML, and return its unit 
!!          number
!-----------------------------------------------!
  function input_unit (nml, input_unit_no)
    implicit none
    character(*), intent (in) :: nml
    integer, intent(in) :: input_unit_no
    integer :: input_unit, iostat
    character(500) :: line
    intrinsic adjustl, trim

    input_unit = input_unit_no
    if (input_unit_no > 0) then
      rewind (unit=input_unit_no)
      do
        read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
        if (iostat /= 0) then
          rewind (unit=input_unit_no)
          exit
        end if
        if (trim(adjustl(line)) == "&"//nml) then
          backspace (unit=input_unit_no)
          return
        end if
      end do
    end if
  end function input_unit

end module file


