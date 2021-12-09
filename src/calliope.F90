program calliope
  use mp, only: init_mp, finish_mp, proc0
  use params, only: init_params, nwrite, nwrite_fld_section, nwrite_fld_3D, nwrite_SF2, nsave_restart
  use terminate, only: init_terminate, monitor_terminate, check_terminate, terminated
  use grid, only: init_grid
  use time, only: init_time, tt, nstep, dt
  use fields, only: init_fields, finish_fields
  use dealias, only: init_filter
  use diagnostics, only: init_diagnostics, finish_diagnostics, loop_diagnostics, &
                         loop_diagnostics_fields_secion, loop_diagnostics_SF2
  use time_stamp, only: put_time_stamp, &
                        timer_total, timer_init, timer_diagnostics, timer_diagnostics_SF2, &
                        timer_advance, timer_nonlinear_terms, & 
                        timer_fft
  use model_specific, only :init_model_specific
  use io, only: loop_io_fields_3D, save_restart
  use advance, only :solve
  implicit none

  integer :: istep

  call init_params
  call init_mp
  call banner
  call init_terminate
  call init_grid
  call init_time

  if (proc0) then
    write(*, '("Solving ")', advance='no')
    write(*, '(a)') trim(_MODEL_)
  endif

  if (proc0) call put_time_stamp(timer_total)

  if (proc0) call put_time_stamp(timer_init)
  call init_fields
  call init_model_specific
  call init_diagnostics
  call init_filter
  if (proc0) call put_time_stamp(timer_init)

  if(tt == 0.d0) then ! skip when restarting
    call loop_diagnostics
    call loop_diagnostics_fields_secion
    call loop_io_fields_3D
    call save_restart
  endif

  do istep = 1, nstep
    call solve
    if(mod(istep, nwrite            ) == 0) call loop_diagnostics               ! output for time history & spectra
    if(mod(istep, nwrite_fld_section) == 0) call loop_diagnostics_fields_secion ! output for cross seciton of fields
    if(mod(istep, nwrite_fld_3D     ) == 0) call loop_io_fields_3D              ! output for full 3D fields
    if(mod(istep, nwrite_SF2        ) == 0) call loop_diagnostics_SF2           ! output for 2nd order structure function
    if(mod(istep, nsave_restart     ) == 0) call save_restart
    if(mod(istep, nwrite) == 0 .and. proc0) write(*, "('step = ', I10, '/', I10, ',  time = ', f15.8)") istep, nstep, tt
    if(monitor_terminate) then
      call check_terminate
      if(terminated) exit
    endif
  enddo

  call loop_diagnostics
  call loop_diagnostics_fields_secion
  call loop_io_fields_3D
  call save_restart

  call finish_diagnostics
  call finish_fields

  if (proc0) call put_time_stamp(timer_total)

  if (proc0) then
    print '(/,'' Initialization'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &'' Advance steps'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    nonlinear terms'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    FFT'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &'' diagnostics'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %)'',/, &
          &''    SF2'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %)'',/, &
          &'' total from timer is:'',T24,0pf10.3,'' min'',/)', &
          timer_init(1)/60.,timer_init(1)/timer_total(1), &
          timer_advance(1)/60.,timer_advance(1)/timer_total(1), &
          timer_nonlinear_terms(1)/60.,timer_nonlinear_terms(1)/timer_total(1), &
          timer_fft(1)/60.,timer_fft(1)/timer_total(1), &
          timer_diagnostics(1)/60.,timer_diagnostics(1)/timer_total(1), &
          timer_diagnostics_SF2(1)/60.,timer_diagnostics_SF2(1)/timer_total(1), &
          timer_total(1)/60.

    print *, '# of steps advanced', nstep
    print '(12X, "Final dt", es15.3)', dt
  endif

  call finish_mp
contains

  subroutine banner
    implicit none
    integer,dimension(8) :: datetime

    if(proc0) then
      call date_and_time(VALUES=datetime)
print *
write(*, "('  -------------------------------------------------------------------- ')", advance='no'); print *
write(*, "(' |       _____          _      _      _____ ____  _____  ______       |')", advance='no'); print *
write(*, "(' |      / ____|   /\   | |    | |    |_   _/ __ \|  __ \|  ____|      |')", advance='no'); print *
write(*, "(' |     | |       /  \  | |    | |      | || |  | | |__) | |__         |')", advance='no'); print *
write(*, "(' |     | |      / /\ \ | |    | |      | || |  | |  ___/|  __|        |')", advance='no'); print *
write(*, "(' |     | |____ / ____ \| |____| |____ _| || |__| | |    | |____       |')", advance='no'); print *
write(*, "(' |      \_____/_/    \_\______|______|_____\____/|_|    |______|      |')", advance='no'); print *
write(*, "(' |                                                                    |')", advance='no'); print *
write(*, "(' |', 36X, ' executed on ', i4, '-', i2, '-', i2, ' ', i2, ':', i2, 2X, ' |')", advance='no') &
                      datetime(1), datetime(2), datetime(3), datetime(5), datetime(6); print *
write(*, "('  -------------------------------------------------------------------- ')", advance='no'); print *
print *
    endif
  end subroutine banner

end program calliope
