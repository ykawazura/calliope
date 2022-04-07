program calliope
  use mp, only: init_mp, finish_mp, proc0
  use params, only: init_params, write_intvl, write_intvl_2D, write_intvl_3D
  use params, only: write_intvl_kpar, write_intvl_SF2, save_restart_intvl
  use terminate, only: init_terminate, monitor_terminate, check_terminate, terminated
  use grid, only: init_grid
  use time, only: init_time, tt, nstep, dt
  use fields, only: init_fields, finish_fields
  use dealias, only: init_filter
  use diagnostics, only: init_diagnostics, finish_diagnostics, loop_diagnostics, &
                         loop_diagnostics_2D, loop_diagnostics_kpar, loop_diagnostics_SF2
  use time_stamp, only: put_time_stamp, &
                        timer_total, timer_init, &
                        timer_diagnostics_total, timer_diagnostics_SF2, timer_diagnostics_kpar, &
                        timer_io_total, timer_io_2D, timer_io_3D, timer_save_restart, &
                        timer_advance, timer_nonlinear_terms, & 
                        timer_fft
  use model_specific, only :init_model_specific
  use io, only: loop_io_3D, save_restart
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
    call loop_diagnostics_2D
    call loop_io_3D
    call loop_diagnostics_kpar
    call save_restart
  endif

  do istep = 1, nstep
    call solve
    if(check_write_now(write_intvl       )) call loop_diagnostics      ! output for time history & spectra
    if(check_write_now(write_intvl_2D    )) call loop_diagnostics_2D   ! output for cross seciton of fields
    if(check_write_now(write_intvl_3D    )) call loop_io_3D            ! output for full 3D fields
    if(check_write_now(write_intvl_kpar  )) call loop_diagnostics_kpar ! output for kpar & delta b/b0 using k-filtering
    if(check_write_now(write_intvl_SF2   )) call loop_diagnostics_SF2  ! output for 2nd order structure function
    if(check_write_now(save_restart_intvl)) call save_restart
    if(check_write_now(write_intvl) .and. proc0) write(*, "('step = ', I10, '/', I10, ',  time = ', f15.8)") istep, nstep, tt
    if(monitor_terminate) then
      call check_terminate
      if(terminated) exit
    endif
  enddo

  call loop_diagnostics
  call loop_diagnostics_2D
  call loop_io_3D
  call loop_diagnostics_kpar
  call save_restart

  call finish_diagnostics
  call finish_fields

  if (proc0) call put_time_stamp(timer_total)

  if (proc0) then
    print '(/,'' Initialization'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &'' Advance steps'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    nonlinear terms'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    FFT'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &'' Diagnostics'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    SF2'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    kpar'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &'' IO'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    2D'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &''    3D'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/, &
          &'' Save restart'',T25,0pf9.3,'' min'',T40,2pf5.1,'' %'',/,/, &
          &'' total from timer is:'',T24,0pf10.3,'' min'',/, &
          &'' (init + adv + diag + io + restart:'',T40,2pf5.1,'' %)'',/)', &
          timer_init(1)/60.,timer_init(1)/timer_total(1), &
          timer_advance(1)/60.,timer_advance(1)/timer_total(1), &
            timer_nonlinear_terms(1)/60.,timer_nonlinear_terms(1)/timer_total(1), &
            timer_fft(1)/60.,timer_fft(1)/timer_total(1), &
          timer_diagnostics_total(1)/60.,timer_diagnostics_total(1)/timer_total(1), &
            timer_diagnostics_SF2(1)/60.,timer_diagnostics_SF2(1)/timer_total(1), &
            timer_diagnostics_kpar(1)/60.,timer_diagnostics_kpar(1)/timer_total(1), &
          timer_io_total(1)/60.,timer_io_total(1)/timer_total(1), &
            timer_io_2D(1)/60.,timer_io_2D(1)/timer_total(1), &
            timer_io_3D(1)/60.,timer_io_3D(1)/timer_total(1), &
          timer_save_restart(1)/60.,timer_save_restart(1)/timer_total(1), &
          timer_total(1)/60., &
          (timer_init(1) + timer_advance(1) + timer_diagnostics_total(1) &
           + timer_io_total(1) + timer_save_restart(1))/timer_total(1)

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


  function check_write_now(intvl)
    use p3dfft
    use time, only: tt, dt
    implicit none
    real(r8), intent(in) :: intvl
    logical  :: check_write_now
    integer, parameter :: i16 = selected_int_kind(16)

    check_write_now = floor(tt/intvl, i16) /= floor((tt+dt)/intvl, i16)

  end function check_write_now

end program calliope
