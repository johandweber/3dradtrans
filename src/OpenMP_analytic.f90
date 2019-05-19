!==============================================================================
! **  MODULES defining global variables and arrays for 3d_rad_trans_metals   **
!==============================================================================
!
! content of  absorption.f90:
! ==========================
! module M_absorption
!
! 
! content of data_input.f90:
! ===========================
! module M_data_input
!
!
! content of data_types.f90:
! ============================
! module M_data_types
!
!
! content of definitions.f90:
! ===========================
! module M_definitions
!
!
! content of geometry.f90:
! ========================
! module M_raysave
! module M_grid_memory
! module M_geometry
! module M_cosmology
!
!
! content of hydrogen.f90:
! ========================
! module M_hydrogen
!
!
! content of min_healpix.f90:
! ===========================
! module M_min_healpix
!
!
! content of module natural_constants.f90:
! =========================================
! module M_natural_constants.f90
!
!
! content of output.f90:
! ======================
! module M_output
!
! contents of write_ifrit.f90:
! ============================
! M_output_ifrit

!===============================================================================
PROGRAM radtrans_3d
!
  use M_data_types
  use M_definitions          , ONLY:  load_parameters,&
                                      read_from_file, resumed, diffuse,&
                                      periodic, trace_expansion,&
                                      ifrit_filename,&
                                      cpu_start, cpu_finish,&
                                      wallclock_start, wallclock_finish,&
                                      last_step, o_external,&
                                      progwccount, progwcrate, progwccountlast,&
                                      progcputime, progcputimelast,&
                                      t_start, t_end, t_all,t_step,&
                                      compute_temperature, include_metals,&
                                      num_sources, redshift,&
                                      command_line_options,&
                                      spectra, o_thermal_pressure,&
                                      add_heating, extra_heating,&
                                      verbose_stdout, use_mask,&
                                      binary_output, include_He
  
  use M_data_input            , ONLY: data_transform
  use M_raysave               , ONLY: allocate_rays
  use M_grid_memory           , ONLY: allocate_grid, mask, NHeI, nHeII, nHeIII, nHe_complete
  use M_geometry              , ONLY: init_geometry
  use M_absorption            , ONLY: SET_RAY_DIRECTIONS,&
                                      ALLOCATE_PHOTON_ARRAYS, PHOTON_ZERO,&
                                      GET_J, GET_J_PERIODIC, DIFFUSE_FIELD
  use M_hydrogen              , ONLY: INITIALIZE_ABUNDANCES,&
                                      INITIALIZE_HYDROGEN_FROM_FILE,&
                                      INITIALIZE_HYDROGEN, INITIALIZE_MASK,&
                                      RESUME_COMPUTATION,&
                                      RESUME_COMPUTATION_FROM_BINARY,&
                                      UPDATE_CHI, UPDATE_ETA, CALC_HYDROGEN,&
                                      CALC_HELIUM2, CALC_METALS, &
                                      UPDATE_ELECTRON_DENSITY,&
                                      INITIALIZE_HEATING_POLY,&
                                      APPLY_ADD_HEATING,&
                                      APPLY_EXTRA_HEATING_COOLING,&
                                      UPDATE_THERMAL_PRESSURE,&
                                      UPDATE_TEMPERATURE, EXPAND_SPACE,&
                                      rep
  use M_output                 ,ONLY: WRITE_IFRIT, CREATE_LOGFILE,&
                                      LOG_ITERATION_STEP,&
                                      CALL_EXTERNAL_HC_PROGRAM,&
                                      WRITE_IFRIT_BIN
!  
  implicit none
!  
  integer(i4b) :: i
!==============================================================================  
!
  call cpu_time(cpu_start)
  call system_clock(wallclock_start,progwcrate)
  if (command_argument_count() .eq. 2) then
     call get_command_argument(1,command_line_options(1))
     call get_command_argument(2, command_line_options(2))
     if (command_line_options(1) .eq. '-resume') then
        resumed=.true.
     else
        stop 'Invalid ccommand line options'
     end if
  end if
!
  call LOAD_PARAMETERS
  call ALLOCATE_PHOTON_ARRAYS
!  
  do i=1, num_sources
     call DATA_TRANSFORM(spectra(i))
     call PHOTON_ZERO(i)
  end do
!
  call ALLOCATE_RAYS
  call ALLOCATE_GRID
  call INIT_GEOMETRY
  call SET_RAY_DIRECTIONS

!
  call INITIALIZE_ABUNDANCES
!  
  if (read_from_file) then
     call INITIALIZE_HYDROGEN_FROM_FILE(ifrit_filename)
  else 
     call INITIALIZE_HYDROGEN
  end if
!
  if (resumed) then
     if (binary_output) then
        call RESUME_COMPUTATION_FROM_BINARY
     else
        call RESUME_COMPUTATION
     end if
  end if
!
  if (use_mask) then
     call INITIALIZE_MASK
  end if
!
  if (add_heating) then
     call INITIALIZE_HEATING_POLY
  end if
! 
!for output of initial hydrogen density
  last_step = .true.
  call WRITE_IFRIT
  last_step = .false.
  call system_clock(progwccountlast)
  open(unit=1, file='time.txt')
  open(unit=33, file='coordinates.txt')
  t_step=t_start 
!  
  do while(t_all/31557600._dp .lt. t_end)
!
     write(*,*)
     write(*,*)
     write(*,*) "TIMESTEP #",rep+1
     write(*,*) "elapsed time in yr. :", t_all/31557600._dp
     rep=rep+1
!
     call UPDATE_CHI
     write(*,*) " update_chi"
!
     if (diffuse) then
        call UPDATE_ETA
        write(*,*) "update_eta"
     end if
!
     print*, "before get_j"
     if (periodic) then
        call GET_J_PERIODIC
     else
        call GET_J
     end if
     write(*,*) "get_j"
!
     if (diffuse) then 
        call DIFFUSE_FIELD
        write(*,*) "diffuse_ field"
     end if
! 
     call CALC_HYDROGEN
     write(*,*) "calc_hydrogen"
!
     if (include_He) then
        call CALC_HELIUM2
        write(*,*) "calc_helium2"
     end if
!
     if (include_metals) then
          call CALC_METALS
          write(*,*) "calc metals"
       end if
!
     call UPDATE_ELECTRON_DENSITY
     write(*,*) "update electron density"
!
     if (add_heating) then
        call APPLY_ADD_HEATING
        write(*,*) "apply_add_heating"
     end if

     if (extra_heating) then
        call CALL_EXTERNAL_HC_PROGRAM
        call APPLY_EXTRA_HEATING_COOLING
     end if
!
    if (compute_temperature) then
        call UPDATE_TEMPERATURE
        write(*,*) "update_temperature"
     end if

     if (o_thermal_pressure) then
        call UPDATE_THERMAL_PRESSURE
        write(*,*) "update_thermal_pressure"
     end if

 !
     if (binary_output) then
        call WRITE_IFRIT_BIN
     else
        call WRITE_IFRIT
     end if
     write(*,*) "write_ifrit"
!
     if (trace_expansion) then
        call EXPAND_SPACE
        write(*,*) "expand_space"
     end if

!
     if (verbose_stdout) then
        call LOG_ITERATION_STEP
        write(*,*) "log_iteration_step"
     end if

     call cpu_time(progcputime)
     call system_clock(progwccount, progwcrate)

     if (trace_expansion) then
        write(1,'(F15.3, F15.3, F10.4, I6, F15.3 ,F15.3)') &
             t_step/31556925.9747_dp, t_all/31556925.9747_dp,&
             redshift,rep , (progwccount-progwccountlast)/(1._dp*progwcrate),&
             progcputime-progcputimelast
        flush(1)
     else
        write(1,'(F15.3, F15.3, I6, F15.3 ,F15.3)') &
             t_step/31556925.9747_dp, t_all/31556925.9747_dp, rep,&
             (progwccount-progwccountlast)/(1._dp*progwcrate),&
             progcputime-progcputimelast
        flush(1)
     end if
!
     progcputimelast=progcputime
     progwccountlast=progwccount
!
  end do
!  
  last_step = .true.
  if (binary_output) then
     call WRITE_IFRIT_BIN
  else
     call WRITE_IFRIT
  end if
  call cpu_time(cpu_finish)
  call system_clock(wallclock_finish,progwcrate)

  call CREATE_LOGFILE
!  
  write(*,*) "Wallclock time in s :", 1._dp*(wallclock_finish-wallclock_start)/&
                                      progwcrate
  write(*,*) "CPU time in s       :", cpu_finish-cpu_start
  close(1)
  close(33)
!------------------------------------------------------------------------------ 
end program radtrans_3d
