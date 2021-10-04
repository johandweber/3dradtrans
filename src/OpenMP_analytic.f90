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
PROGRAM RADTRANS_3D
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
  call load_parameters
  call allocate_photon_arrays
!  
  do i=1, num_sources
     call data_transform(spectra(i))
     call photon_zero(i)
  end do
!
  call allocate_rays
  call allocate_grid
  call init_geometry
  call set_ray_directions

!
  call initialize_abundances
!  
  if (read_from_file) then
     call initialize_hydrogen_from_file(ifrit_filename)
  else 
     call initialize_hydrogen
  end if
!
  if (resumed) then
     if (binary_output) then
        call resume_computation_from_binary
     else
        call resume_computation
     end if
  end if
!
  if (use_mask) then
     call initialize_mask
  end if
!
  if (add_heating) then
     call initialize_heating_poly
  end if
! 
!for output of initial hydrogen density
  last_step = .true.
  call write_ifrit
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
     call update_chi
     write(*,*) " update_chi"
!
     if (diffuse) then
        call update_eta
        write(*,*) "update_eta"
     end if
!
     print*, "before get_j"
     if (periodic) then
        call get_j_periodic
     else
        call get_j
     end if
     write(*,*) "get_j"
!
     if (diffuse) then 
        call diffuse_field
        write(*,*) "diffuse_ field"
     end if
! 
     call calc_hydrogen
     write(*,*) "calc_hydrogen"
!
     if (include_He) then
        call calc_helium2
        write(*,*) "calc_helium2"
     end if
!
     if (include_metals) then
          call calc_metals
          write(*,*) "calc metals"
       end if
!
     call update_electron_density
     write(*,*) "update electron density"
!
     if (add_heating) then
        call apply_add_heating
        write(*,*) "apply_add_heating"
     end if

     if (extra_heating) then
        call call_external_hc_program
        call apply_extra_heating_cooling
     end if
!
    if (compute_temperature) then
        call update_temperature
        write(*,*) "update_temperature"
     end if
!
     if (o_thermal_pressure) then
        call update_thermal_pressure
        write(*,*) "update_thermal_pressure"
     end if
 !
     if (binary_output) then
        call write_ifrit_bin
     else
        call write_ifrit
     end if
     write(*,*) "write_ifrit"
!
     if (trace_expansion) then
        call expand_space
        write(*,*) "expand_space"
     end if

!
     if (verbose_stdout) then
        call log_iteration_step
        write(*,*) "log_iteration_step"
     end if
!
     call cpu_time(progcputime)
     call system_clock(progwccount, progwcrate)
!
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
     call write_ifrit_bin
  else
     call write_ifrit
  end if
  call cpu_time(cpu_finish)
  call system_clock(wallclock_finish,progwcrate)
!
  call create_logfile
!  
  write(*,*) "Wallclock time in s :", 1._dp*(wallclock_finish-wallclock_start)/&
                                      progwcrate
  write(*,*) "CPU time in s       :", cpu_finish-cpu_start
  close(1)
  close(33)
!------------------------------------------------------------------------------ 
END PROGRAM RADTRANS_3D
