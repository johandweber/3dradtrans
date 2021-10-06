!-------------------------------------------------------------------------------
!                           ****MODULE M_definitions****
! Here a number of values which stay constant during the run, but might be 
! changed in different runs are stored.
!-------------------------------------------------------------------------------

module M_definitions
!--- routines:
!    LOAD_PARAMETERS   
!
use M_data_types
implicit none

!frequency points
integer,parameter  :: points=100                        ! frequency points
                                                        ! number of considered elements

!Number of CPU-Threads
integer(i4b):: num_threads=1                            !Number of CPU-Threads

! maximal value for x (total range -x_max:+x_max)
integer(i4b):: x_max=15
!
! maximal value for y (total range -y_max:+y_max)
integer(i4b):: y_max=15
!
! maximal value for z (total range -z_max:+z_max)
integer(i4b):: z_max=15
!
! In the spherical case, the maximal number of radius
! points is given by the variable r_max.
! 
! To be able to reuse the code that computes the
! temperature and ionization structure, we than map r_max to
! z_max and set x_max and y_max to 0.
integer(i4b):: rp_max=100

! minimal value for x (total range -x_max:+x_max)
integer(i4b):: x_min=15
!
! minimal value for y (total range -y_max:+y_max)
integer(i4b):: y_min=15
!
! minimal value for z (total range -z_max:+z_max)
integer(i4b):: z_min=15
!
! number of radiation sources 
integer(i4b):: num_sources
!
! coordinates of radiation sources  
integer(i4b), dimension(:), allocatable:: source_x
integer(i4b), dimension(:), allocatable:: source_y
integer(i4b), dimension(:), allocatable:: source_z
!
! Radius of emitting source in solar radii
real(dp),     dimension(:), allocatable:: R_RS
!
! time, when source starts to emit
real(dp),     dimension(:), allocatable:: start_activity
!
! lifetime of emitting sources
real(dp),     dimension(:), allocatable:: lifetime
!
! filename from which source-spectrum is read
character (len=3000), dimension(:), allocatable:: spectra
!
! command-line options
character (len=3000), dimension(2) :: command_line_options
!
! filename trunk when execution is resumed
character (len=3000) :: resume_string
!
! is execution of program resumed ?
logical             :: resumed
!                       
! is set .true. when time_end is reached
logical             :: last_step = .false.
!
! if set .true., correction for speed of light is acitvated
logical             :: c_correct_on=.true.
!
! read file for density structure
logical             :: read_from_file=.false.
!
! start with ionized gas instead of neutral gas
logical             :: start_ionized=.false.
!
! running number for the source
integer(i4b)        :: source_counter
!
! variables for CPU_time measurement
real(dp)            :: cpu_start, cpu_finish
real(dp)            :: progcputime, progcputimelast=0
                       
!
integer(i8b)        :: wallclock_start, wallclock_finish
integer(i8b)        :: progwccount,progwcrate,progwccountlast=0
!
! initial fraction of neutral hydroge, number density of hydrogen
! (in homogeneous case)
real(dp)     :: nHI_fraction=0.9990_dp, nH_complete_scalar=10.0_dp
!
! "Security factor" as described in Abel et al. 2002
! The complete diagonal of the volume is considered to be the radius
real(dp)     :: f_a=0.1                         
                                                                                  
real(dp)     ::  l_cell          = 1.0_dp, &     ! lenght of a cell (in parsec)
                 t_start         = 1000000._dp,& ! time for first time_step
                 t_end           = 1E+5_dp, &    ! time, after which 
                                                 ! simulation ends   
                 t_all, &                        ! total passed time
                 worldage,&                      ! total paased time 
                                                 ! + begin of simulation
                 t_step, &                       ! time step size
                 t_step_max        = 4E+12_dp, & ! upper limit for time step size
                 factor            = 25._dp,&    ! factor to determine
                                                 ! the size of the timesteps
                 H_0               = 71._dp,&    ! Hubble constant
                 omega_m           = 0.269_dp,&  ! rho_matter/rho_crit
                 initial_redshift,&
                 redshift,&
                 redshift_old

real(dp)            :: lower_energy_limit=83363.7_dp/109678.8_dp


character (len=3000):: add_heating_filename=''               
character (len=3000):: ifrit_filename=''       
character (len=3000):: external_hc_program=''
character (len=3000):: external_chieta_program=''
character (len=3000):: mask_filename=''
character (len=3000):: source_filename=''
character (len=3000):: initial_file_spherical=''

integer(i4b)        :: rep_output=1               ! every "rep_output" repetition
                                                  ! an IFRIT-file is written

! abundance relative to hydrogen instead ofrelative to the solar metallicity
logical             :: abund_rel_to_H               = .false.
logical             :: add_heating                  = .false.                                           
logical             :: artificially_fix_temperature = .false.
logical             :: binary_output                = .false.
logical             :: charge_transfer              = .false.
logical             :: compute_temperature          = .false.
logical             :: diffuse                      = .false. 
logical             :: diffuse_random               = .true.
logical             :: equidistant                  = .true.
logical             :: extra_heating                = .false.
logical             :: extra_radiation_field        = .false.
logical             :: include_Ar                   = .true.
logical             :: include_C                    = .true.
logical             :: include_He                   = .true.
logical             :: include_metals               = .true.
logical             :: include_N                    = .true.
logical             :: include_Ne                   = .true.
logical             :: include_O                    = .true.
logical             :: include_S                    = .true.
logical             :: log_heating_cooling          = .false.
logical             :: mergespec                    = .false.
logical             :: no_source_filename           = .false.
logical             :: o_abundan                    = .false.
logical             :: o_external                   = .false.
logical             :: old_solar_abundances         = .false.
logical             :: optabu_3d                    = .false.
logical             :: o_rad_force                  = .false.
logical             :: o_thermal_pressure           = .false.
logical             :: o_write_escape_fraction      = .false.
logical             :: outward_only                 = .false.
logical             :: periodic                     = .false.
logical             :: photoheat_metals             = .false.
logical             :: spherical_case_A             = .false.
logical             :: spherical_case_B             = .false.
logical             :: planar_case_A                = .false.
logical             :: planar_case_B                = .false.
logical             :: trace_expansion              = .false.
logical             :: verbose_stdout               = .false.
logical             :: use_mask                     = .false.
!
real(dp)            :: Heabundance, Cabundance, Nabundance, Oabundance,&
                       Neabundance,Sabundance, Arabundance
integer(i4b)        :: Heions=3, Cions=4, Nions=4, Oions=4, Neions=4,&
                       Sions=4, Arions=4
!
real(dp)            :: Hetosolar = 1._dp,&
                       Litosolar = 1._dp,&
                       Betosolar = 1._dp,&
                       Btosolar  = 1._dp,&
                       Ctosolar  = 1._dp,&
                       Ntosolar  = 1._dp,&
                       Otosolar  = 1._dp,&
                       Ftosolar  = 1._dp,&
                       Netosolar = 1._dp,&
                       Natosolar = 1._dp,&
                       Mgtosolar = 1._dp,&
                       Altosolar = 1._dp,&
                       Sitosolar = 1._dp,&
                       Ptosolar  = 1._dp,&
                       Stosolar  = 1._dp,&
                       Cltosolar = 1._dp,&
                       Artosolar = 1._dp,&
                       Ktosolar  = 1._dp,&
                       Catosolar = 1._dp,&
                       Titosolar = 1._dp,&
                       Vtosolar  = 1._dp,&
                       Crtosolar = 1._dp,&
                       Mntosolar = 1._dp,&
                       Fetosolar = 1._dp,&
                       Cotosolar = 1._dp,&
                       Nitosolar = 1._dp,&
                       Cutosolar = 1._dp,&
                       Zntosolar = 1._dp,&
                       Gatosolar = 1._dp,&
                       Getosolar = 1._dp,&
                       Astosolar = 1._dp,&
                       Setosolar = 1._dp,&
                       Brtosolar = 1._dp,&
                       Krtosolar = 1._dp,&
                       Rbtosolar = 1._dp,&
                       Srtosolar = 1._dp,&
                       Ytosolar  = 1._dp,&
                       Zrtosolar = 1._dp,&
                       Nbtosolar = 1._dp,&
                       Motosolar = 1._dp,&
                       Tctosolar = 1._dp,&
                       Rutosolar = 1._dp,&
                       Rhtosolar = 1._dp,&
                       Pdtosolar = 1._dp,&
                       Agtosolar = 1._dp,&
                       Cdtosolar = 1._dp,&
                       Intosolar = 1._dp,&
                       Sntosolar = 1._dp,&
                       Sbtosolar = 1._dp,&
                       Tetosolar = 1._dp,&
                       Itosolar  = 1._dp,&
                       Xetosolar = 1._dp,&
                       Cstosolar = 1._dp,&
                       Batosolar = 1._dp,&
                       Latosolar = 1._dp,&
                       Cetosolar = 1._dp,&
                       Prtosolar = 1._dp,&
                       Ndtosolar = 1._dp,&
                       Pmtosolar = 1._dp,&
                       Satosolar = 1._dp,&
                       Eutosolar = 1._dp,&
                       Gdtosolar = 1._dp,&
                       Tbtosolar = 1._dp,&
                       Dytosolar = 1._dp,&
                       Hotosolar = 1._dp,&
                       Ertosolar = 1._dp,&
                       Tmtosolar = 1._dp,&
                       Ybtosolar = 1._dp,&
                       Lutosolar = 1._dp,&
                       Hftosolar = 1._dp,&
                       Tatosolar = 1._dp,&
                       Wtosolar  = 1._dp,&
                       Retosolar = 1._dp,&
                       Ostosolar = 1._dp,&
                       Irtosolar = 1._dp,&
                       Pttosolar = 1._dp,&
                       Autosolar = 1._dp,&
                       Hgtosolar = 1._dp,&
                       Tltosolar = 1._dp,&
                       Pbtosolar = 1._dp,&
                       Bitosolar = 1._dp,&
                       Potosolar = 1._dp,&
                       Attosolar = 1._dp,&
                       Rntosolar = 1._dp,&
                       Frtosolar = 1._dp,&
                       Ratosolar = 1._dp,&
                       Actosolar = 1._dp,&
                       Thtosolar = 1._dp,&
                       Patosolar = 1._dp,&
                       Utosolar  = 1._dp,&
                       Nptosolar = 1._dp,&
                       Putosolar = 1._dp,&
                       Amtosolar = 1._dp,&
                       Cmtosolar = 1._dp

! Note that in the case abund_rel_to_H = .true.
! the default abundances do _not_ correspond 
! to the solar avlues.
! Instead, they are set to 1.E-12 times the
! abundancde of hydrogen, which in most cases 
! should be negligibly small

real(dp)            :: HetoH  = 1.E-12_dp,&
                       LitoH  = 1.E-12_dp,&
                       BetoH  = 1.E-12_dp,&
                       BtoH   = 1.E-12_dp,&
                       CtoH   = 1.E-12_dp,&
                       NtoH   = 1.E-12_dp,&
                       OtoH   = 1.E-12_dp,&
                       FtoH   = 1.E-12_dp,&
                       NetoH  = 1.E-12_dp,&
                       NatoH  = 1.E-12_dp,&
                       MgtoH  = 1.E-12_dp,&
                       AltoH  = 1.E-12_dp,&
                       SitoH  = 1.E-12_dp,&
                       PtoH   = 1.E-12_dp,&
                       StoH   = 1.E-12_dp,&
                       CltoH  = 1.E-12_dp,&
                       ArtoH  = 1.E-12_dp,&
                       KtoH   = 1.E-12_dp,&
                       CatoH  = 1.E-12_dp,&
                       TitoH  = 1.E-12_dp,&
                       VtoH   = 1.E-12_dp,&
                       CrtoH  = 1.E-12_dp,&
                       MntoH  = 1.E-12_dp,&
                       FetoH  = 1.E-12_dp,&
                       CotoH  = 1.E-12_dp,&
                       NitoH  = 1.E-12_dp,&
                       CutoH  = 1.E-12_dp,&
                       ZntoH  = 1.E-12_dp,&
                       GatoH  = 1.E-12_dp,&
                       GetoH  = 1.E-12_dp,&
                       AstoH  = 1.E-12_dp,&
                       SetoH  = 1.E-12_dp,&
                       BrtoH  = 1.E-12_dp,&
                       KrtoH  = 1.E-12_dp,&
                       RbtoH  = 1.E-12_dp,&
                       SrtoH  = 1.E-12_dp,&
                       YtoH   = 1.E-12_dp,&
                       ZrtoH  = 1.E-12_dp,&
                       NbtoH  = 1.E-12_dp,&
                       MotoH  = 1.E-12_dp,&
                       TctoH  = 1.E-12_dp,&
                       RutoH  = 1.E-12_dp,&
                       RhtoH  = 1.E-12_dp,&
                       PdtoH  = 1.E-12_dp,&
                       AgtoH  = 1.E-12_dp,&
                       CdtoH  = 1.E-12_dp,&
                       IntoH  = 1.E-12_dp,&
                       SntoH  = 1.E-12_dp,&
                       SbtoH  = 1.E-12_dp,&
                       TetoH  = 1.E-12_dp,&
                       ItoH   = 1.E-12_dp,&
                       XetoH  = 1.E-12_dp,&
                       CstoH  = 1.E-12_dp,&
                       BatoH  = 1.E-12_dp,&
                       LatoH  = 1.E-12_dp,&
                       CetoH  = 1.E-12_dp,&
                       PrtoH  = 1.E-12_dp,&
                       NdtoH  = 1.E-12_dp,&
                       PmtoH  = 1.E-12_dp,&
                       SatoH  = 1.E-12_dp,&
                       EutoH  = 1.E-12_dp,&
                       GdtoH  = 1.E-12_dp,&
                       TbtoH  = 1.E-12_dp,&
                       DytoH  = 1.E-12_dp,&
                       HotoH  = 1.E-12_dp,&
                       ErtoH  = 1.E-12_dp,&
                       TmtoH  = 1.E-12_dp,&
                       YbtoH  = 1.E-12_dp,&
                       LutoH  = 1.E-12_dp,&
                       HftoH  = 1.E-12_dp,&
                       TatoH  = 1.E-12_dp,&
                       WtoH   = 1.E-12_dp,&
                       RetoH  = 1.E-12_dp,&
                       OstoH  = 1.E-12_dp,&
                       IrtoH  = 1.E-12_dp,&
                       PttoH  = 1.E-12_dp,&
                       AutoH  = 1.E-12_dp,&
                       HgtoH  = 1.E-12_dp,&
                       TltoH  = 1.E-12_dp,&
                       PbtoH  = 1.E-12_dp,&
                       BitoH  = 1.E-12_dp,&
                       PotoH  = 1.E-12_dp,&
                       AttoH  = 1.E-12_dp,&
                       RntoH  = 1.E-12_dp,&
                       FrtoH  = 1.E-12_dp,&
                       RatoH  = 1.E-12_dp,&
                       ActoH  = 1.E-12_dp,&
                       ThtoH  = 1.E-12_dp,&
                       PatoH  = 1.E-12_dp,&
                       UtoH   = 1.E-12_dp,&
                       NptoH  = 1.E-12_dp,&
                       PutoH  = 1.E-12_dp,&
                       AmtoH  = 1.E-12_dp,&
                       CmtoH  = 1.E-12_dp

!
real(dp)            :: diffuse_threshold   = 1E-3_dp
integer(dp)         :: diffuse_random_rays = 1
logical             :: log_diffuse         = .false.
!
real(dp)            :: initial_temperature    =  7500._dp
real(dp)            :: artificial_temperature = 10000._dp
real(dp)            :: filling_factor         = 1.0_dp

namelist/input_values/&
     abund_rel_to_h, add_heating,add_heating_filename,&
     artificially_fix_temperature, Arions, Artosolar,ArtoH,& 
     artificial_temperature,binary_output,&
     c_correct_on, charge_transfer, compute_temperature,Cions,Ctosolar,CtoH,&
     diffuse,diffuse_random, diffuse_random_rays,diffuse_threshold,&
     extra_heating, extra_radiation_field, external_chieta_program,&
     equidistant, external_hc_program, f_a, factor, filling_factor,&
     H_0, Heions, HetoH, Hetosolar,ifrit_filename, initial_file_spherical,&
     include_Ar, include_C, include_He, include_metals, include_N, include_Ne,&
     include_O, include_S, ifrit_filename,initial_redshift,&
     initial_temperature,l_cell,last_step,log_diffuse,&
     log_heating_cooling,lower_energy_limit,mergespec,mask_filename,&
     Neions, Netosolar,NetoH,nH_complete_scalar,nHI_fraction,no_source_filename,&
     Nions, Ntosolar,NtoH,num_threads,o_abundan,o_external, o_rad_force,&
     o_thermal_pressure,o_write_escape_fraction,old_solar_abundances,&
     omega_m,optabu_3d, Oions, Otosolar,OtoH, outward_only,& 
     periodic,photoheat_metals, planar_case_A, planar_case_B,&
     read_from_file, rep_output,&
     rp_max, start_ionized,Sions,source_filename,&
     spherical_case_A, spherical_case_B,&
     Stosolar,StoH, t_start,t_end, t_step_max,trace_expansion,use_mask,&
     verbose_stdout,x_max,y_max,z_max

 
contains

!================================================================================
  subroutine LOAD_PARAMETERS
!
!    Reads the namelist that controls the program
!    and does some initialization
!
! called by main program
!
  implicit none
!
  integer(i4b):: errstat
!===============================================================================
!  
  read(*,input_values)
!  
  x_min = -x_max
  y_min = -y_max
  z_min = -z_max
!
  if (abund_rel_to_h) then
     Heabundance = HetoH
     Cabundance  = CtoH
     Nabundance  = NtoH
     Oabundance  = OtoH
     Neabundance = NetoH
     Sabundance  = StoH
     Arabundance = ArtoH
  else
     if (old_solar_abundances) then
        Heabundance =  10._dp**(11.00_dp-12._dp)*Hetosolar
        Cabundance  =  10._dp**( 8.52_dp-12._dp)*Ctosolar
        Nabundance  =  10._dp**( 7.92_dp-12._dp)*Ntosolar
        Oabundance  =  10._dp**( 8.83_dp-12._dp)*Otosolar
        Neabundance =  10._dp**( 8.08_dp-12._dp)*Netosolar
        Sabundance  =  10._dp**( 7.33_dp-12._dp)*Stosolar
        Arabundance =  10._dp**( 6.40_dp-12._dp)*Artosolar    
     else
        Heabundance =  10._dp**(10.93_dp-12._dp)*Hetosolar
        Cabundance  =  10._dp**( 8.43_dp-12._dp)*Ctosolar
        Nabundance  =  10._dp**( 7.83_dp-12._dp)*Ntosolar
        Oabundance  =  10._dp**( 8.69_dp-12._dp)*Otosolar
        Neabundance =  10._dp**( 7.93_dp-12._dp)*Netosolar
        Sabundance  =  10._dp**( 7.12_dp-12._dp)*Stosolar
        Arabundance =  10._dp**( 6.83_dp-12._dp)*Artosolar
     end if
  end if
!
  if (no_source_filename) then
     read(*,*)
     read(*,*)
     read(*,*) num_sources
  else
     open(2, file=adjustl(trim(source_filename)))
     read(2,*)
     read(2,*) num_sources
  end if
!
! in the spherical case, there is only one source  
  if (spherical_case_A .or. spherical_case_B)&  
       num_sources=1
!
  allocate(source_x(1:num_sources), stat=errstat)  
  if (errstat .ne. 0) then
     stop 'fatal error: could not allocate source_x'
  end if
!  
  allocate(source_y(1:num_sources), stat=errstat)
  if (errstat .ne. 0) then
     stop 'fatal error: could not allocate source_y'
  end if
!  
  allocate(source_z(1:num_sources), stat=errstat)
  if (errstat .ne. 0) then
     stop 'fatal error: could not allocate source_z'
  end if
!  
  if (spherical_case_A .or. spherical_case_B ) then
     source_x(1) = 0
     source_y(1) = 0
     source_z(1) = -z_max
  end if

  allocate(R_RS(1:num_sources), stat=errstat)
  if (errstat .ne. 0) then
     stop 'fatal error: could not allocate R_RS'
  end if
!  
  allocate(spectra(1:num_sources), stat=errstat)
  if (errstat .ne. 0) then
     stop 'fatal error: could not allocate spectra'
  end if
!  
  allocate(start_activity(1:num_sources), stat=errstat)
  if (errstat .ne. 0) then
     stop 'fatal error: could not allocate start_activity'
  end if
!  
  allocate(lifetime(1:num_sources), stat=errstat)
  if (errstat .ne. 0) then
     stop 'fatal error: could not allocate lifetime'
  end if
!  
  if (no_source_filename) then
     do source_counter=1,num_sources
        read(*,*)
        read(*,*) source_x(source_counter), source_y(source_counter),&
             source_z(source_counter), &
             & R_RS(source_counter), spectra(source_counter),&
             start_activity(source_counter),&
             & lifetime(source_counter)
     end do
  else
     do source_counter=1,num_sources
        read(2,*)
        read(2,*) source_x(source_counter), source_y(source_counter),&
             source_z(source_counter), &
             & R_RS(source_counter), spectra(source_counter),&
             start_activity(source_counter),&
             & lifetime(source_counter)
     end do
     close(2)
  end if
!
 write(*,*)
 write(*,*) 'xrange :' , x_min , '...' , x_max
 write(*,*) 'yrange :' , y_min , '...' , y_max
 write(*,*) 'zrange :' , z_min , '...' , z_max
 write(*,*)
 write(*,*) 'Number of sources' , num_sources
 write(*,*)
!
 if (spherical_case_A .or. spherical_case_B) then
    write(*,*) 'R_RS  (solar radii) :', R_RS(1)
    write(*,*) 'spectrum :', trim(spectra(1))
    write(*,*) 'start_activity[yrs.] :', start_activity(1)
    write(*,*) 'lifetime[yrs.] :',lifetime(1)
    write(*,*)  
 else
    do source_counter=1, num_sources
       write(*,*) 'coordinates #',source_counter,' :',&
            source_x(source_counter), &
            &source_y(source_counter),&
            source_z(source_counter)
       write(*,*) 'R_RS  (solar radii) #',source_counter,' :', R_RS(source_counter)
       write(*,*) 'spectra  #',source_counter,' :', trim(spectra(source_counter))
       write(*,*) 'start_activity[yrs.] #',source_counter,' :',&
            start_activity(source_counter)
       write(*,*) 'lifetime[yrs.] #',source_counter,' :',lifetime(source_counter)
       write(*,*)
    end do
 end if
! 
 if (c_correct_on) then
    write(*,*) 'speed of light correction: on'
 else
    write(*,*) 'speed of light correction: off'
 end if
!
 write(*,*)
 write(*,*) 'spatial security factor :', f_a
 write(*,*)
 write(*,*) 'cell length in parsec :', l_cell
 write(*,*)
 write(*,*) 'timestep at the start of the program in s :',t_start
 write(*,*) 'total time of simulation in yrs. :', t_end
 write(*,*) 'temporal security factor :',factor  
 write(*,*)
!
 if (read_from_file) then
    write(*,*) 'ifrit file to be read :', ifrit_filename
 else
    write(*,*) 'initial density fraction of neutral hydrogen :',nHI_fraction
    write(*,*) 'total density of hydrogen :', nH_complete_scalar
 end if
! 
end subroutine LOAD_PARAMETERS
!-------------------------------------------------------------------------------
     
end module M_definitions

!-------------------------------------------------------------------------------
!                       ****END MODULE M_definitions****
!-------------------------------------------------------------------------------
