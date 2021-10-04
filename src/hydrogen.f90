!===============================================================================
!                         ****Module M_hydrogen****
!===============================================================================
!
! Computes the ionization and recombination rates for the considered ions, and
! updates the number densities. Furthermore the heating and cooling rates and 
! the temperature are determined 
!==============================================================================
!
module M_hydrogen
!--- routines:
!    INITIALIZE_HYDROGEN
!    INITIALIZE_HYDROGEN_FROM_FILE
!    RESUME_COMPUTATION
!    UPDATE_CHI
!    UPDATE_ETA
!    CALC_HYDROGEN
!    CALC_HELIUM2
!    CALC_METALS
!    UPDATE_ELECTRON_DENSITY
!    UPDATE_TEMPERATURE
!    UPDATE_THERMAL_PRESSURE
!    FFCOOLING
!    RECHII
!    RECHEII
!    RECHEIII
!    SEATON
!    RADIATIVE_REC0
!    RADIATIVE_REC3
!    COL_FL_COOLING
!    DIELECTRONIC_REC
!    HYDREC
!    HELREC
!    PHOTOINTEGRAL
!    PHOTOHEATINTEGRAL
!    DIELECTRONIC_HOT
!    DIELECTRONIC_COOL
!    EXPAND_SPACE
!      
  use M_data_types
!
  implicit none
! 
  real(dp)     :: nHII_t , x_co, y_co, z_co, nHI_co
  integer(i4b) :: rep=0
!  
  real(dp), save  :: dull_dummy 
contains

!============================================================================
  subroutine INITIALIZE_MASK
!used for masking out
    use M_definitions              ,only: x_max, y_max, z_max, mask_filename
    use M_grid_memory              ,only: mask
    implicit none
    integer(i4b)                       :: x, y, z
    integer(i4b)                       :: dummy1, dummy2, dummy3
!============================================================================
    open (unit=500, file=mask_filename)
    read(500,*) dummy1, dummy2, dummy3
    do z=-z_max, +z_max
       do y=-y_max, +y_max
          do x=-x_max, +x_max
             read(500,*) mask(x,y,z)
          end do
       end do
    end do
    close(500)
  end subroutine INITIALIZE_MASK

!=============================================================================
  subroutine INITIALIZE_HEATING_POLY
    use M_definitions              ,only: x_max, y_max, z_max,&
                                          add_heating_filename
    use M_grid_memory              ,only: add_heating_poly
    implicit none
!
    integer(i4b)                       :: x,y,z  
!=============================================================================
!
    open(unit=600, file=add_heating_filename)
    do x=-x_max,+x_max
       do y=-y_max, +y_max
          do z=-z_max, +z_max
             read(600,*) add_heating_poly(1,x,y,z),&
                         add_heating_poly(2,x,y,z),&
                         add_heating_poly(3,x,y,z),&
                         add_heating_poly(4,x,y,z)
          end do
       end do
    end do
    close(600)
end subroutine INITIALIZE_HEATING_POLY
!

!=============================================================================
  subroutine INITIALIZE_ABUNDANCES
! Initializes the abundanes (in terms of the number density) relative to
! the number density of hydrogen
! The hydrogen number density (and the ionization structure) will be 
! defined in the subroutines INITIALIZE_HYDROGEN and 
! INITIALIZE_HYDROGEN_FROM FILE
! It is _not_ called by RESUME_COMPUTATION

    use M_grid_memory             ,only: nH_complete, nHe_complete,&
                                         nC, nN, nO, nNe, nS,nAr,&
                                         element_abundances
    use M_definitions             ,only: Heabundance, Cabundance,&
                                         Nabundance, Oabundance,&
                                         Neabundance, Sabundance,&
                                         Arabundance,&
                                         o_abundan, optabu_3d,&
                                         include_metals
!
    implicit none
!==============================================================================
    if (o_abundan) then
       call read_abundances
       if (optabu_3d) then
            nHe_complete=nH_complete*element_abundances(2,:,:,:)
            if (include_metals) then
               nC (:,:,:,0)=nH_complete*element_abundances(6,:,:,:)
               nN (:,:,:,0)=nH_complete*element_abundances(7,:,:,:)
               nO (:,:,:,0)=nH_complete*element_abundances(8,:,:,:)
               nNe(:,:,:,0)=nH_complete*element_abundances(10,:,:,:)
               nS (:,:,:,0)=nH_complete*element_abundances(16,:,:,:)
               nAr(:,:,:,0)=nH_complete*element_abundances(18,:,:,:)
            end if
         else
            nHe_complete=nH_complete*element_abundances(2,0,0,0)
            if (include_metals) then
               nC (:,:,:,0)=nH_complete*element_abundances(6,0,0,0)
               nN (:,:,:,0)=nH_complete*element_abundances(7,0,0,0)
               nO (:,:,:,0)=nH_complete*element_abundances(8,0,0,0)
               nNe(:,:,:,0)=nH_complete*element_abundances(10,0,0,0)
               nS (:,:,:,0)=nH_complete*element_abundances(16,0,0,0)
               nAr(:,:,:,0)=nH_complete*element_abundances(18,0,0,0)
               print*, "nAr =", nAr(0,0,0,0)
            end if
         end if
      else
         nHe_complete=nH_complete*Heabundance
         if (include_metals) then
            nC (:,:,:,0)=nH_complete*Cabundance
            nN (:,:,:,0)=nH_complete*Nabundance
            nO (:,:,:,0)=nH_complete*Oabundance
            nNe(:,:,:,0)=nH_complete*Neabundance
            nS (:,:,:,0)=nH_complete*Sabundance
            nAr(:,:,:,0)=nH_complete*Arabundance
         end if
      end if

      return
!
    contains
!      
      subroutine read_abundances
        use M_definitions      , only: x_max, y_max, z_max
        implicit none
! 
        integer (i4b)  :: x,y,z, element, ios
!
        open (20, file='ABUNDAN', form='UNFORMATTED', iostat=ios, status='OLD')
        rewind 20
        if (optabu_3d) then
           do z=-z_max,+z_max
              do y=-y_max, +y_max
                 do x=-x_max,+x_max
                    do element=1, 30
                       read(20, ERR=10)  element_abundances(:,x,y,z)
                       print*, element
                    end do
                 end do
              end do
           end do
        else
           read(20, ERR=10)  element_abundances(:,0,0,0)
           write(*,*) 'Element', 'Abundance'
        end if
        close(20)
        return
10      write (*,*) 'SOMETHING WENT WRONG WHEN READING ABUNDAN'
        write(*,*) 'IOS=', ios
        stop
      end subroutine read_abundances
!
  end subroutine INITIALIZE_ABUNDANCES

!=============================================================================
  subroutine INITIALIZE_HYDROGEN
!
!    Initializes the vaiables defined in M_grid_definitions which describe
!    the initial conditions concerning the number densities of all ions
!    (not only hydrogen - the name is for historical reasons) 
!    as well as the temperature structure. This routine is called, if a
!    a simulation assumes a homogeneous density distribution
!
!   called by th main program
!    
    use M_natural_constants,  only: nc_boltzmann    
    use M_definitions,        only: nH_complete_scalar, start_ionized,&
                                    nHI_fraction, include_metals,&
                                     initial_temperature,&
                                    trace_expansion,&
                                    initial_redshift, redshift_old
    use M_grid_memory,        only: ne, nH_complete,nHI,&
                                    nHII, nHI_old,nHe_complete,&
                                    nHeI, nHeII, nHeIII,&
                                    nC, nN, nO, nNe, nS,nAr,&
                                    temperature, energycontent
!
  implicit none
!============================================================================== 
!
    print*,"initialize_hydrogen"
    nH_complete=nH_complete_scalar
!
    call initialize_abundances
!
    if (start_ionized) then
       nHI     = nHI_fraction*nH_complete_scalar 
       nHII    = (1._dp-nHI_fraction)*nH_complete_scalar
       nHI_old =nH_complete       
!
       nHeIII       = nHe_complete*(1._dp-2.E-10_dp)
       nHeII        = nHe_complete*1E-10_dp
       nHeI         = nHe_complete*1E-10_dp
!       
       if (include_metals) then          
!
          nC(:,:,:,4)=nC(:,:,:,0)*(1._dp-3*1e-10_dp)
          nC(:,:,:,3)=nC(:,:,:,0)*(1e-10_dp)
          nC(:,:,:,2)=nC(:,:,:,0)*(1e-10_dp)
          nC(:,:,:,1)=nC(:,:,:,0)*(1e-10_dp)
!                    
          nN(:,:,:,4)=nN(:,:,:,0)*(1._dp-3*1e-10_dp)
          nN(:,:,:,3)=nN(:,:,:,0)*(1e-10_dp)
          nN(:,:,:,2)=nN(:,:,:,0)*(1e-10_dp)
          nN(:,:,:,1)=nN(:,:,:,0)*(1e-10_dp)
!                  
          nO(:,:,:,4)=nO(:,:,:,0)*(1._dp-3*1e-10_dp)
          nO(:,:,:,3)=nO(:,:,:,0)*(1e-10_dp)
          nO(:,:,:,2)=nO(:,:,:,0)*(1e-10_dp)
          nO(:,:,:,1)=nO(:,:,:,0)*(1e-10_dp)
!          
          nNe(:,:,:,4)=nNe(:,:,:,0)*(1._dp-3*1e-10_dp)
          nNe(:,:,:,3)=nNe(:,:,:,0)*(1e-10_dp)
          nNe(:,:,:,2)=nNe(:,:,:,0)*(1e-10_dp)
          nNe(:,:,:,1)=nNe(:,:,:,0)*(1e-10_dp)
!          
          nS(:,:,:,4)=nS(:,:,:,0)*(1._dp-3*1e-10_dp)
          nS(:,:,:,3)=nS(:,:,:,0)*(1e-10_dp)
          nS(:,:,:,2)=nS(:,:,:,0)*(1e-10_dp)
          nS(:,:,:,1)=nS(:,:,:,0)*(1e-10_dp)
!
          nAr(:,:,:,4)=nAr(:,:,:,0)*(1._dp-3*1e-10_dp)
          nAr(:,:,:,3)=nAr(:,:,:,0)*(1e-10_dp)
          nAr(:,:,:,2)=nAr(:,:,:,0)*(1e-10_dp)
          nAr(:,:,:,1)=nAr(:,:,:,0)*(1e-10_dp)
      end if
    else
!
       nHI     = nHI_fraction*nH_complete_scalar 
       nHII    = (1._dp-nHI_fraction)*nH_complete_scalar
       nHI_old = nH_complete
!
       nHeI         = nHe_complete*(1._dp-2.E-10_dp)
       nHeII        = nHe_complete*1E-10_dp
       nHeIII       = nHe_complete*1E-10_dp
!
       if (include_metals) then
!                    
          nC(:,:,:,1) = nC(:,:,:,0)*(1._dp-3*1E-10_dp)
          nC(:,:,:,2) = nC(:,:,:,0)*(1E-10_dp)
          nC(:,:,:,3) = nC(:,:,:,0)*(1E-10_dp)
          nC(:,:,:,4) = nC(:,:,:,0)*(1E-10_dp)
!
          nN(:,:,:,1) = nN(:,:,:,0)*(1._dp-3*1E-10_dp)
          nN(:,:,:,2) = nN(:,:,:,0)*(1E-10_dp)
          nN(:,:,:,3) = nN(:,:,:,0)*(1E-10_dp)
          nN(:,:,:,4) = nN(:,:,:,0)*(1E-10_dp)
!                    
          nO(:,:,:,1) = nO(:,:,:,0)*(1._dp-3*1E-10_dp)
          nO(:,:,:,2) = nO(:,:,:,0)*(1E-10_dp)
          nO(:,:,:,3) = nO(:,:,:,0)*(1E-10_dp)
          nO(:,:,:,4) = nO(:,:,:,0)*(1E-10_dp)
!          
          nNe(:,:,:,1) = nNe(:,:,:,0)*(1._dp-3*1E-10_dp)
          nNe(:,:,:,2) = nNe(:,:,:,0)*(1E-10_dp)
          nNe(:,:,:,3) = nNe(:,:,:,0)*(1E-10_dp)
          nNe(:,:,:,4) = nNe(:,:,:,0)*(1E-10_dp)
!          
          nS(:,:,:,1) = nS(:,:,:,0)*(1._dp-3*1E-10_dp)
          nS(:,:,:,2) = nS(:,:,:,0)*(1E-10_dp)
          nS(:,:,:,3) = nS(:,:,:,0)*(1E-10_dp)
          nS(:,:,:,4) = nS(:,:,:,0)*(1E-10_dp)
!
          nAr(:,:,:,1) = nAr(:,:,:,0)*(1._dp-3*1E-10_dp)
          nAr(:,:,:,2) = nAr(:,:,:,0)*(1E-10_dp)
          nAr(:,:,:,3) = nAr(:,:,:,0)*(1E-10_dp)
          nAr(:,:,:,4) = nAr(:,:,:,0)*(1E-10_dp)
       end if
    end if
!        
    ne=(nH_complete-nHI)+nHe_complete*(2E-10_dp)
!    
    temperature=initial_temperature
    if (include_metals) then
       energycontent=3._dp/2._dp*nc_boltzmann*temperature*&
            (ne+nH_complete+nHe_complete+&
            nC(:,:,:,0)+ nN(:,:,:,0)+nO(:,:,:,0)+ nNe(:,:,:,0)+nS(:,:,:,0))
    else       
       energycontent=3._dp/2._dp*nc_boltzmann*temperature*&
            (ne+nH_complete+nHe_complete)
     end if  
!    
    if (trace_expansion) then
       redshift_old=initial_redshift
    end if
!
  end subroutine INITIALIZE_HYDROGEN


  
!================================================================================  
  subroutine INITIALIZE_HYDROGEN_FROM_FILE (ifrit_filename)
!
!    Like, INITIALIZE_HYDROGEN, but the density structure is read
!    from the file ifrit_filename. Appropriate for simulations with
!    inhomogeneous density structure
!
! called by the main program    
!    
    use M_natural_constants,  only: nc_boltzmann
    use M_definitions,        only: x_max, y_max, z_max,x_min, y_min, z_min,&
                                    start_ionized,&
                                    include_metals,&
                                    Heabundance, initial_temperature,&
                                    trace_expansion,&
                                    initial_redshift,&
                                    redshift_old, points
    use M_grid_memory,        only: ne, nH_complete,nHI,&
                                    nHII, nHe_complete,&
                                    nHeI, nHeII, nHeIII,&
                                    nC, nN, nO, nNe, nS, nAr,&
                                    temperature, energycontent
!
  implicit none
!
  character (len=3000):: ifrit_filename
  integer:: x_range, y_range, z_range
  integer:: x_counter, y_counter, z_counter
!================================================================================
!
    open(1,file=ifrit_filename)    
    read(1,*) x_range, y_range, z_range
    if (x_range .ne. x_max-x_min+1 .or. &
         y_range .ne. y_max-y_min+1 .or. &
         z_range .ne. z_max-z_min+1) then
       stop 'error: file does not match required array size'
    end if
    
    do z_counter=-(z_range/2),z_range/2
       do y_counter=-(y_range/2), y_range/2
          do x_counter=-(z_range/2), z_range/2
             read(1,*) nH_complete(x_counter, y_counter, z_counter),&
                  nHI(x_counter,y_counter,z_counter)
             nHII(x_counter,y_counter,z_counter)=&
                  nH_complete(x_counter,y_counter,z_counter)-&
                  nHI(x_counter,y_counter,z_counter)
          end do
       end do
    end do
!
    call initialize_abundances
!
    if (start_ionized) then
!
       nHeIII=nHe_complete*(1._dp-2E-10_dp)
       nHeII=nHe_complete*1E-10_dp
       nHeI=nHe_complete- nHeII-nHeIII
!
       ne=nHII+nHeII+2*nHeIII
       temperature=10000._dp
!       
       if (include_metals) then
!          
          nC (:,:,:,4) = nC(:,:,:,0)*(1.-3*1E-10_dp)
          nC (:,:,:,3) = nC (:,:,:,0)*(1E-10_dp)
          nC (:,:,:,2) = nC (:,:,:,0)*(1E-10_dp)
          nC (:,:,:,1) = nC (:,:,:,0)*(1E-10_dp)
!          
          nN (:,:,:,4) = nN (:,:,:,0)*(1.-3*1E-10_dp)
          nN (:,:,:,3) = nN (:,:,:,0)*(1E-10_dp)
          nN (:,:,:,2) = nN (:,:,:,0)*(1E-10_dp)
          nN (:,:,:,1) = nN (:,:,:,0)*(1E-10_dp)
!
          nO (:,:,:,4) = nO (:,:,:,0)*(1._dp-3*1E-10)
          nO (:,:,:,3) = nO (:,:,:,0)*(1E-10_dp)
          nO (:,:,:,2) = nO (:,:,:,0)*(1E-10_dp)
          nO (:,:,:,1) = nO (:,:,:,0)*(1E-10_dp)
!                    
          nNe(:,:,:,4) = nNe(:,:,:,0)*(1._dp-3*1E-10_dp)
          nNe(:,:,:,3) = nNe(:,:,:,0)*(1E-10_dp)
          nNe(:,:,:,2) = nNe(:,:,:,0)*(1E-10_dp)
          nNe(:,:,:,1) = nNe(:,:,:,0)*(1E-10_dp)
!          
          nS (:,:,:,4) = nS (:,:,:,0)*(1._dp-3*1E-10_dp)
          nS (:,:,:,3) = nS (:,:,:,0)*(1E-10_dp)
          nS (:,:,:,2) = nS (:,:,:,0)*(1E-10_dp)
          nS (:,:,:,1) = nS (:,:,:,0)*(1E-10_dp)
!
          nAr(:,:,:,4) = nAr(:,:,:,0)*(1._dp-3*1E-10_dp)
          nAr(:,:,:,3) = nAr(:,:,:,0)*(1E-10_dp)
          nAr(:,:,:,2) = nAr(:,:,:,0)*(1E-10_dp)
          nAr(:,:,:,1) = nAr(:,:,:,0)*(1E-10_dp)
!
       end if
!
    else
       nHe_complete=nH_complete*Heabundance 
       nHeI=nHe_complete*(1._dp-2E-10_dp)
       nHeII=nHe_complete*1E-10_dp
       nHeIII=nHe_complete- nHeI-nHeII
!
       ne=nHII+nHeII+2*nHeIII
       temperature=initial_temperature
!       
       if (include_metals) then
          nC (:,:,:,1) = nC (:,:,:,0)*(1._dp-3*1E-10)
          nC (:,:,:,2) = nC (:,:,:,0)*(1E-10)
          nC (:,:,:,3) = nC (:,:,:,0)*(1E-10)
          nC (:,:,:,4) = nC (:,:,:,0)*(1E-10)
!          
          nN (:,:,:,1) = nN (:,:,:,0)*(1._dp-3*1E-10_dp)
          nN (:,:,:,2) = nN (:,:,:,0)*(1E-10_dp)
          nN (:,:,:,3) = nN (:,:,:,0)*(1E-10_dp)
          nN (:,:,:,4) = nN (:,:,:,0)*(1E-10_dp)
!          
          nO (:,:,:,1) = nO (:,:,:,0)*(1._dp-3*1E-10_dp)
          nO (:,:,:,2) = nO (:,:,:,0)*(1E-10_dp)
          nO (:,:,:,3) = nO (:,:,:,0)*(1E-10_dp)
          nO (:,:,:,4) = nO (:,:,:,0)*(1E-10_dp)
!                    
          nNe(:,:,:,1) = nNe (:,:,:,0)*(1._dp-3*1E-10_dp)
          nNe(:,:,:,2) = nNe(:,:,:,0)*(1E-10_dp)
          nNe(:,:,:,3) = nNe(:,:,:,0)*(1E-10_dp)
          nNe(:,:,:,4) = nNe(:,:,:,0)*(1E-10_dp)
!          
          nS (:,:,:,1) = nS(:,:,:,0)*(1._dp-3*1E-10_dp)
          nS (:,:,:,2) = nS(:,:,:,0)*(1E-10_dp)
          nS (:,:,:,3) = nS(:,:,:,0)*(1E-10_dp)
          nS (:,:,:,4) = nS(:,:,:,0)*(1E-10_dp)
!
          nAr(:,:,:,1) = nAr(:,:,:,0)*(1._dp-3*1E-10_dp)
          nAr(:,:,:,2) = nAr(:,:,:,0)*(1E-10_dp)
          nAr(:,:,:,3) = nAr(:,:,:,0)*(1E-10_dp)
          nAr(:,:,:,4) = nAr(:,:,:,0)*(1E-10_dp)
       end if
!
    end if
!
!    print*, "CHI"
!    do frequency_counter=0, points-1
!       print*, nug(frequency_counter), chi(0,0,0,frequency_counter)
!    end do
!    
    if (trace_expansion) then
       redshift_old=initial_redshift
    end if
    
    temperature=7500._dp
    if (include_metals) then
       energycontent=3._dp/2._dp*nc_boltzmann*temperature*&
            (ne+nH_complete+nHe_complete+nC(:,:,:,0)+&
            nN(:,:,:,0)+nO(:,:,:,0)+nNe(:,:,:,0)+nS(:,:,:,0))
       else
       energycontent=3._dp/2._dp*nc_boltzmann*temperature*&
            (ne+nH_complete+nHe_complete)
       end if
    
    print*,"END INITIALIZE"
  end subroutine INITIALIZE_HYDROGEN_FROM_FILE



!================================================================================ 
  subroutine RESUME_COMPUTATION
!
!  The program can be restarted by providing the output files of a former
!  program run.
!  This is useful in the case of unexpected interruptions like computer crashes 
!
    use M_natural_constants,  only: nc_boltzmann
    use M_definitions,        only: x_max, y_max, z_max, command_line_options,&
                                    t_all, t_step, t_start,&
                                    include_metals, compute_temperature
    use M_grid_memory,        only: nHI, nHII, nH_complete,nHeI, nHeII,&
                                    nHeIII, nHe_complete,ne,nC,nN,nO, nNe,& 
                                    nS,nAr,temperature,energycontent
!      
  implicit none
!
    integer             :: xrange, yrange, zrange
    integer             :: x, y, z
    integer             :: string_length
    real(dp)            :: dummy1, dummy2    
    character (len=3000):: trunk
!================================================================================
!
    trunk=adjustr(trim(adjustl(command_line_options(2))))
    string_length=len_trim(trunk)
    print*,trunk
    print*,trunk(1:string_length-5) 
    read(trunk(1:string_length-5),*) t_all
!
    t_all=t_all*(365.25_dp*24*3600)
    t_step=t_start
    print*, "T_step:", t_step
!
    open(1001, file=adjustr(trim(adjustl(trunk)))//'.H.txt')
    print*,adjustr(trim(adjustl(trunk)))//'.H.txt' 
    read(1001,*) xrange, yrange, zrange
    if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
         zrange .ne. 2*z_max+1) &
         stop 'Data sets are incompatible'
    do z=-z_max, +z_max
       do y=-y_max,+y_max
          do x=-x_max, +x_max
             read(1001,*) nHI(x,y,z), nHII(x,y,z), nH_complete(x,y,z)
          end do
       end do
    end do
    close (1001)
!
    open(1002, file=adjustr(trim(adjustl(trunk)))//'.He.txt')
    read(1002,*) xrange, yrange, zrange
    if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
         zrange  .ne. 2*z_max+1) &
         stop 'Data sets are incompatible'
    do z=-z_max, +z_max
       do y=-y_max,+y_max
          do x=-x_max, +x_max
             read(1002,*) nHeI(x,y,z), nHeII(x,y,z), nHeIII(x,y,z),&
                  nHe_complete(x,y,z)
          end do
       end do
    end do
    close(1002)
!
    do z=-z_max, +z_max
       do y=-y_max, +y_max
          do x=-x_max, +x_max
             ne(x,y,z)=nHII(x,y,z)+nHeII(x,y,z)+2*nHeII(x,y,z)
          end do
       end do
    end do
!
    if (include_metals) then
       open(1003, file=adjustr(trim(adjustl(trunk)))//'.C.txt')
       open(1004, file=adjustr(trim(adjustl(trunk)))//'.N.txt')
       open(1005, file=adjustr(trim(adjustl(trunk)))//'.O.txt')
       open(1006, file=adjustr(trim(adjustl(trunk)))//'.Ne.txt')
       open(1007, file=adjustr(trim(adjustl(trunk)))//'.S.txt')
       open(1008, file=adjustr(trim(adjustl(trunk)))//'.Ar.txt')

!
       read(1003,*) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1004,*) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1005,*) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1006,*) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1007,*) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1008,*) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'

!
       do z=-z_max, +z_max
          do y=-y_max,+y_max
             do x=-x_max, +x_max
                read(1003,*) nC(x,y,z,1), nC(x,y,z,2), nC(x,y,z,3),&
                     nC(x,y,z,4),nC(x,y,z,0)
             end do
          end do
       end do
       do z=-z_max,+z_max
          do y=-y_max,+y_max
             do x=-x_max,+x_max
                read(1004,*) nN(x,y,z,1), nN(x,y,z,2), nN(x,y,z,3),&
                     nN(x,y,z,4), nN(x,y,z,0)
             end do
          end do
       end do
       do z=-z_max,+z_max
          do y=-y_max,+y_max
             do x=-x_max,+x_max
                read(1005,*) nO(x,y,z,1), nO(x,y,z,2), nO(x,y,z,3),&
                     nO(x,y,z,4), nO(x,y,z,0)
             end do
          end do
       end do
       do z=-z_max,+z_max
          do y=-y_max,+y_max
             do x=-x_max,+x_max
                read(1006,*) nNe(x,y,z,1), nNe(x,y,z,2), nNe(x,y,z,3),&
                     nNe(x,y,z,4), nNe(x,y,z,0)
             end do
          end do
       end do
       do z=-z_max,+z_max
          do y=-y_max,+y_max
             do x=-x_max,+x_max
                read(1007,*) nS(x,y,z,1), nS(x,y,z,2), nS(x,y,z,3),&
                     nS(x,y,z,4), nS(x,y,z,0)
             end do
          end do
       end do
       do z=-z_max,+z_max
          do y=-y_max,+y_max
             do x=-x_max,+x_max
                read(1008,*) nAr(x,y,z,1), nAr(x,y,z,2), nAr(x,y,z,3),&
                     nAr(x,y,z,4), nAr(x,y,z,0)
             end do
          end do
       end do
       do z=-z_max,+z_max
          do y=-y_max,+y_max
             do x=-x_max,+x_max
                ne(x,y,z)=ne(x,y,z)+&
                     1*(nC(x,y,z,2)+nN(x,y,z,2)+ nO(x,y,z,2)+&
                     nNe(x,y,z,2)+ nS(x,y,z,2)+nAr(x,y,z,2))+&
                     2*(nC(x,y,z,3)+nN(x,y,z,3)+ nO(x,y,z,3)+&
                     nNe(x,y,z,3)+ nS(x,y,z,3)+nAr(x,y,z,3))+&
                     3*(nC(x,y,z,4)+nN(x,y,z,4)+ nO(x,y,z,4)+&
                     nNe(x,y,z,4)+ nS(x,y,z,4)+nAr(x,y,z,4))
              end do
          end do
       end do
!
       close (1003)
       close (1004)
       close (1005)
       close (1006)
       close (1007)
    end if
!
    if(compute_temperature) then
       print*, trunk
       print*,adjustr(trim(adjustl(trunk)))//'.temp.txt'
       open(1050, file=adjustr(trim(adjustl(trunk)))//'.temp.txt')
       read(1050,*) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       do z= -z_max, +z_max
          do y=-y_max, +y_max
             do x=-x_max, +x_max
                read(1050,*) dummy1,temperature(x,y,z),dummy2
                energycontent(x,y,z)=3._dp/2._dp*nc_boltzmann*temperature(x,y,z)*&
                     (ne(x,y,z)+nH_complete(x,y,z)+nHe_complete(x,y,z)+&
                     nC(x,y,z,0)+nN(x,y,z,0)+nO(x,y,z,0)+&
                     nNe(x,y,z,0)+nS(x,y,z,0))
             end do
          end do
       end do
       close (1050)
    end if     
!
  end subroutine RESUME_COMPUTATION
!  
!================================================================================ 
  subroutine RESUME_COMPUTATION_FROM_BINARY
!
!  The program can be restarted by providing the output files of a former
!  program run.
!  This is useful in the case of unexpected interruptions like computer crashes 
!
    use M_natural_constants,  only: nc_boltzmann
    use M_definitions,        only: x_max, y_max, z_max, command_line_options,&
                                    t_all, t_step, t_start,&
                                    include_metals, compute_temperature
    use M_grid_memory,        only: nHI, nHII, nH_complete,nHeI, nHeII,&
                                    nHeIII, nHe_complete,ne,nC,nN,nO, nNe,& 
                                    nS,nAr,temperature,energycontent,&
                                    nHIf, nHIIf, nH_completef,nHeIf, nHeIIf,&
                                    nHeIIIf, nHe_completef,nef,nCf,nNf,&
                                    nOf, nNef,nSf,nArf,temperaturef,&
                                    energycontentf
!      
  implicit none
!
    integer             :: xrange, yrange, zrange
    integer             :: x, y, z
    integer             :: string_length
!   The values in the file are single-precision floating point numbers
    real(sp), dimension (:,:,:) , allocatable  :: dummy1f, dummy2f
    integer                                    :: allocflag
    character (len=3000):: trunk
!==============================================================================
!
    allocate(dummy1f(-x_max:+x_max,-y_max:+y_max,-z_Max:+z_max),stat=allocflag)
    allocate(dummy2f(-x_max:+x_max,-y_max:+y_max,-z_Max:+z_max),stat=allocflag)
    if (allocflag .ne. 0)&
         stop 'RESUME_COMPUTATION_FROM_BINARY: error allocating memory'

    trunk=adjustr(trim(adjustl(command_line_options(2))))
    string_length=len_trim(trunk)
    print*,trunk
    print*,trunk(1:string_length-5) 
    read(trunk(1:string_length-5),*) t_all
!
    t_all=t_all*(365.25_dp*24*3600)
    t_step=t_start
    print*, "T_step:", t_step
!
    open(1001, file=adjustr(trim(adjustl(trunk)))//'.H.bin',&
         form='unformatted')
    print*,adjustr(trim(adjustl(trunk)))//'.H.bin' 
    read(1001) xrange, yrange, zrange
    if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
         zrange .ne. 2*z_max+1) &
         stop 'Data sets are incompatible'
    read(1001) nHIf
    nHI=real(nHIf,dp)
    read(1001) nHIIf
    nHII=real(nHIIf,dp)
    read(1001) nH_completef
    nH_complete=real(nH_complete,dp)
    close (1001)
!
    open(1002, file=adjustr(trim(adjustl(trunk)))//'.He.bin',&
         form='unformatted')
    read(1002) xrange, yrange, zrange
    if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
         zrange  .ne. 2*z_max+1) &
         stop 'Data sets are incompatible'
    read(1002) nHeIf
    nHeI=real(nHeIf,dp)
    read(1002) nHeIIf
    nHeII=real(nHeIIf,dp)
    read(1002) nHeIIIf
    nHeIII=real(nHeIIIf,dp)
    read(1002) nHe_completef
    nHe_complete=real(nHe_Completef,dp) 
    close(1002)
!
    ne=nHII+nHeII+2*nHeII
!
    if (include_metals) then
       open(1003, file=adjustr(trim(adjustl(trunk)))//'.C.bin',&
            form='unformatted')
       open(1004, file=adjustr(trim(adjustl(trunk)))//'.N.bin',&
            form='unformatted')
       open(1005, file=adjustr(trim(adjustl(trunk)))//'.O.bin',&
            form='unformatted')
       open(1006, file=adjustr(trim(adjustl(trunk)))//'.Ne.bin',&
            form='unformatted')
       open(1007, file=adjustr(trim(adjustl(trunk)))//'.S.bin',&
            form='unformatted')
       open(1008, file=adjustr(trim(adjustl(trunk)))//'.Ar.bin',&
            form='unformatted')

!
       read(1003) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1004) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1005) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1006) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1007) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1008) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'

!
       read(1003) nCf(:,:,:,1)
       read(1003) nCf(:,:,:,2)
       read(1003) nCf(:,:,:,3)
       read(1003) nCf(:,:,:,4)
       read(1003) nCf(:,:,:,0)
       nN=real(nC,dp)
!
       read(1004) nNf(:,:,:,1)
       read(1004) nNf(:,:,:,2)
       read(1004) nNf(:,:,:,3)
       read(1004) nNf(:,:,:,4)
       read(1004) nNf(:,:,:,0)
       nN=real(nNF,dp)
                
       read(1005) nOf(:,:,:,1)
       read(1005) nOf(:,:,:,2)
       read(1005) nOf(:,:,:,3)
       read(1005) nOf(:,:,:,4)
       read(1005) nOf(:,:,:,0)
       nO=real(nOf,dp)
!
       read(1006) nNef(:,:,:,1)
       read(1006) nNef(:,:,:,2)
       read(1006) nNef(:,:,:,3)
       read(1006) nNef(:,:,:,4)
       read(1006) nNef(:,:,:,0)
       nNe=real(nNef,dp)
!
       read(1007) nSf(:,:,:,1)
       read(1007) nSf(:,:,:,2)
       read(1007) nSf(:,:,:,3)
       read(1007) nSf(:,:,:,4)
       read(1007) nSf(:,:,:,0)
       nS=real(nSf,dp)
!
       read(1008) nArf(:,:,:,1)
       read(1008) nArf(:,:,:,2)
       read(1008) nArf(:,:,:,3)
       read(1008) nArf(:,:,:,4)
       read(1008) nArf(:,:,:,0)
       nAr=real(nArf,dp)
!
       ne=ne+&
            1*(nC (:,:,:,2)+nN(:,:,:,2)+ nO(:,:,:,2) +&
              nNe (:,:,:,2)+nS(:,:,:,2)+nAr(:,:,:,2))+&
            2*(nC (:,:,:,3)+nN(:,:,:,3)+ nO(:,:,:,3) +&
               nNe(:,:,:,3)+nS(:,:,:,3)+nAr(:,:,:,3))+&
            3*(nC (:,:,:,4)+nN(:,:,:,4)+ nO(:,:,:,4)+&
               nNe(:,:,:,4)+nS(:,:,:,4)+nAr(:,:,:,4))
!
       close (1003)
       close (1004)
       close (1005)
       close (1006)
       close (1007)
       close (1008)
    end if
    !
    if(compute_temperature) then
       print*, trunk
       print*,adjustr(trim(adjustl(trunk)))//'.temp.bin'
       open(1050, file=adjustr(trim(adjustl(trunk)))//'.temp.bin',&
            form='unformatted')
       read(1050) xrange, yrange, zrange
       if (xrange .ne. 2*x_max+1 .or. yrange .ne. 2*y_max+1 .or.&
            zrange .ne. 2*z_max+1) &
            stop 'Data sets are incompatible'
       read(1050) dummy1f
       read(1050) temperaturef
       read(1050) dummy2f
       energycontent=3._dp/2._dp*nc_boltzmann*temperature*&
            (ne+nH_complete+nHe_complete+&
            nC(:,:,:,0)+nN(:,:,:,0)+nO(:,:,:,0)+&
            nNe(:,:,:,0)+nS(:,:,:,0))
       close (1050)
    end if
!
    deallocate(dummy1f,stat=allocflag)
    deallocate(dummy2f,stat=allocflag)
    if (allocflag .ne. 0)&
         stop 'RESUME_COMPUTATION_FROM_BINARY: error de-allocating memory'
!
  end subroutine RESUME_COMPUTATION_FROM_BINARY


!===============================================================================  
  subroutine UPDATE_CHI
! 
!    This routine updates the frequency_dependent opacitied for each 
!    of the grid points
!
!  called by the main program
!
    use M_natural_constants,  only: nc_light   
    use M_definitions,        only: x_max, y_max, z_max, points,&
                                      include_metals
    use M_grid_memory,        only: chi, nHI,  nHeI, nHeII,&
                                    nC, nN, nO, nNe, nS, nAr
    use M_data_input,         only: nug
!
  implicit none
!
    integer:: x,y,z,frequency_counter
!===============================================================================
!
    !$OMP PARALLEL PRIVATE(x,y, frequency_counter)
    !$OMP DO
    do z=-z_max, +z_max
       do y=-y_max,+y_max
          do x=-x_max, +x_max

             do frequency_counter=0,points-1
!                
                chi(x, y, z, frequency_counter)=&       !HI
                     nHI(x,y,z)*&
                     seaton(nug(frequency_counter),&
                     1.0973731569E+05_dp*2.99792458E+10_dp,&
                     6.3E-18_dp ,1.34_dp, 2.99_dp)    
!               
                chi(x, y, z, frequency_counter)=&       !HeI
                     chi(x, y, z, frequency_counter)+&
                     nHeI(x,y,z)*&
                     seaton(nug(frequency_counter),&
                     198305.0_dp*2.99792458E+10_dp,&
                     7.83E-18_dp ,1.66_dp, 2.05_dp)
!                               
                chi(x, y, z, frequency_counter)=&       !HeII
                     chi(x, y, z, frequency_counter)+&
                     nHeII(x,y,z)*&
                     seaton(nug(frequency_counter),438908.7_dp*2.99792458E+10_dp,&
                     1.58E-18_dp ,1.34_dp, 2.99_dp)    
!                
                if (include_metals) then

! Carbon
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nC(x,y,z,1)*&
                        seaton(nug(frequency_counter),90833.1_dp*nc_light,&
                        17.500E-18_dp ,4.00_dp ,2.50_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nC(x,y,z,2)*&
                        seaton(nug(frequency_counter),196622.4_dp*nc_light,&
                        4.600E-18_dp, 2.00_dp,1.80_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nC(x,y,z,3)*&
                        seaton(nug(frequency_counter),386241.0_dp*nc_light,&
                        1.430E-18_dp, 1.00_dp ,1.50_dp) 
!                   
! Nitrogen
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nN(x,y,z,1)*&
                        seaton(nug(frequency_counter),117314.6_dp*nc_light,&
                        12.720E-18_dp, 4.90_dp, 1.95_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nN(x,y,z,2)*&
                        seaton(nug(frequency_counter),238777.9_dp*nc_light,&
                        8.395E-18_dp, 2.00_dp, 2.15_dp )                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nN(x,y,z,3)*&
                        seaton(nug(frequency_counter),382587.5_dp*nc_light,&
                        2.160E-18_dp, 0.90e-18_dp, 2.00_dp)

 
! Oxygen
                   chi(x, y, z, frequency_counter)=&       
                        chi(x, y, z, frequency_counter)+&
                        nO(x,y,z,1)*&
                        seaton(nug(frequency_counter), 109758.7_dp*nc_light,&
                        2.940E-18_dp,2.66_dp, 1.00_dp)                    
                   chi(x, y, z, frequency_counter)=&      
                        chi(x, y, z, frequency_counter)+&
                        nO(x,y,z,1)*&
                        seaton(nug(frequency_counter), 136575.5_dp*nc_light,&
                        3.850E-18_dp, 4.38_dp, 1.50_dp)                   
                   chi(x, y, z, frequency_counter)=&       
                        chi(x,y, z, frequency_counter)+&
                        nO(x,y,z,1)*&
                        seaton(nug(frequency_counter), 229691.8_dp*nc_light,&
                        2.260E-18_dp, 4.31_dp, 1.50_dp)
                   chi(x, y, z, frequency_counter)=&       !OII 
                        chi(x, y, z, frequency_counter)+&
                        nO(x,y,z,2)*&
                        seaton(nug(frequency_counter), 283758.7_dp*nc_light,&
                        7E-18_dp,2.00_dp, 2.00_dp)                   
                   chi(x, y, z, frequency_counter)=&       !OIII
                        chi(x, y, z, frequency_counter)+&
                        nO(x,y,z,3)*&
                        seaton(nug(frequency_counter),443243.0_dp*nc_light,&
                        4.460E-18_dp, 2.00_dp, 2.30_dp )
!                   
! Neon
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nNe(x,y,z,1)*&
                        seaton(nug(frequency_counter),174192.4_dp*nc_light,&
                        5.350E-18_dp ,3.77_dp ,1.00_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nNe(x,y,z,2)*&
                        seaton(nug(frequency_counter),356930._dp*nc_light,&
                        8.790E-18_dp ,3.00_dp, 1.90_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nNe(x,y,z,3)*&
                        seaton(nug(frequency_counter),513828.3_dp*nc_light,&
                        1.800E-18_dp, 2.28_dp, 2.00_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nNe(x,y,z,3)*&
                        seaton(nug(frequency_counter),559820.7_dp*nc_light,&
                        2.500E-18_dp,2.35_dp ,2.50_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nNe(x,y,z,3)*&
                        seaton(nug(frequency_counter),579390.6_dp*nc_light,&
                        1.480E-18_dp,2.22_dp, 2.50_dp)
!                   
!Sulfur
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nS(x,y,z,1)*&
                        seaton(nug(frequency_counter),83363.7_dp*nc_light,&
                        3.140E-18_dp ,71.00_dp ,2.50_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nS(x,y,z,2)*&
                        seaton(nug(frequency_counter),189386.1_dp*nc_light,&
                        6.780E-18_dp ,2.00_dp,5.00_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nS(x,y,z,3)*&
                        seaton(nug(frequency_counter),282824.7_dp*nc_light,&
                        0.363E-18_dp, 2.00_dp, 2.00_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nS(x,y,z,3)*&
                        seaton(nug(frequency_counter),353934.5_dp*nc_light,&
                        0.364E-18_dp,2.00_dp ,2.00_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nS(x,y,z,3)*&
                        seaton(nug(frequency_counter),405699.7_dp*nc_light,&
                        0.090E-18_dp,4.00_dp,1.75_dp  )
!Argon
                  chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nAr(x,y,z,1)*&
                        seaton(nug(frequency_counter),127587.0_dp*nc_light,&
                        0.755E-18_dp ,2.00_dp ,2.00_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nAr(x,y,z,2)*&
                        seaton(nug(frequency_counter),222888.2_dp*nc_light,&
                        28.000E-18_dp ,9.50_dp,5.00_dp)                   
                   chi(x, y, z, frequency_counter)=&     
                        chi(x, y, z, frequency_counter)+&
                        nAr(x,y,z,3)*&
                        seaton(nug(frequency_counter),329420.5_dp*nc_light,&
                        1.000E-18_dp, 4.00_dp, 2.50_dp)                   
!
                end if
             end do
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine UPDATE_CHI



!================================================================================  
  subroutine UPDATE_ETA

!
!    Computes the emissivity of the gas. This subroutine is called if the
!    option "diffuse" is .true.
!    Warning: not complite. It currently only consideres the emissivity
!    caused by the recombination of hydrogen
! 
! called by the main program 
!
    use M_natural_constants,  only: nc_light, nc_pi, nc_lyman_edge, nc_parsec,&
                                    nc_planck, nc_electronmass, nc_boltzmann
    use M_definitions,        only: x_max, y_max, z_max, points, l_cell,&
                                    log_diffuse, use_mask
    use M_data_input,         only: nug
    use M_grid_memory,        only: eta, temperature, nHII, ne, mask
!
    implicit none
!
    integer(i4b):: x,y,z,frequency_counter
    real(dp)    :: bolometric_luminosity
    real(dp)    :: emitted_photons=0._dp
!================================================================================
!
    !$OMP PARALLEL PRIVATE(x,y,z, frequency_counter)
    !$OMP DO
    do z=-z_max, +z_max
       do y=-y_max,+y_max
          do x=-x_max, +x_max
             if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 0)&
                     cycle
             end if
             do frequency_counter=0,points-1          !HI
                if (nug(frequency_counter).lt.1.0973731569E+05*2.99792458E+10_dp)&
                     cycle 
                eta(x, y, z, frequency_counter)=&
                     (2*nc_planck*nug(frequency_counter)**3/nc_light**2)*&
                     ((nc_planck**2)/(2*nc_pi*nc_electronmass*nc_boltzmann*&
                     temperature(x,y,z)))**(3./2.)*&
                     seaton(nug(frequency_counter),nc_lyman_edge,&
                     6.3e-18_dp ,1.34_dp, 2.99_dp)*&
                     nHII(x,y,z)*ne(x,y,z)*&
                     exp(-(nc_planck*(nug(frequency_counter)-nc_lyman_edge))/&
                     (nc_boltzmann*temperature(x,y,z)))
             end do
         end do
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL
!

   if (log_diffuse) then
      bolometric_luminosity=0._dp
      emitted_photons=0._dp
      do z=1, z_max
         do y=1, y_max
            do x=1,x_max
               if (use_mask) then
                  if (mod(mask(x,y,z),2) .eq. 0)&
                       cycle
               end if
               do frequency_counter=1,points-2          !HI
                  bolometric_luminosity=bolometric_luminosity+&
                       (1./2.)*(&
                       eta(x,y,z,frequency_counter)+&
                       eta(x,y,z,frequency_counter+1))*&
                       abs(nug(frequency_counter+1)-(nug(frequency_counter)))*&
                       (l_cell*nc_parsec)**3*4*nc_pi
                  emitted_photons=emitted_photons+&
                       (1./2.)*(&
                       eta(x,y,z,frequency_counter)/&
                       (nc_planck*nug(frequency_counter))+&
                       eta(x,y,z,frequency_counter+1)/&
                       (nc_planck*nug(frequency_counter+1)))*&
                       abs(nug(frequency_counter+1)-&
                       (nug(frequency_counter)))*&
                       (l_cell*nc_parsec)**3*4*nc_pi                       
               end do
            end do
         end do
      end do
!
      print*, "bolometric diffuse emission in erg s^-1 /cell:",&
           bolometric_luminosity
      print*, "emitted 'diffuse photons' in s^-1", emitted_photons 
!
   end if
 end subroutine UPDATE_ETA

!=============================================================================   
 subroutine SOLVE_STATE(x,y,z,nion, nM, ion_rate_M, recomb_rate_M, heat_cool_M)
!
!   performs the math to solve the new state of the ionizatation and temperature
!   structure
!   currently ony for metals and with nions=4
! 
    use M_definitions,        only: t_step,&
                                    x_max, y_max, z_max,&
                                    compute_temperature,&
                                    use_mask
    use M_grid_memory,        only: energycontent, mask

   implicit none
!
   integer(i4b)                           :: nion
   integer(i4b)                           :: x,y,z
   integer(i4b)                           :: indexcounter ! to build matrix
   real(dp), dimension (nion, nion)       :: metalmatrix, metaleigenvectors,&
                                             metaleigenvectors_copy,&
                                             metalleftvectors
   real(dp), dimension (1:nion)           :: metalre, metalim, metal_start,&
                                             metalpivots
   real(dp), dimension (1:nion)           :: ion_rate_M,recomb_rate_M,&
                                             heat_cool_M
   real(dp), dimension (-x_max:+x_max,-y_max:+y_max,-z_max:+z_max, 0:nion)&
                                             :: nM
   real(dp), dimension (50)               :: work  
   integer(i4b),   dimension(1)           :: maxeig
   integer(i4b)                           :: info, counter  !for lapack
   logical                                :: temperature_masked
!=============================================================================
!
   metalmatrix=0
   if (nion .ne. 4) stop "nions not equal 4"
   do indexcounter=1,4                 ! diagonal elements
      metalmatrix(indexcounter, indexcounter) = -ion_rate_M(indexcounter)-&
                                                recomb_rate_M(indexcounter)
   end do
   do indexcounter=1, nion-1
      ! upper band
      metalmatrix(indexcounter, indexcounter+1) = recomb_rate_M(indexcounter+1)
      !lower band
      metalmatrix(indexcounter+1, indexcounter) = ion_rate_M(indexcounter) 
   end do

   metalmatrix(1,1)=-ion_rate_M(1)
   metalmatrix(1,2)=recomb_rate_M(2)
   metalmatrix(2,1)=ion_rate_M(1)
   metalmatrix(2,2)=-recomb_rate_M(2)-ion_rate_M(2)
   metalmatrix(2,3)=recomb_rate_M(3)
   metalmatrix(3,2)=ion_rate_M(2)
   metalmatrix(3,3)=-recomb_rate_M(3)-ion_rate_M(3)
   metalmatrix(3,4)=recomb_rate_M(4)
   metalmatrix(4,3)=ion_rate_M(3)
   metalmatrix(4,4)=-recomb_rate_M(4)
!
!   print*, "HELLO"
   call dgeev('N','V', nion,metalmatrix,nion,metalre, metalim,&
        metalleftvectors,nion,metaleigenvectors,nion,work, 50, info) 
   MAXEIG=maxloc(METALRE)
!       
   metal_start(1)=nM(x,y,z,1)
   metal_start(2)=nM(x,y,z,2)
   metal_start(3)=nM(x,y,z,3)
   metal_start(4)=nM(x,y,z,4)
!              
   metaleigenvectors_copy=metaleigenvectors
!   print*, "WORLD"
   call dgesv(nion,1,METALEIGENVECTORS_COPY,nion,METALPIVOTS,&
              METAL_START,nion,INFO)
!       
   nM(x,y,z,:)=0._dp
!       
   do counter=1,nion 
      !
      temperature_masked=.false.
      if (use_mask) then
         if (mask(x,y,z) .gt. 0)&
              temperature_masked=.true.  
      end if
!
      if (compute_temperature .and. .not. temperature_masked) then
         if (.true.) then !if (counter .eq. maxeig(1)) then
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(1)*metal_start(counter)*&
                 metaleigenvectors(1,counter)*t_step
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(2)*metal_start(counter)*&
                 metaleigenvectors(2,counter)*t_step
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(3)*metal_start(counter)*&
                 metaleigenvectors(3,counter)*t_step
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(4)*metal_start(counter)*&
                 metaleigenvectors(4,counter)*t_step                   
         else                
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(1)*metal_start(counter)*&
                 metaleigenvectors(1,counter)*&
                 (1./metalre(counter))*&
                 (exp(metalre(counter)*t_step)-1_dp)
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(2)*metal_start(counter)*&
                 metaleigenvectors(2,counter)*&
                 (1./metalre(counter))*&
                 (exp(metalre(counter)*t_step)-1._dp)
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(3)*metal_start(counter)*&
                 metaleigenvectors(3,counter)*&
                 (1./metalre(counter))*&
                 (exp(metalre(counter)*t_step)-1._dp)
            energycontent(x,y,z)=energycontent(x,y,z)+&
                 heat_cool_M(4)*metal_start(counter)*&
                 metaleigenvectors(4,counter)*&
                 (1./metalre(counter))*&
                 (exp(metalre(counter)*t_step)-1._dp)
         end if
      end if
!                
      nM(x,y,z,1)=nM(x,y,z,1)+metal_start(counter)*&
           metaleigenvectors(1,counter)*&
           exp(metalre(counter)*t_step)
      nM(x,y,z,2)=nM(x,y,z,2)+metal_start(counter)*&
           metaleigenvectors(2,counter)*&
           exp(metalre(counter)*t_step)
      nM(x,y,z,3)=nM(x,y,z,3)+metal_start(counter)*&
           metaleigenvectors(3,counter)*&
           exp(metalre(counter)*t_step)
      nM(x,y,z,4)=nM(x,y,z,4)+metal_start(counter)*&
           metaleigenvectors(4,counter)*&
           exp(metalre(counter)*t_step)
   end do
   nM(x,y,z,0)=nM(x,y,z,1)+nM(x,y,z,2)+nM(x,y,z,3)+nM(x,y,z,4)            
 end subroutine SOLVE_STATE
 
 
!================================================================================ 
  subroutine CALC_HYDROGEN
!
!    Computes the ionization structure of hydrogen.
!    As a side effect, thelength of the next timestep is computed.
!    This "dual use" should be changed in the future
!
!   called by the main program
!
    use M_natural_constants,  only: nc_light, nc_gyr, nc_boltzmann, nc_planck    
    use M_definitions,        only: t_step, t_step_max, x_max, y_max, z_max,&
                                    charge_transfer, &
                                    diffuse, filling_factor, points,&
                                    compute_temperature, t_all, worldage,&
                                    factor, initial_redshift, Omega_M, H_0,&
                                    redshift, use_mask
    use M_grid_memory,        only: ne, nHI, nHII, J_nu, nH_complete,& 
                                    temperature, temperature_old,&
                                    energycontent, mask
    use M_data_input,         only: nug
    use M_cosmology,          only: ZTOAGE, AGETOZ
!
  implicit none
!
    integer(i4b)  :: x, y, z
    real(dp)      :: ion_rate, xeq, x0, ne_old, col_ion_rate_HI, rec_rate
    real(dp)      :: t_step_old
    real(dp)      :: heat_cool_HI, heat_cool_HII
    logical       :: temperature_masked
!
    real   (dp), dimension (2,2)    :: HYDROGENMATRIX, HYDROGENEIGENVECTORS
    real   (dp), dimension (2,2)    :: HYDROGENMATRIX_COPY,&
                                       HYDROGENEIGENVECTORS_COPY
    real   (dp), dimension (2,2)    :: HYDROGENLEFTVECTORS
    real   (dp), dimension (2)      :: HYDROGENRE,HYDROGENIM
    real   (dp), dimension(2)       :: HYDROGEN_START, HYDROGEN_START_COPY
    real   (DP), dimension (2)      :: HYDROGEN_OCCUP_OLD
    integer (i4b), dimension(2)     :: HYDROGENPIVOTS
    integer (i4b), dimension(1)     :: MAXEIG
    real   (dp), dimension (50)     :: WORK 
    integer (i4b)                   :: INFO, counter
!================================================================================ 
!   
    write(*,*)
    write(*,*) "started calc hydrogen: determination of timestep"
    t_step_old=t_step
    t_step=t_step_max
!
    do z=-z_max,+z_max
       do y=-y_max,+y_max
          do x=-x_max,+x_max
!
             if (use_mask) then 
               if (mod(mask(x,y,z),2) .eq. 1 )&
                     cycle
             end if
!
             if (diffuse) then
                rec_rate=ne(x,y,z)*HYDREC(temperature(x,y,z),1)*&
                     nc_boltzmann*temperature(x,y,z)/filling_factor
             else             
                rec_rate=ne(x,y,z)*HYDREC(temperature(x,y,z),4)*nc_boltzmann*&
                     temperature(x,y,z)/filling_factor
             end if
             nHII_t= nHI(x,y,z)*&
                  photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*1.0973731569E+05_dp, 6.30E-18_dp,1.34_dp, 2.99_dp) &
                  -(nH_complete(x,y,z)-nHI(x,y,z))*rec_rate
             if (nHII_t .ne. 0) then ! .and. t_step .le. t_step_max) then
                if(t_step .gt. max(abs((factor*(ne(x,y,z)/nHII_t))),&
                     abs((factor*(nHI(x,y,z)/nHII_t))))) then
                   t_step = max(abs((factor*(ne(x,y,z)/nHII_t))),&
                        abs((factor*(nHI(x,y,z)/nHII_t))))
                   t_step=min(t_step,2*t_step_old)
                   x_co=x
                   y_co=y
                   z_co=z
                   nHI_co=nHI(x,y,z)
                end if
             end if
          end do
       end do
    end do
!
    print*, "max_delta_t 0",maxval(abs(temperature-temperature_old))
    print*, "t_step_no_temp:", t_step/(365.25_dp*24*3600)
!
    temperature_masked=.false.
    if (use_mask) then
       if (mask(x,y,z) .gt. 0)&
          temperature_masked=.true.  
    end if
!
    if (compute_temperature .and. .not. temperature_masked) then
       if (maxval(abs(temperature-temperature_old)) .lt. 1e-3) then
          t_step=t_step                             ! do nothing
       else 
          t_step=min(t_step,t_step*500/maxval(abs(temperature-temperature_old)))
      endif
   end if
   write(*,*) "t_step:",t_step/(365.25_dp*24*3600)
   write(*,*)
   if (t_step .lt.  24*3600) then         !minimal length of a timestep: one day
      t_step=24*3600
   end if   
   t_step=min(t_step,3*t_step_old) 

!       
    write(33,*) x_co, y_co, z_co, nHI_co
!
    print*, "BEFORE OPENMP"

    !print*,"determination of ionization structure of hydrogen"
    !$OMP PARALLEL PRIVATE(&
    !$OMP x,y,z, ne_old, rec_rate, ion_rate, col_ion_rate_HI, xeq, x0, &
    !$OMP heat_cool_HI, heat_cool_HII,temperature_masked,&
    !$OMP MAXEIG, HYDROGENMATRIX, HYDROGENEIGENVECTORS,HYDROGENMATRIX_COPY,&
    !$OMP HYDROGENEIGENVECTORS_COPY,HYDROGENLEFTVECTORS, HYDROGENRE,&
    !$OMP HYDROGENIM,HYDROGEN_START,HYDROGEN_START_COPY,HYDROGEN_OCCUP_OLD,&
    !$OMP HYDROGENPIVOTS, WORK, INFO, counter)
    !$OMP DO     
!
    do z=-z_max,+z_max
       do y=-y_max,+y_max
          do x=-x_max,+x_max
!
             if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1) then
                   print*, "cycle", x,y,z
                   cycle
                end if
             end if
 !
             ne_old=ne(x,y,z)
             ne_old=ne(x,y,z)
             ion_rate=photointegral(points-1, j_nu(x,y,z,1:), nug(1:), &
                  nc_light*1.0973731569E+05, 6.30E-18_dp,1.34_dp, 2.99_dp)
             if (diffuse) then
                rec_rate=ne(x,y,z)*HYDREC(temperature(x,y,z),1)/filling_factor
             else
                rec_rate=ne(x,y,z)*HYDREC(temperature(x,y,z),4)/filling_factor
             end if
             if (charge_transfer) &
                  call hydrogen_charge_transfer(x,y,z)
             col_ion_rate_HI=  ne(x,y,z)*1.55E+13/sqrt(temperature(x,y,z))*&
                  0.1_dp*6.256E-18_dp*&
                  exp(-(nc_planck*nc_light*109678.8)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*109678.8_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
             ion_rate=ion_rate+col_ion_rate_HI
!
!
             temperature_masked=.false.

             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
!                          
             if (compute_temperature .and. .not. temperature_masked) then  
                heat_cool_HI=0._dp
                heat_cool_HII=0._dp
!             print*, "x,y,z",x,y,z
                !            print*, "temperature:", temperature
                heat_cool_HI=photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                     nc_light*1.0973731569E+05_dp, 6.30E-18_dp,1.34_dp, 2.99_dp)
                heat_cool_HI=heat_cool_HI-&                       ! for 20000 K
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),1216._dp,&
                     1.04_dp,2._dp, 1E+9_dp)
                heat_cool_HI=heat_cool_HI-&                       ! for 20000 K
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),1026._dp,&
                     0.87_dp,2._dp, 1E+9_dp)
                if (diffuse) then 
                   heat_cool_HII=-ne(x,y,z)*HYDREC(temperature(x,y,z),2)*&
                   nc_boltzmann*temperature(x,y,z)
                else
                 heat_cool_HII=-ne(x,y,z)*HYDREC(temperature(x,y,z),5)*&
                     nc_boltzmann*temperature(x,y,z)
                end if
                
!!$                if (log_heating_cooling) then
!!$                   heat_H(x,y,z)=0._dp
!!$                   cool_H(x,y,z)=0._dp
!!$                   heat_H(x,y,z)=heat_H(x,y,z)+photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:), &
!!$                        nc_light*1.0973731569E+05, 6.30e-18_dp,1.34_dp, 2.99_dp)*nHI(x,y,z)
!!$                   cool_H(x,y,z)=cool_H(x,y,z)+&
!!$                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),1216._dp, 1.04_dp,2._dp, 1e+9_dp)*nHI(x,y,z)
!!$                   cool_H(x,y,z)=cool_H(x,y,z)+&
!!$                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),1026._dp, 0.87_dp,2._dp, 1e+9_dp)
!!$                   cool_H(x,y,z)=cool_H(x,y,z)+HYDREC(temperature(x,y,z),5)*&
!!$                        nc_boltzmann*temperature(x,y,z)*nHII(x,y,z)*ne(x,y,z)                        
!!$                end if
!                print*, "hxcdebug hydrogen"
!                print*, x,y,z,heat_cool_HI, heat_cool_hii 
             end if
             hydrogenmatrix=0._dp
!
             hydrogenmatrix(1,1)=-ion_rate
             hydrogenmatrix(1,2)= rec_rate
             hydrogenmatrix(2,1)= ion_rate
             hydrogenmatrix(2,2)= -rec_rate
!             
             call dgeev('N','V',2,hydrogenmatrix,2,hydrogenre,hydrogenim,&
                  hydrogenleftvectors,2,hydrogeneigenvectors,2,work,50,info)
             maxeig=maxloc(HYDROGENRE)
!
             hydrogen_start(1)=nHI(x,y,z)
             hydrogen_start(2)=nHII(x,y,z)
             hydrogen_start_copy=hydrogen_start
             hydrogeneigenvectors_copy=hydrogeneigenvectors
             call dgesv(2,1,hydrogeneigenvectors_copy,2,hydrogenpivots,hydrogen_start,2,info)
!             
             nHI(x,y,z)=0._dp
             nHII(x,y,z)=0._dp
!             
             do counter=1,2
                nHI(x,y,z)=nHI(x,y,z)+hydrogen_start(counter)*&
                     hydrogeneigenvectors(1,counter)*&
                     exp(hydrogenre(counter)*t_step)
                nHII(x,y,z)=nHII(x,y,z)+hydrogen_start(counter)*&
                     hydrogeneigenvectors(2,counter)*&
                     exp(hydrogenre(counter)*t_step)
             end do
!             
             if (compute_temperature .and. .not. temperature_masked) then
!
                if (.true.) then 
                   energycontent(x,y,z)=energycontent(x,y,z)+&
                        heat_cool_HI*hydrogen_start_copy(1)*t_step
                   energycontent(x,y,z)=energycontent(x,y,z)+&
                        heat_cool_HII*hydrogen_start_copy(2)*t_step
                else             !currently never called !            
                       energycontent(x,y,z)=energycontent(x,y,z)+&
                           heat_cool_HI*hydrogen_start(counter)*&
                           hydrogeneigenvectors(1,counter)*&
                           (1./hydrogenre(counter))*&
                           (exp(hydrogenre(counter)*t_step)-1._dp)
                      energycontent(x,y,z)=energycontent(x,y,z)+&
                           heat_cool_HII*hydrogen_start(counter)*&
                           hydrogeneigenvectors(2,counter)*&
                           (1./hydrogenre(counter))*&
                           (exp(hydrogenre(counter)*t_step)-1._dp)
                end if
             end if
          end do             
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL 
    print*, "end parallel"
!
    t_all=t_all+t_step
    worldage=t_all+ztoage(initial_redshift,Omega_M, H_0)
    redshift=agetoz(worldage/nc_gyr, Omega_M, H_0)

  contains

    subroutine hydrogen_charge_transfer(x,y,z)
      use M_grid_memory,      only: nHeII, nHeIII, nC, nN, nO, nNe, nS, nAr,&
                                    temperature 
      implicit none
      ! Helium
      integer(i4b)               :: x,y,z
      ion_rate = ion_rate + nHeII(x,y,z)/&
           filling_factor*chargetransfercoeff(&
           7.47e-6_dp,2.06_dp,9.93_dp,-3.89_dp,temperature(x,y,z))
      ion_rate = ion_rate +  nHeIII(x,y,z)/&
           filling_factor*chargetransfercoeff(&
           1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp, temperature(x,y,z))
      ! Carbon
      rec_rate = rec_rate + nC(x,y,z,1)/&
           filling_factor*chargetransfercoeff(&
           1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp,temperature(x,y,z))
      ion_rate = ion_rate + nC(x,y,z,2)/&
           filling_factor*chargetransfercoeff(&
           1.76e-9_dp,8.33_dp,4278.78_dp,-6.41_dp,temperature(x,y,z))
      ion_rate = ion_rate + nC(x,y,z,3)/&
           filling_factor*chargetransfercoeff(&
                    1.67e-4_dp,2.79_dp,304.72_dp,-4.07_dp, temperature(x,y,z))
      ion_rate = ion_rate + nC(x,y,z,4)/&
           filling_factor*chargetransfercoeff(&
           3.25e-4_dp,0.21_dp,0.19_dp,3.90_dp,temperature(x,y,z))
      ! Nitrogen
      rec_rate = rec_rate + nN(x,y,z,1)/&
           filling_factor*chargetransfercoeff(&
           4.55e-3_dp,-0.29_dp,-0.92_dp,-8.38_dp,temperature(x,y,z))
      ion_rate = ion_rate + nN(x,y,z,2)/&
           filling_factor*chargetransfercoeff(&
           1.01e-3_dp,-0.29_dp,-0.92_dp,-8.38_dp,temperature(x,y,z))
      ion_rate = ion_rate + nN(x,y,z,3)/&
           filling_factor*chargetransfercoeff(&
           3.05e-1_dp,0.60_dp,2.65_dp,-0.93_dp, temperature(x,y,z))
      ion_rate = ion_rate + nN(x,y,z,4)/&
           filling_factor*chargetransfercoeff(&
           4.54_dp,0.57_dp,-0.65_dp,-0.89_dp,temperature(x,y,z))  
      ! Oxygen
      rec_rate = rec_rate + nO(x,y,z,1)/&
           filling_factor*chargetransfercoeff(&
           7.40e-2_dp,0.47_dp,24.37_dp,-0.74_dp,temperature(x,y,z))
      ion_rate = ion_rate + nO(x,y,z,2)/&
           filling_factor*chargetransfercoeff(&
           1.04_dp,3.15e-2_dp,-0.61_dp,-9.73_dp,temperature(x,y,z))
      ion_rate = ion_rate + nO(x,y,z,3)/&
           filling_factor*chargetransfercoeff(&
           1.04_dp,0.27_dp,2.02_dp,-5.92_dp, temperature(x,y,z))
      ion_rate = ion_rate + nO(x,y,z,4)/&
           filling_factor*chargetransfercoeff(&
           3.98_dp,0.26_dp,0.56_dp,-2.62_dp,temperature(x,y,z))  
      ! Neon
      ion_rate = ion_rate + nNe(x,y,z,3)/&
           filling_factor*chargetransfercoeff(&
           1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp, temperature(x,y,z))
      ion_rate = ion_rate+ nNe(x,y,z,4)/&
           filling_factor*chargetransfercoeff(&
           14.73_dp,4.52e-2_dp,-0.84_dp,-0.31_dp,temperature(x,y,z))  
      ! Sulfur
      rec_rate = rec_rate + nS(x,y,z,1)/&
           filling_factor*chargetransfercoeff(&
           1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp,temperature(x,y,z))
      ion_rate = ion_rate + nS(x,y,z,2)/&
           filling_factor*chargetransfercoeff(&
           3.82e-7_dp,11.10_dp,2.57e+4_dp,-8.22_dp,temperature(x,y,z))
      ion_rate = ion_rate + nS(x,y,z,3)/&
           filling_factor*chargetransfercoeff(&
           1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp,temperature(x,y,z))
      ion_rate = ion_rate + nS(x,y,z,4)/&
           filling_factor*chargetransfercoeff(&
           2.29_dp,4.02e-2_dp,1.59_dp,-6.06_dp,temperature(x,y,z))  
      ! Argon
      ion_rate = ion_rate + nAr(x,y,z,3)/&
           filling_factor*chargetransfercoeff(&
           1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp, temperature(x,y,z))
      ion_rate = ion_rate + nAr(x,y,z,4)/&
           filling_factor*chargetransfercoeff(&
           4.57_dp,0.27_dp,-0.18_dp,-1.57_dp,temperature(x,y,z))  
    end subroutine hydrogen_charge_transfer
 
 end subroutine CALC_HYDROGEN
  
  

!================================================================================  
  subroutine CALC_HELIUM2
!
!    Computes the ionization structure of helium.
!    There is no longer a subroutine "CALC_HELIUM"
! 
! called by main program
!
    use M_natural_constants,  only: nc_light, nc_boltzmann, nc_planck
    use M_definitions,        only: x_max, y_max, z_max, points, filling_factor,&
                                    diffuse, compute_temperature,&
                                    log_heating_cooling, t_step,&
                                    use_mask, charge_transfer
    use M_grid_memory,        only: nHI, nHeI, nHeII, nHeIII,&
                                    ne, J_nu, temperature, energycontent, &
                                    heat_He, cool_He, mask
    use M_data_input,         only: nug
!
    implicit none
!    
    integer(i4b) :: x, y, z
    real(dp)     :: ion_rate_HeI, ion_rate_HeII, &
                    col_ion_rate_HeI, col_ion_rate_HeII,&
                    recomb_rate_HeII, recomb_rate_HeIII,&
                    ne_old, heat_cool_HeI, heat_cool_HeII, heat_cool_HeIII
!
    real   (dp), dimension (3,3)   :: HELIUMMATRIX, HELIUMEIGENVECTORS
    real   (dp), dimension (3,3)   :: HELIUMMATRIX_COPY, HELIUMEIGENVECTORS_COPY
    real   (dp), dimension (3,3)   :: HELIUMLEFTVECTORS
    real   (dp), dimension   (3)   :: HELIUMRE,HELIUMIM,HELIUM_START,&
                                      HELIUM_START_OLD
    real   (dp), dimension   (3)   :: HELIUM_OCCUP_OLD
    integer (i4b), dimension (3)   :: HELIUMPIVOTS
    integer (i4b), dimension (1)   :: MAXEIG
    real   (dp), dimension  (50)   :: WORK 
    integer (i4b)                  :: INFO, counter
    logical                        :: temperature_masked
!===============================================================================
!  
  print*, "calc helium"
!    
    !$OMP PARALLEL private (y,x, ion_rate_HeI, col_ion_rate_HeI, ion_rate_HeII,&
    !$OMP col_ion_rate_HeII,&
    !$OMP recomb_rate_HeII, recomb_rate_HeIII, maxeig, HELIUMMATRIX,&
    !$OMP HELIUMEIGENVECTORS,&
    !$OMP HELIUMMATRIX_COPY, HELIUMEIGENVECTORS_COPY, HELIUMLEFTVECTORS,&
    !$OMP HELIUMRE,HELIUMIM,HELIUM_START,HELIUM_START_OLD,HELIUM_OCCUP_OLD,&
    !$OMP HELIUMPIVOTS,WORK,INFO, counter,&
    !$OMP heat_cool_HeI, heat_cool_HeII, heat_cool_HeIII, temperature_masked)  
    !$OMP DO 
    do z=-z_max,+z_max
       do y=-y_max,+y_max
          do x=-x_max,+x_max
!
             if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1)&
                     cycle
             end if
            ne_old=ne(x,y,z)
!
             ion_rate_HeI=photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*198305.0_dp , 7.466E-18_dp,1.20_dp, 2.00_dp)    
             col_ion_rate_HeI=  ne(x,y,z)*1.55E+13/sqrt(temperature(x,y,z))*0.1*&
                  7.466E-18*&
                  exp(-(nc_planck*nc_light*198305.0_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*198305.0_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
             ion_rate_HeI=ion_rate_HeI+col_ion_rate_HeI             
             ion_rate_HeII=photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*438908.7_dp , 1.573e-18_dp,1.30_dp, 3.00_dp)
             col_ion_rate_HeII= ne(x,y,z)*1.55e+13_dp/sqrt(temperature(x,y,z))*&
                  0.2_dp*1.573e-18_dp*&
                  exp(-(nc_planck*nc_light*438908.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*438908.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
             ion_rate_HeII=ion_rate_HeII+col_ion_rate_HeII
!
             if (diffuse) then
                recomb_rate_HeII=HELREC(temperature(x,y,z),1)*ne(x,y,z)/&
                     filling_factor             
                recomb_rate_HeIII=2*HYDREC(temperature(x,y,z)/4,1)*ne(x,y,z)/&
                     filling_factor    
             else             
                recomb_rate_HeII=HELREC(temperature(x,y,z),4)*ne(x,y,z)/&
                     filling_factor    
                recomb_rate_HeIII=2*HYDREC(temperature(x,y,z)/4,4)*ne(x,y,z)/&
                     filling_factor    
             end if
!
             if(charge_transfer) then
                recomb_rate_HeII = recomb_rate_HeII+ nHI(x,y,z)/&
                     filling_factor*chargetransfercoeff(&
                     7.47e-6_dp,2.06_dp,9.93_dp,-3.89_dp,temperature(x,y,z))
                recomb_rate_HeIII = recomb_rate_HeIII+ nHI(x,y,z)/&
                     filling_factor*chargetransfercoeff(&
                     1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp, temperature(x,y,z))
             endif
!                          
             temperature_masked=.false.
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
!
             if (compute_temperature .and. .not. temperature_masked) then
                heat_cool_HeI=&
                     photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                     nc_light*198305.0_dp, 7.466E-18_dp,1.20_dp, 2.00_dp)
                heat_cool_HeII=&
                     photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                     nc_light*438908.7_dp , 1.573E-18_dp,1.30_dp, 3.00_dp)
                heat_cool_HeII=heat_cool_HeII-&
                     HELREC(temperature(x,y,z),5)*ne(x,y,z)*&
                     nc_boltzmann*temperature(x,y,z)
                if (diffuse) then
                      heat_cool_HeIII=-2*HYDREC(temperature(x,y,z)/4.0_dp,2)*&
                           ne(x,y,z)*nc_boltzmann*temperature(x,y,z)
                else
                heat_cool_HeIII=-2*HYDREC(temperature(x,y,z)/4.0_dp,5)*&
                     ne(x,y,z)*nc_boltzmann*temperature(x,y,z)
                end if
!
               if (log_heating_cooling) then
                  heat_He(x,y,z)=0._dp
                  cool_He(x,y,z)=0._dp
                  heat_He(x,y,z)= heat_He(x,y,z)+&
                       photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                       nc_light*198305.0_dp, 7.466E-18_dp,1.20_dp, 2.00_dp)*&
                       nHeI(x,y,z)
                  heat_He(x,y,z)=heat_He(x,y,z)+&
                       photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                       nc_light*438908.7_dp , 1.573E-18_dp,1.30_dp, 3.00_dp)*&
                       nHeII(x,y,z)
                  if (diffuse) then
                     cool_He(x,y,z)=cool_He(x,y,z)+&
                          HELREC(temperature(x,y,z),2)*ne(x,y,z)*&
                          nc_boltzmann*temperature(x,y,z)*nHeII(x,y,z)
                  else
                     cool_He(x,y,z)=cool_He(x,y,z)+&
                          HELREC(temperature(x,y,z),5)*ne(x,y,z)*&
                          nc_boltzmann*temperature(x,y,z)*nHeII(x,y,z)
                  end if
                  if (diffuse) then
                     cool_He(x,y,z)=cool_He(x,y,z)+&
                          2*HYDREC(temperature(x,y,z)/4.0_dp,2)*ne(x,y,z)*&
                          nc_boltzmann*temperature(x,y,z)*nHeIII(x,y,z)
                  else
                     cool_He(x,y,z)=cool_He(x,y,z)+&
                          2*HYDREC(temperature(x,y,z)/4.0_dp,5)*ne(x,y,z)*&
                          nc_boltzmann*temperature(x,y,z)*nHeIII(x,y,z)
                  end if
               end if
            end if
!                                     
            heliummatrix=0._dp
            heliummatrix(1,1)=-ion_rate_HeI
            heliummatrix(1,2)=recomb_rate_HeII
            heliummatrix(2,1)=ion_rate_HeI
            heliummatrix(2,2)=-recomb_rate_HeII-ion_rate_HeII
            heliummatrix(2,3)=recomb_rate_HeIII
            heliummatrix(3,2)=ion_rate_HeII
            heliummatrix(3,3)=-recomb_rate_HeIII
            !                         
            call dgeev('N','V', 3,heliummatrix,3,heliumre, heliumim,&
                 heliumleftvectors, 3, heliumeigenvectors,3,work, 50, info) 
            MAXEIG=maxloc(HELIUMRE)
            helium_start(1)=nHeI(x,y,z)
            helium_start(2)=nHeII(x,y,z)
            helium_start(3)=nHeIII(x,y,z)
!             
            helium_start_old(1)=helium_start(1)
            helium_start_old(2)=helium_start(2)
            helium_start_old(3)=helium_start(3)
!                          
            heliumeigenvectors_copy=heliumeigenvectors
            call dgesv(3,1,HELIUMEIGENVECTORS_COPY,3,HELIUMPIVOTS,&
                 HELIUM_START,3,INFO)
!             
            nHeI(x,y,z)=0._dp
            nHeII(x,y,z)=0._dp
            nHeIII(x,y,z)=0._dp
!             
            do counter=1,3
               nHeI(x,y,z)=nHeI(x,y,z)+helium_start(counter)*&
                    heliumeigenvectors(1,counter)*&
                    exp(heliumre(counter)*t_step)
               nHeII(x,y,z)=nHeII(x,y,z)+helium_start(counter)*&
                    heliumeigenvectors(2,counter)*&
                    exp(heliumre(counter)*t_step)
               nHeIII(x,y,z)=nHeIII(x,y,z)+helium_start(counter)*&
                    heliumeigenvectors(3,counter)*&
                    exp(heliumre(counter)*t_step)
            end do
!
            if (compute_temperature) then
               if (.true.) then
                  energycontent(x,y,z)=energycontent(x,y,z)+&
                       heat_cool_HeI*helium_start_old(1)*t_step
                  energycontent(x,y,z)=energycontent(x,y,z)+&
                       heat_cool_HeII*helium_start_old(2)*t_step
                  energycontent(x,y,z)=energycontent(x,y,z)+&
                       heat_cool_HeIII*helium_start_old(3)*t_step
                !else currently never called 
               else      
                  energycontent(x,y,z)=energycontent(x,y,z)+&
                       heat_cool_HeI*helium_start(counter)* &
                       heliumeigenvectors(1,counter)*&
                       (1./heliumre(counter))*exp(heliumre(counter)*t_step)
                  energycontent(x,y,z)=energycontent(x,y,z)+&
                       heat_cool_HeII*helium_start(counter)*&
                       heliumeigenvectors(2,counter)*&
                       (1./heliumre(counter))*exp(heliumre(counter)*t_step)
                  energycontent(x,y,z)=energycontent(x,y,z)+&
                       heat_cool_HeIII*helium_start(counter)*&
                       heliumeigenvectors(3,counter)*&
                       (1./heliumre(counter))*exp(heliumre(counter)*t_step)
               end if
             end if
!             
             if (nHeI(x,y,z) .lt. 0._dp) then     ! this should never happen
                print*, "Particle conservation did not work properly !"
                print*,x,y,z
                print*,nHeI(x,y,z) 
                stop
                nHeI(x,y,z) = 0._dp
             end if
 !            
          end do
       end do
    end do    
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine CALC_HELIUM2
  


!================================================================================
  subroutine CALC_METALS
!  
!     Computes the ionization structures of C, N, O, Ne ans S. 
!     Considered are the ionization stages I to IV   
!
    use M_natural_constants, only: nc_light, nc_boltzmann, nc_planck
    use M_definitions,       only: x_max, y_max, z_max,&
                                  charge_transfer,&
                                   compute_temperature, photoheat_metals,&
                                   filling_factor,points,&
                                   log_heating_cooling, use_mask 
    use M_data_input,        only: nug
    use M_grid_memory,       only: ne, nHI, nHII,  nC, nN, nO, nNe, nS, nAr,&
                                   temperature, j_nu,&
                                   cool_C, cool_N, cool_O, cool_Ne,&
                                   cool_S, cool_Ar, mask
  implicit none
!
  real(dp):: ne_old
!===============================================================================   
!    
    print*, "calc_metals"
!    
    write(*,*) "calc carbon"
    call calc_carbon()
    write(*,*) "calc nitrogen"
    call calc_nitrogen()
    write(*,*) "calc oxygen"
    call calc_oxygen()
    write(*,*) "calc neon"
    call calc_neon()
    write(*,*) "calc sulfur"
    call calc_sulfur()
    write(*,*) "calc argon"
    call calc_argon()
!-------------------------------------------------------------------------------

  contains

!==============================================================================


!================================================================================
     subroutine calc_carbon()
!
!      computes ionization structure and of carbon.
!     
! 
     implicit none
!
    integer(i4b), parameter                  :: nion=4
    integer(i4b)                             :: x,y,z,iteration_counter
!
    real(dp), dimension(1:nion) :: ion_rate_C=0,&
                                   col_ion_rate_C=0,&
                                   recomb_rate_C=0,&
                                   heat_cool_C=0
    logical                     :: temperature_masked
!        
!================================================================================
!
!$OMP PARALLEL PRIVATE(x, y, z, iteration_counter,&
!$OMP   ion_rate_C,&
!$OMP   col_ion_rate_C,&
!$OMP   recomb_rate_C,&
!$OMP   heat_cool_C, temperature_masked) 
!$OMP DO 
    do z=-z_max, +z_max
       do y= -y_max,+y_max
          do x=-x_max, x_max
!
            if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1)&
                     cycle
             end if
!
             dull_dummy=1+(x+x_max)+21*(y+y_max)+21*21*(z+z_max)
             ne_old=ne(x,y,z)           
! collisional_ionization rates
             col_ion_rate_C(1)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.1_dp*12.20E-18_dp*&
                  exp(-(nc_planck*nc_light*9.09E+4_dp)/&
                     (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*9.09E+4_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_C(2)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.2_dp*4.60E-18_dp*&
                  exp(-(nc_planck*nc_light*1.97E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*1.97E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_C(3)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.3*3.86e-18*&
                  exp(-(nc_planck*nc_light*3.86E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*3.86E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
!ion_rate= collisional ionization rates + photoionization rates
             ion_rate_C(1) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*9.09E+4_dp, &
                  12.20E-18_dp, 3.32_dp, 2.00_dp)
             ion_rate_C(1)=ion_rate_C(1)+col_ion_rate_C(1)
!       
             ion_rate_C(2) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*1.97e+5, &
                  4.60E-18_dp, 1.95_dp, 3.00_dp)
             ion_rate_C(2)=ion_rate_C(2)+col_ion_rate_C(2)
!       
             ion_rate_C(3) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*3.86e+5, &
                  1.60E-18_dp, 2.60_dp, 3.00_dp)
             ion_rate_C(3)=ion_rate_C(3)+col_ion_rate_C(3)
!       
             recomb_rate_C(2)= radiative_rec3(4.70E-13_dp,0.624_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_C(2)=recomb_rate_C(2)+&
                  dielectronic_hot(2.54E-3_dp,4.42E-2_dp,1.57E+5_dp,3.74E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z) 
             recomb_rate_C(2)=recomb_rate_C(2)+& 
                  dielectronic_cool(0.0108_dp,-0.1075_dp,0.2810_dp,-0.0193_dp,&
                  -0.1127_dp,temperature(x,y,z))*ne(x,y,z)
             recomb_rate_C(2)=recomb_rate_C(2)/filling_factor
!       
             recomb_rate_C(3)=radiative_rec3(2.30E-12_dp,0.645_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_C(3)=recomb_rate_C(3)+&
                  dielectronic_hot(6.15e-3_dp,5.88E-2_dp,1.41E+5_dp,1.41E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_C(3)=recomb_rate_C(3)+& 
                  dielectronic_cool(1.8267_dp,4.1012_dp,4.8443_dp,0.2261_dp,&
                  0.5960_dp,temperature(x,y,z))*ne(x,y,z)
             recomb_rate_C(3)=recomb_rate_C(3)/filling_factor
!       
             recomb_rate_C(4)=radiative_rec3(4.90E-12_dp,0.803_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_C(4)=recomb_rate_C(4)+&
                  dielectronic_hot(1.62E-3_dp,3.43E-1_dp,8.19E+4_dp,1.59E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_C(4)=recomb_rate_C(4)+& 
                  dielectronic_cool(2.3196_dp,10.7328_dp,6.8830_dp,-0.1824_dp,&
                  0.4101_dp,temperature(x,y,z))*ne(x,y,z)      
             recomb_rate_C(4)=recomb_rate_C(4)/filling_factor
!
             if (charge_transfer) then
                ion_rate_C(1) = ion_rate_C(1) + nHII(x,y,z)/&
                     filling_factor*chargetransfercoeff(&
                     1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp,temperature(x,y,z))
                recomb_rate_C(2) = recomb_rate_C(2)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.76e-9_dp,8.33_dp,4278.78_dp,-6.41_dp,temperature(x,y,z))
                recomb_rate_C(3) = recomb_rate_C(3)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.67e-4_dp,2.79_dp,304.72_dp,-4.07_dp, temperature(x,y,z))
                recomb_rate_C(4) = recomb_rate_C(4)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    3.25_dp,0.21_dp,0.19_dp,3.90_dp,temperature(x,y,z))  
             end if
!
             temperature_masked=.false.
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
!
             if (compute_temperature .and. .not. temperature_masked) then
                heat_cool_C=0._dp
!              
                if (photoheat_metals) then
                   heat_cool_C(1)=heat_cool_C(1)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*9.09E+4_dp,& 
                        12.20E-18_dp, 3.32_dp, 2.00_dp)/filling_factor
                   heat_cool_C(2)=heat_cool_C(2)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*1.97e+5, &
                        4.60E-18_dp, 1.95_dp, 3.00_dp)/filling_factor
                   heat_cool_C(3)=heat_cool_C(3)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*3.86e+5, &
                        1.60E-18_dp, 2.60_dp, 3.00_dp)/filling_factor 
                end if
!
                if (log_heating_cooling)  then
                   cool_C(x,y,z)=0._dp
                end if
 
             end if
             call solve_state(x,y,z,nion, nC, ion_rate_C, recomb_rate_C,&
                              heat_cool_C)
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine calc_carbon

!================================================================================ 
  subroutine calc_nitrogen()
!
!              Computes ionization structure of nitrogen and its influence
!              on the temperature structure of the gas
!
    implicit none
!
    integer(i4b), parameter :: nion=4
    integer(i4b):: x, y, z, iteration_counter
!
    real(dp), dimension(1:nion) :: ion_rate_N=0,&
                                   col_ion_rate_N=0,&
                                   recomb_rate_N=0,&
                                   heat_cool_N=0
    logical                     :: temperature_masked
!        
!================================================================================
!$OMP PARALLEL PRIVATE(x, y, z, iteration_counter,&
!$OMP   ion_rate_N, col_ion_rate_N,&
!$OMP   recomb_rate_N, heat_cool_N, temperature_masked) 
!$OMP DO 
    do z=-z_max, +z_max
       do y= -y_max,+y_max
          do x=-x_max, x_max
!
            if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1)&
                     cycle
             end if
!
             ne_old=ne(x,y,z)
!
             col_ion_rate_N(1)=  ne(x,y,z)*1.55E+13/sqrt(temperature(x,y,z))*&
                  0.1_dp*12.720E-18_dp*&
                  exp(-(nc_planck*nc_light*117314.6_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*117314.6_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
             
             col_ion_rate_N(2)=  ne(x,y,z)*1.55E+13_dp/sqrt(temperature(x,y,z))*&
                  0.2_dp*8.395e-18_dp*&
                  exp(-(nc_planck*nc_light*238777.9_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*238777.9_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
             
             col_ion_rate_N(3)=  ne(x,y,z)*1.55E+13_dp/sqrt(temperature(x,y,z))*&
                  0.3_dp*2.160E-18_dp*&
                  exp(-(nc_planck*nc_light*382587.5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*382587.5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!             
             ion_rate_N(1) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*117314.6, 11.40E-18_dp, 4.29_dp, 2.00_dp)
             ion_rate_N(1)=ion_rate_N(1)+col_ion_rate_N(1)
!                          
             ion_rate_N(2) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*238777.9, &
                  6.65E-18_dp, 2.86_dp, 3.00_dp)
             ion_rate_N(2)=ion_rate_N(2)+col_ion_rate_N(2)
!             
             ion_rate_N(3) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*382587.5, &
                  2.060E-18_dp, 0.90_dp, 2.00_dp)
             ion_rate_N(3)=ion_rate_N(3)+col_ion_rate_N(3)
!             
             recomb_rate_N(2)=radiative_rec3(4.10E-13_dp,0.608_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_N(2)= recomb_rate_N(2)+&
                  dielectronic_hot(2.98E-3_dp,0._dp,2.20E+5_dp,1E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z) 
             recomb_rate_N(2)=recomb_rate_N(2)/filling_factor
!             
             recomb_rate_N(3)=radiative_rec3(2.20E-12_dp,0.639_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_N(3)= recomb_rate_N(3)+&
                  dielectronic_hot(7.41E-3_dp, 7.64E-2_dp, 2.01E+5_dp,&
                  7.37E+4_dp,temperature(x,y,z)) * ne(x,y,z)
             recomb_rate_N(3)=recomb_rate_N(3)+&
                  dielectronic_cool(0.0320_dp,-0.6624_dp,4.3191_dp,0.0003_dp,&
                  0.5946_dp,temperature(x,y,z))*ne(x,y,z)
             recomb_rate_N(3)=recomb_rate_N(3)/filling_factor
!             
             recomb_rate_N(4)=radiative_rec3(5.00E-12_dp,0.676_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_N(4)=recomb_rate_N(4)+&
                  dielectronic_hot(1.13E-2_dp,1.64E-1_dp,1.73E+5_dp,2.25E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_N(4)=recomb_rate_N(4)+&
                  dielectronic_cool(-0.8806_dp,11.2406_dp,30.7066_dp,-1.1721_dp,&
                  0.6127_dp,temperature(x,y,z))*ne(x,y,z)
             recomb_rate_N(4)=recomb_rate_N(4)/filling_factor
! 
             if (charge_transfer) then
                ion_rate_N(1) = ion_rate_N(1) + nHII(x,y,z)/&
                     filling_factor*chargetransfercoeff(&
                     4.55e-3_dp,-0.29_dp,-0.92_dp,-8.38_dp,temperature(x,y,z))
                recomb_rate_N(2) = recomb_rate_N(2)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.01e-3_dp,-0.29_dp,-0.92_dp,-8.38_dp,temperature(x,y,z))
                recomb_rate_N(3) = recomb_rate_N(3)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    3.05e-1_dp,0.60_dp,2.65_dp,-0.93_dp, temperature(x,y,z))
                recomb_rate_N(4) = recomb_rate_N(4)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    4.54_dp,0.57_dp,-0.65_dp,-0.89_dp,temperature(x,y,z))  
             end if                         
!
             temperature_masked=.false.
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
             if (compute_temperature .and. .not. temperature_masked) then
                heat_cool_N=0._dp
                if (photoheat_metals) then
                   heat_cool_N(1)=heat_cool_N(1)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*117314.6,&
                        11.40E-18_dp, 4.29_dp, 2.00_dp)
                   heat_cool_N(2)=heat_cool_N(2)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*238777.9, &
                        6.65E-18_dp, 2.86_dp, 3.00_dp)
                   heat_cool_N(3)=heat_cool_N(3)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*382587.5, &
                        2.060E-18_dp, 0.90_dp, 2.00_dp)                     
                end if
!
                heat_cool_N(2)=heat_cool_N(2)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     6570._dp, 2.64_dp,9._dp, 6.6E+4_dp)/filling_factor
                heat_cool_N(2)=heat_cool_N(2)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     2055000._dp, 0.41_dp,1._dp, 8E+1_dp)/filling_factor
                heat_cool_N(2)=heat_cool_N(2)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     765000._dp, 0.27_dp,1._dp, 3.1E+2_dp)/filling_factor
                heat_cool_N(3)=-heat_cool_N(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     573430._dp, 1.45_dp,2._dp, 1.5E+3_dp)/filling_factor
                heat_cool_N(4)=0._dp
                !             
                if (log_heating_cooling)  then
                   cool_N(x,y,z)=0._dp
                   cool_N(x,y,z)=cool_N(x,y,z)+nN(x,y,z,2)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        6570._dp, 2.64_dp, 9._dp, 6.6E+4_dp)/filling_factor
                   cool_N(x,y,z)=cool_N(x,y,z)*nN(x,y,z,2)*&
                        (-col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        6570._dp,2.64_dp,9._dp, 6.6E+4_dp))/filling_factor  
                   cool_N(x,y,z)=cool_N(x,y,z)+nN(x,y,z,2)*&
                        (col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        765000._dp, 0.27_dp,1._dp, 3.1E+2_dp))/filling_factor  
                   cool_N(x,y,z)=cool_N(x,y,z)+nN(x,y,z,3)*&
                        (col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        573430._dp, 1.45_dp,2._dp, 1.5E+3_dp))/filling_factor
                end if
             end if
!
             call solve_state(x,y,z,nion, nN, ion_rate_N, recomb_rate_N,&
                              heat_cool_N)                          
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine calc_nitrogen

!===============================================================================  
  subroutine calc_oxygen()
!
!              Computes ionization structure of nitrogen and its influence
!              on the temperature structure of the gas.
!
             implicit none
!
             integer(i4b), parameter             :: nion=4
!
             integer(i4b)                        :: x, y, z, iteration_counter
!
             real(dp),      dimension (1:nion)   :: ion_rate_O=0,&
                                                    col_ion_rate_O=0,&
                                                    recomb_rate_O=0,&
                                                    heat_cool_O=0
             logical                             :: temperature_masked
!        
!===============================================================================
!
!$OMP PARALLEL PRIVATE(x, y, z, iteration_counter,&
!$OMP   ion_rate_O, col_ion_rate_O,&
!$OMP   recomb_rate_O, heat_cool_O, temperature_masked)
!$OMP DO 
    do z=-z_max, +z_max
       do y= -y_max,+y_max
          do x=-x_max, x_max
!
             if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1)&
                     cycle
             end if
!
             ne_old=ne(x,y,z)
             col_ion_rate_O(1)=  ne(x,y,z)*1.55E+13/&
                  sqrt(temperature(x,y,z))*0.1_dp*2.940E-18_dp*&
                  exp(-(nc_planck*nc_light*109758.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*109758.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor  
             !             
             col_ion_rate_O(2)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*0.2*10.400E-18_dp*&
                  exp(-(nc_planck*nc_light*283758.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*283758.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!             
             col_ion_rate_O(3)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*0.3*4.460E-18*&
                  exp(-(nc_planck*nc_light*443243.0)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*443243.0)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!                          
             ion_rate_O(1) = photointegral(points-1, j_nu(x,y,z,1:),&
                  nug(1:), nc_light*109758.7 , &
                  2.940E-18_dp,2.66_dp, 1.00_dp)
             ion_rate_O(1)=ion_rate_O(1)+col_ion_rate_O(1)
!                          
             ion_rate_O(1) = ion_rate_O(1)+photointegral(points-1, j_nu(x,y,z,1:),&
                  nug(1:), nc_light*136300_dp , &
                  3.85E-18_dp,4.38_dp, 1.50_dp)
             ion_rate_O(1)=ion_rate_O(1)+col_ion_rate_O(1)
!             
             ion_rate_O(1) = ion_rate_O(1)+photointegral(points-1, j_nu(x,y,z,1:),&
                  nug(1:), nc_light*1.500E+5_dp , &
                  2.26E-18_dp,4.31_dp, 1.50_dp)
             ion_rate_O(1)=ion_rate_O(1)+col_ion_rate_O(1)
!                         
             ion_rate_O(2)= photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*283758.7_dp , &
                  7.32E-18_dp,3.84_dp, 2.50_dp)
             ion_rate_O(2)=ion_rate_O(2)+col_ion_rate_O(2)
!             
             ion_rate_O(3)= photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*443243.0_dp, &
                  4.320E-18_dp,2.01_dp, 3.00_dp)
             ion_rate_O(3)=ion_rate_O(3)+col_ion_rate_O(3)
!             
             recomb_rate_O(2)=radiative_rec3(3.10E-13_dp,0.678_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_O(2)=recomb_rate_O(2)+&
                  dielectronic_hot(1.11E-3_dp,9.25E-2_dp,1.75E+5_dp,1.45E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z) 
             if (temperature(x,y,z) .lt. 2.0E+4_dp) then
                recomb_rate_O(2)=recomb_rate_O(2)+&
                     dielectronic_cool(-0.0001_dp,0.0001_dp,0.0956_dp,0.0193_dp,&
                     0.4106_dp,temperature(x,y,z))*ne(x,y,z)
             else
                recomb_rate_O(2)=recomb_rate_O(2)+&
                     dielectronic_cool(0.3715_dp,-0.0239_dp,-0.0597_dp,0.0678_dp,&
                     0.7993_dp,temperature(x,y,z))*ne(x,y,z)
             end if
 
             recomb_rate_O(2)=recomb_rate_O(2)/filling_factor             
             recomb_rate_O(3)=radiative_rec3(2.00E-12_dp,0.646_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_O(3)=recomb_rate_O(3)+&
                  dielectronic_hot(5.07E-3_dp,1.81E-1_dp,1.98E+5_dp,3.35E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_O(3)=recomb_rate_O(3)+&
                  dielectronic_cool(-0.0036_dp,0.7519_dp,1.5252_dp,-0.0838_dp,&
                  0.2769_dp,temperature(x,y,z))*ne(x,y,z)             
             recomb_rate_O(3)=recomb_rate_O(3)/filling_factor             
             recomb_rate_O(4)=radiative_rec3(5.10E-12_dp,0.66_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_O(4)= recomb_rate_O(4)+&
                  dielectronic_hot(1.48E-2_dp,3.05E-1_dp,2.41E+5_dp,2.83E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_O(4)=recomb_rate_O(4)+&
                  dielectronic_cool(0._dp,21.879_dp,16.273_dp,-0.702_dp,&
                  1.1899_dp,temperature(x,y,z))*ne(x,y,z)
             recomb_rate_O(4)=recomb_rate_O(4)/filling_factor
!
             if (charge_transfer) then
                ion_rate_O(1) = ion_rate_O(1) + nHII(x,y,z)/&
                     filling_factor*chargetransfercoeff(&
                     7.40e-2_dp,0.47_dp,24.37_dp,-0.74_dp,temperature(x,y,z))
                recomb_rate_O(2) = recomb_rate_O(2)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.04_dp,3.15e-2_dp,-0.61_dp,-9.73_dp,temperature(x,y,z))
                recomb_rate_O(3) = recomb_rate_O(3)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.04_dp,0.27_dp,2.02_dp,-5.92_dp, temperature(x,y,z))
                recomb_rate_O(4) = recomb_rate_O(4)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    3.98_dp,0.26_dp,0.56_dp,-2.62_dp,temperature(x,y,z))  
             end if
!
             temperature_masked=.false.
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
             if (compute_temperature .and. .not. temperature_masked) then
!                          
                heat_cool_O=0._dp
!             
                if (photoheat_metals) then
                   heat_cool_O(1)=heat_cool_O(1)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:),&
                        nug(1:), nc_light*1.500E+5_dp , &
                        2.26E-18_dp,4.31_dp, 1.50_dp)
                   heat_cool_O(2)=heat_cool_O(2)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*283758.7_dp , &
                        7.32E-18_dp,3.84_dp, 2.50_dp)
                   heat_cool_O(3)=heat_cool_O(3)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*443243.0_dp, &
                        4.320E-18_dp,2.01_dp, 3.00_dp)                     
                end if
!
                heat_cool_O(2)=heat_cool_O(2)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     3728.8_dp, 1.34_dp,4._dp, 3.4E+3_dp)/&
                     filling_factor
                heat_cool_O(3)=heat_cool_O(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     5000._dp, 2.29_dp,9._dp, 6.8E+5_dp)/&
                     filling_factor
                heat_cool_O(3)=heat_cool_O(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),883560._dp,&
                     0.55_dp,1._dp, 5.1E+2_dp)/&
                     filling_factor
                heat_cool_O(3)=heat_cool_O(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),326610._dp,&
                     0.27_dp,1._dp, 3.6E+3_dp)/&
                     filling_factor
                heat_cool_O(4)=heat_cool_O(4)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     259130._dp, 2.34_dp,2._dp, 1.5E+3_dp)/&
                     filling_factor
                !
                if (log_heating_cooling) then
                   cool_O(x,y,z)=0._dp
                   cool_O(x,y,z)=cool_O(x,y,z)+nO(x,y,z,2)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        3728.8_dp, 1.34_dp,4._dp, 3.4E+3_dp)/filling_factor
                   cool_O(x,y,z)=cool_O(x,y,z)+nO(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        5000._dp, 2.29_dp,9._dp, 6.8E+5_dp)/filling_factor
                   !
                   cool_O(x,y,z)=cool_O(x,y,z)+nO(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        883560._dp, 0.55_dp,1._dp, 5.1E+2_dp)/filling_factor
                   !                
                   cool_O(x,y,z)=cool_O(x,y,z)+nO(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        326610._dp, 0.27_dp,1._dp, 3.6E+3_dp)/filling_factor
                   !
                   cool_O(x,y,z)=cool_O(x,y,z)+nO(x,y,z,4)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        259130._dp, 2.34_dp,2._dp, 1.5E+3_dp)/filling_factor
                end if
             end if
!
             call solve_state(x,y,z,nion, nO, ion_rate_O, recomb_rate_O, heat_cool_O)                          
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
end subroutine calc_oxygen

!===============================================================================
     subroutine calc_neon()
!
             ! Computes ionization structure of neon and its influence
!              on the temperature structure of the gas.
!
     implicit none
!
       integer(i4b),  parameter              :: nion=4
       integer(i4b)                          :: x, y, z, iteration_counter
       real   (dp ),  dimension(1:nion)      :: ion_rate_Ne=0,&
                                                col_ion_rate_Ne=0,&
                                                recomb_rate_Ne=0,&
                                                heat_cool_Ne=0
       logical                               :: temperature_masked
!        
!===============================================================================
!
!$OMP PARALLEL PRIVATE(x, y, z, iteration_counter,&
!$OMP   ion_rate_Ne, col_ion_rate_Ne,&
!$OMP   recomb_rate_Ne, heat_cool_Ne, temperature_masked)
!$OMP DO 
    do z=-z_max, +z_max
       do y= -y_max,+y_max
          do x=-x_max, x_max
!
             if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1)&
                     cycle
             end if
!
             ne_old=ne(x,y,z)
             col_ion_rate_Ne(1)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.1_dp*5.35e-18_dp*&
                  exp(-(nc_planck*nc_light*1.739E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*1.739E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_Ne(2)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.2_dp*4.16e-18_dp*&
                  exp(-(nc_planck*nc_light*3.314E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*3.314E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_Ne(3)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.3_dp*3.11e-18_dp*&
                  exp(-(nc_planck*nc_light*5.141E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*5.141E+5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             ion_rate_Ne(1) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*1.739E+5_dp, &
                  5.35E-18_dp, 3.77_dp, 1.00_dp)
             ion_rate_Ne(1)=ion_rate_Ne(1)+col_ion_rate_Ne(1)
!       
             ion_rate_Ne(2) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*3.314E+5_dp, &
                  4.16E-18_dp, 2.72_dp, 1.50_dp)
!       
             ion_rate_Ne(2) = ion_rate_Ne(2)+photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:),&
                  nc_light*3.572E+5_dp, &
                  2.17E-18_dp, 2.15_dp, 1.50_dp)
!       
             ion_rate_Ne(2) = ion_rate_Ne(2)+photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:),&
                  nc_light*3.871E+5_dp, &
                  0.52E-18_dp, 2.13_dp, 1.50_dp)
             ion_rate_Ne(2)=ion_rate_Ne(2)+col_ion_rate_Ne(2)
!       
             ion_rate_Ne(3) = photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:), nc_light*5.141E+5_dp, &
                  1.80E-18_dp, 2.28_dp, 2.00_dp)
!       
             ion_rate_Ne(3) = ion_rate_Ne(3)+photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:), nc_light*5.551E+5_dp, &
                  2.50E-18_dp, 2.35_dp, 2.50_dp)
!             
             ion_rate_Ne(3) = ion_rate_Ne(3)+photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:), nc_light*5.763E+5_dp, &
                  1.48E-18_dp, 2.23_dp, 2.50_dp)
             ion_rate_Ne(3)=ion_rate_Ne(3)+col_ion_rate_Ne(3)
!                          
             recomb_rate_Ne(2)=radiative_rec3(2.20E-13_dp,0.759_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ne(2)=recomb_rate_Ne(2)+&
                  dielectronic_hot(9.77E-4_dp,7.30E-2_dp,3.11E+5_dp,2.06E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ne(2)=recomb_rate_Ne(2)/filling_factor
!       
             recomb_rate_Ne(3)=radiative_rec3(1.50E-12_dp,0.693_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ne(3)=recomb_rate_Ne(3)+&
                  dielectronic_hot(2.65E-3_dp,2.42E-1_dp,2.84E+5_dp,3.07E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ne(3)=recomb_rate_Ne(3)+& 
                  dielectronic_cool(0.0129_dp,-0.1779_dp,0.9353_dp,-0.0682_dp,&
                  0.4516_dp,temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ne(3)=recomb_rate_Ne(3)/filling_factor
!             
             recomb_rate_Ne(4)=radiative_rec3(4.40E-12_dp,0.675_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ne(4)=recomb_rate_Ne(4)+&
                  dielectronic_hot(3.69E-3_dp,1.01E+0_dp,2.24E+5_dp,2.94E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)  
             recomb_rate_Ne(4)=recomb_rate_Ne(4)+& 
                  dielectronic_cool(3.6781_dp,14.1481_dp,17.1175_dp,-0.5017_dp,&
                  0.2313_dp,temperature(x,y,z))*ne(x,y,z)            
             recomb_rate_Ne(4)=recomb_rate_Ne(4)/filling_factor
!
             if (charge_transfer) then
! Ro coefficients given for Ne+ in Kindon and Ferland (1997)
!                recomb_rate_Ne(2) = recomb_rate_Ne(2)+ nHI(x,y,z)/&
!                    filling_factor*chargetransfercoeff(&
!                    ,,,,temperature(x,y,z))
                recomb_rate_Ne(3) = recomb_rate_Ne(3)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp, temperature(x,y,z))
                recomb_rate_Ne(4) = recomb_rate_Ne(4)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    14.73_dp,4.52e-2_dp,-0.84_dp,-0.31_dp,temperature(x,y,z))  
             end if
!

             temperature_masked=.false.
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
!
             if (compute_temperature .and. .not. temperature_masked) then             
                heat_cool_Ne=0._dp
!
                if (photoheat_metals) then
                   heat_cool_Ne(1)=heat_cool_Ne(1)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*1.739E+5_dp, &
                        5.35E-18_dp, 3.77_dp, 1.00_dp)
                   heat_cool_Ne(2)=heat_cool_Ne(2)+&
                        photoheatintegral(points-1,j_nu(x,y,z,1:), nug(1:),&
                        nc_light*3.572E+5_dp, &
                        2.17E-18_dp, 2.15_dp, 1.50_dp)
                   heat_cool_Ne(3)=heat_cool_Ne(3)+&
                        photoheatintegral(points-1,j_nu(x,y,z,1:), nug(1:),&
                        nc_light*5.551E+5_dp, &
                        2.50E-18_dp, 2.35_dp, 2.50_dp)    
                end if
!           
                heat_cool_Ne(2)=heat_cool_Ne(2)-&
                     col_fs_cooling(temperature(x,y,z),&
                     ne(x,y,z),128140._dp, 0.28_dp,2._dp, 7.1e+5_dp)/filling_factor
                heat_cool_Ne(3)=heat_cool_Ne(3)-&
                     col_fs_cooling(temperature(x,y,z),&
                     ne(x,y,z),4000._dp, 1.36_dp,9._dp, 6.9e+6_dp)/&
                     filling_factor
                heat_cool_Ne(3)=heat_cool_Ne(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     360200._dp, 0.24_dp,5._dp, 2.1e+5_dp)/filling_factor
                heat_cool_Ne(3)=heat_cool_Ne(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     108600._dp, 0.21_dp,5._dp, 3.1e+4_dp)/filling_factor
                heat_cool_Ne(4)=0._dp
!       
                if (log_heating_cooling) then
                   cool_Ne(x,y,z)=0._dp
                   cool_Ne(x,y,z)=cool_Ne(x,y,z)+nNe(x,y,z,2)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        128140._dp, 0.28_dp,4._dp, 7.1e+5_dp)/filling_factor  
                   cool_Ne(x,y,z)=cool_Ne(x,y,z)+nNe(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        4000._dp, 1.36_dp,9._dp, 6.9e+6_dp)/&
                        filling_factor   
                   cool_Ne(x,y,z)=cool_Ne(x,y,z)+nNe(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),360200._dp,&
                        0.24_dp,5._dp, 2.1e+5_dp)/filling_factor                
                   cool_Ne(x,y,z)=cool_Ne(x,y,z)+nNe(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),108600._dp,&
                        0.21_dp,5._dp, 3.1e+4_dp)/filling_factor
                end if
             end if
! 
            call solve_state(x,y,z,nion, nNe, ion_rate_Ne, recomb_rate_Ne, heat_cool_Ne)                          
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL    
  end subroutine calc_neon

!===============================================================================



!===============================================================================
     subroutine calc_sulfur()
!
!    ! Computes ionization structure of nitrogen and its influence
!      on the temperature structure of the gas.
!
      implicit none
!
       integer(i4b) , parameter              :: nion=4
       integer(i4b):: x, y, z, iteration_counter
!
       real(dp),      dimension (1:nion)     :: ion_rate_S=0,&
                                                col_ion_rate_S=0,&
                                                recomb_rate_S=0,&
                                                heat_cool_S
       logical                               :: temperature_masked
!
!===============================================================================
!
!$OMP PARALLEL PRIVATE(x, y, z, iteration_counter,&
!$OMP   ion_rate_S, col_ion_rate_S,&
!$OMP   recomb_rate_S, heat_cool_S, temperature_masked)
!$OMP DO 
    do z=-z_max, +z_max
       do y= -y_max,+y_max
          do x=-x_max, x_max
!
            if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1)&
                     cycle
             end if
!
             ne_old=ne(x,y,z)
             col_ion_rate_S(1)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.1_dp*3.140E-18_dp*&
                  exp(-(nc_planck*nc_light*83363.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*83363.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_S(2)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.2_dp*6.780E-18_dp*&
                  exp(-(nc_planck*nc_light*189386.1_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*189386.1_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_S(3)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.3*0.363E-18_dp*&
                  exp(-(nc_planck*nc_light*282824.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*282824.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_S(3)=  col_ion_rate_S(3)+ne(x,y,z)*1.55e+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.3*0.364E-18_dp*&
                  exp(-(nc_planck*nc_light*353934.5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*353934.5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_S(3)=  col_ion_rate_S(3)+&
                  ne(x,y,z)*1.55e+13_dp/sqrt(temperature(x,y,z))*&
                  0.3_dp*0.090E-18_dp*&
                  exp(-(nc_planck*nc_light*405699.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*405699.7_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!                    
             ion_rate_S(1) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*83363.7_dp, &
                  3.140E-18_dp,71.00_dp, 2.5_dp)
             ion_rate_S(1)=ion_rate_S(1)+col_ion_rate_S(1)
!       
             ion_rate_S(2) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*189386.1_dp, &
                  6.780E-18_dp, 2.00_dp, 5.00_dp)
             ion_rate_S(2)=ion_rate_S(2)+col_ion_rate_S(2)
!       
             ion_rate_S(3) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*282824.7_dp, &
                  0.363E-18_dp, 2.00_dp, 2.00_dp)
!       
             ion_rate_S(3) = ion_rate_S(3)+photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:),&
                  nc_light*353934.5_dp, &
                  0.364E-18_dp, 2.00_dp, 2.00_dp)
!       
             ion_rate_S(3) = ion_rate_S(3)+photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:), nc_light*405699.7, &
                  0.090E-18_dp, 4.00_dp, 1.75_dp)
             ion_rate_S(3)=ion_rate_S(3)+col_ion_rate_S(3)
!                    
             recomb_rate_S(2)=radiative_rec3(4.10E-13_dp,0.630_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_S(2)=recomb_rate_S(2)+&
                  dielectronic_hot(1.62E-3_dp,0._dp,1.25E+5_dp,1.0E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_S(2)=recomb_rate_S(2)/filling_factor
!       
             recomb_rate_S(3)=radiative_rec3(1.80E-12_dp,0.686_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_S(3)=recomb_rate_S(3)+&
                  dielectronic_hot(1.09E-2_dp,1.20E-2_dp,1.92E+5_dp,1.80E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z) 
             recomb_rate_S(3)=recomb_rate_S(3)/filling_factor
!       
             recomb_rate_S(4)=radiative_rec3(2.70E-12_dp,0.745_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_S(4)=recomb_rate_S(4)+&
                  dielectronic_hot(3.35E-2_dp,6.59E-2_dp,1.89E+5_dp,1.59E+5_dp,&
                  temperature(x,y,z))*ne(x,y,z)              
             recomb_rate_S(4)=recomb_rate_S(4)/filling_factor
!
             if (charge_transfer) then
                ion_rate_S(1) = ion_rate_S(1) + nHII(x,y,z)/&
                     filling_factor*chargetransfercoeff(&
                     1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp,temperature(x,y,z))
                recomb_rate_S(2) = recomb_rate_S(2)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    3.82e-7_dp,11.10_dp,2.57e+4_dp,-8.22_dp,temperature(x,y,z))
                recomb_rate_S(3) = recomb_rate_S(3)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp,temperature(x,y,z))
                recomb_rate_S(4) = recomb_rate_S(4)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    2.29_dp,4.02e-2_dp,1.59_dp,-6.06_dp,temperature(x,y,z))  
             end if
!
             temperature_masked=.false.
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
!
             if (compute_temperature .and. .not. temperature_masked) then             
!
                heat_cool_S=0._dp
!
                if (photoheat_metals) then
                   heat_cool_S(1)=heat_cool_S(1)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*83363.7_dp,&
                        3.140E-18_dp,71.00_dp, 2.5_dp)
                   heat_cool_S(2)=heat_cool_S(2)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*189386.1_dp,&
                        6.780E-18_dp, 2.00_dp, 5.00_dp)
                   heat_cool_S(3)=heat_cool_S(3)+&
                        photoheatintegral(points-1,j_nu(x,y,z,1:), nug(1:),&
                        nc_light*405699.7,&
                        0.090E-18_dp, 4.00_dp, 1.75_dp)  
                end if
!             
                heat_cool_S(1)=0._dp
                heat_cool_S(2)=heat_cool_S(2)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     6720._dp, 6.9_dp,4._dp, 1E+9_dp)/filling_factor
                heat_cool_S(3)=heat_cool_S(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     9600._dp, 6.95_dp,9._dp, 1E+9_dp)/filling_factor
                heat_cool_S(3)=heat_cool_S(3)-&
                     col_fs_cooling(temperature(x,y,z),&
                     ne(x,y,z),334700._dp, 3.98_dp,1._dp, 1E+9_dp)/filling_factor
                heat_cool_S(3)=heat_cool_S(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     120000._dp, 0.27_dp,1._dp, 1E+9_dp)/filling_factor
                heat_cool_S(4)=heat_cool_S(4)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     105150._dp, 6.9_dp,2._dp, 1E+9_dp)/filling_factor
!       
                if (log_heating_cooling) then
                   cool_S(x,y,z)=0._dp
                   cool_S(x,y,z)=cool_S(x,y,z)+nS(x,y,z,2)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        6720._dp,6.9_dp,4._dp, 1E+9_dp)/filling_factor
                   cool_S(x,y,z)=cool_S(x,y,z)+nS(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        9600._dp, 6.95_dp,9._dp, 1E+9_dp)/&
                        filling_factor 
                   cool_S(x,y,z)=cool_S(x,y,z)+nS(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        334700._dp, 3.98_dp,1._dp, 1E+9_dp)/filling_factor
                   cool_S(x,y,z)=cool_S(x,y,z)+nS(x,y,z,4)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        105150._dp, 6.9_dp,2._dp, 1E+9_dp)/filling_factor
                end if
             end if
!
             call solve_state(x,y,z,nion, nS, ion_rate_S, recomb_rate_S, heat_cool_S)                                       
!       
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
     end subroutine calc_sulfur
!

!===============================================================================
     subroutine calc_argon()
!
             ! Computes ionization structure of neon and its influence
!              on the temperature structure of the gas.
!
     implicit none
!
       integer(i4b),  parameter              :: nion=4
       integer(i4b)                          :: x, y, z, iteration_counter
       real   (dp ),  dimension(1:nion)      :: ion_rate_Ar=0,&
                                                col_ion_rate_Ar=0,&
                                                recomb_rate_Ar=0,&
                                                heat_cool_Ar=0
       logical                               :: temperature_masked
!        
!===============================================================================
!
!$OMP PARALLEL PRIVATE(x, y, z, iteration_counter,&
!$OMP   ion_rate_Ar, col_ion_rate_Ar,&
!$OMP   recomb_rate_Ar, heat_cool_Ar, temperature_masked)
!$OMP DO 
    do z=-z_max, +z_max
       do y= -y_max,+y_max
          do x=-x_max, x_max
!
             if (use_mask) then
                if (mod(mask(x,y,z),2) .eq. 1)&
                     cycle
             end if
!
             ne_old=ne(x,y,z)
             col_ion_rate_Ar(1)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.1_dp*0.755E-18_dp*&
                  exp(-(nc_planck*nc_light*127587.0_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*127587.0_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_Ar(2)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.2_dp*28.000E-18_dp*&
                  exp(-(nc_planck*nc_light*222888.2_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*222888.2_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             col_ion_rate_Ar(3)=  ne(x,y,z)*1.55E+13_dp/&
                  sqrt(temperature(x,y,z))*&
                  0.3_dp*1.00E-18_dp*&
                  exp(-(nc_planck*nc_light*329420.5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/&
                  ((nc_planck*nc_light*329420.5_dp)/&
                  (nc_boltzmann*temperature(x,y,z)))/filling_factor 
!       
             ion_rate_Ar(1) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*127587.0, &
                  0.755E-18_dp, 2.00_dp, 2.00_dp)
             ion_rate_Ar(1)=ion_rate_Ar(1)+col_ion_rate_Ar(1)
!       
             ion_rate_Ar(2) = photointegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                  nc_light*3.314E+5_dp, &
                  28.00E-18_dp, 9.50_dp, 5.0_dp)
!       !       
             ion_rate_Ar(3) = photointegral(points-1,&
                  j_nu(x,y,z,1:), nug(1:), nc_light*329420.5_dp, &
                  1.00E-18_dp, 4.00_dp, 2.50_dp)
!                          
             recomb_rate_Ar(2)=radiative_rec3(3.77E-13_dp,0.651_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ar(2)=recomb_rate_Ar(2)+&
                  dielectronic_hot(1.00E-03_dp,5.00E-03_dp,3.20E+05_dp,3.10E+05_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ar(2)=recomb_rate_Ar(2)/filling_factor
!       
             recomb_rate_Ar(3)=radiative_rec3(1.95E-12_dp,0.752_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ar(3)=recomb_rate_Ar(3)+&
                  dielectronic_hot(1.10E-02_dp,4.50E-02_dp,2.90E+05_dp,5.50E+05_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ar(3)=recomb_rate_Ar(3)/filling_factor
!             
             recomb_rate_Ar(4)=radiative_rec3(3.23E-12_dp,0.869_dp,&
                  temperature(x,y,z))*ne(x,y,z)
             recomb_rate_Ar(4)=recomb_rate_Ar(4)+&
                  dielectronic_hot(3.40E-02_dp,5.70E-02_dp,2.39E+05_dp,6.00E+05_dp,&
                  temperature(x,y,z))*ne(x,y,z)  
             recomb_rate_Ar(4)=recomb_rate_Ar(4)/filling_factor
!
             if (charge_transfer) then
!               No coefficients for Ar+ given in Kingddon & Ferland (1997) 
!               recomb_rate_Ar(2) = recomb_rate_C(2)+ nHI(x,y,z)/&
!                    filling_factor*chargetransfercoeff(&
!                    _dp,_dp,_dp,_dp,temperature(x,y,z))
                recomb_rate_Ar(3) = recomb_rate_Ar(3)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    1.00e-5_dp,0.00_dp,0.00_dp,0.00_dp, temperature(x,y,z))
                recomb_rate_Ar(4) = recomb_rate_Ar(4)+ nHI(x,y,z)/&
                    filling_factor*chargetransfercoeff(&
                    4.57_dp,0.27_dp,-0.18_dp,-1.57_dp,temperature(x,y,z))  
             end if
!
             temperature_masked=.false.
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
!
             if (compute_temperature .and. .not. temperature_masked) then             
                heat_cool_Ar=0._dp
!
                if (photoheat_metals) then
                   heat_cool_Ar(1)=heat_cool_Ar(1)+&
                        photoheatintegral(points-1, j_nu(x,y,z,1:), nug(1:),&
                        nc_light*1.739E+5_dp, &
                        5.35E-18_dp, 3.77_dp, 1.00_dp)
                   heat_cool_Ar(2)=heat_cool_Ar(2)+&
                        photoheatintegral(points-1,j_nu(x,y,z,1:), nug(1:),&
                        nc_light*3.572E+5_dp, &
                        2.17E-18_dp, 2.15_dp, 1.50_dp)
                   heat_cool_Ar(3)=heat_cool_Ar(3)+&
                        photoheatintegral(points-1,j_nu(x,y,z,1:), nug(1:),&
                        nc_light*5.551E+5_dp, &
                        2.50E-18_dp, 2.35_dp, 2.50_dp)    
                end if
!           
                heat_cool_Ar(2)=heat_cool_Ar(2)-&
                     col_fs_cooling(temperature(x,y,z),&
                     ne(x,y,z),128140._dp, 0.28_dp,2._dp, 7.1e+5_dp)/filling_factor
                heat_cool_Ar(3)=heat_cool_Ar(3)-&
                     col_fs_cooling(temperature(x,y,z),&
                     ne(x,y,z),4000._dp, 1.36_dp,9._dp, 6.9e+6_dp)/&
                     filling_factor
                heat_cool_Ar(3)=heat_cool_Ar(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     360200._dp, 0.24_dp,5._dp, 2.1e+5_dp)/filling_factor
                heat_cool_Ar(3)=heat_cool_Ar(3)-&
                     col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                     108600._dp, 0.21_dp,5._dp, 3.1e+4_dp)/filling_factor
                heat_cool_Ar(4)=0._dp
!       
                if (log_heating_cooling) then
                   cool_Ar(x,y,z)=0._dp
                   cool_Ar(x,y,z)=cool_Ar(x,y,z)+nAr(x,y,z,2)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        128140._dp, 0.28_dp,4._dp, 7.1e+5_dp)/filling_factor  
                   cool_Ar(x,y,z)=cool_Ar(x,y,z)+nAr(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),&
                        4000._dp, 1.36_dp,9._dp, 6.9e+6_dp)/&
                        filling_factor   
                   cool_Ar(x,y,z)=cool_Ar(x,y,z)+nAr(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),360200._dp,&
                        0.24_dp,5._dp, 2.1e+5_dp)/filling_factor                
                   cool_Ar(x,y,z)=cool_Ar(x,y,z)+nAr(x,y,z,3)*&
                        col_fs_cooling(temperature(x,y,z), ne(x,y,z),108600._dp,&
                        0.21_dp,5._dp, 3.1e+4_dp)/filling_factor
                end if
             end if
! 
            call solve_state(x,y,z,nion, nAr, ion_rate_Ar, recomb_rate_Ar, heat_cool_Ar)                          
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL    
  end subroutine calc_argon
!     
  end subroutine CALC_METALS

!===============================================================================
  subroutine APPLY_ADD_HEATING
!  
! Applies the additional heating rates as a function of the temperature
!
    use M_definitions,       only: x_max, y_max, z_max, use_mask, t_step
    use M_grid_memory,       only: add_heating_poly, energycontent, mask, ne
!
    implicit none
    integer(i4b)                :: x,y,z
    logical                     :: temperature_masked
!===============================================================================
!
!$OMP PARALLEL PRIVATE (x,y,z, temperature_masked)
!$OMP DO
    do z=-z_max,+z_max
       do y=-y_max,+y_max
          do x=-x_max, +x_max             
             temperature_masked=.false.

             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     temperature_masked=.true.  
             end if
!
             if (.not. temperature_masked) &
                  energycontent(x,y,z)=energycontent(x,y,z)+&
                 (ne(x,y,z)**3*add_heating_poly(1,x,y,z)+&
                  ne(x,y,z)**2*add_heating_poly(2,x,y,z)+&
                  ne(x,y,z)   *add_heating_poly(3,x,y,z)+&
                  add_heating_poly(4,x,y,z))*&
                  t_step
          end do
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL
  end subroutine APPLY_ADD_HEATING

!================================================================================
  subroutine APPLY_EXTRA_HEATING_COOLING
!
! Reads the heating (and/or cooling rates) form the file extraheating.bin    
!
    use M_definitions,       only: x_max, y_max, z_max, use_mask, binary_output
    use M_grid_memory,       only: energycontent
!
    implicit none
!
    integer(i4b)                :: x,y,z
    real(dp)                    :: heatvalue
    real(sp)                    :: heatvaluef   ! For reaing the singe-precision
                                                ! floating point binary files 
    if (binary_output) then
       open(unit=100, file='extra_hc.bin', form='unformatted')
       read(100) x,y,z
       if (x .ne. 2*x_max+1 .or. y .ne. 2*y_max+1 .or. z .ne. 2*z_max+1)&
            stop 'size of extra_hc.txt does not match'
       read(100) heatvaluef
       energycontent(x,y,z)=energycontent(x,y,z)+real(heatvalue,8)
       close(100)    
    else
       open(unit=100, file='extra_hc.txt')
       read(100,*) x,y,z
       if (x .ne. 2*x_max+1 .or. y .ne. 2*y_max+1 .or. z .ne. 2*z_max+1)&
            stop 'size of extra_hc.txt does not match'
       do z=-z_max,+z_max
          do y=-y_max,+y_max
             do x=-x_max, +x_max
                read(100,*) heatvalue
                energycontent(x,y,z)=energycontent(x,y,z)+heatvalue
             end do
          end do
       end do
       close(100)
    end if
!================================================================================    
  end subroutine APPLY_EXTRA_HEATING_COOLING


!=============================================================================== 
  subroutine UPDATE_ELECTRON_DENSITY
!   
!    Updates the array ne
!
!   called by the main program
!
    use M_definitions,       only: include_metals
    use M_grid_memory,       only: ne,nHII,nHeII,nHEIII,nC,nN,nO,nNe,nS, nAr
!
  implicit none
!
    ne=nHII+nHeII+2*nHeIII
    if (include_metals) then
       ne=ne+nC(:,:,:,2)+nN(:,:,:,2)+nO(:,:,:,2)+nNe(:,:,:,2)+nS(:,:,:,2)+&
            nAr(:,:,:,2)
       ne=ne+2*(nC(:,:,:,3)+nN(:,:,:,3)+nO(:,:,:,3)+nNe(:,:,:,3)+nS(:,:,:,3)+&
            nAr(:,:,:,3))
       ne=ne+3*(nC(:,:,:,4)+nN(:,:,:,4)+nO(:,:,:,4)+nNe(:,:,:,4)+nS(:,:,:,4)+&
            nAr(:,:,:,4))
    end if
  end subroutine UPDATE_ELECTRON_DENSITY



!================================================================================  
  subroutine UPDATE_TEMPERATURE
! 
!    Updates the temperature and the energy content of the grid cells
!
! called by the main program
! 
    use M_natural_constants, only: nc_boltzmann
    use M_definitions,       only: x_max, y_max, z_max,&
                                   artificially_fix_temperature, &
                                   artificial_temperature, t_step,&
                                   use_mask
    use M_grid_memory,       only: temperature, temperature_old,&
                                   energycontent,ne, nH_complete,&
                                   nHe_complete, nC, nN, nO, nNe, nS,&
                                   mask
!
  implicit none
!
    integer(i4b):: x,y,z
!================================================================================
    write(*,*) "update temperature called"
    if (artificially_fix_temperature) then
       do z=-z_max, +z_max
          do y=-y_max,+y_max
             do x=-x_max, +x_max
                temperature(x,y,z)=artificial_temperature
                energycontent(x,y,z)= 3._dp/2._dp*nc_boltzmann*&
                     artificial_temperature*&
                     (nH_complete(x,y,z)+nHe_complete(x,y,z)+ne(x,y,z)+&
                     nC(x,y,z,0)+nN(x,y,z,0)+nO(x,y,z,0)+nNe(x,y,z,0)+nS(x,y,z,0))
             end do
          end do
       end do
    end if
    temperature_old=temperature
    do z=-z_max, +z_max
       do y=-y_max,+y_max
          do x=-x_max, +x_max
             if (use_mask) then
                if (mask(x,y,z) .gt. 0)&
                     cycle
             end if
! 
             energycontent(x,y,z)=energycontent(x,y,z)-ffcooling(x,y,z)*t_step
!
!  Currenty, the temperature is forced to be in the interval 3000 K ... 50000K
             temperature(x,y,z)=max(energycontent(x,y,z)*(2._dp/3._dp)/&
                  (nc_boltzmann*(ne(x,y,z)+nH_complete(x,y,z)+&
                  nHe_complete(x,y,z)+&
                  nC(x,y,z,0)+nN(x,y,z,0)+nO(x,y,z,0)+nNe(x,y,z,0)+nS(x,y,z,0))),&
                  3000._dp)
             temperature(x,y,z)=min(temperature(x,y,z),50000._dp)
!
             if (temperature(x,y,z) .lt. (3000._dp+1E-6_dp)) then
                energycontent(x,y,z)= 3._dp/2._dp*nc_boltzmann*3000._dp*&
                     (nH_complete(x,y,z)+nHe_complete(x,y,z)+ne(x,y,z)+&
                     nC(x,y,z,0)+nN(x,y,z,0)+nO(x,y,z,0)+nNe(x,y,z,0)+&
                     nS(x,y,z,0))               
             end if
             if (temperature(x,y,z) .gt. (50000._dp+1E-3_dp)) then
                energycontent(x,y,z)= 3._dp/2._dp*nc_boltzmann*50000._dp*&
                     (nH_complete(x,y,z)+nHe_complete(x,y,z)+ne(x,y,z)+&
                     nC(x,y,z,0)+nN(x,y,z,0)+nO(x,y,z,0)+nNe(x,y,z,0)+nS(x,y,z,0))
             end if
          end do
       end do
    end do
    !  read(*,*)
  end subroutine  UPDATE_TEMPERATURE

!===============================================================================
  subroutine UPDATE_THERMAL_PRESSURE
!
!   computes the thermal pressure (in units ofdyn/cm^-2 = 0.1 Pa) 
!   for each grid cell
!
! called by the main program
! 
    use M_natural_constants,    only: nc_boltzmann, nc_parsec
    use M_definitions,          only: x_max, y_max, z_max
    use M_grid_memory,          only: ne, nH_complete, nHe_complete,&
                                      nC, nN, nO, nNe, nS,&
                                      temperature, thermal_pressure
! 
  implicit none
!
    integer(i4b) :: x, y, z
    real   ( dp) :: total_particle_density
!    
!================================================================================
    do z=-z_max,+z_max
       do y=-y_max,+y_max
          do x=-x_max, +x_max
             total_particle_density=ne(x, y, z)+&
                  nH_complete(x,y,z)+nHe_complete(x,y,z)+&
                  nC(x,y,z,0)+nN(x,y,z,0)+nO(x,y,z,0)+nNe(x,y,z,0)+nS(x,y,z,0)
             thermal_pressure(x,y,z)=total_particle_density*&
                  nc_boltzmann*temperature(x,y,z)
          end do
       end do
    end do
  end subroutine UPDATE_THERMAL_PRESSURE

!================================================================================  
  real(dp) function FFCOOLING(x,y,z)
!
!    computes the cooling by free-free radiation, considered is the interaction 
!    between hydrogen and helium ions on the one hand and electons on the other 
!    hand
!
! called by UPDATE_TEMPERATURE 
!
    use M_definitions,    only: filling_factor, log_heating_cooling
    use M_grid_memory,    only: ne, nHII, nHeII, nHeIII, temperature, cool_ff
!
  implicit none
!
    integer(i4b)::x,y,z
!================================================================================    
    ffcooling=1.42E-27_dp*sqrt(temperature(x,y,z))*1.3_dp*ne(x,y,z)*nHII(x,y,z)+&
              1.42E-27_dp*sqrt(temperature(x,y,z))*1.3_dp*ne(x,y,z)*nHeII(x,y,z)+&
            4*1.42E-27_dp *sqrt(temperature(x,y,z))*1.3_dp*ne(x,y,z)*nHeIII(x,y,z)
    ffcooling=ffcooling/filling_factor
    if (log_heating_cooling)&
         cool_ff(x,y,z)=ffcooling
  end function FFCOOLING



!================================================================================
  real(dp) function RECHII(temperature)                     ! deprectiated
!
!    Computes the recombination of HII.
!
! Currently not used, HYDREC is used instead
!
  implicit none
    real(dp):: temperature
!================================================================================
    recHII=exp(-21.31752118_dp-0.833338575_dp*log(temperature))
  end function RECHII



!================================================================================  
  real(dp) function RECHEII(temperature)                    ! depreciated
! 
!    Computes the recombination of He II
!
! Curently not used, HELREC ist used instead
!
   implicit none
!
    real(dp):: temperature
!================================================================================
    recHeII=exp(-21.66825983_dp-0.769391963_dp*log(temperature))
  end function recHeII



!================================================================================  
  real(dp) function RECHEIII(temperature)                   ! depreciated
! 
!    Computes the recombination of He II
!
! Curently not used, HELREC ist used instead
!
  implicit none
    real(dp):: temperature
!================================================================================
    recHeIII=exp(-16.89734798_dp-1.138766988_dp*log(temperature))
  END FUNCTION recHeIII



!================================================================================  
  real(dp) function SEATON (nu,nu_0, alpha_0,beta, s)
!
!    Computes the frequency-dependent Seaton photoionization cross-sections
!
!  called by PHOTOINTEGRAL, PHOTOHEATINTEGRAL, UPDATE_CHI
!    
  implicit none
    real(dp) :: nu, nu_0, alpha_0, beta, s
!================================================================================
    if (nu .ge. nu_0) then
       seaton = alpha_0*(beta*(nu/nu_0)**(-s)+ (1._dp-beta)*(nu/nu_0)**(-s-1))
    else
       seaton=0._dp
    end if
    !  print*, "SEATON", nu, seaton
  end function SEATON
  
!================================================================================
  real(dp) function RADIATIVE_REC0 (alphaR10000,T)          ! depreciated
!
!    Computes radiative recombination rates
!
! Currently not used. Instead, RADIATIVE_REC3 is used
!
  implicit none
    real (dp) :: alphaR10000,T
!================================================================================
!
    radiative_rec0=alphaR10000*(T/10000)**(-0.8)
  END FUNCTION RADIATIVE_REC0



!================================================================================  
  real(dp) function RADIATIVE_REC2(alphaR10000,f,T)         ! depreciated
!
!    Computes radiative recombination rates
!
! Currently not used. Instead, RADIATIVE_REC3 is used
!
  implicit none
    real(dp):: alphaR10000,f,T
    real(dp):: phi
!================================================================================
    phi=(T/10000._dp)**(-0.4_dp)
    radiative_rec2=alphaR10000*sqrt(T/10000._dp)*(f+(1._dp-f)*phi)
  END FUNCTION RADIATIVE_REC2



!===============================================================================  
  real(dp) function RADIATIVE_REC3(arad,eta,t)
!
!    Computes radiative recombination rates 
!    (cf. www.pa.uky.edu/~verner/rec.html and the references therein)
!
! used by CALC_METALS
!
  implicit none
    real(dp):: arad, eta, t
!===============================================================================
    radiative_rec3=arad*(t/1E+4_dp)**(-eta)
  end function RADIATIVE_REC3



!================================================================================
  real(dp) function COL_FS_COOLING  &
       (temperature, electrondensity, lambda, ypsilon, omega1, crit)
!
!    Computes the cooling by the radiative decay of collisionally excited
!    states 
!
! called by CALC_HYDROGEN, CALC_METALS
!
    use M_natural_constants,  only: nc_planck, nc_light, nc_boltzmann
!
  implicit none

    real(dp) :: electrondensity, temperature, lambda, ypsilon, omega1, crit
!================================================================================
!
    col_fs_cooling=&
         ELECTRONDENSITY*8.629E-6_dp/sqrt(temperature)*YPSILON/OMEGA1*&
         exp(-(nc_planck*nc_light/(LAMBDA*1E-8))  / (nc_boltzmann*TEMPERATURE))*&
         (nc_planck*nc_light/(lambda*1E-8))*&
         (crit)/(electrondensity+crit)
  end function COL_FS_COOLING


  
!===============================================================================  
  real(dp) function DIELECTRONIC_REC(a,b,T)               ! deprecated
!
!    Computes dielectronic recombination rates.
!
! currently not used, Instead: DIELETRONIC_HOT, DIELECTRONIC_COOL
!
  implicit none
!
    real(dp):: a,b,T
!===============================================================================
!
    dielectronic_rec=10._dp**(a+b*log10(T))
  end function DIELECTRONIC_REC
  


!===============================================================================
  real(dp) function DIELECTRONIC_HOT(a,b,t0,t1,t)
!
!    Computes high-temperature dielectronic recombination rates
!    (cf. www.pa.uky.edu/~verner/rec.html and the references therein)
!
! called by CALC_METALS 
!
  implicit none
!
    real(dp) :: a,b,t0,t1,t
!==============================================================================
!
    dielectronic_hot=a*t**(-3._dp/2._dp)*exp(-t0/t)*(1+B*exp(-t1/t))
  end function DIELECTRONIC_HOT



!===============================================================================
  real(dp) function DIELECTRONIC_COOL(a,b,c,d,f,t)
!
!   Computes high-temperature dielectronic recombination rates    
!   (cf. www.pa.uky.edu/~verner/rec.html and the references therein)
!
! called by CALC_METALS
!
  implicit none
!
  real(dp) :: a,b,c,d,f,t, t10000
!================================================================================
!
    t10000=t/10000.0_dp
    dielectronic_cool=1E-12_dp*(a/t10000+b+c*t10000+d*t10000**2)*&
         t10000**(-3._dp/2._dp)*exp(-f/t10000)
  end function DIELECTRONIC_COOL



!===============================================================================  
  real(dp) function PHOTOINTEGRAL (N, J, nu, nu_0, alpha_0,beta, s)
!
!   computes the photoionization rate (cf. Hoffmann et al. 2012)
!
! called by CALC_HYDROGEN, CALC_HELIUM2, CALC_METALS
!
    use M_natural_constants,  only:  nc_pi, nc_planck
    use M_data_input,         only : nug
!
  implicit none
!
    integer(i4b) ::N, frequencycounter
    real(dp)     :: nu_0 , alpha_0 , beta, s
!
    real(dp), dimension (:) :: J, nu
!===============================================================================
!
    photointegral=0._dp
    do frequencycounter=1, N-1
       photointegral = photointegral+(1._dp/2._dp)*&
            (J(frequencycounter)*&
            seaton(nug(frequencycounter), nu_0,alpha_0, beta,s)&
            /(nc_planck*nu(frequencycounter))+&
            J(frequencycounter+1)*&
            seaton(nug(frequencycounter+1), nu_0,alpha_0, beta,s)/&
            (nc_planck*nu(frequencycounter+1)))*&
            (nu(frequencycounter)-nu(frequencycounter+1))
    end do
    photointegral=photointegral*4*nc_pi
  end function PHOTOINTEGRAL
 
!================================================================================
!!$  real(dp) function CHARGETRANSFERION(iontype,x,y,z)
!!$!
!!$! Computes the ionization rate coefficients due to charge exchange processes
!!$! The argument is the considered ion
!!$!
!!$    use M_grid_memory,        only: nHI, nHII, nO
!!$!
!!$    implicit none

!!$!
!!$    character (len=*) :: iontype
!!$    integer(i4b)      :: x,y,z
!!$!
!!$    if (iontype .eq. "HI") then
!!$       chargetransferion=1.01E-9_dp*nO(x,y,z,2)
!!$    end if
!!$    if (iontype .eq. "OI") then
!!$       chargetransferion=9.10E-10_dp*nHII(x,y,z)
!!$    end if
!!$  end function CHARGETRANSFERION
!!$!
!!$!================================================================================  
!!$  real(dp) function CHARGETRANSFERRECOMB(iontype,x,y,z)
!!$!
!!$! Computes the recombination rate coefficients due to charge exchange processes
!!$! The argument is the considered ion
!!$!
!!$    use M_grid_memory,        only: nHI, nHII, nO
!!$!
!!$    implicit none
!!$!
!!$    character (len=*) :: iontype
!!$    integer(i4b)      :: x,y,z
!!$!
!!$    if (iontype .eq. "HII") then
!!$       chargetransferrecomb=9.10E-10_dp*nO(x,y,z,1)
!!$    end if
!!$    if (iontype .eq. "OII") then
!!$       chargetransferrecomb=1.01E-9_dp*nHI(x,y,z)
!!$    end if 
!!$ end function CHARGETRANSFERRECOMB
!
!===============================================================================
 real(dp) function CHARGETRANSFERCOEFF(a,b,c,d,t)
!
! Computes the charge transfer ionization and recombination rate coefficients 
! accoring to the fit formulae in J.B Kingdon and G.J. Ferland 
! (1996ApJS...100..205K) 
!
   implicit none
!
   real(dp)           ::  a,b,c,d,t      ! name of the coefficients in the paper
   real(dp)           ::  t4             ! t in units of 10000 K
!
   t4=t/1e+4_dp
   chargetransfercoeff = (a*1e-9_dp)*t4**b*(1._dp+c*exp(d*t4))
 end function CHARGETRANSFERCOEFF
!
!================================================================================ 
  real(dp) function PHOTOHEATINTEGRAL (N, J, nu, nu_0, alpha_0,beta, s)
!  
!    Computes the energy input by photoheating.
!
! called_by CALC_HYDROGEN, CALC_HELIUM2
!
    use M_natural_constants,  only: nc_planck, nc_pi                      
    use M_data_input,         only: nug
!
  implicit none
!
    integer(i4b):: N, frequencycounter
    real(dp)    :: nu_0 , alpha_0 , beta, s
!
    real(dp), dimension (:) :: J, nu

!================================================================================
!
    photoheatintegral=0._dp
    do frequencycounter=1, N-1
       photoheatintegral = photoheatintegral+(1._dp/2._dp)*&
            (J(frequencycounter)*&
            (nc_planck*nu(frequencycounter)-nc_planck*nu_0)*&
            seaton(nug(frequencycounter), nu_0,alpha_0, beta,s)/&
            (nc_planck*nu(frequencycounter))+&
            J(frequencycounter)*&
            (nc_planck*nu(frequencycounter+1)-nc_planck*nu_0)*&
            seaton(nug(frequencycounter+1), nu_0,alpha_0, beta,s)&
            /(nc_planck*nu(frequencycounter+1)))*&
            (nu(frequencycounter)-nu(frequencycounter+1))
    end do
    photoheatintegral=photoheatintegral*4*nc_pi
  end function PHOTOHEATINTEGRAL


  
!===============================================================================
  real(dp) function HYDREC(T,N)
    !
    ! Computes total recombination rates and total recombination cooling
    ! coefficients for hydrogen by interpolating the table of
    ! D.G.Hummer, MNRAS 268, 109-112 (1994)
    ! Arguments: T is temperature
    !            N is 1 for alpha_A
    !                 2 for beta_A
    !                 3 for beta_ff
    !                 4 for alpha_b
    !                 5 for beta_B 
    ! Table entries are:
    ! Log T, Q*alpha_1, Q*alpha_B, Q*beta_1, Q*beta_B, Q*beta_ff, Q*beta_B_tot
    ! where Q=Sqrt(T) and beta_B_tot=beta_B+beta_ff
    ! M 3.05 ADI
    !
    !  called by CALC_HYDROGEN
    !
    !input
  IMPLICIT NONE
    !
    REAL   (dp )  :: T
    INTEGER(i4b)  :: N
    !local
    INTEGER(i4b)  :: I,J
    REAL   (dp )  :: Q,S,A,X
    REAL   (dp ), DIMENSION(31,7) :: C
!===============================================================================
!
    DATA ((C(I,J),J=1,7),I=1,31)/ &
         1.0_dp,1.646E-11_dp,9.283E-11_dp,1.646E-11_dp,8.287E-11_dp,1.061E-11_dp,9.348E-11_dp,&
         1.2_dp,1.646E-11_dp,8.823E-11_dp,1.646E-11_dp,7.821E-11_dp,1.068E-11_dp,8.889E-11_dp,&
         1.4_dp,1.646E-11_dp,8.361E-11_dp,1.646E-11_dp,7.356E-11_dp,1.076E-11_dp,8.432E-11_dp,&
         1.6_dp,1.646E-11_dp,7.898E-11_dp,1.646E-11_dp,6.892E-11_dp,1.085E-11_dp,7.977E-11_dp,&
         1.8_dp,1.646E-11_dp,7.435E-11_dp,1.646E-11_dp,6.430E-11_dp,1.095E-11_dp,7.525E-11_dp,&
         2.0_dp,1.646E-11_dp,6.973E-11_dp,1.645E-11_dp,5.971E-11_dp,1.106E-11_dp,7.077E-11_dp,&
         2.2_dp,1.645E-11_dp,6.512E-11_dp,1.644E-11_dp,5.515E-11_dp,1.118E-11_dp,6.633E-11_dp,&
         2.4_dp,1.645E-11_dp,6.054E-11_dp,1.643E-11_dp,5.062E-11_dp,1.132E-11_dp,6.194E-11_dp,&
         2.6_dp,1.644E-11_dp,5.599E-11_dp,1.641E-11_dp,4.614E-11_dp,1.145E-11_dp,5.758E-11_dp,&
         2.8_dp,1.642E-11_dp,5.147E-11_dp,1.638E-11_dp,4.170E-11_dp,1.161E-11_dp,5.332E-11_dp,&
         3.0_dp,1.640E-11_dp,4.700E-11_dp,1.633E-11_dp,3.734E-11_dp,1.181E-11_dp,4.915E-11_dp,&
         3.2_dp,1.636E-11_dp,4.258E-11_dp,1.625E-11_dp,3.306E-11_dp,1.202E-11_dp,4.508E-11_dp,&
         3.4_dp,1.629E-11_dp,3.823E-11_dp,1.613E-11_dp,2.888E-11_dp,1.224E-11_dp,4.112E-11_dp,&
         3.6_dp,1.620E-11_dp,3.397E-11_dp,1.594E-11_dp,2.484E-11_dp,1.248E-11_dp,3.733E-11_dp,&
         3.8_dp,1.605E-11_dp,2.983E-11_dp,1.565E-11_dp,2.098E-11_dp,1.274E-11_dp,3.373E-11_dp,&
         4.0_dp,1.582E-11_dp,2.584E-11_dp,1.522E-11_dp,1.736E-11_dp,1.303E-11_dp,3.039E-11_dp,&
         4.2_dp,1.548E-11_dp,2.204E-11_dp,1.460E-11_dp,1.402E-11_dp,1.335E-11_dp,2.737E-11_dp,&
         4.4_dp,1.499E-11_dp,1.846E-11_dp,1.374E-11_dp,1.103E-11_dp,1.369E-11_dp,2.472E-11_dp,&
         4.6_dp,1.431E-11_dp,1.520E-11_dp,1.260E-11_dp,8.442E-12_dp,1.403E-11_dp,2.247E-11_dp,&
         4.8_dp,1.341E-11_dp,1.226E-11_dp,1.119E-11_dp,6.279E-12_dp,1.434E-11_dp,2.062E-11_dp,&
         5.0_dp,1.227E-11_dp,9.696E-12_dp,9.571E-12_dp,4.539E-12_dp,1.460E-11_dp,1.914E-11_dp,&
         5.2_dp,1.093E-11_dp,7.514E-12_dp,7.844E-12_dp,3.192E-12_dp,1.478E-11_dp,1.797E-11_dp,&
         5.4_dp,9.454E-12_dp,5.710E-12_dp,6.146E-12_dp,2.185E-12_dp,1.485E-11_dp,1.704E-11_dp,&
         5.6_dp,7.920E-12_dp,4.257E-12_dp,4.601E-12_dp,1.458E-12_dp,1.482E-11_dp,1.628E-11_dp,&
         5.8_dp,6.427E-12_dp,3.117E-12_dp,3.295E-12_dp,9.484E-13_dp,1.468E-11_dp,1.563E-11_dp,&
         6.0_dp,5.058E-12_dp,2.244E-12_dp,2.262E-12_dp,6.023E-13_dp,1.444E-11_dp,1.505E-11_dp,&
         6.2_dp,3.866E-12_dp,1.590E-12_dp,1.494E-12_dp,3.738E-13_dp,1.414E-11_dp,1.451E-11_dp,&
         6.4_dp,2.877E-12_dp,1.110E-12_dp,9.520E-13_dp,2.268E-13_dp,1.380E-11_dp,1.402E-11_dp,&
         6.6_dp,2.089E-12_dp,7.642E-13_dp,5.878E-13_dp,1.348E-13_dp,1.344E-11_dp,1.358E-11_dp,&
         6.8_dp,1.485E-12_dp,5.199E-13_dp,3.528E-13_dp,7.859E-14_dp,1.311E-11_dp,1.318E-11_dp,&
         7.0_dp,1.036E-12_dp,3.498E-13_dp,2.066E-13_dp,4.499E-14_dp,1.280E-11_dp,1.285E-11_dp/
    !
    !--------------------------------------------------------------------------------------
    ! statement function for table interpolation
    X(A,I,J)=C(I,J)+A*(C(I+1,J)-C(I,J))
    !-------------------------------------------------------------------------
    !
    Q=SQRT (T)
    S=LOG10(T)
    IF(S.GT.7.0.OR.S.LT.1.0) THEN
       print*,T,S,log10(T)
       STOP 'HYDREC: Temperature out of range'
    END IF
    DO I=1,30
       IF(S.GE.C(I,1).AND.S.LE.C(I+1,1)) GOTO 1
    END DO
    STOP 'HYDREC: Unexpected error'
1   CONTINUE
    A=(S-C(I,1))/(C(I+1,1)-C(I,1))
    IF    (N.EQ.1) THEN
       HYDREC=(X(A,I,2)+X(A,I,3))/Q
    ELSEIF(N.EQ.2) THEN
       HYDREC=(X(A,I,4)+X(A,I,5))/Q
    ELSEIF(N.EQ.3) THEN
       HYDREC=X(A,I,6)/Q
    ELSEIF(N .EQ. 4) THEN
       HYDREC=X(A,I,3)/Q          
    ELSEIF(N .EQ. 5) THEN
       HYDREC=X(A,I,5)/Q
    ELSE
       STOP 'HYDREC: Bad parameter'
    ENDIF
    RETURN
  end function HYDREC
  


!================================================================================
  real(dp) function HELREC(T,N)
    !
    ! Computes total recombination rates and total recombination cooling
    ! coefficients for hydrogen by interpolating the table of
    ! D.G.Hummer, MNRAS 268, 109-112 (1994)
    ! Arguments: T is temperature
    !            N is 1 for alpha_A
    !                 2 for beta_A
    !                 3 for beta_ff
    !                 4 for alpha_b
    !                 5 for beta_B 
    
    ! Table entries are:
    ! Log T, Q*alpha_1, Q*alpha_B, Q*beta_1, Q*beta_B, Q*beta_ff, Q*beta_B_tot
    ! where Q=Sqrt(T) and beta_B_tot=beta_B+beta_ff
    ! 
    ! called by CALC_HELIUM2
    !
   IMPLICIT NONE
    !
    !input
    REAL   (dp )  :: T
    INTEGER(i4b)  :: N
    !local
    INTEGER(i4b)  :: I,J
    REAL   (dp )  :: Q,S,A,X
    REAL   (dp ), DIMENSION(18,7) :: C
!===============================================================================
    !
    DATA ((C(I,J),J=1,7),I=1,18)/ &
         !        alpha_1       alpha_B       beta_1        beta_b         beta_ff        beta_b_tot
         1.0_dp, 1.569E-11_dp, 9.284E-11_dp, 1.569e-11_dp, 8.347E-11_dp,  1.061E-11_dp,  9.408e-11_dp,&
         1.2_dp, 1.569E-11_dp, 8.847E-11_dp, 1.569E-11_dp, 7.889E-11_dp,  1.068E-11_dp,  8.957e-11_dp,&
         1.4_dp, 1.569E-11_dp, 8.403E-11_dp, 1.569E-11_dp, 7.430E-11_dp,  1.076E-11_dp,  8.506e-11_dp,&
         1.6_dp, 1.599E-11_dp, 7.952E-11_dp, 1.569E-11_dp, 6.971E-11_dp,  1.085E-11_dp,  8.056e-11_dp,&
         1.8_dp, 1.569E-11_dp, 7.499E-11_dp, 1.569E-11_dp, 6.512E-11_dp,  1.095E-11_dp,  7.607e-11_dp,&
         2.0_dp, 1.569E-11_dp, 7.044E-11_dp, 1.569E-11_dp, 6.056E-11_dp,  1.106E-11_dp,  7.162e-11_dp,&
         2.2_dp, 1.569E-11_dp, 6.589E-11_dp, 1.570E-11_dp, 5.603E-11_dp,  1.118E-11_dp,  6.721e-11_dp,&
         2.4_dp, 1.560E-11_dp, 6.136E-11_dp, 1.570E-11_dp, 5.154E-11_dp,  1.132E-11_dp,  6.286e-11_dp,&
         2.6_dp, 1.570E-11_dp, 5.685E-11_dp, 1.571E-11_dp, 4.710E-11_dp,  1.145E-11_dp,  5.855e-11_dp,&
         2.8_dp, 1.571E-11_dp, 5.238E-11_dp, 1.572E-11_dp, 4.274E-11_dp,  1.162E-11_dp,  5.436e-11_dp,&
         3.0_dp, 1.572E-11_dp, 4.797E-11_dp, 1.574E-11_dp, 3.847E-11_dp,  1.181E-11_dp,  5.028e-11_dp,&
         3.2_dp, 1.573E-11_dp, 4.364E-11_dp, 1.578E-11_dp, 3.431E-11_dp,  1.202E-11_dp,  4.633e-11_dp,&
         3.4_dp, 1.576E-11_dp, 3.940E-11_dp, 1.583E-11_dp, 3.031E-11_dp,  1.224E-11_dp,  4.255e-11_dp,&
         3.6_dp, 1.580E-11_dp, 3.528E-11_dp, 1.591E-11_dp, 2.650E-11_dp,  1.248E-11_dp,  3.898e-11_dp,&
         3.8_dp, 1.586E-11_dp, 3.132E-11_dp, 1.602E-11_dp, 2.291E-11_dp,  1.274E-11_dp,  3.565e-11_dp,&
         4.0_dp, 1.595E-11_dp, 2.755E-11_dp, 1.619E-11_dp, 1.960E-11_dp,  1.336E-11_dp,  3.296e-11_dp,&
         4.2_dp, 1.608E-11_dp, 2.401E-11_dp, 1.641E-11_dp, 1.660E-11_dp,  1.336E-11_dp,  2.996e-11_dp,&
         4.4_dp, 1.626E-11_dp, 2.073E-11_dp, 1.670E-11_dp, 1.394E-11_dp,  1.369E-11_dp,  2.763e-11_dp/
    
    !--------------------------------------------------------------------------------------
    ! statement function for table interpolation
    X(A,I,J)=C(I,J)+A*(C(I+1,J)-C(I,J))
    !-------------------------------------------------------------------------
    !
    Q=SQRT (T)
    S=LOG10(T)
    
    
    IF(S .LT. 4.4_dp .AND. S.GT.1.0_dp) THEN
       
       DO I=1,17
          IF(S.GE.C(I,1).AND.S.LE.C(I+1,1)) GOTO 1
       END DO
       STOP 'HELREC: Unexpected error'
1      CONTINUE
       A=(S-C(I,1))/(C(I+1,1)-C(I,1))
       
       
       IF    (N.EQ.1) THEN
          HELREC=(X(A,I,2)+X(A,I,3))/Q
       ELSEIF(N.EQ.2) THEN
          HELREC=(X(A,I,4)+X(A,I,5))/Q
       ELSEIF(N.EQ.3) THEN
          HELREC=X(A,I,6)/Q
       ELSEIF(N .EQ. 4) THEN
          HELREC=X(A,I,3)/Q
       ELSEIF(N .EQ. 5) THEN
          HELREC=X(A,I,5)/Q
       ELSE
          STOP 'HYDREC: Bad parameter'
       ENDIF
    ELSEIF (S .LE. 7 .and. S .gt. 1.0) THEN
       IF    (N.EQ.1) THEN
          HELREC=HYDREC(T,1)*((1.626E-11_dp+2.073E-11_dp)/&
               (1.499E-11_dp+1.847E-11_dp))
       ELSEIF(N.EQ.2) THEN
          HELREC=HYDREC(T,2)*((1.670E-11_dp+1.394E-11_dp)/&
               (1.374E-11_dp+1.103E-11_dp))
       ELSEIF(N.EQ.3) THEN
          HELREC=HYDREC(T,3)* 1.0_dp                              ! free-free ist the same
       ELSEIF(N .EQ. 4) THEN
          HELREC=HYDREC(T,4)*(2.073_dp/1.847_dp)
       ELSEIF(N .EQ. 5) THEN
          HELREC=HYDREC(T,5)*(1.369E-11_dp/1.103E-11_dp)
       ELSE
          STOP 'HELREC: Bad parameter'
       ENDIF
    ELSE
       STOP 'HELREC: Temperature out of range'
    ENDIF
!!$      iF (HELREC .lt. 1)&
!!$           stop "mess with helrec, neagtive rates !!!"
    RETURN
!-------------------------------------------------------------------------
  end function HELREC
!------------------------------------------------------------------------------

!================================================================================
  subroutine EXPAND_SPACE
!
!    Considers the expansion of space: increases l_cell, 
!    reduces particle densiteis.
!
! called by main program
!
    use M_definitions,        only: redshift, redshift_old, l_cell,&
                                    include_metals
    use M_grid_memory,        only: ne, nH_complete, nHI, nHII,nHe_complete,&
                                    nHeI, nHeII, nHeIII, nC, nN, nO, nNe, nS
! 
  implicit none
!
    real(dp) :: expansion
!===============================================================================
!    
    expansion=(redshift_old+1)/(redshift+1)
!    
    l_cell=l_cell*expansion
!   
    nH_complete=nH_complete/expansion**3
    ne=ne/expansion**3
    nHI=nHI/expansion**3
    nHII=nHII/expansion**3
    nHeI=nHeI/expansion**3
    nHeII=nHeII/expansion**3
    nHeIII=nHeIII/expansion**3
    nHe_complete=nHe_complete/expansion**3
!    
    if (include_metals)then
       nC=nC/expansion**3
       nN=nN/expansion**3
       nO=nO/expansion**3
       nNe=nNe/expansion**3
       nS=nS/expansion**3
    end if
!    
    redshift_old=redshift
  end subroutine  EXPAND_SPACE
  
end module M_hydrogen

