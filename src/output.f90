!===============================================================================
!                              ****Module M_output****
! This module writes the values of the nHI array (containing the density of 
! neutral hydrogen)  in a textfile readable by ifrit. Additionally, there is 
! a subroutine that writes a logfile. 
!===============================================================================
!
module M_output
!--- routines:
!    WRITE_IFRIT
!    WRITE_IFRIT2
!    CREATE_LOGFILE
!    PARTICLES

  use M_data_types
  
  implicit none
  contains
  
!===============================================================================
 subroutine WRITE_IFRIT
!
!    Writes the abundances of the considered ions, the electron density and 
!    the temperature in Textfiles that can be read by the IFRiT visualization
!    program
!
!  called by the main program
!
   use M_natural_constants, only: nc_light 
   use M_definitions,       only: rep_output, last_step, t_all, include_metals,&
                                  compute_temperature, log_heating_cooling,&
                                  extra_heating, extra_radiation_field,&
                                  o_thermal_pressure,&
                                  x_max, y_max, z_max, t_step,&
                                  o_write_escape_fraction, num_sources, points,&
                                  include_He, include_C, include_N, include_O,&
                                  include_Ne, include_S, include_Ar,&
                                  Heions, Cions, Nions, Oions, Neions,&
                                  Sions, Arions
   use M_data_input,        only: nug, nug_spect
   use M_grid_memory,       only: nHI, nHeI, nH_complete,&
                                  nHeI, nHeII, nHeIII, nHe_complete,&
                                  nC, nN, nO, nNe, nS,nAr,&
                                  temperature, temperature_old, ne,&
                                  heat_H, cool_H, heat_He,&
                                  energycontent, thermal_pressure,&
                                  cool_He, cool_C, cool_N, cool_O, cool_Ne,&
                                  cool_S, cool_Ar, cool_ff
   use M_raysave,           only: escapefraction 
   use M_hydrogen,          only: rep
! 
    implicit none
!
    integer(i4b)        :: file_counter=1000
    integer(i4b)        :: x, y, z, m=0, f, s
    character (len=3000):: string
    character (len=3000):: prefix_string
    character (len=4)   :: file_counter_string
!===============================================================================
!  
    if((mod(rep,rep_output) .eq. 0) .or. last_step .or.&
         extra_heating .or. extra_radiation_field )  then
       if (extra_heating .or. extra_radiation_field) then
          prefix_string='dump'
       else if (last_step) then
          prefix_string='final'
       else
          write(string,*) t_all/31557600._dp
          write(file_counter_string,'(I4)') file_counter
          prefix_string=adjustr(trim(adjustl(string)))//'_'//file_counter_string
       end if
!
!    re-entry point in the case both "dump" and "and the time-files  
100    continue 
       if (adjustr(trim(adjustl(prefix_string))) .ne. 'dump') then
          file_counter=file_counter+1
       end if
          
       write(file_counter_string,'(I4)') file_counter
       open(105,file=adjustr(trim(adjustl(prefix_string)))//'.H.txt')
       open(106,file=adjustr(trim(adjustl(prefix_string)))//'.He.txt')
       open(107,file=adjustr(trim(adjustl(prefix_string)))//'.e.txt')
       if (include_metals) then
          open(108,file=adjustr(trim(adjustl(prefix_string)))//'.O.txt')
          open(109,file=adjustr(trim(adjustl(prefix_string)))//'.N.txt')
          open(110,file=adjustr(trim(adjustl(prefix_string)))//'.C.txt')
          open(111,file=adjustr(trim(adjustl(prefix_string)))//'.Ne.txt')
          open(112,file=adjustr(trim(adjustl(prefix_string)))//'.S.txt')
          open(113,file=adjustr(trim(adjustl(prefix_string)))//'.Ar.txt')
       end if
!
       if (compute_temperature) then
          open(150,file=adjustr(trim(adjustl(prefix_string)))//'.temp.txt')
          open(151,file=adjustr(trim(adjustl(prefix_string)))//'.energy.txt')
       end if
!
       if (log_heating_cooling) then
          open(200, file=adjustr(trim(adjustl(prefix_string)))//'.hheatcool.txt')
          open(201, file=adjustr(trim(adjustl(prefix_string)))//'.heheatcool.txt')
          open(202, file=adjustr(trim(adjustl(prefix_string)))//'.cheatcool.txt')
          open(203, file=adjustr(trim(adjustl(prefix_string)))//'.nheatcool.txt')
          open(204, file=adjustr(trim(adjustl(prefix_string)))//'.oheatcool.txt')
          open(205, file=adjustr(trim(adjustl(prefix_string)))//'.neheatcool.txt')
          open(206, file=adjustr(trim(adjustl(prefix_string)))//'.sheatcool.txt')
          open(208, file=adjustr(trim(adjustl(prefix_string)))//'.arheatcool.txt')
! TODO change filenumber for free-free radiation
          open(207, file=adjustr(trim(adjustl(prefix_string)))//'.ffheatcool.txt')
       end if
!
       if (o_thermal_pressure) then
          open(300, file=adjustr(trim(adjustl(prefix_string)))//'.tp.txt')
       end if
!
       write(105,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       write(106,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       write(107,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
!
       if (include_metals) then
          write(108,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(109,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(110,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(111,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(112,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(113,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       end if
! 
       if (compute_temperature) then
          write(150,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(151,*)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       end if
!
       if (log_heating_cooling) then
          write(200,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(201,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(202,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(203,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(204,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(205,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(206,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(206,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(207,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(208,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
       end if
!
       if (o_thermal_pressure)&
            write(300,*) (2*x_max+1),(2*y_max+1), (2*z_max+1)
!     
       do z=-z_max,z_max
          do y=-y_max,y_max
             do x=-x_max,x_max
!
                write(105,'(3ES15.6E3)') nHI(x, y, z),&
                     nH_complete(x, y,z)-nHI(x, y, z),&
                     nH_complete(x, y, z)
                if (include_He)&
                     write(106,'(4ES15.6E3)') nHeI(x, y, z),&
                     nHeII(x, y, z),&
                     nHeIII(x, y, z),&
                     nHe_complete(x, y, z)
                write(107,*) ne(x, y, z)
!              
                if (include_metals) then
!
                   write(108,ifrit_formatstring(Oions)) nO(x, y, z,1:Oions),&
                        nO(x, y, z,0)
!                
                   write(109,ifrit_formatstring(Nions)) nN(x, y, z,1:Nions),&
                        nN(x, y, z,0)
!
                   write(110,ifrit_formatstring(Cions)) nC(x, y, z,1:Cions),&
                        nC(x, y, z,0)
!                 
                   write(111,ifrit_formatstring(Neions)) nNe(x, y, z,1:Neions),&
                        nNe(x, y, z,0)
!                 
                   write(112,ifrit_formatstring(Sions)) nS(x, y, z,1:Sions),&
                        nS(x, y, z,0)
!                  
                   write(113,ifrit_formatstring(Arions)) nAr(x, y, z,1:Arions),&
                        nAr(x, y, z,0)
                end if
!              
                if (compute_temperature .and. rep .ge. 2) then
                   write(150,'(3ES15.6E3)') temperature(x, y, z),&
                        temperature_old(x, y, z),&
                        ( temperature(x, y, z)-&
                        temperature_old(x, y, z))/&
                        t_step*365.25*24*3600
                   write(151,*) energycontent(x,y,z)
                end if
!
                if (log_heating_cooling) then
                   write(200,*) heat_H(x, y, z),&
                        cool_H(x, y, z),&
                        heat_H(x, y, z)-&
                        cool_H(x, y, z)
                   write(201,*) heat_He(x, y, z),&
                        cool_He(x, y, z),&
                        heat_He(x, y, z)-cool_He(x, y, z)
                   write(202,*) 0, cool_C(x, y, z),&
                        0- cool_C(x, y, z)
                   write(203,*) 0, cool_N(x, y, z),&
                        0- cool_N(y, y, z)
                   write(204,*) 0, cool_O(x, y, z),&
                        0- cool_O(x, y, z)
                   write(205,*) 0, cool_Ne(x, y, z),&
                        0- cool_Ne(x, y, z)
                   write(206,*) 0, cool_S(x, y, z),&
                        0- cool_S(x, y, z)    
                   write(207,*) 0, cool_ff(x, y, z),&
                        0- cool_ff(x, y, z)
                   write(208,*) 0, cool_Ar(x, y, z),&
                        0- cool_ff(x, y, z)
                end if

                if (o_thermal_pressure)&
                     write(300,*) thermal_pressure(x, y, z)   
             end do
        end do
     end do
     close(105)
     close(106)
     close(107)
     if (include_metals) then
        close(108)
        close(109)
        close(110)
        close(111)
        close(112)
        close(113)
     end if
     close(150)
     close(151)
     if (log_heating_cooling) then
        close(200)
        close(201)
        close(202)
        close(203)
        close(204)
        close(205)
        close(206)
        close(207)
        close(208)
     end if
!    
     if (o_thermal_pressure)&
          close(300)
!
     if (prefix_string .eq. 'dump') then       !write both dump files & "snapshots"
        if (last_step) then
           prefix_string='final'
           goto 100
        else if (mod(rep,rep_output) .eq. 0) then 
           write(string,*) t_all/31557600._dp
           write(file_counter_string,'(I4)') file_counter
           write(*,*) 'file-counter_string:',file_counter_string
           prefix_string=adjustr(trim(adjustl(string)))//'_'//file_counter_string
           goto 100
        endif
       end if
!
     call particles(prefix_string)     
!
     if (o_write_escape_fraction) then
        open(120,file=adjustr(trim(adjustl(prefix_string)))//'.escape.txt')
        do s=1, num_sources
           write(120,*) '#SOURCE #,',s,':'
           do f=1, points-1
              write(120,*) 1e+8*nc_light/nug(f), nug_spect(f), escapefraction (s,f)
           end do
           write(*,*)
           write(*,*)
        end do
        close(120)
     end if
  end if

contains

  character (len=11) function ifrit_formatstring (nion)
    implicit none
    integer(i4b)::  nion
    select case (nion)
       case (2)
          ifrit_formatstring = '(3ES15.6E3)'
       case (3)
          ifrit_formatstring = '(4ES15.6E3)'
       case (4)
          ifrit_formatstring = '(5ES15.6E3)'          
       case (5)
          ifrit_formatstring = '(6ES15.6E3)'
       case (6)
          ifrit_formatstring = '(7ES15.6E3)'          
    end select
  end function ifrit_formatstring
end subroutine WRITE_IFRIT

!===============================================================================
 subroutine WRITE_IFRIT_BIN
!
!    Writes the abundances of the considered ions, the electron density and 
!    the temperature in binary files that can be read by the IFRiT visualization
!    program
!
!  called by the main program
!
   use M_natural_constants, only: nc_light 
   use M_definitions,       only: rep_output, last_step, t_all, include_metals,&
                                  compute_temperature, log_heating_cooling,&
                                  binary_output,&
                                  extra_heating, extra_radiation_field,&
                                  o_thermal_pressure,&
                                  x_max, y_max, z_max, t_step,&
                                  o_write_escape_fraction, num_sources, points,&
                                  Heions, Cions, Nions, Oions, Neions,&
                                  Sions, Arions
   use M_data_input,        only: nug, nug_spect
   use M_grid_memory,       only: nHI, nHeI, nH_complete,&
                                  nHeI, nHeII, nHeIII, nHe_complete,&
                                  nC, nN, nO, nNe, nS,nAr,&
                                  temperature, temperature_old, ne,&
                                  heat_H, cool_H, heat_He,&
                                  energycontent, thermal_pressure,&
                                  cool_He, cool_C, cool_N, cool_O, cool_Ne,&
                                  cool_S, cool_Ar, cool_ff,&
                                  nHIf, nHeIf, nH_completef,&
                                  nHeIf, nHeIIf, nHeIIIf, nHe_completef,&
                                  nCf, nNf, nOf, nNef, nSf,nArf,&
                                  temperaturef, temperature_oldf, nef,&
                                  heat_Hf, cool_Hf, heat_Hef,&
                                  energycontentf, thermal_pressuref,&
                                  cool_Hef, cool_Cf, cool_Nf, cool_Of, cool_Nef,&
                                  cool_Sf, cool_Arf, cool_fff
   use M_raysave,           only: escapefraction 
   use M_hydrogen,          only: rep
! 
    implicit none
!
    integer(i4b)        :: file_counter=1000
    integer(i4b)        :: x, y, z, m=0, f, s
    character (len=3000):: string
    character (len=3000):: prefix_string
    character (len=4)   :: file_counter_string
!
!===============================================================================
!
    if((mod(rep,rep_output) .eq. 0) .or. last_step .or.&
         extra_heating .or. extra_radiation_field )  then
       call convert_to_float
       if ( prefix_string .ne. 'dump')&
            file_counter=file_counter+1       
       if (extra_heating .or. extra_radiation_field) then
          prefix_string='dump'
       else if (last_step) then
          prefix_string='final'
       else
          write(string,*) t_all/31557600._dp
          write(file_counter_string,'(I4)') file_counter
          write(*,*) 'file_counter_string:',file_counter_string
          prefix_string=adjustr(trim(adjustl(string)))//'_'//file_counter_string
       end if
!
!    re-entry point in the case both "dump" and "and the time-files  
100    continue                      
       write(file_counter_string,'(I4)') file_counter
       open(105,file=adjustr(trim(adjustl(prefix_string)))//'.H.bin',&
            form='unformatted')
       open(106,file=adjustr(trim(adjustl(prefix_string)))//'.He.bin',&
            form='unformatted')
       open(107,file=adjustr(trim(adjustl(prefix_string)))//'.e.bin',&
            form='unformatted')
       if (include_metals) then
          open(108,file=adjustr(trim(adjustl(prefix_string)))//'.O.bin',&
               form='unformatted')
          open(109,file=adjustr(trim(adjustl(prefix_string)))//'.N.bin',&
               form='unformatted')
          open(110,file=adjustr(trim(adjustl(prefix_string)))//'.C.bin',&
               form='unformatted')
          open(111,file=adjustr(trim(adjustl(prefix_string)))//'.Ne.bin',&
               form='unformatted')
          open(112,file=adjustr(trim(adjustl(prefix_string)))//'.S.bin',&
               form='unformatted')
          open(113,file=adjustr(trim(adjustl(prefix_string)))//'.Ar.bin',&
               form='unformatted')
       end if
!
       if (compute_temperature) then
          open(150,file=adjustr(trim(adjustl(prefix_string)))//'.temp.bin',&
               form='unformatted')
          open(151,file=adjustr(trim(adjustl(prefix_string)))//'.energy.bin',&
               form='unformatted')
       end if
!
       if (log_heating_cooling) then
          open(200, file=adjustr(trim(adjustl(prefix_string)))//'.hheatcool.bin',&
               form='unformatted')
          open(201, file=adjustr(trim(adjustl(prefix_string)))//'.heheatcool.bin',&
               form='unformatted')
          open(202, file=adjustr(trim(adjustl(prefix_string)))//'.cheatcool.bin',&
               form='unformatted')
          open(203, file=adjustr(trim(adjustl(prefix_string)))//'.nheatcool.bin',&
               form='unformatted')
          open(204, file=adjustr(trim(adjustl(prefix_string)))//'.oheatcool.bin',&
               form='unformatted')
          open(205, file=adjustr(trim(adjustl(prefix_string)))//'.neheatcool.bin',&
               form='unformatted')
          open(206, file=adjustr(trim(adjustl(prefix_string)))//'.sheatcool.bin',&
               form='unformatted')
          open(208, file=adjustr(trim(adjustl(prefix_string)))//'.arheatcool.bin',&
               form='unformatted')
! TODO change filenumber for free-free radiation
          open(207, file=adjustr(trim(adjustl(prefix_string)))//'.ffheatcool.bin',&
               form='unformatted')
       end if
!
       if (o_thermal_pressure) then
          open(300, file=adjustr(trim(adjustl(prefix_string)))//'.tp.bin',&
               form='unformatted')
       end if
!
       write(105)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       write(106)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       write(107)(2*x_max+1),(2*y_max+1),(2*z_max+1)
!
       if (include_metals) then
          write(108)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(109)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(110)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(111)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(112)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(113)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       end if
! 
       if (compute_temperature) then
          write(150)(2*x_max+1),(2*y_max+1),(2*z_max+1)
          write(151)(2*x_max+1),(2*y_max+1),(2*z_max+1)
       end if

       if (log_heating_cooling) then
          write(200) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(201) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(202) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(203) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(204) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(205) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(206) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(206) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(207) (2*x_max+1),(2*y_max+1), (2*z_max+1)
          write(208) (2*x_max+1),(2*y_max+1), (2*z_max+1)
       end if

       if (o_thermal_pressure)&
            write(300) (2*x_max+1),(2*y_max+1), (2*z_max+1)       
!
       write(105) nHIf
       write(105) nH_completef-nHIf
       write(105) nH_completef              
       write(106) nHeIf
       write(106) nHeIIf
       write(106) nHeIIIf
       write(106) nHe_completef
!
       write(107) nef
!              
       if (include_metals) then
!
          write(108) nOf(:,:,:,1:Oions),   nOf(:,:,:,0)
!                
          write(109) nNf(:,:,:,1:Nions),   nNf(:,:,:,0)
!
          write(110) nCf(:,:,:,1:Cions),   nCf(:,:,:,0)
!                 
          write(111) nNef(:,:,:,1:Neions), nNef(:,:,:,0)
!                       
          write(112) nSf(:,:,:,1:Sions),   nSf(:,:,:,0)
!
          write(113) nArf(:,:,:,1:Arions), nArf(:,:,:,0)
       end if
!              
       if (compute_temperature .and. rep .ge. 2) then
          write(150) temperaturef
          write(150) temperature_oldf
          write(150) (temperaturef-&
               temperature_oldf)/&
               t_step*365.25*24*3600
          write(151)  energycontentf
       end if
!
       if (log_heating_cooling) then
          write(200) heat_Hf
          write(200) cool_Hf
          write(200) heat_Hf-cool_Hf
!          
          write(201) heat_Hef
          write(201) cool_Hef
          write(201) heat_Hef-cool_Hef
!          
          write(202) 0._sp*cool_Cf
          write(202) cool_Cf
          write(202) 0._sp*cool_Cf- cool_Cf
!          
          write(203) 0._sp*cool_Nf
          write(203) cool_Nf
          write(203) 0._sp*cool_Nf- cool_Nf
!          
          write(204) 0._sp*cool_Of
          write(204) cool_Of
          write(204) 0._sp*cool_Of- cool_Of
!
          write(205) 0._sp*cool_Nef
          write(205) cool_Nef
          write(205) 0._sp*cool_Nef- cool_Nef
!          
          write(206) 0._sp*cool_Sf
          write(206) cool_Sf
          write(206) 0._sp*cool_Sf- cool_Sf    
!
          write(207) 0._sp*cool_fff
          write(207) cool_fff
          write(207) 0._sp*cool_fff- cool_fff
!          
          write(208) 0._sp*cool_Arf
          write(208) cool_Arf
          write(208) 0._sp*cool_Arf- cool_fff
       end if
!       
       if (o_thermal_pressure)&
            write(300) thermal_pressuref   
       close(105)
       close(106)
       close(107)
       if (include_metals) then
          close(108)
          close(109)
          close(110)
          close(111)
          close(112)
          close(113)
       end if
       close(150)
       close(151)
       if (log_heating_cooling) then
          close(200)
          close(201)
          close(202)
          close(203)
          close(204)
          close(205)
          close(206)
          close(207)
          close(208)
       end if
 !    
       if (o_thermal_pressure)&
            close(300)
!
       if (prefix_string .eq. 'dump') then       !write both dump files & "snapshots"
          if (last_step) then
             prefix_string='final'
             goto 100
          else if (mod(rep,rep_output) .eq. 0) then
             write(string,*) t_all/31557600._dp
             write(file_counter_string,'(I4)') file_counter
             write(*,*) 'file-counter_string:',file_counter_string
             prefix_string=adjustr(trim(adjustl(string)))//'_'//file_counter_string
             goto 100
          end if
       end if
       !
       call particles(prefix_string)     
       !
       if (o_write_escape_fraction) then
          open(120,file=adjustr(trim(adjustl(prefix_string)))//'.escape.txt')
          do s=1, num_sources
             write(120,*) '#SOURCE #,',s,':'
             do f=1, points-1
                write(120,*) 1e+8*nc_light/nug(f), nug_spect(f), escapefraction (s,f)
             end do
             write(*,*)
             write(*,*)
          end do
          close(120)
       end if
!
       m=m+1
!
    end if
!    
  contains
!
    subroutine convert_to_float
      implicit none
!      
      nef=real(ne,sp)
!      
      nHIf=real(nHI,sp)
      nH_completef=real(nH_complete,sp)
      nHeIf=real(nHeI,sp)
      nHeIIf=real(nHeII,sp)
      nHeIIIf=real(nHeIII,sp)
      nHe_completef=real(nHe_complete,sp)
      temperaturef=real(temperature,sp)
      temperature_oldf=real(temperature_old,sp)
      energycontentf=real(energycontent,sp)
      !               
      if (include_metals) then
         nCf=real(nC,sp)
         nNf=real(nN,sp)
         nOf=real(nO,sp)
         nNef=real(nNe,sp)
         nSf=real(nS,sp)
         nArf=real(nAr,sp)
      end if
!
      if (log_heating_cooling) then
         cool_Hf  = real(cool_H,sp)
         cool_Hef = real(cool_Hf,sp)
         cool_fff = real(cool_ff,sp)
         heat_Hf  = real(heat_H,sp)
         heat_Hef = real(heat_He,sp)
         if (include_metals) then
            cool_Cf=real(cool_C,sp)
            cool_Nf=real(cool_N,sp)
            cool_Of=real(cool_O,sp)
            cool_Nef=real(cool_Ne,sp)
            cool_Sf=real(cool_S,sp)
            cool_Arf=real(cool_Ar,sp)
         end if
      end if
!
      if(o_thermal_pressure) then
           thermal_pressuref=real(thermal_pressure,sp)
        end if
!     
    end subroutine convert_to_float
    
end subroutine WRITE_IFRIT_BIN

 
!================================================================================   
subroutine CREATE_LOGFILE
!
! Writes a logfile with the input parameters
! WARNING: This function is not up to date. As the parameters are already
! contained in the input file, this function may be removed in the future.
!
! called by the main program.
!
  use M_definitions,    only: x_max, y_max, z_max, x_min, y_min, z_min,&
                                t_all, num_sources, source_x, source_y,&
                                source_z, R_Rs, spectra, start_activity,&
                                lifetime, c_correct_on, f_a,l_cell, t_start,&
                                t_end, factor, read_from_file, nHI_fraction,&
                                nH_complete_scalar,cpu_start, cpu_finish,&
                                ifrit_filename
  use M_raysave,        only: l_max
  use M_hydrogen,       only: rep
  implicit none
!  
  integer(i4b) :: source_counter
!===============================================================================  
!  
  open(3, file='logfile.txt')
  write(3,*) 'Computation time in s          :', cpu_finish-cpu_start
  write(3,*) 'Number of repetitions          :', rep
  write(3,*) 'Elapsed simulated time in yrs. :', t_all, 'years'
  write(3,*) '***********************************************************'
  write(3,*) 'Number of rays', 12*4**l_max
  
  
  write(3,*)
  write(3,*) 'xrange :' , x_min , '...' , x_max
  write(3,*) 'yrange :' , y_min , '...' , y_max
  write(3,*) 'zrange :' , z_min , '...' , z_max
  write(3,*)
  write(3,*) 'Number of sources' , num_sources
  write(3,*)
!
  do source_counter=1, num_sources
     write(3,*) 'coordinates #',source_counter,' :', source_x(source_counter), &
          &source_y(source_counter), source_z(source_counter)     
     write(3,*) 'R_RS  (solar radii) #',source_counter,' :', R_RS(source_counter)
     write(3,*) 'spectra  #',source_counter,' :', spectra(source_counter)
     write(3,*) 'start_activity[yrs.] #',source_counter,' :',&
                start_activity(source_counter)
     write(3,*) 'lifetime[yrs.] #',source_counter,' :',lifetime(source_counter)
     write(3,*)
  end do
!  
  if (c_correct_on) then
     write(3,*) 'speed of light correction: on'
  else
     write(3,*) 'speed of light correction: off'
  end if
!  
  write(3,*)  
  write(3,*) 'spatial security factor :', f_a
  write(3,*)
  write(3,*) 'cell length in parsec :', l_cell
  write(3,*)
  write(3,*) 'timestep at the start of the program in s :',t_start
  write(3,*) 'total time of simulation in yrs. :', t_end
  write(3,*) 'temporal security factor :',factor
  write(3,*)
!  
  if (read_from_file) then
    write(3,*) 'ifrit file to be read ', ifrit_filename
  else
    write(3,*) 'initial density fraction of neutral hydrogen :',nHI_fraction
    write(3,*) 'total density of hydrogen :', nH_complete_scalar
  end if
!
  close(3)
  end subroutine CREATE_LOGFILE


  
!==============================================================================
  subroutine LOG_ITERATION_STEP
!
!   Prints the ionization fractions and the temperature (weighted by n_e^2
!   for each timestep, if the option verbose_stdout is set to .true.
!
!  called by the main program
!    
    use M_definitions,       only: include_metals,&
                                   Heions, Cions, Nions, Oions, Neions,&
                                   Sions, Arions
                                   
    use M_grid_memory,       only: nHI, nHII, nHeI, nH_complete,&
                                   nHeI, nHeII, nHeIII, nHe_complete,&
                                   nC, nN, nO, nNe, nS,nAr,&
                                   temperature, ne
!
    implicit none
    integer(i4b)    :: counter
!
!=============================================================================
    write(*,*)
!    
    if (max (Heions,Cions,Nions, Oions, Neions, Sions, Arions) .eq. 3)& 
         write(*,"(A10,A8,A8,A8)") "Element", "I","II", "III"    
    if (max (Heions,Cions,Nions, Oions, Neions, Sions, Arions) .eq. 4)& 
         write(*,"(A10,A8,A8,A8,A8)") "Element", "I","II", "III", "IV"
    if (max (Heions,Cions,Nions, Oions, Neions, Sions, Arions) .eq. 5)& 
         write(*,"(A10,A8,A8,A8,A8,A8)") "Element", "I","II", "III", "IV","V"
!
    write(*,log_formatstring(2))&
         "H",&    
         100._dp*sum(nHI )/sum(nH_complete),&
         100._dp*sum(nHII)/sum(nH_complete)
   write(*,log_formatstring(Heions))&
        "He",&    
        100_dp*sum(nHeI)/sum(nHe_complete),&
        100_dp*sum(nHeII)/sum(nHe_complete),&
        100_dp*sum(nHeIII)/sum(nHe_complete)
!
   if (include_metals) then
      write(*,log_formatstring(Cions))& 
           "C",&
           (100_dp*sum(nC(:,:,:,counter))/sum(nC(:,:,:,0)), counter=1,Cions)
      write(*,log_formatstring(Nions))&
           "N",&
           (100_dp*sum(nN(:,:,:,counter))/sum(nN(:,:,:,0)), counter=1,Nions)
      write(*,log_formatstring(Oions))&
           "O",&
           (100_dp*sum(nO(:,:,:,counter))/sum(nO(:,:,:,0)), counter=1,Oions)
      write(*,log_formatstring(Neions))&
           "Ne",&
           (100_dp*sum(nNe(:,:,:,counter))/sum(nNe(:,:,:,0)), counter=1,Neions)
      write(*,log_formatstring(Sions))&
           "S",&
           (100_dp*sum(nS(:,:,:,counter))/sum(nS(:,:,:,0)), counter=1,Sions)
      write(*,log_formatstring(Arions))&
           "Ar",&
           (100_dp*sum(nAr(:,:,:,counter))/sum(nAr(:,:,:,0)), counter=1,Arions)
   end if
!
   write(*,*)
   write(*,*) "Maximal gas Temperature in K:", maxval(temperature)
   write(*,*) "Minimal gas temperature in K:", minval(temperature)
   write(*,*) "Temperature of the gas, weighted by the square of the electron"
   write(*,*) "number density:", sum (ne**2*temperature)/sum(ne**2)
   write(*,*)
!
 contains
!
   character(len=11 ) function log_formatstring (nion)
     implicit none
     integer(i4b)     :: nion

     select case (nion)
     case (2)
        log_formatstring = '(A10,2F8.1)'
     case (3)
        log_formatstring = '(A10,3F8.1)'
     case (4)
        log_formatstring = '(A10,4F8.1)'
     case (5)
        log_formatstring = '(A10,5F8.1)'        
     end select
    
   end function log_formatstring
 end subroutine LOG_ITERATION_STEP


 
  
!===============================================================================
  subroutine PARTICLES(string)
!
!    The properties of the sources that are active at the current time in a file
!
!   called by WRITE_IFRIT!
!   
    use M_definitions,      only: x_max, y_max, z_max, x_min, y_min, z_min,&
                                  source_x, source_y, source_z,&
                                  num_sources, t_all, t_step, start_activity,&
                                  lifetime
!
    implicit none
!
    integer(i4b)               :: so, so_number
    character (len=3000)       :: string
!===============================================================================
!
    open(7,file='particles_'//adjustl(trim(string))//'.txt')
    
    so_number=0
    do so=1,num_sources
       if((t_all+t_step)/31557600._dp .ge. start_activity(so)&
         .and. (t_all+t_step)/31557600._dp .le. start_activity(so)+lifetime(so))&
            so_number=so_number+1 
    end do
!    
    write(7,*) so_number
    write(7,*) x_min, y_min, z_min, x_max, y_max, z_max
    do so=1,num_sources
       if((t_all+t_step)/31557600._dp .ge. start_activity(so)&
            .and. (t_all+t_step)/31557600._dp .le.&
            start_activity(so)+lifetime(so)) then
!         last number is attribute
          write(7,*) source_x(so), source_y(so), source_z(so), 1
       end if
    end do
!
    close(7)
  end subroutine PARTICLES

!================================================================================
  subroutine CALL_EXTERNAL_HC_PROGRAM
    use M_definitions,      only:  external_hc_program
    implicit none
!================================================================================
   call SYSTEM (external_hc_program) !Legacy, for example first generation flang
  ! call EXECUTE_COMMAND_LINE(external_hc_program) (corresponding to Fortran 2008)
  end subroutine CALL_EXTERNAL_HC_PROGRAM
!
end module M_output

!================================================================================
!================================================================================

!================================================================================
!                              ****Module M_output****
!================================================================================
