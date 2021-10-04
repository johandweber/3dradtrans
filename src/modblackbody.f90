!=============================================================================
!                            MODULE M_blackbody-mod
!=============================================================================
!
! Module that contains subroutines to deal with the properties of
! blackbody radiatiors
! MOD jweber 17-12
!-----------------------------------------------------------------------------
!
!
module M_blackbody
!--- routines
!    write_bb
!
  use M_data_types
  use M_natural_constants,           only: nc_light, nc_planck, nc_boltzmann
!
  implicit none
!
contains

!==============================================================================  
   subroutine WRITE_BB(temperature, radius, templatefile, outputfile)
!    writes the a blackbody spectrum into a file that resembles a
!    "MERGESPEC" file of WM-Basic that can be read by occup_3d_ext
!    if the "mergespec" switch in the input file is set to ".true.".
     implicit none

     real(dp)           :: temperature  ! temperature of the bb in K
     real(dp)           :: radius       ! radius in solar radii
     character(len=*)   :: templatefile !existing "MERGESPEC" file, from
                                        ! which the frequengy grid is read
     character(len=*)   :: outputfile   ! output file in the "MERGESPEC" format
      
     integer            :: start_int, dummy_int, lf_m, ifred_m, ifblue_m
     real(dp)           :: corr_m, hflux_m, rstar_m, tefflum_m
     real(dp)           :: wavelength, hedd, conthedd, heddratio
     integer(i4b)       :: wavelengthcounter, ntcounter
     real(dp)           :: frequency, hnue
!==============================================================================     
!
     open(file=trim(templatefile), unit=101)
     open(file=trim(outputfile), unit=102)
     read(101,*) dummy_int, lf_m, corr_m, ifred_m, ifblue_m, hflux_m,&
            rstar_m, tefflum_m
     write(102,'(I7,I7,ES15.6,I7,I7,ES15.6,ES15.6,ES15.6)')&
          dummy_int, lf_m, corr_m, ifred_m,&
          ifblue_m, hflux_m,  radius, temperature
     write(102,*)
     write(102,*)
     do wavelengthcounter=1, lf_m
        read(101,*) wavelength , hedd, conthedd, heddratio
        frequency=nc_light/(wavelength*1E-8)                  !1A =1E-8 cm
        hnue=(1._dp/4._dp) * (2*nc_planck*frequency**3)/(nc_light**2)*&
             (1._dp/(exp((nc_planck*frequency)/(nc_boltzmann*temperature))&
             -1._dp))
        write(102,'(ES15.6,ES15.6,ES15.6,ES15.6)') wavelength, hnue, hnue, 1._dp
     end do
     write(102,*)
     write(102,*)
     write(102,*) 1, 41
! Create dummy values for rosseland opacities and flux conservation,
! which do not make any physical sense, but are not required for the
! occup_3d_ext program
     do ntcounter=1, 41                           
        write(102,*) 1e+6,1 
     end do
     close(101)
     close(102)
   end subroutine WRITE_BB
    
end module M_blackbody
