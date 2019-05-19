!==============================================================================
program sb2ms
!
! converts the output spectra of Starburst99 into the MERGESPEC format
!
  use M_data_types
  use M_natural_constants          ,only: nc_light, nc_pi, nc_sol_rad
!
  implicit none
!
  integer(i4b)                         :: no_of_args
  integer(i4b)                         :: line_counter=0, lines, dummy_counter
  character (len=255)                  :: infile_name, outfile_name
  real(dp), dimension(:), allocatable  :: age,wavelength,&
                                          total,stellar,nebular,&
                                          total_nu,stellar_nu,nebular_nu
  integer(i4b)                         :: ntz_m=41
!==============================================================================

  no_of_args = command_argument_count() 
  if (no_of_args .ne. 2) then
     write(*,*) "incorrect number of arguments"
     write(*,*)
     write(*,*) "correct usage of the program:"
     write(*,*) "sb2ms starburst_file output_file"
     write(*,*)
     stop
  end if
  call get_command_argument (1,infile_name)
  open(1000, file=infile_name)!infile_name)
  do while (.true.)
     read(1000,*, end=100)
     line_counter=line_counter+1
  end do
100 rewind(1000)
  lines=line_counter
  allocate(age(1:lines))
  allocate(wavelength(1:lines))
  allocate(total(1:lines))
  allocate(stellar(1:lines))
  allocate(nebular(1:lines))
  allocate(total_nu(1:lines))
  allocate(stellar_nu(1:lines))
  allocate(nebular_nu(1:lines))
  do line_counter=1, lines
     read(1000,*) age(line_counter),&
                  wavelength(line_counter),&
                  total(line_counter),&
                  stellar(line_counter),&
                  nebular(line_counter)
!    years in s
     age(line_counter)=age(line_counter)*24._dp*3600_dp*365._dp
!    wavelengths in cm
     wavelength(line_counter)=wavelength(line_counter)*1E-8
!    The energyies are written as log10s     
!
     total(line_counter)  = 1E+8*10._dp**total(line_counter)
     stellar(line_counter)= 1E+8*10._dp**stellar(line_counter)
     nebular(line_counter)= 1e+8*10._dp**nebular(line_counter)
!
     total_nu(line_counter)   = wavelength(line_counter)**2/nc_light*&
                                total(line_counter)/(16*nc_pi**2*nc_sol_rad**2)
     stellar_nu(line_counter) = wavelength(line_counter)**2/nc_light*&
                                stellar(line_counter)/(16*nc_pi**2*nc_sol_rad**2)
     nebular_nu(line_counter) = wavelength(line_counter)**2/nc_light*&
                                nebular(line_counter)/(16*nc_pi**2*nc_sol_rad**2)
  end do
  close(1000)
!
  call get_command_argument(2, outfile_name)
  open(unit=1001, file=outfile_name)
  write(1001,"('  1 ',I6,ES16.8,2I6,3ES16.8//)")  lines, 0._dp,0,0,&
                                                  0._dp, 1.0_dp, 1E+4_dp
  do line_counter=lines,1,-1
     write(1001,"(4ES16.8)") wavelength(line_counter)*1E+8_dp,&
                             total_nu(line_counter),&
                             total_nu(line_counter), 1.0_dp
  end do
  write(1001,"(//'  1 ',I6)")ntz_m
  do dummy_counter=1,ntz_m
     write(1001,"(3ES16.8)")  1._dp, 1._dp
  end do
  close(1001)
end program sb2ms
