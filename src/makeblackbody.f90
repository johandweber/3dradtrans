!=============================================================================
program MAKEBLACKBODY
! creates a fake "MERGESPEC" file that can be read from occup_3d_ext
!
  use M_data_types
  use M_blackbody
!
  implicit none
!
  integer            :: no_of_args
  integer            :: x_index, y_index, z_index, num_lines
  integer            :: linecounter
  character(len=255) :: temperature_s,radius_s, templatefile_s, outputfile_s
  character(len=255) :: inputfile_s
  character(len=255) :: time_filename, mergespec_filename
  real(dp)           :: lifetime, temperature, radius
  real(dp)           :: lifetime_old=0.0_dp

!  
!=============================================================================
!
  no_of_args=command_argument_count()
  if ((no_of_args .ne. 2) .and. (no_of_args .ne. 4)) then
     write(*,*) 'wrong number of arguments'
     write(*,*)
     write(*,*) 'correct usage:'
     write(*,*)
     write(*,*)&
        'makeblackbody <radius/R_sun> <temperature/K> <templatefile> <outputfile>'
     write(*,*)
     write(*,*) 'where <templatefile> is a MERGESPEC output file of WM-Basic that'
     write(*,*) 'defines thze wavenlength grid'
     write(*,*)
     write(*,*) 'or'
     write(*,*)
     write(*,*) 'makeblackbody <inputfile> <templatefile>'
     write(*,*) 'where the structure of the input file is as follows:'
     write(*,*)
     write(*,*) '# Header: first line, can contain  comments'
     write(*,*) '<x_index> <y_index> <z_index> <no_of_timepoints>'
     write(*,*) 'time_1 radius_1 t_eff_1'
     write(*,*) 'time_2 radius_2 t_eff_2'
     write(*,*) '...'
     write(*,*) 'time_not radius_not t_eff_not'
     write(*,*)
     write(*,*) 'x_index, y_index and z_index are the grid coordinates'
     write(*,*) 'of the source'
     stop
!
  elseif (no_of_args .eq. 2) then
     call get_command_argument(1, inputfile_s)
     call get_command_argument(2, templatefile_s)
     write(*,*) inputfile_s
     open(unit=11, file=inputfile_s)
     open(unit=21, file='sources.txt')
     write(*,*) inputfile_s
     read(11,*)                                   ! line with comments
     read(11,*)  x_index, y_index, z_index, num_lines
     write(21,*) "#", "file generated automatically from ",trim(inputfile_s)
     write(21,*) num_lines
     do linecounter=1, num_lines
        read(11,*) lifetime, radius, temperature
        write(time_filename,*) lifetime_old
        time_filename=adjustr(trim(adjustl(time_filename)))//'.MERGESPEC'
        write(21,*)
        write(21,'(I4,I4,I4,E25.10, 3X, A50,E20.10, E20.10)') x_index, y_index,&
             z_index, radius, time_filename, lifetime_old, lifetime-lifetime_old
        call write_bb(temperature, radius, templatefile_s,time_filename)
        lifetime_old=lifetime
     end do
     close(11)
     close(21)
 !    
  elseif (no_of_args .eq. 4) then
     call get_command_argument(1, radius_s)
     call get_command_argument(2, temperature_s)
     call get_command_argument(3, templatefile_s)
     call get_command_argument(4, outputfile_s)
     read(temperature_s,*) temperature
     read(radius_s,*) radius
     call write_bb(temperature, radius, templatefile_s, outputfile_s)
     stop
  endif

end program MAKEBLACKBODY
