program makelexington
! Creates the inittial conditions for some of the tests presened in 
! Pequignot et al. 2001  (2001ASPC..247..533P)
  use M_data_types
  use M_blackbody,                 only: Write_bb
  implicit none
!
! variables for batch mode
  logical       :: is_batch=.false.
  integer(i4b)  :: no_of_args
  integer(i4b)         :: mainmenuid, x_max_global
  character(len=20)   ::  selectionstring,  x_max_string
  no_of_args=command_argument_count()
  

  if (no_of_args .eq. 2) then
     is_batch=.true.
     call get_command_argument(1,selectionstring)
     read(selectionstring,*) mainmenuid
     call get_command_argument(2,x_max_string)
     read (x_max_string,*) x_max_global
  else
     write(*,*) 'Please choose which test to create'
     write(*,*)
     write(*,*) "HII40...........................................(1)"
     write(*,*) "HII20...........................................(2)"
     write(*,*) "PN150...........................................(3)"
     write(*,*) "PN75............................................(4)"
     write(*,*)
     write(*,*) "Abort...........................................(0)"
     read(*,*) mainmenuid
  end if
  write(*,*)
  write(*,*)
  select case(mainmenuid)
  case(1)
     call makehii40
  case(2)
     call makehii20
  case (3)
     call makepn150
  case(4)
     call makepn75
  end select

  stop
 
  contains
    
   subroutine makehii40
     implicit none
!
     integer(i4b)  :: x_max
! 
     if (is_batch) then
        x_max=x_max_global
        else
           write(*,*) "The resolution is(2*x_max+1)**3 grid cells, where the source"
           write(*,*) "is located in the middle of the volume. Please input x_max."
           read(*,*) x_max
           write(*,*)
        end if

     call write_abundances&
          ("0.00_1001.H.txt" , x_max, 0.97_dp, 5._dp, 100._dp, 2)
     call write_abundances&
          ("0.00_1001.He.txt", x_max, 0.97_dp, 5._dp, 10._dp, 3)
     call write_abundances&
          ("0.00_1001.C.txt", x_max, 0.97_dp, 5._dp , 0.022_dp, 4)
     call write_abundances&
          ("0.00_1001.N.txt", x_max, 0.97_dp, 5._dp , 0.005_dp, 4)
     call write_abundances&
          ("0.00_1001.O.txt", x_max, 0.97_dp, 5._dp , 0.033_dp, 4)
     call write_abundances&
          ("0.00_1001.Ne.txt", x_max, 0.97_dp, 5._dp, 0.005_dp, 4)
     call write_abundances&
          ("0.00_1001.S.txt", x_max, 0.97_dp, 5._dp, 0.0009_dp, 4)

! There is no Argon included in the HII-40 Lexington testcase, but occup_3d_ext
!  needs a file containing the density dructure of argon. Thus it is set to
!  a (hopefully) negligible value
     call write_abundances&
          ("0.00_1001.Ar.txt", x_max, 0.97_dp, 5._dp , 1e-10_dp, 4)
     call write_abundances&
          ("0.00_1001.e.txt", x_max, 0.97_dp, 5._dp , 0._dp, 1)

     call write_temperature ("0.00_1001.temp.txt", x_max, 7500._dp)
     call write_inputfile(x_max,5._dp)
     ! copy a mergespec (possibly for an effective temperature of 
     ! approx. 40kK to the workung directory and rename it MERGESPEC-40
     call write_bb(40000._dp, 19.6_dp,"TEMPLATE40kK", "LEXINGTON-HII40")
     call write_source("source.txt", "LEXINGTON-HII40",18.68_dp, 1E+9_dp)
   end subroutine makehii40
!
!
!   
   subroutine makehii20
     implicit none
!
     integer(i4b)  :: x_max
! 
     if (is_batch) then
        x_max=x_max_global
        else
           write(*,*) "The resolution is(2*x_max+1)**3 grid cells, where the source"
           write(*,*) "is located in the middle of the volume. Please input x_max."
           read(*,*) x_max
           write(*,*)
        end if

     call write_abundances&
          ("0.00_1001.H.txt" , x_max, 0.97_dp, 3._dp, 100._dp, 2)
     call write_abundances&
          ("0.00_1001.He.txt", x_max, 0.97_dp, 3._dp, 10._dp, 3)
     call write_abundances&
          ("0.00_1001.C.txt", x_max, 0.97_dp, 3._dp , 0.022_dp, 4)
     call write_abundances&
          ("0.00_1001.N.txt", x_max, 0.97_dp, 3._dp , 0.005_dp, 4)
     call write_abundances&
          ("0.00_1001.O.txt", x_max, 0.97_dp, 3._dp , 0.033_dp, 4)
     call write_abundances&
          ("0.00_1001.Ne.txt", x_max, 0.97_dp,3._dp, 0.005_dp, 4)
     call write_abundances&
          ("0.00_1001.S.txt", x_max, 0.97_dp, 3._dp, 0.0009_dp, 4)

! There is no Argon included in the HII-20 Lexington testcase, but occup_3d_ext
!  needs a file containing the density dructure of argon. Thus it is set to
!  a (hopefully) negligible value
     call write_abundances&
          ("0.00_1001.Ar.txt", x_max, 0.97_dp, 3._dp , 1e-10_dp, 4)
     call write_abundances&
          ("0.00_1001.e.txt", x_max, 0.97_dp, 3._dp , 0._dp, 1)

     call write_temperature ("0.00_1001.temp.txt", x_max, 7500._dp)
     call write_inputfile(x_max,3._dp)
     ! copy a mergespec (possibly for an effective temperature of 
     ! approx. 40kK to the workung directory and rename it MERGESPEC-40
     call write_bb(20000._dp, 19.6_dp,"TEMPLATE40kK", "LEXINGTON-HII20")
     call write_source("source.txt", "LEXINGTON-HII20",103.4_dp, 1E+9_dp) 
   end subroutine makehii20
   
   subroutine makepn150
     implicit none
     write(*,*) "Not implemented yet"
     stop
   end subroutine makepn150

   subroutine makepn75
     implicit none
     write(*,*) "Not implemented yet"
     stop
   end subroutine makepn75
   
   subroutine write_abundances(filename, x_max, inner, outer, abundance, nions)
     implicit none 
     character (len=*) :: filename
     integer(i4b)      :: x_max, nions
!    inner     : inner radius of the nebula in pc
!    outer     : distande between the source and the end intersection of the
!                coordinate axes with the border of the simulated volume
!    abundance : total abundance of an element in cm^-3 
     real(dp)          :: inner, outer, abundance
!
     integer(i4b)      :: xcounter, ycounter, zcounter, ioncounter
     real(dp)          :: l_cell, pseudoemptyfactor
!
     l_cell=(2*outer)/(2*x_max+1)
     open (unit=101, file=trim(filename))
     write(101,*) 2*x_max+1, 2*x_max+1, 2*x_max+1
     do zcounter=-x_max,+x_max      ! The volume is cubic, thus no z_max required
        do ycounter=-x_max, +x_max
           do xcounter=-x_max, +x_max
              if (l_cell*sqrt(1.0_dp*xcounter**2+ycounter**2+zcounter**2)&
                   .lt. inner) then
                 pseudoemptyfactor=1e-3_dp
              else
                 pseudoemptyfactor=1._dp
              end if
              write(101,'(ES16.9)',advance='no')&
                   pseudoemptyfactor*abundance
              do ioncounter=2,nions
                 write(101,'(ES16.9)',advance='no') 0._dp
              end do
              if (nions.gt.1) then
                 write(101,'(ES16.9)',advance='yes')&    ! line break
                      pseudoemptyfactor*abundance
              else
                 write(101,*)                            ! line break
              end if
           end do
        end do
     end do
     close (101)
   end subroutine write_abundances

   subroutine write_temperature(filename,x_max,temperature)
     implicit none
!
      character(len=*) :: filename
      integer(i4b)     :: x_max
      real(dp)         :: temperature
!
      integer(i4b)     :: cellcounter
!
     open(unit=101, file=trim(filename))
     write(101,*) 2*x_max+1,2*x_max+1, 2*x_max+1
     do cellcounter=1,(2*x_max+1)**3
        write(101,'(3ES16.9)') temperature, temperature, 0._dp 
     end do
     close(101)
   end subroutine write_temperature

   subroutine write_inputfile(x_max,outer)
     implicit none
!
     integer(i4b)  :: x_max
     real(dp)      :: outer
!     
     open(unit=101, file="input_values");
     write(101,*) '&input_values'
     write(101,*) 'c_correct_on=F,'
     write(101,*) 'charge_transfer=T,'
     write(101,*) 'compute_temperature=.true.,'
     write(101,*) 'Ctosolar=0.0,'     ! Metallicities are read via the
                                      ! "-resume" option 
     write(101,*) 'diffuse=F,'
     write(101,*) 'f_a=0.05,'
     write(101,*) 'factor=10.15,'
     write(101,*) 'filling_factor=1.0,'
     write(101,*) 'ifrit_filename="",'
     write(101,*) 'include_metals=T,'
     write(101,*) 'initial_redshift=100.,'
     write(101,*) 'l_cell=',(2*outer)/(2*x_max+1),','
     write(101,*) 'last_step=T,'
     write(101,*) 'mergespec=T,'
     write(101,*) 'Netosolar=0.0,'   
     write(101,*) 'nH_complete_scalar=0.0,'
     write(101,*) 'nHI_fraction=0.0,'
     write(101,*) 'Ntosolar=0.0,'
     write(101,*) 'num_threads=12,'
     write(101,*) 'Otosolar=0.0,'
     write(101,*) 'read_from_file=F,'
     write(101,*) 'rep_output=100,'
     write(101,*) 'source_filename="source.txt",'
     write(101,*) 'Stosolar=0.0,'
     write(101,*) 't_end=1.5e+4,'
     write(101,*) 't_start=1e+8,'
     write(101,*) 't_step_max=1e+10,'
     write(101,*) 'trace_expansion=F,'
     write(101,*) 'verbose_stdout=T,'
     write(101,*) 'x_max=',x_max,','
     write(101,*) 'y_max=',x_max,','
     write(101,*) 'z_max=',x_max,','
     write(101,*) '/'
     close(101)
   end subroutine write_inputfile
   
   subroutine write_source(sourcefilename, mergespecname,radius, tmax)
     implicit none
     character (len=*) :: sourcefilename, mergespecname
     real(dp)          :: radius, tmax
     open(unit=101, file=trim(sourcefilename))
     write(101,*) '#Number of sources'
     write(101,*) 1
     write(101,*)
     write(101,'(I4,I4,I4,ES15.6,A50, ES15.6,ES15.6)') 0,0,0,radius, trim(mergespecname),0., tmax
     write(101,*)
     close(101)
   end subroutine write_source

end program makelexington
