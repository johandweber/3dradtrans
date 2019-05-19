! This program cuts out a part of the volume 
! of an occup_3d_ext output file.
! Still the output files represent the same 
! dimensions. The rest of the volume is "filled"
! with an almost perfect vaccuum (n_ion=10^-30 cm^-3)
! for all considered ion species
! The intent of the residduum is to make sure that
! the data can be postprocessed by programs that can not
! handle "true" vacuum
!
! Mod jweber 2019-04

program restrict_volume
!
  use M_data_types
  use M_natural_constants
  use M_definitions
  implicit none
!
  character(len=1024) :: xcenterstring,&
                         ycenterstring,&
                         zcenterstring,&
                         radiusstring, &
                         filenamestring

  integer(i4b)::         no_of_args

  real(dp)::             xcenter,&
                         ycenter,&
                         zcenter,&
                         radius

  logical, dimension(:,:,:), allocatable :: restrict_mask

  no_of_args=command_argument_count()

  if (no_of_args .ne. 5) then
     call write_errormessage
  end if

  call get_command_argument (1, xcenterstring)
  call get_command_argument (2, ycenterstring)
  call get_command_argument (3, zcenterstring)
  call get_command_argument (4, radiusstring)
  call get_command_argument (5, filenamestring)

  read (xcenterstring,*) xcenter
  read (ycenterstring,*) ycenter
  read (zcenterstring,*) zcenter
  read (radiusstring, *) radius 

  read(*,input_values) 

  allocate(restrict_mask(-x_max:+x_max, -y_max:+y_max, -z_max:+z_max))

  call generate_mask(xcenter,ycenter, zcenter, radius)



  call restrict_file(trim(filenamestring)//".e.txt",1,&
                     trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                     trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                     trim(zcenterstring)//"_"//trim(radiusstring)//&
                     ".e.txt", 1E-30_dp*(1+2+sumto(Cions-1)+sumto(Nions-1)+&
                     sumto(Oions-1)+sumto(Neions-1)+sumto(Sions-1)+&
                     sumto(Arions-1)))

  if (compute_temperature) then
     call restrict_file(trim(filenamestring)//".temp.txt",1,&
                        trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                        trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                        trim(zcenterstring)//"_"//trim(radiusstring)//&
                        ".temp.txt", 5000._dp)
  end if

  call restrict_file(trim(filenamestring)//".H.txt",3,&
                     trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                     trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                     trim(zcenterstring)//"_"//trim(radiusstring)//&
                     ".H.txt", 1E-30_dp)

  call restrict_file(trim(filenamestring)//".He.txt",4,&
                     trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                     trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                     trim(zcenterstring)//"_"//trim(radiusstring)//&
                     ".He.txt", 1E-30_dp)
  
  
  if (include_metals) then
     call restrict_file(trim(filenamestring)//".C.txt",Cions+1,&
                        trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                        trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                        trim(zcenterstring)//"_"//trim(radiusstring)//&
                        ".C.txt", 1E-30_dp)

     call restrict_file(trim(filenamestring)//".N.txt",Nions+1,&
                        trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                        trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                        trim(zcenterstring)//"_"//trim(radiusstring)//&
                        ".N.txt", 1E-30_dp)


     call restrict_file(trim(filenamestring)//".O.txt",Oions+1,&
                        trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                        trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                        trim(zcenterstring)//"_"//trim(radiusstring)//&
                        ".O.txt", 1E-30_dp)

     call restrict_file(trim(filenamestring)//".Ne.txt",Neions+1,&
                        trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                        trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                        trim(zcenterstring)//"_"//trim(radiusstring)//&
                        ".Ne.txt", 1E-30_dp)

     call restrict_file(trim(filenamestring)//".S.txt",Sions+1,&
                        trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                        trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                        trim(zcenterstring)//"_"//trim(radiusstring)//&
                        ".S.txt", 1E-30_dp)


     call restrict_file(trim(filenamestring)//".Ar.txt",Arions+1,&
                        trim(filenamestring)//"_"//trim(xcenterstring)//"_"//&
                        trim(xcenterstring)//"_"//trim(ycenterstring)//"_"//&
                        trim(zcenterstring)//"_"//trim(radiusstring)//&
                       ".Ar.txt", 1E-30_dp)
  end if

  stop

contains
  
  subroutine write_errormessage
    implicit none
     write(*,*)
     write(*,*) "Wrong number of command line arguments"
     write(*,*)
     write(*,*) "Correct usage:"
     write(*,*) &          
          "restrict_volume xcent ycent zcent rad name_without_suff < input_file"
     write(*,*)
     write(*,*) "xcent: x-position of the center of the included sphere in pc in"
     write(*,*) "a coordinate system where the coordinates of the center of"
     write(*,*) "the simulated volume are (0|0|0). The length of an edge of a"
     write(*,*) "simulation cell in pc is provided by the variable l_cell"
     write(*,*) "that is read from the name list contained in input_file" 
     write(*,*) "(see below)."
     write(*,*)
     write(*,*) "ycent: y-position of the center of the included sphere in pc"
     write(*,*) "in the coordinate system described above."
     write(*,*)
     write(*,*) "zcent: z-position of the center of the included sphere in pc"
     write(*,*) "in the coordinate system described above."
     write(*,*) 
     write(*,*) "rad: radius of the included sphere in pc"
     write(*,*)
     write(*,*) "For xcent, ycent, zcent, and rad a floating-pint value can"
     write(*,*) "be used."
     write(*,*)
     write(*,*) "name_without_suffix: Corresponds to the names of the"
     write(*,*) "files containing the output data, which are modified,"
     write(*,*) "without the suffix"
     write(*,*) "For example, if the simaulation results are written in"
     write(*,*) "files with the names 12.356_1007.xxxx.txt, where xxxx specifies"
     write(*,*) "the type of the output data (e.g, the ion number densities of"
     write(*,*) "the considered elements, the electron number densities, the" 
     write(*,*) "temperature ect.), name_without suff is 12.356_1007."
     write(*,*) "The output files are then called"
     write(*,*) "123456_1007.xcent_ycent_zcent_rad.xxxx.txt ."
     write(*,*)
     write(*,*) "input_file: name of the file that contains the parameters of"
     write(*,*) "the simulation, whose output data are used, in a"
     write(*,*) "Fortran namelist."
     write(*,*) "To avoid wrong results, make sure to use the dase file that"
     write(*,*) "has been used by the corresponding run of occup_3d_ext."
     write(*,*)
     stop
  end subroutine write_errormessage

  subroutine restrict_file (inputfilename, no_of_columns, outputfilename,&
                            replacement_value)
    implicit none
    integer(i4b)                               :: no_of_columns
    integer(i4b)                               :: xrange, yrange, zrange
    integer(i4b)                               :: x,y,z, colcounter
    real(dp)                                   :: replacement_value
    real(dp), dimension (1:1024), save         :: columns 
    character(len=*)                           :: inputfilename, outputfilename
    open(unit=101, file=inputfilename)
    open (unit=102,file=outputfilename)
    read(101,*)  xrange, yrange, zrange
    write(102,*) xrange, yrange, zrange
    do z=-zrange/2,+zrange/2
       do y=-yrange/2,+yrange/2
          do x=-xrange/2,+xrange/2
             read(101,*) columns(1:no_of_columns)
             if (restrict_mask(x,y,z)) then
                do colcounter=1, no_of_columns-1
                   write(102,'(ES16.7E3)', advance='no') columns(colcounter)
                end do
                write(102,'(ES16.7E3)') columns(no_of_columns) 
             else
                do colcounter=1, no_of_columns-1
                   write(102,'(ES16.7E3)', advance='no') replacement_value
                end do                
                write(102,'(ES16.7E3)') replacement_value
             end if
          end do
       end do
    end do
    close(102)
    close(101)
  end subroutine restrict_file


  subroutine generate_mask(xcenter, ycenter, zcenter, radius)
    implicit none
    real(dp) :: xcenter, ycenter, zcenter, radius
    integer  :: x, y, z
    do z=-z_max, +z_max
       do y=-y_max, +y_max
          do x=-x_max, +x_max
             if (sqrt(&
                  (l_cell*x-xcenter)**2+&       ! here we rely on automatic 
                  (l_cell*y-ycenter)**2+&       ! type conversion
                  (l_cell*z-zcenter)**2)&
                 .le. radius) then
                restrict_mask(x,y,z)=.true.
             else
                restrict_mask(x,y,z)=.false.
             end if
          end do
       end do
    end do
  end subroutine generate_mask


! The function sumto computes the sum
! 1,..,maximal value. It is required to make the number
! density of the electors in the pseudo-"vacuum" consistent
! with the number densities of the ions

  integer(i4b) function sumto(maximal_value)
    implicit none
    integer(i4b)   :: maximal_value
    sumto=(maximal_value)*(maximal_value+1)/2
  end function sumto

end program restrict_volume
