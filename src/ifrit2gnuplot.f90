program IFRIT2GNUPLOT
!
! Converts infrit text data into 2D data readable by
! gnuplot. 3D values will be added up along directions that can 
! be specified in the second command line argument.
 
  use iso_fortran_env
  implicit none
  integer:: xrange, yrange, zrange
  integer:: xcounter, ycounter, zcounter
  integer:: no_of_command_lines, no_cols
  integer:: columncounter
  character(len=2) :: direction
  character(len=2) :: no_cols_string
  double precision, dimension (:,:,:,:), allocatable :: ifrit_values
  integer:: written_elements=0
!
  no_of_command_lines=command_argument_count()
!
  call get_command_argument(1,direction)
  call get_command_argument(2,no_cols_string)
!
  if (no_of_command_lines .ne. 2 .or.&
      ((direction .ne. "-x") .and.(direction .ne. "-y") .and.&
      (direction .ne. "-z"))) then
     write(ERROR_UNIT,*) "Number of command-line arguments", no_of_command_lines
     write(ERROR_UNIT,*) "correct usage:"
     write(ERROR_UNIT,*) "cat <ifritfile> ifrit2gnuplot -x|-y|-z <no_of_columns> > gnuplotfile"
     write(ERROR_UNIT,*)
     write(ERROR_UNIT,*) "For more details see the manual page"
     write(ERROR_UNIT,*)
     stop
  end if
 ! 
  read(no_cols_string,*) no_cols
!
  read(*,*) xrange, yrange, zrange
!
  allocate (ifrit_values(xrange, yrange, zrange, no_cols)) 
  do zcounter=1,zrange
     do ycounter=1, yrange
        do xcounter=1, xrange
           read(*,*) (ifrit_values(xcounter, ycounter, zcounter,columncounter),columncounter=1, no_cols)
        end do
     end do
  end do
!
  if (direction .eq. "-x") then
     do columncounter=1, no_cols
        do zcounter=1, zrange
           do ycounter=1, yrange
              write(*,'(ES18.8)', advance='no') sum(ifrit_values(:, ycounter, zcounter, columncounter))
           end do
           write(*,*)
        end do
        write(*,*)
        write(*,*)
     end do
  end if
!
  if (direction .eq. "-y") then
     do columncounter=1, no_cols
        do zcounter=1, zrange
           do xcounter=1, xrange
              write(*,'(ES18.8)', advance='no') sum(ifrit_values(xcounter, :, zcounter, columncounter))
           end do
           write(*,*)
        end do
        write(*,*)
        write(*,*)
     end do
  end if
!
  if (direction .eq. "-z") then
     do columncounter=1, no_cols
        do ycounter=1, yrange
           do xcounter=1, xrange
              write(*,'(ES18.8)', advance='no') sum(ifrit_values(xcounter, ycounter,:, columncounter))
              written_elements=written_elements+1
           end do
           write(*,*)
        end do
        write(*,*)
        write(*,*)
     end do
  end if


end program ifrit2gnuplot
