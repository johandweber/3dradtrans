	program cutzeq0
	  use iso_fortran_env
	  implicit none
	  integer:: x,y,z,x_max, y_max, z_max
	  character(len=255) :: totalline
	  read(*,*) x_max, y_max, z_max
	  write(*,*) x_max, y_max, 1
	  do z=-z_max/2,+z_max/2
	     do y=-y_max/2,+y_max/2
	       do x=-x_max/2,+x_max/2
	        write(ERROR_UNIT,*) x,y,z
	         read(*,'(A)') totalline
	        if (z .eq. 0)&
	          write(*,'(A)') totalline
              end do
           end do
        end do
	end program cutzeq0
