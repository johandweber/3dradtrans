program faketemp
  implicit none
  integer:: x,xmax,y,ymax, z, zmax
  real:: temp
  read(*,*) xmax, ymax, zmax, temp
  xmax=2*xmax+1
  ymax=2*ymax+1
  zmax=2*zmax+1
  write(*,*) xmax, ymax, zmax
  do z=1, zmax
     do y=1, ymax
        do x=1, xmax
           write(*,*) temp, temp, 0.0
        end do
     end do
  end do
end program faketemp
 
    
