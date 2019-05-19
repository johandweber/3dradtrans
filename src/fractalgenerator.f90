program fractal
  implicit none
  real, dimension(:,:,:), allocatable :: density
  integer, dimension(:,:,:), allocatable :: numbers_in_grid
  
  real:: xpos,ypos,zpos
  
  real:: fractal_dimension
  
  integer::  x_max, y_max, z_max
  integer::  tiers, tier1, tier_counter
  integer:: x, y, z
  real:: emptyradius, smooth, clumpy
  integer:: generation_counter
  integer(8),save  :: number_of_points, void_points
  integer, dimension(20):: number_in_tier 
  real, dimension(20)::ranges
  real :: nop
  real:: ion_fraction
  
  write(*,*) 'Please intput x_max (xrange: [-x_max:+x_max],'//&
       ' i.e., 2*x_max+1 cells in x-direction) :'
  read(*,*)   x_max
  write(*,*)
  write(*,*) 'Please input y_max (yrange: [-y_max:+y_max),'//&
       ' i.e., 2*y_max+1 cells in y-direction) :'       
  read(*,*)   y_max
  write(*,*)
  write(*,*) 'Please input z_max (zrange: [-z_max, +z_max],'//&
       ' i.e., 2*z_max+1 cells in z-direction) :'              
  read(*,*) z_max
  allocate (density(-x_max:+x_max,-y_max:+y_max,-z_max:z_max))
  allocate (numbers_in_grid(-x_max:+x_max,-y_max:+y_max,-z_max:z_max))
  density=0
  numbers_in_grid=0
  number_of_points=0
  void_points=0
  write(*,*)
  write(*,*) 'Please input number of tiers :'
  read(*,*) tiers
  write(*,*)
  write(*,*) 'Please input number of tier1 values :'
  read(*,*) tier1
  write(*,*)
  write(*,*) 'Please input fractal dimension :'
  read(*,*) fractal_dimension 
  write(*,*)
  write(*,*) '"Smooth" density fracion in hydrogen atoms per cm^-3'
  read(*,*) smooth
  write(*,*)
  write(*,*) 'Mean density of "clumpy" density fraction in atoms per cm^-3'
  read(*,*) clumpy
  write(*,*)
  write(*,*) 'Please input initial ionization fraction'
  read(*,*) ion_fraction
  write(*,*)
  write(*,*) 'Please input radius of central hole (in grid cells):'
  read(*,*) emptyradius
  open(12, file="particles.txt")
  write(12,*) tier1
  write(12,*) 0,0,0,1,1,1
  do generation_counter=tiers,1,-1
    if (generation_counter .eq. tiers) then
      number_in_tier(generation_counter)= 64 !tier1
    else
      number_in_tier(generation_counter)=32
    end if  
  end do
  
  do generation_counter=5,1,-1
    if (generation_counter .eq. 5) then
      ranges(generation_counter)=1.
    else
      ranges(generation_counter)=ranges(generation_counter+1)/&
      (number_in_tier(generation_counter+1))**(1./fractal_dimension)
    end if
  end do
  
  write(*,*) 'generation','number_in_tier','ranges'
  
  do generation_counter = 5,1,-1
    write(*,*) generation_counter, number_in_tier(generation_counter),&
               ranges(generation_counter)
  end do
  write(*,*)
  do tier_counter=1, tier1 
    call random_number(xpos)
    call random_number(ypos)
    call random_number(zpos)
!    print*, 'tier_counter :', tier_counter
    print*, tier_counter,xpos,ypos,zpos
!    write(12,*) x,y,z,1
    call level(xpos,ypos,zpos, tiers,32,tier1)
  end do
    
  do x=-x_max,+x_max
     do y=-y_max,+y_max
        do z=-z_max,+z_max
           if (numbers_in_grid(x, y,z) .eq. 0) then
              void_points=void_points+1
           end if
        end do
     end do
  end do
  
  do x=-x_max,+x_max
    do y=-y_max,+y_max
       do z=-z_max,+z_max
          if (sqrt((1.*x)**2+(1.*y)**2+(1.*z)**2)&
              .lt. emptyradius) then
             numbers_in_grid(x,y, z)=0
          end if
       end do
    end do
 end do

  
  write(*,*) 'Number of void points :', void_points
  write(*,*) 'Percentage of void volume :',&
              real(void_points)/&
              real((2*x_max+1)*(2*y_max+1)*(2*z_max+1))*100., "%"

  nop=tier1
  do generation_counter=tiers-1,1,-1
    nop=nop*number_in_tier(generation_counter)
  end do
  print*, 'nop :', nop
    
!  nop=(2.*2./3.)*(2*x_max+1)*(2*y_max+1)*(2*z_max+1)/nop
    
  open(11,file='fractaltest.txt')
  write(11,*) 2*x_max+1, 2*y_max+1, 2*z_max+1
  do z=-x_max, +x_max
    do y=-y_max, +y_max
      do x=-x_max,+x_max
        write(11,*)& 
        smooth+numbers_in_grid(x,y,z)*clumpy*&
        (2*x_max+1)*(2*y_max+1)*(2*z_max+1)/nop,&
        smooth+numbers_in_grid(x,y,z)*clumpy*&
        (2*x_max+1)*(2*y_max+1)*(2*z_max+1)/nop
      end do
    end do
  end do
  close(11)
  close(12)
  write (*,*) 'Number of particles :', number_of_points
  write (*,*) 'Mean number of "particles" per grid cell :', 1.0*number_of_points/((2*x_max+1)*(2*y_max+1)*(2*z_max+1))
  stop 
  
contains

recursive subroutine level (x_0,y_0, z_0,generation,N,N_old)
  implicit none
  real:: x_0, y_0, z_0
  real:: x_local,y_local,z_local
  integer:: generation,N, N_counter,N_old
  
  if (generation .eq. 1) then
     if(.true.) then
!   if(x_0 .gt. 0. .and. x_0 .lt. 1. .and. y_0 .gt. 0. .and. y_0 .lt.1 .and. z_0 .gt.0 .and. z_0 .lt. 1) then

      numbers_in_grid(floor(-x_max+periodic(x_0,1.)*(2*x_max+0.99999)),&
                      floor(-y_max+periodic(y_0,1.)*(2*y_max+0.99999)),&
                      floor(-z_max+periodic(z_0,1.)*(2*z_max+0.99999)))=&
      numbers_in_grid(floor(-x_max+periodic(x_0,1.)*(2*x_max+0.99999)),&
                      floor(-y_max+periodic(y_0,1.)*(2*y_max+0.99999)),&
                      floor(-z_max+periodic(z_0,1.)*(2*z_max+0.99999)))+1
!      write(*,*) generation, N_counter, x_0, y_0, z_0
!      write(*,*)
      number_of_points=number_of_points+1
    end if
  else
    do N_counter=1,N
      call random_number(x_local)
      x_local=x_local-0.5
      call random_number(y_local)
      y_local=y_local-0.5
      call random_number(z_local)
      z_local=z_local-0.5
!       x_local=periodic(x_0 + x_local*ranges(generation),1.)
!       y_local=periodic(y_0 + y_local*ranges(generation),1.)
!       z_local=periodic(z_0 + z_local*ranges(generation),1.)
      x_local=x_0 + x_local*ranges(generation-1)
      y_local=y_0 + y_local*ranges(generation-1)
      
      z_local=z_0 + z_local*ranges(generation-1)
!      write(*,*) generation, N_counter, x_0, y_0, z_0
!      read(*,*)
!      print*, generation, N_counter,1.0*number_of_points/((1.*32)**tiers)*100,'%'
      call level(x_local,y_local,z_local,generation-1,32,32)
    end do
  end if
  return
end subroutine level 

 real function periodic (a,b)
   implicit none
   real:: a,b, tmp
   if (a.lt. -b .or. a .gt. 2*b) &
     stop 'range error'
   if (a .ge. 0. .and. a .le. b) then
     tmp=a
   else if (a .lt. 0) then
     tmp=b+a
   else if (a .gt. b) then
     tmp=a-b
   end if
   if (tmp .le.0 .or. tmp .ge. 2)&
     stop 'periodic out of range'  
   periodic=tmp 
 end function periodic


end program fractal


