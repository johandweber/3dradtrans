program threedhydro
  use M_data_types
  use M_natural_constants
!
  implicit none
!
  character (len=1024)     :: mode, full_line, ifrit_filename="STWITHD3DR"
  character (len=1)        :: plane
  real(dp)                 :: nhyd=10._dp, nhyd_fill=1E-2_dp,&
                              nHI_fraction=0.9999_dp,empty_radius,&
                              frac_dimension, macro_filling_factor=1,&
                              dens_scale,&
                              prefact, r_0, z_0, radius, &
                              x_half_axis, y_half_axis, z_half_axis, cutoff,&
                              dens_exponent, sigma
  integer(i4b)             :: x_max=0, y_max=0, z_max=0,&
                              x_center=0, y_center=0, z_center=0,&
                              frac_tiers, frac_tier_1
  logical                  :: random_mode=.false., mode_set=.false.,&
                              ensure_minimum=.true.
  real(dp), dimension(:,:,:), allocatable  :: grid
!
! just counter variables
  integer (i4b)            :: x,y,z

!
  namelist/hydro_params/&
       mode, nhyd, nhyd_fill, nHI_fraction,&
       frac_tiers,frac_dimension, frac_tier_1,&
       empty_radius, macro_filling_factor,&
       dens_exponent, dens_scale,&
       prefact, r_0, z_0,sigma,&
       radius, x_half_axis, y_half_axis,z_half_axis,&
       plane,&
       x_max, y_max, z_max, random_mode,&
       x_center, y_center, z_center,&
       ensure_minimum,&
       ifrit_filename
!  
  read(*,'(A)') full_line
  read(*,*) x_max, y_max, z_max
  allocate(grid(-x_max:+x_max, -y_max:+y_max, -z_max:+z_max))
  grid=nhyd_fill
!
  if (index (full_line, "HOMOGENEOUS") .ne. 0) then
     read(full_line,*) mode, nhyd
     call homogen
  end if
!  
  if (index (full_line, "FRACTAL") .ne. 0) then
     read(full_line,*) mode
     call fractal
  end if
!
  if (index (full_line, "EXPONENTIALRAD") .ne. 0) then
     read(full_line,*) mode, x_center, y_center, z_center, nhyd, dens_scale
     call exponentialrad
  end if
!
  if (index (full_line, "EXPONENTIALPLANE") .ne. 0) then
     read(full_line,*) mode, plane, x_center, y_center, z_center
     call exponentialplane
  end if
!
  if (index (full_line, "POWERRADCUTOFF") .ne. 0) then
     read(full_line,*) mode, x_center, y_center, z_center, nhyd, cutoff, r_0,&
                       dens_exponent
     call powerradcutoff
  end if
!
  if (index (full_line, "POWERPLANECUTOFF") .ne. 0) then
     read(full_line,*) mode, plane, x_center, y_center, z_center, nhyd, cutoff,&
                       z_0, dens_exponent
     call powerplanecutoff
    end if
  !
  if (index (full_line, "SPHERE") .ne. 0) then
     read(full_line,*) mode, x_center, y_center, z_center, nhyd, radius
     call sphere
  end if
!
  if (index (full_line, "ELLIPSOID") .ne. 0) then
     read(full_line,*) mode, x_center, y_center, z_center, nhyd,&
                       x_half_axis, y_half_axis, z_half_axis
     call  ellipsoid
  end if
!
  if (index (full_line, "GAUSSIAN") .ne. 0) then
     read(full_line,*) mode, x_center, y_center, z_center, nhyd,&
                       nhyd, sigma
     call  gaussian
  end if
!
  if (.not. mode_set)&
       stop 'STOP: DENSITY DISTRIBUTION NOT SPECIFIED'

!
  if (random_mode) then
     call random_remove
  end if

  write(*,'(3I6)') (2*x_max+1), (2*y_max+1), (2*z_max+1)
  do z=-z_max,+z_max
     do y=-y_max,+y_max
        do x=-x_max,+x_max
           if (ensure_minimum) then
              if (grid(x,y,z) .lt. nhyd_fill)&
                   grid(x,y,z)=nhyd_fill
           end if
           write(*,'(2ES20.10)') grid(x,y,z), grid(x,y,z)*nHI_fraction
        end do
     end do
  end do
  stop
!-----------------------------------------------------------------------------

  contains

    subroutine homogen
      implicit none
      grid=nhyd
      mode_set=.true.
    end subroutine homogen
  
    subroutine fractal
      implicit none
      ! source code still to be added
      mode_set=.true.
    end subroutine fractal
  
    subroutine exponentialrad
      implicit none
      integer(i4b):: x, y, z
      real   ( dp):: distance 
      do z=-z_max, +z_max
         do y=-y_max, +y_max
            do x=-x_max, +x_max
               distance=sqrt(1._dp*(x-x_center)**2+(y-y_center)**2+&
                                    (z-z_center)**2)
               grid(x,y,z)=nhyd*exp(-distance/dens_scale)
            end do
         end do
      end do
      mode_set=.true.
    end subroutine exponentialrad

    subroutine exponentialplane
      implicit none
      integer(i4b):: x,y,z
      real   ( dp):: distance 
         do z=-z_max, +z_max
            do y=-y_max, +y_max
               do x=-x_max, +x_max
                  if (plane .eq. 'x' .or. plane .eq. 'X')&
                       distance = abs(x-x_center)
                  if (plane .eq. 'y' .or. plane .eq. 'Y')&
                       distance = abs(y-y_center)
                  if (plane .eq. 'z' .or. plane .eq. 'Z')&
                       distance = abs(z-z_center)
                  grid(x,y,z)=nhyd*exp(-distance/dens_scale)                 
               end do
            end do
         end do
         mode_set=.true.
    end subroutine exponentialplane

    subroutine powerradcutoff
      implicit none
      integer(i4b):: x,y,z
      real   ( dp):: distance
      do z=-z_max, +z_max
         do y=-y_max, +y_max
            do x=-x_max, +x_max
               distance=sqrt(1._dp*(x-x_center)**2+(y-y_center)**2+&
                                   (z-z_center)**2)    
               if (x .eq. x_center .and. y .eq. y_center .and. &
                   z .eq. z_center) then
                  grid(x,y,z)= nhyd 
               else
                  if (prefact*(distance/r_0)**dens_exponent .gt. cutoff) then
                     grid(x,y,z)=cutoff
                  else 
                     grid(x,y,z)= nhyd*(distance/r_0)**(-dens_exponent)
                  end if
               end if
            end do
         end do
      end do
      mode_set=.true.
    end subroutine powerradcutoff

    subroutine powerplanecutoff
      implicit none
      integer(i4b):: x,y,z
      real   ( dp):: distance
      do z=-z_max, +z_max
         do y=-y_max, +y_max
            do x=-x_max, +x_max
               if (plane .eq. 'x' .or. plane .eq. 'X')&
                    distance = abs(x-x_center)
               if (plane .eq. 'y' .or. plane .eq. 'Y')&
                    distance = abs(y-y_center)
               if (plane .eq. 'z' .or. plane .eq. 'Z')&
                  distance = abs(z-z_center)
!
               if (distance .lt. 0.1) then
                  grid(x,y,z) = nhyd
               else
                  if (nhyd*(distance/z_0)**dens_exponent&
                       .gt. cutoff) then
                     grid(x,y,z)=cutoff
                  else 
                     grid(x,y,z)=nhyd*(distance/z_0)**(-dens_exponent)
                  end if
               end if
            end do
         end do
      end do
      mode_set=.true.
    end subroutine powerplanecutoff

    subroutine sphere
      implicit none
      integer(i4b):: x,y,z
      real   ( dp):: xdist, ydist, zdist, distance
      do z=-z_max,+z_max
         do y=-y_max,+y_max
            do x=-x_max,+x_max
               xdist=(x-x_center)
               ydist=(y-y_center)
               zdist=(z-z_center)
               distance = sqrt(xdist**2+ydist**2+zdist**2)
               if (distance .le. radius) then
                  grid(x,y,z)=nhyd
               end if
            end do
         end do
      end do
      mode_set=.true.
    end subroutine sphere

    subroutine ellipsoid
      implicit none
      integer(i4b)  :: x,y,z
      real(dp) :: xdist, ydist, zdist
      do z=-z_max,+z_max
         do y =-y_max, +y_max
            do x=-x_max, +x_max
               xdist=(x-x_center)
               ydist=(y-y_center)
               zdist=(z-z_center)
               if ((real(xdist**2)/real(x_half_axis**2))+&
                   (real(ydist**2)/real(y_half_axis**2))+&
                   (real(zdist**2)/real(z_half_axis**2)) .le. 1._dp) then
                  grid(x,y,z)=nhyd
               end if
            end do
         end do
      end do
      mode_set=.true.
    end subroutine ellipsoid

    subroutine random_remove
      implicit none
      integer (i4b):: x,y,z
      real    (dp) :: harvest    ! as opposed to (random_)seed ...
      call random_seed()
      do z=-z_max,+z_max
         do y=-y_max,+y_max
            do x=-x_max,+x_max
               call random_number(harvest)
               if (harvest .gt. macro_filling_factor)&
                    grid(x,y,z) = nhyd_fill
            end do
         end do
      end do
    end subroutine random_remove

    subroutine gaussian
      implicit none
      integer (i4b):: x,y,z
      real    ( dp):: xdist, ydist,zdist,distance
      print*, "GAUSSIAN"
      do z=-z_max,+z_max
         do y=-y_max,+y_max
            do x=-x_max, +x_max
               xdist=abs(x-x_center)
               ydist=abs(y-y_center)
               zdist=abs(z-z_center)
               distance=sqrt(xdist**2+ydist**2+zdist**2)
               grid(x,y,z)=nhyd*exp(-distance**2/sigma**2)
            end do
         end do
      end do
      mode_set=.true.
    end subroutine gaussian

end program
