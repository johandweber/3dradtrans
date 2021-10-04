!===============================================================================
!                                ****MODULE M_absorption****
! Calculates number of photons emitted by the radiation source, performs the
! radiative transfer and save the men intensity 
!===============================================================================

module M_absorption
!---routines:
!    ALLOCATE_PHOTON_ARRAYS
!    PHOTON_ZERO    
!    SET_RAY_DIRECTIONS
!    DIFFUSE_FIELD
!    GET_J
!    GET_J_PERIODIC
!    
 
  use M_data_types


  implicit none

  real(dp), dimension(1:8):: photons, absorbed
  real(dp), dimension(1:8):: absorbed_H, absorbed_HeI, absorbed_HeII
!
  real(dp), dimension(:,:,:), allocatable :: photons_freq_source,&
                                             photons_freq_source1 
!
  real(dp), dimension(:,:), allocatable   :: photons_freq, photons_freq1,&
                                             photons_old, photons_old1,&
                                             absorbed_freq, absorbed_freq1,&
                                             absorbed_freq_H, absorbed_freq1_H,&
                                             absorbed_freq_HeI, absorbed_freq1_HeI,&
                                             absorbed_freq_HeII, absorbed_freq1_HeII
  real(dp), dimension(:,:,:), allocatable :: photons_freqt, photons_freq1t,&
                                             photons_oldt, photons_old1t,&
                                             absorbed_freqt, absorbed_freq1t,&
                                             absorbed_freq_Ht, absorbed_freq1_Ht,&
                                             absorbed_freq_HeIt, absorbed_freq1_HeIt,&
                                             absorbed_freq_HeIIt, absorbed_freq1_HeIIt
!
contains
!
!================================================================================
  subroutine ALLOCATE_PHOTON_ARRAYS
!
!     Allocated the arrays required for the description of the
!     energy distributions of the sources
!
! called by main program
! 
   use M_definitions,            only: num_sources, points
!
   implicit none
   !
   integer:: errstat
!================================================================================
!   
   allocate(photons_freq_source(1:num_sources, 1:points, 1:8),stat=errstat)
   if (errstat .ne. 0) stop 'error allocating photons_freq_source'
!   
   allocate(photons_freq_source1(1:num_sources, 1:points, 1:8),stat=errstat)
   if (errstat .ne. 0) stop 'error allocating photons_freq_source1'
!
 end subroutine ALLOCATE_PHOTON_ARRAYS


 
!================================================================================ 
 subroutine PHOTON_ZERO(source_index)
!
!    Calculation of number of initial photons for source per frequency
!
!  called by main program
!   
   use M_natural_constants,      only: nc_light, nc_pi, nc_planck
   use M_definitions,            only: r_rs, points
   use M_data_input,             only: nug, nug_spect, weights_nu
!
 implicit none
!
   integer(i4b) :: f, c ,source_index
!=================================================================================
! 
 photons=0
!   
  do f=1,points-1
     do c =1,8
        photons_freq_source(source_index,f,c) = 4._dp*4._dp*nc_pi**2*&
             (R_RS(source_index)*6.9599E+10_dp)**2*nug_spect(f)*&
             nc_planck*nug(f)
        photons_freq_source1(source_index,f+1, c) = 4._dp*4._dp*nc_pi**2&
             *(R_RS(source_index)*6.9599E+10_dp)**2*nug_spect(f)*&
             nc_planck*nug(f+1)
        
      print*, "photon_zero",nc_light/nug_spect(f)*1E+8,4._dp*&
           (4._dp*nc_pi**2*(R_RS(source_index)*6.9599E+10_dp)**2*&
           nug_spect(f))*nc_planck*nug_spect(f)
      photons(c)=photons(c)+(photons_freq_source(source_index,f,c)&
           +photons_freq_source1(source_index,f+1,c))/2._dp*weights_nu(f) 
    end do
   end do
 end subroutine PHOTON_ZERO



!================================================================================
 subroutine SET_RAY_DIRECTIONS
!
!    Selects a random ray direction (required for MC description of diffuse
!    radiation field.
!
!  called by main program
!   
   use M_definitions,          only: x_max, y_max, z_max, diffuse_random_rays
   use M_grid_memory,          only: ray_directions,ray_octants
   use M_raysave,              only: rays
!
 implicit none
!
   integer(i4b):: x,y,z, ray_counter
   real(dp)    :: rnd, rnd2
!================================================================================
!
   call random_seed()
   do z=-z_max, z_max
      do y=-y_max, +y_max
         do x=-x_max,+x_max
            do ray_counter=1, diffuse_random_rays
               call random_number(rnd)
               call random_number(rnd2)
               ray_directions(x,y,z,ray_counter)=int(rnd*(rays-1))
               ray_octants(x,y,z, ray_counter)=int(rnd2*8)+1
            end do
         end do
      end do
   end do
 end subroutine SET_RAY_DIRECTIONS



!================================================================================
 subroutine DIFFUSE_FIELD
!
!    Performs the radiative transfer for the diffuse radiatiopn field
!
!  called by the main program
!
   use M_natural_constants,       only: nc_pi, nc_parsec
   use M_data_input,              only: nug
   use M_definitions,             only: x_max, y_max, z_max, diffuse_threshold,&
                                        l_cell, diffuse_random_rays,&
                                        x_min, y_min, z_min,points
                                        
   use M_grid_memory,             only: ne, nh_complete,eta, j_nu, chi,&
                                        ray_directions, ray_octants
   use M_raysave,                 only: ray_cells, coord_sign
!
 implicit none
!
   integer(i4b) :: x,y,z,f,c,e,r,rr_counter
   real(dp)     :: cell_zero
   real(dp)     :: tau, length_calc
   real(dp)     :: rnd, rnd2
   integer      :: cell_counter
!
   real(dp), dimension(64) :: photons_freq_diffuse, absorbed_freq_diffuse

!================================================================================
!
   cell_counter=0
!   
   !$OMP PARALLEL PRIVATE(cell_zero,tau,absorbed_freq_diffuse,& 
   !$OMP photons_freq_diffuse,&
   !$OMP z,y,x,c,e,r,rr_counter,rnd, rnd2,  length_calc)
   !$OMP DO
   do z=-z_max, +z_max
      do y=-y_max,+y_max
         do x=-x_max, +x_max
            if (ne(x,y,z) .lt. diffuse_threshold*nH_complete(x,y,z)) cycle 
            cell_counter=cell_counter+1
            do f=1, points-1
               cell_zero= (l_cell*nc_parsec)**3*eta(x,y,z,f)*4*nc_pi 
               photons_freq_diffuse=cell_zero/diffuse_random_rays
               if (x.eq. 0 .and. y .eq. 0 .and. z.eq. 0) then
                  print*, nug(f)/(1.0973731569E+05*2.99792458E+10_dp)
               end if
               do rr_counter=1,diffuse_random_rays
                  e=0
                  call random_seed()
                  call random_number(rnd)
                  call random_number(rnd2)
                  r=ray_directions(x,y,z,rr_counter)
                  c=ray_octants(x,y,z, rr_counter)
                  do while (&
                       abs(x+coord_sign(1,c)*(ray_cells(e,r)%x-x_min)) .le.&
                       x_max .and. &
                       abs(y+coord_sign(2,c)*(ray_cells(e,r)%y-y_min)) .le.&
                       y_max .and. &
                       abs(z+coord_sign(3,c)*(ray_cells(e,r)%z-z_min)) .le.&
                       z_max)
                     length_calc=ray_cells(e,r)%length
                     tau=chi(x+coord_sign(1,c)*(ray_cells(e,r)%x-x_min),&
                          y+coord_sign(2,c)*(ray_cells(e,r)%y-y_min),&
                          z+coord_sign(3,c)*(ray_cells(e,r)%z-z_min),f)*&
                          nc_parsec*length_calc
                     absorbed_freq_diffuse(rr_counter)=&
                          photons_freq_diffuse(rr_counter)*(1-exp(-tau))
                     photons_freq_diffuse(rr_counter)=&
                          photons_freq_diffuse(rr_counter)*exp(-tau)
                     if((chi(x+coord_sign(1,c)*(ray_cells(e,r)%x-x_min),&
                         y+coord_sign(2,c)*(ray_cells(e,r)%y-y_min),&
                         z+coord_sign(3,c)*(ray_cells(e,r)%z-z_min),f)).ne.0) then
                        !$OMP ATOMIC
                        j_nu(x+coord_sign(1,c)*(ray_cells(e,r)%x-x_min),&
                             y+coord_sign(2,c)*(ray_cells(e,r)%y-y_min),&
                             z+coord_sign(3,c)*(ray_cells(e,r)%z-z_min),f)=&
                             j_nu(x+coord_sign(1,c)*(ray_cells(e,r)%x-x_min),&
                             y+coord_sign(2,c)*(ray_cells(e,r)%y-y_min),&
                             z+coord_sign(3,c)*(ray_cells(e,r)%z-z_min),f)+&
                             absorbed_freq_diffuse(rr_counter)/&
                             (chi(x+coord_sign(1,c)*(ray_cells(e,r)%x-x_min),&
                             y+coord_sign(2,c)*(ray_cells(e,r)%y-y_min),&
                             z+coord_sign(3,c)*(ray_cells(e,r)%z-z_min),f) *&
                             4*nc_pi*(nc_parsec*l_cell)**3)
                     else
                        j_nu(x+coord_sign(1,c)*(ray_cells(e,r)%x-x_min),&
                             y+coord_sign(2,c)*(ray_cells(e,r)%y-y_min),&
                             z+coord_sign(3,c)*(ray_cells(e,r)%z-z_min),f)=0                      
                     end if
                     e=e+1
                  end do
               end do
            end do
         end do
      end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
   print*, "No. of cells considered for diffuse radiation field:", cell_counter
 end subroutine DIFFUSE_FIELD
 


!================================================================================
 subroutine GET_J
!
!    Performs the compuation of the radiative transfer,
!    not (explicitly) considering the diffuse radiation field.
!    If there is no case-B assumption, DIFFUISE_FIELD has to be called
!    separately
!
!  called by main program
!
   use M_natural_constants,   only: nc_pi, nc_light, nc_parsec
   use M_definitions,         only: x_max, y_max, z_max,x_min, y_min, z_min,&
                                    num_sources, source_x, source_y, source_z,&
                                    l_cell,&
                                    c_correct_on, o_write_escape_fraction,&
                                    t_all, t_step, start_activity, lifetime,&
                                    num_threads, points
   use M_grid_memory,         only: j_nu, xj_abs, chi
   use M_raysave,             only: coord_sign, ray_cells, escapefraction,&
                                    totaltau, rays, ray_split, l_max
!                                           
 implicit none
!
   integer(i4b) :: so,t,errstat,f,e,r,c
   integer(i4b) :: local_x, local_y, local_z
   real(dp)     :: length_sum, length_calc
   logical      :: exit_marker
!
   real(dp), dimension(0:points-1,1:8):: tau

!===============================================================================
!
   j_nu=0._dp
   xj_abs=0._dp
   if(allocated(escapefraction))&
        escapefraction=0._dp
!
   do so=1, num_sources
      if ((t_all+t_step)/31557600._dp .lt. start_activity(so) .or.&
           (t_all+t_step)/31557600._dp .gt. start_activity(so)+lifetime(so)) cycle
      if (allocated(totaltau))&
           totaltau=0._dp
!
      allocate(photons_freqt(0:points,1:8,0:num_threads-1), stat=errstat)
      if (errstat .ne. 0) stop 'error allocating photons_freq'
      allocate(absorbed_freqt(0:points,1:8,0:num_threads-1))
      if (errstat .ne. 0) stop 'error allocating absorbned_freq'
!
      !$OMP PARALLEL DEFAULT(SHARED)&
      !$OMP PRIVATE (errstat,r,c,f,e,length_sum,tau,&
      !$OMP          local_x, local_y, local_z, length_calc, exit_marker)
      !$OMP DO
      do t=0, num_threads-1
         do r=t,rays-1, num_threads
             do c=1,8
               do f=1,points-1
                  e=0
                  length_sum=0
                  exit_marker = .false.
                  photons_freqt(f,c,t)=photons_freq_source(so,f,c)/&
                       (12*4**L_MAX*RAY_SPLIT(R))
                  do while (.true.)   ! The exit condition follows in the next if clause
                     local_x=source_x(so)+coord_sign(1,c)*&
                          (ray_cells(e,r)%x-x_min)
                     local_y=source_y(so)+coord_sign(2,c)*&
                          (ray_cells(e,r)%y-y_min)
                     local_z=source_z(so)+coord_sign(3,c)*&
                          (ray_cells(e,r)%z-z_min)
!
                     if(abs(local_x) .gt. x_max .or. &
                          abs(local_y) .gt. y_max .or. &
                          abs(local_z).gt. z_max)&
                          exit
 !                    
                     exit_marker=.false.
                     if (c_correct_on) then
                        length_sum=length_sum+ray_cells(e,r)%length*nc_parsec
                        if (length_sum .gt. nc_light*(t_all-start_activity(so)))&
                             then
                           exit
                        else
                           length_calc=ray_cells(e,r)%length
                        end if
                     else
                        continue
                        length_sum=length_sum+ray_cells(e,r)%length*nc_parsec
                        length_calc=ray_cells(e,r)%length
                     end if
                     tau(f,c)=chi(local_x,local_y,local_z,f)*nc_parsec*length_calc
                     if (o_write_escape_fraction)&
                          totaltau(r,c,f)=totaltau(r,c,f)+tau(f,c)
!                     
                     absorbed_freqt(f,c,t)=photons_freqt(f,c,t)*(1-exp(-tau(f,c)))
                     photons_freqt(f,c,t)= photons_freqt(f,c,t)*exp(-tau(f,c))
                     if ((chi(local_x,local_y, local_z,f)) .ne. 0) then
                        !                        !$OMP ATOMIC  
                        j_nu(local_x, local_y, local_z,f)=&
                             j_nu(local_x,local_y,local_z,f)+&
                             absorbed_freqt(f,c,t)/&
                             (chi(local_x, local_y, local_z,f) *&
                             4*nc_pi*(nc_parsec*l_cell)**3)
                     else
                        j_nu(local_x,local_y, local_z,f)=0
                     end if
                     e=e+1
                  end do             ! end of raysegments
               end do                ! end of frequency points
            end do                   ! end of octants
         end do                      ! end of rays         
      end do                         ! end of threads
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate (photons_freqt)
      deallocate (absorbed_freqt)
!     
      if (o_write_escape_fraction) then   
         !$OMP PARALLEL PRIVATE(c,r)
         !$OMP DO
         do f=1,points-1
            do c=1,8
               do r=0, rays-1
                  escapefraction(so,f)=escapefraction(so,f)+&
                       exp(-totaltau(r,c,f))/(8*rays*ray_split(r))
               end do
            end do
         end do
         !$OMP END DO
         !$OMP END PARALLEL
      end if
!
   end do                            ! end of sources
 end subroutine GET_J
 
!================================================================================ 
 integer(i4b) function FLIP_X (xnotflipped)
!
!    Helper function for periodic boundary_conditions
!    x-coordinates for periodic boundary conditions 
!
!  called by GET_J_PERIODIC
!
   use M_definitions,       only: x_max
!
   implicit none
!
   integer::  xnotflipped
!================================================================================
!
   if (xnotflipped .ge. -x_max) then
      flip_x= -x_max+mod(xnotflipped+x_max, 2*x_max+1)
   else
      flip_x =  mod(xnotflipped+x_max+1, 2*x_max+1)+x_max
   end if
 end function FLIP_X


 
!================================================================================ 
 integer(i4b) function FLIP_Y (ynotflipped)
!
!    Helper function for periodic boundary_conditions
!    y -coordinates for periodic boundary conditions 
!
!  called by GET_J_PERIODIC
!
   use M_definitions,       only: y_max
!
 implicit none
!
   integer::  ynotflipped
!================================================================================

   if (ynotflipped .ge. -y_max) then
      flip_y= -y_max+mod(ynotflipped+y_max, 2*y_max+1)
   else
      flip_y =  mod(ynotflipped+y_max+1, 2*y_max+1)+y_max
   end if
 end function FLIP_Y
 


!================================================================================ 
 integer function FLIP_Z (znotflipped)
!
!    Helper function for periodic boundary_conditions
!    z-coordinates for periodic boundary conditions 
!
!  called by GET_J_PERIODIC
!
   use M_definitions,       only: z_max
!
 implicit none
!
   integer::  znotflipped
!================================================================================
!
   if (znotflipped .ge. -z_max) then
      flip_z= -z_max+mod(znotflipped+z_max, 2*z_max+1)
   else
      flip_z =  mod(znotflipped+z_max+1, 2*z_max+1)+z_max
   end if
 end function FLIP_Z
 
 
!================================================================================ 
 subroutine GET_J_PERIODIC
!
! Radiative transfer for periodic boundary conditions
!
! called by main program
!
   use M_definitions,           only: x_max, y_max, z_max, x_min, y_min, z_min,&
                                      source_x, source_y, source_z,&
                                      c_correct_on, l_cell,&
                                      num_sources, t_all, t_step,&
                                      start_activity, lifetime, num_threads,&
                                      points, o_write_escape_fraction
   
   use M_natural_constants,     only: nc_pi, nc_light, nc_parsec
   use M_grid_memory,           only: j_nu, xj_abs, chi
   use M_raysave,               only: escapefraction,totaltau, rays,&
                                      ray_split, ray_cells, coord_sign,&
                                      l_max
!
   implicit none
!
   integer(i4b) :: so,t,errstat,f,e,r,c
   real(dp)     :: length_sum, length_calc
   logical      :: exit_marker
!
   real(dp), dimension(0:points-1,1:8):: tau
!================================================================================
!
   j_nu=0._dp
   xj_abs=0._dp
   if(allocated(escapefraction))&
        escapefraction=0._dp
   do so=1, num_sources
      if ((t_all+t_step)/31557600._dp .lt. start_activity(so) .or.&
           (t_all+t_step)/31557600._dp .gt. start_activity(so)+lifetime(so)) cycle
      if (allocated(totaltau))&
           totaltau=0._dp
      
      !$OMP PARALLEL DEFAULT(SHARED)  PRIVATE&
      !$OMP (errstat,r,c,f,e,length_sum,tau, absorbed_freq, photons_freq,&
      !$OMP                                    length_calc, exit_marker)
      !$OMP DO
      do t=0, num_threads-1
         allocate(photons_freq(0:points,1:8), stat=errstat)
         if (errstat .ne. 0) stop 'error allocating photons_freq'
         allocate(absorbed_freq(0:points,1:8))
         if (errstat .ne. 0) stop 'error allocating absorbned_freq'
         do r=t,rays-1, num_threads
            do c=1,8
               do f=1,points-1
                  e=0
                  length_sum=0
                  exit_marker = .false.
                  photons_freq(f,c)=photons_freq_source(so,f,c)/&
                       (12*4**L_MAX*RAY_SPLIT(R))                 
                  do while (&
                       ray_cells(e,r)%x .le. (2*x_max) .and. &
                       ray_cells(e,r)%y .le. (2*y_max) .and. &
                       ray_cells(e,r)%z .le. (2*z_max) .and.&
                       e .lt. 3*max(2*x_max+1,2*y_max+1,2*z_max+1))                     
!
                     exit_marker=.false.
                     if (c_correct_on) then
                        length_sum=length_sum+ray_cells(e,r)%length*nc_parsec
                        if (length_sum .gt. nc_light*&
                             (t_all-start_activity(so))) then
                           exit
                        else
                           length_calc=ray_cells(e,r)%length
                        end if
                     else
                        continue
                        length_sum=length_sum+ray_cells(e,r)%length*nc_parsec
                        length_calc=ray_cells(e,r)%length
                     end if
                     tau(f,c)=chi( flip_x(source_x(so)+coord_sign(1,c)*&
                          (ray_cells(e,r)%x-x_min)),&
                          flip_y(source_y(so)+coord_sign(2,c)*&
                          (ray_cells(e,r)%y-y_min)),&
                          flip_z(source_z(so)+coord_sign(3,c)*&
                          (ray_cells(e,r)%z-z_min)),f)*nc_parsec*length_calc
                     if (o_write_escape_fraction)&
                          totaltau(r,c,f)=totaltau(r,c,f)+tau(f,c)
                     
                     absorbed_freq(f,c)=photons_freq(f,c)*(1-exp(-tau(f,c)))
                     photons_freq(f,c)= photons_freq(f,c)*exp(-tau(f,c))
                     if ((chi(flip_x( source_x(so)+coord_sign(1,c)*&
                          (ray_cells(e,r)%x-x_min)),&
                          flip_y(source_y(so)+coord_sign(2,c)*&
                          (ray_cells(e,r)%y-y_min)),&
                          flip_z(source_z(so)+coord_sign(3,c)*&
                          (ray_cells(e,r)%z-z_min)),f)) .ne. 0) then                                !               
                        j_nu(flip_x(source_x(so)+coord_sign(1,c)*&
                             (ray_cells(e,r)%x-x_min)),&
                             flip_y(source_y(so)+coord_sign(2,c)*&
                             (ray_cells(e,r)%y-y_min)),&
                             flip_z(source_z(so)+coord_sign(3,c)*&
                             (ray_cells(e,r)%z-z_min)),f)=&
                             j_nu(flip_x(source_x(so)+coord_sign(1,c)*&
                             (ray_cells(e,r)%x-x_min)),&
                             flip_y(source_y(so)+coord_sign(2,c)*&
                             (ray_cells(e,r)%y-y_min)),&
                             flip_z(source_z(so)+coord_sign(3,c)*&
                             (ray_cells(e,r)%z-z_min)),f)+&
                             absorbed_freq(f,c)/&
                             (chi( flip_x(source_x(so)+coord_sign(1,c)*&
                             (ray_cells(e,r)%x-x_min)),&
                             flip_y(source_y(so)+coord_sign(2,c)*&
                             (ray_cells(e,r)%y-y_min)),&
                             flip_z(source_z(so)+coord_sign(3,c)*&
                             (ray_cells(e,r)%z-z_min)),f) *&
                             4*nc_pi*(nc_parsec*l_cell)**3)
                     else
                        j_nu(flip_x(source_x(so)+coord_sign(1,c)*&
                             (ray_cells(e,r)%x-x_min)),&
                             flip_y(source_y(so)+coord_sign(2,c)*&
                             (ray_cells(e,r)%y-y_min)),&
                             flip_z(source_z(so)+coord_sign(3,c)*&
                             (ray_cells(e,r)%z-z_min)),f)=0
                     end if
                     e=e+1
                  end do             ! end of raysegments
               end do                ! end of frequency points
            end do                   ! end of octants
         end do                      ! end of rays
         deallocate (photons_freq)
         deallocate (absorbed_freq)
      end do                         ! end of threads
      !$OMP END DO
      !$OMP BARRIER
      !$OMP END PARALLEL
!
      if (o_write_escape_fraction) then   
         !$OMP PARALLEL PRIVATE(c,r)
         !$OMP DO
         do f=1,points-1
            do c=1,8
               do r=0, rays-1
                  escapefraction(so,f)=escapefraction(so,f)+&
                       exp(-totaltau(r,c,f))/(8*rays*ray_split(r))
               end do
            end do
         end do
         !$OMP END DO
         !$OMP END PARALLEL
      end if
   end do                            !end of sources
 end subroutine GET_J_PERIODIC


 
!==============================================================================
 subroutine GET_J_SPHERICAL
!    Computes J for the spherically-symmetric case
!    Has yet to be implemented... 
   use M_definitions,       only: points
   implicit none
   integer(i4b)                :: f ! local frequency counter
!==============================================================================
!
   do f=1, points-1
      
   end do   
!
 end subroutine GET_J_SPHERICAL

!==============================================================================
 real(dp) function spherical_shell(rinner,router)
!    Compoutes the volume of a spherical shell
!    Has yet to be implemented
   use M_natural_constants, only: nc_pi, nc_parsec
!
   implicit none
   real (dp) :: rinner, router               ! in parsec
!==============================================================================
!  return value in cm^-3
   spherical_shell = (4._dp/3._dp)*nc_pi*nc_parsec**3*(router**3-rinner**3)
 end function spherical_shell

end module M_absorption

!==============================================================================
!                            ****END MODULE M_absorption****
!==============================================================================
