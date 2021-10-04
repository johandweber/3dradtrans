!===============================================================================
!                               ****MODULE M_RAYSAVE****
! In this module the variables needed for the geometrical part of the raytracing 
! are stored and allocated. It is needed by the  module traceray 
!===============================================================================

module M_raysave
!--- routines
!    ALLOCATE_RAYS 
!    INIT_GEOMETRY
!    SIGNUM

  use M_data_types
  
  implicit none
  
  type cell                            ! geometrical properties of a cell
    integer(i2b)::x
    integer(i2b)::y
    integer(i2b)::z
    real(dp):: length
  end type cell
!  
  ! all cells
  type(cell), dimension(:,:), allocatable  :: ray_cells      
!
! "weaken" the ray to compensate mirroring
  integer(i4b), dimension(:), allocatable  :: ray_split       
!  
! number of 3D-Pixels used by one ray (not used so far)  
  integer(i2b),dimension(:), allocatable   :: total_length       
!
  real(dp), dimension (:,:,:), allocatable :: totaltau
!  
  real(dp), dimension (:,:), allocatable   :: escapefraction
!
! cells per ray (possibly not necessary)
  integer, dimension(:), allocatable       :: ray_number_of_cells   
!
  real(dp), dimension(1:3):: directional_vector
!  
  integer:: number_of_subvolumes
  integer:: subvolume_counter
  integer:: rays  = 0                                           
  integer:: ray_counter
  integer:: ray_max                 ! number of "valid" rays
  integer:: l_max
  
  integer, dimension(1:3,1:8):: coord_sign
  
  contains
 
!================================================================================
    SUBROUTINE ALLOCATE_RAYS
!
!    Allocates the arrays describing the rays
!
!   called by main program
!
!
      use M_natural_constants ,only: nc_pi
      use M_definitions       ,only: f_a, num_sources, o_write_escape_fraction,&
                                     points, x_max, y_max, z_max
!   
    implicit none
!
      integer(i4b):: errstat
      integer(i4b):: dummy_counter1, dummy_counter2
!================================================================================ 
!
      coord_sign(1,1)=+1 ; coord_sign(2,1)=+1; coord_sign(3,1)=+1
      coord_sign(1,2)=+1 ; coord_sign(2,2)=+1; coord_sign(3,2)=-1
      coord_sign(1,3)=+1 ; coord_sign(2,3)=-1; coord_sign(3,3)=+1
      coord_sign(1,4)=+1 ; coord_sign(2,4)=-1; coord_sign(3,4)=-1
      coord_sign(1,5)=-1 ; coord_sign(2,5)=+1; coord_sign(3,5)=+1
      coord_sign(1,6)=-1 ; coord_sign(2,6)=+1; coord_sign(3,6)=-1
      coord_sign(1,7)=-1 ; coord_sign(2,7)=-1; coord_sign(3,7)=+1
      coord_sign(1,8)=-1 ; coord_sign(2,8)=-1; coord_sign(3,8)=-1
      
 
      l_max=ceiling(log(f_a*4._dp*nc_pi*&
           & ((2._dp*1._dp*x_max+1._dp)**2+(2*1._dp*y_max*2+1._dp)**2+&
           (2*1._dp*z_max+1)**2))/log(4._dp))
!      
      if (o_write_escape_fraction) then
         allocate(totaltau(0:12*4**l_max/6,1:8,0:points))
         allocate(escapefraction(1:num_sources,0:points))
      end if
!            
      write(*,*) '# l_max:', l_max
      write(*,*) '# number of rays: ',12*4**l_max
!      
      errstat=0
      allocate(ray_cells(0:3*max(2*x_max+3, 2*y_max+3, 2*z_max+3),&
           0:12*4**l_max/6), stat=errstat)
      allocate(ray_split(0:12*4**l_max/6), stat=errstat)
      allocate(ray_number_of_cells(0:12*4**l_max/6), stat=errstat)
!      
      if (errstat .ne. 0) then
         stop 'memory allocation error'
      end if
!      
!     initialize ray_cells with (absurd) initial conditions
      do dummy_counter1=0,12*4**l_max/6 
         do dummy_counter2=0,3*max(2*x_max+1,2*y_max+1,2*z_max+1)
            !print*, dummy_counter1, dummy_counter2
            ray_cells(dummy_counter2, dummy_counter1)%x=20000
            ray_cells(dummy_counter2, dummy_counter1)%y=20000
            ray_cells(dummy_counter2, dummy_counter1)%z=20000
         end do
      end do
    end subroutine ALLOCATE_RAYS
!
end module M_raysave
!===============================================================================
!                             ****END MODULE M_RAYSAVE****
!===============================================================================
!
!===============================================================================
!                             ****MODULE M_grid_memory****
! In this module the arrays for the hydrogen density and the density of 
! absorbed photons are defined and allocated.
!================================================================================

module M_grid_memory
!--- routines:
!    ALLOCATE_GRID
!
  use M_data_types
  implicit none
  real( dp), dimension(:,:,:), allocatable     :: nHI, nHII, absorbed_sum,&
                                                  nH_complete, nHI_old,&
                                                  absorbed_sum_H,&
                                                  absorbed_sum_HeI,&
                                                  absorbed_sum_HeII,&
                                                  temperature,&
                                                  temperature_old,&
                                                  energycontent, thermal_pressure
  real( sp), dimension(:,:,:), allocatable     :: nHIf, nHIIf, absorbed_sumf,&
                                                  nH_completef, nHI_oldf,&
                                                  absorbed_sum_Hf,&
                                                  absorbed_sum_HeIf,&
                                                  absorbed_sum_HeIIf,&
                                                  temperaturef,&
                                                  temperature_oldf,&
                                                  energycontentf, thermal_pressuref

  real( dp), dimension(:,:,:), allocatable     :: nHeI, nHeI_old, nHeII,nHeIII,&
                                                  nHeII_old, nHe_complete, ne
  real( sp), dimension(:,:,:), allocatable     :: nHeIf, nHeI_oldf, nHeIIf,nHeIIIf,&
                                                  nHeII_oldf, nHe_completef, nef

!
  real( dp), dimension(:,:,:), allocatable     :: cool_H, cool_He, cool_C,&
                                                  cool_N,cool_O,&
                                                  cool_Ne, cool_S,&
                                                  cool_Ar, cool_ff
  real( sp), dimension(:,:,:), allocatable     :: cool_Hf, cool_Hef, cool_Cf,&
                                                  cool_Nf,cool_Of,&
                                                  cool_Nef, cool_Sf,&
                                                  cool_Arf, cool_fff

!
  real( dp), dimension(:,:,:), allocatable     :: heat_H, heat_He
  real( sp), dimension(:,:,:), allocatable     :: heat_Hf, heat_Hef

!
  real(dp), dimension(:,:,:), allocatable     ::  add_heating_cooling
!
  integer(i4b), dimension (:,:,:), allocatable :: mask
!
  real( dp), dimension (:,:,:,:), allocatable  :: nC, nO, nN , nNe, nS, nAr
  real( sp), dimension (:,:,:,:), allocatable  :: nCf, nOf, nNf , nNef, nSf, nArf

!
  real(dp), dimension (:,:,:,:), allocatable   :: element_abundances
!
  real( dp), dimension(:,:,:,:), allocatable   :: absorbed_sum_copy,&
                                                  absorbed_sum_copy_H,&
                                                  absorbed_sum_copy_HeI,&
                                                  absorbed_sum_copy_HeII
!
  real( dp)   , dimension (:,:,:,:),   allocatable :: chi, j_nu, xj_abs, eta
  real( dp)   , dimension (:,:,:,:,:), allocatable :: j_nu_copy
  integer(i4b), dimension (:,:,:,:)  , allocatable :: ray_directions, ray_octants
  real( dp)   , dimension (:,:,:,:)  , allocatable :: add_heating_values,&
                                                      add_heating_poly

!

  contains
!
!================================================================================  
    subroutine ALLOCATE_GRID
!
!    Allocates the arrays describing the proberties of the grid cells
!
!  called by main program
!
      use M_definitions   , only: spherical, x_max, y_max, z_max, rp_max,&
                                  num_threads, points, log_heating_cooling,&
                                  include_metals,&
                                  diffuse, diffuse_random_rays,&
                                  o_thermal_pressure,&
                                  o_abundan, optabu_3d, use_mask,&
                                  add_heating, extra_heating,&
                                  binary_output,&
                                  Heions, Cions, Nions, Oions, Neions,&
                                  Sions, Arions
!
    implicit none
!
      integer:: errstat=0
!===============================================================================      
!
      if (spherical) then        ! remapping the spherical coordinates to the
                                 ! Cartesian ones 
         x_max=0
         y_max=0
         z_max=rp_max
      end if
!
      allocate (ne(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
!      
      allocate (nHI(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
      allocate (nHII(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
      allocate (nH_complete(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                stat=errstat)
      allocate (nHeI(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
      allocate (nHeII(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
      allocate (nHeIII(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
      allocate (nHe_complete(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (temperature(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (temperature_old(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (ray_directions(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,&
           1:diffuse_random_rays), stat=errstat)
      allocate (ray_octants(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,&
           1:diffuse_random_rays), stat=errstat)
      allocate (energycontent(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
!
      if (errstat .ne. 0) stop "error allocating basic files"
!               
      if (include_metals) then
         if (Cions .ne. 3 .and. Cions .ne. 4)&
              stop "There have to be 3 or 4 ionzation stages of C"
         allocate (nC(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Cions),&
              stat=errstat)
         if (Nions .ne. 3 .and. Nions .ne. 4)&
              stop "There have to be 3 or 4 ionzation stages of N"
         allocate (nN(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Nions),&
              stat=errstat)
         if (Oions .ne. 3 .and. Oions .ne. 4)&
              stop "There have to be 3 or 4 ionzation stages of O"         
         allocate (nO(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Oions),&
              stat=errstat)
          if (Neions .ne. 3 .and. Neions .ne. 4)&
              stop "There have to be 3 or 4 ionzation stages of Ne"                 
         allocate (nNe(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max, 0:Neions),&
              stat=errstat)
          if (Sions .ne. 3 .and. Sions .ne. 4)&
              stop "There have to be 3 or 4 ionzation stages of S"                          
         allocate (nS(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Sions),&
              stat=errstat)
          if (Arions .ne. 3 .and. Arions .ne. 4)&
              stop "There have to be 3 or 4 ionzation stages of Ar"                                   
         allocate (nAr(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Arions),&
              stat=errstat)
      end if
!
      if (log_heating_cooling) then
         allocate(cool_H(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(cool_He(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(cool_ff(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(heat_H(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(heat_He(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         if (include_metals) then
            allocate(cool_C(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_N(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_O(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_Ne(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_S(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_Ar(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
         end if
      end if
!
      if (add_heating) then
         allocate(add_heating_poly(1:4,-x_max:+x_max, -y_max:+y_max,&
              -z_max:+z_max))
      end if
!
      if (extra_heating) then
         allocate(add_heating_cooling(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max))
      end if
!
      if (use_mask)&
           allocate(mask(-x_max:+x_max, -y_max:+y_max, -z_max:+z_max))
!
      allocate (absorbed_sum(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (absorbed_sum_H(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (absorbed_sum_HeI(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (absorbed_sum_HeII(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (nHI_old(-x_max-1:+x_max,-y_max:+y_max,-z_max:+z_max),&
           stat=errstat)
      allocate (chi(-x_max:x_max,-y_max:y_max,-z_max:z_max,0:points-1),&
           stat=errstat)
! 
      
! 
     if (diffuse) &
           allocate (eta(-x_max:x_max,-y_max:y_max,-z_max:z_max,0:points-1),&
           stat=errstat)
!
      allocate (j_nu(-x_max:x_max,-y_max:y_max,-z_max:z_max,0:points-1),&
           stat=errstat)
      allocate (j_nu_copy(-x_max:x_max,-y_max:y_max,-z_max:z_max,&
           0:points-1,0:num_threads-1), stat=errstat)
      allocate (xj_abs(-x_max:x_max,-y_max:y_max,-z_max:z_max, 0:points-1),&
           stat=errstat)
!           
!
      if(o_thermal_pressure) &
           allocate(thermal_pressuref(-x_max:x_max,-y_max:y_max,-z_max:z_max),&
           stat=errstat)
!      
      if (o_abundan) then
         if (optabu_3d) then 
            allocate(element_abundances(1:30,&
                 -x_max:+x_max, -y_max:+y_max, -z_max+z_max))
         else
            allocate(element_abundances(1:30,0:0, 0:0, 0:0)) !only one cell
         end if
      end if
!      
      if (errstat .ne. 0) then
         stop 'error allocating memory'
      end if
!
      if (binary_output) then
         call allocate_float_values
      end if
!
    contains
!
      subroutine allocate_float_values
        implicit none
        allocate (nef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
!      
        allocate (nHIf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
        allocate (nHIIf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
        allocate (nH_completef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
             stat=errstat)
        allocate (nHeIf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
        allocate (nHeIIf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
        allocate (nHeIIIf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max), stat=errstat)
        allocate (nHe_completef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
             stat=errstat)
        allocate (temperaturef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
             stat=errstat)
        allocate (temperature_oldf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
             stat=errstat)
        allocate (energycontentf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
             stat=errstat)
        !
        if (errstat .ne. 0) stop "error allocating basic files"
!               
        if (include_metals) then
           allocate (nCf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Cions),&
                stat=errstat)
           allocate (nNf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Nions),&
                stat=errstat)
           allocate (nOf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Oions),&
                stat=errstat)
           allocate (nNef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max, 0:Neions),&
                stat=errstat)
           allocate (nSf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Sions),&
                stat=errstat)
           allocate (nArf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max,  0:Arions),&
                stat=errstat)
      end if
!
      if (log_heating_cooling) then
         allocate(cool_Hf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(cool_Hef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(cool_fff(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(heat_Hf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         allocate(heat_Hef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
              stat=errstat)
         if (include_metals) then
            allocate(cool_Cf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_Nf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_Of(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_Nef(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_Sf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
            allocate(cool_Arf(-x_max:+x_max,-y_max:+y_max,-z_max:+z_max),&
                 stat=errstat)
         end if
      end if
!
      if (add_heating) then
         allocate(add_heating_poly(1:4,-x_max:+x_max, -y_max:+y_max,&
                  -z_max:+z_max))
      end if
!
!
      if(o_thermal_pressure) &
           allocate(thermal_pressuref(-x_max:x_max,-y_max:y_max,-z_max:z_max),&
           stat=errstat)
!      
!      
      if (errstat .ne. 0) then
         stop 'error allocating memory in allocate_float_values'
      end if
!        
    end subroutine allocate_float_values
      
    end subroutine ALLOCATE_GRID
end module M_grid_memory
!===============================================================================
!                          ****END MODULE M_grid_memory****
!===============================================================================
!
!===============================================================================
!                          ****MODULE M_geometry****
!==============================================================================
! Here the calculations for the geometrical part of the raytracing are performed. 
! The result is saved in  an array that can be read by other routines or the main 
! program
!--------------------------------------------------------------------------------
!
module M_geometry
!--- routines TRACERAY

  use M_data_types  
! 
  implicit none
!  
  contains

!===============================================================================  
    subroutine TRACERAY
!
!    Computes that geometrical path of a ray through the grid.
!    Note: It does not perform the radiative transfer
!    This is done in the routines GET_J or GET_J_PERIODIC (and possible DIFFUSE)
!
!  called by INIT_GEOMETRY
!
      use M_definitions,     only: x_max, y_max, z_max, l_cell
      use M_raysave,         only: directional_vector, ray_split, rays,&
                                   ray_cells, rays
!
    implicit none
!
! savety parameter to prevent divison by zero
      real(dp), parameter :: eps=1E-8               
!      
! see AMN 99, Formula 13
      real(dp)     :: r_x, r_y, r_z                  
!
! "old" ray lenth to be subtracted from the
      real(dp)     :: grid_old             
!
! indexes for the trace of the ray
      integer(i4b) :: il, jl, kl  
!
! sign of the direction in x-, y- or z- dirction                   
      integer(i4b) :: s_mu, s_gamma, s_zeta     
!
! loop index for subvolumes     
      integer(i4b) :: element_counter                
!================================================================================
!      
      grid_old=0._dp                                 ! initialize grid_old
      il=-x_max; jl=-y_max; kl=-z_max                ! initialize il, jl and kl
      element_counter=0                              ! set element_counter to 0 
!      
      s_mu=signum(directional_vector(1)+eps)
      s_gamma=signum(directional_vector(2)+eps)
      s_zeta=signum(directional_vector(3)+eps)
!      
      ray_split(rays)=1
!            
      if ( (s_mu .lt. 0)  .or.  (s_gamma .lt. 0)  .or. ( s_zeta .lt. 0)) return 
!      
      do while ((abs(il).le.x_max+1) .and. (abs(jl).le.y_max+1) .and.&
           ( abs(kl) .le. z_max+1) )
         
         r_x=abs((il+(s_mu+1)/2._dp)-(-x_max+0.5_dp))/&
              max(eps, abs(directional_vector(1)))
         r_y=abs((jl+(s_gamma+1)/2._dp)-(-y_max+0.5_dp))/&
              max (eps, abs(directional_vector(2)))
         r_z=abs((kl+(s_zeta+1)/2._dp)-(-z_max+0.5_dp))/&
              max(eps, abs(directional_vector(3)))
!         
         if (r_x .le. min (r_y,r_z)) then
            if (abs(il) .le. x_max+1) then
               ray_cells(element_counter, rays)%length=(r_x-grid_old)*l_cell
               ray_cells(element_counter, rays)%x=int(il, kind=2)
               ray_cells(element_counter, rays)%y=int(jl, kind=2)
               ray_cells(element_counter, rays)%z=int(kl, kind=2)
               grid_old=r_x
            else 
               exit
            end if
            il=il+s_mu
            element_counter=element_counter+1        
            cycle      
         end if
!         
         if (r_y .le. min (r_x,r_z)) then
            if (abs(jl) .le. y_max+1) then
               ray_cells(element_counter, rays)%length=(r_y-grid_old)*l_cell
               ray_cells(element_counter, rays)%x=int(il, kind=2)
               ray_cellS(element_counter, rays)%y=int(jl, kind=2)
               ray_cells(element_counter, rays)%z=int(kl, kind=2)
               grid_old=r_y
            else 
               exit
            end if
            jl=jl+s_gamma
            element_counter=element_counter+1
            cycle      
         end if
!         
         if (r_z .le. min (r_x,r_y)) then
            if (abs(kl) .le. z_max+1) then
               ray_cells(element_counter, rays)%length=(r_z-grid_old)*l_cell
               ray_cells(element_counter, rays)%x=int(il, kind=2)
               ray_cells(element_counter, rays)%y=int(jl, kind=2)
               ray_cells(element_counter, rays)%z=int(kl, kind=2)
               grid_old=r_z
            else 
               exit
            end if
            kl=kl+s_zeta
            element_counter=element_counter+1
            cycle      
         end if
      end do
      
      if (directional_vector(1) .lt. eps) then
         ray_split(rays)=ray_split(rays)*2
      end if
!      
      if (directional_vector(2) .lt. eps) then
         ray_split(rays)=ray_split(rays)*2
      end if
!      
      if (directional_vector(3) .lt. eps) then
         ray_split(rays)=ray_split(rays)*2
      end if
      rays=rays+1    
    end subroutine TRACERAY



!================================================================================    
  subroutine INIT_GEOMETRY
!
!    Precomputes the rays (with resect to the origin , for different sources there 
!    is a different offset
!
!  called by main program
!
    use M_raysave,           only: l_max, directional_vector
    use M_min_healpix,       only: PIX2VEC_NEST
!
  implicit none
!
    integer(i4b):: ray_counter
!================================================================================
!
    do ray_counter=0, 12*4**l_max-1
      call pix2vec_nest(2**l_max,ray_counter,directional_vector)
      call traceray
    end do
  end subroutine INIT_GEOMETRY 



!===============================================================================
  integer(i4b) function SIGNUM(a)
!
!    Provides a version of the sign-function
!
!   called by TRACERAY
!  
  implicit none
!
    real(dp)    :: a
    integer(i4b):: fn
!================================================================================
!
    fn=1
    if(a < 0) then
      fn=-1
    end if
    signum=fn
  end function SIGNUM
  
end module M_geometry 
!=============================================================================== 
!                      ****END MODULE M_geometry****
!===============================================================================
!
!
!===============================================================================
!                        M_geometrical_utilities
!===============================================================================
!  Helper functions concerning geometrical aspects of the grid
!  M 17.08 jweber
!
module M_geo_utils
!
  use M_data_types
!
  implicit none

 contains

!============================================================================== 
   subroutine consec2coord(consec, x,y,z)
!  translates a consecutive number (with respect to the beginning of the
!  allocated array) into grix coordinates x,y, and z
!
     use M_definitions,                   only: x_max,y_max, z_max
!
     implicit none 
!
     integer(i4b)   :: consec, x,y,z
     x=mod(consec-1, (2*x_max+1))-x_max
     y=(mod(consec-1, (2*x_max+1)*(2*y_max+1))/(2*x_max+1))-y_max
     z=(consec-1)/((2*y_max+1)*(2*x_max+1))-z_max
   end subroutine consec2coord

 end module M_geo_utils
!===============================================================================
!                  +***END MODULE M_geometrical_utilities****
!===============================================================================
!
!
!===============================================================================
!                       *****Module M_cosmology****
!===============================================================================
!"converts" redshift into scale facor and the other way round
!===============================================================================
module M_cosmology

  use M_data_types

  implicit none

contains

!================================================================================  
  real(dp) function ZTOAGE(z, omega_m, H_0)
! 
!    Converts redshift into worldage (since Big Bang).
!
!   called by CALC_HYDROGEN
!
    use M_natural_constants,      only: nc_km, nc_megaparsec
!
  implicit none
    real(dp):: z, omega_m, H_0, H_0_copy
    real(dp):: omega_lambda, scalefactor
    real(dp):: age
!================================================================================
!
    omega_lambda=1.0_sp-omega_m
    scalefactor=1.0_dp/(z+1.0e+00)
    H_0_copy=H_0*(nc_km/nc_megaparsec)
    age=2._dp/3._dp*(1._dp/(H_0_copy*sqrt(omega_lambda)))*&
         log(&
         sqrt(omega_lambda/omega_m*scalefactor**3) +&
         sqrt(1.+omega_lambda/omega_m*scalefactor**3)&
         )
    ztoage=age
  end function ZTOAGE
 


!=============================================================================== 
  double precision function AGETOZ (age, omega_m, H_0)
!
!   Converts worldage (since Big Bang) into redshift
!
!   called by CALC_HYDROGEn
!
    use M_natural_constants,        only: nc_km, nc_megaparsec, nc_gyr
!
  implicit none
    double precision:: age, age_copy, omega_m, H_0, H_0_copy
    double precision:: omega_lambda, scalefactor
    double precision:: z
!================================================================================
!
    omega_lambda=1.0d0-omega_m
    H_0_copy=H_0*(nc_km/nc_megaparsec)
    age_copy=age*nc_gyr
!    print*, age
    scalefactor=(omega_m/omega_lambda)**(1./3.)*&
         sinh(3./2.*H_0_copy*age_copy*sqrt(omega_lambda)) **(2./3.)
!    print*, "sf:",scalefactor
    z=1./scalefactor-1.0
!    print*, "z:",z
    agetoz=z    
  end function AGETOZ

end module M_cosmology

!===============================================================================
!                      ****END M_cosmology****                                
!===============================================================================

