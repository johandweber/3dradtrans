!module natural_constants
!  use data_types
!  implicit none
!  real(dp), parameter :: nc_lyman_edge=1.0973731569E+05_8*2.99792458E+10_8
!  real(dp), parameter :: nc_light =2.99792458E+10_8
!  real(dp), parameter :: nc_planck = 6.6260755E-27_8
!  real(dp), parameter :: nc_boltzmann = 1.380658E-16_8
!  real(dp), parameter::  nc_sol_lum = 3.826E+33_8
!  real(dp), parameter :: nc_sol_rad = 6.9599E+10_8
!  real(dp), parameter::  nc_parsec= 3.08567808E+18_8
!  real(dp), parameter::  nc_h_ion_cross  = 6.3E-18_8
!  real(dp), parameter::  nc_pi=3.141592653589793116_8
!  real(dp), parameter::  nc_sphere=4*nc_pi
!end module natural_constants

module M_elements
  use M_data_types
  implicit none
  integer, parameter:: HYDROGEN   =  1
  integer, parameter:: HELIUM     =  2
  integer, parameter:: CARBON     =  3
  integer, parameter:: NITROGEN   =  4
  integer, parameter:: OXYGEN     =  5
  integer, parameter:: NEON       =  6
  integer, parameter:: SULFUR     =  7 
  integer, parameter:: ARGON      =  8
  integer, parameter:: ELEMENTS_MAX=ARGON
end module M_elements

module M_mathematics
  use M_data_types
  use M_natural_constants
contains
!
  real(dp) function col_fs_cooling  &
       (TEMPERATURE, ELECTRONDENSITY, LAMBDA, YPSILON, OMEGA2, CRIT)
    IMPLICIT NONE
    real(dp) :: ELECTRONDENSITY, TEMPERATURE,LAMBDA, YPSILON, OMEGA2, CRIT
    col_fs_cooling=&
         ELECTRONDENSITY*8.629e-6/sqrt(temperature)*YPSILON/OMEGA2*&
         EXP(-(nc_planck*nc_light/(LAMBDA*1e-8))  / (nc_boltzmann*TEMPERATURE))/&
         ((ELECTRONDENSITY+CRIT)/CRIT)*(nc_planck*nc_light/(LAMBDA*1e-8))
  end function col_fs_cooling
!
  real(dp) function recombination_line (TEMPERATURE, ELECTRONDENSITY,A,B)
    implicit none
    real(dp)::temperature, electrondensity,a,b
    recombination_line=electrondensity*10**(a+b*log10(TEMPERATURE))
  end function recombination_line
!
  real(dp) function maskedvalue(inputvalue,testvalue,lowerlimit)
! only relevant values (above a threshold for the emission measure)
! should be plotted   
    implicit none
    real(dp):: inputvalue, testvalue, lowerlimit
!
    if (testvalue .lt. lowerlimit) then
       maskedvalue=-1E-30_dp              ! Very small negaive value
                                          ! The exact value does not have any
                                          ! meaning.
                                          
    else
       maskedvalue=inputvalue
    end if
  end function maskedvalue

end module M_mathematics



program nebular_statistics
  use M_data_types
  use M_natural_constants
  use M_elements
  use M_mathematics
  implicit none

  real(dp), dimension (:,:,:,:,:), allocatable :: ions
  real(dp), dimension (:,:,:) , allocatable    :: electrons, density,&
                                                  temperature
!
  real(dp), dimension(:,:,:), allocatable  :: emissivity_h_alpha
  real(dp) :: max_emissivity_h_alpha_x=0,& 
                      max_emissivity_h_alpha_y=0,&
                      max_emissivity_h_alpha_z=0   
!
  real(dp), dimension (:,:,:), allocatable ::  emissivity_h_beta
  real(dp) :: max_emissivity_h_beta_x=0,&
                      max_emissivity_h_beta_y=0,&
                      max_emissivity_h_beta_z=0   
!
  real(dp), dimension(:,:,:), allocatable  :: emissivity_h_gamma 
  real(dp) :: max_emissivity_h_gamma_x=0,&
                      max_emissivity_h_gamma_y=0,&
                      max_emissivity_h_gamma_z=0   
!
  real(dp), dimension(:,:,:), allocatable  ::  emissivity_h_delta
  real(dp) :: max_emissivity_h_delta_x=0,&
                      max_emissivity_h_delta_y=0,&
                      max_emissivity_h_delta_z=0   
!
  real(dp), dimension(:,:,:), allocatable  ::  emissivity_pasch_alpha
  real(dp) :: max_emissivity_pasch_alpha_x=0,&
                      max_emissivity_pasch_alpha_y=0,&
                      max_emissivity_pasch_alpha_z=0   
!
  real(dp), dimension(:,:,:), allocatable  ::  emissivity_pasch_beta
  real(dp) :: max_emissivity_pasch_beta_x=0,&
                      max_emissivity_pasch_beta_y=0,&
                      max_emissivity_pasch_beta_z=0   
!
  real(dp), dimension(:,:,:), allocatable  :: emissivity_brack_alpha
  real(dp) :: max_emissivity_brack_alpha_x=0,&
                      max_emissivity_brack_alpha_y=0,&
                      max_emissivity_brack_alpha_z=0   
!
 real(dp), dimension(:,:,:), allocatable  ::  emissivity_brack_beta
  real(dp) :: max_emissivity_brack_beta_x=0,&
                      max_emissivity_brack_beta_y=0,&
                      max_emissivity_brack_beta_z=0   
!  
  real(dp), dimension(:,:,:),allocatable ::  emissivity_hei_optical
  real(dp) :: max_emissivity_hei_optical_x=0,&
                      max_emissivity_hei_optical_y=0,&
                      max_emissivity_hei_optical_z=0
!
  real(dp), dimension(:,:,:),allocatable ::  emissivity_cii_1335
  real(dp) :: max_emissivity_cii_1335_x=0,&
                      max_emissivity_cii_1335_y=0,&
                      max_emissivity_cii_1335_z=0
!
  real(dp), dimension(:,:,:),allocatable ::  emissivity_cii_2325
  real(dp) :: max_emissivity_cii_2325_x=0,&
                      max_emissivity_cii_2325_y=0,&
                      max_emissivity_cii_2325_z=0
!  
  real(dp), dimension(:,:,:),allocatable ::  emissivity_ciii_1907
  real(dp) :: max_emissivity_ciii_1907_x=0,&
                      max_emissivity_ciii_1907_y=0,&
                      max_emissivity_ciii_1907_z=0
!  
  real(dp), dimension(:,:,:),allocatable ::  emissivity_nii_optical
  real(dp) :: max_emissivity_nii_optical_x=0,&
                      max_emissivity_nii_optical_y=0,&
                      max_emissivity_nii_optical_z=0   
!
  real(dp), dimension(:,:,:),allocatable ::  emissivity_nii_auroral
  real(dp) :: max_emissivity_nii_auroral_x=0,&
                      max_emissivity_nii_auroral_y=0,&
                      max_emissivity_nii_auroral_z=0   
!
  real(dp), dimension(:,:,:),allocatable ::  emissivity_nii_infrared122
  real(dp) :: max_emissivity_nii_infrared122_x=0,&
                      max_emissivity_nii_infrared122_y=0,&
                      max_emissivity_nii_infrared122_z=0   

!
  real(dp), dimension(:,:,:),allocatable ::  emissivity_niii_infrared57
  real(dp) :: max_emissivity_niii_infrared57_x=0,&
                      max_emissivity_niii_infrared57_y=0,&
                      max_emissivity_niii_infrared57_z=0   
!
  real(dp), dimension(:,:,:),allocatable ::   emissivity_oi_optical
  real(dp) :: max_emissivity_oi_optical_x=0,&
                      max_emissivity_oi_optical_y=0,&
                      max_emissivity_oi_optical_z=0   
!
  real(dp), dimension(:,:,:),allocatable ::   emissivity_oii_optical
  real(dp) :: max_emissivity_oii_optical_x=0,&
                      max_emissivity_oii_optical_y=0,&
                      max_emissivity_oii_optical_z=0   
!
  real(dp), dimension(:,:,:),allocatable ::   emissivity_oii_7320
  real(dp) :: max_emissivity_oii_7320_x=0,&
                      max_emissivity_oii_7320_y=0,&
                      max_emissivity_oii_7320_z=0   
!
  real(dp), dimension(:,:,:),allocatable :: emissivity_oiii_optical
  real(dp) :: max_emissivity_oiii_optical_x=0,&
                      max_emissivity_oiii_optical_y=0,&
                      max_emissivity_oiii_optical_z=0
!
  real(dp), dimension(:,:,:),allocatable :: emissivity_oiii_auroral
  real(dp) :: max_emissivity_oiii_auroral_x=0,&
                      max_emissivity_oiii_auroral_y=0,&
                      max_emissivity_oiii_auroral_z=0
!  
  real(dp), dimension(:,:,:), allocatable :: emissivity_oiii_infrared52
  real(dp) :: max_emissivity_oiii_infrared52_x=0,&
                      max_emissivity_oiii_infrared52_y=0,&
                      max_emissivity_oiii_infrared52_z=0
!
  real(dp), dimension(:,:,:), allocatable :: emissivity_oiii_infrared88
  real(dp) :: max_emissivity_oiii_infrared88_x=0,&
                      max_emissivity_oiii_infrared88_y=0,&
                      max_emissivity_oiii_infrared88_z=0
!
  real(dp), dimension(:,:,:), allocatable :: emissivity_oiv_infrared
  real(dp) :: max_emissivity_oiv_infrared_x=0,&
                      max_emissivity_oiv_infrared_y=0,&
                      max_emissivity_oiv_infrared_z=0
!
  real(dp), dimension(:,:,:), allocatable :: emissivity_neii_infrared13
  real(dp) :: max_emissivity_neii_infrared13_x=0,&
                      max_emissivity_neii_infrared13_y=0,&
                      max_emissivity_neii_infrared13_z=0
!
  real(dp), dimension(:,:,:), allocatable :: emissivity_neiii_optical
  real(dp) :: max_emissivity_neiii_optical_x=0,&
                      max_emissivity_neiii_optical_y=0,&
                      max_emissivity_neiii_optical_z=0  
  
  real(dp), dimension(:,:,:), allocatable :: emissivity_neiii_infrared16
  real(dp) :: max_emissivity_neiii_infrared16_x=0,&
                     max_emissivity_neiii_infrared16_y=0,&
                     max_emissivity_neiii_infrared16_z=0
!
  real(dp), dimension(:,:,:),allocatable :: emissivity_sii_4068
  real(dp) :: max_emissivity_sii_4068_x=0,&
                      max_emissivity_sii_4068_y=0,&
                      max_emissivity_sii_4068_z=0  
!  
  real(dp), dimension(:,:,:),allocatable :: emissivity_sii_6716
  real(dp) :: max_emissivity_sii_6716_x=0,&
                      max_emissivity_sii_6716_y=0,&
                      max_emissivity_sii_6716_z=0
!  
  real(dp), dimension(:,:,:),allocatable :: emissivity_siii_optical
  real(dp) :: max_emissivity_siii_optical_x=0,&
                      max_emissivity_siii_optical_y=0,&
                      max_emissivity_siii_optical_z=0
!
  real(dp), dimension(:,:,:),allocatable :: emissivity_siii_infrared19
  real(dp) :: max_emissivity_siii_infrared19_x=0,&
                      max_emissivity_siii_infrared19_y=0,&
                      max_emissivity_siii_infrared19_z=0
! 
  real(dp), dimension(:,:,:),allocatable :: emissivity_siii_infrared34
  real(dp) :: max_emissivity_siii_infrared34_x=0,&
                      max_emissivity_siii_infrared34_y=0,&
                      max_emissivity_siii_infrared34_z=0
!
  real(dp), dimension(:,:,:),allocatable :: emissivity_siv_infrared11
  real(dp) :: max_emissivity_siv_infrared11_x=0,&
                      max_emissivity_siv_infrared11_y=0,&
                      max_emissivity_siv_infrared11_z=0
! 
  real(dp) :: max_emissionmeasure_x, max_emissionmeasure_y, max_emissionmeasure_z
!
  real(dp), dimension(:,:,:),allocatable :: EMISSIVITY_OIII_RADIO
 
  real(dp)                             :: l_cell !length of a grid cell in pc

  integer, dimension (ELEMENTS_MAX,4) :: CELL_IONS=0

  real(dp):: nenpsum=0D0, tnenpsum=0D0

  real(dp):: mean_of_squares=0.
  real(dp):: square_of_mean=0.
  real(dp):: xHbeta, yHbeta, zHbeta
  real(dp):: em_x, em_y, em_z, em_min=1.0_dp
  integer                  :: number_of_arguments
  integer                  :: x, y, z
  integer                  :: x_max, y_max, z_max 
  character(len=3000)      :: trunkstring,l_cell_string
!
  number_of_arguments=command_argument_count()
  call get_command_argument(1,trunkstring) 
  call get_command_argument(2,l_cell_string)
  read(l_cell_string,*) l_cell
  call open_input_files
! 
  read(1,*) x_max, y_max, z_max
!
  call allocate_cell_arrays
!
  read(2,*)
!  read(3,*)

  ! ANALYZE HYDROGEN

  call read_grid_data

  call compute_weighted_temperatures

  call write_global_abundance_data

  call determine_most_abundant_stages

  call write_most_abundant_stages
 
  call calc_emissivities
  
  call write_global_emission_data
 
!write spectral_maps 

! ... in x-direction

 call open_x_output_files

 call write_maps_x
 
 call close_output_files
 
 
!... in y direction

 call open_y_output_files

 call write_maps_y

 call close_output_files

!... in z-direction

call open_z_output_files

call write_maps_z

call close_output_files
! 
!
 call close_input_files

 stop

 contains 

   subroutine open_input_files
     implicit none
     open(unit=1, file=trim(trunkstring)//'.temp.txt')
     open(unit=2, file=trim(trunkstring)//'.e.txt')
     open(unit=100+HYDROGEN    , file=trim(trunkstring)//'.H.txt')
     open(unit=100+HELIUM      , file=trim(trunkstring)//'.He.txt')
     open(unit=100+CARBON      , file=trim(trunkstring)//'.C.txt')
     open(unit=100+NITROGEN    , file=trim(trunkstring)//'.N.txt')
     open(unit=100+OXYGEN      , file=trim(trunkstring)//'.O.txt')
     open(unit=100+NEON        , file=trim(trunkstring)//'.Ne.txt')
     open(unit=100+SULFUR      , file=trim(trunkstring)//'.S.txt')
     open(unit=100+ARGON       , file=trim(trunkstring)//'.Ar.txt')     
   end subroutine open_input_files

   subroutine read_grid_data
     implicit none
!
     read(100+HYDROGEN,*)
     read(100+HELIUM,*)
     read(100+CARBON,*)
     read(100+NITROGEN,*)
     read(100+OXYGEN,*)
     read(100+NEON,*)
     read(100+SULFUR,*)
     read(100+ARGON,*)
!     
     do z=1,z_max
        do y=1, y_max
           do x=1, x_max
              read(1,*) temperature(x, y, z)
              read(2,*) electrons(x,y, z)
              !           read(3,*) density(x, y, z)
              read(100+HYDROGEN,*) ions(x, y, z,HYDROGEN,1),&
                   ions(x, y, z,HYDROGEN,2),&
                   density(x, y, z)
              read(100+HELIUM,*) ions(x, y, z,HELIUM,1),&
                   ions(x, y, z,HELIUM,2),&
                   ions(x, y, z,HELIUM,3)
              !           print*, x, y, Z_counter, ions(x, y, z,HELIUM,3) 
              read(100+CARBON,*) ions(x, y, z,CARBON,1),&
                   ions(x, y, z,CARBON,2),&
                   ions(x, y, z,CARBON,3),&
                   ions(x, y, z,CARBON,4)
              read(100+NITROGEN,*) ions(x, y, z,NITROGEN,1),&
                   ions(x, y, z,NITROGEN,2),&
                   ions(x, y, z,NITROGEN,3),&
                   ions(x, y, z,NITROGEN,4)
              read(100+OXYGEN,*) ions(x, y, z,OXYGEN,1),&
                   ions(x, y, z,OXYGEN,2),&
                   ions(x, y, z,OXYGEN,3),&
                   ions(x, y, z,OXYGEN,4)
              read(100+NEON,*) ions(x, y, z,NEON,1),&
                   ions(x, y, z,NEON,2),&
                   ions(x, y, z,NEON,3),&
                   ions(x, y, z,NEON,4)
              read(100+SULFUR,*) ions(x, y, z,SULFUR,1),&
                   ions(x, y, z,SULFUR,2),&
                   ions(x, y, z,SULFUR,3),&
                   ions(x, y, z,SULFUR,4)
              read(100+ARGON,*) ions(x, y, z,ARGON,1),&
                   ions(x, y, z,ARGON,2),&
                   ions(x, y, z,ARGON,3),&
                   ions(x, y, z,ARGON,4)           
           end do
        end do
     end do
   end subroutine read_grid_data

   subroutine compute_weighted_temperatures
     implicit none
!
     do z=1,z_max
        do y=1,y_max
           do x=1,x_max
              nenpsum=nenpsum+electrons(x,y, z)*&
                   ions(x, y, z,HYDROGEN,2)
              tnenpsum=tnenpsum+temperature(x, y, z)*&
                   electrons(x, y, z)*&
                   ions(x, y, z, HYDROGEN,2)
           end do
        end do
     end do     
   end subroutine compute_weighted_temperatures


   subroutine determine_most_abundant_stages
     implicit none
     do z=1,z_max
        do y=1,y_max
           do x=1,x_max
              CELL_IONS(HYDROGEN,MAXLOC(ions(x,y,z,HYDROGEN,:)))=&
                   CELL_IONS(HYDROGEN,MAXLOC(ions(x,y,z,HYDROGEN,:)))+1
              CELL_IONS(HELIUM,MAXLOC(ions(x,y,z,HELIUM,:)))=&
                   CELL_IONS(HELIUM,MAXLOC(ions(x,y,z,HELIUM,:)))+1
              CELL_IONS(CARBON,MAXLOC(ions(x,y,z,CARBON,:)))=&
                   CELL_IONS(CARBON,MAXLOC(ions(x,y,z,CARBON,:)))+1
              CELL_IONS(NITROGEN,MAXLOC(ions(x,y,z,NITROGEN,:)))=&
                   CELL_IONS(NITROGEN,MAXLOC(ions(x,y,z,NITROGEN,:)))+1
              CELL_IONS(OXYGEN,MAXLOC(ions(x,y,z,OXYGEN,:)))=&
                   CELL_IONS(OXYGEN,MAXLOC(ions(x,y,z,OXYGEN,:)))+1
              CELL_IONS(NEON,MAXLOC(ions(x,y,z,NEON,:)))=&
                   CELL_IONS(NEON,MAXLOC(ions(x,y,z,NEON,:)))+1
              CELL_IONS(SULFUR,MAXLOC(ions(x,y,z,SULFUR,:)))=&
                   CELL_IONS(SULFUR,MAXLOC(ions(x,y,z,SULFUR,:)))+1
              CELL_IONS(ARGON,MAXLOC(ions(x,y,z,ARGON,:)))=&
                   CELL_IONS(ARGON,MAXLOC(ions(x,y,z,ARGON,:)))+1           
           end do
        end do
     end do
   end subroutine determine_most_abundant_stages

   subroutine write_most_abundant_stages
     implicit none
!
     write(*,*)
     write(*,*)
     write(*,*) 'Fraction of cells with the most abundant ionization stage:'
     write(*,*) '=========================================================='
     write(*,'(A17, 2ES15.5)') 'H',&
          CELL_IONS(HYDROGEN,1)*1.0/sum(CELL_IONS(HYDROGEN,:)),&
          CELL_IONS(HYDROGEN,2)*1.0/sum(CELL_IONS(HYDROGEN,:))
     write(*,'(A17, 3ES15.5)') 'He',&
          CELL_IONS(HELIUM,1)*1.0/sum(CELL_IONS(HELIUM,:)),&
          CELL_IONS(HELIUM,2)*1.0/sum(CELL_IONS(HELIUM,:)),&
          CELL_IONS(HELIUM,3)*1.0/sum(CELL_IONS(HELIUM,:))
     write(*,'(A17, 4ES15.5)') 'C',&
          CELL_IONS(CARBON,1)*1.0/sum(CELL_IONS(CARBON,:)),&
          CELL_IONS(CARBON,2)*1.0/sum(CELL_IONS(CARBON,:)),&
          CELL_IONS(CARBON,3)*1.0/sum(CELL_IONS(CARBON,:)),&
          CELL_IONS(CARBON,4)*1.0/sum(CELL_IONS(CARBON,:))
     write(*,'(A17, 4ES15.5)') 'N',&
          CELL_IONS(NITROGEN,1)*1.0/sum(CELL_IONS(NITROGEN,:)),&
          CELL_IONS(NITROGEN,2)*1.0/sum(CELL_IONS(NITROGEN,:)),&
          CELL_IONS(NITROGEN,3)*1.0/sum(CELL_IONS(NITROGEN,:)),&
          CELL_IONS(NITROGEN,4)*1.0/sum(CELL_IONS(NITROGEN,:))
     write(*,'(A17, 4ES15.5)') 'O',&
          CELL_IONS(OXYGEN,1)*1.0/sum(CELL_IONS(OXYGEN,:)),&
          CELL_IONS(OXYGEN,2)*1.0/sum(CELL_IONS(OXYGEN,:)),&
          CELL_IONS(OXYGEN,3)*1.0/sum(CELL_IONS(OXYGEN,:)),&
          CELL_IONS(OXYGEN,4)*1.0/sum(CELL_IONS(OXYGEN,:))
     write(*,'(A17, 4ES15.5)') 'Ne',&
          CELL_IONS(NEON,1)*1.0/sum(CELL_IONS(NEON,:)),&
          CELL_IONS(NEON,2)*1.0/sum(CELL_IONS(NEON,:)),&
          CELL_IONS(NEON,3)*1.0/sum(CELL_IONS(NEON,:)),&
          CELL_IONS(NEON,4)*1.0/sum(CELL_IONS(NEON,:))
     write(*,'(A17, 4ES15.5)') 'S',&
          CELL_IONS(SULFUR,1)*1.0/sum(CELL_IONS(SULFUR,:)),&
          CELL_IONS(SULFUR,2)*1.0/sum(CELL_IONS(SULFUR,:)),&
          CELL_IONS(SULFUR,3)*1.0/sum(CELL_IONS(SULFUR,:)),&
          CELL_IONS(SULFUR,4)*1.0/sum(CELL_IONS(SULFUR,:))
     write(*,'(A17, 4ES15.5)') 'Ar',&
          CELL_IONS(ARGON,1)*1.0/sum(CELL_IONS(ARGON,:)),&
          CELL_IONS(ARGON,2)*1.0/sum(CELL_IONS(ARGON,:)),&
          CELL_IONS(ARGON,3)*1.0/sum(CELL_IONS(ARGON,:)),&
          CELL_IONS(ARGON,4)*1.0/sum(CELL_IONS(ARGON,:))  
     write(*,*)
   end subroutine write_most_abundant_stages

   subroutine calc_line_emssivities
 
  end subroutine calc_line_emssivities
     
   subroutine close_input_files
     implicit none
     close(1)
     close(2)
     ! close(3)
     close(100+HYDROGEN)
     close(100+HELIUM)
     close(100+CARBON)
     close(100+NITROGEN)
     close(100+OXYGEN)
     close(100+NEON)
     close(100+SULFUR)
     close(100+ARGON)     
   end subroutine close_input_files

   subroutine allocate_cell_arrays
     implicit none
     allocate (ions(x_max, y_max, z_max,ELEMENTS_MAX,4))
     ions=0
     allocate(electrons(x_max,y_max,z_max))
     electrons=0
     allocate(density(x_max, y_max, z_max))
     density=0
     allocate(temperature(x_max, y_max, z_max))
     temperature=0
     allocate(emissivity_h_alpha(x_max, y_max, z_max))
     emissivity_h_alpha=0
     allocate(emissivity_h_beta(x_max, y_max, z_max))
     emissivity_h_beta=0
     allocate(emissivity_h_gamma(x_max, y_max, z_max))
     emissivity_h_gamma=0
     allocate(emissivity_h_delta(x_max, y_max, z_max))
     emissivity_h_delta=0
     allocate(emissivity_pasch_alpha(x_max, y_max, z_max))
     emissivity_pasch_alpha=0
     allocate(emissivity_pasch_beta(x_max, y_max, z_max))
     emissivity_pasch_beta=0
     allocate(emissivity_brack_alpha(x_max, y_max, z_max))
     emissivity_brack_alpha=0
     allocate(emissivity_brack_beta(x_max, y_max, z_max))
     emissivity_brack_beta=0
     allocate(emissivity_hei_optical(x_max, y_max, z_max))
     emissivity_hei_optical=0
     allocate(emissivity_cii_1335(x_max,y_max,z_max))
     emissivity_cii_1335=0
     allocate(emissivity_cii_2325(x_max,y_max,z_max))
     emissivity_cii_2325=0
     allocate(emissivity_ciii_1907(x_max,y_max,z_max))
     emissivity_ciii_1907=0     
     allocate(emissivity_nii_optical(x_max,y_max, z_max))
     emissivity_nii_optical=0
     allocate(emissivity_nii_auroral(x_max,y_max, z_max))
     emissivity_nii_auroral=0
     allocate(emissivity_nii_infrared122(x_max,y_max, z_max))
     emissivity_nii_infrared122=0
     allocate(emissivity_niii_infrared57(x_max,y_max, z_max))
     emissivity_niii_infrared57=0
     allocate(emissivity_oi_optical(x_max,y_max, z_max))
     emissivity_oi_optical=0     
     allocate(emissivity_oii_optical(x_max,y_max, z_max))
     emissivity_oii_optical=0
     allocate(emissivity_oii_7320(x_max,y_max, z_max))
     emissivity_oii_7320=0     
     allocate(emissivity_oiii_optical(x_max,y_max, z_max))
     emissivity_oiii_optical=0
     allocate(emissivity_oiii_auroral(x_max,y_max, z_max))
     emissivity_oiii_auroral=0
     allocate(emissivity_oiii_infrared52(x_max,y_max, z_max))
     emissivity_oiii_infrared52=0     
     allocate(emissivity_oiii_infrared88(x_max,y_max, z_max))
     emissivity_oiii_infrared88=0
     allocate(emissivity_oiii_radio(x_max,y_max, z_max))
     emissivity_oiii_radio=0
     allocate(emissivity_oiv_infrared(x_max, y_max, z_max))
     emissivity_oiv_infrared=0
     allocate(emissivity_neii_infrared13(x_max,y_max, z_max))
     emissivity_neii_infrared13=0     
     allocate(emissivity_neiii_optical(x_max,y_max, z_max))
     emissivity_neiii_optical=0
     allocate(emissivity_neiii_infrared16(x_max,y_max, z_max))
     emissivity_neiii_infrared16=0
     allocate(emissivity_sii_4068(x_max,y_max, z_max))
     emissivity_sii_4068=0
     allocate(emissivity_sii_6716(x_max,y_max, z_max))
     emissivity_sii_6716=0
     allocate(emissivity_siii_optical(x_max,y_max, z_max))
     emissivity_siii_optical=0
     allocate(emissivity_siii_infrared19(x_max,y_max, z_max))
     emissivity_siii_infrared19=0
     allocate(emissivity_siii_infrared34(x_max,y_max, z_max))
     emissivity_siii_infrared34=0
     allocate(emissivity_siv_infrared11(x_max,y_max, z_max))
     emissivity_siv_infrared11=0

   end subroutine allocate_cell_arrays


   subroutine open_x_output_files
     implicit none
     open(unit=1000+100*HYDROGEN+1,&
          file=trim(trunkstring)//".H_alpha_x.txt")
     write(1000+100*HYDROGEN+1,*) 1, y_max, z_max
     open(unit=1000+100*HYDROGEN+2,&
          file=trim(trunkstring)//".H_beta_x.txt")
     write(1000+100*HYDROGEN+2,*) 1, y_max, z_max
     open(unit=1000+100*HYDROGEN+3,&
          file=trim(trunkstring)//".H_gamma_x.txt")
     write(1000+100*HYDROGEN+3,*) 1, y_max, z_max
     open(unit=1000+100*HYDROGEN+4,&
          file=trim(trunkstring)//".H_delta_x.txt")
     write(1000+100*HYDROGEN+4,*) 1, y_max, z_max
     open(unit=1000+100*HYDROGEN+5,&
          file=trim(trunkstring)//".Pasch_alpha_x.txt")
     write(1000+100*HYDROGEN+5,*) 1, y_max, z_max
     open(unit=1000+100*HYDROGEN+6,&
          file=trim(trunkstring)//".Pasch_beta_x.txt")
     write(1000+100*HYDROGEN+6,*) 1, y_max, z_max
     open(unit=1000+100*HYDROGEN+7,&
          file=trim(trunkstring)//".Brack_alpha_x.txt")
     write(1000+100*HYDROGEN+7,*) 1, y_max, z_max
     open(unit=1000+100*HYDROGEN+8,&
          file=trim(trunkstring)//".Brack_beta_x.txt")
     write(1000+100*HYDROGEN+8,*) 1, y_max, z_max
!
     open(unit=1000+100*HELIUM+1,&
          file=trim(trunkstring)//".HeI-optical_x.txt")
     write(1000+100*HELIUM+1,*) 1, y_max, z_max     
!
     open(unit=1000+100*CARBON+1,&
          file=trim(trunkstring)//".CII-1335_x.txt")
     write(1000+100*CARBON+1,*) 1, y_max, z_max 
     open(unit=1000+100*CARBON+2,&
          file=trim(trunkstring)//".CII-2325_x.txt")
     write(1000+100*CARBON+2,*) 1, y_max, z_max
     open(unit=1000+100*CARBON+3,&
          file=trim(trunkstring)//".CIII-1907_x.txt")
     write(1000+100*CARBON+3,*) 1, y_max, z_max     
!
     open(unit=1000+100*NITROGEN+1,&
          file=trim(trunkstring)//".NII-optical_x.txt")
     write(1000+100*NITROGEN+1,*) 1, y_max, z_max
     open(unit=1000+100*NITROGEN+2,&
          file=trim(trunkstring)//".NII-auroral_x.txt")
     write(1000+100*NITROGEN+2,*) 1, y_max, z_max
     open(unit=1000+100*NITROGEN+3,&
          file=trim(trunkstring)//".NII-infrared122_x.txt")
     write(1000+100*NITROGEN+3,*) 1, y_max, z_max
     open(unit=1000+100*NITROGEN+4,&
          file=trim(trunkstring)//".NIII-infrared57_x.txt")
     write(1000+100*NITROGEN+4,*) 1, y_max, z_max
!
     open(unit=1000+100*OXYGEN+1,&
          file=trim(trunkstring)//".OI-optical_x.txt")
     write(1000+100*OXYGEN+1,*) 1, y_max, z_max
     open(unit=1000+100*OXYGEN+2,&
          file=trim(trunkstring)//".OII-optical_x.txt")
     write(1000+100*OXYGEN+2,*) 1, y_max, z_max
     open(unit=1000+100*OXYGEN+3,&
          file=trim(trunkstring)//".OII-7320_x.txt")
     write(1000+100*OXYGEN+3,*) 1, y_max, z_max
     open(unit=1000+100*OXYGEN+4,&
          file=trim(trunkstring)//".OIII-optical_x.txt")
     write(1000+100*OXYGEN+4, *) 1, y_max, z_max
     open(unit=1000+100*OXYGEN+5,&
          file=trim(trunkstring)//".OIII-auroral_x.txt")
     write(1000+100*OXYGEN+5, *) 1, y_max, z_max
     open(unit=1000+100*OXYGEN+6,&
          file=trim(trunkstring)//".OIII-infrared52_x.txt")
     write(1000+100*OXYGEN+6, *) 1, y_max, z_max
     open(unit=1000+100*OXYGEN+7,&
          file=trim(trunkstring)//".OIII-infrared88_x.txt")
     write(1000+100*OXYGEN+7, *) 1, y_max, z_max
     open(unit=1000+100*OXYGEN+8,&
          file=trim(trunkstring)//".OIV-infrared_x.txt")
     write(1000+100*OXYGEN+8, *) 1, y_max, z_max
!
     open(unit=1000+100*NEON+1,&
          file=trim(trunkstring)//".NeII-infrared13_x.txt")
     write(1000+100*NEON+1, *) 1, y_max, z_max
     open(unit=1000+100*NEON+2,&
          file=trim(trunkstring)//".NeIII-optical_x.txt")
     write(1000+100*NEON+2, *) 1, y_max, z_max
     open(unit=1000+100*NEON+3,&
          file=trim(trunkstring)//".NeIII-infrared16_x.txt")
     write(1000+100*NEON+3, *) 1, y_max, z_max
!
     open(unit=1000+100*SULFUR+1,&
          file=trim(trunkstring)//".SII_4068_x.txt")
     write(1000+100*SULFUR+1, *) 1, y_max, z_max
     open(unit=1000+100*SULFUR+2,&
          file=trim(trunkstring)//".SII_6716_x.txt")
     write(1000+100*SULFUR+2, *) 1, y_max, z_max
     open(unit=1000+100*SULFUR+3,&
          file=trim(trunkstring)//".SIII-optical_x.txt")
     write(1000+100*SULFUR+3, *) 1, y_max, z_max
     open(unit=1000+100*SULFUR+4,&
          file=trim(trunkstring)//".SIII-infrared19_x.txt")
     write(1000+100*SULFUR+4, *) 1, y_max, z_max
     open(unit=1000+100*SULFUR+5,&
          file=trim(trunkstring)//".SIII-infrared34_x.txt")
     write(1000+100*SULFUR+5, *) 1, y_max, z_max
     open(unit=1000+100*SULFUR+6,&
          file=trim(trunkstring)//".SIV-infrared11_x.txt")
     write(1000+100*SULFUR+6, *) 1, y_max, z_max
!
     open(unit=10000, file=trim(trunkstring)//".emissionmeasure_x.txt")
     write(10000, *) 1, y_max, z_max     
   end subroutine open_x_output_files
   
   subroutine open_y_output_files
     implicit none
     open(unit=1000+100*HYDROGEN+1,&
          file=trim(trunkstring)//".H_alpha_y.txt")
     write(1000+100*HYDROGEN+1,*) x_max, 1, z_max
     open(unit=1000+100*HYDROGEN+2,&
          file=trim(trunkstring)//".H_beta_y.txt")
     write(1000+100*HYDROGEN+2,*) x_max, 1, z_max
     open(unit=1000+100*HYDROGEN+3,&
          file=trim(trunkstring)//".H_gamma_y.txt")
     write(1000+100*HYDROGEN+3,*) x_max, 1, z_max
     open(unit=1000+100*HYDROGEN+4,&
          file=trim(trunkstring)//".H_delta_y.txt")
     write(1000+100*HYDROGEN+4,*) x_max, 1, z_max
     open(unit=1000+100*HYDROGEN+5,&
          file=trim(trunkstring)//".Pasch_alpha_y.txt")
     write(1000+100*HYDROGEN+5,*) x_max, 1, z_max
     open(unit=1000+100*HYDROGEN+6,&
          file=trim(trunkstring)//".Pasch_beta_y.txt")
     write(1000+100*HYDROGEN+6,*) x_max, 1, z_max
     open(unit=1000+100*HYDROGEN+7,&
          file=trim(trunkstring)//".Brack_alpha_y.txt")
     write(1000+100*HYDROGEN+7,*) x_max, 1, z_max
     open(unit=1000+100*HYDROGEN+8,&
          file=trim(trunkstring)//".Brack_beta_y.txt")
     write(1000+100*HYDROGEN+8,*) x_max, 1, z_max
!
     open(unit=1000+100*HELIUM+1,&
          file=trim(trunkstring)//".HeI-optical_y.txt")
     write(1000+100*HELIUM+1,*) x_max, 1, z_max     
!
     open(unit=1000+100*CARBON+1,&
          file=trim(trunkstring)//".CII-1335_y.txt")
     write(1000+100*CARBON+1,*) x_max, 1, z_max 
     open(unit=1000+100*CARBON+2,&
          file=trim(trunkstring)//".CII-2325_y.txt")
     write(1000+100*CARBON+2,*) x_max, 1, z_max
     open(unit=1000+100*CARBON+3,&
          file=trim(trunkstring)//".CIII-1907_y.txt")
     write(1000+100*CARBON+3,*) x_max, 1, z_max     
!
     open(unit=1000+100*NITROGEN+1,&
          file=trim(trunkstring)//".NII-optical_y.txt")
     write(1000+100*NITROGEN+1,*) x_max, 1, z_max
     open(unit=1000+100*NITROGEN+2,&
          file=trim(trunkstring)//".NII-auroral_y.txt")
     write(1000+100*NITROGEN+2,*) x_max, 1, z_max
     open(unit=1000+100*NITROGEN+3,&
          file=trim(trunkstring)//".NII-infrared122_y.txt")
     write(1000+100*NITROGEN+3,*) x_max, 1, z_max
     open(unit=1000+100*NITROGEN+4,&
          file=trim(trunkstring)//".NIII-infrared57_y.txt")
     write(1000+100*NITROGEN+4,*) x_max, 1, z_max
     
!
     open(unit=1000+100*OXYGEN+1,&
          file=trim(trunkstring)//".OI-optical_y.txt")
     write(1000+100*OXYGEN+1,*) x_max, 1, z_max
     open(unit=1000+100*OXYGEN+2,&
          file=trim(trunkstring)//".OII-optical_y.txt")
     write(1000+100*OXYGEN+2,*) x_max, 1, z_max
     open(unit=1000+100*OXYGEN+3,&
          file=trim(trunkstring)//".OII-7320_y.txt")
     write(1000+100*OXYGEN+3,*) x_max, 1, z_max
     open(unit=1000+100*OXYGEN+4,&
          file=trim(trunkstring)//".OIII-optical_y.txt")
     write(1000+100*OXYGEN+4, *) x_max, 1, z_max
     open(unit=1000+100*OXYGEN+5,&
          file=trim(trunkstring)//".OIII-auroral_y.txt")
     write(1000+100*OXYGEN+5, *) x_max, 1, z_max
     open(unit=1000+100*OXYGEN+6,&
          file=trim(trunkstring)//".OIII-infrared52_y.txt")
     write(1000+100*OXYGEN+6, *) x_max, 1, z_max
     open(unit=1000+100*OXYGEN+7,&
          file=trim(trunkstring)//".OIII-infrared88_y.txt")
     write(1000+100*OXYGEN+7, *) x_max, 1, z_max
     open(unit=1000+100*OXYGEN+8,&
          file=trim(trunkstring)//".OIV-infrared_y.txt")
     write(1000+100*OXYGEN+8, *) x_max, 1, z_max
!
     open(unit=1000+100*NEON+1,&
          file=trim(trunkstring)//".NeII-infrared13_y.txt")
     write(1000+100*NEON+1, *) x_max, 1, z_max
     open(unit=1000+100*NEON+2,&
          file=trim(trunkstring)//".NeIII-optical_y.txt")
     write(1000+100*NEON+2, *) x_max, 1, z_max
     open(unit=1000+100*NEON+3,&
          file=trim(trunkstring)//".NeIII-infrared16_y.txt")
     write(1000+100*NEON+3, *) x_max, 1, z_max
!
     open(unit=1000+100*SULFUR+1,&
          file=trim(trunkstring)//".SII_4068_y.txt")
     write(1000+100*SULFUR+1, *) x_max, 1, z_max
     open(unit=1000+100*SULFUR+2,&
          file=trim(trunkstring)//".SII_6716_y.txt")
     write(1000+100*SULFUR+2, *) x_max, 1, z_max
     open(unit=1000+100*SULFUR+3,&
          file=trim(trunkstring)//".SIII-optical_y.txt")
     write(1000+100*SULFUR+3, *) x_max, 1, z_max
     open(unit=1000+100*SULFUR+4,&
          file=trim(trunkstring)//".SIII-infrared19_y.txt")
     write(1000+100*SULFUR+4, *) x_max, 1, z_max
     open(unit=1000+100*SULFUR+5,&
          file=trim(trunkstring)//".SIII-infrared34_y.txt")
     write(1000+100*SULFUR+5, *) x_max, 1, z_max
     open(unit=1000+100*SULFUR+6,&
          file=trim(trunkstring)//".SIV-infrared11_y.txt")
     write(1000+100*SULFUR+6, *) x_max, 1, z_max
!
     open(unit=10000, file=trim(trunkstring)//".emissionmeasure_y.txt")
     write(10000, *) x_max, 1, z_max     
   end subroutine open_y_output_files

   subroutine open_z_output_files
     implicit none
     open(unit=1000+100*HYDROGEN+1,&
          file=trim(trunkstring)//".H_alpha_z.txt")
     write(1000+100*HYDROGEN+1,*) x_max, y_max, 1
     open(unit=1000+100*HYDROGEN+2,&
          file=trim(trunkstring)//".H_beta_z.txt")
     write(1000+100*HYDROGEN+2,*) x_max, y_max, 1
     open(unit=1000+100*HYDROGEN+3,&
          file=trim(trunkstring)//".H_gamma_z.txt")
     write(1000+100*HYDROGEN+3,*) x_max, y_max,1
     open(unit=1000+100*HYDROGEN+4,&
          file=trim(trunkstring)//".H_delta_z.txt")
     write(1000+100*HYDROGEN+4,*) x_max, y_max, 1
     open(unit=1000+100*HYDROGEN+5,&
          file=trim(trunkstring)//".Pasch_alpha_z.txt")
     write(1000+100*HYDROGEN+5,*) x_max, y_max, 1
     open(unit=1000+100*HYDROGEN+6,&
          file=trim(trunkstring)//".Pasch_beta_z.txt")
     write(1000+100*HYDROGEN+6,*) x_max, y_max, 1
     open(unit=1000+100*HYDROGEN+7,&
          file=trim(trunkstring)//".Brack_alpha_z.txt")
     write(1000+100*HYDROGEN+7,*) x_max, y_max, 1
     open(unit=1000+100*HYDROGEN+8,&
          file=trim(trunkstring)//".Brack_beta_z.txt")
     write(1000+100*HYDROGEN+8,*) x_max, y_max, 1
!
     open(unit=1000+100*HELIUM+1,&
          file=trim(trunkstring)//".HeI-optical_z.txt")
     write(1000+100*HELIUM+1,*) x_max, y_max, 1     
!
     open(unit=1000+100*CARBON+1,&
          file=trim(trunkstring)//".CII-1335_z.txt")
     write(1000+100*CARBON+1,*) x_max, y_max, 1 
     open(unit=1000+100*CARBON+2,&
          file=trim(trunkstring)//".CII-2325_z.txt")
     write(1000+100*CARBON+2,*) x_max, y_max, 1
     open(unit=1000+100*CARBON+3,&
          file=trim(trunkstring)//".CIII-1907_z.txt")
     write(1000+100*CARBON+3,*) x_max, y_max, 1     
!
     open(unit=1000+100*NITROGEN+1,&
          file=trim(trunkstring)//".NII-optical_z.txt")
     write(1000+100*NITROGEN+1,*) x_max, y_max,1
     open(unit=1000+100*NITROGEN+2,&
          file=trim(trunkstring)//".NII-auroral_z.txt")
     write(1000+100*NITROGEN+2,*) x_max, y_max,1
     open(unit=1000+100*NITROGEN+3,&
          file=trim(trunkstring)//".NII-infrared122_z.txt")
     write(1000+100*NITROGEN+3,*) x_max, y_max, 1
     open(unit=1000+100*NITROGEN+4,&
          file=trim(trunkstring)//".NIII-infrared57_z.txt")
     write(1000+100*NITROGEN+4,*) x_max, y_max, 1
!
     open(unit=1000+100*OXYGEN+1,&
          file=trim(trunkstring)//".OI-optical_z.txt")
     write(1000+100*OXYGEN+1,*) x_max, y_max, 1
     open(unit=1000+100*OXYGEN+2,&
          file=trim(trunkstring)//".OII-optical_z.txt")
     write(1000+100*OXYGEN+2,*) x_max, y_max, 1
     open(unit=1000+100*OXYGEN+3,&
          file=trim(trunkstring)//".OII-7320_z.txt")
     write(1000+100*OXYGEN+3,*) x_max, y_max, 1
     open(unit=1000+100*OXYGEN+4,&
          file=trim(trunkstring)//".OIII-optical_z.txt")
     write(1000+100*OXYGEN+4, *) x_max, y_max, 1
     open(unit=1000+100*OXYGEN+5,&
          file=trim(trunkstring)//".OIII-auroral_z.txt")
     write(1000+100*OXYGEN+5, *) x_max, y_max, 1
     open(unit=1000+100*OXYGEN+6,&
          file=trim(trunkstring)//".OIII-infrared52_z.txt")
     write(1000+100*OXYGEN+6, *) x_max, y_max, 1
     open(unit=1000+100*OXYGEN+7,&
          file=trim(trunkstring)//".OIII-infrared88_z.txt")
     write(1000+100*OXYGEN+7, *) x_max, y_max, 1
     open(unit=1000+100*OXYGEN+8,&
          file=trim(trunkstring)//".OIV-infrared_z.txt")
     write(1000+100*OXYGEN+8, *) x_max, y_max, 1
!
     open(unit=1000+100*NEON+1,&
          file=trim(trunkstring)//".NeII-infrared13_z.txt")
     write(1000+100*NEON+1, *) x_max, y_max, 1
     open(unit=1000+100*NEON+2,&
          file=trim(trunkstring)//".NeIII-optical_z.txt")
     write(1000+100*NEON+2, *) x_max, y_max, 1
     open(unit=1000+100*NEON+3,&
          file=trim(trunkstring)//".NeIII-infrared16_z.txt")
     write(1000+100*NEON+3, *) x_max, y_max, 1
!
     open(unit=1000+100*SULFUR+1,&
          file=trim(trunkstring)//".SII_4068_z.txt")
     write(1000+100*SULFUR+1, *) x_max, y_max, 1
     open(unit=1000+100*SULFUR+2,&
          file=trim(trunkstring)//".SII_6716_z.txt")
     write(1000+100*SULFUR+2, *) x_max, y_max, 1
     open(unit=1000+100*SULFUR+3,&
          file=trim(trunkstring)//".SIII-optical_z.txt")
     write(1000+100*SULFUR+3, *) x_max, y_max, 1
     open(unit=1000+100*SULFUR+4,&
          file=trim(trunkstring)//".SIII-infrared19_z.txt")
     write(1000+100*SULFUR+4, *) x_max, y_max, 1
     open(unit=1000+100*SULFUR+5,&
          file=trim(trunkstring)//".SIII-infrared34_z.txt")
     write(1000+100*SULFUR+5, *) x_max, y_max, 1
     open(unit=1000+100*SULFUR+6,&
          file=trim(trunkstring)//".SIV-infrared11_z.txt")
     write(1000+100*SULFUR+6, *) x_max, y_max, 1
!
     open(unit=10000, file=trim(trunkstring)//".emissionmeasure_z.txt")
     write(10000, *) x_max, y_max,1     
   end subroutine open_z_output_files

   subroutine close_output_files
     implicit none
     close(1000+100*HYDROGEN+1)
     close(1000+100*HYDROGEN+2)
     close(1000+100*HYDROGEN+3)
     close(1000+100*HYDROGEN+4)
     close(1000+100*HYDROGEN+5)
     close(1000+100*HYDROGEN+6)
     close(1000+100*HYDROGEN+7)
     close(1000+100*HYDROGEN+8)
     close(1000+100*HELIUM+1)
     close(1000+100*CARBON+1)
     close(1000+100*CARBON+2)
     close(1000+100*CARBON+3)
     close(1000+100*NITROGEN+1)
     close(1000+100*NITROGEN+2)
     close(1000+100*NITROGEN+3)
     close(1000+100*NITROGEN+4)
     close(1000+100*OXYGEN+1)
     close(1000+100*OXYGEN+2)
     close(1000+100*OXYGEN+3)
     close(1000+100*OXYGEN+4)
     close(1000+100*OXYGEN+5)
     close(1000+100*OXYGEN+6)
     close(1000+100*OXYGEN+7)
     close(1000+100*OXYGEN+8)
     close(1000+100*NEON+1)
     close(1000+100*NEON+2)
     close(1000+100*NEON+3)
     close(1000+100*SULFUR+1)
     close(1000+100*SULFUR+2)
     close(1000+100*SULFUR+3)
     close(1000+100*SULFUR+4)
     close(1000+100*SULFUR+5)
     close(1000+100*SULFUR+6)
     close(10000) 
   end subroutine close_output_files

   subroutine write_global_abundance_data
     implicit none
     write(*,*) 'Minimal hydrogen density in cm^-3:',minval(density)
     write(*,*) 'Mean hydrogen density in cm^-3:', sum (density)/(x_max*y_max*z_max) 
     write(*,*) 'Maximal hydrogen density in cm^-3:', maxval(density)
     write(*,*)
     write(*,*) 'Minimal electron density in cm^-3:',minval(electrons)
     write(*,*) 'Mean electron density in cm^-3:', sum(electrons)/(x_max*y_max*z_max)
     write(*,*) 'Maximal electron density in cm^-3:', maxval(electrons)
     write(*,*)
     write(*,*) 'Minimal temperature in K:', minval(temperature)
     write(*,*) 'Maximal temperature in K:', maxval(temperature)
     write(*,*) 'Temperature weigthed with ne X np:', tnenpsum/nenpsum 
     write(*,*)
     
     do z=1,z_max
        do y=1,y_max
           do x=1,x_max
              mean_of_squares=mean_of_squares+density(x, y, z)**2           
           end do
        end do
     end do
     mean_of_squares=mean_of_squares/(X_max*y_max*z_max)
     square_of_mean=(sum (density)/(x_max*y_max*z_max))**2
     write(*,*) 'Generalized clumping factor:',mean_of_squares/square_of_mean
     write(*,*) 'Standard deviation:', sqrt(mean_of_squares-square_of_mean)
     write(*,*)
     write(*,*) 'Number fractions:'
     write(*,*) '================='
     write(*,'(A17, 2ES15.5)') 'H',&
          sum(ions(:,:,:,HYDROGEN,1))/&
          (sum(ions(:,:,:,HYDROGEN,1))+sum(ions(:,:,:,HYDROGEN,2))),&
          sum(ions(:,:,:,HYDROGEN,2))/&
          (sum(ions(:,:,:,HYDROGEN,1))+sum(ions(:,:,:,HYDROGEN,2))) 
     write(*,'(A17, 3ES15.5)') 'He',&
          sum(ions(:,:,:,HELIUM,1))/&
          (sum(ions(:,:,:,HELIUM,1))+sum(ions(:,:,:,HELIUM,2))+sum(ions(:,:,:,HELIUM,3))),&
          sum(ions(:,:,:,HELIUM,2))/&
          (sum(ions(:,:,:,HELIUM,1))+sum(ions(:,:,:,HELIUM,2))+sum(ions(:,:,:,HELIUM,3))),& 
          sum(ions(:,:,:,HELIUM,3))/&
          (sum(ions(:,:,:,HELIUM,1))+sum(ions(:,:,:,HELIUM,2))+sum(ions(:,:,:,HELIUM,3))) 
     write(*,'(A17, 4ES15.5)') 'C',&
          sum(ions(:,:,:,CARBON,1))/ sum(ions(:,:,:,CARBON,:)),&
          sum(ions(:,:,:,CARBON,2))/ sum(ions(:,:,:,CARBON,:)),&
          sum(ions(:,:,:,CARBON,3))/ sum(ions(:,:,:,CARBON,:)),&
          sum(ions(:,:,:,CARBON,4))/ sum(ions(:,:,:,CARBON,:))
     write(*,'(A17, 4ES15.5)') 'N',&
          sum(ions(:,:,:,NITROGEN,1))/ sum(ions(:,:,:,NITROGEN,:)),&
          sum(ions(:,:,:,NITROGEN,2))/ sum(ions(:,:,:,NITROGEN,:)),&
          sum(ions(:,:,:,NITROGEN,3))/ sum(ions(:,:,:,NITROGEN,:)),&
          sum(ions(:,:,:,NITROGEN,4))/ sum(ions(:,:,:,NITROGEN,:))
     write(*,'(A17, 4ES15.5)') 'O',&
          sum(ions(:,:,:,OXYGEN,1))/ sum(ions(:,:,:,OXYGEN,:)),&
          sum(ions(:,:,:,OXYGEN,2))/ sum(ions(:,:,:,OXYGEN,:)),&
          sum(ions(:,:,:,OXYGEN,3))/ sum(ions(:,:,:,OXYGEN,:)),&
          sum(ions(:,:,:,OXYGEN,4))/ sum(ions(:,:,:,OXYGEN,:))
     write(*,'(A17, 4ES15.5)') 'Ne',&
          sum(ions(:,:,:,NEON,1))/ sum(ions(:,:,:,NEON,:)),&
          sum(ions(:,:,:,NEON,2))/ sum(ions(:,:,:,NEON,:)),&
          sum(ions(:,:,:,NEON,3))/ sum(ions(:,:,:,NEON,:)),&
          sum(ions(:,:,:,NEON,4))/ sum(ions(:,:,:,NEON,:))
     write(*,'(A17, 4ES15.5)') 'S',&
          sum(ions(:,:,:,SULFUR,1))/ sum(ions(:,:,:,SULFUR,:)),&
          sum(ions(:,:,:,SULFUR,2))/ sum(ions(:,:,:,SULFUR,:)),&
          sum(ions(:,:,:,SULFUR,3))/ sum(ions(:,:,:,SULFUR,:)),&
          sum(ions(:,:,:,SULFUR,4))/ sum(ions(:,:,:,SULFUR,:))
     write(*,'(A17, 4ES15.5)') 'Ar',&
          sum(ions(:,:,:,ARGON,1))/ sum(ions(:,:,:,ARGON,:)),&
          sum(ions(:,:,:,ARGON,2))/ sum(ions(:,:,:,ARGON,:)),&
          sum(ions(:,:,:,ARGON,3))/ sum(ions(:,:,:,ARGON,:)),&
          sum(ions(:,:,:,ARGON,4))/ sum(ions(:,:,:,ARGON,:))     
   end subroutine write_global_abundance_data

   subroutine calc_emissivities
     implicit none
!
     do z=1, z_max
        do y=1,y_max
           do x=1,x_max
              
              emissivity_cii_1335(x,y,z)=&
                   ions(x,y,z,CARBON,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   1335D0, 1.37D0+0.610D0,2.D0,5.1D+5)   ! two transitions with
                                                         ! very similar energies
              
              emissivity_cii_2325(x,y,z)=&
                   ions(x,y,z,CARBON,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   2324.7D0, 0.282D0,2.D0,5.0D+2)+&
                   ions(x,y,z,CARBON,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   2323.6D0, 0.419D0,2.D0,5.0D+2)+&
                   ions(x,y,z,CARBON,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   2.3226D0, 0.265D0,2.D0,5.0D+2)+&
                   emissivity_cii_1335(x,y,z)
              
              emissivity_ciii_1907(x,y,z)=&
                   ions(x,y,z,CARBON,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   1907D0, 1.05D0,1.D0,5.1D+5)
              
              emissivity_nii_optical(x, y, z)=&
                   ions(x,y,z,NITROGEN,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   6535.D0,2.64D0,9.D0,6.6D+4)
              
              emissivity_nii_auroral(x, y, z)=&
                   ions(x,y,z,NITROGEN,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   3063D0,0.29D0,9.D0,1.D+9)*&
                   (3063.D0/5755.D0)
              
              emissivity_nii_infrared122(x, y, z)=&
                   ions(x,y,z,NITROGEN,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   765000D0,0.27D0,1.D0,1.D+9)*&
                   (765000D0/1218900D0)
              
              emissivity_niii_infrared57(x, y, z)=&
                   ions(x,y,z,NITROGEN,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   573430D0,1.45D0,2.D0,1D+9)
              
              
              emissivity_oi_optical(x, y, z)=&
                   ions(x,y,z,OXYGEN,1)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   (6300.D0),&
                   0.27D0,9.D0,1D+9)
              
              emissivity_oii_optical(x, y, z)=&
                   ions(x,y,z,OXYGEN,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   3727D0,1.34D0,4.D0,6D+3)
              
              emissivity_oiii_optical(x, y, z)=&
                   ions(x,y,z,OXYGEN,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   5000D0,2.29D0,9.D0,6.8D+5)
              
              
              emissivity_oiii_auroral(x, y, z)=&
                   ions(x,y,z,OXYGEN,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   2321D0,0.29D0,9.D0,1.D+9)*&
                   (2321.D0/4363.D0)
              
              emissivity_oiii_infrared52(x, y, z)=&
                   ions(x,y,z,OXYGEN,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   326610.D0,0.27D0,1.D0,3.6D+3)*&
                   (326610.D0/518140.D0)
              
              
              emissivity_oiii_infrared88(x, y, z)=&
                   ions(x,y,z,OXYGEN,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   883560.D0,0.55D0,1.D0,5D+2)+&
                   ions(x,y,z,OXYGEN,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   326610.D0,0.27D0,1.D0,3.6D+3)*&
                   (326610.D0/883560.D0)
              
              emissivity_oiv_infrared(x, y, z)=&
                   ions(x,y,z,OXYGEN,4)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   259130D0,2.34D0,2.D0,1D+9)
              
              emissivity_neii_infrared13(x,y,z)=&
                   ions(x,y,z,NEON,2)*&
                   col_fs_cooling(&
                   temperature(x,y,z),&
                   electrons(x,y,z),&
                   128140.D0,0.28D0,4.D0,7.1D+5)
              
              emissivity_neiii_optical(x,y,z)=&
                   ions(x,y,z,NEON,3)*&
                   col_fs_cooling(&
                   temperature(x,y,z),&
                   electrons(x,y,z),&
                   3869.D0,1.36D0,9.D0,7.1D+5)
              
              emissivity_neiii_infrared16(x,y,z)=&
                   ions(x,y,z,NEON,3)*&
                   col_fs_cooling(&
                   temperature(x,y,z),&
                   electrons(x,y,z),&
                   155550.D0,0.77D0,5.D0,2.1D+5)+&
                   ions(x,y,z,NEON,3)*&
                   col_fs_cooling(&
                   temperature(x,y,z),&
                   electrons(x,y,z),&
                   108600.D0,0.21D0,5.D0,3.1D+4)*&
                   10.86D0/15.555D0
              
              emissivity_sii_4068(x, y, z)=&
                   ions(x,y,z,SULFUR,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   4068D0,3.53D0,4.D0,1D+9)*&    ! there is not only decay into 
                   ((0.078D0+0.16D0)/&           ! the ground-, but also in the
                   (0.078D0+0.16D0+0.091D0))     ! first excited state              
              
              emissivity_sii_6716(x, y, z)=&
                   ions(x,y,z,SULFUR,2)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   6716D0,6.90D0,4.D0,1D+9)
              
              emissivity_siii_optical(x, y, z)=&
                   ions(x,y,z,SULFUR,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   9000D0,6.95D0,9.D0,1D+9)
              
              emissivity_siii_infrared19(x, y, z)=&
                   ions(x,y,z,SULFUR,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   120000D0,1.31D0,1.D0,1D+9)*&
                   (120000D0/187130D0)
              
              emissivity_siii_infrared34(x, y, z)=&
                   ions(x,y,z,SULFUR,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   334700D0,3.98D0,1.D0,1D+9)+&
                   ions(x,y,z,SULFUR,3)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   120000D0,1.31D0,1.D0,1D+9)*&
                   (120000D0/ 334700D0)
              

              emissivity_siv_infrared11(x, y, z)=&
                   ions(x,y,z,SULFUR,4)*&
                   col_fs_cooling(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   105140D0,6.33D0,2.D0,1D+9)
              
              emissivity_h_beta(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -21.599D0,-0.830D0)
              
              emissivity_h_alpha(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -20.7918D0,-0.915951D0)
              
              emissivity_h_gamma(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&                
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -22.064D0,-0.797D0)
              
              emissivity_h_delta(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -22.352D0,-0.790D0)

              emissivity_pasch_alpha(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&
                   temperature(x, y, z),&           
                   electrons(x, y, z),&
                   -20.919D0,-1.113D0)
              
              emissivity_pasch_beta(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -21.680D0,-1.005D0)

              emissivity_brack_alpha(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -21.114D0,-1.220D0)

              emissivity_brack_beta(x, y, z)=&
                   ions(x, y, z,HYDROGEN,2)*&
                   recombination_line(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -21.730D0,-1.129D0)
              
              emissivity_hei_optical(x,y,z)=&
                   ions(x,y,z,HELIUM,2)*&
                   recombination_line(&
                   temperature(x, y, z),&
                   electrons(x, y, z),&
                   -20.5969D0,-1.04757D0)
           end do
        end do
     end do
   end subroutine calc_emissivities

   subroutine write_global_emission_data
     implicit none
     write(*,'(1A17,2A28)')&
          'Line(s)','Total Emission in erg s^-1','Emission relative to H_beta'
     write(*,'(1A17,2A28)') '=======','==========================',&
                            '==========================='
     write(*,*)
     write(*,'(A17,2ES28.5)') 'H_alpha',&
          sum(emissivity_h_alpha)*(l_cell*nc_parsec)**3,&
          sum(emissivity_h_alpha)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'H_beta',&
          sum(emissivity_h_beta)*(l_cell*nc_parsec)**3,&
          sum(emissivity_h_beta)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'H_gamma',&
          sum(emissivity_h_gamma)*(l_cell*nc_parsec)**3,&
          sum(emissivity_h_gamma)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'H_delta',&
          sum(emissivity_h_delta)*(l_cell*nc_parsec)**3,&
          sum(emissivity_h_delta)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'pasch_alpha',&
          sum(emissivity_pasch_alpha)*(l_cell*nc_parsec)**3,&
          sum(emissivity_pasch_alpha)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'pasch_beta',&
          sum(emissivity_pasch_beta)*(l_cell*nc_parsec)**3,&
          sum(emissivity_pasch_beta)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'brack_alpha',&
          sum(emissivity_brack_alpha)*(l_cell*nc_parsec)**3,&
          sum(emissivity_brack_alpha)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'brack_beta',&
          sum(emissivity_brack_beta)*(l_cell*nc_parsec)**3,&
          sum(emissivity_brack_beta)/sum(emissivity_h_beta)
!
     write(*,*)
     write(*,'(A17,2ES28.5)') 'HeI 5876',&
          sum(emissivity_hei_optical)*(l_cell*nc_parsec)**3,&
          sum(emissivity_hei_optical)/sum(emissivity_h_beta)
!
     write(*,*)
     write(*,'(A17,2ES28.5)') 'CII 1335',&
          sum(emissivity_cii_1335)*(l_cell*nc_parsec)**3,&
          sum(emissivity_cii_1335)/sum(emissivity_h_beta)      
     write(*,'(A17,2ES28.5)') 'CII] 2325+',&
          sum(emissivity_cii_2325)*(l_cell*nc_parsec)**3,&
          sum(emissivity_cii_2325)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') 'CIII] 1907+1909',&
          sum(emissivity_ciii_1907)*(l_cell*nc_parsec)**3,&
          sum(emissivity_ciii_1907)/sum(emissivity_h_beta)      
!
     write(*,*)
     write(*,'(A17,2ES28.5)') '[NII] 6584+6548',&
          sum(emissivity_nii_optical)*(l_cell*nc_parsec)**3,&
          sum(emissivity_nii_optical)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[NII] 5755',&
          sum(emissivity_nii_auroral)*(l_cell*nc_parsec)**3,&
          sum(emissivity_nii_auroral)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[NII] 122 mu',&
          sum(emissivity_nii_infrared122)*(l_cell*nc_parsec)**3,&
          sum(emissivity_nii_infrared122)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[NIII] 57.3 mu',&
          sum(emissivity_niii_infrared57)*(l_cell*nc_parsec)**3,&
          sum(emissivity_niii_infrared57)/sum(emissivity_h_beta) 
!
     write(*,*)
     write(*,'(A17,2ES28.5)') '[0I] 6300+6363',&
          sum(emissivity_oi_optical)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oi_optical)/sum(emissivity_h_beta) 
     write(*,'(A17,2ES28.5)') '[0II] 3726+3729',&
          sum(emissivity_oii_optical)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oii_optical)/sum(emissivity_h_beta) 
     write(*,'(A17,2ES28.5)') '[0II] 7320+7330',&
          sum(emissivity_oii_7320)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oii_7320)/sum(emissivity_h_beta) 
     write(*,'(A17,2ES28.5)') '[0III] 5007+4959',&
          sum(emissivity_oiii_optical)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oiii_optical)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[0III] 4363',&
          sum(emissivity_oiii_auroral)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oiii_auroral)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[0III] 51.8 mu',&
          sum(emissivity_oiii_infrared52)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oiii_infrared52)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[0III] 88.3 mu',&
          sum(emissivity_oiii_infrared88)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oiii_infrared88)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[0IV] 25.9 mu',&
          sum(emissivity_oiv_infrared)*(l_cell*nc_parsec)**3,&
          sum(emissivity_oiv_infrared)/sum(emissivity_h_beta)
!
     write(*,*)
     write(*,'(A17,2ES28.5)') '[NeII] 12.8 mu',&
          sum(emissivity_neii_infrared13)*(l_cell*nc_parsec)**3,&
          sum(emissivity_neii_infrared13)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[NeIII] 3869+3968',&
          sum(emissivity_neiii_optical)*(l_cell*nc_parsec)**3,&
          sum(emissivity_neiii_optical)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[NeIII] 15.5 mu',&
          sum(emissivity_neiii_infrared16)*(l_cell*nc_parsec)**3,&
          sum(emissivity_neiii_infrared16)/sum(emissivity_h_beta)
!
     write(*,*)
     write(*,'(A17,2ES28.5)') '[SII] 4068+4076',&
          sum(emissivity_sii_4068)*(l_cell*nc_parsec)**3,&
          sum(emissivity_sii_4068)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[SII] 6716+6731',&
          sum(emissivity_sii_6716)*(l_cell*nc_parsec)**3,&
          sum(emissivity_sii_6716)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[SIII] 9531+5069',&
          sum(emissivity_siii_optical)*(l_cell*nc_parsec)**3,&
          sum(emissivity_siii_optical)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[SIII] 18.7 mu',&
          sum(emissivity_siii_infrared19)*(l_cell*nc_parsec)**3,&
          sum(emissivity_siii_infrared19)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[SIII] 33.6 mu',&
          sum(emissivity_siii_infrared34)*(l_cell*nc_parsec)**3,&
          sum(emissivity_siii_infrared34)/sum(emissivity_h_beta)
     write(*,'(A17,2ES28.5)') '[SIV] 10.5 um',&
          sum(emissivity_siv_infrared11)*(l_cell*nc_parsec)**3,&
          sum(emissivity_siv_infrared11)/sum(emissivity_h_beta)
!     
     write(*,*)
     write(*,*)
     write(*,'(1A17,A42)') 'Line(s)','max. specific  Emission in erg s^-1 cm^-3'
     write(*,'(1A17,A42)') '=======','========================================='
     write(*,*)
     write(*,'(A17,ES42.5)') 'H_alpha', maxval(emissivity_h_alpha)
     write(*,'(A17,ES42.5)') 'H_beta', maxval(emissivity_h_beta)
     write(*,'(A17,ES42.5)') 'H_gamma', maxval(emissivity_h_gamma)
     write(*,'(A17,ES42.5)') 'H_delta', maxval(emissivity_h_delta)
     write(*,'(A17,ES42.5)') 'pasch_alpha', maxval(emissivity_pasch_alpha)
     write(*,'(A17,ES42.5)') 'pasch_beta', maxval(emissivity_pasch_beta)
     write(*,'(A17,ES42.5)') 'brack_alpha', maxval(emissivity_brack_alpha)
     write(*,'(A17,ES42.5)') 'brack_beta', maxval(emissivity_brack_beta)
     write(*,*)
     write(*,'(A17,ES42.5)') 'HeI 5876', maxval(emissivity_hei_optical)
     write(*,*)
     write(*,'(A17,ES42.5)') 'CII 1335', maxval(emissivity_cii_1335)
     write(*,'(A17,ES42.5)') 'CII] 2325+', maxval(emissivity_cii_2325)
     write(*,'(A17,ES42.5)') 'CIII] 1907+1909', maxval(emissivity_ciii_1907)   
     write(*,*)
     write(*,'(A17,ES42.5)') '[NII] 6584+6548', maxval(emissivity_nii_optical)
     write(*,'(A17,ES42.5)') '[NII] 5755', maxval(emissivity_nii_auroral)
     write(*,'(A17,ES42.5)') '[NII] 122 mu', maxval(emissivity_nii_infrared122)
     write(*,'(A17,ES42.5)') '[NIII] 57 mu', maxval(emissivity_niii_infrared57)
     write(*,*)
     write(*,'(A17,ES42.5)') '[0I] 6300+6363', maxval(emissivity_oi_optical)
     write(*,'(A17,ES42.5)') '[0II] 3726+3729', maxval(emissivity_oii_optical)
     write(*,'(A17,ES42.5)') '[0III] 5007+4959', maxval(emissivity_oiii_optical)
     write(*,'(A17,ES42.5)') '[0III] 4363', maxval(emissivity_oiii_auroral)
     write(*,'(A17,ES42.5)') '[0III] 51.8 mu', maxval(emissivity_oiii_infrared52)
     write(*,'(A17,ES42.5)') '[0III] 88.3 mu', maxval(emissivity_oiii_infrared88)
     write(*,'(A17,ES42.5)') '[OIV] 25.9 um', maxval(emissivity_oiv_infrared)
     write(*,*)
     write(*,'(A17,ES42.5)') '[NeII] 12.8 mu', maxval(emissivity_neii_infrared13)
     write(*,'(A17,ES42.5)') '[NeIII] 3869+3968',&
                              maxval(emissivity_neiii_optical)
     write(*,'(A17,ES42.5)') '[NeIII] 15.5 mu',&
                             maxval(emissivity_neiii_infrared16)     
     write(*,*)
     write(*,'(A17,ES42.5)') '[SII] 4068+4076', maxval(emissivity_sii_4068)
     write(*,'(A17,ES42.5)') '[SII] 6716+6731', maxval(emissivity_sii_6716)
     write(*,'(A17,ES42.5)') '[SIII] 9531+9069', maxval(emissivity_siii_optical)
     write(*,'(A17,ES42.5)') '[SIII] 18.7 um', maxval(emissivity_siii_infrared19)
     write(*,'(A17,ES42.5)') '[SIII] 33.6 um', maxval(emissivity_siii_infrared34)
     write(*,'(A17,ES42.5)') '[SIV] 10.5 um', maxval(emissivity_siv_infrared11)   
   end subroutine write_global_emission_data

   subroutine write_maps_x
     implicit none
     do z=1, z_max
        do y=1, y_max
           xHbeta=sum(emissivity_h_beta(:,y, z))*l_cell*nc_parsec
           em_x=sum (electrons(:, y, z)**2)*l_cell
           write(10000,*) em_x
          
           write(1000+100*HYDROGEN+1,*)&
                sum(emissivity_h_alpha(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_h_alpha(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(6563*1E-8))/nc_sphere,&
                sum(emissivity_h_alpha(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_h_alpha(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_alpha_x)&
                max_emissivity_h_alpha_x=&
                sum(emissivity_h_alpha(:,y, z))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+2,*)&
                sum(emissivity_h_beta(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_h_beta(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4861*1E-8))/nc_sphere,&
                sum(emissivity_h_beta(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
                           
           if (sum(emissivity_h_beta(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_beta_x)&
                max_emissivity_h_beta_x=&
                sum(emissivity_h_beta(:,y, z))*l_cell*nc_parsec
       
           write(1000+100*HYDROGEN+3,*)&
                sum(emissivity_h_gamma(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_h_gamma(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4340*1E-8))/nc_sphere,&
                sum(emissivity_h_gamma(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_h_gamma(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_gamma_x)&
                max_emissivity_h_gamma_x=&
                sum(emissivity_h_gamma(:,y, z))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+4,*)&
                sum(emissivity_h_delta(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_h_delta(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4102*1E-8))/nc_sphere,&
                sum(emissivity_h_delta(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_h_delta(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_delta_x)&
                max_emissivity_h_delta_x=&
                sum(emissivity_h_delta(:,y, z))*l_cell*nc_parsec 
                      
           write(1000+100*HYDROGEN+5,*)&
                sum(emissivity_pasch_alpha(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_pasch_alpha(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(18745*1E-8))/nc_sphere,&
                sum(emissivity_pasch_alpha(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_pasch_alpha(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_pasch_alpha_x)&
                max_emissivity_pasch_alpha_x=&
                sum(emissivity_pasch_alpha(:,y, z))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+6,*)&
                sum(emissivity_pasch_beta(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_pasch_beta(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(12814*1E-8))/nc_sphere,&
                sum(emissivity_pasch_beta(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_pasch_beta(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_pasch_beta_x)&
                max_emissivity_pasch_beta_x=&
                sum(emissivity_pasch_beta(:,y, z))*l_cell*nc_parsec 
           
           
           write(1000+100*HYDROGEN+7,*)&
                sum(emissivity_brack_alpha(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_brack_alpha(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(40525*1E-8))/nc_sphere,&
                sum(emissivity_brack_alpha(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
                
           
           if (sum(emissivity_brack_alpha(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_brack_alpha_x)&
                max_emissivity_brack_alpha_x=&
                sum(emissivity_brack_alpha(:,y, z))*l_cell*nc_parsec 
                      
           write(1000+100*HYDROGEN+8,*)&
                sum(emissivity_brack_beta(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_brack_beta(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(26259*1E-8))/nc_sphere,&
                sum(emissivity_brack_beta(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_brack_beta(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_brack_beta_x)&
                max_emissivity_brack_beta_x=&
                sum(emissivity_brack_beta(:,y, z))*l_cell*nc_parsec 
!
!
           write(1000+100*HELIUM+1,*)&
                sum(emissivity_hei_optical(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_hei_optical(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(5876*1E-8))/nc_sphere,&
                sum(emissivity_hei_optical(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_hei_optical(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_hei_optical_x)&
                max_emissivity_hei_optical_x=&
                sum(emissivity_hei_optical(:,y, z))*l_cell*nc_parsec 
!
!
           write(1000+100*CARBON+1,*)&
                sum(emissivity_cii_1335(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_cii_1335(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(2325*1E-8))/nc_sphere,&
                sum(emissivity_cii_1335(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
                
           
           if (sum(emissivity_cii_1335(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_cii_1335_x)&
                max_emissivity_cii_1335_x=&
                sum(emissivity_cii_1335(:,y, z))*l_cell*nc_parsec 

           write(1000+100*CARBON+2,*)&
                sum(emissivity_cii_2325(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_cii_2325(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(1335*1E-8))/nc_sphere,&
                sum(emissivity_cii_2325(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_cii_2325(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_cii_2325_x)&
                max_emissivity_cii_2325_x=&
                sum(emissivity_cii_2325(:,y, z))*l_cell*nc_parsec 

           write(1000+100*CARBON+3,*)&
                sum(emissivity_ciii_1907(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_ciii_1907(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(1908*1E-8))/nc_sphere,&
                sum(emissivity_ciii_1907(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_ciii_1907(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_ciii_1907_x)&
                max_emissivity_ciii_1907_x=&
                sum(emissivity_ciii_1907(:,y, z))*l_cell*nc_parsec 
!
!          
           write(1000+100*NITROGEN+1,*)&
                sum(emissivity_nii_optical(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_nii_optical(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(6548+6584)*1E-8))/nc_sphere,&
                sum(emissivity_nii_optical(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_nii_optical(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_optical_x)&
                max_emissivity_nii_optical_x=&
                sum(emissivity_nii_optical(:,y, z))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+2,*)&
                sum(emissivity_nii_auroral(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_nii_auroral(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(5755*1E-8))/nc_sphere,&
                sum(emissivity_nii_auroral(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_nii_auroral(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_auroral_x)&
                max_emissivity_nii_auroral_x=&
                sum(emissivity_nii_auroral(:,y, z))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+3,*)&
                sum(emissivity_nii_infrared122(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_nii_infrared122(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(121.89*1E-4))/nc_sphere,&
                sum(emissivity_nii_infrared122(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_nii_infrared122(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_infrared122_x)&
                max_emissivity_nii_infrared122_x=&
                sum(emissivity_nii_infrared122(:,y, z))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+4,*)&
                sum(emissivity_niii_infrared57(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_niii_infrared57(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(57.343*1E-4))/nc_sphere,&
                sum(emissivity_niii_infrared57(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_niii_infrared57(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_niii_infrared57_x)&
                max_emissivity_niii_infrared57_x=&
                sum(emissivity_niii_infrared57(:,y, z))*l_cell*nc_parsec 
!           
!           
           write(1000+100*OXYGEN+1,*)&
                sum(emissivity_oi_optical(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oi_optical(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( (1./2.)*(6300+6363)*1E-8))/nc_sphere,&
                sum(emissivity_oi_optical(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oi_optical(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oi_optical_x)&
                max_emissivity_oi_optical_x=&
                sum(emissivity_oi_optical(:,y, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+2,*)&
                sum(emissivity_oii_optical(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oii_optical(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 3727*1E-8))/nc_sphere,&
                sum(emissivity_oii_optical(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oii_optical(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oii_optical_x)&
                max_emissivity_oii_optical_x=&
                sum(emissivity_oii_optical(:,y, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+3,*)&
                sum(emissivity_oii_7320(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oii_7320(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 3727*1E-8))/nc_sphere,&
                sum(emissivity_oii_7320(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oii_7320(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oii_7320_x)&
                max_emissivity_oii_7320_x=&
                sum(emissivity_oii_7320(:,y, z))*l_cell*nc_parsec 
                      
           write(1000+100*OXYGEN+4,*)&
                sum(emissivity_oiii_optical(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_optical(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(5007+4959)*1E-8))/nc_sphere,&
                sum(emissivity_oiii_optical(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oiii_optical(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_optical_x)&
                max_emissivity_oiii_optical_x=&
                sum(emissivity_oiii_optical(:,y, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+5,*)&
                sum(emissivity_oiii_auroral(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_auroral(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4363*1E-8))/nc_sphere,&
                sum(emissivity_oiii_auroral(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oiii_auroral(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_auroral_x)&
                max_emissivity_oiii_auroral_x=&
                sum(emissivity_oiii_auroral(:,y, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+6,*)&
                sum(emissivity_oiii_infrared52(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_infrared52(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(51.814*1E-4))/nc_sphere,&
                sum(emissivity_oiii_infrared52(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oiii_infrared52(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_infrared52_x)&
                max_emissivity_oiii_infrared52_x=&
                sum(emissivity_oiii_infrared52(:,y, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+7,*)&
                sum(emissivity_oiii_infrared88(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_infrared88(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(88.356*1E-4))/nc_sphere,&
                sum(emissivity_oiii_infrared88(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oiii_infrared88(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_infrared88_x)&
                max_emissivity_oiii_infrared88_x=&
                sum(emissivity_oiii_infrared88(:,y, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+8,*)&
                sum(emissivity_oiv_infrared(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_oiv_infrared(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(25.916*1E-4))/nc_sphere,&
                sum(emissivity_oiv_infrared(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_oiv_infrared(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiv_infrared_x)&
                max_emissivity_oiv_infrared_x=&
                sum(emissivity_oiv_infrared(:,y, z))*l_cell*nc_parsec 
!
!
           write(1000+100*NEON+1,*)&
                sum(emissivity_neii_infrared13(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_neii_infrared13(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 12.814*1E-4))/nc_sphere,&
                sum(emissivity_neii_infrared13(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_neii_infrared13(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_neii_infrared13_x)&
                max_emissivity_neii_infrared13_x=&
                sum(emissivity_neii_infrared13(:,y, z))*l_cell*nc_parsec 
                      
           write(1000+100*NEON+2,*)&
                sum(emissivity_neiii_optical(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_neiii_optical(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(3869+3968)*1E-8))/nc_sphere,&
                sum(emissivity_neiii_optical(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_neiii_optical(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_neiii_optical_x)&
                max_emissivity_neiii_optical_x=&
                sum(emissivity_neiii_optical(:,y, z))*l_cell*nc_parsec 

           write(1000+100*NEON+3,*)&
                sum(emissivity_neiii_infrared16(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_neiii_infrared16(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(15.555*1E-4))/nc_sphere,&
                sum(emissivity_neiii_infrared16(:,y, z))*l_cell*nc_parsec/&
                xHbeta
           
           if (sum(emissivity_neiii_infrared16(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_neiii_infrared16_x)&
                max_emissivity_neiii_infrared16_x=&
                sum(emissivity_neiii_infrared16(:,y, z))*l_cell*nc_parsec 
!
!
           write(1000+100*SULFUR+1,*)&
                sum(emissivity_sii_4068(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_sii_4068(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(4068+4076)*1E-8))/nc_sphere,&
                sum(emissivity_sii_4068(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_sii_4068(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_sii_4068_x)&
                max_emissivity_sii_4068_x=&
                sum(emissivity_sii_4068(:,y, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+2,*)&
                sum(emissivity_sii_6716(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_sii_6716(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(6716+6731)*1E-8))/nc_sphere,&
                sum(emissivity_sii_6716(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_sii_6716(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_sii_6716_x)&
                max_emissivity_sii_6716_x=&
                sum(emissivity_sii_6716(:,y, z))*l_cell*nc_parsec 
                      
           write(1000+100*SULFUR+3,*)&
                sum(emissivity_siii_optical(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_siii_optical(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(9532+9069)*1E-8))/nc_sphere,&
                sum(emissivity_siii_optical(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_siii_optical(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_optical_x)&
                max_emissivity_siii_optical_x=&
                sum(emissivity_siii_optical(:,y, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+4,*)&
                sum(emissivity_siii_infrared19(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_siii_infrared19(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(18.713*1E-4))/nc_sphere,&
                sum(emissivity_siii_infrared19(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_siii_infrared19(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_infrared19_x)&
                max_emissivity_siii_infrared19_x=&
                sum(emissivity_siii_infrared19(:,y, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+5,*)&
                sum(emissivity_siii_infrared34(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_siii_infrared34(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(33.47*1E-4))/nc_sphere,&
                sum(emissivity_siii_infrared34(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_siii_infrared34(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_infrared34_x)&
                max_emissivity_siii_infrared34_x=&
                sum(emissivity_siii_infrared34(:,y, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+6,*)&
                sum(emissivity_siv_infrared11(:,y, z))*l_cell*nc_parsec,&
                sum(emissivity_siv_infrared11(:,y, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(10.514*1E-4))/nc_sphere,&
                sum(emissivity_siv_infrared11(:,y, z))*l_cell*nc_parsec/&
                xHbeta*maskedvalue(1.0_dp,em_x,em_min)
           
           if (sum(emissivity_siv_infrared11(:,y, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siv_infrared11_x)&
                max_emissivity_siv_infrared11_x=&
                sum(emissivity_siv_infrared11(:,y, z))*l_cell*nc_parsec 
           
        end do
     end do
   end subroutine write_maps_x


   subroutine write_maps_y
     implicit none
     do z=1, z_max
        do x=1, x_max
!           
           yHbeta=sum(emissivity_h_beta(x,:, z))*l_cell*nc_parsec
           em_y=sum (electrons(x, :, z)**2)*l_cell
           write(10000,*) em_y
!
           write(1000+100*HYDROGEN+1,*)&
                sum(emissivity_h_alpha(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_h_alpha(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(6563*1E-8))/nc_sphere,&
                sum(emissivity_h_beta(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_h_alpha(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_alpha_y)&
                max_emissivity_h_alpha_y=&
                sum(emissivity_h_alpha(x,:, z))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+2,*)&
                sum(emissivity_h_beta(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_h_beta(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4861*1E-8))/nc_sphere,&
                sum(emissivity_h_beta(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_h_beta(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_beta_y)&
                max_emissivity_h_beta_y=&
                sum(emissivity_h_beta(x,:, z))*l_cell*nc_parsec
       
           write(1000+100*HYDROGEN+3,*)&
                sum(emissivity_h_gamma(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_h_gamma(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4340*1E-8))/nc_sphere,&
                sum(emissivity_h_gamma(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_h_gamma(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_gamma_y)&
                max_emissivity_h_gamma_y=&
                sum(emissivity_h_gamma(x,:, z))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+4,*)&
                sum(emissivity_h_delta(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_h_delta(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4102*1E-8))/nc_sphere,&
                sum(emissivity_h_delta(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_h_delta(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_h_delta_y)&
                max_emissivity_h_delta_y=&
                sum(emissivity_h_delta(x,:, z))*l_cell*nc_parsec 
                      
           write(1000+100*HYDROGEN+5,*)&
                sum(emissivity_pasch_alpha(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_pasch_alpha(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(18745*1E-8))/nc_sphere,&
                sum(emissivity_pasch_alpha(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_pasch_alpha(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_pasch_alpha_y)&
                max_emissivity_pasch_alpha_y=&
                sum(emissivity_pasch_alpha(x,:, z))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+6,*)&
                sum(emissivity_pasch_beta(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_pasch_beta(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(12814*1E-8))/nc_sphere,&
                sum(emissivity_pasch_beta(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_pasch_beta(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_pasch_beta_y)&
                max_emissivity_pasch_beta_y=&
                sum(emissivity_pasch_beta(x,:, z))*l_cell*nc_parsec 
           
           
           write(1000+100*HYDROGEN+7,*)&
                sum(emissivity_brack_alpha(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_brack_alpha(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(40525*1E-8))/nc_sphere,&
                sum(emissivity_brack_alpha(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_brack_alpha(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_brack_alpha_y)&
                max_emissivity_brack_alpha_y=&
                sum(emissivity_brack_alpha(x,:, z))*l_cell*nc_parsec 
                      
           write(1000+100*HYDROGEN+8,*)&
                sum(emissivity_brack_beta(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_brack_beta(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(26259*1E-8))/nc_sphere,&
                sum(emissivity_brack_beta(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_brack_beta(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_brack_beta_y)&
                max_emissivity_brack_beta_y=&
                sum(emissivity_brack_beta(x,:, z))*l_cell*nc_parsec 
!
!
           write(1000+100*HELIUM+1,*)&
                sum(emissivity_hei_optical(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_hei_optical(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(5876*1E-8))/nc_sphere,&
                sum(emissivity_hei_optical(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_hei_optical(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_hei_optical_y)&
                max_emissivity_hei_optical_y=&
                sum(emissivity_hei_optical(x,:, z))*l_cell*nc_parsec 
!
!
           write(1000+100*CARBON+1,*)&
                sum(emissivity_cii_1335(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_cii_1335(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(2325*1E-8))/nc_sphere,&
                sum(emissivity_cii_1335(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
 
           if (sum(emissivity_cii_1335(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_cii_1335_y)&
                max_emissivity_cii_1335_y=&
                sum(emissivity_cii_1335(x,:, z))*l_cell*nc_parsec 

           write(1000+100*CARBON+2,*)&
                sum(emissivity_cii_2325(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_cii_2325(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(1335*1E-8))/nc_sphere,&
                sum(emissivity_cii_2325(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_cii_2325(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_cii_2325_y)&
                max_emissivity_cii_2325_y=&
                sum(emissivity_cii_2325(x,:, z))*l_cell*nc_parsec 

           write(1000+100*CARBON+3,*)&
                sum(emissivity_ciii_1907(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_ciii_1907(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(1908*1E-8))/nc_sphere,&
                sum(emissivity_ciii_1907(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_ciii_1907(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_ciii_1907_y)&
                max_emissivity_ciii_1907_y=&
                sum(emissivity_ciii_1907(x,:, z))*l_cell*nc_parsec 
!
!          
           write(1000+100*NITROGEN+1,*)&
                sum(emissivity_nii_optical(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_nii_optical(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(6548+6584)*1E-8))/nc_sphere,&
                sum(emissivity_nii_optical(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_nii_optical(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_optical_y)&
                max_emissivity_nii_optical_y=&
                sum(emissivity_nii_optical(x,:, z))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+2,*)&
                sum(emissivity_nii_auroral(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_nii_auroral(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(5755*1E-8))/nc_sphere,&
                sum(emissivity_nii_auroral(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_nii_auroral(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_auroral_y)&
                max_emissivity_nii_auroral_y=&
                sum(emissivity_nii_auroral(x,:, z))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+3,*)&
                sum(emissivity_nii_infrared122(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_nii_infrared122(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(121.89*1E-4))/nc_sphere,&
                sum(emissivity_nii_infrared122(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_nii_infrared122(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_infrared122_y)&
                max_emissivity_nii_infrared122_y=&
                sum(emissivity_nii_infrared122(x,:, z))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+4,*)&
                sum(emissivity_niii_infrared57(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_niii_infrared57(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(57.343*1E-4))/nc_sphere,&
                sum(emissivity_niii_infrared57(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_niii_infrared57(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_niii_infrared57_y)&
                max_emissivity_niii_infrared57_y=&
                sum(emissivity_niii_infrared57(x,:, z))*l_cell*nc_parsec 
!           
!           
           write(1000+100*OXYGEN+1,*)&
                sum(emissivity_oi_optical(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oi_optical(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( (1./2.)*(6300+6363)*1E-8))/nc_sphere,&
                sum(emissivity_oi_optical(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oi_optical(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oi_optical_y)&
                max_emissivity_oi_optical_y=&
                sum(emissivity_oi_optical(x,:, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+2,*)&
                sum(emissivity_oii_optical(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oii_optical(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 3727*1E-8))/nc_sphere,&
                sum(emissivity_oii_optical(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oii_optical(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oii_optical_y)&
                max_emissivity_oii_optical_y=&
                sum(emissivity_oii_optical(x,:, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+3,*)&
                sum(emissivity_oii_7320(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oii_7320(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 3727*1E-8))/nc_sphere,&
                sum(emissivity_oii_7320(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oii_7320(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oii_7320_y)&
                max_emissivity_oii_7320_y=&
                sum(emissivity_oii_7320(x,:, z))*l_cell*nc_parsec 
                      
           write(1000+100*OXYGEN+4,*)&
                sum(emissivity_oiii_optical(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_optical(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(5007+4959)*1E-8))/nc_sphere,&
                sum(emissivity_oiii_optical(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oiii_optical(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_optical_y)&
                max_emissivity_oiii_optical_y=&
                sum(emissivity_oiii_optical(x,:, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+5,*)&
                sum(emissivity_oiii_auroral(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_auroral(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4363*1E-8))/nc_sphere,&
                sum(emissivity_oiii_auroral(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oiii_auroral(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_auroral_y)&
                max_emissivity_oiii_auroral_y=&
                sum(emissivity_oiii_auroral(x,:, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+6,*)&
                sum(emissivity_oiii_infrared52(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_infrared52(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(51.814*1E-4))/nc_sphere,&
                sum(emissivity_oiii_infrared52(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oiii_infrared52(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_infrared52_y)&
                max_emissivity_oiii_infrared52_y=&
                sum(emissivity_oiii_infrared52(x,:, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+7,*)&
                sum(emissivity_oiii_infrared88(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oiii_infrared88(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(88.356*1E-4))/nc_sphere,&
                sum(emissivity_oiii_infrared88(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oiii_infrared88(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_infrared88_y)&
                max_emissivity_oiii_infrared88_y=&
                sum(emissivity_oiii_infrared88(x,:, z))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+8,*)&
                sum(emissivity_oiv_infrared(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_oiv_infrared(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(25.916*1E-4))/nc_sphere,&
                sum(emissivity_oiv_infrared(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_oiv_infrared(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_oiv_infrared_y)&
                max_emissivity_oiv_infrared_y=&
                sum(emissivity_oiv_infrared(x,:, z))*l_cell*nc_parsec 
!
!
           write(1000+100*NEON+1,*)&
                sum(emissivity_neii_infrared13(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_neii_infrared13(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 12.814*1E-4))/nc_sphere,&
                sum(emissivity_neii_infrared13(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_neii_infrared13(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_neii_infrared13_y)&
                max_emissivity_neii_infrared13_y=&
                sum(emissivity_neii_infrared13(x,:, z))*l_cell*nc_parsec 
                      
           write(1000+100*NEON+2,*)&
                sum(emissivity_neiii_optical(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_neiii_optical(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(3869+3968)*1E-8))/nc_sphere,&
                sum(emissivity_neiii_optical(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_neiii_optical(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_neiii_optical_y)&
                max_emissivity_neiii_optical_y=&
                sum(emissivity_neiii_optical(x,:, z))*l_cell*nc_parsec 

           write(1000+100*NEON+3,*)&
                sum(emissivity_neiii_infrared16(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_neiii_infrared16(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(15.555*1E-4))/nc_sphere,&
                sum(emissivity_neiii_infrared16(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_neiii_infrared16(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_neiii_infrared16_y)&
                max_emissivity_neiii_infrared16_y=&
                sum(emissivity_neiii_infrared16(x,:, z))*l_cell*nc_parsec 
!
!
           write(1000+100*SULFUR+1,*)&
                sum(emissivity_sii_4068(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_sii_4068(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(4068+4076)*1E-8))/nc_sphere,&
                sum(emissivity_sii_4068(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_sii_4068(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_sii_4068_y)&
                max_emissivity_sii_4068_y=&
                sum(emissivity_sii_4068(x,:, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+2,*)&
                sum(emissivity_sii_6716(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_sii_6716(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(6716+6731)*1E-8))/nc_sphere,&
                sum(emissivity_sii_6716(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_sii_6716(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_sii_6716_y)&
                max_emissivity_sii_6716_y=&
                sum(emissivity_sii_6716(x,:, z))*l_cell*nc_parsec 
                      
           write(1000+100*SULFUR+3,*)&
                sum(emissivity_siii_optical(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_siii_optical(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(9532+9069)*1E-8))/nc_sphere,&
                sum(emissivity_siii_optical(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_siii_optical(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_optical_y)&
                max_emissivity_siii_optical_y=&
                sum(emissivity_siii_optical(x,:, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+4,*)&
                sum(emissivity_siii_infrared19(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_siii_infrared19(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(18.713*1E-4))/nc_sphere,&
                sum(emissivity_siii_infrared19(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_siii_infrared19(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_infrared19_y)&
                max_emissivity_siii_infrared19_y=&
                sum(emissivity_siii_infrared19(x,:, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+5,*)&
                sum(emissivity_siii_infrared34(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_siii_infrared34(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(33.47*1E-4))/nc_sphere,&
                sum(emissivity_siii_infrared34(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_siii_infrared34(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_infrared34_y)&
                max_emissivity_siii_infrared34_y=&
                sum(emissivity_siii_infrared34(x,:, z))*l_cell*nc_parsec 

           write(1000+100*SULFUR+6,*)&
                sum(emissivity_siv_infrared11(x,:, z))*l_cell*nc_parsec,&
                sum(emissivity_siv_infrared11(x,:, z))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(10.514*1E-4))/nc_sphere,&
                sum(emissivity_siv_infrared11(x,:, z))*l_cell*nc_parsec/&
                yHbeta*maskedvalue(1.0_dp,em_y,em_min)
           
           if (sum(emissivity_siv_infrared11(x,:, z))*l_cell*nc_parsec .gt.&
                max_emissivity_siv_infrared11_y)&
                max_emissivity_siv_infrared11_y=&
                sum(emissivity_siv_infrared11(x,:, z))*l_cell*nc_parsec 
!
!           
       end do
     end do
   end subroutine write_maps_y

   subroutine write_maps_z
     implicit none
     do y=1, y_max
        do x=1, x_max
!
           zHbeta=sum(emissivity_h_beta(x,y,:))*l_cell*nc_parsec
           em_z=sum (electrons(x, y, :)**2)*l_cell
           write(10000,*) em_z
           
!
           write(1000+100*HYDROGEN+1,*)&
                sum(emissivity_h_alpha(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_h_alpha(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(6563*1E-8))/nc_sphere,&
                sum(emissivity_h_alpha(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
                           
           if (sum(emissivity_h_alpha(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_h_alpha_z)&
                max_emissivity_h_alpha_z=&
                sum(emissivity_h_alpha(x,y,:))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+2,*)&
                sum(emissivity_h_beta(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_h_beta(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4861*1E-8))/nc_sphere,&
                sum(emissivity_h_beta(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
                           
           if (sum(emissivity_h_beta(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_h_beta_z)&
                max_emissivity_h_beta_z=&
                sum(emissivity_h_beta(x,y,:))*l_cell*nc_parsec
       
           write(1000+100*HYDROGEN+3,*)&
                sum(emissivity_h_gamma(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_h_gamma(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4340*1E-8))/nc_sphere,&
                sum(emissivity_h_gamma(x,y,:))*l_cell*nc_parsec/&
                zHbeta
           
           if (sum(emissivity_h_gamma(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_h_gamma_z)&
                max_emissivity_h_gamma_z=&
                sum(emissivity_h_gamma(x,y,:))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+4,*)&
                sum(emissivity_h_delta(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_h_delta(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4102*1E-8))/nc_sphere,&
                sum(emissivity_h_delta(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_h_delta(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_h_delta_z)&
                max_emissivity_h_delta_z=&
                sum(emissivity_h_delta(x,y,:))*l_cell*nc_parsec 
                      
           write(1000+100*HYDROGEN+5,*)&
                sum(emissivity_pasch_alpha(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_pasch_alpha(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(18745*1E-8))/nc_sphere,&
                sum(emissivity_pasch_alpha(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_pasch_alpha(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_pasch_alpha_z)&
                max_emissivity_pasch_alpha_z=&
                sum(emissivity_pasch_alpha(x,y,:))*l_cell*nc_parsec 
           
           write(1000+100*HYDROGEN+6,*)&
                sum(emissivity_pasch_beta(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_pasch_beta(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(12814*1E-8))/nc_sphere,&
                sum(emissivity_pasch_beta(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_pasch_beta(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_pasch_beta_z)&
                max_emissivity_pasch_beta_z=&
                sum(emissivity_pasch_beta(x,y,:))*l_cell*nc_parsec 
           
           
           write(1000+100*HYDROGEN+7,*)&
                sum(emissivity_brack_alpha(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_brack_alpha(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(40525*1E-8))/nc_sphere,&
                sum(emissivity_brack_alpha(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_brack_alpha(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_brack_alpha_z)&
                max_emissivity_brack_alpha_z=&
                sum(emissivity_brack_alpha(x,y,:))*l_cell*nc_parsec 
                      
           write(1000+100*HYDROGEN+8,*)&
                sum(emissivity_brack_beta(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_brack_beta(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(26259*1E-8))/nc_sphere,&
                sum(emissivity_brack_beta(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_brack_beta(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_brack_beta_z)&
                max_emissivity_brack_beta_z=&
                sum(emissivity_brack_beta(x,y,:))*l_cell*nc_parsec 
!
!
           write(1000+100*HELIUM+1,*)&
                sum(emissivity_hei_optical(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_hei_optical(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(5876*1E-8))/nc_sphere,&
                sum(emissivity_hei_optical(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_hei_optical(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_hei_optical_z)&
                max_emissivity_hei_optical_z=&
                sum(emissivity_hei_optical(x,y,:))*l_cell*nc_parsec 
!
!
           write(1000+100*CARBON+1,*)&
                sum(emissivity_cii_1335(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_cii_1335(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(2325*1E-8))/nc_sphere,&
                sum(emissivity_cii_1335(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
 
           if (sum(emissivity_cii_1335(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_cii_1335_z)&
                max_emissivity_cii_1335_z=&
                sum(emissivity_cii_1335(x,y,:))*l_cell*nc_parsec 

           write(1000+100*CARBON+2,*)&
                sum(emissivity_cii_2325(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_cii_2325(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(1335*1E-8))/nc_sphere,&
                sum(emissivity_cii_2325(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_cii_2325(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_cii_2325_z)&
                max_emissivity_cii_2325_z=&
                sum(emissivity_cii_2325(x,y,:))*l_cell*nc_parsec 

           write(1000+100*CARBON+3,*)&
                sum(emissivity_ciii_1907(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_ciii_1907(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(1908*1E-8))/nc_sphere,&
                sum(emissivity_ciii_1907(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_ciii_1907(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_ciii_1907_z)&
                max_emissivity_ciii_1907_z=&
                sum(emissivity_ciii_1907(x,y,:))*l_cell*nc_parsec 
!
!          
           write(1000+100*NITROGEN+1,*)&
                sum(emissivity_nii_optical(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_nii_optical(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(6548+6584)*1E-8))/nc_sphere,&
                sum(emissivity_nii_optical(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_nii_optical(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_optical_z)&
                max_emissivity_nii_optical_z=&
                sum(emissivity_nii_optical(x,y,:))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+2,*)&
                sum(emissivity_nii_auroral(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_nii_auroral(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(5755*1E-8))/nc_sphere,&
                sum(emissivity_nii_auroral(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_nii_auroral(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_auroral_z)&
                max_emissivity_nii_auroral_z=&
                sum(emissivity_nii_auroral(x,y,:))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+3,*)&
                sum(emissivity_nii_infrared122(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_nii_infrared122(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(121.89*1E-4))/nc_sphere,&
                sum(emissivity_nii_infrared122(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_nii_infrared122(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_nii_infrared122_z)&
                max_emissivity_nii_infrared122_z=&
                sum(emissivity_nii_infrared122(x,y,:))*l_cell*nc_parsec 

           write(1000+100*NITROGEN+4,*)&
                sum(emissivity_niii_infrared57(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_niii_infrared57(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(57.343*1E-4))/nc_sphere,&
                sum(emissivity_niii_infrared57(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_niii_infrared57(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_niii_infrared57_z)&
                max_emissivity_niii_infrared57_z=&
                sum(emissivity_niii_infrared57(x,y,:))*l_cell*nc_parsec 
!           
!           
           write(1000+100*OXYGEN+1,*)&
                sum(emissivity_oi_optical(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oi_optical(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( (1./2.)*(6300+6363)*1E-8))/nc_sphere,&
                sum(emissivity_oi_optical(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oi_optical(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oi_optical_z)&
                max_emissivity_oi_optical_z=&
                sum(emissivity_oi_optical(x,y,:))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+2,*)&
                sum(emissivity_oii_optical(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oii_optical(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 3727*1E-8))/nc_sphere,&
                sum(emissivity_oii_optical(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oii_optical(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oii_optical_z)&
                max_emissivity_oii_optical_z=&
                sum(emissivity_oii_optical(x,y,:))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+3,*)&
                sum(emissivity_oii_7320(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oii_7320(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 3727*1E-8))/nc_sphere,&
                sum(emissivity_oii_7320(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oii_7320(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oii_7320_z)&
                max_emissivity_oii_7320_z=&
                sum(emissivity_oii_7320(x,y,:))*l_cell*nc_parsec 
                      
           write(1000+100*OXYGEN+4,*)&
                sum(emissivity_oiii_optical(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oiii_optical(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(5007+4959)*1E-8))/nc_sphere,&
                sum(emissivity_oiii_optical(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oiii_optical(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_optical_z)&
                max_emissivity_oiii_optical_z=&
                sum(emissivity_oiii_optical(x,y,:))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+5,*)&
                sum(emissivity_oiii_auroral(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oiii_auroral(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(4363*1E-8))/nc_sphere,&
                sum(emissivity_oiii_auroral(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oiii_auroral(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_auroral_z)&
                max_emissivity_oiii_auroral_z=&
                sum(emissivity_oiii_auroral(x,y,:))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+6,*)&
                sum(emissivity_oiii_infrared52(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oiii_infrared52(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(51.814*1E-4))/nc_sphere,&
                sum(emissivity_oiii_infrared52(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oiii_infrared52(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_infrared52_z)&
                max_emissivity_oiii_infrared52_z=&
                sum(emissivity_oiii_infrared52(x,y,:))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+7,*)&
                sum(emissivity_oiii_infrared88(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oiii_infrared88(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(88.356*1E-4))/nc_sphere,&
                sum(emissivity_oiii_infrared88(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oiii_infrared88(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oiii_infrared88_z)&
                max_emissivity_oiii_infrared88_z=&
                sum(emissivity_oiii_infrared88(x,y,:))*l_cell*nc_parsec 

           write(1000+100*OXYGEN+8,*)&
                sum(emissivity_oiv_infrared(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_oiv_infrared(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(25.916*1E-4))/nc_sphere,&
                sum(emissivity_oiv_infrared(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_oiv_infrared(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_oiv_infrared_z)&
                max_emissivity_oiv_infrared_z=&
                sum(emissivity_oiv_infrared(x,y,:))*l_cell*nc_parsec 
!
!
           write(1000+100*NEON+1,*)&
                sum(emissivity_neii_infrared13(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_neii_infrared13(x,y, :))*l_cell*nc_parsec/&
                (nc_planck*nc_light/( 12.814*1E-4))/nc_sphere,&
                sum(emissivity_neii_infrared13(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_neii_infrared13(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_neii_infrared13_z)&
                max_emissivity_neii_infrared13_z=&
                sum(emissivity_neii_infrared13(x,y,:))*l_cell*nc_parsec 
                      
           write(1000+100*NEON+2,*)&
                sum(emissivity_neiii_optical(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_neiii_optical(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(3869+3968)*1E-8))/nc_sphere,&
                sum(emissivity_neiii_optical(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_neiii_optical(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_neiii_optical_z)&
                max_emissivity_neiii_optical_z=&
                sum(emissivity_neiii_optical(x,y,:))*l_cell*nc_parsec 

           write(1000+100*NEON+3,*)&
                sum(emissivity_neiii_infrared16(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_neiii_infrared16(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(15.555*1E-4))/nc_sphere,&
                sum(emissivity_neiii_infrared16(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_neiii_infrared16(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_neiii_infrared16_z)&
                max_emissivity_neiii_infrared16_z=&
                sum(emissivity_neiii_infrared16(x,y,:))*l_cell*nc_parsec 
!
!
           write(1000+100*SULFUR+1,*)&
                sum(emissivity_sii_4068(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_sii_4068(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(4068+4076)*1E-8))/nc_sphere,&
                sum(emissivity_sii_4068(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_sii_4068(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_sii_4068_z)&
                max_emissivity_sii_4068_z=&
                sum(emissivity_sii_4068(x,y,:))*l_cell*nc_parsec 

           write(1000+100*SULFUR+2,*)&
                sum(emissivity_sii_6716(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_sii_6716(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(6716+6731)*1E-8))/nc_sphere,&
                sum(emissivity_sii_6716(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_sii_6716(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_sii_6716_z)&
                max_emissivity_sii_6716_z=&
                sum(emissivity_sii_6716(x,y,:))*l_cell*nc_parsec 
                      
           write(1000+100*SULFUR+3,*)&
                sum(emissivity_siii_optical(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_siii_optical(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/((1./2.)*(9532+9069)*1E-8))/nc_sphere,&
                sum(emissivity_siii_optical(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_siii_optical(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_optical_z)&
                max_emissivity_siii_optical_z=&
                sum(emissivity_siii_optical(x,y,:))*l_cell*nc_parsec 

           write(1000+100*SULFUR+4,*)&
                sum(emissivity_siii_infrared19(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_siii_infrared19(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(18.713*1E-4))/nc_sphere,&
                sum(emissivity_siii_infrared19(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_siii_infrared19(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_infrared19_z)&
                max_emissivity_siii_infrared19_z=&
                sum(emissivity_siii_infrared19(x,y,:))*l_cell*nc_parsec 

           write(1000+100*SULFUR+5,*)&
                sum(emissivity_siii_infrared34(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_siii_infrared34(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(33.47*1E-4))/nc_sphere,&
                sum(emissivity_siii_infrared34(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_siii_infrared34(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_siii_infrared34_z)&
                max_emissivity_siii_infrared34_z=&
                sum(emissivity_siii_infrared34(x,y,:))*l_cell*nc_parsec 

           write(1000+100*SULFUR+6,*)&
                sum(emissivity_siv_infrared11(x,y,:))*l_cell*nc_parsec,&
                sum(emissivity_siv_infrared11(x,y,:))*l_cell*nc_parsec/&
                (nc_planck*nc_light/(10.514*1E-4))/nc_sphere,&
                sum(emissivity_siv_infrared11(x,y,:))*l_cell*nc_parsec/&
                zHbeta*maskedvalue(1.0_dp,em_z,em_min)
           
           if (sum(emissivity_siv_infrared11(x,y,:))*l_cell*nc_parsec .gt.&
                max_emissivity_siv_infrared11_z)&
                max_emissivity_siv_infrared11_z=&
                sum(emissivity_siv_infrared11(x,y,:))*l_cell*nc_parsec 
!          
           write(10000,*) sum (electrons(x, y, :)**2)*l_cell
        end do
     end do
   end subroutine write_maps_z

!
!  CAUTION: THESE ROUTINES CURRENTLY ARE NOT CALLED 
!           FURTHERMORE NOT THE ENTIRE NUMBER OF LINES 
!           HAS BEEN INCLUDED
!

   subroutine write_max_emission_data_x
     implicit none
     write(*,*)
     write(*,*) 'x-direction:'
     
     write(*,'(A17,2A30 )')&
          'Line', 'max. Linestrength in', 'max. photons s^-1 cm^-2 sr^-1'
     write(*,'(A17,2A30)')&
          '',      'erg s^-1 cm^-2 sr^-2',''
     write(*,'(A17,2A30 )')&
          '====', '====================', '============================='
     write(*,'(A17,2ES30.5)')&
          'H_alpha', max_emissivity_h_alpha_x/nc_sphere, &
          max_emissivity_h_alpha_x/(nc_planck*nc_light/(6563*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_beta', max_emissivity_h_beta_x/nc_sphere, &
          max_emissivity_h_beta_x/(nc_planck*nc_light/(4861*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_gamma', max_emissivity_h_gamma_x/nc_sphere, &
          max_emissivity_h_gamma_x/(nc_planck*nc_light/(4340*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_delta', max_emissivity_h_delta_x/nc_sphere, &
          max_emissivity_h_delta_x/(nc_planck*nc_light/(4102*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'pasch_alpha', max_emissivity_pasch_alpha_x/nc_sphere, &
          max_emissivity_pasch_alpha_x/(nc_planck*nc_light/(18745*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'pasch_beta', max_emissivity_pasch_alpha_x/nc_sphere, &
          max_emissivity_pasch_beta_x/(nc_planck*nc_light/(12814*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'brack_alpha', max_emissivity_brack_alpha_x/nc_sphere, &
          max_emissivity_brack_alpha_x/(nc_planck*nc_light/(40525*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'brack_beta', max_emissivity_brack_beta_x/nc_sphere, &
          max_emissivity_brack_beta_x/(nc_planck*nc_light/(26259*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[OII] 3726+3729', max_emissivity_oii_optical_x/nc_sphere, &
          max_emissivity_oii_optical_x/(nc_planck*nc_light/( 3727*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[OIII] 5007+4959', max_emissivity_oiii_optical_x/nc_sphere, &
          max_emissivity_oiii_optical_x/(nc_planck*nc_light/(5000*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[NII] 6584+6548', max_emissivity_nii_optical_x/nc_sphere, &
          max_emissivity_nii_optical_x/(nc_planck*nc_light/( 6535*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[SIII] 9531+9069', max_emissivity_siii_optical_x/nc_sphere, &
          max_emissivity_siii_optical_x/(nc_planck*nc_light/(9400*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[OIV] 25.913 u', max_emissivity_oiv_infrared_x/nc_sphere, &
          max_emissivity_oiv_infrared_x/(nc_planck*nc_light/(259130*1e-8))/nc_sphere
     
   end subroutine write_max_emission_data_x

   subroutine write_max_emission_data_y
     implicit none
     write(*,*)

     write(*,*) 'y-direction:'
     
     write(*,'(A17,2A30 )')&
          'Line', 'max. Linestrength in', 'max. photons s^-1 cm^-2 sr^-1'
     write(*,'(A17,2A30)')&
          '',      'erg s^-1 cm^-2 sr^-2',''
     write(*,'(A17,2A30 )')&
          '====', '====================', '============================='
     write(*,'(A17,2ES30.5)')&
          'H_alpha', max_emissivity_h_alpha_y/nc_sphere, &
          max_emissivity_h_alpha_y/(nc_planck*nc_light/(6563*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_beta', max_emissivity_h_beta_y/nc_sphere, &
          max_emissivity_h_beta_y/(nc_planck*nc_light/(4861*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_gamma', max_emissivity_h_gamma_y/nc_sphere, &
          max_emissivity_h_gamma_y/(nc_planck*nc_light/(4340*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_delta', max_emissivity_h_delta_y/nc_sphere, &
          max_emissivity_h_delta_y/(nc_planck*nc_light/(4102*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'pasch_alpha', max_emissivity_pasch_alpha_y/nc_sphere, &
          max_emissivity_pasch_alpha_y/(nc_planck*nc_light/(18745*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'pasch_beta', max_emissivity_pasch_alpha_y/nc_sphere, &
          max_emissivity_pasch_beta_y/(nc_planck*nc_light/(12814*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'brack_alpha', max_emissivity_brack_alpha_y/nc_sphere, &
          max_emissivity_brack_alpha_y/(nc_planck*nc_light/(40525*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'brack_beta', max_emissivity_brack_beta_y/nc_sphere, &
          max_emissivity_brack_beta_y/(nc_planck*nc_light/(26259*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[OII] 3726+3729', max_emissivity_oii_optical_y/nc_sphere, &
          max_emissivity_oii_optical_y/(nc_planck*nc_light/( 3727*1e-8))&
          /nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[OIII] 5007+4959', max_emissivity_oiii_optical_y/nc_sphere, &
          max_emissivity_oiii_optical_y/(nc_planck*nc_light/(5000*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[NII] 6584+6548', max_emissivity_nii_optical_y/nc_sphere, &
          max_emissivity_nii_optical_y/(nc_planck*nc_light/( 6535*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[SIII] 9531+9069', max_emissivity_siii_optical_y/nc_sphere, &
          max_emissivity_siii_optical_y/(nc_planck*nc_light/(9400*1e-8))/&
          nc_sphere
     

     write(*,'(A17,2ES30.5)')&
          '[OIV] 25.9 u', max_emissivity_oiv_infrared_y/nc_sphere, &
          max_emissivity_oiv_infrared_y/(nc_planck*nc_light/(259130*1e-8))/&
          nc_sphere
     
   end subroutine write_max_emission_data_y

   subroutine write_max_emission_data_z
     implicit none
     write(*,*)

     write(*,*) 'z-direction:'
     
     write(*,'(A17,2A30 )')&
          'Line', 'max. Linestrength in', 'max. photons s^-1 cm^-2 sr^-1'
     write(*,'(A17,2A30)')&
          '',      'erg s^-1 cm^-2 sr^-2',''
     write(*,'(A17,2A30 )')&
          '====', '====================', '============================='
     write(*,'(A17,2ES30.5)')&
          'H_alpha', max_emissivity_h_alpha_z/nc_sphere, &
          max_emissivity_h_alpha_z/(nc_planck*nc_light/(6563*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_beta', max_emissivity_h_beta_z/nc_sphere, &
          max_emissivity_h_beta_z/(nc_planck*nc_light/(4861*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_gamma', max_emissivity_h_gamma_z/nc_sphere, &
          max_emissivity_h_gamma_z/(nc_planck*nc_light/(4340*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'H_delta', max_emissivity_h_delta_z/nc_sphere, &
          max_emissivity_h_delta_z/(nc_planck*nc_light/(4102*1e-8))/nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'pasch_alpha', max_emissivity_pasch_alpha_z/nc_sphere, &
          max_emissivity_pasch_alpha_z/(nc_planck*nc_light/(18745*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'pasch_beta', max_emissivity_pasch_alpha_z/nc_sphere, &
          max_emissivity_pasch_beta_z/(nc_planck*nc_light/(12814*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'brack_alpha', max_emissivity_brack_alpha_z/nc_sphere, &
          max_emissivity_brack_alpha_z/(nc_planck*nc_light/(40525*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          'brack_beta',max_emissivity_brack_beta_z/nc_sphere, &
          max_emissivity_brack_beta_z/(nc_planck*nc_light/(26259*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[OII] 3726+3729', max_emissivity_oii_optical_z/nc_sphere, &
          max_emissivity_oii_optical_z/(nc_planck*nc_light/( 3727*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[OIII] 5007+4959', max_emissivity_oiii_optical_z/nc_sphere, &
          max_emissivity_oiii_optical_z/(nc_planck*nc_light/(5000*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[NII] 6584+6548', max_emissivity_nii_optical_z/nc_sphere, &
          max_emissivity_nii_optical_z/(nc_planck*nc_light/( 6535*1e-8))/&
          nc_sphere
     
     write(*,'(A17,2ES30.5)')&
          '[SIII] 9531+9069', max_emissivity_siii_optical_z/nc_sphere, &
          max_emissivity_siii_optical_z/(nc_planck*nc_light/(9400*1e-8))/&
          nc_sphere
     

     write(*,'(A17,2ES30.5)')&
          '[OIV] 25.9 u', max_emissivity_oiv_infrared_z/nc_sphere, &
          max_emissivity_oiv_infrared_z/(nc_planck*nc_light/(259130*1e-8))/&
          nc_sphere
   end subroutine write_max_emission_data_z


end program nebular_statistics
