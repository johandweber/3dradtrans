!=============================================================================
!                         ****Module M_subsystems
! Contains routines to the occupation numbers of different excited levels
! of the same ion and the resulting line ratios.
! This is not intended to be a full replacement of the time-dependent,
! solution based on the ionization stages as it is implemented 
! (several times...) in hydrogen.f90
! The routines provided here might be used in occup_3d_ext or
! in other routines.
! 
!=============================================================================
!
module M_subsystems
!---- routines:
!     SUBSYSTEM_OII
!
! MOD JWEBER 07.16
  use M_data_types
!
  implicit none
  
!
contains
!
!=============================================================================
  subroutine SOLVE_MATRIX  (n, rate_matrix, initial,&
                           occupation, t)
!
! Computes the occupation numbers for for a n-state system, where the sum of
! all occupation numbers is conserved.
! This is a mathematical function and does not compute the rate
! coefficients, which have to be computed by anotther part of the
! program that calls this subroutine.
!
    implicit none
!   input
    integer(i4b)                    :: n
    real(dp), dimension (n,n)       :: rate_matrix
    real(dp), dimension(n)          :: initial
    real(dp)                        :: t
!   output
    real(dp), dimension(n)         :: occupation
!   local 
    integer(i4b)                    :: diagonal_counter, inner_counter,&
                                       outer_counter
    integer                         :: info=0
    real(dp), dimension(:,:),allocatable :: leftvectors,rightvectors,&
                                            rightvectors_copy
    real(dp), dimension(:),allocatable   :: eigvalre, eigvalim, work, start
    integer(i4b), dimension(:), allocatable      :: pivots

    real(dp), dimension(1)          :: maxeig   
!=============================================================================
!
    allocate(leftvectors(n,n), rightvectors(n,n), rightvectors_copy(n,n))
    allocate(eigvalre(n), eigvalim(n), work(50*n), start(n), pivots(n))
!    
    eigvalre=0._dp
    eigvalim=0._dp
    leftvectors=0._dp
    rightvectors=0._dp
    rightvectors_copy=0._dp
!
    call dgeev('N','V', n,rate_matrix,n,eigvalre, eigvalim,&
         leftvectors,n,rightvectors,n,work, 50, info)
!
    rightvectors_copy=rightvectors
!
    call dgesv(n,1,rightvectors_copy,n,pivots,&
              initial,n,info)
    occupation=0._dp
    do outer_counter=1, n
       do inner_counter=1,n 
          occupation(outer_counter)=occupation(outer_counter)+&
               initial(inner_counter)*&
               rightvectors(outer_counter, inner_counter)*&
               exp(eigvalre(inner_counter)*t)   
       end do
    end do
!
    deallocate(leftvectors, rightvectors, rightvectors_copy)
    deallocate(eigvalre, eigvalim, work, start)

 end subroutine SOLVE_MATRIX

!==============================================================================
 subroutine GENERATE_MATRIX(n, einstein, collisional, rate_matrix,&
                            statweight, energy, ne, temp)
   !
   ! Computesthe rate matrix for a system of the ground state and several
   ! excited states of a model atom
   ! It currently considers collisional excitation and de-excitation rates
   ! and as well as the spontaneous radiative decay
!
   use M_data_types
   use M_natural_constants,   only: nc_pi, nc_planck, nc_boltzmann, nc_light
!
   implicit none
!
!  input
   integer(i4b)               :: n
   real(dp)                   :: ne, temp
!     statistical weights of the levels, energy above the ground level
   real(dp), dimension(n)     :: statweight, energy
!      einstein contains the einstein coefficients, collisional the collision, 
!      strengths
   real(dp), dimension(n,n)   :: einstein, collisional
!  output
   real(dp), dimension(n,n)   :: rate_matrix
!  local
   integer(i4b)               :: columncounter, linecounter
   real(dp)                   :: deltaE, weightupper, weightlower
!==============================================================================
   rate_matrix=einstein
!  upward collisional rates   
   do linecounter=2, n
      do columncounter=1, linecounter-1
         deltaE= nc_planck*(energy(linecounter)-energy(columncounter))*nc_light
         weightlower=statweight(columncounter)
         rate_matrix(linecounter,columncounter)=&
              rate_matrix(linecounter, columncounter)+&
              ne*8.629e-6/sqrt(temp)*collisional(columncounter, linecounter)/&
              weightlower*exp(-deltaE/(nc_boltzmann*temp))        
      end do
   end do
!  downward collisional rates
   do linecounter=1,n-1
      do columncounter=linecounter+1, n
         weightupper=statweight(columncounter)
         rate_matrix(linecounter, columncounter)=&
              rate_matrix(linecounter, columncounter)+&
              ne*8.629e-6/sqrt(temp)*collisional(columncounter, linecounter)/&
              weightlower                     
      end do
   end do

   do linecounter=1,n
      rate_matrix(linecounter, linecounter)=&
      -sum (rate_matrix(:,linecounter))
   end do

 end subroutine GENERATE_MATRIX


!
!==============================================================================
  subroutine COMPARE_EMISSION_PHOTONS (n, occupation, einstein, &
                                       emission_pho)
    implicit none
!   input
    integer(i4b)              :: n
    real(dp), dimension(n)    :: occupation
    real(dp), dimension(n,1)  :: occupation_matrix
    real(dp), dimension(n, n) :: einstein
!   output
    real(dp), dimension(n,n)  :: emission_pho  
!   local
    integer(i4b)              :: statecounter
!==============================================================================    
!
    do statecounter=1,n
       emission_pho(statecounter,:) = einstein(statecounter,:)*&
                                      occupation                                    
    end do
!
  end subroutine COMPARE_EMISSION_PHOTONS
!
!==============================================================================
  subroutine COMPARE_EMISSION_ENERGY (n, occupation, einstein, energy,&
                                       emission_en)
    implicit none
!   input
    integer(i4b)              :: n
    real(dp), dimension(n,n)  :: einstein
    real(dp), dimension(n)    :: occupation, energy
!   output
    real(dp), dimension(n,n)  :: emission_en
!   local
    integer(i4b)              :: linecounter, columncounter, statecounter
    real(dp), dimension(n,n)  :: energymatrix
!==============================================================================
!    
    do linecounter=1,n
       do columncounter=1,n
          energymatrix(linecounter, columncounter) =&
               energy(columncounter)-energy(linecounter)
       end do
    end do
    do statecounter=1,n
       emission_en(statecounter,:) = einstein(statecounter,:)*&
            energymatrix(statecounter,:)*occupation
    end do
  end subroutine COMPARE_EMISSION_ENERGY
!
end module M_subsystems

!=============================================================================
