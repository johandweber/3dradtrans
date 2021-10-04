!=============================================================================
program DISKCUT
! Computes the star formation rate within an sqare "cut" of the plane of a 
! galactic disk (including the space "above" and "below" the disk)
! It assumes an exponential decay of the star fromation density of a galaxy as
! a function of radius and does not explicitly consider substructure like 
! spiral arms.
!
  use M_data_types
!
  implicit none
  integer :: x,y 
  integer::  resol=20000
  real(dp):: scale_length, total_sfr, size_of_cut, center_of_cut
  real(dp):: cut_sfr=0._dp
!=============================================================================
!
  write(*,*) 'Please input characteristic scale length of galactic dist in kpc :'
  read( *,*) scale_length
  write(*,*)
  write(*,*) 'Please input total str formation rate of galaxy in M_sun/yr. :'
  read (*,*) total_sfr
  write(*,*)
  write(*,*) 'Please input distance from center of the cut to'
  write(*,*) 'galactic center in kpc'
  read (*,*) center_of_cut 
  write(*,*) 
  write(*,*) 'Please input the size of the extract in kpc'
  read (*,*)  size_of_cut
  write(*,*)
  write(*,*) 'Computing the SFR within the cut, please wait...'
  do y=-resol,+resol
     do x=-resol, +resol
        cut_sfr=cut_sfr+sfr_density(&
             sqrt((center_of_cut+x*size_of_cut/(2*resol+1))**2+&
             (y*size_of_cut/(2*resol+1))**2),&
             scale_length,total_sfr)*&
             (size_of_cut/(2*resol+1))**2
     end do
  end do
  write(*,*) 'The star formation rate within the cut in M_sun/yr is'
  write(*,*) cut_sfr
!
  contains 
!
    real(dp) function sfr_density(distance, sc_len,total)
      implicit none
      real(dp) :: distance, sc_len, total
      sfr_density=total/(sc_len**2*(8._dp*atan(1._dp)))*exp(-distance/sc_len)
    end function sfr_density
!
end program DISKCUT
