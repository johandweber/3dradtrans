!===============================================================================
!                          *****MODULE NATURAL_CONSTANTS****
!===============================================================================
! This is a collection of some important natural constants (given in the cgs 
! system)

module M_natural_constants
    
  use M_data_types
    
  real(dp), parameter :: nc_lyman_edge=1.0973731569E+05_dp*2.99792458E+10_dp
  real(dp), parameter :: nc_light =2.99792458E+10_dp
  real(dp), parameter :: nc_planck = 6.6260755E-27_dp
  real(dp), parameter :: nc_boltzmann = 1.380658E-16_dp
  real(dp), parameter :: nc_sol_lum = 3.826E+33_dp
  real(dp), parameter :: nc_sol_rad = 6.9599E+10_dp
  real(dp), parameter :: nc_parsec= 3.08567808E+18_dp
  real(dp), parameter :: nc_h_ion_cross  = 6.3E-18_dp
  real(dp), parameter :: nc_pi=3.141592653589793116_dp
  real(dp), parameter :: nc_sphere=4*nc_pi
  real(dp), parameter :: nc_electronmass=9.10938291e-28_dp
  real(dp), parameter :: nc_megaparsec=nc_parsec*1E+6_dp
  real(dp), parameter :: nc_km=1E+5_dp
  real(dp), parameter :: nc_yr=365*24*3600
  real(dp), parameter :: nc_myr=1e+6_dp*nc_yr
  real(dp), parameter :: nc_gyr=1e+9_dp*nc_yr
   
end module M_natural_constants
!===============================================================================
!                      *****END MODULE NATURAL_CONSTANTS****
!===============================================================================
