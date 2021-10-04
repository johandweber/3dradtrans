!=============================================================================
! Continua - This program computes the stellar and nebular continuua based
! on the output from occup_3d_ext 
!=============================================================================
!
! MOD - jweber 2018-12
!

program continua
  use M_data_types
  use M_natural_constants   , only: nc_pi, nc_sol_rad
  use M_definitions         , only: l_cell,source_filename,&
                                    num_sources,&
                                    source_x, source_y, source_z, R_RS,&
                                    start_activity, lifetime,  spectra,&
                                    load_parameters
  use M_data_input          , only: input, wlength, spectrum, row_end,nu
  implicit none
  integer(i4b)                   :: source_counter
  real(dp)                       :: stellar_luminosity, total_luminosity
  integer(i4b)                   :: wavelenght_counter
  integer(i4b), parameter        :: length_of_wavelenght_list=4
  real(dp), dimension (1:length_of_wavelenght_list) :: wavelenght_list=&
                                   (/6563._dp, 4861._dp, 4340._dp, 4101._dp/)
!=============================================================================
!
  call load_parameters
!
  do source_counter=1,num_sources
     call input(spectra(source_counter))
     write(*,*) "source #", source_counter
     write(*,*) "x-coordinate :", source_x (source_counter)
     write(*,*) "y-coordinate :", source_x (source_counter)
     write(*,*) "z-coordinate :", source_z (source_counter)
     write(*,*) "radius/R-sun :", R_RS     (source_counter)
     write(*,*) "spectrum:",      spectra  (source_counter)
     write(*,*) 
     do wavelenght_counter=1,length_of_wavelenght_list
        call write_stellar_continua&
             (wavelenght_list(wavelenght_counter),&
             row_end,&
             wlength(1:row_end),&
             spectrum, R_RS(source_counter))
     end do
  end do
!  
  contains

    
    subroutine write_stellar_continua(lambda, no_of_wl, lambda_grid ,H, radius)
      implicit none
      integer(i4b) :: no_of_wl
      real(dp), dimension (1:no_of_wl) :: lambda_grid, H
      real(dp) :: lambda, radius
!
      radius=radius*nc_sol_rad

      write(*,*) 'stellar continuum luminosity in erg s^-1 Hz^-1 :'
      write(*,*)
      write(*,*) 'stellar continuum luminosity in erg s^-1 A^-1  :'
    end subroutine write_stellar_continua

    
    real(dp) function calc_stellar_continuum&
            (lambda, lambda_grid, no_of_wl, H, radius)
      implicit none
      integer(i4b) :: no_of_wl
      real(dp), dimension (1:no_of_wl) :: lambda_grid, H
      real(dp)     :: lambda, radius
    end function calc_stellar_continuum
!  
end program continua
