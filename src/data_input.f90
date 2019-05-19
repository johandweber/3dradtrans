!===============================================================================
!                                 ****MODULE DATA_INPUT****
! reads file with spectrum and transforms it to new grid under conservation of 
! flux start with subroutine "data_transform"
!==============================================================================

module M_data_input
!--- routines:
!    DATA_TRANSFORM
!    INPUT
!    NU_GRID
!    INPUT_TRANSFORM
!
  use M_data_types 
  use M_definitions,                only: points
!
  implicit none
!  
  real   (dp ), dimension (:), allocatable   :: wlength, spectrum, nu 
  real   (dp ), dimension (0:points+5)       :: weights_nu, nug,nugh,nugh_spect,&
                                                nug_spect
  integer(i4b) :: fp, row_end=2900 
!  
contains
    
!================================================================================ 
  subroutine DATA_TRANSFORM(filename)
!
!  Is responsible for reading the descriptions of the SEDs and transfroming them
!  to the internal frequency grid of the program (with a reduced number of 
!  frequency points)
!
! called by main program 
!
  implicit none
!
    character (len=3000) :: filename
!================================================================================
!
    call input(filename)
    call nu_grid
    call input_transform
  end subroutine DATA_TRANSFORM

      

!================================================================================ 
  subroutine INPUT(filename)
!
!    Reads the spectra
!
!   called by DATA_TRANSFORM
!       
    
  use M_definitions,              only: mergespec
  use M_natural_constants,        only: nc_planck, nc_light
   
!
  implicit none
!
  integer (i4b)           :: indx           ! "index" is an intrinsic function
  character (len=3000)    :: filename  
!  
  integer,  dimension(10) :: dump_int_arguments
  real(dp), dimension(10) :: dump_real_arguments
!================================================================================
!    
    print*, "start input"
!        
    if (mergespec) then
       print*, "using mergespec"
       open(2, file=adjustl(trim(filename)))
       read(2,*) dump_int_arguments(1),row_end, dump_real_arguments(1),&
            dump_int_arguments(2), dump_int_arguments(3),&
            dump_real_arguments(2), dump_real_arguments(3), dump_real_arguments(4)
!
       if (allocated(wlength)) then
          deallocate(wlength)
          deallocate(spectrum)
          deallocate(nu)
       end if
!       
       if (.not. allocated(wlength)) then
          allocate (wlength (0:row_end))
          allocate (spectrum (0: row_end))
          allocate (nu (0:row_end))
       end if
 !      
       wlength=0
       nu=0
       spectrum=0
       indx=1
       do fp=0,1
          read(2,*)                          ! empty lines
       end do
!       
       do fp=1, row_end
          read(2,*) wlength(fp),spectrum(fp), dump_real_arguments(5),&
               dump_real_arguments(6)
          if(wlength(fp).lt.0.5E+5_dp) then
             nu(indx)=nc_light*1E+8_dp/wlength(fp)
             spectrum(indx)=spectrum(fp)/(nc_planck*nu(indx))
             indx=indx+1
          end if
       end do
!       
       row_end=indx-1
       close(2)
!      
    else
!
       row_end=2900  
!
       if (.not. allocated(wlength)) then
          allocate (wlength (0:row_end))
          allocate (spectrum (0: row_end))
          allocate (nu (0:row_end))
          wlength=0
          nu=0
          spectrum=0
          indx=1
       end if
!       
       open(2,file=adjustl(trim(filename)))
       do fp=0,15
          read(2,*)                         ! empty lines
       end do
!       
       indx=1
       do fp=0,row_end
          read(2,*) wlength(fp), spectrum(fp)
          if(wlength(fp).lt.5.1E+5_dp) then
             nu(indx)=nc_light*1E+8_dp/wlength(fp)
             spectrum(indx)=spectrum(fp)/(nc_planck*nu(indx))
             indx=indx+1
          end if
       end do
!       
       row_end=indx-1
       close(2)
    endif
!
  end subroutine INPUT
!    
!===============================================================================  
  subroutine NU_GRID
!
!    Creates the frequency grid.
!
!   called by DATA_TRANSFORM
!
    use M_natural_constants,    only: nc_light, nc_lyman_edge
    use M_definitions,          only: points, lower_energy_limit
!
  implicit none
!
    integer (i4b)  :: p
!===============================================================================
!    
    ! print*, "nu_grid started"
    nug=0
    nugh=0
    nug_spect=0
    nugh_spect=0
    do p=1,points+5
       nug(p)=nc_light/((nc_light/(lower_energy_limit*nc_lyman_edge))/points*p)
    end do
 !   
    !half grid nugh
    nugh(1)=nug(1)
    do p=2, points
       nugh(p)=(nug(p-1)+nug(p))/2._dp
    end do
    nugh(points+1)=nug(points)
! 
 end subroutine NU_GRID



  
!================================================================================ 
 subroutine INPUT_TRANSFORM
!
!    Remaps the read spectra to the internal frequency grid
!
! called by DATA_TRANSFORm
!
   use M_natural_constants, only: nc_light
   use M_definitions,       only: points
!
 implicit none 
!   
   real   (dp )                          :: nugh_under, nugh_over 
   integer(i4b)                          :: p
!
   real   ( dp), dimension (0:points)    :: interval
   integer(i4b), dimension (0:points)    :: interval_fill
!================================================================================
!  
   nugh_spect(1)=spectrum(row_end)
 !  
   do p=2,points+1
      do fp=0,row_end-1
         if(nu(fp).lt.nugh(p)) then
            nugh_over=spectrum(fp)
            nugh_under=spectrum(fp+1)
         end if
         if(nu(fp).eq.nugh(p)) then
            nugh_over=spectrum(fp)
            nugh_under=spectrum(fp)
         end if
      end do
!      
      nugh_spect(p)=(nugh_under+nugh_over)/2._dp
      !print*,p,nugh(p),nugh_spect(p)
   end do

!      
   fp=row_end-1
   interval_fill=0
   interval=0
   do p=1,points
      do while(nu(fp).gt.nugh(p+1))
         interval_fill(p) = interval_fill(p)+1
         fp=fp-1
      end do
!      print*,'interval_fill', p,interval_fill(p)
   end do
!   
   fp=row_end-1
   p=0
   do while(p.lt.points)
      p=p+1
      if(interval_fill(p).eq.0) then
         interval(p) = (nugh_spect(p)+nugh_spect(p+1))/2._dp*abs(nugh(p)-nugh(p+1))
         print*
         print*, "a"
      else
         if(interval_fill(p).eq.1) then
            interval(p) = (nugh_spect(p)+spectrum(fp))/2._dp*abs(nugh(p)-nu(fp))
            interval(p) = interval(p)+(spectrum(fp)+nugh_spect(p+1))/2._dp*&
                 abs(nu(fp)-nugh(p+1))
            print*
            print*, "b"
         else
            interval(p) = (nugh_spect(p)+spectrum(fp))/2._dp*abs(nugh(p)-nu(fp))
            do while(nu(fp-1).gt.nugh(p+1))
               interval(p) = interval(p)+((spectrum(fp)+spectrum(fp-1))/2._dp)*&
                    abs(nu(fp)-nu(fp-1))
               fp=fp-1
            end do
            interval(p) = interval(p)+(spectrum(fp)+nugh_spect(p+1))/2._dp*&
                 abs(nu(fp)-nugh(p+1))  
            fp=fp-1
            print*
            print*, "c"
         end if
      end if
!      
      nug_spect(p) = interval(p)/(nugh(p)-nugh(p+1))
      if (nug_spect(p) .lt. 0) then
         print*
         print*, "NEGATIVE FLUX!"
         print*, "Wavelength:"   ,1e+8_dp*nc_light/nug(p)
         print*, "Flux:", nug_spect(p)
         print*, "p:",p
         print*, "nugh(p)",nugh(p)
         print*, "nugh(p+1)",nugh(p+1)
         print*, "interval_fill", interval_fill(p)
         print*, "interval(p)", interval(p)
         print*, "interval(p)", interval(p)
         stop
      else
         print*
         print*, "Wavelength:"   ,1E+8_dp*nc_light/nug(p)
         print*, "Flux:", nug_spect(p)
         print*, "p:",p
         print*, "nugh(p)",nugh(p)
         print*, "nugh(p+1)",nugh(p+1)
         print*, "interval_fill", interval_fill(p)
         print*, "interval(p)", interval(p)         
      end if
   end do
!   
! calculate weights
   weights_nu(1)=nug(1)-nug(2)
   do p=2,points-1
      weights_nu(p)=nug(p)-nug(p+1)
   end do
!   
 end subroutine INPUT_TRANSFORM
!-----------------------------------------------------------------------
      
end module M_data_input
    
    
