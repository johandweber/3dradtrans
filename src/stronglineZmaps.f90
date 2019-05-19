program stronglinezmaps
  use iso_fortran_env
  use M_data_types
  use M_natural_constants
  use M_strongline
!
  implicit none
!

  real(dp) :: lower_em_limit =1.0_dp

  real(dp) :: totalH, totalO

  character(len=200):: argument
  real(dp) :: dummy1, dummy2               ! variables that have to be read from files,
                                           ! but are never used by the program
  real(dp), dimension(:,:), allocatable :: nii_x, nii_y, nii_z,&
                                           oii_x, oii_y, oii_z,&
                                           oiii_x, oiii_y, oiii_z,&
                                           halpha_x, halpha_y, halpha_z,&
                                           hbeta_x, hbeta_y, hbeta_z,&
                                           em_x, em_y, em_z
  integer(i4b)                          :: x, y, z,&
                                           xmax, ymax, zmax

  if (command_argument_count() .ne. 1) then
     write(ERROR_UNIT,*) "Wrong number of command line options !"
     write(ERROR_UNIT,*)
     write(ERROR_UNIT,*) "correct usage"
     write(ERROR_UNIT,*) "stronglineZmaps <filenamewithaoutsuffix>"
     write(ERROR_UNIT,*)
     stop
  end if

  call get_command_argument(1,argument)

  open(unit=101, file=trim(argument)//".emissionmeasure_x.txt")
  open(unit=102, file=trim(argument)//".emissionmeasure_y.txt")
  open(unit=103, file=trim(argument)//".emissionmeasure_z.txt")

  open(unit=1101, file=trim(argument)//".H_alpha_x.txt")
  open(unit=1102, file=trim(argument)//".H_alpha_y.txt")
  open(unit=1103, file=trim(argument)//".H_alpha_z.txt")

  open(unit=1111, file=trim(argument)//".H_beta_x.txt")
  open(unit=1112, file=trim(argument)//".H_beta_y.txt")
  open(unit=1113, file=trim(argument)//".H_beta_z.txt")

  open(unit=7201, file=trim(argument)//".NII-optical_x.txt")
  open(unit=7202, file=trim(argument)//".NII-optical_y.txt")
  open(unit=7203, file=trim(argument)//".NII-optical_z.txt")

  open(unit=8201, file=trim(argument)//".OII-optical_x.txt")
  open(unit=8202, file=trim(argument)//".OII-optical_y.txt")
  open(unit=8203, file=trim(argument)//".OII-optical_z.txt")

  open(unit=8301, file=trim(argument)//".OIII-optical_x.txt")
  open(unit=8302, file=trim(argument)//".OIII-optical_y.txt")
  open(unit=8303, file=trim(argument)//".OIII-optical_z.txt")


!============================================================================!
!                       integration along the x-axis                         !
!============================================================================!

  read(101,*) xmax, ymax, zmax
  ! skip first lines of the other files for projections 
  ! along x-axis
  read(1101,*)
  read(1111,*)
  read(7201,*)
  read(8201,*)
  read(8301,*)


  allocate (em_x(1:ymax, 1:zmax))
  allocate (halpha_x(1:ymax, 1:zmax))
  allocate (hbeta_x(1:ymax, 1:zmax))
  allocate (nii_x(1:ymax, 1:zmax))
  allocate (oii_x(1:ymax, 1:zmax))
  allocate (oiii_x(1:ymax, 1:zmax))
  
  do z=1,zmax
     do y=1,ymax
        read(101,*) em_x(y,z)
        read(1101,*) halpha_x(y,z), dummy1, dummy2
        read(1111,*) hbeta_x(y,z), dummy1, dummy2
        read(7201,*) nii_x(y,z), dummy1, dummy2
        read(8201,*) oii_x(y,z) , dummy1, dummy2
        read(8301,*) oiii_x(y,z), dummy1, dummy2
     end do
  end do

  open (unit=201,file=trim(argument)//".M91_x.txt")
  write(201,*) 1, ymax, zmax
  do z=1,zmax
     do y=1,ymax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_x(y,z) .lt. Lower_em_limit .or. oiii_x(y,z) .eq. 0._dp) then
           write(201,*) -1._dp
        else
           write(201,*) M91(nii_x(y,z), oii_x(y,z), oiii_x(y,z), hbeta_x(y,z))
        end if
     end do
  end do
  close(201)
  
  open (unit=202,file=trim(argument)//".KD02_x.txt")
  write(202,*) 1, ymax, zmax
  do z=1,zmax
     do y=1,ymax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_x(y,z) .lt. Lower_em_limit .or. oiii_x(y,z) .eq. 0._dp) then
           write(202,*) -1._dp
        else
           write(202,*) KD02(nii_x(y,z), oii_x(y,z), oiii_x(y,z), hbeta_x(y,z))
        end if
     end do
  end do
  close(202)

 
  open (unit=203,file=trim(argument)//".KK04_x.txt")
  write(203,*) 1, ymax, zmax
  do z=1,zmax
     do y=1,ymax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_x(y,z) .lt. Lower_em_limit .or. oiii_x(y,z) .eq. 0._dp) then
           write(203,*) -1._dp
        else
           write(203,*) KK04(nii_x(y,z), oii_x(y,z), oiii_x(y,z), hbeta_x(y,z))
        end if
     end do
  end do
  close(203)

  open (unit=204,file=trim(argument)//".Z94_x.txt")
  write(204,*) 1, ymax, zmax
  do z=1,zmax
     do y=1,ymax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_x(y,z) .lt. Lower_em_limit .or. oiii_x(y,z) .eq. 0._dp) then
           write(204,*) -1._dp
        else
           write(204,*) Z94(nii_x(y,z), oii_x(y,z), oiii_x(y,z), hbeta_x(y,z))
        end if
     end do
  end do
  close(204)

  open (unit=205,file=trim(argument)//".P05_x.txt")
  write(205,*) 1, ymax, zmax
  do z=1,zmax
     do y=1,ymax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_x(y,z) .lt. Lower_em_limit .or. oiii_x(y,z) .eq. 0._dp) then
           write(205,*) -1._dp
        else
           write(205,*) P05(nii_x(y,z), oii_x(y,z), oiii_x(y,z), hbeta_x(y,z))
        end if
     end do
  end do
  close(205)


  open (unit=206,file=trim(argument)//".D02_x.txt")
  write(206,*) 1, ymax, zmax
  do z=1,zmax
     do y=1,ymax
        ! In this callibrtation the [OIII] line is not considered, so the case where
        ! there is [OIII] emission needs not to be handled
        if (em_x(y,z) .lt. Lower_em_limit) then
           write(206,*) -1._dp
        else
           write(206,*) D02(nii_x(y,z), halpha_x(y,z))
        end if
     end do
  end do
  close(206)

  write(*,*) "Line ratios:"
  write(*,'(A30,F10.4)') 'H_alpha/H_beta :',sum(halpha_x(:,:))/sum(hbeta_x(:,:))
  write(*,'(A30,F10.4)') '[NII]_nebular/H_beta :',&
                                           sum(nii_x(:,:))/sum(hbeta_x(:,:))
  write(*,'(A30,F10.4)') '[OII]_nearUV/H_beta :',&
                                           sum(oii_x(:,:))/sum(hbeta_x(:,:))
  write(*,'(A30,F10.4)') '[OIII]_nebular/H_beta :',&
                                           sum(oiii_x(:,:))/sum(hbeta_x(:,:))
  write(*,*)

  write(*,*)
  write(*,*) "Summary considering the entire fluxes, integrated along "//&
             "the x-axis:"
  write(*,*) "(The results should be independent of the direction,"//&
               "else there is an error !)"
  write(*,'(A20,A20, A20)') 'Callibration', 'n(O)/n(H)', '12+log10(n(O)/n(H))'
  write(*,'(A20,ES20.4, ES20.4)') 'M91',&
                                  Ztoratio(M91(sum(nii_x(:,:)),&
                                               sum(oii_x(:,:)),&
                                               sum(oiii_x(:,:)),&
                                               sum(hbeta_x(:,:)))),&
                                            M91(sum(nii_x(:,:)),&
                                               sum(oii_x(:,:)),&
                                               sum(oiii_x(:,:)),&
                                               sum(hbeta_x(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'KD02',&
                                  Ztoratio(KD02(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))),&
                                           KD02(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'KK04',&
                                  Ztoratio(KK04(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))),&
                                           KK04(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'Z94',&
                                  Ztoratio( Z94(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))),&
                                            Z94(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'P05',&
                                  Ztoratio( P05(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))),&
                                            P05(sum(nii_x(:,:)),&
                                                sum(oii_x(:,:)),&
                                                sum(oiii_x(:,:)),&
                                                sum(hbeta_x(:,:)))
  write(*,'(A20,ES20.4, ES20.4)')'D02',&
                                  Ztoratio( D02(sum(nii_x(:,:)),&
                                                sum(halpha_x(:,:)))),&
                                            D02(sum(nii_x(:,:)),&
                                                sum(halpha_x(:,:)))

!!$  deallocate (em_x)
!!$  deallocate (halpha_x)
!!$  deallocate (hbeta_x)
!!$  deallocate (nii_x)
!!$  deallocate (oii_x)
!!$  deallocate (oiii_x)


!=============================================================================!
!                       integration along the y-axis                          !
!=============================================================================!

  read(102,*) xmax, ymax, zmax
  ! skip first lines of the other files for projections 
  ! along x-axis
  read(1102,*)
  read(1112,*)
  read(7202,*)
  read(8202,*)
  read(8302,*)

  allocate (em_y(1:xmax, 1:zmax))
  allocate (halpha_y(1:xmax, 1:zmax))
  allocate (hbeta_y(1:xmax, 1:zmax))
  allocate (nii_y(1:xmax, 1:zmax))
  allocate (oii_y(1:xmax, 1:zmax))
  allocate (oiii_y(1:xmax, 1:zmax))
  
  do z=1,zmax
     do x=1,xmax
        read(102,*) em_y(x,z)
        read(1102,*) halpha_y(x,z), dummy1, dummy2
        read(1112,*) hbeta_y(x,z), dummy1, dummy2
        read(7202,*) nii_y(x,z), dummy1, dummy2
        read(8202,*) oii_y(x,z) , dummy1, dummy2
        read(8302,*) oiii_y(x,z), dummy1, dummy2
     end do
  end do

  open (unit=201,file=trim(argument)//".M91_y.txt")
  write(201,*) xmax,1, zmax
  do z=1,zmax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_y(x,z) .lt. Lower_em_limit .or. oiii_y(x,z) .eq. 0._dp) then
           write(201,*) -1._dp
        else
           write(201,*) M91(nii_y(x,z), oii_y(x,z), oiii_y(x,z), hbeta_y(x,z))
        end if
     end do
  end do
  close(201)
  
  open (unit=202,file=trim(argument)//".KD02_y.txt")
  write(202,*) xmax,1, zmax
  do z=1,zmax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_y(x,z) .lt. Lower_em_limit .or. oiii_y(x,z) .eq. 0._dp) then
           write(202,*) -1._dp
        else
           write(202,*) KD02(nii_y(x,z), oii_y(x,z), oiii_y(x,z), hbeta_y(x,z))
        end if
     end do
  end do
  close(202)

  open (unit=203,file=trim(argument)//".KK04_y.txt")
  write(203,*) xmax, 1,  zmax
  do z=1,zmax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_y(x,z) .lt. Lower_em_limit .or. oiii_y(x,z) .eq. 0._dp) then
           write(203,*) -1._dp
        else
           write(203,*) KK04(nii_y(x,z), oii_y(x,z), oiii_y(x,z), hbeta_y(x,z))
        end if
     end do
  end do
  close(203)

  open (unit=204,file=trim(argument)//".Z94_y.txt")
  write(204,*) xmax, 1, zmax
  do z=1,zmax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_y(x,z) .lt. Lower_em_limit .or. oiii_y(x,z) .eq. 0._dp) then
           write(204,*) -1._dp
        else
           write(204,*) Z94(nii_y(x,z), oii_y(x,z), oiii_y(x,z), hbeta_y(x,z))
        end if
     end do
  end do
  close(204)
  
  open (unit=205,file=trim(argument)//".P05_y.txt")
  write(205,*) xmax, 1, zmax
  do z=1,zmax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_y(x,z) .lt. Lower_em_limit .or. oiii_y(x,z) .eq. 0._dp) then
           write(205,*) -1._dp
        else
           write(205,*) P05(nii_y(x,z), oii_y(x,z), oiii_y(x,z), hbeta_y(x,z))
        end if
     end do
  end do
  close(205)


  open (unit=206,file=trim(argument)//".D02_y.txt")
  write(206,*) xmax, y, zmax
  do z=1,zmax
     do x=1,xmax
        ! In this callibrtation the [OIII] line is not considered, so the case where
        ! there is [OIII] emission needs not to be handled
        if (em_y(x,z) .lt. Lower_em_limit) then
           write(206,*) -1._dp
        else
           write(206,*) D02(nii_y(x,z), halpha_y(x,z))
        end if
     end do
  end do
  close(206)

  write(*,*)
  write(*,*) "Summary considering the entire fluxes, integrated along "//&
             "the y-axis:"
  write(*,*) "(The results should be independent of the direction,"//&
               "else there is an error !)"
  write(*,'(A20,A20, A20)') 'Callibration', 'n(O)/n(H)', '12+log10(n(O)/n(H))'
  write(*,'(A20,ES20.4, ES20.4)') 'M91',&
                                  Ztoratio(M91(sum(nii_y(:,:)),&
                                               sum(oii_y(:,:)),&
                                               sum(oiii_y(:,:)),&
                                               sum(hbeta_y(:,:)))),&
                                            M91(sum(nii_y(:,:)),&
                                               sum(oii_y(:,:)),&
                                               sum(oiii_y(:,:)),&
                                               sum(hbeta_y(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'KD02',&
                                  Ztoratio(KD02(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))),&
                                           KD02(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'KK04',&
                                  Ztoratio(KK04(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))),&
                                           KK04(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'Z94',&
                                  Ztoratio( Z94(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))),&
                                            Z94(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'P05',&
                                  Ztoratio( P05(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))),&
                                            P05(sum(nii_y(:,:)),&
                                                sum(oii_y(:,:)),&
                                                sum(oiii_y(:,:)),&
                                                sum(hbeta_y(:,:)))
  write(*,'(A20,ES20.4, ES20.4)')'D02',&
                                  Ztoratio( D02(sum(nii_y(:,:)),&
                                                sum(halpha_y(:,:)))),&
                                            D02(sum(nii_y(:,:)),&
                                                sum(halpha_y(:,:)))

!!$  deallocate (em_y)
!!$  deallocate (halpha_y)
!!$  deallocate (hbeta_y)
!!$  deallocate (nii_y)
!!$  deallocate (oii_y)
!!$  deallocate (oiii_y)


!=============================================================================!
!                       integration along the z-axis                          !
!=============================================================================!

  read(103,*) xmax, ymax, zmax
  ! skip first lines of the other files for projections 
  ! along x-axis
  read(1103,*)
  read(1113,*)
  read(7203,*)
  read(8203,*)
  read(8303,*)

  allocate (em_z(1:xmax, 1:ymax))
  allocate (halpha_z(1:xmax, 1:ymax))
  allocate (hbeta_z(1:xmax, 1:ymax))
  allocate (nii_z(1:xmax, 1:ymax))
  allocate (oii_z(1:xmax, 1:ymax))
  allocate (oiii_z(1:xmax, 1:ymax))
  
  do y=1,ymax
     do x=1,xmax
        read(103,*) em_z(x,y)
        read(1103,*) halpha_z(x,y), dummy1, dummy2
        read(1113,*) hbeta_z(x,y), dummy1, dummy2
        read(7203,*) nii_z(x,y), dummy1, dummy2
        read(8203,*) oii_z(x,y) , dummy1, dummy2
        read(8303,*) oiii_z(x,y), dummy1, dummy2
     end do
  end do

  open (unit=201,file=trim(argument)//".M91_z.txt")
  write(201,*) xmax, ymax, 1
  do y=1,ymax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_z(x,y) .lt. Lower_em_limit .or. oiii_z(x,y) .eq. 0._dp) then
           write(201,*) -1._dp
        else
           write(201,*) M91(nii_z(x,y), oii_z(x,y), oiii_z(x,y), hbeta_z(x,y))
        end if
     end do
  end do
  close(201)
  
  open (unit=202,file=trim(argument)//".KD02_z.txt")
  write(202,*) xmax, ymax, 1
  do y=1,ymax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_z(x,y) .lt. Lower_em_limit .or. oiii_z(x,y) .eq. 0._dp) then
           write(202,*) -1._dp
        else
           write(202,*) KD02(nii_z(x,y), oii_z(x,y), oiii_z(x,y), hbeta_z(x,y))
        end if
     end do
  end do
  close(202)

 
  open (unit=203,file=trim(argument)//".KK04_z.txt")
  write(203,*) xmax, ymax, 1
  do y=1,ymax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_z(x,y) .lt. Lower_em_limit .or. oiii_z(x,y) .eq. 0._dp) then
           write(203,*) -1._dp
        else
           write(203,*) KK04(nii_z(x,y), oii_z(x,y), oiii_z(x,y), hbeta_z(x,y))
        end if
     end do
  end do
  close(203)

  open (unit=204,file=trim(argument)//".Z94_z.txt")
  write(204,*) xmax, ymax, 1
  do y=1,ymax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_z(x,y) .lt. Lower_em_limit .or. oiii_z(x,y) .eq. 0._dp) then
           write(204,*) -1._dp
        else
           write(204,*) Z94(nii_z(x,y), oii_z(x,y), oiii_z(x,y), hbeta_z(x,y))
        end if
     end do
  end do
  close(204)

  open (unit=205,file=trim(argument)//".P05_z.txt")
  write(205,*) xmax, ymax, 1
  do y=1,ymax
     do x=1,xmax
        ! CAUTION: this routine tests for oiii_x .eq. 0 (floating pint comparison !)
        if (em_z(x,y) .lt. Lower_em_limit .or. oiii_z(x,y) .eq. 0._dp) then
           write(205,*) -1._dp
        else
           write(205,*) P05(nii_z(x,y), oii_z(x,y), oiii_z(x,y), hbeta_z(x,y))
        end if
     end do
  end do
  close(205)


  open (unit=206,file=trim(argument)//".D02_z.txt")
  write(206,*) xmax, ymax, 1
  do y=1,ymax
     do x=1,xmax
        ! In this callibrtation the [OIII] line is not considered, so the case where
        ! there is [OIII] emission needs not to be handled
        if (em_y(x,y) .lt. Lower_em_limit) then
           write(206,*) -1._dp
        else
           write(206,*) D02(nii_z(x,y), halpha_z(x,y))
        end if
     end do
  end do
  close(206)

  write(*,*) "Summary considering the entire fluxes, integrated along "//&
             "the z-axis:"
  write(*,*) "(The results should be independent of the direction,"//&
               "else there is an error !)"
!  write(*,'(A20,A20, A20)') 'Callibration', 'n(O)/n(H)', '12+log10(n(O)/n(H))'
  write(*,'(A20,ES20.4, ES20.4)') 'M91',&
                                  Ztoratio(M91(sum(nii_z(:,:)),&
                                               sum(oii_z(:,:)),&
                                               sum(oiii_z(:,:)),&
                                               sum(hbeta_z(:,:)))),&
                                            M91(sum(nii_z(:,:)),&
                                               sum(oii_z(:,:)),&
                                               sum(oiii_z(:,:)),&
                                               sum(hbeta_z(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'KD02',&
                                  Ztoratio(KD02(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))),&
                                           KD02(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'KK04',&
                                  Ztoratio(KK04(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))),&
                                           KK04(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'Z94',&
                                  Ztoratio( Z94(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))),&
                                            Z94(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))
  write(*,'(A20,ES20.4, ES20.4)') 'P05',&
                                  Ztoratio( P05(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))),&
                                            P05(sum(nii_z(:,:)),&
                                                sum(oii_z(:,:)),&
                                                sum(oiii_z(:,:)),&
                                                sum(hbeta_z(:,:)))
  write(*,'(A20,ES20.4, ES20.4)')'D02',&
                                  Ztoratio( D02(sum(nii_z(:,:)),&
                                                sum(halpha_z(:,:)))),&
                                            D02(sum(nii_z(:,:)),&
                                                sum(halpha_z(:,:)))

  deallocate (em_z)
  deallocate (halpha_z)
  deallocate (hbeta_z)
  deallocate (nii_z)
  deallocate (oii_z)
  deallocate (oiii_z)

!============================================================================!
!                                   clean-up work                            !
!============================================================================

  close(101)   ! emission measure along x-axis
  close(102)   ! emission measure along y-axis
  close(103)   ! emission measure along z-axis

  close(1101)   ! H-alpha projected along x-axis
  close(1102)   ! H_alpha projected along y-axis
  close(1103)   ! h_alpha projected along z-axis   

  close(1111)   ! H-beta projected along x-axis
  close(1112)   ! H-beta projected along y-axis
  close(1113)   ! H-beta projected along z-axis   

  close(7201)   ! "nebular" [NII] projected along x-axis
  close(7202)   ! "nebular" [NII] projected along y-axis
  close(7203)   ! "nebular" [NII] projected along z-axis   

  close(8201)   ! near-UV [OII] projected along x-axis
  close(8202)   ! near-UV [OII] projected along y-axis
  close(8203)   ! near-UV [OII] projected along z-axis   

  close(8301)   ! "nebular" [OIII] projected along x-axis
  close(8302)   ! "nebular" [OIII] projected along y-axis
  close(8303)   ! "nebular" [OIII] projected along z-axis   


end program stronglinezmaps
