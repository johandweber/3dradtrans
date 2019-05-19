program get_rec_rate
!
! Determines the entruire recombination rate of pure hydrogen
! if it doess not match the ionization rate, you are in trouble...
! 
  implicit none
  integer, parameter :: i4b=selected_int_kind(8)
  integer, parameter :: dp = selected_real_kind(13,300)

  double precision, parameter :: nc_parsec=3.0856780800D+18

  integer:: xrange, yrange, zrange
  integer:: x,y,z
  integer:: argc
  character (len=3000):: filenameprefix
  character (len=80)  :: l_cell_string
  double precision    :: l_cell
  double precision, dimension (:,:,:), allocatable :: H,e, temp
  double precision, dimension (:,:,:), allocatable :: cellrecrate
  double precision    :: dummy1, dummy2
!
! 

  argc=command_argument_count()

  if (argc .ne. 2) then
     write(*,*) 'Invalid number of arguments !'
     write(*,*)
     write(*,*) 'correct usage of the program:'
     write(*,*) 'get_rec_rate <filename_prefix> <l_cell in parsec>'
     write(*,*)
     write(*,*) ' <filename_prefix> is the part of the filename before the'
     write(*,*) 'first dot'
     stop
  end if

  call get_command_argument(1, filenameprefix)
  call get_command_argument(2, l_cell_string)
  read(l_cell_string,*) l_cell

  write(*,*) "opening ",adjustr(trim(adjustl(filenameprefix)))&
       //".H.txt"
  open(unit=101, name=adjustr(trim(adjustl(filenameprefix)))&
       //".H.txt")
  write(*,*) "opening ",adjustr(trim(adjustl(filenameprefix)))&
       //".e.txt"
  open(unit=102, name=adjustr(trim(adjustl(filenameprefix)))&
       //".e.txt")
  write(*,*) "opening ",adjustr(trim(adjustl(filenameprefix)))&
       //".temp.txt"
  open(unit=103, name=adjustr(trim(adjustl(filenameprefix)))&
       //".temp.txt")
  read(101,*) xrange, yrange, zrange
  read(102,*) xrange, yrange, zrange
  read(103,*) xrange, yrange, zrange
  print*, "xrange=", xrange
  print*, "yrange=", yrange
  print*, "ztange=", zrange
  read(*,*)
  allocate(H(xrange, yrange, zrange))
  allocate(e(xrange, yrange, zrange))
  allocate(temp(xrange, yrange, zrange))
  allocate(cellrecrate(xrange, yrange, zrange))
  do z=1, zrange
     do y=1, yrange
        do x=1,xrange
           read(101,*) dummy1,H(x,y,z), dummy2
           read(102,*) e(x,y,z)
           read(103,*) dummy1, temp(x,y,z), dummy2
           cellrecrate(x,y,z) = H(x,y,z)*e(x,y,z)*(l_cell*nc_parsec)**3*&
                HYDREC(temp(x,y,z),4)
           print*,x,y,z, cellrecrate(x,y,z)
        end do
     end do
  end do
  write(*,*)
  write(*,*) "Total number of recombination processes per second:"
  write(*,*) sum(cellrecrate)
  write(*,*) 
  close(101)
  close(102)
  close(103)

  contains

!===============================================================================
  real(dp) function HYDREC(T,N)
    !
    ! Computes total recombination rates and total recombination cooling
    ! coefficients for hydrogen by interpolating the table of
    ! D.G.Hummer, MNRAS 268, 109-112 (1994)
    ! Arguments: T is temperature
    !            N is 1 for alpha_A
    !                 2 for beta_A
    !                 3 for beta_ff
    !                 4 for alpha_b
    !                 5 for beta_B 
    ! Table entries are:
    ! Log T, Q*alpha_1, Q*alpha_B, Q*beta_1, Q*beta_B, Q*beta_ff, Q*beta_B_tot
    ! where Q=Sqrt(T) and beta_B_tot=beta_B+beta_ff
    ! M 3.05 ADI
    !
    !  called by CALC_HYDROGEN
    !
    !input
  IMPLICIT NONE
    !
    REAL   (dp )  :: T
    INTEGER(i4b)  :: N
    !local
    INTEGER(i4b)  :: I,J
    REAL   (dp )  :: Q,S,A,X
    REAL   (dp ), DIMENSION(31,7) :: C
!===============================================================================
!
    DATA ((C(I,J),J=1,7),I=1,31)/ &
         1.0_dp,1.646E-11_dp,9.283E-11_dp,1.646E-11_dp,8.287E-11_dp,1.061E-11_dp,9.348E-11_dp,&
         1.2_dp,1.646E-11_dp,8.823E-11_dp,1.646E-11_dp,7.821E-11_dp,1.068E-11_dp,8.889E-11_dp,&
         1.4_dp,1.646E-11_dp,8.361E-11_dp,1.646E-11_dp,7.356E-11_dp,1.076E-11_dp,8.432E-11_dp,&
         1.6_dp,1.646E-11_dp,7.898E-11_dp,1.646E-11_dp,6.892E-11_dp,1.085E-11_dp,7.977E-11_dp,&
         1.8_dp,1.646E-11_dp,7.435E-11_dp,1.646E-11_dp,6.430E-11_dp,1.095E-11_dp,7.525E-11_dp,&
         2.0_dp,1.646E-11_dp,6.973E-11_dp,1.645E-11_dp,5.971E-11_dp,1.106E-11_dp,7.077E-11_dp,&
         2.2_dp,1.645E-11_dp,6.512E-11_dp,1.644E-11_dp,5.515E-11_dp,1.118E-11_dp,6.633E-11_dp,&
         2.4_dp,1.645E-11_dp,6.054E-11_dp,1.643E-11_dp,5.062E-11_dp,1.132E-11_dp,6.194E-11_dp,&
         2.6_dp,1.644E-11_dp,5.599E-11_dp,1.641E-11_dp,4.614E-11_dp,1.145E-11_dp,5.758E-11_dp,&
         2.8_dp,1.642E-11_dp,5.147E-11_dp,1.638E-11_dp,4.170E-11_dp,1.161E-11_dp,5.332E-11_dp,&
         3.0_dp,1.640E-11_dp,4.700E-11_dp,1.633E-11_dp,3.734E-11_dp,1.181E-11_dp,4.915E-11_dp,&
         3.2_dp,1.636E-11_dp,4.258E-11_dp,1.625E-11_dp,3.306E-11_dp,1.202E-11_dp,4.508E-11_dp,&
         3.4_dp,1.629E-11_dp,3.823E-11_dp,1.613E-11_dp,2.888E-11_dp,1.224E-11_dp,4.112E-11_dp,&
         3.6_dp,1.620E-11_dp,3.397E-11_dp,1.594E-11_dp,2.484E-11_dp,1.248E-11_dp,3.733E-11_dp,&
         3.8_dp,1.605E-11_dp,2.983E-11_dp,1.565E-11_dp,2.098E-11_dp,1.274E-11_dp,3.373E-11_dp,&
         4.0_dp,1.582E-11_dp,2.584E-11_dp,1.522E-11_dp,1.736E-11_dp,1.303E-11_dp,3.039E-11_dp,&
         4.2_dp,1.548E-11_dp,2.204E-11_dp,1.460E-11_dp,1.402E-11_dp,1.335E-11_dp,2.737E-11_dp,&
         4.4_dp,1.499E-11_dp,1.846E-11_dp,1.374E-11_dp,1.103E-11_dp,1.369E-11_dp,2.472E-11_dp,&
         4.6_dp,1.431E-11_dp,1.520E-11_dp,1.260E-11_dp,8.442E-12_dp,1.403E-11_dp,2.247E-11_dp,&
         4.8_dp,1.341E-11_dp,1.226E-11_dp,1.119E-11_dp,6.279E-12_dp,1.434E-11_dp,2.062E-11_dp,&
         5.0_dp,1.227E-11_dp,9.696E-12_dp,9.571E-12_dp,4.539E-12_dp,1.460E-11_dp,1.914E-11_dp,&
         5.2_dp,1.093E-11_dp,7.514E-12_dp,7.844E-12_dp,3.192E-12_dp,1.478E-11_dp,1.797E-11_dp,&
         5.4_dp,9.454E-12_dp,5.710E-12_dp,6.146E-12_dp,2.185E-12_dp,1.485E-11_dp,1.704E-11_dp,&
         5.6_dp,7.920E-12_dp,4.257E-12_dp,4.601E-12_dp,1.458E-12_dp,1.482E-11_dp,1.628E-11_dp,&
         5.8_dp,6.427E-12_dp,3.117E-12_dp,3.295E-12_dp,9.484E-13_dp,1.468E-11_dp,1.563E-11_dp,&
         6.0_dp,5.058E-12_dp,2.244E-12_dp,2.262E-12_dp,6.023E-13_dp,1.444E-11_dp,1.505E-11_dp,&
         6.2_dp,3.866E-12_dp,1.590E-12_dp,1.494E-12_dp,3.738E-13_dp,1.414E-11_dp,1.451E-11_dp,&
         6.4_dp,2.877E-12_dp,1.110E-12_dp,9.520E-13_dp,2.268E-13_dp,1.380E-11_dp,1.402E-11_dp,&
         6.6_dp,2.089E-12_dp,7.642E-13_dp,5.878E-13_dp,1.348E-13_dp,1.344E-11_dp,1.358E-11_dp,&
         6.8_dp,1.485E-12_dp,5.199E-13_dp,3.528E-13_dp,7.859E-14_dp,1.311E-11_dp,1.318E-11_dp,&
         7.0_dp,1.036E-12_dp,3.498E-13_dp,2.066E-13_dp,4.499E-14_dp,1.280E-11_dp,1.285E-11_dp/
    !
    !--------------------------------------------------------------------------------------
    ! statement function for table interpolation
    X(A,I,J)=C(I,J)+A*(C(I+1,J)-C(I,J))
    !-------------------------------------------------------------------------
    !
    Q=SQRT (T)
    S=LOG10(T)
    IF(S.GT.7.0.OR.S.LT.1.0) THEN
       print*,T,S,log10(T)
       STOP 'HYDREC: Temperature out of range'
    END IF
    DO I=1,30
       IF(S.GE.C(I,1).AND.S.LE.C(I+1,1)) GOTO 1
    END DO
    STOP 'HYDREC: Unexpected error'
1   CONTINUE
    A=(S-C(I,1))/(C(I+1,1)-C(I,1))
    IF    (N.EQ.1) THEN
       HYDREC=(X(A,I,2)+X(A,I,3))/Q
    ELSEIF(N.EQ.2) THEN
       HYDREC=(X(A,I,4)+X(A,I,5))/Q
    ELSEIF(N.EQ.3) THEN
       HYDREC=X(A,I,6)/Q
    ELSEIF(N .EQ. 4) THEN
       HYDREC=X(A,I,3)/Q          
    ELSEIF(N .EQ. 5) THEN
       HYDREC=X(A,I,5)/Q
    ELSE
       STOP 'HYDREC: Bad parameter'
    ENDIF
    RETURN
  end function HYDREC


end program get_rec_rate
