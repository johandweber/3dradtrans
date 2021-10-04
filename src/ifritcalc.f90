!=============================================================================
program IFRITCALC
!
! Allows calculations with IFRiT files in a non-interactive way 
!
!
  use M_data_types
!
  implicit none
  real(dp), dimension(:,:,:), allocatable :: op1, op2, result
  real(dp), dimension(1:1000)             :: readline
  real(dp)                                :: rconst, lconst
  integer(i4b)                            :: no_of_lines, linesel1, linesel2
  integer(i4b)                            :: no_of_args
!
! is there a valid operator ?
  logical                                 :: matched_option=.false.  
!
! Attention: The definition of x_max is different from the definition in
! occup_3d_ext. The reason is  that the index range in occup_3d_ext
! is -x_max...+x_max while in the ifrit coordinate_system it is
! 1... x_max (i.e. x_max(ifrit)=2*x_max(occup_3d_ext)+1
  integer(i4b)                            :: xcounter, ycounter, zcounter,&
                                             x_max, y_max, z_max,&
                                             x_old, y_old, z_old
  integer(i4b)                            :: argcounter
  real(dp)                                :: leftconst, rightconst, rightconst1,&
                                             rightconst2, rightconst3, rightconst4
  character(len=255)                      :: action
!------------------------------------------------------------------------------  
!
  call get_command_argument(1,action)
!
  if (trim(action) .eq. "add") then     
     call readinputbinaryop()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) op1(xcounter,ycounter, zcounter)+op2(xcounter,ycounter, zcounter)
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "sub") then     
     call readinputbinaryop()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) op1(xcounter,ycounter, zcounter)-op2(xcounter,ycounter, zcounter)
           end do
        end do
     end do
     matched_option=.true.
  end if

  if (trim(action) .eq. "mul") then     
     call readinputbinaryop()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) op1(xcounter,ycounter, zcounter)*op2(xcounter,ycounter, zcounter)
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "div") then     
     call readinputbinaryop()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) op1(xcounter,ycounter, zcounter)/op2(xcounter,ycounter, zcounter)
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "pow") then     
     call readinputbinaryop()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) op1(xcounter,ycounter, zcounter)**op2(xcounter,ycounter, zcounter)
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "log") then     
     call readinputfunction()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) log(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "log10") then     
     call readinputfunction()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) log10(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "exp") then     
     call readinputfunction()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) exp(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "exp10") then     
     call readinputfunction()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) 10._dp**(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "addc") then     
     call readinputrightconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) (op1(xcounter,ycounter, zcounter))+rightconst
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "subc") then     
     call readinputrightconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) (op1(xcounter,ycounter, zcounter))-rightconst
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "mulc") then     
     call readinputrightconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) (op1(xcounter,ycounter, zcounter))*rightconst
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "divc") then     
     call readinputrightconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) (op1(xcounter,ycounter, zcounter))/rightconst
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "powc") then     
     call readinputrightconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) (op1(xcounter,ycounter, zcounter))**rightconst
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "cadd") then     
     call readinputleftconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) leftconst+(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "csub") then     
     call readinputleftconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) leftconst-(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "cmul") then     
     call readinputleftconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) leftconst*(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "cdiv") then     
     call readinputleftconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) leftconst/(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "cpow") then     
     call readinputleftconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*) leftconst**(op1(xcounter,ycounter, zcounter))
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "minmax") then     
     call readinput2rightconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              write(*,*)&
                   min(&
                   max((op1(xcounter,ycounter, zcounter)), rightconst1),&
                   rightconst2)
           end do
        end do
     end do
     matched_option=.true.
  end if
!
  if (trim(action) .eq. "mask") then     
     call readinput2rightconst()
     write(*,*) x_max, y_max, z_max
     do zcounter=1,z_max
        do ycounter=1, y_max
           do xcounter=1, x_max
              if (op1(xcounter,ycounter, zcounter) .lt. rightconst1&
                   .or. op1(xcounter,ycounter, zcounter) .gt. rightconst2) then
                 write(*,*) 0
              else
                 write(*,*) 1
              end if
           end do
        end do
     end do
     matched_option=.true.
  end if

!------------------------------------------------------------------------------
contains

  subroutine readinputbinaryop ()
    implicit none
    integer(i4b)                          :: no_of_params
    integer(i4b)                          :: col1, col2
    integer(i4b)                          :: dummycounter
    character(len=255)                    :: col1name, col2name
    character(len=255)                    :: filename1, filename2
!
    call get_command_argument(2,filename1)
    call get_command_argument(3,col1name)
    read(col1name,*) col1
    call get_command_argument(4,filename2)
    call get_command_argument(5,col2name)
    read(col2name,*) col2
!
    open(unit=1001, file =trim(filename1))
    read(1001,*) x_max, y_max, z_max
    x_old=x_max; y_old=y_max; z_old=z_max
    allocate(op1(1:x_max, 1:y_max, 1:z_max))
    do zcounter=1, z_max
       do ycounter=1,y_max
          do xcounter=1, x_max
             read(1001,*)(op1(xcounter, ycounter, zcounter), dummycounter=1, col1)
          end do
       end do
    end do
    close(1001)
!
    open(unit=1002, file =trim(filename2))
    read(1002,*) x_max, y_max, z_max
    if (x_old .ne. x_max .or. y_old .ne.y_max .or.  z_old .ne. z_max) then
       stop 'The sizes of the files do not fit !'
    end if
    allocate(op2(1:x_max, 1:y_max, 1:z_max))
    do zcounter=1, z_max
       do ycounter=1,y_max
          do xcounter=1, x_max
             read(1002,*)(op2(xcounter, ycounter, zcounter), dummycounter=1, col2)
          end do
       end do
    end do
    close(1002)
!    
  end subroutine readinputbinaryop

  subroutine readinputfunction ()
    implicit none
    integer(i4b)                          :: no_of_params
    integer(i4b)                          :: col1, col2
    integer(i4b)                          :: dummycounter
    character(len=255)                    :: col1name
    character(len=255)                    :: filename1
!
    call get_command_argument(2,filename1)
    call get_command_argument(3,col1name)
    read(col1name,*) col1
!
    print*, trim(filename1)
    open(unit=1001, file =trim(filename1))
    read(1001,*) x_max, y_max, z_max
    x_old=x_max; y_old=y_max; z_old=z_max
    allocate(op1(1:x_max, 1:y_max, 1:z_max))
    do zcounter=1, z_max
       do ycounter=1,y_max
          do xcounter=1, x_max
             read(1001,*)(op1(xcounter, ycounter, zcounter),dummycounter=1,col1)
          end do
       end do
    end do
    close(1001)
!    
  end subroutine readinputfunction

  subroutine readinputrightconst ()
    implicit none
    integer(i4b)                          :: no_of_params
    integer(i4b)                          :: col1, col2
    integer(i4b)                          :: dummycounter
    character(len=255)                    :: col1name, rightconstname
    character (len=255)                   :: rightconstname1,&
                                             rightconstname2,&
                                             rightconstname3,&
                                             rightconstname4
    character(len=255)                    :: filename1
!
    call get_command_argument(2,filename1)
    call get_command_argument(3,col1name)
    read(col1name,*) col1
    call get_command_argument(4, rightconstname)
    read(rightconstname,*) rightconst                  
!
    print*, trim(filename1)
    open(unit=1001, file =trim(filename1))
    read(1001,*) x_max, y_max, z_max
    x_old=x_max; y_old=y_max; z_old=z_max
    allocate(op1(1:x_max, 1:y_max, 1:z_max))
    do zcounter=1, z_max
       do ycounter=1,y_max
          do xcounter=1, x_max
             read(1001,*)(op1(xcounter, ycounter, zcounter),dummycounter=1,col1)
          end do
       end do
    end do
    close(1001)
!    
  end subroutine readinputrightconst
!
  subroutine readinputleftconst ()
    implicit none
    integer(i4b)                          :: no_of_params
    integer(i4b)                          :: col1, col2
    integer(i4b)                          :: dummycounter
    character(len=255)                    :: col1name, leftconstname
    character(len=255)                    :: filename1
!
    call get_command_argument(2, leftconstname)
    read(leftconstname,*) leftconst                  
    call get_command_argument(3,filename1)
    call get_command_argument(4,col1name)
    read(col1name,*) col1
!
    print*, trim(filename1)
    open(unit=1001, file =trim(filename1))
    read(1001,*) x_max, y_max, z_max
    x_old=x_max; y_old=y_max; z_old=z_max
    allocate(op1(1:x_max, 1:y_max, 1:z_max))
    do zcounter=1, z_max
       do ycounter=1,y_max
          do xcounter=1, x_max
             read(1001,*)(op1(xcounter, ycounter, zcounter),dummycounter=1,col1)
          end do
       end do
    end do
    close(1001)
!    
  end subroutine readinputleftconst

  subroutine readinput2rightconst ()
    implicit none
    integer(i4b)                          :: no_of_params
    integer(i4b)                          :: col1, col2
    integer(i4b)                          :: dummycounter
    character(len=255)                    :: col1name, rightconstname1,&
                                             rightconstname2
    character(len=255)                    :: filename1
!
    call get_command_argument(2,filename1)
    call get_command_argument(3,col1name)
    read(col1name,*) col1
    call get_command_argument(4, rightconstname1)
    read(rightconstname1,*) rightconst1                  
    call get_command_argument(5, rightconstname2)
    read(rightconstname2,*) rightconst2                  
    open(unit=1001, file =trim(filename1))
    read(1001,*) x_max, y_max, z_max
    x_old=x_max; y_old=y_max; z_old=z_max
    allocate(op1(1:x_max, 1:y_max, 1:z_max))
    do zcounter=1, z_max
       do ycounter=1,y_max
          do xcounter=1, x_max
             read(1001,*)(op1(xcounter, ycounter, zcounter),dummycounter=1,col1)
          end do
       end do
    end do
    close(1001)
  end subroutine readinput2rightconst


end program IFRITCALC
!=============================================================================
