PROGRAM ABSOLUTETORELATIVE 

 implicit none
  integer:: no_of_ions, argc
  integer:: xrange, yrange, zrange, xcounter, x,y,z
  double precision:: HeI_abs, HeII_abs, HeIII_abs, totalhelium
  double precision:: ion1, ion2, ion3, ion4, ion5, total
  double precision:: meandensity
  character(len=3000):: arg1,arg2

  double precision :: min=-1., max=-1., mean=0.

  argc=command_argument_count()

  call get_command_argument(1,arg1)

  read(arg1,*) no_of_ions

  read(*,*) xrange, yrange, zrange
  write(*,*) xrange, yrange, zrange

  if (no_of_ions .eq. 2) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2
              total=ion1+ion2
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(2ES18.8E3)') ion1/total, ion2/total
           end do
        end do
     end do
  end if

  
  if (no_of_ions .eq. 3) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2, ion3
              total=ion1+ion2+ion3
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(3ES18.8E3)') ion1/total, ion2/total, ion3/total
           end do
        end do
     end do
  end if

  
  if (no_of_ions .eq. 4) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2, ion3, ion4
              total=ion1+ion2+ion3+ion4
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(4ES18.8E3)') ion1/total, ion2/total, ion3/total, ion4/total
           end do
        end do
     end do
  end if

  if (no_of_ions .eq. 5) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2, ion3, ion4,ion5
              total=ion1+ion2+ion3+ion4+ion5
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              total=ion1+ion2+ion3+ion4+ion5
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(5ES18.8E3)') ion1/total, ion2/total, ion3/total, ion4/total, ion5/total
           end do
        end do
     end do
  end if

  if (no_of_ions .eq. -2) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2
              total=ion1+ion2
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(2ES18.8E3)') total
           end do
        end do
     end do
  end if

  
  if (no_of_ions .eq. -3) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2, ion3
              total=ion1+ion2+ion3
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(3ES18.8E3)') total
           end do
        end do
     end do
  end if

  
  if (no_of_ions .eq. -4) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2, ion3, ion4
              total=ion1+ion2+ion3+ion4
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(4ES18.8E3)') total
           end do
        end do
     end do
  end if

  if (no_of_ions .eq. -5) then
     do z=1, zrange
        do y=1, yrange
           do x=1, xrange
              read(*,*) ion1, ion2, ion3, ion4,ion5
              total=ion1+ion2+ion3+ion4+ion5
              if (min .lt. 0 .or. total .lt. min) min=total
              if (total .gt. max) max=total
              total=ion1+ion2+ion3+ion4+ion5
              mean=mean+total/(xrange*yrange*zrange)
              write(*,'(5ES18.8E3)') total
           end do
        end do
     end do
  end if


  write(*,*) 'min :', min
  write(*,*) 'mean :', mean
  write(*,*) 'max:', max
END PROGRAM ABSOLUTETORELATIVE
