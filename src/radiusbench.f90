PROGRAM radiusbench

  implicit none

  integer:: xrange, yrange, zrange
  integer:: midx, midy, midz
  integer:: x, y, z
  integer:: xsource, ysource, zsource
  integer:: argument_counter
  character(len=3000) :: infile, outfile, distancename,datacolsname, xs, ys, zs  
  double precision:: distance
  double precision, dimension (1:5) ::data
  integer:: datacols

  argument_counter=command_argument_count()

  if (argument_counter .ne. 7) then
     print*, "Incorrect number of command line arguments !"
     print*, "correct usage:"
     print*, "radiusbench <cellsize> <data_cols> <x_source> <y_source> <z_source> <infile> <outfile>"
     stop
  end if
  call get_command_argument(1, distancename)
  call get_command_argument(2, datacolsname)
  call get_command_argument(3, xs)
  call get_command_argument(4, ys)
  call get_command_argument(5, zs)
  call get_command_argument(6, infile)
  call get_command_argument(7, outfile)
  
  
  open(unit=1, file=trim(infile))
  open(unit=2, file=trim(outfile))
  read(distancename,*) distance
  read(datacolsname,*) datacols
  read(xs,*) xsource
  read(ys,*) ysource
  read(zs,*) zsource


  read(1,*) xrange, yrange, zrange
  midx=(xrange+1)/2
  midy=(yrange+1)/2
  midz=(zrange+1)/2


  select case (datacols)
     case (0)
        do z=1, zrange
           do y=1, yrange
              do x=1, xrange
                 write(2,'(ES20.10)') get_distance(distance,&
                                                   (x-midx-xsource),&
                                                   (y-midy-ysource),&
                                                   (z-midz-zsource))
              end do
           end do
        end do
     case (1)
        do z=1, zrange
           do y=1, yrange
              do x=1, xrange
                 read(1,*) data(1)
                 write(2,'(2ES20.10)') get_distance(distance,&
                                                    (x-midx-xsource),&
                                                    (y-midy-ysource),&
                                                    (z-midz-zsource)) ,&
                                                    data(1)
              end do
           end do
        end do
        case(2)
           do z=1, zrange
              do y=1, yrange
                 do x=1, xrange
                    read(1,*) data(1), data(2)
                    write(2,'(2ES20.10)') get_distance(distance,&
                         (x-midx-xsource),&
                         (y-midy-ysource),&
                         (z-midz-zsource)),&
                         data(1),data(2)
                 end do
              end do
           end do
           case(3)
              do z=1,zrange
                 do y=1,yrange
                    do x=1,xrange
                   read(1,*) data(1), data(2),data(3)
                    write(2,'(4ES20.10)') get_distance(distance,&
                         (x-midx-xsource),&
                         (y-midy-ysource),&
                         (z-midz-zsource)),&
                         data(1),data(2),data(3)
                    end do
                 end do
              end do
           case (4)
              do z=1,zrange
                 do y=1, yrange
                    do x=1, xrange
                   read(1,*) data(1), data(2),data(3),data(4)
                    write(2,'(5ES20.10)') get_distance(distance,&
                         (x-midx-xsource),&
                         (y-midy-ysource),& 
                         (z-midz-zsource)),&
                         data(1),data(2),data(3),&
                         data(4)
                    end do
                 end do
              end do
           case(5)
              do z=1,zrange
                 do y=1, yrange
                    do x=1, xrange
                       read(1,*) data(1), data(2),data(3),data(4), data(5)
                       write(2,'(6ES20.10)') get_distance(distance,&
                            (x-midx-xsource),&
                            (y-midy-ysource),&
                            (z-midz-zsource)),&
                            data(1),data(2),data(3),&
                            data(4), data(5)
                    end do
                 end do
              end do
           end select
  close(1)
  close(2)
  
  stop

  contains
    function get_distance(cellsize, deltax, deltay, deltaz)
      implicit none
      double precision:: get_distance
      double precision:: cellsize
      integer :: deltax, deltay, deltaz
      get_distance=sqrt((cellsize*deltax)**2+(cellsize*deltay)**2+&
                        (cellsize*deltaz)**2)
    end function get_distance

END PROGRAM radiusbench
