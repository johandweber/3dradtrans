!===============================================================================
!                          ****MODULE M_min_healpix****
! This routine was taken from the Healpix project ((c) Gorski et. al.) and is 
! licensed by GPL  It gives the components of an unit vector according to the 
! Healpix Scheme 
!===============================================================================

module M_min_healpix
!--- routines: 
!    MK_PIX2XY
!    PIX2VEC_NEST
!        
    use M_data_types
    use M_natural_constants         , only: nc_pi 
    implicit none
!
    integer(i4b), parameter                        :: LGT=4    
!   2^13 : largest nside available
    INTEGER(i4b), private, parameter               :: ns_max=8192 

!   initialise arrays pix2x and pix2y used in several routines
    integer(i4b), private, save, dimension(0:1023) :: pix2x=0, pix2y=0

    real(dp), parameter:: halfpi=nc_pi/2._dp
    
    contains
     
!================================================================================ 
  subroutine MK_PIX2XY()
!
!     constructs the array giving x and y in the face from pixel number
!     for the nested (quad-cube like) ordering of pixels
!
!     the bits corresponding to x and y are interleaved in the pixel number
!     one breaks up the pixel number by even and odd bits
!
! called by PIX2VEC_NEST
!
  implicit none
!    
    integer(i4b) ::  kpix, jpix, ix, iy, ip, id
!================================================================================
!
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31
    enddo
!
    return
  end subroutine MK_PIX2XY
  


!================================================================================
  subroutine PIX2VEC_NEST(nside, ipix, vector, vertex)
!
!    renders vector (x,y,z) coordinates of the nominal pixel center
!    for the pixel number ipix (NESTED scheme)
!    given the map resolution parameter nside
!    also returns the (x,y,z) position of the 4 pixel vertices (=corners)
!    in the order N,W,S,E
!
! called by INIT_GEOMETRY
!
  implicit none
!
    integer(i4b), INTENT(IN) :: nside, ipix
    real(dp), INTENT(OUT), dimension(1:) :: vector
    real(dp),     INTENT(OUT),dimension(1:,1:), optional :: vertex

    integer(i4b) :: npix, npface, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4
    real(dp) :: z, fn, fact1, fact2, sth, phi

    integer(i4b) ::  ix, iy, face_num
    ! common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    ! coordinate of the lowest corner of each face
    INTEGER(i4b), dimension(1:12) :: jrll =&
         (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(i4b), dimension(1:12) :: jpll =&
         (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

    real(dp) :: z_nv, z_sv
    logical(lgt) :: do_vertex
!================================================================================
!
    if (nside<1 .or. nside>ns_max) stop "nside out of range"
    npix = 12 * nside**2
    if (ipix <0 .or. ipix>npix-1) stop "ipix out of range"
!
! initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    fn = real(nside,kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)
    nl4   = 4*nside
!
    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          stop " pix2vec_ring : vertex array has wrong size "
       endif
    endif
!
! finds the face, and the number in the face
    npface = nside**2
!
    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}
!
! finds the x,y on the face (starting from the lowest corner)
! from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits
!
    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
!
! transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}
!
!     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}
!
    nr = nside                  ! equatorial region (the most frequent)
    z  = (2*nside-jr)*fact2
    kshift = MODULO(jr - nside, 2)
    if (do_vertex) then
       z_nv = (2*nside-jr+1)*fact2
       z_sv = (2*nside-jr-1)*fact2
       if (jr == nside) then ! northern transition
          z_nv =  1.0_dp - (nside-1)**2 * fact1
       elseif (jr == 3*nside) then  ! southern transition
          z_sv = -1.0_dp + (nside-1)**2 * fact1
       endif
    endif
    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = 1.0_dp - (nr-1)**2*fact1
          z_sv = 1.0_dp - (nr+1)**2*fact1
       endif
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_dp + nr*nr*fact1
       kshift = 0
       if (do_vertex) then
          z_nv = - 1.0_dp + (nr+1)**2*fact1
          z_sv = - 1.0_dp + (nr-1)**2*fact1
       endif
    endif
!
! computes the phi coordinate on the sphere, in [0,2Pi]
! 'phi' number in the ring in {1,4*nr}
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2 
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4
!
    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)
    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z
!
    return
  end subroutine PIX2VEC_NEST
  
  
  !=======================================================================

end module M_min_healpix
!===============================================================================
!                       ****END MODULE M_min_healpix****
!===============================================================================
