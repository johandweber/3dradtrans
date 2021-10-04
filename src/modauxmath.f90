!=============================================================================
! M_auxmath provides auxiliary mathematical functions
!=============================================================================

module M_auxmath 
  use  M_data_types
  implicit none

! So far the type polynomialinfo contains just the order of the polynomial
! as an integer value

  type polynomialinfo
     integer(i4b) :: polynomial_degree
     logical :: onlyrealcoefficients
     
  end type polynomialinfo
  
  interface solvel
     module procedure solvelc, solvelr
  end interface solvel

  interface solveq
     module procedure solveqc, solveqr
  end interface solveq

  interface solvec
     module procedure solvecc, solvecr
  end interface solvec

  interface solveb
     module procedure solvebc, solvebr
  end interface solveb

contains

!=============================================================================
  function  SOLVELC (a1, a0, info)
! solver for linear equation
! solves a1 * x + a2 = 0
    implicit none
    complex(dp)                    :: solvelc
    complex(dp)                    :: a1, a0
    type(polynomialinfo), optional :: info
!=============================================================================
    if (a1 .eq. cmplx(0._dp, 0._dp)) then
       if (a1 .eq. 0 .and. present(info)) then
          info%polynomial_degree=0
       end if
       if (a0 .eq. cmplx(0._dp, 0._dp)) then
          stop '0 * x + 0 = 0 is fullfiffed for all real/complex values of x'
       else
          stop "0 * x + C = 0 is never fulfilled for any value of x"
       end if
    end if
    if (present(info))&
         info%polynomial_degree=1
    solvelc= -a0/a1 
  end function SOLVELC


!==============================================================================
  function SOLVEQC(a2, a1, a0,info)
! solver for quadratic equation
! solves a2 * x**2 + a1 * x + a0 = 0
    implicit none
    complex(dp)                    :: a2, a1, a0
    complex, dimension (1:2)       :: solveqc
    type(polynomialinfo), optional :: info
!==============================================================================
    if (a2 .eq.  cmplx(0._dp, 0._dp)) then
       if (present(info)) then
          info%polynomial_degree=1
          solveqc(1)=solvel(a1,a0,info)
          solveqc(2)=solvel(a1,a0,info)
       else
          solveqc(1)=solvel(a1,a0)
          solveqc(2)=solvel(a1,a0)
       end if
    else
       if (present(info))&
            info%polynomial_degree=2
       solveqc(1) = (-a1 - sqrt(a1**2-4*a2*a0))/(2*a2)
       solveqc(2) = (-a1 + sqrt(a1**2-4*a2*a0))/(2*a2)
    end if
 end function SOLVEQC

!==============================================================================
  function SOLVECC(a3, a2, a1, a0,info)
! solver for cubic equation
! solves a3 * x **3 + a2 * x**2 + a1 * x + a0 = 0
    implicit none
!
    complex(dp), parameter         :: zeta=cmplx(-0.5_dp, +0.5_dp*sqrt(3._dp))
!
    complex(dp)                    :: a3, a2, a1, a0
    type(polynomialinfo), optional :: info
    complex(dp), dimension (1:3)   :: solvecc
    complex(dp), dimension (0:2)   :: subquadratic
    complex(dp), dimension (1:2)   :: quadratic
!
    complex(dp)                    :: offset, p, q, discriminant, u1, v1
!==============================================================================
!
    if (a3 .eq.  cmplx(0._dp, 0._dp)) then
       if (present(info)) then
          info%polynomial_degree=2
          quadratic=solveq(a2,a1,a0,info)
       else
          quadratic=solveq(a2,a1,a0)
       end if
       solvecc(1)=quadratic(1)
       solvecc(2)=quadratic(2)
       solvecc(3)=quadratic(1)       
    else
       offset = a2/(3._dp*a3)
!       print*, "offset:",offset
       p= (3._dp*a3*a1-a2**2)/(9._dp*a3**2)
!       print*, "p:",p
       q= (2._dp*a2**3)/(54._dp*a3**3) - (a2*a1)/(6._dp*a3**2) + a0/(2._dp*a3)
!       print*, "q:",q
       discriminant=q**2+p**3
!       print*, discriminant
       
       u1=(-q+sqrt(discriminant))**cmplx((1._dp/3._dp),0._dp,kind=dp)
!       print*, "u1:",u1
       v1=(-q-sqrt(discriminant))**cmplx((1._dp/3._dp),0._dp,kind=dp)
!       v1=-p/u1
!       print*, "v1:",v1
       solvecc(1)=u1+v1-offset
       subquadratic(2)=a3                    !  polynomial division
       subquadratic(1)=a2+subquadratic(2)*solvecc(1)
       subquadratic(0)=a1+subquadratic(1)*solvecc(1)
       if (present(info)) then
          info%polynomial_degree=2
          quadratic=solveq(subquadratic(2),subquadratic(1),subquadratic(0),info)
       else
          quadratic=solveq(subquadratic(2),subquadratic(1),subquadratic(0))
       end if
       solvecc(2)=quadratic(1)
       solvecc(3)=quadratic(2)              
    end if
 end function SOLVECC

!==============================================================================
  function SOLVEBC(a4, a3, a2, a1, a0, info)
! solver for biquadratic equation
! solves a4 * x**4 +  a3 * x **3 + a2 * x**2 + a1 * x + a0 = 0
    implicit none
!
    complex(dp)                    :: a4, a3, a2, a1, a0
    type(polynomialinfo), optional :: info
    complex, dimension (1:4)       :: solvebc
    complex, dimension (1:3)       :: cubic
!
    complex(dp)                    :: alpha, beta, gamma, P, Q, U, y, w, z
!==============================================================================
!
    if (a4 .eq.  cmplx(0._dp, 0._dp)) then
       if (present(info)) then
          info%polynomial_degree=3
          cubic=solvec(a3,a2,a1,a0,info)
       end if
       solvebc(1)=cubic(1)
       solvebc(2)=cubic(2)
       solvebc(3)=cubic(3)
       solvebc(4)=cubic(1)       
    else
       if (aimag(a4) .eq. 0._dp .and. aimag(a3) .eq. 0._dp .and. &
            aimag(a2) .eq. 0._dp .and. aimag(a1) .eq. 0._dp .and. &
            aimag(a0) .eq. 0._dp) then
          if (present(info))&
               info%onlyrealcoefficients=.true.
       else
          if (present(info))&
               info%onlyrealcoefficients=.false.
       end if
       if (present(info))&
            info%polynomial_degree=4
       alpha= - (3._dp*a3**2)/(8._dp * a4**2)+a2/a4
       beta= a3**3/(8._dp*a4**3)- (a3*a2)/(2*a4**2) +a1/a4
       gamma= -(3*a3**4)/(256*a4**4) + (a3**2*a2)/(16*a4**3) -&
               (a3*a1)/(4*a4**2) + a0/a4
       P=-alpha**2/12._dp-gamma
       Q=-alpha**3/108._dp+alpha*gamma/3._dp-beta**2/8._dp
       U=(-Q/2._dp + sqrt(Q**2/4._dp + P**3/27._dp))**(1._dp/3._dp)
       if (P .eq. 0) then
          y = -(5._dp/6._dp)*alpha-Q**(1._dp/3._dp)
       else
          y= -(5._dp/6._dp)*alpha+U-P/(3._dp*U)
       end if
       w=sqrt(alpha+2._dp*y)
       z=beta/(2._dp*w)
       solvebc(1)= -a3/(4._dp*a4)+&
                  0.5_dp*(-w-sqrt(-(alpha+2._dp*y)-2*(alpha-beta/w)))
       solvebc(2)= -a3/(4._dp*a4)+&
                  0.5_dp*(-w+sqrt(-(alpha+2._dp*y)-2*(alpha-beta/w)))
       solvebc(3)= -a3/(4._dp*a4)+&
                  0.5_dp*(+w-sqrt(-(alpha+2._dp*y)-2*(alpha+beta/w)))
       solvebc(4)= -a3/(4._dp*a4)+&
                  0.5_dp*(+w+sqrt(-(alpha+2._dp*y)-2*(alpha+beta/w)))
    end if
 end function SOLVEBC

!=============================================================================
 function SOLVELR(a1,a0,info)
! As above, but for real arguments. Return value remains complex
   implicit none
   real(dp)                        :: a1, a0
   type(polynomialinfo), optional  :: info
   complex(dp)                     :: solvelr
!=============================================================================  

   if (present(info)) then
      solvelr=solvelc(cmplx(a1,0._dp,dp), cmplx(a0, 0._dp,dp), info) 
   else
      solvelr=solvelc(cmplx(a1,0._dp,dp), cmplx(a0, 0._dp,dp))
   end if
      
 end function SOLVELR
  

!=============================================================================
 function SOLVEQR(a2,a1,a0,info)
! As above, but for real arguments. Return values remain complex
   implicit none
   real(dp)                        :: a2,a1, a0
   type(polynomialinfo), optional  :: info
   complex(dp), dimension (1:2)    :: solveqr
!=============================================================================  
!
   if (present(info)) then
      solveqr=solveqc(cmplx(a2, 0._dp,dp), cmplx(a1,0._dp,dp),&
           cmplx(a0, 0._dp,dp), info)
   else
       solveqr=solveqc(cmplx(a2, 0._dp,dp), cmplx(a1,0._dp,dp),&
                       cmplx(a0, 0._dp,dp))     
   end if
!
 end function SOLVEQR


!=============================================================================
 function SOLVECR(a3,a2, a1,a0,info)
! As above, but for real arguments. Return values remain complex
   implicit none
   real(dp)                        :: a3,a2,a1, a0
   type(polynomialinfo), optional  :: info
   complex(dp), dimension (1:3)    :: solvecr
!=============================================================================  
!
   if (present(info)) then
      solvecr=solvecc(cmplx(a3,0._dp,dp), cmplx(a2, 0._dp,dp),&
                      cmplx(a1,0._dp,dp), cmplx(a0, 0._dp,dp), info)
   else
      solvecr=solvecc(cmplx(a3,0._dp,dp), cmplx(a2, 0._dp,dp),&
                      cmplx(a1,0._dp,dp), cmplx(a0, 0._dp,dp))
   end if
!
 end function SOLVECR


!=============================================================================
 function SOLVEBR(a4,a3,a2,a1,a0,info)
! As above, but for real arguments. Return values remain complex
   implicit none
   real(dp)                        :: a4,a3,a2,a1, a0
   type(polynomialinfo), optional  :: info
   complex(dp), dimension (1:4)    :: solvebr
!=============================================================================  

   if (present(info)) then
      solvebr=solvebc(cmplx(a4, 0._dp,dp), cmplx(a3,0._dp,dp),&
                      cmplx(a2, 0._dp,dp), cmplx(a1,0._dp,dp),&
                      cmplx(a0, 0._dp,dp), info)
   else
      solvebr=solvebc(cmplx(a4, 0._dp,dp), cmplx(a3,0._dp,dp),&
                      cmplx(a2, 0._dp,dp), cmplx(a1,0._dp,dp),&
                      cmplx(a0, 0._dp,dp))
   end if
 end function SOLVEBR

!=============================================================================
 function FACULTYR(n)
! faculty of n as a real value (a real value can represent a larger number
! than an integer)
   implicit none
   real(dp)                     :: a4,a3,a2,a1, a0
   type(polynomialinfo)         :: info
   real(dp)                     :: facultyr
   integer(i4b)                 :: n           !argument of the function
   integer                      :: counter

!=============================================================================  
!
   facultyr=1._dp
!
   if (n .lt. 0) then
      stop 'error: negative arguments for the faculty function are invalid'
   else if (n .eq. 0) then 
      facultyr=1._dp
      return
   else
      do counter=1,n
         facultyr=facultyr*real(counter, KIND=dp)
      end do
   end if
 end function FACULTYR

!============================================================================
function EXPM1(x)
! The result of exp(x)-1 leads to a "small difference of large numbers"  
! problem if x is close to 0. This function fixes the problem
  implicit none
  real(dp)                     :: x, expm1
  integer(i4b)                 :: counter
!
  real(dp) , dimension (1:19), parameter :: taylorcoeff=&
  (/   1._dp/1.00000000000000_dp,&
       1._dp/2.00000000000000_dp,&
       1._dp/6.00000000000000_dp,&
       1._dp/24.0000000000000_dp,&
       1._dp/120.000000000000_dp,&
       1._dp/720.000000000000_dp,&
       1._dp/5040.00000000000_dp,&
       1._dp/40320.0000000000_dp,&
       1._dp/362880.000000000_dp,&
       1._dp/3628800.00000000_dp,&
       1._dp/39916800.0000000_dp,&
       1._dp/479001600.000000_dp,&
       1._dp/6227020800.00000_dp,&
       1._dp/87178291200.0000_dp,&
       1._dp/1307674368000.00_dp,&
       1._dp/20922789888000.0_dp,&
       1._dp/355687428096000._dp,&
       1._dp/6.402373705728000E+015_dp,&
       1._dp/1.216451004088320E+017&       
  /)
  real, dimension(1:19)         :: xpot
!============================================================================  
  expm1=0._dp
  if (abs(x) .ge. 1) then 
     expm1=exp(x)-1._dp
  else
     xpot(1)=x
     do counter=2,19
        xpot(counter)=xpot(counter-1)*x  
     end do
! We start at high powers (i.e. small summands) to avoid rounding errors.
! The upper limit of 17 habe been cosen, such that exp(1._dp)-1 = expm1(1._dp)
! for 64 bit real numbers, i.e. real(dp), as defined in M_data_types
     do counter=19,1,-1
        expm1=expm1+taylorcoeff(counter)*xpot(counter)
    end do          
  end if
end function EXPM1

!============================================================================
function EXPM1H(x)
! The result of exp(x)-1 leads to a "small difference of large numbers"  
! problem if x is close to 0. This function fixes the problem
! This is an alternative implementation that is based on the
! tangens hyperbolicus. Its quality may depend strongly
! on the implementation of the numerical implementation of the tangens 
! hyperbolicus, but turned out to be very good in the case of the
! version supplied with current (2018) versions of libc
! 
! The method has been taken from
! https://www.plunk.org/~hatch/rightway.php
! (retrieved 2018-10-18)
  implicit none
  real(dp)                     :: x, expm1h
!===========================================================================
expm1h= tanh(x/2._dp)*(exp(x)+1._dp)
end function EXPM1H

end module M_auxmath
