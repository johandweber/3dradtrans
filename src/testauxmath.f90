!=============================================================================
program testauxilmath
  use M_data_types
  use M_auxmath
  implicit none
  type (polynomialinfo) :: info
  integer(i4b) :: counter
  real(dp)     :: test, start, finish
!=============================================================================
  write(*,*) "scalar arguments, ith optional info argument :"
  write(*,*) "=============================================="
  write(*,*)
  write(*,*) "solution of (x-1) = 0"
  write(*,*) solvel(1._dp,-1._dp, info)
  write(*,*) "solution of x=0"
  write(*,*) solvel(1._dp, 0._dp,info)
  write(*,*) "solution of 5*x-1-i=0"
  write(*,*) solvel (cmplx(5._dp,0._dp,kind=dp),cmplx(1._dp,1._dp,kind=dp),info)
  write(*,*) "solution of x = i" 
  write(*,*) solvel( cmplx(1._dp,0._dp,kind=dp),& 
                    -cmplx(0._dp,1._dp, kind=dp), info)
  write(*,*)
  write(*,*) "solution of (x-1) * (x-2) = 0"
  write(*,*) solveq(1._dp, -3._dp, +2._dp, info)
  write(*,*) "solution of x ** 2 = 0"
  write(*,*) solveq(1._dp, 0._dp, +0._dp, info)
  write(*,*) "solution of x ** 2 =1"
  write(*,*) solveq(1._dp, 0._dp, -1._dp, info)
  write(*,*) "solution of x**2=i"
  write(*,*) solveq(cmplx(1._dp,0._dp,kind=dp),cmplx(0._dp,0._dp,kind=dp),&
                   -cmplx(0._dp,1._dp,kind=dp), info)
  write(*,*)
! three different real roots
  write(*,*) "solution of (x-1) * (x-2) * (x-3) = 0"
  write(*,*) solvec(1._dp, -6._dp, + 11._dp, -6._dp, info)
  write(*,*)
! three different solutions method is known to work
! (computed with pencil & paper)
  write(*,*) "solution of x * (x-1) * (x-2) = 0"
  write(*,*) solvec(1._dp, -3._dp, + 2._dp, 0._dp, info)
! one real tripple root
  write(*,*)  "solution of x ** 3 = 0"
  write(*,*) solvec(1._dp, 0._dp, 0._dp, 0._dp, info)
 write(*,*) "solution of x**3 = 1"
  write(*,*) solvec(1._dp, 0._dp,0._dp,-1._dp, info)
! one double real and one single real solution
  write(*,*) "solution of x * (x-1) ** 2 = 0"
  write(*,*) solvec (1._dp, -2._dp, +1._dp, 0._dp, info)
  write(*,*) "solution of x**3 = i"
  write(*,*) solvec(cmplx(1._dp,0._dp,kind=dp), cmplx(0._dp,0._dp,kind=dp),&
                    cmplx(0._dp,0._dp,kind=dp),-cmplx(0._dp,1._dp,kind=dp),&
                    info)
  write(*,*)
  write(*,*) "solution of (x-1) * (x-2) * (x-3) * (x-4) = 0"
  write(*,*) solveb(1._dp, -10._dp, +35._dp, -50._dp, +24._dp, info)
  write(*,*) "solution of x**4 =-1"
  write(*,*) solveb(1._dp,0._dp, 0._dp, 0._dp,1._dp, info)
  write(*,*) "solution of x**4=i"
  write(*,*) solveb(cmplx(1._dp,0._dp,kind=dp),cmplx(0._dp,0._dp, kind=dp),&
                    cmplx(1._dp,0._dp,kind=dp),cmplx(0._dp,0._dp, kind=dp),&
                   -cmplx(0._dp,1._dp,kind=dp),info)
  write(*,*)
  write(*,*)
  write(*,*) "Without optional info argument :"
  write(*,*) "================================"
  write(*,*)
  write(*,*)
  write(*,*) "solution of (x-1) = 0"
  write(*,*) solvel(1._dp,-1._dp)
  write(*,*) "solution of x=0"
  write(*,*) solvel(1._dp, 0._dp)
  write(*,*)
  write(*,*) "solution of (x-1) * (x-2) = 0"
  write(*,*) solveq(1._dp, -3._dp, +2._dp)
  write(*,*) "solution of x ** 2 = 0"
  write(*,*) solveq(1._dp, 0._dp, +0._dp)
  write(*,*)
! three different real roots
  write(*,*) "solution of (x-1) * (x-2) * (x-3) = 0"
  write(*,*) solvec(1._dp, -6._dp, + 11._dp, -6._dp)
  write(*,*)
! three different solutions method is known to work
! (computed with pencil & paper)
  write(*,*) "solution of x * (x-1) * (x-2) = 0"
  write(*,*) solvec(1._dp, -3._dp, + 2._dp, 0._dp)
  write(*,*) 
! one real tripple root
  write(*,*)  "solution of x ** 3 = 0"
  write(*,*) solvec(1._dp, 0._dp, 0._dp, 0._dp)
! one double real and one single real solution
  write(*,*) "solution of x * (x-1) ** 2 = 0"
  write(*,*) solvec (1._dp, -2._dp, +1._dp, 0._dp)
  write(*,*)
  write(*,*) "solution of (x-1) * (x-2) * (x-3) * (x-4) = 0"
  write(*,*) solveb(1._dp, -10._dp, +35._dp, -50._dp, +24._dp)
  write(*,*)
  write(*,*)
!
!  Table of faculties
!
  do counter=1,170
     write(*,*) counter,"! =", facultyr(counter) 
  end do
!
  write(*,*)
  write(*,'(4A25)') "x", "exp(x)-1._dp", "expm1(x)","(exp(x)-1._dp)-expm1(x)"
  do counter=0, -30,-1
     write(*,'(4ES25.15)') 10._dp**counter, exp(10._dp**counter)-1._dp,&
          expm1(10._dp**counter),&
          exp(10._dp**counter)-1._dp-expm1(10._dp**counter) 
  end do
  
! some benchmarks conxerning the performance of the self-implemented expm(1)
! note this implemantation of expm1 is optimized for precision, not for 
! performance

  write(*,*)
  write(*,*)
  write(*,*) "EXMPM1 TEST 1:"
  write(*,*) " do counter=1,100000000"
  write(*,*) "    test=test+exp(real(counter, kind=dp)*1E-9_dp)-1._dp"
  write(*,*) " end do"
  test=0._dp
  call cpu_time(start)
  do counter=1,100000000
     test=test+exp(real(counter, kind=dp)*1E-9_dp)-1._dp
  end do
  call cpu_time(finish)
  write(*,*) "result value of benchmark: ", test
  write(*,*) finish-start, " seconds"
  write(*,*)
  write(*,*)
  write(*,*) " EXPM1 TEST 2:"
  write(*,*) " do counter=1,100000000"
  write(*,*) "    test=test+expm1(real(counter, kind=dp)*1E-9_dp)"
  write(*,*) " end do"
  test=0._dp
  call cpu_time(start)
  do counter=1,100000000
     test=test+expm1(real(counter, kind=dp)*1E-9_dp)
  end do
  call cpu_time(finish)
  write(*,*) "result value of benchmark: ", test
  write(*,*) "run time:", finish-start, " seconds"
  write(*,*)
  write(*,*)
  write(*,*) " EXPM1 TEST 3:"
  write(*,*) " do counter=1,100000000"
  write(*,*) "    test=test+expm1h(real(counter, kind=dp)*1E-9_dp)"
  write(*,*) " end do"
  test=0._dp
  call cpu_time(start)
  do counter=1,100000000
     test=test+expm1h(real(counter, kind=dp)*1E-9_dp)
  end do
  call cpu_time(finish)
  write(*,*) "result value of benchmark: ", test
  write(*,*) "run time:", finish-start, " seconds"
  write(*,*)
  write(*,*)

end program testauxilmath
