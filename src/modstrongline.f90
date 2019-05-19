!=================================Module M_R23=============================================!
!                                                                                          !
! This module provides a subset of metallicity callibrations based on stron emission lines !
! used in the paper by Kewlex&Ellison(2008)                                                !
!                                                                                          !
!==========================================================================================!

module M_strongline
!
  use M_data_types
  use M_natural_constants
!
  implicit none
!
contains
    logical function isupperbranch(nii, oii)
!   else it is the lower branch
      implicit none
      real(dp) :: nii, oii
      if (log10(nii/oii) .gt. -1.2_dp) then
         isupperbranch=.true.
      else
         isupperbranch=.false.
      end if
    end function isupperbranch
    
    
    real(dp) function Ztoratio(Z)
      implicit none
      real(dp):: Z
!      Ztoratio=Z
      Ztoratio=10._dp**(Z-12._dp)
    end function Ztoratio

    real(dp) function ratiotoZ(ratio)
      implicit none
      real(dp):: ratio
      ratiotoZ=12._dp+log10(ratio)
    end function ratiotoZ
    
    
    real(dp) function M91(nii, oii, oiii, hbeta)
      implicit none
      real(dp) :: nii, oii, oiii, hbeta
      real(dp) :: lower, upper, x, y

      x=log10((oii+oiii)/hbeta)
      y=log10(oiii/hbeta)

      lower =12._dp - 4.944_dp+0.767_dp*x+0.602_dp*x**2+&
             y*(0.29_dp+0.332_dp*x-0.331_dp*x**2)
      upper= 12._dp-2.939-0.2*x-0.237*x**2-0.305*x**3-0.0283*x**4-&
             y*(0.0047_dp-0.0221_dp*x-0.102_dp*x**2-&
             0.0817_dp*x**3-0.00717_dp*x**4)
      if (isupperbranch(nii, oii)) then
         m91=upper
      else
         m91=lower
      end if
    end function M91
  
    real(dp) function KD02(nii, oii, oiii, hbeta)
!   Solution of a quartic equation:
!   log([NII]/[OII])=kd02_fit(nii,oii)
!    However, in the
!   considered range the relation (Kewley & Ellison 2008, Eq. A3))
!   the relation 12+log(O/H) <==> log([NII]/[OII]) in monotoneous.
!   Thus, we can simply apply the method of nested intervals
      implicit none
      real(dp)     :: nii, oii, oiii, hbeta
      real(dp)     :: metallower, metalmiddle, metalupper
      integer(i4b) :: nestedcounter
!
      metalupper=10._dp
      metallower=8.4_dp
      KD02=0._dp
      if (log10(nii/oii) .lt. -1.2_dp) then
         KD02=0.5_dp*(M91(nii,oii, oiii, hbeta)+KK04(nii, oii, oiii, hbeta))
      else
         metalmiddle=(metalupper+metallower)/2._dp
         do nestedcounter=1,30
            if (kd02_fit(metalmiddle) .gt. log10(nii/oii)) then
               metalupper=metalmiddle
            else
               metallower=metalmiddle
            end if
            metalmiddle=0.5_dp*(metalupper+metallower)
!!$         write(ERROR_UNIT,'(4ES15.4)') log10(nii/oii), kd02_fit(metalmiddle),&
!!$                             metalmiddle, 10._dp**(metalmiddle-12._dp)

         end do
         KD02=metalmiddle
      end if
    end function KD02

    real(dp) function kd02_fit(Z)
!   fit polynomial (Eq. A3 of Kewley&Ellison 2008)  
      implicit none
      real(dp) :: Z
      kd02_fit  =  1106.8660_dp               &
                  - 532.15451_dp    * Z       &
                  +  96.373260_dp   * Z ** 2  &
                  -   7.8106123_dp  * Z ** 3  &
                  +   0.23928247_dp * Z ** 4
    end function kd02_fit

    real(dp) function  KK04(nii, oii,oiii, hbeta)
      implicit none
      real(dp)     :: nii, oii, oiii, hbeta
      real(dp)     :: Zlower=8.2_dp, Zupper=8.7_dp
      real(dp)     :: logqlower, logqupper
      real(dp)     :: x,y
      integer(i4b) :: iteration_steps
      y=log10(oiii/oii)
      x=log10((oii+oiii)/hbeta)
      if (log10(nii/oii) .lt. -1.2_dp) then
         do iteration_steps=1,5
            logqlower = (32.81_dp-1.153_dp*y**2+&
                 Zlower*(-3.396_dp-0.025_dp*y+0.1444_dp*y**2))/&
                 (4.603_dp-0.3119_dp*y-0.163*y**2+&
                 Zlower*(-0.48_dp+0.00271_dp*y+0.02037_dp*y**2))
            Zlower = 9.40_dp+0.777_dp*x+0.951_dp*x**2-&
                 0.072_dp*x**3-0.811_dp*x**4-&
                 logqlower*(0.272_dp+0.547_dp*x-0.513_dp*x**2)
            KK04=Zlower
         end do
      else
         do iteration_steps=1,5
            logqupper = (32.81_dp-1.153_dp*y**2+&
                 Zupper*(-3.396_dp-0.025_dp*y+0.1444_dp*y**2))/&
                 (4.603_dp-0.3119_dp*y-0.163_dp*y**2+&
                 Zupper*(-0.48_dp+0.00271*y+0.02037*y**2))
            Zupper = 9.72_dp-0.777_dp*x-0.951_dp*x**2&
                 -0.072*x**3-0.8111_dp*x**4-&
                 logqupper*&
                 (0.0737_dp-0.0713_dp*x-0.141_dp*x**2+&
                 0.0373_dp*x**3-0.058_dp*x**4)
            KK04=Zupper
         end do
      end if
    end function KK04

    real(dp) function Z94(nii,oii,oiii,hbeta)
!   only valid for upper branch
      implicit none
      real(dp):: nii,oii, oiii, hbeta
      real(dp):: x
      x=log10((oii+oiii)/hbeta)
      if (nii/oii .lt. -1.2_dp) then
         Z94=-1_dp    ! Fit not defined in this metallicity range 
      else  
         Z94 = 9.265_dp&
               -0.33_dp*x-0.202_dp*x**2&
               -0.207_dp*x**3-0.333_dp*x**4
      end if
    end function Z94

    real(dp) function P05(nii,oii,oiii,hbeta)
      implicit none
      real(dp):: nii, oii, oiii, hbeta
      real(dp):: P, R23
      R23=(oii+oiii)/hbeta
      P=oiii/R23
      if (log10(nii/oiii) .lt. -1.2) then
         P05 = (R23-106.4_dp-106.8_dp*P+3.40_dp*P**2) /&
               (17.72_dp+6.60*P+6.95*P**2-0.302*R23)
      else
         P05 =(R23+726.1_dp+842.2_dp*P+337.5_dp*P**2) /&
              (85.96_dp+82.76_dp*P+43.98_dp*P**2+1.793_dp*R23)
      end if
    end function P05

    real(dp) function D02 (nii, halpha)
      implicit none
      real(dp):: nii, halpha
      real(dp):: n2
      n2=log10(nii/halpha)
      D02=9.12_dp+0.73_dp*n2
    end function D02

end module M_strongline
