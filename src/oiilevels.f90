!=============================================================================
  PROGRAM OIILEVELS
!
!   Computes the ratio of the OII emission lines in
!   ionized gas in dependence of the electron density 
!   and the gas temperature
!
    use M_data_types
    use M_natural_constants   ,only : nc_pi, nc_boltzmann, nc_planck
    use M_subsystems          ,only : generate_matrix, solve_matrix,&
                                      compare_emission_energy
    implicit none
    integer(i4b), parameter        :: n=5
    integer(i4b)                   :: linecounter
    integer(i4b)                   :: densitycounter
    integer(i4b)                   :: temperaturecounter
    real(dp), dimension (n)        :: energy =&
                                         (/    0.00_dp,&  ! 4So3/2 
                                           26810.55_dp,&  ! 2Do5/2
                                           26830.57_dp,&  ! 2Do3/2
                                           40468.01_dp,&  ! 2Po3/2
                                           40470.00_dp/)  ! 1Po1/2
!
    real(dp), dimension (n)        :: weight =&
                                         (/   4.00_dp,&  ! 4So3/2
                                              6.00_dp,&  ! 2Do5/2
                                              4.00_dp,&  ! 2Do3/2
                                              4.00_dp,&  ! 2Po3/2
                                              2.00_dp/)  ! 1Po1/2
!
    real(dp), dimension (n)        :: initial =&
                                         (/   1.00_dp,&  ! 4So3/2
                                              0.00_dp,&  ! 2Do5/2
                                              0.00_dp,&  ! 2Do3/2
                                              0.00_dp,&  ! 2Po3/2
                                              0.00_dp/)  ! 1Po1/2
    real(dp), dimension (n)        :: occupation=&
                                         (/   0.00_dp,&
                                              0.00_dp,&
                                              0.00_dp,&
                                              0.00_dp,&
                                              0.00_dp/)
         
!
    real(dp), dimension(n,n)      :: einstein, collision
    real(dp), dimension(n,n)      :: rate_matrix, emission_en
    real(dp)                      :: ne
!=============================================================================

    collision=0._dp
    einstein=0._dp
    occupation=0._dp
    rate_matrix=0._dp
    emission_en=0._dp
    
    collision(1,2) = 1.34_dp*6._dp/10._dp   !  4So3/2 <--> 2Do5/2
    collision(2,1) = collision(1,2)

    collision(1,3) = 1.34_dp*4._dp/10._dp   !  4So3/2 <--> 2Do3/2
    collision(3,1) = collision (1,3)

    collision(1,4) = 0.40_dp*4._dp/6._dp    !  4So3/2 <--> 2Po3/2
    collision(4,1) = collision(1,4)

    collision(1,5) = 0.40_dp*2._dp/6._dp    !  4So3/2 <--> 2Po1/2
    collision(5,1) = collision(1,5)

    collision(2,3) = 1.17_dp                !  2Do5/2 <--> 2Do3/2
    collision(3,2) = collision(2,3)
           
    collision(2,4) = 1.71_dp*0.428          !  2Do5/2 <--> 2Po3/2
    collision(4,2) = collision(2,3)

    collision(2,5) = 1.71_dp*0.173_dp       !  2Do5/2 <--> 2Po1/3
    collision(5,2) = collision(2,5)
         
    collision(3,4) = 1.71_dp*0.239_dp       !  2Do3/2 <--> 2Po3/2
    collision(4,3) = collision(3,4)

    collision(3,5) = 1.71_dp*0.161_dp       !  2Do3/2 <--> 2Po1/2
    collision(5,3) = collision (3,5)

    collision(4,5) = 0.29_dp                !  2Po3/2 <--> 2Po1/2


    einstein(1,2)  = 3.6e-05_dp             !  2Do5/2 --> 4So3/2
    einstein(1,3)  = 1.6e-04_dp             !  2Do3/2 --> 4So3/2 
    einstein(1,4)  = 5.7e-02_dp             !  2Po3/2 --> 4So3/2
    einstein(1,5)  = 2.3e-02_dp             !  2Po5/2 --> 4So3/2
    einstein(2,3)  = 1.3e-07_dp             !  2Do3/2 --> 2Do5/2
    einstein(2,4)  = 1.1e-01_dp             !  2Po3/2 --> 2Do5/2 
    einstein(2,5)  = 5.6e-02_dp             !  2Po1/2 --> 2Do5/2
    einstein(3,4)  = 5.8e-02_dp             !  2Po3/2 --> 2Do3/2
    einstein(3,5)  = 9.4e-02_dp             !  2Po1/2 --> 2Do3/2
    einstein(4,5)  = 1.3e-10_dp             !  2Po1/2 --> 2Po3/2

    do temperaturecounter=5000,20000,5000
       write(*,*) "# Temperature in K:" , temperaturecounter
       write (*,'("#",2A25)') "*n_e in cm⁻3","lambda3729/lambda3726"
       do densitycounter=0,1000
          ne=(1e+5**0.001)**densitycounter
          call GENERATE_MATRIX (n,einstein, collision,rate_matrix,&
               weight, energy, ne, 1._dp*temperaturecounter)
          initial =&
                                         (/   1.00_dp,&  ! 4So3/2
                                              0.00_dp,&  ! 2Do5/2
                                              0.00_dp,&  ! 2Do3/2
                                              0.00_dp,&  ! 2Po3/2
                                              0.00_dp/)  ! 1Po1/2

          call SOLVE_MATRIX    (n, rate_matrix, initial, occupation, 1e+10_dp)
          
          call COMPARE_EMISSION_ENERGY (n, occupation, einstein,&
                                     energy, emission_en)
    
          write(*,'(2ES25.5)') ne , emission_en(1,2)/emission_en(1,3)
       end do
    write(*,*)
    write(*,*)
 end do

    
END PROGRAM OIILEVELS
