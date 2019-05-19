!=============================================================================!
!
!" PROGRAM R23
!
! This program will use the output of occup_3d_ext.
! Then it will use the simulated image to retrieve the
! metallicity again.
! This is done to study problemns both withn the methodn and with the 3d code.
!=============================================================================!

PROGRAM R23
  use iso_fortran_env
  use M_data_types
  use M_natural_constants
  use M_strongline
  implicit none
  integer(i4b) :: counter
!-----------------------------------------------------------------------------

  write(*,*) "Lexington HII-40 test"
  write(*,*) "Input O/H ratio:"
  write(*,'(ES15.3)') (33._dp*1e-5_dp)
  write(*,*)

  write(*,'(7A15)') '', 'M91', 'KD02', 'KK04', 'Z94', 'P05','D02'
  write(*,'(A15, 6ES15.4)') "Med",&
                            Ztoratio(M91(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(KD02(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(KK04(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(Z94(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(P05(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(D02(0.725_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "GF",&
                            Ztoratio(M91(0.669_dp,1.94_dp,2.21_dp,1.0_dp)),&
                            Ztoratio(KD02(0.669_dp,1.94_dp,2.21_dp,1.0_dp)),&
                            Ztoratio(KK04(0.669_dp,1.94_dp,2.21_dp,1.0_dp)),&
                            Ztoratio(Z94(0.669_dp,1.94_dp,2.21_dp,1.0_dp)),&
                            Ztoratio(P05(0.669_dp,1.94_dp,2.21_dp,1.0_dp)),&
                            Ztoratio(D02(0.669_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "HN",&
                            Ztoratio(M91(0.817_dp,1.94_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(KD02(0.817_dp,1.94_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(KK04(0.817_dp,1.94_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(Z94(0.817_dp,1.94_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(P05(0.817_dp,1.94_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(D02(0.817_dp,2.8_dp))


 write(*,'(A15, 6ES15.4)')  "DP",&
                            Ztoratio(M91(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(KD02(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(KK04(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(Z94(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(P05(0.725_dp,2.12_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(D02(0.725_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "TK",&
                            Ztoratio(M91(0.69_dp,1.6_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(KD02(0.69_dp,1.6_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(KK04(0.69_dp,1.6_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(Z94(0.69_dp,1.6_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(P05(0.69_dp,1.6_dp,2.20_dp,1.0_dp)),&
                            Ztoratio(D02(0.69_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "PH",&
                            Ztoratio(M91(0.736_dp,2.19_dp,1.93_dp,1.0_dp)),&
                            Ztoratio(KD02(0.736_dp,2.19_dp,1.93_dp,1.0_dp)),&
                            Ztoratio(KK04(0.736_dp,2.19_dp,1.93_dp,1.0_dp)),&
                            Ztoratio(Z94(0.736_dp,2.19_dp,1.93_dp,1.0_dp)),&
                            Ztoratio(P05(0.736_dp,2.19_dp,1.93_dp,1.0_dp)),&
                            Ztoratio(D02(0.736_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "RS",&
                            Ztoratio(M91(0.723_dp,1.88_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(KD02(0.723_dp,1.88_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(KK04(0.723_dp,1.88_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(Z94(0.723_dp,1.88_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(P05(0.723_dp,1.88_dp,2.17_dp,1.0_dp)),&
                            Ztoratio(D02(0.723_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "RR",&
                            Ztoratio(M91(0.807_dp,2.26_dp,2.08_dp,1.0_dp)),&
                            Ztoratio(KD02(0.807_dp,2.26_dp,2.08_dp,1.0_dp)),&
                            Ztoratio(KK04(0.807_dp,2.26_dp,2.08_dp,1.0_dp)),&
                            Ztoratio(Z94(0.807_dp,2.26_dp,2.08_dp,1.0_dp)),&
                            Ztoratio(P05(0.807_dp,2.26_dp,2.08_dp,1.0_dp)),&
                            Ztoratio(D02(0.807_dp,2.8_dp))


 write(*,'(A15, 6ES15.4)')  "BE",&
                            Ztoratio(M91(0.56_dp,2.03_dp,2.13_dp,1.0_dp)),&
                            Ztoratio(KD02(0.56_dp,2.03_dp,2.13_dp,1.0_dp)),&
                            Ztoratio(KK04(0.56_dp,2.03_dp,2.13_dp,1.0_dp)),&
                            Ztoratio(Z94(0.56_dp,2.03_dp,2.13_dp,1.0_dp)),&
                            Ztoratio(P05(0.56_dp,2.03_dp,2.13_dp,1.0_dp)),&
                            Ztoratio(D02(0.56_dp,2.8_dp))


  write(*,*)
  write(*,*) "Lexington HII-20 test"
  write(*,*) "Input O/H ratio:"
  write(*,'(ES15.4)') (33._dp*1e-5_dp)
  write(*,*)

  write(*,'(7A15)') '', 'M91', 'KD02', 'KK04', 'Z94', 'P05','D02'
  write(*,'(A15, 6ES15.4)') "Med",&
                            Ztoratio(M91(0.803_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(KD02(0.803_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(KK04(0.803_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(Z94(0.803_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(P05(0.803_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(D02(0.803_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "GF",&
                            Ztoratio(M91(0.745_dp,1.01_dp,0.0021_dp,1.0_dp)),&
                            Ztoratio(KD02(0.745_dp,1.01_dp,0.0021_dp,1.0_dp)),&
                            Ztoratio(KK04(0.745_dp,1.01_dp,0.0012_dp,1.0_dp)),&
                            Ztoratio(Z94(0.745_dp,1.01_dp,0.0021_dp,1.0_dp)),&
                            Ztoratio(P05(0.745_dp,1.01_dp,0.0021_dp,1.0_dp)),&
                            Ztoratio(D02(0.745_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "HN",&
                            Ztoratio(M91(0.786_dp,1.13_dp,0.0016_dp,1.0_dp)),&
                            Ztoratio(KD02(0.786_dp,1.13_dp,0.0016_dp,1.0_dp)),&
                            Ztoratio(KK04(0.786_dp,1.13_dp,0.0016_dp,1.0_dp)),&
                            Ztoratio(Z94(0.786_dp,1.13_dp,0.0016_dp,1.0_dp)),&
                            Ztoratio(P05(0.786_dp,1.13_dp,0.0016_dp,1.0_dp)),&
                            Ztoratio(D02(0.786_dp,2.8_dp))


 write(*,'(A15, 6ES15.4)')  "DP",&
                            Ztoratio(M91(0.785_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(KD02(0.785_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(KK04(0.785_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(Z94(0.785_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(P05(0.785_dp,1.10_dp,0.0015_dp,1.0_dp)),&
                            Ztoratio(D02(0.785_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "TK",&
                            Ztoratio(M91(0.925_dp,1.04_dp,0.0024_dp,1.0_dp)),&
                            Ztoratio(KD02(0.925_dp,1.04_dp,0.0024_dp,1.0_dp)),&
                            Ztoratio(KK04(0.925_dp,1.04_dp,0.0024_dp,1.0_dp)),&
                            Ztoratio(Z94(0.925_dp,1.04_dp,0.0024_dp,1.0_dp)),&
                            Ztoratio(P05(0.925_dp,1.04_dp,0.0024_dp,1.0_dp)),&
                            Ztoratio(D02(0.925_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "PH",&
                            Ztoratio(M91(0.843_dp,1.22_dp,0.0014_dp,1.0_dp)),&
                            Ztoratio(KD02(0.843_dp,1.22_dp,0.0014_dp,1.0_dp)),&
                            Ztoratio(KK04(0.843_dp,1.22_dp,0.0014_dp,1.0_dp)),&
                            Ztoratio(Z94(0.843_dp,1.22_dp,0.0014_dp,1.0_dp)),&
                            Ztoratio(P05(0.843_dp,1.22_dp,0.0014_dp,1.0_dp)),&
                            Ztoratio(D02(0.843_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "RS",&
                            Ztoratio(M91(0.803_dp,1.08_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(KD02(0.803_dp,1.08_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(KK04(0.803_dp,1.08_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(Z94(0.803_dp,1.08_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(P05(0.803_dp,1.08_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(D02(0.803_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "RR",&
                            Ztoratio(M91(0.915_dp,1.17_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(KD02(0.915_dp,1.17_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(KK04(0.915_dp,1.17_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(Z94(0.915_dp,1.17_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(P05(0.915_dp,1.17_dp,0.0010_dp,1.0_dp)),&
                            Ztoratio(D02(0.915_dp,2.8_dp))



  write(*,*)
  write(*,*) "Lexington PN-75 test"
  write(*,*) "Input O/H ratio:"
  write(*,'(ES15.4)') (30._dp*1e-5_dp)
  write(*,*)

  write(*,'(7A15)')  '', 'M91', 'KD02', 'KK04', 'Z94', 'P05','D02'
  write(*,'(A15, 6ES15.4)') "Med",&
                            Ztoratio(M91(0.097_dp,0.262_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(KD02(0.097_dp,0.262_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(KK04(0.097_dp,0.262_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(Z94(0.097_dp,0.262_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(P05(0.097_dp,0.262_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(D02(0.097_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "GF",&
                            Ztoratio(M91(0.069_dp,0.178_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(KD02(0.069_dp,0.178_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(KK04(0.069_dp,0.178_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(Z94(0.069_dp,0.178_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(P05(0.069_dp,0.178_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(D02(0.069_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "HN",&
                            Ztoratio(M91(0.097_dp,0.262_dp,13.2_dp,1.0_dp)),&
                            Ztoratio(KD02(0.097_dp,0.262_dp,13.2_dp,1.0_dp)),&
                            Ztoratio(KK04(0.097_dp,0.262_dp,13.2_dp,1.0_dp)),&
                            Ztoratio(Z94(0.097_dp,0.262_dp,13.2_dp,1.0_dp)),&
                            Ztoratio(P05(0.097_dp,0.262_dp,13.2_dp,1.0_dp)),&
                            Ztoratio(D02(0.097_dp,2.8_dp))


 write(*,'(A15, 6ES15.4)')  "DP",&
                            Ztoratio(M91(0.089_dp,0.266_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(KD02(0.089_dp,0.266_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(KK04(0.089_dp,0.266_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(Z94(0.089_dp,0.266_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(P05(0.089_dp,0.266_dp,11.7_dp,1.0_dp)),&
                            Ztoratio(D02(0.089_dp,2.8_dp))


 write(*,'(A15, 6ES15.4)')  "PH",&
                            Ztoratio(M91(0.108_dp,0.262_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(KD02(0.108_dp,0.262_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(KK04(0.108_dp,0.262_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(Z94(0.108_dp,0.262_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(P05(0.108_dp,0.262_dp,10.1_dp,1.0_dp)),&
                            Ztoratio(D02(0.108_dp,2.8_dp))

 write(*,'(A15, 6ES15.4)')  "RS",&
                            Ztoratio(M91(0.119_dp,0.311_dp,11.8_dp,1.0_dp)),&
                            Ztoratio(KD02(0.119_dp,0.311_dp,11.8_dp,1.0_dp)),&
                            Ztoratio(KK04(0.119_dp,0.311_dp,11.8_dp,1.0_dp)),&
                            Ztoratio(Z94(0.119_dp,0.311_dp,11.8_dp,1.0_dp)),&
                            Ztoratio(P05(0.119_dp,0.311_dp,11.8_dp,1.0_dp)),&
                            Ztoratio(D02(0.119_dp,2.8_dp))
!
END PROGRAM R23
