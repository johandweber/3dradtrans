module M_chemical_data
  use M_data_types
  
  integer, parameter :: kel =30                           ! number of considered elements
  
  real(dp)           :: Hetosolar = 1._dp,&
                        Litosolar = 1._dp,&
                        Betosolar = 1._dp,&
                        Btosolar  = 1._dp,&
                        Ctosolar  = 1._dp,&
                        Ntosolar  = 1._dp,&
                        Otosolar  = 1._dp,&
                        Ftosolar  = 1._dp,&
                        Netosolar = 1._dp,&
                        Natosolar = 1._dp,&
                        Mgtosolar = 1._dp,&
                        Altosolar = 1._dp,&
                        Sitosolar = 1._dp,&
                        Ptosolar  = 1._dp,&
                        Stosolar  = 1._dp,&
                        Cltosolar = 1._dp,&
                        Artosolar = 1._dp,&
                        Ktosolar  = 1._dp,&
                        Catosolar = 1._dp,&
                        Titosolar = 1._dp,&
                        Vtosolar  = 1._dp,&
                        Crtosolar = 1._dp,&
                        Mntosolar = 1._dp,&
                        Fetosolar = 1._dp,&
                        Cotosolar = 1._dp,&
                        Nitosolar = 1._dp,&
                        Cutosolar = 1._dp,&
                        Zntosolar = 1._dp,&
                        Gatosolar = 1._dp,&
                        Getosolar = 1._dp,&
                        Astosolar = 1._dp,&
                        Setosolar = 1._dp,&
                        Brtosolar = 1._dp,&
                        Krtosolar = 1._dp,&
                        Rbtosolar = 1._dp,&
                        Srtosolar = 1._dp,&
                        Ytosolar  = 1._dp,&
                        Zrtosolar = 1._dp,&
                        Nbtosolar = 1._dp,&
                        Motosolar = 1._dp,&
                        Tctosolar = 1._dp,&
                        Rutosolar = 1._dp,&
                        Rhtosolar = 1._dp,&
                        Pdtosolar = 1._dp,&
                        Agtosolar = 1._dp,&
                        Cdtosolar = 1._dp,&
                        Intosolar = 1._dp,&
                        Sntosolar = 1._dp,&
                        Sbtosolar = 1._dp,&
                        Tetosolar = 1._dp,&
                        Itosolar  = 1._dp,&
                        Xetosolar = 1._dp,&
                        Cstosolar = 1._dp,&
                        Batosolar = 1._dp,&
                        Latosolar = 1._dp,&
                        Cetosolar = 1._dp,&
                        Prtosolar = 1._dp,&
                        Ndtosolar = 1._dp,&
                        Pmtosolar = 1._dp,&
                        Smtosolar = 1._dp,&
                        Eutosolar = 1._dp,&
                        Gdtosolar = 1._dp,&
                        Tbtosolar = 1._dp,&
                        Dytosolar = 1._dp,&
                        Hotosolar = 1._dp,&
                        Ertosolar = 1._dp,&
                        Tmtosolar = 1._dp,&
                        Ybtosolar = 1._dp,&
                        Lutosolar = 1._dp,&
                        Hftosolar = 1._dp,&
                        Tatosolar = 1._dp,&
                        Wtosolar  = 1._dp,&
                        Retosolar = 1._dp,&
                        Ostosolar = 1._dp,&
                        Irtosolar = 1._dp,&
                        Pttosolar = 1._dp,&
                        Autosolar = 1._dp,&
                        Hgtosolar = 1._dp,&
                        Tltosolar = 1._dp,&
                        Pbtosolar = 1._dp,&
                        Bitosolar = 1._dp,&
                        Potosolar = 1._dp,&
                        Attosolar = 1._dp,&
                        Rntosolar = 1._dp,&
                        Frtosolar = 1._dp,&
                        Ratosolar = 1._dp,&
                        Actosolar = 1._dp,&
                        Thtosolar = 1._dp,&
                        Patosolar = 1._dp,&
                        Utosolar  = 1._dp,&
                        Nptosolar = 1._dp,&
                        Putosolar = 1._dp,&
                        Amtosolar = 1._dp,&
                        Cmtosolar = 1._dp

  logical, dimension(1:96)  :: is_considered = (/&
                       .TRUE. , & ! H 
                       .TRUE. , & ! He 
                       .FALSE., & ! Li 
                       .FALSE., & ! Be 
                       .FALSE., & ! B
                       .TRUE. , & ! C
                       .TRUE. , & ! N
                       .TRUE. , & ! O
                       .FALSE., & ! F
                       .TRUE. , & ! Ne
                       .FALSE., & ! Na 
                       .FALSE., & ! Mg
                       .FALSE., & ! Al
                       .FALSE., & ! Si
                       .FALSE., & ! P
                       .TRUE. , & ! S
                       .FALSE. ,& ! Cl 
                       .FALSE. ,& ! Ar
                       .FALSE. ,& ! K
                       .FALSE. ,& ! Ca
                       .FALSE. ,& ! Sc
                       .FALSE. ,& ! Ti 
                       .FALSE. ,& ! V
                       .FALSE. ,& ! Cr
                       .FALSE. ,& ! Mn
                       .FALSE. ,& ! Fe
                       .FALSE. ,& ! Co 
                       .FALSE. ,& ! Ni
                       .FALSE. ,& ! Cu
                       .FALSE. ,& ! Zn
                       .FALSE. ,& ! Ga
                       .FALSE. ,& ! Ge 
                       .FALSE. ,& ! As
                       .FALSE. ,& ! Se
                       .FALSE. ,& ! Br
                       .FALSE. ,& ! Kr
                       .FALSE. ,& ! Rb
                       .FALSE. ,& ! Sr
                       .FALSE. ,& ! Y
                       .FALSE. ,& ! Zr
                       .FALSE. ,& ! Nb
                       .FALSE. ,& ! Mo
                       .FALSE. ,& ! Tc
                       .FALSE. ,& ! Ru
                       .FALSE. ,& ! Rh
                       .FALSE. ,& ! Pd
                       .FALSE. ,& ! Ag
                       .FALSE. ,& ! Cd
                       .FALSE. ,& ! In
                       .FALSE. ,& ! Sn
                       .FALSE. ,& ! Sb
                       .FALSE. ,& ! Te
                       .FALSE. ,& ! I
                       .FALSE. ,& ! Xe 
                       .FALSE. ,& ! Cs
                       .FALSE. ,& ! Ba
                       .FALSE. ,& ! La
                       .FALSE. ,& ! Ce
                       .FALSE. ,& ! Pr 
                       .FALSE. ,& ! Nd
                       .FALSE. ,& ! Pm
                       .FALSE. ,& ! Sm
                       .FALSE. ,& ! Eu
                       .FALSE. ,& ! Gd 
                       .FALSE. ,& ! Tb
                       .FALSE. ,& ! Dy
                       .FALSE. ,& ! Ho
                       .FALSE. ,& ! Er
                       .FALSE. ,& ! Tm 
                       .FALSE. ,& ! Yb
                       .FALSE. ,& ! Lu
                       .FALSE. ,& ! Hf
                       .FALSE. ,& ! Ta
                       .FALSE. ,& ! W 
                       .FALSE. ,& ! Re
                       .FALSE. ,& ! Os
                       .FALSE. ,& ! Ir
                       .FALSE. ,& ! Pt
                       .FALSE. ,& ! Au 
                       .FALSE. ,& ! Hg
                       .FALSE. ,& ! Tl
                       .FALSE. ,& ! Pb
                       .FALSE. ,& ! Bi
                       .FALSE. ,& ! Po 
                       .FALSE. ,& ! At
                       .FALSE. ,& ! Rn
                       .FALSE. ,& ! Fr
                       .FALSE. ,& ! Ra
                       .FALSE. ,& ! Ac 
                       .FALSE. ,& ! Th
                       .FALSE. ,& ! Pa
                       .FALSE. ,& ! U
                       .FALSE. ,& ! Np
                       .FALSE. ,& ! Pu 
                       .FALSE. ,& ! Am
                       .FALSE.  & ! Cm
                       /)

  integer, dimension(1:96)  :: lower_ion= (/&
                       1 ,& ! H 
                       1 ,& ! He 
                       0 ,& ! Li 
                       0 ,& ! Be 
                       0 ,& ! B
                       1 ,& ! C
                       1 ,& ! N
                       1 ,& ! O
                       0 ,& ! F
                       1 ,& ! Ne
                       0 ,& ! Na 
                       0 ,& ! Mg
                       0 ,& ! Al
                       0 ,& ! Si
                       0 ,& ! P
                       1 ,& ! S
                       0 ,& ! Cl 
                       0 ,& ! Ar
                       0 ,& ! K
                       0 ,& ! Ca
                       0 ,& ! Sc
                       0 ,& ! Ti 
                       0 ,& ! V
                       0 ,& ! Cr
                       0 ,& ! Mn
                       0 ,& ! Fe
                       0 ,& ! Co 
                       0 ,& ! Ni
                       0 ,& ! Cu
                       0 ,& ! Zn
                       0 ,& ! Ga
                       0 ,& ! Ge 
                       0 ,& ! As
                       0 ,& ! Se
                       0 ,& ! Br
                       0 ,& ! Kr
                       0 ,& ! Rb
                       0 ,& ! Sr
                       0 ,& ! Y
                       0 ,& ! Zr
                       0 ,& ! Nb
                       0 ,& ! Mo
                       0 ,& ! Tc
                       0 ,& ! Ru
                       0 ,& ! Rh
                       0, & ! Pd
                       0, & ! Ag
                       0, & ! Cd
                       0, & ! In
                       0, & ! Sn
                       0, & ! Sb
                       0 ,& ! Te
                       0 ,& ! I
                       0 ,& ! Xe 
                       0 ,& ! Cs
                       0 ,& ! Ba
                       0 ,& ! La
                       0 ,& ! Ce
                       0 ,& ! Pr 
                       0 ,& ! Nd
                       0 ,& ! Pm
                       0 ,& ! Sm
                       0 ,& ! Eu
                       0 ,& ! Gd 
                       0 ,& ! Tb
                       0 ,& ! Dy
                       0 ,& ! Ho
                       0 ,& ! Er
                       0 ,& ! Tm 
                       0 ,& ! Yb
                       0 ,& ! Lu
                       0 ,& ! Hf
                       0 ,& ! Ta
                       0 ,& ! W 
                       0 ,& ! Re
                       0 ,& ! Os
                       0 ,& ! Ir
                       0 ,& ! Pt
                       0 ,& ! Au 
                       0 ,& ! Hg
                       0 ,& ! Tl
                       0 ,& ! Pb
                       0 ,& ! Bi
                       0 ,& ! Po 
                       0 ,& ! At
                       0 ,& ! Rn
                       0 ,& ! Fr
                       0 ,& ! Ra
                       0 ,& ! Ac 
                       0 ,& ! Th
                       0 ,& ! Pa
                       0 ,& ! U
                       0 ,& ! Np
                       0 ,& ! Pu 
                       0 ,& ! Am
                       0  & ! Cm
                       /)

  integer, dimension(1:96)  :: upper_ion= (/&
                       2 ,& ! H 
                       3 ,& ! He 
                       0 ,& ! Li 
                       0 ,& ! Be 
                       0 ,& ! B
                       4 ,& ! C
                       4 ,& ! N
                       4 ,& ! O
                       0 ,& ! F
                       4 ,& ! Ne
                       0 ,& ! Na 
                       0 ,& ! Mg
                       0 ,& ! Al
                       0 ,& ! Si
                       0 ,& ! P
                       4 ,& ! S
                       0 ,& ! Cl 
                       0 ,& ! Ar
                       0 ,& ! K
                       0 ,& ! Ca
                       0 ,& ! Sc
                       0 ,& ! Ti 
                       0 ,& ! V
                       0 ,& ! Cr
                       0 ,& ! Mn
                       0 ,& ! Fe
                       0 ,& ! Co 
                       0 ,& ! Ni
                       0 ,& ! Cu
                       0 ,& ! Zn
                       0 ,& ! Ga
                       0 ,& ! Ge 
                       0 ,& ! As
                       0 ,& ! Se
                       0 ,& ! Br
                       0 ,& ! Kr
                       0 ,& ! Rb
                       0 ,& ! Sr
                       0 ,& ! Y
                       0 ,& ! Zr
                       0 ,& ! Nb
                       0 ,& ! Mo
                       0 ,& ! Tc
                       0 ,& ! Ru
                       0 ,& ! Rh
                       0 ,& ! Pd
                       0 ,& ! Ag
                       0, & ! Cd
                       0, & ! In
                       0, & ! Sn
                       0 ,& ! Sb
                       0 ,& ! Te
                       0 ,& ! I
                       0 ,& ! Xe 
                       0 ,& ! Cs
                       0 ,& ! Ba
                       0 ,& ! La
                       0 ,& ! Ce
                       0 ,& ! Pr 
                       0 ,& ! Nd
                       0 ,& ! Pm
                       0 ,& ! Sm
                       0 ,& ! Eu
                       0 ,& ! Gd 
                       0 ,& ! Tb
                       0 ,& ! Dy
                       0 ,& ! Ho
                       0 ,& ! Er
                       0 ,& ! Tm 
                       0 ,& ! Yb
                       0 ,& ! Lu
                       0 ,& ! Hf
                       0 ,& ! Ta
                       0 ,& ! W 
                       0 ,& ! Re
                       0 ,& ! Os
                       0 ,& ! Ir
                       0 ,& ! Pt
                       0 ,& ! Au 
                       0 ,& ! Hg
                       0 ,& ! Tl
                       0 ,& ! Pb
                       0 ,& ! Bi
                       0 ,& ! Po 
                       0 ,& ! At
                       0 ,& ! Rn
                       0 ,& ! Fr
                       0 ,& ! Ra
                       0 ,& ! Ac 
                       0 ,& ! Th
                       0 ,& ! Pa
                       0 ,& ! U
                       0 ,& ! Np
                       0 ,& ! Pu 
                       0 ,& ! Am
                       0  & ! Cm
          /)

 

      REAL(dp), DIMENSION (96), PARAMETER :: ALL_AWEIGHT = &
          (/   1.00794_dp  , & !  1  H   Hydrogen
               4.002602_dp , & !  2  He  Helium
               6.941_dp    , & !  3  Li  Lithium
               9.012182_dp , & !  4  Be  Beryllium
              10.811_dp    , & !  5  B   Boron
              12.0107_dp   , & !  6  C   Carbon
              14.0067_dp   , & !  7  N   Nitrogen
              15.9994_dp   , & !  8  O   Oxygen
              18.9984032_dp, & !  9  F   Fluorine
              20.1797_dp   , & ! 10  Ne  Neon
              22.989770_dp , & ! 11  Na  Sodium
              24.3050_dp   , & ! 12  Mg  Magnesium
              26.981538_dp , & ! 13  Al  Aluminium
              28.0855_dp   , & ! 14  Si  Silicon
              30.973761_dp , & ! 15  P   Phosphorus
              32.065_dp    , & ! 16  S   Sulfur
              35.453_dp    , & ! 17  Cl  Chlorine
              39.948_dp    , & ! 18  Ar  Argon
              39.0983_dp   , & ! 19  K   Potassium
              40.078_dp    , & ! 20  Ca  Calcium
              44.955910_dp , & ! 21  Sc  Scandium
              47.867_dp    , & ! 22  Ti  Titanium
              50.9415_dp   , & ! 23  V   Vanadium
              51.9961_dp   , & ! 24  Cr  Chromium
              54.938049_dp , & ! 25  Mn  Manganese
              55.845_dp    , & ! 26  Fe  Iron
              58.933200_dp , & ! 27  Co  Cobalt
              58.6934_dp   , & ! 28  Ni  Nickel
              63.546_dp    , & ! 29  Cu  Copper
              65.409_dp    , & ! 30  Zn  Zinc
              69.723_dp    , & ! 31  Ga  Gallium
              72.64_dp     , & ! 32  Ge  Germanium
              74.92160_dp  , & ! 33  As  Arsenic
              78.96_dp     , & ! 34  Se  Selenium
              79.904_dp    , & ! 35  Br  Bromine
              83.798_dp    , & ! 36  Kr  Krypton
              85.4678_dp   , & ! 37  Rb  Rubidium
              87.62_dp     , & ! 38  Sr  Strontium
              88.90585_dp  , & ! 39  Y   Yttrium
              91.224_dp    , & ! 40  Zr  Zirconium
              92.90638_dp  , & ! 41  Nb  Niobium
              95.94_dp     , & ! 42  Mo  Molybdenum
              98.0_dp      , & ! 43  Tc  Technetium
             101.07_dp     , & ! 44  Ru  Ruthenium
             102.90550_dp  , & ! 45  Rh  Rhodium
             106.42_dp     , & ! 46  Pd  Palladium
             107.8682_dp   , & ! 47  Ag  Silver
             112.411_dp    , & ! 48  Cd  Cadmium
             114.818_dp    , & ! 49  In  Indium
             118.710_dp    , & ! 50  Sn  Tin
             121.760_dp    , & ! 51  Sb  Antimony
             127.60_dp     , & ! 52  Te  Tellurium
             126.90447_dp  , & ! 53  I   Iodine
             131.293_dp    , & ! 54  Xe  Xenon
             132.90545_dp  , & ! 55  Cs  Caesium
             137.327_dp    , & ! 56  Ba  Barium
             138.9055_dp   , & ! 57  La  Lanthanum
             140.116_dp    , & ! 58  Ce  Cerium
             140.90765_dp  , & ! 59  Pr  Praseodymium
             144.24_dp     , & ! 60  Nd  Neodymium
             145.0_dp      , & ! 61  Pm  Promethium
             150.36_dp     , & ! 62  Sm  Samarium
             151.964_dp    , & ! 63  Eu  Europium
             157.25_dp     , & ! 64  Gd  Gadolinium
             158.92534_dp  , & ! 65  Tb  Terbium
             162.500_dp    , & ! 66  Dy  Dysprosium
             164.93032_dp  , & ! 67  Ho  Holmium
             167.259_dp    , & ! 68  Er  Erbium
             168.93421_dp  , & ! 69  Tm  Thulium
             173.04_dp     , & ! 70  Yb  Ytterbium
             174.967_dp    , & ! 71  Lu  Lutetium
             178.49_dp     , & ! 72  Hf  Hafnium
             180.9479_dp   , & ! 73  Ta  Tantalum
             183.84_dp     , & ! 74  W   Tungsten
             186.207_dp    , & ! 75  Re  Rhenium
             190.23_dp     , & ! 76  Os  Osmium
             192.217_dp    , & ! 77  Ir  Iridium
             195.078_dp    , & ! 78  Pt  Platinum
             196.96655_dp  , & ! 79  Au  Gold
             200.59_dp     , & ! 80  Hg  Mercury
             204.3833_dp   , & ! 81  Tl  Thallium
             207.2_dp      , & ! 82  Pb  Lead
             208.98038_dp  , & ! 83  Bi  Bismuth
             209.0_dp      , & ! 84  Po  Polonium
             210.0_dp      , & ! 85  At  Astatine
             222.0_dp      , & ! 86  Rn  Radon
             223.0_dp      , & ! 87  Fr  Francium
             226.0_dp      , & ! 88  Ra  Radium
             227.0_dp      , & ! 89  Ac  Actinium
             232.0381_dp   , & ! 90  Th  Thorium
             231.03588_dp  , & ! 91  Pa  Protactinium
             238.02891_dp  , & ! 92  U   Uranium
             0._dp         , & ! 93  Np  Neptunium
             0._dp         , & ! 94  Pu  Plutonium
             0._dp         , & ! 95  Am  Americium
             0._dp           & ! 96  Cm  Curium 
      /)
!
         ! SOLAR ABUNDANCES OF THE ELEMENTS - 1998.
         ! Standard solar composition from N. Grevesse and A. Sauval
         ! in: Space Science Reviews 85, 161, 1998.
         ! The values shown as comments on the right hand side
         ! are from Holweger, 1979.
         ! The values of the element abundances EL_sol are scaled to
         ! H_sol = 12 in the following way EL_sol=10^(EL_sol-H_sol).
         !
         ! List of Elements in Atomic Number Order:
         !=========================================
      real(dp), dimension (96), parameter :: all_solar_abund_old = &
          (/   12.00_dp , & !  1  H   Hydrogen
               10.93_dp , & !  2  He  Helium
                1.10_dp , & !  3  Li  Lithium
                1.40_dp , & !  4  Be  Beryllium
                2.55_dp , & !  5  B   Boron
                8.52_dp , & !  6  C   Carbon
                7.92_dp , & !  7  N   Nitrogen
                8.83_dp , & !  8  O   Oxygen
                4.56_dp , & !  9  F   Fluorine
                8.08_dp , & ! 10  Ne  Neon
                6.33_dp , & ! 11  Na  Sodium
                7.58_dp , & ! 12  Mg  Magnesium
                6.47_dp , & ! 13  Al  Aluminium
                7.55_dp , & ! 14  Si  Silicon
                5.45_dp , & ! 15  P   Phosphorus
                7.33_dp , & ! 16  S   Sulfur
                5.50_dp , & ! 17  Cl  Chlorine
                6.40_dp , & ! 18  Ar  Argon
                5.12_dp , & ! 19  K   Potassium
                6.36_dp , & ! 20  Ca  Calcium
                3.17_dp , & ! 21  Sc  Scandium
                5.02_dp , & ! 22  Ti  Titanium
                4.00_dp , & ! 23  V   Vanadium
                5.67_dp , & ! 24  Cr  Chromium
                5.39_dp , & ! 25  Mn  Manganese
                7.50_dp , & ! 26  Fe  Iron
                4.92_dp , & ! 27  Co  Cobalt
                6.25_dp , & ! 28  Ni  Nickel
                4.21_dp , & ! 29  Cu  Copper
                4.60_dp , & ! 30  Zn  Zinc
                2.88_dp , & ! 31  Ga  Gallium
                3.41_dp , & ! 32  Ge  Germanium
              -10.00_dp , & ! 33  As  Arsenic
              -10.00_dp , & ! 34  Se  Selenium
              -10.00_dp , & ! 35  Br  Bromium
              -10.00_dp , & ! 36  Kr  Krypton
                2.60_dp , & ! 37  Rb  Rubidium
                2.97_dp , & ! 38  Sr  Strontium
                2.24_dp , & ! 39  Y   Yttrium
                2.60_dp , & ! 40  Zr  Zirconium
                1.42_dp , & ! 41  Nb  Niobium
                1.92_dp , & ! 42  Mo  Molybdenum
                1.84_dp , & ! 43  Tc  Technetium
                1.12_dp , & ! 44  Ru  Ruthenium
                1.69_dp , & ! 45  Rh  Rhodium
                0.94_dp , & ! 46  Pd  Palladium
                1.77_dp , & ! 47  Ag  Silver
                1.77_dp , & ! 48  Cd  Cadmium
                1.66_dp , & ! 49  In  Indium
                2.00_dp , & ! 50  Sn  Tin
                1.00_dp , & ! 51  Sb  Antimony
              -10.00_dp , & ! 52  Te  Tellurium 
              -10.00_dp , & ! 53  I   Iodium
              -10.00_dp , & ! 54  Xe  Xenon
              -10.00_dp , & ! 55  Cs  Caesium
                2.13_dp , & ! 56  Ba  Barium
                1.17_dp , & ! 57  La  Lanthanum
                1.58_dp , & ! 58  Ce  Cerium
                0.71_dp , & ! 59  Pr  Praesodynium
                1.50_dp , & ! 60  Nd  Neodynium
              -10.00_dp , & ! 61  Pm  Promethium
                1.01_dp , & ! 62  Sm  Samarium
                0.51_dp , & ! 63  Eu  Europium
                1.12_dp , & ! 64  Gd  Gadlinium
               -0.10_dp , & ! 65  Tb  Terbium
                1.14_dp , & ! 66  Dy  Dysprosium
                0.26_dp , & ! 67  Ho  Holmium
                0.93_dp , & ! 68  Er  Erbium
                0.00_dp , & ! 69  Tm  Thulium
                1.08_dp , & ! 70  Yb  Ytterbium
                0.06_dp , & ! 71  Lu  Luthetium
                0.88_dp , & ! 72  Hf  Haffnium
              -10.00_dp , & ! 73  Ta  Tanatalum    
                1.11_dp , & ! 74  W   Tungsten
              -10.00_dp , & ! 75  Re  Rhenium
                1.45_dp , & ! 76  Os  Osmium
                1.35_dp , & ! 77  Ir  Iridium
                1.80_dp , & ! 78  Pt  Platinum
                1.01_dp , & ! 79  Au  Gold
              -10.00_dp , & ! 80  Hg  Mercury
                0.90_dp , & ! 81  Tl  Thallium
                1.95_dp , & ! 82  Pb  Lead
              -10.00_dp , & ! 83  Bi  Bismuthum
              -10.00_dp , & ! 84  Po  Polonium
              -10.00_dp , & ! 85  At  Astatine
              -10.00_dp , & ! 86  Rn  Radon
              -10.00_dp , & ! 87  Fr  Francium
              -10.00_dp , & ! 88  Ra  Radium
              -10.00_dp , & ! 89  Ac  Actinium
              -10.00_dp , & ! 90  Th  Thorium
              -10.00_dp , & ! 91  Pa  Protactinium
               -0.47_dp , & ! 92  U   Uranium
              -10.00_dp,  & ! 93  Np  Neptunium
              -10.00_dp , & ! 94  Pu  Plutonium
              -10.00_dp , & ! 95  Am  Americium
              -10.00_dp   & ! 96  Cm  Curium 
                /)
                

         ! SOLAR ABUNDANCES OF THE ELEMENTS - 2009.
         ! Standard solar composition from M. Asplund, N. Grevesse, 
         ! A.J. Sauval, and  P. Scott
         ! "The Chemical Composition of the Sun"
         ! Annu. Rev. Astron. Astrophys. 2009, 47:481-522.
         ! To maintain compatibility between older WM-Basic test runs
         ! and newer ones the data from 1998 can still be used.
         ! The middle column contains the data from
         ! N. Grevesse and A. Sauval
         ! in: Space Science Reviews 85, 161, 1998,
         ! and the right-hand column the data from Hollweger 1979. 
         ! The values of the element abundances EL_sol are scaled to
         ! H_sol = 12 in the following way EL_sol=10^(EL_sol-H_sol).
         !
         ! List of Elements in Atomic Number Order:
         !=========================================
      real(dp), dimension (96), parameter :: all_solar_abund = &
          (/   12.00_dp , & !  1  H   Hydrogen
               10.93_dp , & !  2  He  Helium
                1.05_dp , & !  3  Li  Lithium
                1.38_dp , & !  4  Be  Beryllium
                2.70_dp , & !  5  B   Boron
                8.43_dp , & !  6  C   Carbon
                7.83_dp , & !  7  N   Nitrogen
                8.69_dp , & !  8  O   Oxygen
                4.56_dp , & !  9  F   Fluorine
                7.93_dp , & ! 10  Ne  Neon
                6.24_dp , & ! 11  Na  Sodium
                7.60_dp , & ! 12  Mg  Magnesium
                6.45_dp , & ! 13  Al  Aluminium
                7.51_dp , & ! 14  Si  Silicon
                5.41_dp , & ! 15  P   Phosphorus
                7.12_dp , & ! 16  S   Sulfur
                5.50_dp , & ! 17  Cl  Chlorine
                6.40_dp , & ! 18  Ar  Argon
                5.03_dp , & ! 19  K   Potassium
                6.34_dp , & ! 20  Ca  Calcium
                3.15_dp , & ! 21  Sc  Scandium
                4.95_dp , & ! 22  Ti  Titanium
                3.93_dp , & ! 23  V   Vanadium
                5.46_dp , & ! 24  Cr  Chromium
                5.43_dp , & ! 25  Mn  Manganese
                7.50_dp , & ! 26  Fe  Iron
                4.99_dp , & ! 27  Co  Cobalt
                6.22_dp , & ! 28  Ni  Nickel
                4.19_dp , & ! 29  Cu  Copper
                4.56_dp , & ! 30  Zn  Zinc
                3.04_dp , & ! 31  Ga  Gallium
                3.65_dp , & ! 32  Ge  Germanium
              -10.00_dp , & ! 33  As  Arsenic
              -10.00_dp , & ! 34  Se  Selenium
              -10.00_dp , & ! 35  Br  Bromium
                3.25_dp , & ! 36  Kr  Krypton
                2.52_dp , & ! 37  Rb  Rubidium
                2.87_dp , & ! 38  Sr  Strontium
                2.21_dp , & ! 39  Y   Yttrium
                2.58_dp , & ! 40  Zr  Zirconium
                1.46_dp , & ! 41  Nb  Niobium
                1.88_dp , & ! 42  Mo  Molybdenum
              -10.00_dp , & ! 43  Tc  Technetium
                1.75_dp , & ! 44  Ru  Ruthenium
                0.91_dp , & ! 45  Rh  Rhodium
                1.57_dp , & ! 46  Pd  Palladium
                0.94_dp , & ! 47  Ag  Silver
              -10.00_dp , & ! 48  Cd  Cadmium
                0.80_dp , & ! 49  In  Indium
                2.04_dp , & ! 50  Sn  Tin
              -10.00_dp , & ! 51  Sb  Antimony
              -10.00_dp , & ! 52  Te  Tellurium 
              -10.00_dp , & ! 53  I   Iodium
                2.24_dp , & ! 54  Xe  Xenon
              -10.00_dp , & ! 55  Cs  Caesium
                2.18_dp , & ! 56  Ba  Barium
                1.10_dp , & ! 57  La  Lanthanum
                1.58_dp , & ! 58  Ce  Cerium
                0.72_dp , & ! 59  Pr  Praesodynium
                1.42_dp , & ! 60  Nd  Neodynium
              -10.00_dp , & ! 61  Pm  Promethium
                0.96_dp , & ! 62  Sm  Samarium
                0.52_dp , & ! 63  Eu  Europium
                1.07_dp , & ! 64  Gd  Gadlinium
                0.30_dp , & ! 65  Tb  Terbium
                1.10_dp , & ! 66  Dy  Dysprosium
                0.48_dp , & ! 67  Ho  Holmium
                0.92_dp , & ! 68  Er  Erbium
                0.10_dp , & ! 69  Tm  Thulium
                0.84_dp , & ! 70  Yb  Ytterbium
                0.10_dp , & ! 71  Lu  Luthetium
                0.85_dp , & ! 72  Hf  Haffnium
              -10.00_dp , & ! 73  Ta  Tanatalum    
                0.85_dp , & ! 74  W   Tungsten
              -10.00_dp , & ! 75  Re  Rhenium
                1.40_dp , & ! 76  Os  Osmium
                1.38_dp , & ! 77  Ir  Iridium
              -10.00_dp , & ! 78  Pt  Platinum
                0.92_dp , & ! 79  Au  Gold
              -10.00_dp , & ! 80  Hg  Mecury
                0.90_dp , & ! 81  Tl  Thallium
                1.75_dp , & ! 82  Pb  Lead
              -10.00_dp , & ! 83  Bi  Bismuthum
              -10.00_dp , & ! 84  Po  Polonium
              -10.00_dp , & ! 85  At  Astatine
              -10.00_dp , & ! 86  Rn  Radon
              -10.00_dp , & ! 87  Fr  Francium
              -10.00_dp , & ! 88  Ra  Radium
              -10.00_dp , & ! 89  Ac  Actinium
              -10.00_dp , & ! 90  Th  Thorium
              -10.00_dp , & ! 91  Pa  Protactinium
              -10.00_dp , & ! 92  U   Uranium
              -10.00_dp , & ! 93  Np  Neptunium
              -10.00_dp , & ! 94  Pu  Plutonium
              -10.00_dp , & ! 95  Am  Americium
              -10.00_dp   & ! 96  Cm  Curium 
                /)
 

!
         ! List of Elements considered in Atomic Number Order:
         !====================================================
      CHARACTER(2), DIMENSION(96), PARAMETER :: ALL_STRING1 = &
              (/   ' H' , & !  1  H   Hydrogen
                   'He' , & !  2  He  Helium
                   'Li' , & !  3  Li  Lithium
                   'Be' , & !  4  Be  Beryllium
                   ' B' , & !  5  B   Boron
                   ' C' , & !  6  C   Carbon
                   ' N' , & !  7  N   Nitrogen
                   ' O' , & !  8  O   Oxygen
                   ' F' , & !  9  F   Fluorine
                   'Ne' , & ! 10  Ne  Neon
                   'Na' , & ! 11  Na  Sodium
                   'Mg' , & ! 12  Mg  Magnesium
                   'Al' , & ! 13  Al  Aluminium
                   'Si' , & ! 14  Si  Silicon
                   ' P' , & ! 15  P   Phosphorus
                   ' S' , & ! 16  S   Sulfur
                   'Cl' , & ! 17  Cl  Chlorine
                   'Ar' , & ! 18  Ar  Argon
                   ' K' , & ! 19  K   Potassium
                   'Ca' , & ! 20  Ca  Calcium
                   'Sc' , & ! 21  Sc  Scandium
                   'Ti' , & ! 22  Ti  Titanium
                   ' V' , & ! 23  V   Vanadium
                   'Cr' , & ! 24  Cr  Chromium
                   'Mn' , & ! 25  Mn  Manganese
                   'Fe' , & ! 26  Fe  Iron
                   'Co' , & ! 27  Co  Cobalt
                   'Ni' , & ! 28  Ni  Nickel
                   'Cu' , & ! 29  Cu  Copper
                   'Zn' , & ! 30  Zn  Zinc
                   'Ga' , & ! 31  Ga  Gallium
                   'Ge' , & ! 32  Ge  Germanium
                   'As' , & ! 33  As  Arsenic
                   'Se' , & ! 34  Se  Selenium
                   'Br' , & ! 35  Br  Bromium
                   'Kr' , & ! 36  Kr  Krypton
                   'Rb' , & ! 37  Rb  Rubidium
                   'Sr' , & ! 38  Sr  Strontium
                   ' Y' , & ! 39  Y   Yttrium
                   'Zr' , & ! 40  Zr  Zirconium
                   'Nb' , & ! 41  Nb  Niobium
                   'Mo' , & ! 42  Mo  Molybdenum
                   'Tc' , & ! 43  Tc  Technetium
                   'Ru' , & ! 44  Ru  Ruthenium
                   'Rh' , & ! 45  Rh  Rhodium
                   'Pd' , & ! 46  Pd  Palladium
                   'Ag' , & ! 47  Ag  Silver
                   'Cd' , & ! 48  Cd  Cadmium
                   'In' , & ! 49  In  Indium
                   'Sn' , & ! 50  Sn  Tin
                   'Sb' , & ! 51  Sb  Antimony
                   'Te' , & ! 52  Te  Tellurium 
                   ' I' , & ! 53  I   Iodium
                   'Xe' , & ! 54  Xe  Xenon
                   'Cs' , & ! 55  Cs  Caesium
                   'Ba' , & ! 56  Ba  Barium
                   'La' , & ! 57  La  Lanthanum
                   'Ce' , & ! 58  Ce  Cerium
                   'Pr' , & ! 59  Pr  Praesodynium
                   'Nd' , & ! 60  Nd  Neodynium
                   'Pm' , & ! 61  Pm  Promethium
                   'Sm' , & ! 62  Sm  Samarium
                   'Eu' , & ! 63  Eu  Europium
                   'Gd' , & ! 64  Gd  Gadlinium
                   'Tb' , & ! 65  Tb  Terbium
                   'Dy' , & ! 66  Dy  Dysprosium
                   'Ho' , & ! 67  Ho  Holmium
                   'Er' , & ! 68  Er  Erbium
                   'Tm' , & ! 69  Tm  Thulium
                   ' Y' , & ! 70  Yb  Ytterbium
                   'Lu' , & ! 71  Lu  Luthetium
                   'Hf' , & ! 72  Hf  Haffnium
                   'Ta' , & ! 73  Ta  Tanatalum    
                   ' W' , & ! 74  W   Tungsten
                   'Re' , & ! 75  Re  Rhenium
                   'Os' , & ! 76  Os  Osmium
                   'Ir' , & ! 77  Ir  Iridium
                   'Pt' , & ! 78  Pt  Platinum
                   'Au' , & ! 79  Au  Gold
                   'Hg' , & ! 80  Hg  Mecury
                   'Tl' , & ! 81  Tl  Thallium
                   'Pb' , & ! 82  Pb  Lead
                   'Bi' , & ! 83  Bi  Bismuthum
                   'Po' , & ! 84  Po  Polonium
                   'At' , & ! 85  At  Astatine
                   'Rn' , & ! 86  Rn  Radon
                   'Fr' , & ! 87  Fr  Francium
                   'Ra' , & ! 88  Ra  Radium
                   'Ac' , & ! 89  Ac  Actinium
                   'Th' , & ! 90  Th  Thorium
                   'Pa' , & ! 91  Pa  Protactinium
                   ' U' , & ! 92  U   Uranium
                   'Np' , & ! 93  Np  Neptunium
                   'Pu' , & ! 94  Pu  Plutonium
                   'Am' , & ! 95  Am  Americium
                   'Cu'   & ! 96  Cm  Curium
                /) 
 

 
             

end module M_chemical_data
