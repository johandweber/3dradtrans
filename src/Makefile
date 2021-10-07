#########################################################################
ifeq ($(PLATFORM), linux_ifort)

# Linux programs compiled with ifort, the Intel Fortran compiler.
# Newer versions of the Intel fortran compiler packages
# include the MKL library that includes the BLAS and
# LAPACK routines. The MKL library is highly optimized for the
# modern x86 compatible PCs and likely to be faster than libraries
# compiled from the more generic code form netlib.org (see
# the descriptions of the platforms below)
# The Makefile was tested using the bash shell. There may or may not be
# problems with other shells.


  OBJEXT   =.o  #     # Extension of object file
  MODEXT   =.mod  #   # Extension of fortran module files
  EXEEXT   = #        # Extension of executable file
  COMPILER = ifort #  # Name of compiler
  OPENMP   = -qopenmp # OpenMP switch(es)
  NOLINK   = -c #     # Just compile, not link
  NAMING   = -o #     # Option for giving the ouput file a certain name
  VERSION  = -V #     # Showing the version
  PROGRAM  = occup_3d_ext$(EXEXT) # name of the program
  ifndef MKL_ROOT
    MKL_ROOT = -mkl=parallel
  endif
  LAPACK = $(MKL_ROOT)
  ifdef debug
    OPTIMIZE = -static -O0 -g # Optimization switch for compiler
    LINKER   = ifort -static -O0
    DEBUG    = -g -fpe0 -warn all -traceback -debug extended -check all
  else
# -march=core-avx2 means that AVX and AVX2 are supported
    OPTIMIZE = -static -O3 #-march=core-avx2   # Optimization switch for compiler
    LINKER   = ifort -static -Ofast
    DEBUG    =
  endif
# if perfana is set (perfana=1), the option "-p" is used
# for compilation and linking. In this case the program
# supports profiling using the gprof tool.
  ifdef perfana
    GPROF   = -p
  else
    GPROF   =
  endif
  ECHO      = echo 
  REMOVE    = rm -f # -f, because else there is an error when there
                    # are no files (of  a certain pattern) to delete,
                    # preventing the clean target to be completed


#########################################################################
else ifeq ($(PLATFORM), linux_ifx)

# Linux programs compiled with ifx, the LLVM-based Intel Fortran compiler.
# Newer versions of the Intel fortran compiler packages
# include the MKL library that includes the BLAS and
# LAPACK routines. The MKL library is highly optimized for the
# modern x86 compatible PCs and likely to be faster than libraries
# compiled from the more generic code form netlib.org (see
# the descriptions of the platforms below)
# The Makefile was tested using the bash shell. There may or may not be
# problems with other shells.
#
# The support for IFX is not stable yet.


  OBJEXT   =.o  #     # Extension of object file
  MODEXT   =.mod  #   # Extension of fortran module files
  EXEEXT   = #        # Extension of executable file
  COMPILER = ifx #  # Name of compiler
  OPENMP   = -qopenmp # OpenMP switch(es)
  NOLINK   = -c #     # Just compile, not link
  NAMING   = -o #     # Option for giving the ouput file a certain name
  VERSION  = -V #     # Showing the version
  PROGRAM  = occup_3d_ext$(EXEXT) # name of the program
  ifndef MKL_ROOT
    MKL_ROOT = -mkl=parallel
  endif
  LAPACK = $(MKL_ROOT)
  ifdef debug
    OPTIMIZE = -static -O0 -g # Optimization switch for compiler
    LINKER   = ifort -static -O0
    DEBUG    = -g -fpe0 -warn all -traceback -debug extended -check all
  else
    OPTIMIZE = -static -O3 #   # Optimization switch for compiler
    LINKER   = ifort -static -O3
    DEBUG    =
  endif
# if perfana is set (perfana=1), the option "-p" is used
# for compilation and linking. In this case the program
# supports profiling using the gprof tool.
  ifdef perfana
    GPROF   = -p
  else
    GPROF   =
  endif
  ECHO      = echo 
  REMOVE    = rm -f # -f, because else there is an error when there
                    # are no files (of  a certain pattern) to delete,
                    # preventing the clean target to be completed


#########################################################################
else ifeq ($(PLATFORM), linux_aocc)

# Linux programs compiled with AMD's flang Fortran compiler.
# It is assumed that the BLIS library from the AOCL library collection
# is present for the BLAS routines. For LAPACK the presence of a
# linbapack is assumed.
#
# The support for this target is not stable yet.

  OBJEXT   =.o  #     # Extension of object file
  MODEXT   =.mod  #   # Extension of fortran module files
  EXEEXT   = #        # Extension of executable file
  COMPILER = flang      # Name of compiler
  OPENMP   = -fopenmp # OpenMP switch(es)
  NOLINK   = -c #     # Just compile, not link
  NAMING   = -o #     # Option for giving the ouput file a certain name
  VERSION  = -V #     # Showing the version
  PROGRAM  = occup_3d_ext$(EXEXT) # name of the program
  ifndef MKL_ROOT
    MKL_ROOT = -llapack -lblis
  endif
  LAPACK = $(MKL_ROOT)
  ifdef debug
    OPTIMIZE = -static -O0 -g # Optimization switch for compiler
    LINKER   = flang -static
    DEBUG    = -g -fpe0 -warn all -traceback -debug extended -check all
  else
    OPTIMIZE =  -O3 #-march=native #   # Optimization switch for compiler
    LINKER   = flang -O3
    DEBUG    =
  endif
# if perfana is set (perfana=1), the option "-p" is used
# for compilation and linking. In this case the program
# supports profiling using the gprof tool.
  ifdef perfana
    GPROF   = -p
  else
    GPROF   =
  endif
  ECHO      = echo 
  REMOVE    = rm -f # -f, because else there is an error when there
                    # are no files (of  a certain pattern) to delete,
                    # preventing the clean target to be completed



#########################################################################

# Linux programs compiled with gfortran. gfortran is a part of the
# gcc suite and is thus available for free in the repositories of all
# maior Linux distributions.
# Furthermore it assumes the presence of the standard BLAS and LAPACK
# libraries. These linear algebra libraries are also likely to be
# avaiable in the standard libraries. Otherwise you can download them from
# http://netlib.org and compile them yourself. To do so you need just
# a Fortran90 and a C compiler.
# The Makefile was tested using the bash shell. There may or may not be
# problems with other shells.

else ifeq ($(PLATFORM), linux_gfortran)
  OBJEXT   =.o #      # Extension of object file
  MODEXT   =.mod #    # Extension of fortran module files
  EXEEXT   = #        # Extension of executable file
  COMPILER = gfortran # Name of compiler
  OPENMP   = -fopenmp # OpenMP switch(es)
  NOLINK   = -c #     # Just compile, not link
  NAMING   = -o #     # Option for giving the ouput file a certain name
  VERSION  = -v #     # Showing the version
  PROGRAM  = occup_3d_ext$(EXEXT) # name of the program

# OpenBLAS is an alternative Library that contains both the BLAS and the LAPACK functions.
  LAPACK =  -llapack -lblas # -lopenblas
  ifdef debug
    OPTIMIZE = -O0  # Optimization switch for compiler
    LINKER   = gfortran -O0
    DEBUG    = -g3 -Wall -pedantic  -fbacktrace  -fcheck=all -ffpe-trap=invalid,zero,overflow
  else
# If -march=native is commented in, the created binary may no longer work on other CPUs    
    OPTIMIZE = -O3   -march=native   # Optimization switch for compiler
    LINKER   = gfortran
    DEBUG    =
  endif
# if perfana is set (perfana=1), the option "-p" is used
# for compilation and linking. In this case the program
# supports profiling using the gprof tool.
  ifdef perfana
    GPROF   = -pg
  else
    GPROF   =
  endif 
  REMOVE    = rm -f # -f, because else there is an error when there
                    # are no files (of  a certain pattern) to delete,
                    # preventing the clean target to be completed

########################################################################

# Native Windows programs, but assuming there exists a bash shell and
# the gcc compiler is available. This Makefile was tested with the
# MSYS2 environment freely available for windows. So the only difference
# between to the linux_gcc case is the .exe suffix for executable
# files and the use of -lopenblas instead of -llapack -lblas.

else ifeq ($(PLATFORM), windows_msys2_gfortran) 
  OBJEXT   =.o #                 # Extension of object file
  MODEXT   =.mod #               # Extension of fortran module files
  EXEEXT   =.exe #               # Extension of executable file
  COMPILER = gfortran -static -s # Name of compiler, the option "-static" has been
                                 # added for cases where a single command triggers
                                 # both the compilig and linking
                                 # The -s option removes unneccessary symbols from
                                 # the executables, thus reducing the file size
                                 # (The impact may, however, be small.)
  OPENMP   = -fopenmp            # OpenMP switch(es)
  NOLINK   = -c #                # Just compile, not link
  NAMING   = -o #                # Option for giving the ouput file a certain name
  VERSION  = -v #                # Showing the version
  PROGRAM  = occup_3d_ext$(EXEXT) # name of the program
  LAPACK = -lopenblas
  ifdef debug
    OPTIMIZE = -O0  # Optimization switch for compiler
    LINKER   = gfortran -O0
    DEBUG    = -g3 -Wall -pedantic  -fbacktrace  -fcheck=all -ffpe-trap=invalid,zero,overflow
  else
    OPTIMIZE = -O3 # Optimization switch for compiler
    LINKER   = gfortran -static -s
    DEBUG    =
  endif
# if perfana is set (perfana=1), the option "-p" is used
# for compilation and linking. In this case the program
# supports profiling using the gprof tool.
  ifdef perfana
    GPROF   = -pg
  else
    GPROF   =
  endif
  ECHO      = echo
  REMOVE    = rm -f # -f, because else there is an error when there
                    # are no files (of  a certain pattern) to delete,
                    # preventing the clean target to be completed


########################################################################

# Here you can add additional platforms by copying and modifying the
# variable definitions from above.
# We recommend the naming convention
# PLATFORM= operting_system[_subsystem]_compiler[special_optimizations]´
# (the parts in square barckets are optional)
# If you have ported the program to a new environment, we would
# appreciate if you contributed to the project and sent us the
# neccessary definitions (info@3dradtrans.de)
# However, this is not a requirement.

else
    $(info Unknown platform or no platform specified.)
    $(info available platforms:)
    $(info PLATFORM=linux_ifort)
    $(info PLATFORM=linux_ifx)
    $(info PLATFORM=linux_gfortran)
    $(info PLATFORM=windows_msys2_gfortran)
    $(error make process terminated)
endif

#########################################################################

# This is the platform-independent partthat describes the dependencies
# within the project

data_types_f90           = data_types.f90
modauxmath_f90           = modauxmath.f90
modr23_f90               = modr23.f90

definitions_f90          = definitions.f90
chemical_data_f90        = chemical_data.f90
natural_constants_f90    = natural_constants.f90
min_healpix_f90          = min_healpix.f90
geometry_f90             = geometry.f90
data_input_f90           = data_input.f90
ifdef legacy
  absorption_f90         = absorption.f90
else
  absorption_f90         = absorption2.f90
endif
ifdef legacy
  hydrogen_f90           = hydrogen.f90
else
  hydrogen_f90           = hydrogen2.f90
endif
output_f90               = output.f90
subsystems_f90           = subsystems.f90
OpenMP_analytic_f90      = OpenMP_analytic.f90

modblackbody_f90         = modblackbody.f90
modstrongline_f90        = modstrongline.f90

absolutetorelative_f90   = absolutetorelative.f90
createtemperature_f90    = createtemperature.f90
diskcut_f90              = diskcut.f90
fractalgenerator_f90     = fractalgenerator.f90
nebstat_f90              = nebstat.f90
radiusbench_f90          = radiusbench.f90
sb2ms_f90                = sb2ms.f90
threedhydro_ext_f90      = 3dhydro_ext.f90
oiilevels_f90            = oiilevels.f90
ifritcalc_f90            = ifritcalc.f90
ifrit2gnuplot_f90        = ifrit2gnuplot.f90
makelexington_f90        = makelexington.f90
makeblackbody_f90        = makeblackbody.f90
r23_f90                  = r23.f90
stronglineZmaps_f90      = stronglineZmaps.f90
testauxmath_f90          = testauxmath.f90            
restrict_volume_f90      = restrict_volume.f90

data_types_obj           = data_types$(OBJEXT)
modauxmath_obj           = modauxmath$(OBJEXT)
definitions_obj          = definitions$(OBJEXT)
chemical_data_obj        = chemical_data$(OBJEXT)
natural_constants_obj    = natural_constants$(OBJEXT)
min_healpix_obj          = min_healpix$(OBJEXT)
geometry_obj             = geometry$(OBJEXT)
data_input_obj           = data_input$(OBJEXT)
absorption_obj           = absorption$(OBJEXT)
hydrogen_obj             = hydrogen$(OBJEXT)
output_obj               = output$(OBJEXT)
subsystems_obj           = subsystems$(OBJEXT)
OpenMP_analytic_obj      = OpenMP_analytic$(OBJEXT)

modblackbody_obj         = modblackbody$(OBJEXT)
modr23_obj               = modr23$(OBJEXT)
modstrongline_obj        = modstrongline$(OBJEXT)

restrict_volume_obj      = restrict_volume$(OBJEXT)

absolutetorelative       = absolutetorelative$(EXEEXT)
createtemperature        = createtemperature$(EXEEXT)
diskcut                  = diskcut$(EXEEXT)
fractalgenerator         = fractalgenerator$(EXEXT)
nebstat                  = nebstat$(EXEEXT)
radiusbench              = radiusbench$(EXEEXT)
sb2ms                    = sb2ms$(EXEEXT)
threedhydro_ext          = 3dhydro_ext$(EXEEXT)
oiilevels                = oiilevels$(EXEEXT)
ifritcalc                = ifritcalc$(EXEEXT)
ifrit2gnuplot            = ifrit2gnuplot$(EXEEXT)
makelexington            = makelexington$(EXEEXT)
makeblackbody            = makeblackbody$(EXEEXT)
r23                      = r23$(EXEEXT)
stronglineZmaps          = stronglineZmaps$(EXEEXT)
testauxmath              = testauxmath$(EXEEXT)
restrict_volume          = restrict_volume$(EXEEXT)


# The dependence to modauxmath_obj is currently artificial to make sure that 
# it is compiled. It is not used yet.
all: binary $(absolutetorelative) $(createtemperature) $(diskcut)\
            $(fractalgenerator) $(nebstat) $(radiusbench)\
            $(sb2ms) $(threedhydro_ext)  $(oiilevels) $(ifritcalc) $(ifrit2gnuplot)\
            $(makelexington) $(makeblackbody) $(r23) $(stronglineZmaps) $(testauxmath)\
	    $(restrict_volume) $(wmbasicexecs)
	$(info "Platform:")
	$(info $(PLATFORM) )
	$(info "Compiler:")
	$(COMPILER) $(VERSION)

$(data_types_obj): $(data_types_f90)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING)\
	$(data_types_obj) $(data_types_f90)

$(modauxmath_obj): $(modauxmath_f90) $(data_types_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING)\
	$(modauxmath_obj) $(modauxmath_f90)

$(definitions_obj): $(definitions_f90) $(data_types_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING)\
	$(definitions_obj) $(definitions_f90)

$(chemical_data_obj): $(chemical_data_f90) $(data_types_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING)\
	$(chemical_data_obj) $(chemical_data_f90)


$(natural_constants_obj): $(natural_constants_f90) $(data_type_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(natural_constants_obj) $(natural_constants_f90)

$(min_healpix_obj): $(min_healpix_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(min_healpix_obj) $(min_healpix_f90)

$(geometry_obj): $(geometry_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(geometry_obj) $(geometry_f90)

$(data_input_obj):  $(data_input_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(data_input_obj) $(data_input_f90)

$(absorption_obj): $(absorption_f90) $(data_types_obj) $(natural_constants_obj) $(geometry_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(absorption_obj) $(absorption_f90)

$(hydrogen_obj): $(hydrogen_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(hydrogen_obj) $(hydrogen_f90)

$(output_obj): $(output_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(output_obj) $(output_f90)

$(subsystems_obj): $(subsystems_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(subsystems_obj) $(subsystems_f90)

$(OpenMP_analytic_obj): $(OpenMP_analytic_f90) $(data_types_obj) $(natural_constants_obj) 
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(OpenMP_analytic_obj) $(OpenMP_analytic_f90)

$(modblackbody_obj): $(modblackbody_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(modblackbody_obj) $(modblackbody_f90)

$(modstrongline_obj): $(modstrongline_f90) $(data_types_obj) $(natural_constants_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(modstrongline_obj) $(modstrongline_f90)

$(Datadefs_obj): $(Datadefs_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(Datadefs_obj) $(Datadefs_f90)

$(Modgen_obj): $(Datadefs_obj) $(Modgen_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(Modgen_obj) $(Modgen_f90)

$(abund_obj): $(Modgen_obj) $(Datadefs_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(abund_obj) $(abund_f90)

$(restrict_volume_obj): $(data_types_obj) $(natural_constants_obj) $(definitions_obj) $(restrict_volume_f90)
	$(COMPILER) $(OPENMP) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NOLINK) $(NAMING) $(restrict_volume_obj) $(restrict_volume_f90)


binary: $(data_types_obj) $(definitions_obj) $(chemical_data_obj) $(natural_constants_obj)\
        $(min_healpix_obj) $(geometry_obj) $(data_input_obj) $(absorption_obj)\
        $(hydrogen_obj) $(output_obj) $(subsystems_obj) $(OpenMP_analytic_obj)
	$(LINKER) $(OPENMP) $(DEBUG) $(GPROF) $(NAMING) $(PROGRAM) $(data_types_obj)\
        $(definitions_obj) $(natural_constants_obj) $(geometry_obj) $(min_healpix_obj) \
	$(data_input_obj) $(absorption_obj) $(hydrogen_obj) $(output_obj) $(OpenMP_analytic_obj) \
	$(LAPACK) 

$(absolutetorelative): $(absolutetorelative_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(absolutetorelative) $(absolutetorelative_f90)


$(createtemperature): $(createtemperature_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(createtemperature) $(createtemperature_f90)

$(diskcut): $(createtemperature_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(diskcut) $(diskcut_f90)

$(fractalgenerator): $(fractalgenerator_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(fractalgenerator) $(fractalgenerator_f90)


$(nebstat): $(nebstat_f90) $(natural_constants_obj) $(data_types_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(nebstat) $(nebstat_f90)  $(natural_constants_obj) $(data_types_obj)

$(radiusbench): $(radiusbench_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(radiusbench) $(radiusbench_f90)

$(sb2ms): $(sb2ms_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(sb2ms) $(sb2ms_f90)

$(spectinfo): $(spectinfo_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(spectinfo)  $(spectinfo_f90)

$(threedhydro_ext): $(threedhydro_ext_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(threedhydro_ext) $(threedhydro_ext_f90)

$(oiilevels): $(oiilevels_f90) $(subsystems_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(subsystems_obj) $(LAPACK) $(NAMING) $(oiilevels)  $(oiilevels_f90) 

$(ifritcalc): $(ifritcalc_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(ifritcalc) $(ifritcalc_f90)

$(ifrit2gnuplot): $(ifrit2gnuplot_f90)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(NAMING) $(ifrit2gnuplot) $(ifrit2gnuplot_f90)

$(makelexington): $(makelexington_f90) $(modblackbody_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(modblackbody_obj)  -o $(makelexington) $(makelexington_f90)

$(makeblackbody): $(makeblackbody_f90) $(modblackbody_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(modblackbody_obj)  -o $(makeblackbody) $(makeblackbody_f90)

$(testauxmath): $(testauxmath_f90) $(modauxmath_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(modauxmath_obj)  -o $(testauxmath) $(testauxmath_f90)

$(r23): $(r23_f90) $(modstrongline_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(modstrongline_obj) $(NAMING) $(r23) $(r23_f90)

$(stronglineZmaps): $(stronglineZmaps_f90) $(modstrongline_obj)
	$(COMPILER) $(OPTIMIZE) $(DEBUG) $(GPROF) $(modstrongline_obj) $(NAMING) $(stronglineZmaps) $(stronglineZmaps_f90)

$(restrict_volume): $(restrict_volume_obj) $(definitions_obj) $(natural_constants_obj) $(data_types_obj)
	$(LINKER)  $(OPENMP) $(DEBUG) $(GPROF) $(NAMING) $(restrict_volume)\
                   $(restrict_volume_obj) $(definitions_obj) $(natural_constants_obj) $(data_types_obj)

clean:
	$(REMOVE) *$(OBJEXT)
	$(REMOVE) *$(MODEXT) 
	$(REMOVE) $(PROGRAM)$(EXEEXT) 
	$(REMOVE) $(absolutetorelative) $(createtemperature) $(diskcut) $(fractalgenerator) $(fsplit)\
                  $(ifrit2gnuplot) $(ifritcalc) $(makelexington) $(makeblackbody)\
                  $(nebstat) $(oiilevels) $(radiusbench)\
                  $(sb2ms)   $(r23) $(stronglineZmaps) $(threedhydro_ext)  $(testauxmath) $(wmbasicexecs)\
                  $(restrict_volume) 

