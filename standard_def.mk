# GNU Make 3.80

####### FLAGS
# moddir=moddir directory of module file
# FC=mpif90, ifort
# BLAS= gotoblas, gotoblas2, gotoblas2mp
# FFT=fftw3
# OMP=on,off
# STACK_CHECK=on,off
# debug=on  : to build debug version
# prof=on : to build profile information
#
################################

#### FLAGS definition

FFLAGS=
LDFLAGS=

#intel compiler
ifeq ($(COMPILER),intel)
	## serial
	ifeq ($(COMPILER_MPI),on)
		FC=mpif90
	else
		FC=ifort
	endif

	## DEFAULT
	FFLAGS = -module $(MODDIR) 

	## OMP
	ifeq ($(OMP),on)
		FFLAGS += -openmp
		ifeq ($(OMPREP),on) 
			FFLAGS += -openmp-report2
		endif
	endif

	## stack check
	ifeq ($(STACK_CHECK),on) 
		FFLAGS += -traceback -fp-stack-check
	endif	

	## debug	
	ifeq ($(debug),on)	
		FFLAGS += -g -O0
	endif
	## optimization
	ifeq ($(optimization),on)
		FFLAGS += -xW -O2
	endif

	## gprof
	ifeq ($(prof),on)
		FFLAGS += -p
	endif
endif

ifeq ($(COMPILER),pgi)
	ifeq ($(COMPILER_MPI),on)
		FC=mpif90
	else
		FC=pgf95
	endif	
	
	## DEFAULT
	FFLAGS = -module $(MODDIR) 

	ifeq ($(OMP),on)
		FFLAGS += -mp
	endif
	## optimization
	ifeq ($(optimization),on)
		FFLAGS += -O3 -tp barcelona-64 -fastsse -Mmovnt
		FFLAGS += -Mpreprocess -Mprefetch=w -Minline=reshape -fPIC
	endif	
endif

##	FFLAGS += -tp barcelona-64 -fast

# BLAS=gotoblas
ifeq ($(BLAS),gotoblas)
	FFLAGS += -I$(TACC_GOTOBLAS_DIR)
	LDFLAGS += -L$(TACC_GOTOBLAS_LIB) -lgoto_lp64
endif	

# BLAS=gotoblas_mp
ifeq ($(BLAS),gotoblas_mp)
	FFLAGS += -I$(TACC_GOTOBLAS_DIR)
	LDFLAGS += -L$(TACC_GOTOBLAS_LIB) -lgoto_lp64_mp
endif	

# BLAS=gotoblas2
ifeq ($(BLAS),gotoblas2)
	FFLAGS += -Wl,-rpath,${TACC_GOTOBLAS2_LIB}
	LDFLAGS += -L$(TACC_GOTOBLAS_LIB) -lgoto_lp64
endif

# BLAS=gotoblas2mp 
# this link gotoblas2 openmp version
ifeq ($(BLAS),gotoblas2mp)
	FFLAGS += -Wl,-rpath,${TACC_GOTOBLAS2_LIB}
	LDFLAGS += -L$(TACC_GOTOBLAS_LIB) -lgoto_lp64_openmp
endif

ifeq ($(FFT), fftw3)
	FFLAGS += -I$(TACC_FFTW3_INC) 
	LDFLAGS += -L$(TACC_FFTW3_LIB) -lfftw3	
endif	

ifeq ($(tau),on)
	FC=tau_f90.sh	
	FFLAGS += -optTauSelectFile=*
endif
######################################################################## RULES
.SUFFIXES:
.SUFFIXES: .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<
