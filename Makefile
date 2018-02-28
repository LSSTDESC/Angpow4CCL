

ifndef MACH
	MACH := $(shell uname)
endif

ifeq ($(MACH), Linux)
#==== Linux system
# choose compiler
	include Linux_g++_make.inc
#	include Linux_icpc_make.inc
#if unusal loction specify below (here is for nersc)
#	FFTWLIB   = -L$(FFTW_DIR)
#	FFTWINC   = -I$(FFTW_INC)
else
#==== Darwin system
# for clang needs OpenMP support tesed with clang 3.9
#	include Darwin_clang_make.inc
	include Darwin_g++_omp_make.inc
	FFTWDIR  = /opt/local
	FFTWLIB   = -L$(FFTWDIR)/lib
	FFTWINC   = -I$(FFTWDIR)/include
endif

# ==== Profiling ====
PROFIL=
ifdef PROFILING
	PROFIL = -DPROFILING
endif


FFTWLIBN  = $(FFTWLIB) -lfftw3 -lfftw3_threads
#FFTWLIBN = -L$(FFTWLIB) -lfftw3 -lfftw3_omp



OBJ = ./Objs/
EXE = ./bin/
INC = ./inc/
INCEXT = $(INC)Angpow/
SRC = ./src/
LIB = ./lib/


BOOSTDIR = $(INC)
BOOSTINC = -I$(BOOSTDIR)


CPPFLAGS  += $(PROFIL) $(BOOSTINC) $(FFTWINC)
LDFLAGS   += $(FFTWLIBN) -lm


#  Define our target list
.PHONY: default
default: makedir lib angpow angcor magnification

.PHONY: all 
all : makedir lib angpow

.PHONY: tidy
tidy : 
	find . -name "*~" | xargs -I {} rm {}

.PHONY: clean
clean :
	rm -rf ./$(OBJ)/* ./$(EXE)/* ./$(LIB)/*

.PHONY: makedir
makedir :
	@mkdir -p $(OBJ)
	@mkdir -p $(EXE)
	@mkdir -p $(LIB)


#C++ common Objects
CXXOBJ = $(OBJ)walltimer.o \
	$(OBJ)walltime_c.o \
	$(OBJ)angpow_bessel.o \
	$(OBJ)angpow_chebyshevInt.o \
	$(OBJ)angpow_radint.o \
	$(OBJ)angpow_kinteg.o \
	$(OBJ)angpow_pk2cl.o \
	$(OBJ)angpow_tools.o \
	$(OBJ)angpow_clbase.o \
	$(OBJ)angpow_powspec.o \
	$(OBJ)angpow_parameters.o \
	$(OBJ)angpow_cosmo.o \
	$(OBJ)angpow_ccl.o \
	$(OBJ)angpow_ctheta.o 
	#$(OBJ)angpow_growth.o


CXXSHOBJ = walltimer.o \
	walltime_c.o \
	angpow_bessel.o \
	angpow_chebyshevInt.o \
	angpow_radint.o \
	angpow_kinteg.o \
	angpow_pk2cl.o \
	angpow_tools.o \
	angpow_clbase.o \
	angpow_powspec.o \
	angpow_parameters.o \
	angpow_cosmo.o \
	angpow_ccl.o \
	angpow_ctheta.o \
	angpow_growth.o


#C++ common Headers
CXXHDR =  $(INCEXT)walltimer.h \
	$(INCEXT)walltime_c.h \
	$(INCEXT)angpow_quadinteg.h \
	$(INCEXT)angpow_numbers.h \
	$(INCEXT)angpow_exceptions.h \
	$(INCEXT)angpow_func.h  \
	$(INCEXT)angpow_bessel.h \
	$(INCEXT)angpow_fft.h \
	$(INCEXT)angpow_chebyshevInt.h\
	$(INCEXT)angpow_integrand_base.h \
	$(INCEXT)angpow_integrand.h \
	$(INCEXT)angpow_powspec_base.h \
	$(INCEXT)angpow_powspec.h \
	$(INCEXT)angpow_radial.h \
	$(INCEXT)angpow_radial_base.h \
	$(INCEXT)angpow_radint.h \
	$(INCEXT)angpow_kinteg.h \
	$(INCEXT)angpow_utils.h \
	$(INCEXT)angpow_cosmo_base.h \
	$(INCEXT)angpow_cosmo.h \
	$(INCEXT)angpow_pk2cl.h \
	$(INCEXT)angpow_tools.h \
	$(INCEXT)angpow_parameters.h \
	$(INCEXT)angpow_clbase.h \
	$(INCEXT)angpow_ctheta.h 
	#$(INCEXT)angpow_growth_base.h \
	#$(INCEXT)angpow_growth.h 
	#$(INCEXT)angpow_filter_base.h \
	#$(INCEXT)angpow_filter.h
	#$(INCEXT)angpow_limber.h \


#C++ rule for compiling
$(OBJ)%.o: $(SRC)%.cc $(CXXHDR)
	echo "compile... $@"
	$(CXXCOMPILE) -I$(INC) -c $< -o $@

######################
.PHONY: sharelib
sharelib : $(CXXOBJ)
	echo "make shared library..."
	cd $(OBJ); \
	$(CMDSHLCXX) -o ../$(LIB)libangpow.$(SLEXT) $(CXXSHOBJ) $(LDFLAGS)

.PHONY: lib
lib : $(LIB)libangpow.a

$(LIB)libangpow.a : $(CXXOBJ)
	$(AR) $(ARFLAGS) $@ $^

###############
.PHONY:  angpow
angpow: $(EXE)angpow
	echo '---- angpow made'

$(OBJ)angpow.o: angpow.cc $(CXXHDR)
	echo "compile... $@"
	$(CXXCOMPILE)  -I$(INC) -c $< -o $@ 

$(EXE)angpow :  $(OBJ)angpow.o $(LIB)libangpow.a
	echo "Link..."
	$(CXXLINK) -o $@ $(OBJ)angpow.o -L$(LIB) -langpow  $(LDFLAGS)

###############
.PHONY:  angcor
angcor: $(EXE)angcor
	echo '---- angcor made'

$(OBJ)angcor.o: angcor.cc $(CXXHDR)
	echo "compile... $@"
	$(CXXCOMPILE)  -I$(INC) -c $< -o $@ 

$(EXE)angcor :  $(OBJ)angcor.o $(LIB)libangpow.a
	echo "Link..."
	$(CXXLINK) -o $@ $(OBJ)angcor.o -L$(LIB) -langpow  $(LDFLAGS)

###############
.PHONY:  magnification
magnification: $(EXE)magnification
	echo '---- magnification made'

$(OBJ)magnification.o: magnification.cc $(CXXHDR)
	echo "compile... $@"
	$(CXXCOMPILE)  -I$(INC) -c $< -o $@ 

$(EXE)magnification :  $(OBJ)magnification.o $(LIB)libangpow.a
	echo "Link..."
	$(CXXLINK) -o $@ $(OBJ)magnification.o -L$(LIB) -langpow  $(LDFLAGS)

###############
.PHONY: debug
debug:
	echo "Test with bench1 init file"
	$(EXE)angpow angpow_bench1.ini

###############
.PHONY: fullcheck
fullcheck:
	echo "Test with the whole set of init files"
	./verif.sh




