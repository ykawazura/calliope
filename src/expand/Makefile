#=======================================================================
# Makefile for calliope modified by Y. Kawazura @ Sep 19 2019
# Basically, don't change this file!
#=======================================================================
include ../../Makefile.in
include ../../arch/$(arch).in

MODEL_DIR = model/$(MODEL)

FFLAGS += -cpp

#vvvvvvv  2decomp option setting  vvvvvvv#
DECOMP_DIR = ./2decomp

ifeq ($(DECOMP_DOUBLE_PREC),yes)
	DECOMP_FLGS += -DDOUBLE_PREC
endif
ifeq ($(DECOMP_SAVE_SINGLE),yes)
	DECOMP_FLGS += -DSAVE_SINGLE
endif
ifeq ($(DECOMP_DEBG),yes)
	DECOMP_FLGS += -DDEBG
endif
ifeq ($(DECOMP_OVERWRITE),yes)
	DECOMP_FLGS += -DOVERWRITE
endif
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

### List of files for 2DECOMP&FFT
SRC = $(DECOMP_DIR)/decomp_2d.f90 $(DECOMP_DIR)/glassman.f90 $(DECOMP_DIR)/fft_fftw3.f90 $(DECOMP_DIR)/io.f90 expand.f90

#######OPTIONS settings###########
OBJ = $(SRC:.f90=.o)


all: expand

expand : $(OBJ)
	$(FC) $(FFLAGS) -o ../../$@ $(OBJ) $(FFTW3_LIB)

%.o : %.f90
	$(FC) $(FFLAGS) $(DECOMP_FLGS) $(FFTW3_INC) -I$(DECOMP_DIR) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.mod $(DECOMP_DIR)/*.o $(DECOMP_DIR)/*.mod expand

.PHONY: print
print:
	@echo DEFS is $(DEFS)
	@echo OBJ is $(OBJ)
	@echo FFT is $(FFT)
