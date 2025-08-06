#=======================================================================
# Makefile for calliope modified by Y. Kawazura @ Sep 19 2019
# Basically, don't change this file!
#=======================================================================
include ../../Makefile.in
include ../../arch/$(arch).in

MODEL_DIR = model/$(MODEL)

# F90FLAGS = -cpp

#vvvvvvv  2decomp option setting  vvvvvvv#
DECOMP_DIR = ./2decomp
DECOMP_FLGS = -DDOUBLE_PREC
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

### List of files for 2DECOMP&FFT
SRC = $(DECOMP_DIR)/decomp_2d.f90 $(DECOMP_DIR)/glassman.f90 $(DECOMP_DIR)/fft_fftw3.f90 $(DECOMP_DIR)/io.f90 expand.f90

#######OPTIONS settings###########
OBJ = $(SRC:.f90=.o)


all: expand

expand : $(OBJ)
	$(FC) $(F90FLAGS) -o ../../$@ $(OBJ) $(FFTW3_LIB)

%.o : %.f90
	$(FC) $(F90FLAGS) $(DECOMP_FLGS) $(FFTW3_INC) -I$(DECOMP_DIR) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.mod $(DECOMP_DIR)/*.o $(DECOMP_DIR)/*.mod expand

.PHONY: print
print:
	@echo DEFS is $(DEFS)
	@echo OBJ is $(OBJ)
	@echo FFT is $(FFT)
