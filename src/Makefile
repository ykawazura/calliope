#=======================================================================
# Makefile for calliope modified by Y. Kawazura @ Feb 15 2020
# Basically, don't change this file!
#=======================================================================
include ../Makefile.in
include ../arch/$(arch).in

MODEL_DIR = model/$(MODEL)

INCLUDES = $(P3DFFT_INC) $(NETCDF_INC)
LIBS += $(MPI_LIB) $(P3DFFT_LIB) $(FFTW3_LIB) $(NETCDF_LIB)

F90FLAGS += -DMEASURE -DSTRIDE1 -DFFTW 
LDFLAGS += -DMEASURE -DSTRIDE1

ifeq ($(DEBG),yes)
	F90FLAGS += -DDEBG
endif

#vvvvvvv      openmp setting      vvvvvvv#
ifeq ($(USE_OPENMP),yes)
	LIBS += $(LIB_OPENMP)
	F90FLAGS += $(F90FLAGS_OPENMP) -DOPENMP
	LDFLAGS += $(F90FLAGS_OPENMP) -DOPENMP
endif
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

### List of files for 'RMHD'
SRC_RMHD = \
			time_stamp.F90 file.F90 \
			$(MODEL_DIR)/params.F90 mp.F90 mpiio.F90 \
			terminate.F90 grid.F90 time.F90 utils.F90 dealias.F90 \
			shearing_box.F90 \
			$(MODEL_DIR)/fields.F90 \
			$(MODEL_DIR)/force.F90 \
			$(MODEL_DIR)/$(MODEL)_specific.F90 \
			$(MODEL_DIR)/io.F90 $(MODEL_DIR)/diagnostics.F90 \
			$(MODEL_DIR)/advance.F90 calliope.F90

### List of files for 'RRMHD'
SRC_RRMHD = $(SRC_RMHD)

### List of files for 'MHD_INCOMP'
SRC_MHD_INCOMP = \
			time_stamp.F90 file.F90 \
			$(MODEL_DIR)/params.F90 mp.F90 mpiio.F90 \
			terminate.F90 grid.F90 time.F90 utils.F90 dealias.F90 \
			shearing_box.F90 \
			$(MODEL_DIR)/fields.F90 \
			$(MODEL_DIR)/force.F90 \
			$(MODEL_DIR)/$(MODEL)_specific.F90 \
			$(MODEL_DIR)/io.F90 $(MODEL_DIR)/diagnostics.F90 \
			$(MODEL_DIR)/advance.F90 calliope.F90

### List of files for 'MHD_COMP_ISOTH'
SRC_MHD_COMP_ISOTH = \
			time_stamp.F90 file.F90 \
			$(MODEL_DIR)/params.F90 mp.F90 mpiio.F90 \
			terminate.F90 grid.F90 time.F90 utils.F90 dealias.F90 \
			shearing_box.F90 \
			$(MODEL_DIR)/fields.F90 \
			$(MODEL_DIR)/force.F90 \
			$(MODEL_DIR)/$(MODEL)_specific.F90 \
			$(MODEL_DIR)/io.F90 $(MODEL_DIR)/diagnostics.F90 \
			$(MODEL_DIR)/advance.F90 calliope.F90

### List of files for 'MHD_COMP_ADIAB'
SRC_MHD_COMP_ADIAB = \
			time_stamp.F90 file.F90 \
			$(MODEL_DIR)/params.F90 mp.F90 mpiio.F90 \
			terminate.F90 grid.F90 time.F90 utils.F90 dealias.F90 \
			shearing_box.F90 \
			$(MODEL_DIR)/fields.F90 \
			$(MODEL_DIR)/force.F90 \
			$(MODEL_DIR)/$(MODEL)_specific.F90 \
			$(MODEL_DIR)/io.F90 $(MODEL_DIR)/diagnostics.F90 \
			$(MODEL_DIR)/advance.F90 calliope.F90

### List of files for 'EMHD'
SRC_EMHD = \
			time_stamp.F90 file.F90 \
			$(MODEL_DIR)/params.F90 mp.F90 mpiio.F90 \
			terminate.F90 grid.F90 time.F90 utils.F90 dealias.F90 \
			shearing_box.F90 \
			$(MODEL_DIR)/fields.F90 \
			$(MODEL_DIR)/force.F90 \
			$(MODEL_DIR)/$(MODEL)_specific.F90 \
			$(MODEL_DIR)/io.F90 $(MODEL_DIR)/diagnostics.F90 \
			$(MODEL_DIR)/advance.F90 calliope.F90

OBJ = $(SRC_$(MODEL):.F90=.o)

calliope : $(OBJ)
	$(FC) $(LDFLAGS) -o ../$@ $^ $(LIBS)

$(TARGET) : $(TARGET).o
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o : %.F90
	$(FC) $(INCLUDES) $(F90FLAGS) -D$(MODEL) -D_MODEL_=\"$(MODEL)\" -c $< -o $@


.PHONY: print clean

print:
	@echo OBJ = $(OBJ)
	@echo LIBS = $(LIBS)
	@echo F90FLAGS = $(F90FLAGS)
	@echo LDFLAGS = $(LDFLAGS)

clean:
	rm -f *.o *.mod $(MODEL_DIR)/*.o $(MODEL_DIR)/*.mod calliope
