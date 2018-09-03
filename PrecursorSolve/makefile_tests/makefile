######################
#  Make file
#####################

SHELL = /bin/sh

# Compiler options 
F90 = gfortran
FCLAGS = -g -c -Wall -fbounds-check -fdefault-real-8
FLFLAGS = -g
LIBDIR = ../lapack-3.8.0/ 
LIB = -llapack -lblas

# Directory containing the main makefile, used by "VPATH":
MAKEFILE_DIR = .

# Write paths relative to MAKEFILE_DIR, separated by blanks (do not
# include ${MAKEFILE_DIR} in the list):
VPATH := ${MAKEFILE_DIR} $(addprefix ${MAKEFILE_DIR}/src,bin)

# Finds all source files
SRCS = $(wildcard ${MAKEFILE_DIR}/src/*.f90)
OBJS=$(patsubst %.f90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard ${MAKEFILE_DIR}/src/*_M.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

PROGRAM = msr1d
PRG_OBJ = ${PROGRAM}.o

# Clean the suffixes
.SUFFIXES:

# Set the suffixes we are interested in
.SUFFIXES: .f90 .o

# make without parameters will make first target found.
default : $(PROGRAM)

${MAKEFILE_DIR}/src/parameters_fe_M.o : ${MAKEFILE_DIR}/src/parameters_fe_M.f90
	$(F90) $(FCLAGS) -c ${MAKEFILE_DIR}/src/parameters_fe_M.f90

#$(MOD_OBJS) : %.o : %.f90
#	$(F90) $(FCFLAGS) -o $@ $< 

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	$(F90) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(F90) $(FLFLAGS) -o $@ $^

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Print out some stuff to make sure things are cool
debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

# Clean up the directory
clean: 
	rm *.mod *.o

.PHONY. : debug clean
