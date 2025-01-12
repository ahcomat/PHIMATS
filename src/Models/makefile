# Include makefile definitions
include ../MakeFile.defs

# Objects dir
OBJDIR := ../obj

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}

# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

all: ${OBJECTS}

# in PETSc we need "mv *.o ../obj"
${OBJDIR}/%.o: %.cxx
	#@echo ${OBJDIR}
	$(COMPILE.cc) ${$(SETTINGS)} $< -I ${PETSC_DIR}/include/ -I $(EIGEN) -I $(PHIMATSINCLUDES) -L ${PETSC_LIB} -I $(H5ID) -L $(H5LD) -lhdf5 -lhdf5_cpp
	mv *.o ../obj
	
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules