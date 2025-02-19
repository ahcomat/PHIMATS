# Objects dir
OBJDIR := ../obj

# Include makefile definitions
include ../MakeFile.defs

# Finds .cxx files in current folder
SOURCES := ${wildcard *.cxx}

# Name of objects
OBJECTS := ${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

all: ${OBJECTS}

${OBJDIR}/%.o: %.cxx
	@echo ${OBJDIR}
	$(COMPILE.cc) $(CXXFLAGS) $< -I ${PETSC_DIR}/include/ -I$(EIGEN) -I$(PHIMATSINCLUDES) -L${PETSC_LIB}
	mv *.o ../obj
	

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
