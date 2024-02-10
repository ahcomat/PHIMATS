# Include makefile definitions
include ../MakeFile.defs

# include dir
INCLUDES = ../include

# Objects dir
OBJDIR := ../obj

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}

# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

# all: LinearElastic.o	
# %.o: %.cxx
# 	$(COMPILE.cc) $(CPPFLAGS) $< -I ${PETSC_DIR}/include/ -I $(H5ID) -L $(H5LD) -lhdf5 -lhdf5_cpp $(PETSC_LIB)

all: ${OBJECTS}

${OBJDIR}/%.o: %.cxx
	@echo ${OBJDIR}
	$(COMPILE.cc) ${$(SETTINGS)} $< -I ${PETSC_DIR}/include/ -I $(EIGEN) -I $(INCLUDES) -L ${PETSC_LIB}
	mv *.o ../obj
	

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
