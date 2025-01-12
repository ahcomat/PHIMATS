# Include makefile definitions
include MakeFile.defs

# Objects dir
OBJDIR := obj

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}

# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

export EXCLUDE = obj/ include/ FEM_utils/ #Solvers/ FiniteElements/ Materials/ 

SUBDIRS := $(wildcard */)
SUBDIRS :=$(sort $(dir $(wildcard */)))
# Exlude entries from list above
SUBDIRS := $(filter-out $(EXCLUDE), $(SUBDIRS))

all: $(SUBDIRS) ${OBJECTS}

obj/%.o: %.cxx
	# ${CXX} $(CPPFLAGS) -c $< -o $@ -I $(H5ID) -I $(EIGEN) -L $(H5LD) -lhdf5 -lhdf5_cpp
	$(COMPILE.cc) ${$(SETTINGS)} $< -I ${PETSC_DIR}/include/ -I $(EIGEN) -I $(PHIMATSINCLUDES) -L ${PETSC_LIB} -I $(H5ID) -L $(H5LD) -lhdf5 -lhdf5_cpp
	mv *.o obj

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: all $(SUBDIRS)
	
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


