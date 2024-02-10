# Include makefile definitions
include MakeFile.defs

# include dir
INCLUDES = include

# Objects dir
OBJDIR := obj

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}

# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

export EXCLUDE = obj/ include/ FEM_utils/ LinearAlgebra/ #FiniteElements/ Materials/ # Solvers/

SUBDIRS := $(wildcard */)
SUBDIRS :=$(sort $(dir $(wildcard */)))
# Exlude entries from list above
SUBDIRS := $(filter-out $(EXCLUDE), $(SUBDIRS))

all: $(SUBDIRS) ${OBJECTS}

obj/%.o: %.cxx
	${CXX} $(CPPFLAGS) -c $< -o $@ -I $(INCLUDES) -I $(H5ID) -I $(EIGEN) -L $(H5LD) -lhdf5 -lhdf5_cpp

$(SUBDIRS):
	mkdir -p obj
	$(MAKE) -C $@

.PHONY: all $(SUBDIRS)
