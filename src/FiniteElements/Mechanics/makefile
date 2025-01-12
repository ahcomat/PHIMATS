# Include makefile definitions
include ../../MakeFile.defs

# Objects dir
OBJDIR := ../../obj

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}

# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

all: ${OBJECTS}

all: ${OBJECTS}
	@echo "${PHIMATSINCLUDES}"

../../obj/%.o: %.cxx
	${CXX} $(CPPFLAGS) -c $< -o $@ -I $(PHIMATSINCLUDES) -I $(EIGEN) -I $(H5ID) -L $(H5LD) -lhdf5 -lhdf5_cpp