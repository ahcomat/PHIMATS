# Objects dir
OBJDIR := ../../obj

# Include makefile definitions
include ../../MakeFile.defs

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}

# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

all: ${OBJECTS}

../../obj/%.o: %.cxx
	${CXX} $(CXXFLAGS) -c $< -o $@ -I$(PHIMATSINCLUDES) -I$(EIGEN) -I $(H5ID) $(LDFLAGS)