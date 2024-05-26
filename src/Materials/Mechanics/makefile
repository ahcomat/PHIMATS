# Include makefile definitions
include ../../MakeFile.defs

# include dir
INCLUDES = ../../include

# Objects dir
OBJDIR := ../../obj

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}

# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

all: ${OBJECTS}
	# @echo "${OBJECTS}"

../../obj/%.o: %.cxx
	${CXX} $(CFLAGS) -c $< -o $@ -I $(INCLUDES) -I $(EIGEN) -I $(H5ID) -L $(H5LD) -lhdf5 -lhdf5_cpp
