# Include makefile definitions
include MakeFile.defs

# Objects dir
OBJDIR := obj
LIBDIR := lib

# Finds .cxx files in current folder
SOURCES:=${wildcard *.cxx}
# Name of objects
OBJECTS:=${patsubst %.cxx,${OBJDIR}/%.o,${SOURCES}}

export EXCLUDE = obj/ lib/ include/ FEM_utils/ #Solvers/ FiniteElements/ Materials/ 

SUBDIRS := $(wildcard */)
SUBDIRS :=$(sort $(dir $(wildcard */)))
# Exlude entries from list above
SUBDIRS := $(filter-out $(EXCLUDE), $(SUBDIRS))


# Shared library
PHIMATS_LIB := $(LIBDIR)/libphimats.so

all: $(SUBDIRS) ${OBJECTS} $(PHIMATS_LIB)

obj/%.o: %.cxx
	$(COMPILE.cc) $(CXXFLAGS) $< -I${PETSC_DIR}/include/ -I$(EIGEN) -I$(PHIMATSINCLUDES) -L${PETSC_LIB} -I$(H5ID) $(LDFLAGS)
	mv *.o obj

$(SUBDIRS):
	make -C $@

# Create library
$(PHIMATS_LIB): $(OBJECTS)
	# @echo $(wildcard $(OBJDIR)/*.o)  # Debug: List all object files in obj/
	#ar rcs $@ $(wildcard $(OBJDIR)/*.o)
	$(CXX) -shared -fPIC -o $@ $(wildcard $(OBJDIR)/*.o) $(LDFLAGS) ${PETSC_LIB}


.PHONY: all $(SUBDIRS) $(PHIMATS_LIB)
	
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


