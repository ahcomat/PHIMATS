# Include shared definitions and PETSc configurations
include ../MakeFile.defs

# Compile the driver program
$(PROG): $(OBJ)
	@echo ""
	@echo "Linking $(PROG) with libphimats.a"
	@echo ""
	${CXXLINKER} ${$(SETTINGS)} -o $@ ${OBJ} -L$(PHIMATSLIBDIR) -lphimats -I${PETSC_DIR}/include/ -L${PETSC_LIB} -I$(H5ID) $(LDFLAGS)
	${RM} ${OBJ}

${OBJ}: ${SRC}
	@echo $(PHIMATSINCLUDES)
	${COMPILE.cc} ${$(SETTINGS)} $< -o $@ -I$(EIGEN) -I${PETSC_DIR}/include/ -I$(PHIMATSINCLUDES) -L ${PETSC_LIB} -I$(H5ID) $(LDFLAGS)

# Clean target
.PHONY: clean2
clean2:
	rm -f $(PROG)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules