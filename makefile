
.PHONY : clean all options

all : 
	@echo "---------------------------------"
	@echo "Compiling PhiMATSFEM ... "
	@echo "Start: " $(shell date)
	@echo "---------------------------------"
	make -C src/
	@echo "---------------------------------"
	@echo "End: " $(shell date)
	@echo "Compiling PhiMATSFEM was successful "
	@echo "---------------------------------"

options: 
	@echo " "
	@echo "Options for make SETTIGNS="
	@echo " (empty) -> No flags"
	@echo "  DEBUG  -> -g -Wall -DDEBUG"
	@echo "  OPTIM  -> -O3"
	@echo "  PETSC  -> with petsc compiler and linker"
	
clean : 
	rm -rf src/obj
