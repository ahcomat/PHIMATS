
.PHONY : clean all options

# Compiler
export CXX = mpic++

OPTIM = -O3 -std=c++20 
DEBUG = -O0 -g -Wall -DDEBUG -std=c++20 

SETTINGS ?= OPTIM

ifeq ($(SETTINGS), DEBUG)
	CXXFLAGS = $(DEBUG)
else ifeq ($(SETTINGS), OPTIM)
	CXXFLAGS =  $(OPTIM)
endif

all : 
	@echo ""
	@echo "---------------------------------"
	@echo "Compiling PHIMATS ... "
	@echo "Start: " $(shell date)
	@echo "---------------------------------"
	@echo "Settings:"
	@echo "  Compiler:         $(CXX)"
	@echo "  Compilation Mode: $(SETTINGS)"
	@echo "  CXXFLAGS:         $(CXXFLAGS)"
	@echo "---------------------------------"
	@echo ""
	@echo "Create obj directory"
	mkdir -p src/obj
	@echo "Create lib directory"
	mkdir -p src/lib
	@echo ""
	make -C src/
	@echo ""
	@echo "---------------------------------"
	@echo "End: " $(shell date)
	@echo "Compiling PHIMATS was successful "
	@echo "---------------------------------"
	@echo ""

options: 
	@echo " "
	@echo "Options for make SETTIGNS="
	@echo "  	(empty) -> OPTIM"
	@echo "  	OPTIM   -> -O3 -std=c++20 "
	@echo "  	DEBUG   -> -O0 -g -Wall -DDEBUG -std=c++20 "
	@echo " "
	
clean : 
	rm -rf src/obj src/lib