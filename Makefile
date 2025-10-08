# Makefile for INCA Fortran project

.DEFAULT_GOAL := roda.exe

SRC_DIR = src
OBJ_DIR = obj

FC = gfortran

FFLAGS = -Wall -g -msse4.2 -fopenmp -fcheck=all \
         -Waliasing -Wampersand -Wconversion -Wsurprising \
         -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant \
         -J$(OBJ_DIR)

LIB = -llapack -lblas

SRCF90 = $(wildcard $(SRC_DIR)/*.f90)
SRCF   = $(wildcard $(SRC_DIR)/*.f)
OBJF90 = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCF90))
OBJF   = $(patsubst $(SRC_DIR)/%.f,$(OBJ_DIR)/%.o,$(SRCF))
OBJ    = $(OBJF90) $(OBJF)

$(shell mkdir -p $(OBJ_DIR))

# Pattern rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_DIR) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) -I$(OBJ_DIR) -c $< -o $@

# Final executable
roda.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIB)

# Generate dependencies
deps:
	python3 generate_deps.py

# Clean
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod roda.exe test.log

# Test
test:
	./test.sh 

# Dependencies (generated)
include deps.mk


