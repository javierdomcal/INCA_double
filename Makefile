#
#Makefile for RODA program, find a better name (;
#
ODIR=obj
FC=gfortran
FFLAGS= -Wall -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
##FFLAGS= -o3
LIB= -llapack -lblas
SRCF90=$(wildcard *.f90)
SRCF=$(wildcard *.f)

#OBJ=$(SRCF90:.f90=.o) $(SRCF:.f=.o)
OBJ=$(patsubst %.f90,$(ODIR)/%.o,$(SRCF90)) $(patsubst %.f,$(ODIR)/%.o,$(SRCF))


$(ODIR)/%.o: %.f90
	$(FC) $(FFLAGS) -o $(ODIR)/mod5.o -c mod5.f90
	$(FC) $(FFLAGS) -o $(ODIR)/intrastuff.o -c intrastuff.f90
	$(FC) $(FFLAGS) -o $@ -c $< 

$(ODIR)/%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $<


roda.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIB)
	rm -f *.mod

clean:
	@rm -f *.mod $(ODIR)/*.o roda	

#Testing
test: 
	./test.sh >& test.log


