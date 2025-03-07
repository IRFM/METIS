include ../../../arch.inc
include make.inc

#MUMPS_DVER = DROIDE_DBL_SEQ
#MUMPS_ZVER = DROIDE_cmplx16_SEQ


INCMUMPS = -I../include 
LIBMUMPS = -L../lib -ldmumps -lzmumps 

INCTOT = $(IORDERINGS) $(INCSEQ) $(INCMUMPS)
LIBTOT = $(LIBMUMPS) $(LORDERINGS) $(LIBSEQ) $(LIBFORT) -lpthread

all: dmumpsmex.$(EXT)

dmumpsmex.$(EXT): dmumpsmex.c
ifeq ($(MATLABVERSION),-DMATLAB86)
	$(MEX) $(OPTC) -largeArrayDims -c $(INCTOT) -I$(MATLABHOME)/extern/include dmumpsmex.c 
	$(MEX_F90) -largeArrayDims -fortran $(OPTC)  dmumpsmex.o $(INCTOT) $(LIBTOT)
else ifeq ($(MATLABVERSION),-DMATLAB73)
	$(MEX) $(OPTC) -largeArrayDims -c $(INCTOT) -I$(MATLABHOME)/extern/include dmumpsmex.c 
	$(MEX_F90) -largeArrayDims -fortran $(OPTC) -o dmumpsmex.$(EXT) dmumpsmex.o $(INCTOT) $(LIBTOT)
else
	$(MEX) $(OPTC) -c $(INCTOT) -I$(MATLABHOME)/extern/include dmumpsmex.c
	$(MEX_F90) -fortran $(OPTC) -o dmumpsmex.$(EXT) dmumpsmex.o $(INCTOT) $(LIBTOT)
endif
clean:
	rm -f *.o *.$(EXT)
