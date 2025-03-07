include ../arch.inc

all: mexpde1dsolver.$(EXT) interpos_old.$(EXT) interpos.$(EXT)

MUMPSHOME=../libs/MUMPS_4.7.3
SKITHOME=../libs/SPARSKIT2

mexpde1dsolver.$(EXT):  mexpde1dsolver.o pde1solver.o 
	@#$(MEX) -fortran $^ $(MUMPSHOME)/lib/libdmumps.a $(MUMPSHOME)/PORD/lib/libpord.a $(MUMPSHOME)/libseq/libmpiseq.a $(SKITHOME)/libskit.a
ifeq ($(MATLABVERSION),-DMATLAB86)
	$(MEX_F90)  $^ $(MUMPSHOME)/lib/libdmumps.a $(MUMPSHOME)/PORD/lib/libpord.a $(MUMPSHOME)/../metis/libmetis.a  $(MUMPSHOME)/libseq/libmpiseq.a
else
	$(MEX_F90) -fortran $^ $(MUMPSHOME)/lib/libdmumps.a $(MUMPSHOME)/PORD/lib/libpord.a $(MUMPSHOME)/../metis/libmetis.a  $(MUMPSHOME)/libseq/libmpiseq.a
endif

interpos_old.$(EXT): interpos_old.f90 $(includefile)
	$(MEX) $<
	
interpos.$(EXT): interpos.f90
	mex $<

pde1solver.o: pde1solver.f90
	$(F90) -c $(FFLAGS_MUMPS) -I$(MUMPSHOME)/libseq -I$(MUMPSHOME)/include $<

interpos_old.o: interpos_old.f90
	$(MEX) -c $<
	
mexpde1dsolver.o: mexpde1dsolver.c
	$(MEX) -c $<


.PHONY: clean 

clean:
	rm -f *.o *.mod *.$(EXT)

