include ../../arch.inc

SRC_OBJ=cos_pchip.o dchfev.o  dpchfe.o  dpchim.o  dpchst.o fdump.o  i1mach.o  j4save.o  xercnt.o  xerhlt.o  xermsg.o  xerprn.o  xersve.o  xgetua.o

LIBS_LOC=-L. -lpchip 

all: libpchip.a test tpchip.$(EXT)

%.o:%.f90
	$(F90) -c $(F90FLAGS) $<

%.o:%.f
	$(F90) -c $(F90FLAGS) $<	

libpchip.a:$(SRC_OBJ)
	ar -cr libpchip.a $(SRC_OBJ)

tpchip.$(EXT): tpchip.f libpchip.a
	$(MEX) $^

test: test.f90 libpchip.a
	$(F90) $(F90FLAGS) -o test test.f90 $(LIBS_LOC) $(LIBS)

clean:
	rm -f $(SRC_OBJ) test libpchip.a tpchip.$(EXT)
