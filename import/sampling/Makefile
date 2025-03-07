#
# We introduced additional dependencies in order to make the Makefile run in parallel
#

include ../../arch.inc

OBJ=ctable1.o cround.o cspline.o

all: tsplinet.$(EXT) tsample.$(EXT) pchipslopest_mex.$(EXT)

tsplinet.$(EXT): tsample.$(EXT) tsplinet.c
	$(MEX) tsplinet.c -output tsplinet.$(EXT)

tsample.$(EXT): $(OBJ) tsample.c
	$(MEX) $^ -output tsample.$(EXT)

ctable1.o:ctable1.c
	$(CC) -c $(CFLAGS_SAMPLING) $^

cround.o:cround.c
	$(CC) -c $(CFLAGS_SAMPLING) $^

cspline.o:cspline.c
	$(CC) -c $(CFLAGS_SAMPLING) $^

pchipslopest_mex.$(EXT): pchipslopest_mex.c
	$(MEX) $^ -output pchipslopest_mex.$(EXT)


.PHONY: clean

clean:
	rm -rf *.$(EXT) $(OBJ)


