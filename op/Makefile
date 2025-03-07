include ../arch.inc

all: zvalfit.$(EXT) zinterp2d.$(EXT) zinterp2d_hel.$(EXT) zfittime.$(EXT) zfitdata1.$(EXT) zfitdata.$(EXT) lignemex.$(EXT) lignemex2.$(EXT) ellipfastmex.$(EXT)

zvalfit.$(EXT): zvalfit.f
	$(MEX) $^

zinterp2d.$(EXT): zinterp2d.f
	$(MEX) $^

zinterp2d_hel.$(EXT): zinterp2d_hel.f
	$(MEX) $^

zfittime.$(EXT): zfittime.f
	$(MEX) $^

zfitdata1.$(EXT): zfitdata1.f
	$(MEX) $^

zfitdata.$(EXT): zfitdata.f
	$(MEX) $^
lignemex.$(EXT): lignemex.f
	 $(MEX) $^
lignemex2.$(EXT): lignemex2.f
	$(MEX) $^

ellipfastmex.$(EXT): ellipfastmex.c
	$(MEX) $^


.PHONY: clean

clean:
	rm -rf *.$(EXT) 

