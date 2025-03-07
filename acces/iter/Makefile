include ../../arch.inc

all: separatrice.$(EXT)

separatrice.$(EXT): separatrice.f dasinh.f initds.f dcsevl.f d1mach.f
	$(MEX) $^ 

clean:
	rm -rf *.$(EXT)

