# This file is part of MUMPS VERSION 4.7.3
# This Version was built on Fri May 4 15:54:01 2007
#
#
default: double cmplx16 matlab

#all: double simple cmplx cmplx16 matlab
all: double cmplx16 matlab

# Is Makefile.inc available ?
Makefile.inc:
	@echo "##################################################################"
	@echo BEFORE COMPILING MUMPS, YOU SHOULD HAVE AN APPROPRIATE FILE
	@echo Makefile.inc AVALAIBLE. PLEASE LOOK IN THE DIRECTORY ./Make.inc FOR
	@echo EXAMPLES OF Makefile.inc FILES, AT Make.inc/Makefile.inc.generic
	@echo IN CASE YOU NEED TO BUILD A NEW ONE AND READ THE MAIN README FILE
	@echo "##################################################################"
	@exit 1

.PHONY: default requiredobj all double simple cmplx cmplx16 libseqneeded

include ../../arch.inc
include Makefile.inc

requiredobj: Makefile.inc $(LIBSEQNEEDED) ./lib/libpord$(PLAT).a

# dummy MPI library (sequential version)

libseqneeded:
	echo r $(CC) r
	(cd libseq; make)

# Build the libpord.a library and copy it into ./lib
./lib/libpord$(PLAT).a:
	if [ "$(LPORDDIR)" != "" ] ; then \
	cd src; cd $(LPORDDIR); make CC="$(CC)" CFLAGS="$(OPTC)" AR="$(AR)" ARFUNCT= RANLIB="$(RANLIB)"; fi;
	if [ "$(LPORDDIR)" != "" ] ; then \
	cd src; cp $(LPORDDIR)/libpord.a ../lib/libpord$(PLAT).a; fi;

double: requiredobj
	(cd src ; make double)
	(cd test ; make double)
simple: requiredobj
	(cd src ; make simple)
	(cd test ; make simple)
cmplx: requiredobj
	(cd src ; make cmplx)
	(cd test ; make cmplx)
cmplx16: requiredobj
	(cd src ; make cmplx16)
	(cd test ; make cmplx16)
matlab:
	(cd MATLAB ; make all)	
clean:
	(cd src; make clean)
	(cd test; make clean)
	(cd lib; $(RM) *$(PLAT).a)
	(cd libseq; make clean)
	if [ $(LPORDDIR) != "" ] ; then \
	cd src; cd $(LPORDDIR); make realclean; fi;
	(cd MATLAB ; make clean)
