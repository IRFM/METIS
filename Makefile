############################################################################
#
#  Program:         cronos
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   29/08/2006
#
#  Modified:         
#
#  Send bug reports, comments or suggestions to 
#
############################################################################

include arch.inc

TARGETS=separatrice tsample op libs solver 


CLEANTARGETS=cleanseparatrice cleantsample cleanop cleanlibs cleansolver

ifeq ($(CRONOSJET),yes)
        TARGETSJET=jams jet
endif

all: $(TARGETS) $(TARGETSJET)
	@echo " ################"
	@echo " #COMPILATION OK#"
	@echo " ################"
	chmod ugo+x ./certification/test.sh
	chmod ugo+x ./metis
	chmod ugo+x ./op/ztailgrep

separatrice:
	( cd ./acces/iter; $(MAKE) all )
cleanseparatrice:
	( cd ./acces/iter; $(MAKE) clean )
tsample:
	( cd ./import/sampling; $(MAKE) all )
cleantsample:
	( cd ./import/sampling; $(MAKE) clean )
op:
	(cd ./op/xml/xmltree/@xmltree/private; $(MAKE) all )
cleanop:
	(cd ./op/xml/xmltree/@xmltree/private; $(MAKE) clean )
libs:
	(cd ./libs/metis; $(MAKE) all )
	(cd ./libs/MUMPS_4.7.3; $(MAKE) all )
cleanlibs:
	(cd ./libs/metis; $(MAKE) clean )
	(cd ./libs/MUMPS_4.7.3; $(MAKE) clean )
solver:
	(cd ./solver; $(MAKE) all )
cleansolver:
	(cd ./solver; $(MAKE) clean )

.PHONY: clean test libs solver

clean: $(CLEANTARGETS)

testmetis:
	echo "METIS TESTS";rm -f ./log; time -p certification/test.sh metis | tee -a log
testimasinmetis:
	echo "METIS TESTS (IMAS I/O)";rm -f ./log; time -p certification/test.sh imasinmetis | tee -a log
