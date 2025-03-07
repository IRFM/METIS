# This file is part of MUMPS VERSION 4.7.3
# This Version was built on Fri May 4 15:54:01 2007
#
#
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. If you would like to
# use other orderings, you need to obtain the corresponding package and modify
# the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis is now available as an internal ordering for MUMPS.
#

#LSCOTCHDIR = $(HOME)/JY/emilio/bin/generic
#LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps  -lfax -lorder -lscotch -lsymbol -ldof -lgraph -lcommon -lm

LPORDDIR = ../PORD/lib/
IPORD    = -I../PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

#LMETISDIR = ../../metis
#IMETIS    = # Metis doesn't need include files (Fortran interface avail.)
#LMETIS    = -L$(LMETISDIR) -lmetis

# The following variables will be used in the compilation process.
#ORDERINGSF = -Dscotch -Dmetis -Dpord
#ORDERINGSF = -Dmetis -Dpord
ORDERINGSF  = -Dpord
ORDERINGSC  = $(ORDERINGSF)
LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGS = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################
################################################################################

RM      = /bin/rm -f
#CC = 
FC = $(F90) 
FL = $(CC) $(FCFLAGS_MUMPS)
FL = $(F90) $(DEF_MUMPS)

AR = ar vr
RANLIB  = echo


INCSEQ = -I../libseq
LIBSEQ  =  -L../libseq -lmpiseq

LIBBLAS  = $(BLASLIB)
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimization options
OPTF   = $(FFLAGS_MUMPS)
OPTL   = $(FFLAGS_MUMPS)
ifdef MUMPS_NO_MAIN
	OPTC = $(FCFLAGS_MUMPS)
else
	OPTC   = $(FCFLAGS_MUMPS) -DMAIN_COMP 
endif
#End Optimization options

INC = $(INCSEQ)
LIB = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
