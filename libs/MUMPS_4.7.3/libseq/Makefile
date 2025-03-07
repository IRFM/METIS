# This file is part of MUMPS VERSION 4.7.3
# This Version was built on Fri May 4 15:54:01 2007
#
#
.PHONY: libmpiseq

all: libmpiseq

include ../../../arch.inc
include ../Makefile.inc

libmpiseq: libmpiseq$(PLAT).a

libmpiseq$(PLAT).a: mpi.o mpic.o elapse.o
	$(AR) libmpiseq$(PLAT).a mpi.o mpic.o elapse.o
.f.o:
	$(FC) $(OPTF) -c $*.f
.c.o:
	$(CC) $(OPTC) -I. -c $*.c

clean:
	$(RM) *.o *.a
