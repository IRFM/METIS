# This file is part of MUMPS VERSION 4.7.3
# This Version was built on Fri May 4 15:54:01 2007
#
#


default:double

all: simple double cmplx cmplx16

simple : ssimpletest
double : dsimpletest c_example
cmplx  : csimpletest
cmplx16: zsimpletest

include ../../../arch.inc
include ../Makefile.inc

.SUFFIXES: .F .o
.F.o:
	$(FC) $(OPTF) $(INC) -I. -I../include -c $*.F

.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) -I. -I../include -c $*.c

clean:
	$(RM) *.o ssimpletest dsimpletest csimpletest zsimpletest c_example

dsimpletest: ../lib/libdmumps$(PLAT).a dsimpletest.o
	$(FL) -o dsimpletest $(OPTL) dsimpletest.o ../lib/libdmumps$(PLAT).a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)
c_example: ../lib/libdmumps$(PLAT).a c_example.o
	$(FL) -o c_example $(OPTL) c_example.o ../lib/libdmumps$(PLAT).a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)

ssimpletest: ../lib/libsmumps$(PLAT).a ssimpletest.o
	$(FL) -o ssimpletest $(OPTL) ssimpletest.o ../lib/libsmumps$(PLAT).a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)

zsimpletest: ../lib/libzmumps$(PLAT).a zsimpletest.o
	$(FL) -o zsimpletest $(OPTL) zsimpletest.o ../lib/libzmumps$(PLAT).a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)

csimpletest: ../lib/libcmumps$(PLAT).a csimpletest.o
	$(FL) -o csimpletest $(OPTL) csimpletest.o ../lib/libcmumps$(PLAT).a $(LIB) $(LIBBLAS) $(LIBOTHERS) $(LORDERINGS)


../lib/libsmumps$(PLAT).a:
	@ echo 'Error: you should build the library ../lib/libsmumps$(PLAT).a first'
	exit -1
../lib/libdmumps$(PLAT).a:
	@ echo 'Error: you should build the library ../lib/libdmumps$(PLAT).a first'
	exit -1
../lib/libcmumps$(PLAT).a:
	@ echo 'Error: you should build the library ../lib/libcmumps$(PLAT).a first'
	exit -1
../lib/libzmumps$(PLAT).a:
	@ echo 'Error: you should build the library ../lib/libzmumps$(PLAT).a first'
	exit -1


