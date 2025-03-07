# This file is part of MUMPS VERSION 4.7.3
# This Version was built on Fri May 4 15:54:01 2007
#
#
.PHONY: double simple cmplx cmplx16

include ../../../arch.inc
include ../Makefile.inc

double:../lib/libdmumps$(PLAT).a
simple:../lib/libsmumps$(PLAT).a
cmplx16:../lib/libzmumps$(PLAT).a
cmplx:../lib/libcmumps$(PLAT).a

.SUFFIXES: .F .o
.F.o:
	$(FC) $(OPTF) $(INC) $(ORDERINGSF) -I. -I../include -c $*.F
.SUFFIXES: .c .o
.c.o:
	$(CC) $(OPTC) $(INC) $(CDEFS) $(IORDERINGS) $(ORDERINGSC) -I../include -c $*.c
clean:
	$(RM) *.o *.mod
################################################################################
#                            Double precision
################################################################################
../lib/libdmumps$(PLAT).a:dmumps_part1.o dmumps_part2.o dmumps_part3.o dmumps_part4.o dmumps_part5.o dmumps_part6.o dmumps_part7.o dmumps_part8.o dmumps_part9.o dmumps_buffer.o dmumps_bloc2.o dmumps_cv_load.o dmumps_cv.o dmumps_c.o dmumps_orderings.o dmumps_ooc.o dmumps_io.o dmumps_io_basic.o dmumps_io_thread.o dmumps_io_err.o dmumps_struc_def.o
	$(AR) $@ $?
	$(RANLIB) $@

dmumps_cv_load.o:dmumps_buffer.o dmumps_struc_def.o
dmumps_ooc.o : dmumps_struc_def.o
dmumps_buffer.o :dmumps_bloc2.o
dmumps_part1.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o
dmumps_part2.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o
dmumps_part3.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o
dmumps_part4.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o
dmumps_part5.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_cv.o dmumps_ooc.o
dmumps_part6.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o
dmumps_part7.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o
dmumps_part8.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o
dmumps_part9.o:dmumps_bloc2.o dmumps_buffer.o dmumps_cv_load.o dmumps_ooc.o


################################################################################
#                            Simple precision
################################################################################
../lib/libsmumps$(PLAT).a:smumps_part1.o smumps_part2.o smumps_part3.o smumps_part4.o smumps_part5.o smumps_part6.o smumps_part7.o smumps_part8.o smumps_part9.o smumps_buffer.o smumps_bloc2.o smumps_cv_load.o smumps_cv.o smumps_c.o smumps_orderings.o smumps_ooc.o smumps_io.o smumps_io_basic.o smumps_io_thread.o smumps_io_err.o smumps_struc_def.o
	$(AR) $@ $?
	$(RANLIB) $@

smumps_cv_load.o  : smumps_buffer.o smumps_struc_def.o
smumps_ooc.o : smumps_struc_def.o
smumps_buffer.o: smumps_bloc2.o
smumps_part1.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o
smumps_part2.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o
smumps_part3.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o
smumps_part4.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o
smumps_part5.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_cv.o smumps_ooc.o
smumps_part6.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o
smumps_part7.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o
smumps_part8.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o
smumps_part9.o:smumps_bloc2.o smumps_buffer.o smumps_cv_load.o smumps_ooc.o


################################################################################
#                            Complex precision
################################################################################
../lib/libcmumps$(PLAT).a:cmumps_part1.o cmumps_part2.o cmumps_part3.o cmumps_part4.o cmumps_part5.o cmumps_part6.o cmumps_part7.o cmumps_part8.o cmumps_part9.o cmumps_buffer.o cmumps_bloc2.o cmumps_cv_load.o cmumps_cv.o cmumps_c.o cmumps_orderings.o cmumps_ooc.o cmumps_io.o cmumps_io_basic.o cmumps_io_thread.o cmumps_io_err.o cmumps_struc_def.o
	$(AR) $@ $?
	$(RANLIB) $@

cmumps_cv_load.o  : cmumps_buffer.o cmumps_struc_def.o
cmumps_ooc.o : cmumps_struc_def.o
cmumps_buffer.o: cmumps_bloc2.o
cmumps_part1.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o
cmumps_part2.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o
cmumps_part3.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o
cmumps_part4.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o
cmumps_part5.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_cv.o cmumps_ooc.o
cmumps_part6.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o
cmumps_part7.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o
cmumps_part8.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o
cmumps_part9.o:cmumps_bloc2.o cmumps_buffer.o cmumps_cv_load.o cmumps_ooc.o


################################################################################
#                            Double complex precision
################################################################################
../lib/libzmumps$(PLAT).a:zmumps_part1.o zmumps_part2.o zmumps_part3.o zmumps_part4.o zmumps_part5.o zmumps_part6.o zmumps_part7.o zmumps_part8.o zmumps_part9.o zmumps_buffer.o zmumps_bloc2.o zmumps_cv_load.o zmumps_cv.o zmumps_c.o zmumps_orderings.o zmumps_ooc.o zmumps_io.o zmumps_io_basic.o zmumps_io_thread.o zmumps_io_err.o zmumps_struc_def.o
	$(AR) $@ $?
	$(RANLIB) $@

zmumps_cv_load.o  : zmumps_buffer.o zmumps_struc_def.o
zmumps_ooc.o : zmumps_struc_def.o
zmumps_buffer.o: zmumps_bloc2.o
zmumps_part1.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o
zmumps_part2.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o
zmumps_part3.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o
zmumps_part4.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o
zmumps_part5.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_cv.o zmumps_ooc.o
zmumps_part6.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o
zmumps_part7.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o
zmumps_part8.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o
zmumps_part9.o:zmumps_bloc2.o zmumps_buffer.o zmumps_cv_load.o zmumps_ooc.o


