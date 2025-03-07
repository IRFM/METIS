#!/bin/tcsh -fe
cp -fp ./architecture/arch.inc.gateway.gfortran.matlab2013b ./arch.inc
make clean
make
exit
