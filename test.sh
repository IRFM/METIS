#!/bin/tcsh -fe
make testmetis
grep CRONOSTEST log  | grep ERROR | wc -l | exit

