#!/usr/bin/env bash
# Compile the METIS C mex files used by the solver for GNU Octave.
# Requires mkoctfile (Debian/Ubuntu package: octave-dev).
set -euo pipefail
cd "$(dirname "$0")/../import/sampling"
mkoctfile --mex -O2 tsample.c cround.c cspline.c ctable1.c
mkoctfile --mex -O2 tsplinet.c
mkoctfile --mex -O2 pchipslopest_mex.c
rm -f ./*.o
ls -1 ./*.mex
