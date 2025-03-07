mex -setup:'C:\Program Files\MATLAB\R2022b\bin\win64\mexopts\msvc2022.xml' C
cd import\sampling
mex pchipslopest_mex.c
mex tsplinet.c 
mex tsample.c cround.c cspline.c ctable1.c

