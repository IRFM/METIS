mex -setup C++
% put here the rigth command to setup the local compiler as written by the
% command mex -setup C++:
% mex -setup ...
local = pwd;
cd(fullfile(fileparts(which('zineb_path')),'import','sampling'));
disp('============================================')
disp('compiling pchipslopest_mex')
mex -v  pchipslopest_mex.c
disp('============================================')
disp('compiling tsplinet')
mex -v  tsplinet.c
disp('============================================')
disp('compiling tsample')
mex -v  tsample.c cround.c cspline.c ctable1.c
cd(local);

