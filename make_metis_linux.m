local = pwd;
cd(fullfile(local,'import','sampling'));
disp('============================================')
disp('compiling pchipslopest_mex')
mex -v -cxx pchipslopest_mex.c
disp('============================================')
disp('compiling tsplinet')
mex -v -cxx tsplinet.c
disp('============================================')
disp('compiling tsample')
mex -v -cxx tsample.c cround.c cspline.c ctable1.c
cd(local);

% try to install cocos tools
if isempty(which('init_cocos_test_in_metis'))
    zineb_path;
end
try 
    init_cocos_test_in_metis
catch
    disp('Error during installation of COCOS tools');
    lasterr
end
