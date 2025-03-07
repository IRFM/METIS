% fonction d'auto test de metis4imas
function metis4imasautotest(shot,run)

%[error_flag,xsd] = metis4imas(shot,run,occurrence,methode,time,codeparam_filename,interpolation_methode)
% simple test
disp('simple test, creation of xml and xsd')
metis4imas;
disp('simple test, run of test case')
[error_flag,output] = metis4imas(shot,run,'','test');
% re lecture
[error_flag,output] = metis4imas(shot,run,'','read');
% test calcul complet avec entree
disp('run in fast mode')
metis4imas(shot,run,'','fast')
disp('run in full mode')
metis4imas(shot,run,'','full')
disp('initialisation evolution mode')
metis4imas(shot,run,'','init',1);
for k=2:10
	disp('one time evolution mode')
	k
	% k is time
	metis4imas(shot,run,'','one_time',k);
end
% save with restart
disp('test restart');
option.restart='test_restart';
metis4imas(shot,run,'','one_time',k+1,option);
metis4imas(shot,run,'','test_restart',k-2);
metis4imas(shot,run,'','one_time',k-1);
metis4imas(shot,run,'','one_time',k);
% same with reset cpo
disp('test init output IDSs');
option.init_output_imas =1;
metis4imas(shot,run,'','init',1,option);
for k=2:10
	disp('one time evolution mode')
	k
	% k is time
	metis4imas(shot,run,'','one_time',k,option);
end
