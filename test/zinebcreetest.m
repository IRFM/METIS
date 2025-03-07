% ZINEBCREETEST  makes a test file for CRONOS and/or the modules
%------------------------------------------------------------------------------- 
% file :  zinebcreetest.m  ->   zinebcreetest 
% 
% 
% function Matlab  7: 
% 
% This fucntion makes a test file for CRONOS and/or the modules.  
% The test file can be a simple module test, a multi modules test
% or a complete CRONOS feature test( modules, init phase and solver).
% This function can be also used to perform simple module call.
% 
% syntaxe :  
%
%    [{test,info_test,logfile}] = zinebcreetest(data,param,time,testfilename,{typemodule,varargin})
%  
% input :  
%  
%	data         = structure data of CRONOS
%       param        = structure param of CRONOS
%       time         = time selected for the test (s)
%       testfilename = name of the output file used by zineb_test. If empty, just 
%                      perform the execution and not save the results.
%       typemodule   = optionnal argument used to choose what is the test
%	varargin     = complementary argument for one module test
%
%   typemodule :
%	'all'           = prepare tests for each CRONOS module and for CRONOS solver
%       'modules'       = prepare tests for each CRONOS module
%	'module_name'   = prepare test for the module 'module_name'. 
%
%        the list of modul_name is given by the command :
%                list = zineb1test
%
%   varargin : optionnal arguments for test on one module
%
%     ,'reprise',integer  = force the recovery mode of CRONOS 
%     ,'option',structure parameter of the module = set new parameters for the module
%     ,'nomf',module name = change the module connected to CRONOS
%     ,'plotonoff, 1  = make .png file of plots that compare data
%     ,'plotonoff, 0  = don't make .png file
%     ,'updateparam',0 = no automatic update of module parameters
%     ,'updateparam',1 = automatic update of module parameters
%     ,'debug',0       = no error trapping
%     ,'debug',1       = stop in module if error
%     ,'timestep',duration = optionnal time step for the onestep "module_name" (s).
%     ,'clean',0       = do not delete temporary files after test run
%     ,'clean',1       = delete temporary files after test run
%     ,'tolerance',1e-3  = relative error tolerance
%
%  
% output :  
%
%      test          = data structure of the test
%      info_test     = structure of test environnement
%      logfile       = text of execution logfile 
%
% test structure :
%    test structure have one filed per elementary test. 
%    This sub strucure contain a field data and a field param.
%            
% examples :
% 
%   load CRONOS simulation
%	>> zuiload my_cronos_simulation
%
%   create test for equi module and return the output for other use
%	>> [test,info_test,logfile] = zinebcreetest(data,param,800,'my_equi_test','equi');
%  
%   execute for one time the neoclassical module
%	>> [test,info_test,logfile] = zinebcreetest(data,param,800,'','neo');
%  
%   create test for all CRONOS modules
%	>> zinebcreetest(data,param,800,'my_modules_test','modules');
%
%   create test for CRONOS (modules and solver)
%	>> zinebcreetest(data,param,800,'my_cronos_test','all');
%
%   create a test for an unconnected module 
%	>> zinebcreetest(data,param,800,'my_test','fce','nomf','zremafile');
%
%   execute for one time the ICRH module with new parameters and make plot
%       >> my_param = param.cons.fci;
%       >> my_param.save = 'Yes';   
%	>> [test,info_test,logfile] = zinebcreetest(data,param,800,'','fci','option',my_param,'plotonoff',1);
%
% function write by J-F Artaud, tel 62-15 
% version  3.1  du  25/10/2006  
%  
% CVS version
%
%-------------------------------------------------------------------------------  
%  
function  [test,info_test,logfile] = zinebcreetest(data,param,time,testfilename,typemodule,varargin)


% initilaisation des sorties
test = [];
info_test = [];
logfile = '';

% test des entrees
if nargin < 4
	error('syntaxe : test = zinebcreetest(data,param,time,testfilename,{typemodule})');
elseif isempty(data) | isempty(param)
	error('you have to provided valid data for : data,param');
end
if isempty(time) & (length(data.gene.temps) <= 2)
	time = data.gene.temps(1);
elseif isempty(time)
	error('you have to provided valid time');
end
if nargin < 5
	typemodule ='all';
elseif isempty(typemodule)
	typemodule ='all';	
end

% test du temps
if (time > max(data.gene.temps)) |(time <= min(data.gene.temps))
	error('time is out of bounds');
end

% composition du nom du fichier 
root = getappdata(0,'root');
if isempty(root)
	zineb_path;
	root = getappdata(0,'root');
end
if isempty(testfilename)
	file ='';
elseif (strcmp(typemodule,'all') | strcmp(typemodule,'modules')) & isdir(fullfile(root,'certification','fullruns'))
	testfilename = strcat(testfilename,'_',typemodule);
	file = fullfile(root,'certification','fullruns',sprintf('Cronos_test_%s_%d',testfilename,fix(zinebversion*1000)));
elseif ~isempty(strmatch(typemodule,{'coef_all','sources','first','onestep'},'exact') & isdir(fullfile(root,'certification','cronosfunctions')))
	testfilename = strcat(testfilename,'_',typemodule);
	file = fullfile(root,'certification','cronosfunctions',sprintf('Cronos_test_%s_%d',testfilename,fix(zinebversion*1000)));
elseif isdir(fullfile(root,'certification','modules'))
	testfilename = strcat(testfilename,'_',typemodule);
	file = fullfile(root,'certification','modules',sprintf('Cronos_test_%s_%d',testfilename,fix(zinebversion*1000)));
else
	testfilename = strcat(testfilename,'_',typemodule);
	file = fullfile(root,'certification',sprintf('Cronos_test_%s_%d',testfilename,fix(zinebversion*1000)));
end
if isempty(file)
	% pas de test d'exsitence
elseif exist(strcat(file,'.mat'),'file')
	error(sprintf('the name %s have been already used, choose another name\n(%s)',testfilename,strcat(file,'.mat')));
end 

% ouverture du fichier de log
if  ~isempty(file)
	if isdir(fullfile(root,'certification','logfiles'))
		diary(fullfile(root,'certification','logfiles',sprintf('Cronos_test_%s_%d_logfile.txt',testfilename,fix(zinebversion*1000))));
	else
		diary(strcat(file,'_logfile.txt'))
	end
else
	tp = tempname;
	diary(tp);
end

% appel des modules
switch typemodule
case {'all','modules'}
	test = zineballtest(data,param,time,typemodule);
	if isempty(test)
		error('unable to create test');
	end
otherwise
	test = [];
	if nargin > 5
		[data_out,param_out,tolerance] = zineb1test(data,param,time,typemodule,varargin{:});
	else
		[data_out,param_out,tolerance] = zineb1test(data,param,time,typemodule);
	end
	if isempty(data_out)
		error('unable to create test');
	end
	test = setfield(test,typemodule,'data',data_out);
	test = setfield(test,typemodule,'param',param_out);
	test = setfield(test,typemodule,'tolerance',tolerance);
end

% fermeture du log
diary off


% structure d'information 
info_test.date        = clock;
[s,t] = unix('uname -a');
if s == 0
	info_test.machine = t;
else
	error(sprintf('error executing ''uname -a'' (%s)',t));
end
[s,t] = unix('env');
if s == 0
	info_test.env = t;
else
	error(sprintf('error reading environnement variables  (%s)',t));
end
[info_test.version,info_test.date_version]  = zinebversion;
info_test.matlab_version  = version;
info_test.toolbox_version = ver;
info_test.path = path;
info_test.user = getenv('USER');
info_test.root = root;


% lecture du logfile
if  ~isempty(file)
	[s,logfile] = unix(sprintf('cat %s_logfile.txt',file));	
else
	[s,logfile] = unix(sprintf('cat %s',tp));	
	delete(tp);
end

% sauvegarde des donnees
if  ~isempty(file)
	if  verLessThan('matlab','7.0')
				save(file,'logfile','info_test','test','-V6');
	else
				save(file,'logfile','info_test','test');
	end
	%save(file,'logfile','info_test','test','-V6');
	fprintf('Test file %s.mat created\n',file);
end

