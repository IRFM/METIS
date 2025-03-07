% METIS_TEST execute a test define in a METIS file test
%------------------------------------------------------------------------------- 
% file :  metis_test.m  ->   metis_test 
% 
% 
% function Matlab 7 : 
% 
% This function load a METIS file test created with METIS tool,
% re-run METIS and compare old and new result . 
%
% syntaxe :  
%
%  execution of the test : 
%   
%    [ref,out] = metis_test(testfilename);
%   
%    metis_test(testfilename);
%
%  list of existing test :
%
%     metis_test
%    
% input : 
% 
%   testfilename = file name of the test (all file test are in the directory "certification/metis" 
%                  of the CRONOS tree path).
%
%  
% sorties :
%	
%	ref =   structure that define the test             
%  		ref.zerod:  	      0D initial data
%             	ref.z0dinput:  	      METIS parameters
%             	ref.profil0d:         initital profiles
%     		ref info_test_metis:  environnement definition of the test (version ...)
%  
%		ref.info_test : environnement definition of the test (version ...)
%		ref.test      : data that define the test.
%
%       out = structure of results (same than ref , with new results).
%  
%  
% function wrote by J-F Artaud, tel 62-15  
% version  4.0  du  18/03/2008  
%
% CVS version  
%-------------------------------------------------------------------------------  
%  
function  [ref,out] = metis_test(varargin)

% security with disk cache
!sync


% variable global importante pour les test
% utilisation du mexfile si disponible	 
repmex=which(strcat('mexpde1dsolver.',mexext));
if ~isempty(repmex)
	mexsolver = 1;
else
	mexsolver = 0;	
end
setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);


if nargout > 0
	% donnees chargees
	ref.logfile = '';
	% donnees en sortie
	out.logfile = '';
end
if nargin == 0
 	testfilename ='';
else
	testfilename = varargin{1};
end
if nargin < 2
	plotonoff = 0;
else
	plotonoff = varargin{2};
	if isempty(plotonoff) 
		plotonoff = 0;
	end
end
% racine CRONOS
root = getappdata(0,'root');
if isempty(root)
	zineb_path;
	root = getappdata(0,'root');
end

% creation temporaire de la variable globale CRONOS_WORK_DIR dans le cas
% des tests
filename=sprintf('%s/certification',root);
setappdata(0,'CRONOS_WORK_DIR',filename);


% si pas de nom de fichier
if nargin == 0
	testfilename ='';
end
if isempty(testfilename)
	help('metis_test')
	disp(' ')
	disp('===============================================================')
	disp('list of existing METIS test :');
    dir(fullfile(root,'certification','metis','*.mat'));
	%unix(sprintf('ls %s | grep ".mat"',fullfile(root,'certification','metis')));
	return
end



% met le path cronos si besoin
if isempty(which('zinebcreetest'))
	zineb_path;
end

% ajout des path pour chargement automatique du fichier
if isdir(fullfile(root,'certification','metis'))
	addpath(fullfile(root,'certification','metis'),'-begin');
end




[void,testfilename,extvoid] = fileparts(testfilename);
file     = sprintf('%s.mat',testfilename);

if isdir(fullfile(root,'certification','output'))	
	file_log = fullfile(root,'certification','output',sprintf('%s_logfile_out@%s.txt',testfilename,datestr(now,30)));	
else
	file_log = fullfile(root,'certification',sprintf('%s_logfile_out@%s.txt',testfilename,datestr(now,30)));
end
% chargement des donnees du test
fprintf('loadind file %s\n',which(file))
try 
	ref = load(file);
        ref = ref.post;
catch
	error(sprintf('unable to load testfile %s !\n%s',testfilename,file))
end

% ouverture du journal
diary(file_log)

% separateur
disp(' ')
disp('==================================================')
fprintf('Start of %s @ %s\n',testfilename,time);
fprintf('full name of test file : \n\t%s\n\n',file);


% structure d'information 
out.info_test_metis.date        = clock;
if ispc
        out.info_test_metis.machine = computer;
        out.info_test_metis.env = getenv('OS');
else
    [s,t] = unix('uname -a');
    if s == 0
        out.info_test_metis.machine = t;
    else
        error(sprintf('error executing ''uname -a'' (%s)',t));
    end
    [s,t] = unix('env');
    if s == 0
        out.info_test_metis.env = t;
    else
        error(sprintf('error reading environnement variables  (%s)',t));
    end
end
[out.info_test_metis.version,out.info_test_metis.date_version]  = zinebversion;
out.info_test_metis.matlab_version  = version;
out.info_test_metis.toolbox_version = ver;
out.info_test_metis.path = path;
out.info_test_metis.user = getenv('USER');
out.info_test_metis.root = root;

% separateur
disp('--------------------------------------------------')

% version des machines
fprintf('test created on  :\n\t%s\n',ref.info_test_metis.machine);
fprintf('test executed on :\n\t%s\n',out.info_test_metis.machine);
disp('    ')
disp('--------------------------------------------------')

% version de cronos
if out.info_test_metis.version == ref.info_test_metis.version
	if out.info_test_metis.date_version == ref.info_test_metis.date_version
		fprintf('test created and executed with same version of CRONOS (%g from %d )\n', ...
		ref.info_test_metis.version,ref.info_test_metis.date_version);
	else
		fprintf('test created and executed with same version of CRONOS with different creation date\n');
		fprintf('test created with CRONOS version %g from %d\n', ...
		          ref.info_test_metis.version,ref.info_test_metis.date_version); 
		fprintf('test executed with CRONOS version %g from %d\n', ...
		          out.info_test_metis.version,out.info_test_metis.date_version); 
	
	end
else
	fprintf('test created and executed with different versions of CRONOS with different creation date\n');
	fprintf('test created with CRONOS version %g from %d\n', ...
		 ref.info_test_metis.version,ref.info_test_metis.date_version); 
	fprintf('test executed with CRONOS version %g from %d\n', ...
		 out.info_test_metis.version,out.info_test_metis.date_version); 
	fprintf('It''s maybe better to update the test ...\n');
end
disp('    ')
disp('--------------------------------------------------')

% version de matlab et des toolboxes
if ~strcmp(out.info_test_metis.matlab_version,ref.info_test_metis.matlab_version)
	fprintf('test created and execute with different Matlab version\n');
	fprintf('test created with Matlab version %s \n',ref.info_test_metis.matlab_version);
	fprintf('test executed with Matlab version %s \n',out.info_test_metis.matlab_version); 	
else
	fprintf('test created and execute with same Matlab version (%s)\n', ref.info_test_metis.matlab_version);
end
disp('    ')

% force le solveur matlab si necessaire
if isfield(ref.info_test_metis,'mexsolver')
	mexsolver = ref.info_test_metis.mexsolver;
end
repmex=which(strcat('mexpde1dsolver.',mexext));
if isempty(repmex) && (mexsolver == 1)
        disp('unable to run test with mexpde1dsolver : not compiled for this computer');
        mexsolver  = 0;
elseif mexsolver == 1
		disp('using fast mex PDE solver');
else
		disp('using matlab PDE solver');
end			
% securite
setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);

% appel du test
z0dinput                = ref.z0dinput;
z0dinput.option.machine = z0dinput.machine;
z0dinput
z0dinput.option
out.z0dinput = z0dinput;
if length(ref.zerod.temps) == length(ref.profil0d.temps)
	% calcul complet
	[out.zerod,void,out.profil0d] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d); 
	tol0d = 1e-3;
else
	% calcul rapide
	[out.zerod,void,out.profil0d] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
	tol0d = 1e-2;
end

% mise a disposition des donnees dans le workspace
zassignin('base','post.zerod',out.zerod);
zassignin('base','post.profil0d',out.profil0d);
z0dinput_a = z0dinput;
z0dinput_a.exp0d = ref.zerod;
z0dinput_a.mode_exp = -pi;
zassignin('base','post.z0dinput',z0dinput_a);
zassignin('base','z0dinput',z0dinput);

% la tolerance ne peut etre meilleur que la convergence lors du test
% tol0d = max(max(tol0d,max(max(ref.zerod.dini,ref.zerod.diboot),max(ref.zerod.pfus,ref.zerod.w))));

% on ne test pas les disruption
indd = find((out.zerod.disrup > 0) | (ref.zerod.disrup > 0));
nbt  = size(ref.zerod.temps,1);
out_zerod = out.zerod;
ref_zerod = ref.zerod;
if ~isempty(indd)
	noms  = fieldnames(ref.zerod);
	for k = 1:length(noms)
		var = out.zerod.(noms{k});
		if size(var,1) == nbt;
			out_zerod.(noms{k})(indd)  = [];
		end
		var = ref.zerod.(noms{k});
		if size(var,1) == nbt;
			ref_zerod.(noms{k})(indd)  = [];
		end
	end
end
out_zerod.diboot(:) = 0;
ref_zerod.diboot(:) = 0;
out_zerod.dw(:) = 0;
ref_zerod.dw(:) = 0;
out_zerod.dini(:) = 0;
ref_zerod.dini(:) = 0;
out_zerod.dpfus(:) = 0;
ref_zerod.dpfus(:) = 0;

% la source de neutre de bord est exponentiellement sensible aux variations des parametres plasma
% la source d'echange de charge est elle meme tres sensible a la source de neutres
% et souvent negligeable
out_prof = out.profil0d;
ref_prof = ref.profil0d;
out_prof.n0 = out_prof.n0 + out_prof.n0m;
out_prof.s0 = out_prof.s0 + out_prof.s0m;
ref_prof.n0 = ref_prof.n0 + ref_prof.n0m;
ref_prof.s0 = ref_prof.s0 + ref_prof.s0m;


% comparaison des structures
fprintf('difference in zerod data :\n');
zcompstruct(out_zerod,ref_zerod,tol0d);
fprintf('difference in profiles :\n');
zcompstruct(out_prof,ref_prof,tol0d);

if plotonoff  == 1
	zplotstruct(out.zerod,ref.zerod,'zerod data','METIS');
	zplotstruct(out.profil0d,ref.profil0d,'profiles','METIS');
end

% fermeture du journal
diary off

% lecture du journal 
[s,logfile] = unix(sprintf('cat %s',file_log));	
disp('    ')
disp('--------------------------------------------------')
fprintf('End of test %s\n',testfilename);
out.logfile = logfile;





