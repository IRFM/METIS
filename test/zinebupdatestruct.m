% ZINEBUPDATESTRUCT  mise a niveau les structures de donneees d'un test de CRONOS  
%------------------------------------------------------------------------------- 
% fichier :  zinebupdatestuct.m  ->   zinebupdatestruct 
% 
% 
% fonction Matlab 5 : 
% 
% Cette fonction met a niveau les structures de donneees d'un test de CRONOS
%  
% syntaxe :  
%   
%        zinebupdatestruct testfilename
%        [ref,out] = zinebupdatestruct('testfilename');
% 
% entrees :  
%
%          testfilename = nom du fichier de test
%  
% sorties :  
%  
%          ref   = jeu de donnees initial
%          out   = nouveau jeu de donnees mis a niveau 
%  
% fonction ecrite par Artaud , poste 62.15 
% version  4.1  du  07/08/2008  
%  
% liste des modifications :  version CVS
%  
%-------------------------------------------------------------------------------  
%  
function  [ref,out] = zinebupdatestruct(testfilename)

% donnees chargees
ref.test = [];
ref.info_test = [];
ref.logfile = '';

% donnees en sortie
out.test = [];
out.info_test = [];
out.logfile = '';


% met le path cronos si besoin
if isempty(which('zinebcreetest'))
	zineb_path;
end

% nom complet du fichier de test
root = getappdata(0,'root');
if isempty(root)
	zineb_path;
	root = getappdata(0,'root');
end
% ajout des path pour chargement automatique du fichier
if isdir(fullfile(root,'certification','fullruns'))
	addpath(fullfile(root,'certification','fullruns'),'-begin');
end
if isdir(fullfile(root,'certification','cronosfunctions'))
	addpath(fullfile(root,'certification','cronosfunctions'),'-begin');
end
if isdir(fullfile(root,'certification','modules'))
	addpath(fullfile(root,'certification','modules'),'-begin');
end

[void,testfilename,extvoid] = fileparts(testfilename);
[vertest,filnamenoversion]  = strtok(testfilename(end:-1:1),'_');
filnamenoversion            = filnamenoversion(end:-1:1);
vertest                     = str2num(vertest(end:-1:1))./1000;
ajout = '';

if strfind(testfilename,'Cronos_test_') == 1
	file     = which(sprintf('%s.mat',testfilename));
	file_log = fullfile(root,'certification','logfiles',sprintf('%s_logfile_update.txt',testfilename));
	file_update = file;
else
	file     = which(sprintf('Cronos_test_%s.mat',testfilename));
	file_log = fullfile(root,'certification','logfiles',sprintf('Cronos_test_%s_logfile_update.txt',testfilename));
	file_update = file;
end

% chargement des donnees du test
try 
	ref = load(file);
catch
	error(sprintf('unable to load testfile %s !\n%s',testfilename,file))
end

% ouverture du journal
diary(file_log)
% separateur
disp(' ')
disp('==================================================')
fprintf('Start update structures of %s @ %s\n',testfilename,time);
fprintf('full name of file to be updated is : \n\t%s\n\n',file);
fprintf('full name of file updated will be : \n\t%s.mat\n\n',file_update);

% structure d'information 
info_test = ref.info_test;
[info_test.version,info_test.date_version]  = zinebversion;

% mise a jour des sous structure test
test = ref.test;

% boucle sur les test
namelist = fieldnames(test);
for k=1:length(namelist)
	module = namelist{k};
	fprintf('updating test stuctures for %s\n',module);
	dataold   = getfield(test,module,'data');
	paramold  = getfield(test,module,'param');
	paramold.gene.nbt = length(dataold.gene.temps);
	[data,param] = updateuntest(dataold,paramold);
	if isfield(getfield(test,module),'tolerance')
		tolerance = getfield(test,module,'tolerance');
	else
		tolerance = 1e-3;
	end
	test = setfield(test,module,'data',data);
	test = setfield(test,module,'param',param);
	test = setfield(test,module,'tolerance',tolerance);
end
out.test = test;

% fermeture du journal
diary off
% lecture du journal 
[s,logfile] = unix(sprintf('cat %s',file_log));	

% sauvegarde du nouveau test
if  verLessThan('matlab','7.0')
			save(file_update,'logfile','info_test','test','-V6');
else
			save(file_update,'logfile','info_test','test');
end
%save(file_update,'logfile','info_test','test','-V6');
fprintf('Test file %s.mat updated\n',file_update);


% fin
disp('    ')
disp('--------------------------------------------------')
fprintf('End of update of %s\n',testfilename);
out.logfile = logfile;



