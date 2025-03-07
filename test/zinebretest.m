%ZINEBRETEST  appel de zineb1test a partir des donnnes d'un test standart
%------------------------------------------------------------------------------- 
% fichier :  zinebretest.m  ->   zinebretest 
% 
% 
% fonction Matlab 7 : 
% 
% Cette fonction reexecute un test de cronos en changeant les options
%  
% syntaxe :  
%
%     exection du test :
%
%   		[datak,param,tolerance] = zinebretest(filename,{options});
%
% entree :
%		filename = nom du fichier de test
% options :
%
%     ,'reprise',nombre  = force le mode reprise de CRONOS
%     ,'option',structure option du module = change le parametrage du module
%     ,'nomf',nom du module = change le module connecte a CRONOS 
%     ,'plotonoff, 1  = cree les fichiers .png des plot de comparaison des donnees
%     ,'plotonoff, 0  = sorties textes seulement
%     ,'updateparam',0 = pas de mise a jour automatique des parametres d'un module 
%     ,'updateparam',1 = pas mise a jour automatique des parametres d'un module 
%     ,'debug',0       = pas d'arret si plantage d'un module
%     ,'debug',1       = arret lors du plantage d'un module 
%     ,'timestep',delta_t = pas de temps optionnel pour l'option onestep
%     ,'clean',0       = n'efface pas les fichiers cree par les modules
%     ,'clean',1       = efface les fichiers cree par les modules
%     ,'tolerance',1e-3  = tolerance sur le resultat
%
%
% 
% sorties :
% 
%     datak    =   structure data a 1 temps contenant les donnees modifiees .
%                  
%     param    = structure param modifiee. 
%
%     tolerance = tolerance du test
% 
%  
%
% fonction ecrite par J-F Artaud , poste 62-15  
% version  3.1  du  01/02/06  
%  
%  
% liste des modifications :  version CVS
%-------------------------------------------------------------------------------  
%  
function [datak,param,tolerance] = zinebretest(testfilename,varargin)


% racine CRONOS
root = getappdata(0,'root');
if isempty(root)
	zineb_path;
	root = getappdata(0,'root');
end

[MPISTATUS,MPIRESULT]=unix('grep "\-DMPI" arch.inc | grep -v \#');
if MPISTATUS == 0
   setappdata(0,'MPIFLAG',1);
else
   setappdata(0,'MPIFLAG',0);
end

% met le path cronos si besoin
if isempty(which('zinebcreetest'))
	zineb_path;
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
if strfind(testfilename,'Cronos_test_') == 1
	file     = sprintf('%s.mat',testfilename);
else
	file     = sprintf('Cronos_test_%s.mat',testfilename);
end

% chargement des donnees du test
try 
	ref = load(file);
catch
	error(sprintf('unable to load testfile %s !\n%s',testfilename,file))
end

% boucle sur les test
test = ref.test;
noms = fieldnames(test);
for k=1:length(noms)
	vars = test.(noms{k});
	[datak,param,tolerance] = zineb1test(vars.data,vars.param,[],noms{k},varargin{:});
end




