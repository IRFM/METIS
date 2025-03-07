% cahrge en memoire les ressource pour la version anglaise
function zload_anglais

global Zdictionnaire
global Zmots_reserves

% a completer
disp('loading english dialect ressources ...');
%
% chargement du dictionnaire disponible fr -> us
%
root     = getappdata(0,'root');
file     = fullfile(root,'anglais','deja_anglais.mat');
if exist(file,'file') == 2
     load(file);
else
   deja_fr = {};
   deja_auto = {};
   deja_us = {};
   deja_date   = {};
   deja_afaire = {};
end
dico.fr   = deja_fr;
dico.us   = deja_us;
%
% creation de l'index
%
dico.index  = NaN * ones(1,length(deja_fr));
for k = 1:length(deja_fr)
    dico.index(k) = sum(abs(deja_fr{k}));
end
%setappdata(0,'Zdictionnaire',dico);    
Zdictionnaire = dico;

% charge le dernier dictionnaire de mots reserver de cronos disponible
ww       = which('zinfo');
filename = fullfile(fileparts(ww),'zineb_mots_reserves');
try
	   info     = load(filename);
	   nomn     = fieldnames(info);
	   reserve  = getfield(info,nomn{1});
	   %setappdata(0,'Zmots_reserves',reserve);
	   Zmots_reserves = reserve;
catch
	   disp('pas de liste de mots reserves disponible ...')
	   reserve ={};    
end
