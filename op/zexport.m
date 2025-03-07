% ZEXPORT: script d'exportation des donnees pour le debuggage 
%------------------------------------------------------------
% fichier zexport.m
%
%
% script Matlab 5 :
%
% Ce script exporte les donnees depuis une fonction vers
% l'espace de travail de base afin de travailer sur les 
% donnees (cf. zdataplot). Il exporte soit les strutures
% param et data, soit les structures datak et datkp1
%  
% syntaxe  de la commande:
%  
%     zexport;
% 
%
% script ecrite par J-F Artaud , poste 46-78
% version 1.0, du 23/08/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
% script d'exportation des donnees pour le debuggage 
evalin('base','clear data datak datakp1 param');
try
    assignin('base','param',param);
    assignin('base','data',data);
catch
    assignin('base','datak',datak);
    assignin('base','datakp1',datakp1);
end