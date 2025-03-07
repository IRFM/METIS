% ZUICACHE cache rend invisible un objet
%------------------------------------------------------
% fichier zuicache.m ->  zuicache
%
%
% fonction Matlab 5 :
% 
% Cette fonction rend invisible un objet. 
%
% syntaxe  :
%  
%  zuicache(handle)
%
% entrees :
% 
%  handle       = handle de l'objet
%
% sortie : aucune
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 08/02/2001.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuicache(hui)

set(hui,'visible','off');
