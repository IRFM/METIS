% ZUIENABLE		rend actif un objet
%------------------------------------------------------
% fichier zuienable.m ->  zuienable
%
% fonction Matlab 5 :
%	Cette fonction rend actif un objet. 
%
% syntaxe  :
%	zuienable(handle)
%
% entrees :
%	handle = handle du uicontrol
%
% sortie : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 09/02/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuienable(h)
set(h,'enable','on');
