% ZUIVISIBLE rend visible un objet
%------------------------------------------------------
% fichier zuivisible.m ->  zuivisible
%
% fonction Matlab 5 :
%	Cette fonction rend visible un objet. 
%
% syntaxe  :
%	zuivisible(handle)
%
% entrees :
%  handle = handle de l'objet
%
% sortie : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 08/02/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuivisible(hui)

set(hui,'visible','on');
