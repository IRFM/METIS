% ZUIDISABLE     rend un objet inactif 
%------------------------------------------------------
% fichier zuidisable.m  
%
% fonction Matlab 5 :
%	Cette fonction rend actif un objet. 
%
% syntaxe  : 
%	zuidisable(handle)
%
% entrees :
%	handle : handle du uicontrol (objet)
%
% sortie 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 09/02/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuidisable(h)
set(h,'enable','off');
