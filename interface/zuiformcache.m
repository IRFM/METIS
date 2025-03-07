% ZUIFORMCACHE rend invisible un formulaire
%------------------------------------------------------
% fichier zuiformcache.m ->  zuiformcache
%
%
% fonction Matlab 5 :
% 
% Cette fonction rend invisible un formulaire. 
%
% syntaxe  :
%  
%  zuiformcache(handle)
%
% entrees :
% 
%  handle       = handle du formulaire
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
function zuiformcache(hfig)

if ishandle(hfig)
	set(hfig,'visible','off');
end
