% ZUIFORMVISIBLE rend visible un formulaire
%------------------------------------------------------
% fichier zuiformvisible.m ->  zuiformvisible
%
% fonction Matlab 5 :
%	Cette fonction rend visible un formulaire. 
%
% syntaxe  :
%	zuiformvisible(handle)
%
% entrees :
%	handle       = handle du formulaire
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
function zuiformvisible(hfig)

if ishandle(hfig)
	set(hfig,'visible','on');
	figure(hfig);
end
