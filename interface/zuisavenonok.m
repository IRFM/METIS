% ZUISAVENONOK positionne le flag indiquant des donnees a sauver
%-----------------------------------------------------------------
% fichier zuisavenonok.m ->  zuisavenonok
%
% fonction Matlab 5 :
%	Cette fonction positionne le flag indiquant des donnees a sauver. 
%
% syntaxe  :
%  zuisavenonok
%
% entrees : aucune
%
% sortie  : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.6, du 29/06/2001.
% 
% liste des modifications : 
%	* 29/08/2001 -> ajout selected(on) sur le botton sauve
%
%--------------------------------------------------------------
%
function zuisavenonok

zassignin('base','param.edit.saveok',0);
[hfig,h] = zuiformhandle('direct');
if ishandle(hfig)
	set(h.radio_savefile,'foregroundcolor',[0.5 0 0]);
end
