% ZUIRESET    met a 0 un uicontrol radio ou check
%-------------------------------------------------
% fichier zuireset.m ->  zuireset
%
% fonction Matlab 5 :
%	Cette fonction met a 0 un uicontrol radio ou check. 
%
% syntaxe  :
%  zuireset(handle)
%
% entrees :
%	handle = handle du uicontrol
%
% sortie : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.1, du 20/05/2002.
% 
% liste des modifications : 
%
%  * 20/05/2003  -> protection des popup
%
%--------------------------------------------------------------
%
function zuireset(h)

st = get(h,'style');
if  strcmp(st,'popupmenu')
   % zuireste ne s'applique pas au popup
   warning('ZUIRESET ne s''applique pas au popupmenu');
   return
end	
set(h,'value',0);

