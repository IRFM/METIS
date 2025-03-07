% ZUISET   met a 1 un uicontrol radio ou check
%------------------------------------------------------
% fichier zuiset.m ->  zuiset
%
% fonction Matlab 5 :
%	Cette fonction met a 1 un uicontrol radio ou check. 
%
% syntaxe  :
%  zuiset(handle)
%
% entrees :
%	handle = handle du uicontrol
%
% sortie : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.1, du 20/05/2003.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuiset(h)

if  strcmp(st,'popupmenu')
   % zuireste ne s'applique pas au popup
   warning('ZUISET ne s''applique pas au popupmenu');
   return
end	
set(h,'value',1);
