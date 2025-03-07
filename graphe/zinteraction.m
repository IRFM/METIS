% ZINTERACTION gere la mise a jour des graphiques et interfaces
%---------------------------------------------------------------
% fichier zinteraction.m ->  zinteraction
%
%
% fonction Matlab 5 :
%
% Cette fonction gere la mise a jour des graphiques et interfaces. 
%  
% syntaxe  :
%
%   zinteraction;
%   
% remarque :
% 
%  Elle utilise la propriete 'zinteraction' de root
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 03/08/2000.
%
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zinteraction

persistent on
if isempty(on)
	if isappdata(0,'zinteraction')
		on = getappdata(0,'zinteraction');
	else
		on = 0;
	end
end

if ~isempty(on) 
	if on == 1
		drawnow
		zpause;
	end
end
