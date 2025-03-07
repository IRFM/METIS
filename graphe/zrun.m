% ZRUN gere le bouton run. 
%---------------------------------------------------------------
% fichier zrun.m ->  zrun
%
%
% fonction Matlab 5 :
%
% Cette fonction gere le bouton run. 
%  
% syntaxe  :
%  
%   zpause;
%   
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 03/08/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zrun

info.h.run       = findobj(0,'type','uicontrol','tag','cz_run');
info.h.pause     = findobj(0,'type','uicontrol','tag','cz_pause');
info.h.keyboard  = findobj(0,'type','uicontrol','tag','cz_keyboard');
info.h.fin       = findobj(0,'type','uicontrol','tag','cz_fin');
info.h.step      = findobj(0,'type','uicontrol','tag','cz_step');

if isempty(info.h.run) | isempty(info.h.pause) | isempty(info.h.keyboard) | isempty(info.h.fin) | isempty(info.h.step)
	return
end

set(info.h.pause,'value',0);
set(info.h.step,'value',0);

info=zbouton(info);

