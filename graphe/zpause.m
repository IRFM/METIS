% ZPAUSE gere le bouton pause. 
%---------------------------------------------------------------
% fichier zpause.m ->  zpause
%
%
% fonction Matlab 5 :
%
% Cette fonction gere le bouton pause. 
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
function zpause

info.h.run       = findobj(0,'type','uicontrol','tag','cz_run');
info.h.pause     = findobj(0,'type','uicontrol','tag','cz_pause');
info.h.keyboard  = findobj(0,'type','uicontrol','tag','cz_keyboard');
info.h.fin       = findobj(0,'type','uicontrol','tag','cz_fin');
info.h.step      = findobj(0,'type','uicontrol','tag','cz_step');

if isempty(info.h.run) | isempty(info.h.pause) | isempty(info.h.keyboard) | isempty(info.h.fin) | isempty(info.h.step)
	return
end

if get(info.h.pause,'value') == 1
	set(info.h.run,'value',0);
else
		return		
end	
set(info.h.step,'value',0);

info=zbouton(info);

while (info.run == 0)
	pause(1);
	info=zbouton(info);
	if get(info.h.step,'value') == 1
		break;
	end
	
end

