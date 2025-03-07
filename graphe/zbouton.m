% ZBOUTON gere la coherence entre les boutons 
%---------------------------------------------------------------
% fichier zbouton.m ->  zbouton
%
%
% fonction Matlab 5 :
%
% Cette fonction gere la coherence entre les boutons . 
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
function info=zbouton(info)

% test des handles
if ishandle(info.h.run) & ishandle(info.h.pause) & ishandle(info.h.fin)
	info.run = get(info.h.run,'value');
	if info.run == 0 
		set(info.h.pause,'value',1);
	end
	info.fin = get(info.h.fin,'value');
else
	info.run = 1;
	info.fin = 0;
end	

% etats des boutons
set(info.h.run,'value', info.fin | info.run);
set(info.h.pause,'value', ~info.fin & ~info.run);
