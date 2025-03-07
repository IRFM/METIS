% ZVERBOSE affiche les messages d'etat. 
%---------------------------------------------------------------
% fichier zverbose.m ->  zverbose
%
%
% fonction Matlab 5 :
%
% Cette fonction affiche les messages d'etat. 
%  
% syntaxe  :
%  
%   zverbose(string);
%   
% remarque :
% 
%  Elle utilise la propriete 'zverbose' de root
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 03/08/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zverbose(varargin)

persistent on
persistent oni

if isempty(on)
	on = getappdata(0,'zverbose');
end

if ~isempty(on) 
	if on == 1
		sp=strrep(sprintf(varargin{:}),sprintf('\t'),'  ');
		sp=strrep(sp,sprintf('\n'),'');
		fprintf('%s\n',sp);
		if isempty(oni)
			oni = getappdata(0,'zinteraction');
		end
		if ~isempty(oni) 
			if oni == 1
				h=findobj(0,'type','uicontrol','tag','cz_info');
				if ~isempty(h) 
					set(h,'string',sp);
				end
			end
		end           
	end
end
