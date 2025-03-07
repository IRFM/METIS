% ZASSIGIN comme assignin, mais marche avec les structure
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 23/03/2005.
% 
%--------------------------------------------------------------
%
function zassignin(ws,name,v)


% export la donnees
if isempty(v)
	if any(name == '(')
		if ischar(v)
			evalin(ws,strcat(name,'=''?'';'));
		else
			try
				evalin(ws,sprintf('%s = NaN .* ones(size(%s));',name,name));
			catch
				evalin(ws,sprintf('%s = NaN;',name));
			end
		end
	else
		assignin(ws,'variable_intermediaire_zineb',v);
		evalin(ws,strcat(name,'=variable_intermediaire_zineb;'));
	end
else
	assignin(ws,'variable_intermediaire_zineb',v);
	evalin(ws,strcat(name,'=variable_intermediaire_zineb;'));
end
