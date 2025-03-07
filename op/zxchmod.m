% fonction de gestion des droits (execution)
function zxchmod(root,rfile)

% met les droits a executable
[s,t] = unix(sprintf('chmod a+x %s',fullfile(root,rfile)));

if s~=0
	% en cas d'erreur (si on n'est pas le proprio)
	% test si deja executable
	[s,t] = unix(sprintf('ls -l %s',fullfile(root,rfile)));
	if s~=0
		warning(t);
	else
		[d,t] = strtok(t);
		if sum(find(d=='x')) < 3
			warning(sprintf('unable to change the permission of %',fullfile(root,rfile)));
		end
	end
end

