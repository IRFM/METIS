% fantome de la fonction unix 
% elle appel deneb ou hercule sous linux
% pour acces sql
%
function [status,texte] = zunix(commande)

if strcmp(computer,'GLNX86')
	[s,pp] = unix('ping -c 1 deneb.cad.cea.frj');
	 if s == 0
	 	 [status,texte] = unix(sprintf('rsh deneb.cad.cea.fr -l cgc "%s"',commande));
	 else
	 	 [status,texte] = unix(sprintf('rsh hercule.cad.cea.fr -l cgc "%s"',commande));
	 end
else
	[status,texte] = unix(commande);
end

