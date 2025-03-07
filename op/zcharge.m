% attend que la charge soit < x pour lancer l'execution
% arrete le process a 7h s'il est encore en attente
function zcharge(seuil)

if nargin == 0
	seuil = 1;
elseif isempty(seuil)
	seuil = 1;
end

% pour le mail
quand = datestr(now);

charge = inf;
while (charge > seuil) 
	[s,t]=unix('uptime');
	if ((~isempty(t))&(s==0))
		indice=max(find(t==':'));
	else
		indice=[];
	end
	if ~isempty(indice)
		charge=str2num(['[',t((indice+1):(length(t)-1)),']']);
		charge =charge(3);
	else
		charge=inf;
	end
	
	% time-out
	clk = clock;
	if clk(4)  == 7
		disp('c''est trop tard')
		message = sprintf('le job du %s n''a pas pu etre soumis',quand);
		cr=zmail(getenv('USER'),'Job non soumis',message);
		
		exit
	end
	
	fprintf('zcharge: %s  => %g / %d (%g)\n',datestr(now),charge,seuil,cputime);

	% attente avant mesure
	if charge > seuil
		pause(60);
	else
		break;
	end
	
end
