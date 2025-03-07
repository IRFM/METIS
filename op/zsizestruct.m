function  zsizestruct(s1)
% -------------------------------------------------------------------
% List all variables included into structures and substructures of
% the general workspace and display their size
% -------------------------------------------------------------------


% preparation
data.s1jshjehcq = s1;
clef1 = 's1jshjehcq';


% boucle sur les champs
% liste des champs la structures
champ = fieldnames(data);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('data.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    eval(strcat('test=isstruct(',champc,');'));
    if test
    	% cas d'une sous structure -> ajout de champs
    	eval(strcat('champnew=fieldnames(',champc,');'));   	
    	for k=1:length(champnew)
    			% nom complet pour acceder a la variable
	        	champnew{k}=strcat(champc,'.',champnew{k});
		end
		% ajout a la liste des champs
		if isempty(champ)
	    	champ =champnew;
		else
	    	champ=cat(1,champ,champnew);
		end
    else
	 	% juste pour les tests
     		%disp(champc);
		
		if strmatch(strcat('data.',clef1),champc)
			% test complet
			champs1 = champc;
			
			eval(sprintf('var1 = %s ;',champs1));
			fprintf('size(%s) = ',champs1);
			fprintf('%d ',size(var1));
			fprintf('\n');
		end
		
		%fprintf('%s -> %g\n',champc,real(somme));
    end
end

