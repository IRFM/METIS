function sim_deux  = zserialize(sim_in)

% boucle sur les champs
data = struct(sim_in);
% liste des champs la structures
champ = sort(fieldnames(data));
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
		%eval(sprintf('class(%s)',champc))
		champ_out = strrep(champc,'data.','sim_out.');
		type_out  = strrep(champc,'data.','sim_type.');
                type_loc  = eval(sprintf('class(%s)',champc));
		switch type_loc
		case {'double','java.lang.Double','java.lang.Long','java.lang.Integer'}
			eval(sprintf('%s=double(%s);',champ_out,champc));			
			eval(sprintf('%s=''%s'';',type_out,type_loc));			
		case 'java.lang.String'
			eval(sprintf('%s=char(%s);',champ_out,champc));			
			eval(sprintf('%s=''%s'';',type_out,type_loc));			
		case 'java.util.Hashtable'
		       	rep = eval(sprintf('char(%s.toString())',champc));
		       	rep = strrep(rep,'{','');
		       	rep = strrep(rep,'}','');
			while(~isempty(rep))
				[u,rep] = strtok(rep,',');
				champ_h = strtok(u,'=');
				champ_h(champ_h <= ' ') = [];
				val     = eval(sprintf('zserialize(%s.get(''%s''))',champc,champ_h));
				eval(sprintf('%s.%s=val.value;',champ_out,champ_h));							
				eval(sprintf('%s.%s=val.type;',type_out,champ_h));							
			end
		case {'SimulationFile','SimulationModule'}
			eval(sprintf('val=zserialize(%s);',champc));			
			eval(sprintf('%s=val.value;',champ_out));							
			eval(sprintf('%s=val.type;',type_out));							
		otherwise
			warning(sprintf('ZSERIALIZE : class unhandled (%s)',type_loc))
		end
    end
		
end

% regroupement
sim_deux.value = sim_out;
sim_deux.type  = sim_type;
