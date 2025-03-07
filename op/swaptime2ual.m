% cette fonction deplace le temps a la racine de la structure matlab
% s_in est une structure de type cronos (s_in.data(temps,espace)) 
% s_out est une structure compatible avec l'UAL (s_out(temps).data(espace)) pointe vers un cpo
% cpo_name est le nom du cpo
function s_out = swaptime2ual(s_out,s_in,cpo_name)

% recherche du nombre de temps
if isfield(s_in,'time')
	nbt = length(s_in.time);
	time_cpo =1;
else
	nbt = 1; 
	time_cpo =0;
end

% liste des champs la structures
champ = fieldnames(s_in);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('s_in.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while(~isempty(champ))
    %premier champ de la liste
    champc=champ{1};
    champ(1)=[];
    locvar = eval(champc);
    %eval(strcat('test=isstruct(',champc,');'));
    test= isstruct(locvar);
    len = length(locvar);
    if test & (len ==1) 
    	% cas d'une sous structure -> ajout de champs
    	%eval(strcat('champnew=fieldnames(',champc,');'));   	
    	champnew=fieldnames(locvar);	
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
    elseif test 
%	try
%  		% cas des structure interne a plusieur elements
%  		% pour le moment seul la structure position(k).R &  position(k).Z est traite
%  		infield = fieldnames(locvar);
%  		% chaines utiles
%  		sn = strrep(champc,'s_in.','');
%  		sk = strrep(champc,'s_in.','s_out(kz).');
%  		snk = strrep(champc,'s_in.','s_out.');
%  		ss_in =size(locvar);
%  		vide = isempty(locvar);
%  		isanum = isnumeric(locvar);
%  		isacell = iscell(locvar);
%  		for lz = 1:length(infield)
%  			% boucle sur les noms de champs
%  			if (ss_in(1) == nbt) & (nbt > 1)
%  				for kz = 1:nbt
%  					loc = squeeze(locvar(kz,:,:));
%  					%s_out(kz) = eval(sprintf('allocate_%s(s_out(kz),''%s.%s'',loc.(%s))',cpo_name,sn,infield{lz},infield{lz}));
%  					if length(loc) == 1
%  						eval(sprintf('s_out(kz).%s.%s=loc.%s;',sn,infield{lz},infield{lz}));
%  					else
%  						for lmz = 1:length(loc)
%  							eval(sprintf('s_out(kz).%s(lmz).%s=loc(lmz).%s;',sn,infield{lz},infield{lz}));                    
%  						end
%  					end
%  				end
%  			elseif  (~vide) & (time_cpo == 0) 
%  				loc = squeeze(locvar(1,:,:));
%  				%s_out = eval(sprintf('allocate_%s(s_out,''%s.%s'',loc.(%s))',cpo_name,sn,infield{lz},infield{lz}));
%  				if length(loc) == 1
%  					eval(sprintf('s_out.%s.%s=loc.%s;',sn,infield{lz},infield{lz}));
%  				else
%  					for lmz = 1:length(loc)
%  						eval(sprintf('s_out.%s(lmz).%s=loc(lmz).%s;',sn,infield{lz},infield{lz}));                   
%  					end              
%  				end
%  			else
%  				kz =1; 
%  				loc = squeeze(locvar(kz,:,:));
%  				%s_out(kz) = eval(sprintf('allocate_%s(s_out(kz),''%s.%s'',loc.(%s))',cpo_name,sn,infield{lz},infield{lz}));
%  				if length(loc) == 1
%  					eval(sprintf('s_out(kz).%s.%s=loc.%s;',sn,infield{lz},infield{lz}));
%  				else
%  					for lmz = 1:length(loc)
%  						eval(sprintf('s_out(kz).%s(lmz).%s=loc(lmz).%s;',sn,infield{lz},infield{lz}));
%  			
%  					end
%  				end
%  			end
%  		end
%          catch
		% cas autre de structure indexee
		champloc = fieldnames(locvar);
                champnew = {};	
		for k=1:length(champloc)
			for kzz  = 1:len
				% nom complet pour acceder a la variable
				champnew{end+1}=sprintf('%s(%d).%s',champc,kzz,champloc{k});
			end
		end
		% ajout a la liste des champs
		if isempty(champ)
			champ =champnew';
		else
			champ=cat(1,champ,champnew');
		end
%        end
    else
        % juste pour les tests
        %disp(champc);
		
        % chaines utiles
        sn = strrep(champc,'s_in.','');
        sk = strrep(champc,'s_in.','s_out(kz).');
        snk = strrep(champc,'s_in.','s_out.');
 		
        %eval(strcat('ss_in =size(',champc,');'));
        %eval(strcat('vide = isempty(',champc,');'));
	%eval(sprintf('isanum = isnumeric(%s);',champc));
	%eval(sprintf('isacell = iscell(%s);',champc));
        ss_in =size(locvar);
        vide = isempty(locvar);
	isanum = isnumeric(locvar);
	isacell = iscell(locvar);
        if (ss_in(1) == nbt) & (nbt > 1)
            nbsize = length(ss_in);
            for kz =1:nbt
	    	loc = [];
		if isacell
	    		loc = squeeze(locvar{kz});
	    		%loc = squeeze(eval(sprintf('%s{kz}',champc)));
			%s_out(kz) = eval(sprintf('allocate_%s(s_out(kz),''%s'',loc)',cpo_name,sn));			
			eval(sprintf('s_out(kz).%s = loc',sn));
		else
			switch nbsize
			case {1,2}
				loc = locvar(kz,:);
			case 3
				loc = locvar(kz,:,:);
			case 4
			 	loc = locvar(kz,:,:,:);
			case 5
				loc = locvat(kz,:,:,:,:);
			case 6
				loc = locvar(kz,:,:,:,:,:);
			case 7
				loc = locvar(kz,:,:,:,:,:,:);
			otherwise
				error('this matrix dimension is not yet implemented');
			end
			loc = squeeze(loc);
			if ~isanum
				eval(sprintf('%s = loc;',sk));
			elseif all(size(loc) == 1)
				try
                 			%eval(sprintf('%s = squeeze(%s);',sk,champc));
                 			eval(sprintf('%s = loc;',sk));
				catch
					%s_out(kz) = eval(sprintf('allocate_%s(s_out(kz),''%s'',loc)',cpo_name,sn));
					eval(sprintf('s_out(kz).%s = loc;',sn));
				end
           	        else
                		loc = freenan(loc); % type veteur ayant qu'un seul element
				%s_out(kz) = eval(sprintf('allocate_%s(s_out(kz),''%s'',loc)',cpo_name,sn));
                try
				      eval(sprintf('s_out(kz).%s = loc;',sn));
                catch
                    try
                      indpoint = find(sn =='.',1,'last');
                      snmem    = sn;
                      sn       = strcat(sn(1:indpoint-1),'(1)',sn(indpoint:end));
                      eval(sprintf('s_out(kz).%s = loc;',sn));
                    catch
                       indpointm1 =  find(sn =='.');
                       indpointm1 = max(indpointm1(indpointm1<indpoint));
                       if ~isempty(indpointm1)
                            sn       = strcat(snmem(1:indpointm1-1),'(1)',snmem(indpointm1:end));
                            eval(sprintf('s_out(kz).%s = loc;',sn));                          
                       end
                    end
                end
			end
		end
            end
        elseif (~vide) & (time_cpo == 0) 
	    if isacell
	    	loc = squeeze(locvar{1});
		%s_out = eval(sprintf('allocate_%s(s_out,''%s'',loc)',cpo_name,sn));			
		eval(sprintf('s_out.%s = loc;',sn));
            elseif isanum
	    	loc = squeeze(locvar);
                if all(size(loc) == 1)
			try
                 		eval(sprintf('%s = loc;',snk));
			catch
				%s_out = eval(sprintf('allocate_%s(s_out,''%s'',loc)',cpo_name,sn));
				eval(sprintf('s_out.%s = loc;',sn));
			end
                else
			%s_out = eval(sprintf('allocate_%s(s_out,''%s'',loc)',cpo_name,sn));
			eval(sprintf('s_out.%s = loc;',sn));
		end
               
            else    
 		try
                 	eval(sprintf('%s = locvar;',snk));
		catch
			%s_out = eval(sprintf('allocate_%s(s_out,''%s'',locvar)',cpo_name,sn));
			eval(sprintf('s_out.%s = loc;',sn));
		end
                %eval(sprintf('%s = locvar;',snk));
            end
	elseif (~vide) 
	    kz = 1; 
 	    if isacell
		for kz = 1:length(locvar);
	    		loc = squeeze(locvar{kz});
			%s_out(kz) = eval(sprintf('allocate_%s(s_out(kz),''%s'',loc)',cpo_name,sn));
			eval(sprintf('s_out(kz).%s = loc;',sn));
		end			
            elseif isanum
	    	loc = squeeze(locvar);
                if all(size(loc) == 1)
			try
                 		eval(sprintf('%s = loc;',sk));
			catch
				%s_out(1) = eval(sprintf('allocate_%s(s_out(1),''%s'',loc)',cpo_name,sn));
				eval(sprintf('s_out(1).%s = loc;',sn));
			end
                else
			%s_out(1) = eval(sprintf('allocate_%s(s_out(1),''%s'',loc)',cpo_name,sn));
            try
                eval(sprintf('s_out(1).%s = loc;',sn));
            catch
                snz = strrep(sn,'.','(1).');
                eval(sprintf('s_out(1).%s = loc;',snz));
            end
		end
               
            else    
		try
			try
                 		eval(sprintf('%s = locvar;',sk));
			catch
                      		indpoint = find(sk =='.',1,'last');
                      		sk       = strcat(sk(1:indpoint-1),'(1)',sk(indpoint:end));
                      		eval(sprintf('%s = locvar;',sk));
			end
		catch
			try
				%s_out(1) = eval(sprintf('allocate_%s(s_out(1),''%s'',locvar)',cpo_name,sn));
				eval(sprintf('s_out(1).%s = locvar;',sn));
			catch
                      		indpoint = find(sn =='.',1,'last');
                      		sn       = strcat(sn(1:indpoint-1),'(1)',sn(indpoint:end));
                      		eval(sprintf('s_out(1).%s = locvar;',sn));
			end
		end
                %eval(sprintf('%s = locvar;',sk));
            end
   
        end
    end
end

function loc = freenan(loc)

si =size(loc);
% pour le moment uniquement pour les vecteurs
if (length(si) == 2) && (si(1) == 1)
    loc(isnan(loc)) = [];
end





