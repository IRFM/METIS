% UPDATEUNTEST met a niveau la structure de donnee d'un test elementaire
%------------------------------------------------------------------------------- 
% fichier :  uodateuntest.m  ->   updateuntest 
% 
% 
% fonction Matlab 5 : 
% 
% Cette fonction met a niveau la structure de donnees d'un test elementaire 
%  
% syntaxe :  
%   
%        [data,param] = updateuntest(dataold,paramold)
%  
% fonction ecrite par Artaud , poste 62.15 
% version  3.1  du  20/02/2006  
%  
% liste des modifications :  version CVS
%  
%-------------------------------------------------------------------------------  
%  
function [data,param] = updateuntest(dataold,paramold)

% pas d'arret dans cette fonction
dbclear all

% creation des nouvelles structures
[cr,data,param]=zinit('',dataold.gene.temps,paramold.gene.nbrho, ...
                         paramold.gene.origine,min(dataold.gene.temps), ...
                         max(dataold.gene.temps),paramold.gene.nbfci,paramold.gene.nbfce, ...
                         paramold.gene.nbidn,paramold.gene.nbhyb,paramold.gene.nbglacon, ...
                         paramold.gene.nbg);
                         
% recopie de paramold dans param pour les champs existant
% liste des champs la structures
champ = fieldnames(paramold);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('paramold.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while (~isempty(champ))
	%premier champ de la liste
	champc=champ{1};
	champ(1)=[];
    vector_struct = 0;  % indique que le champ appartient a un vecteur de structure
    try
        eval(strcat('test=isstruct(',champc,');'));
    catch
        test = 0;
        vector_struct = 1;
    end
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
    elseif vector_struct == 1
        indpoint = max(find(champc =='.'));
        ffend    = champc(indpoint+1:end);
        ssini    = champc(1:indpoint-1);
        nbelss   = eval(sprintf('length(%s)',ssini));
        for klz  = 1:nbelss
            locval = eval(sprintf('%s(klz).%s',ssini,ffend));
            if ~isempty(nbelss)
                destination = strrep(ssini,'paramold.','param.');
                eval(sprintf('%s(klz).%s= locval;',destination,ffend));
            end
        end

	else
		% recopie
		source      = champc;
		valc        = [];
		eval(strcat('valc = ',source,';'));
		if ~isempty(valc)
		     destination = strrep(champc,'paramold.','param.');
		     eval(strcat(destination,'=',source,';'));
		end
	end
% fin de la boucle while	   
end

% recopie de dataold dans data pour les champs existant
% liste des champs la structures
champ = fieldnames(dataold);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('dataold.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while (~isempty(champ))
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
		% recopie
		source      = champc;
		valc        = eval(source,'[]');
		if ~isempty(valc)
		    destination = strrep(champc,'dataold.','data.');
                    % ajout du test des dimensions
		    valcc  = eval(destination,'[]');
                    if length(size(valcc)) ~= length(size(valc))
                        % le nombre de dimension n'est pas correct
                        fprintf('La donnees "%s" a un nombre de dimension incorrect, elle conservera sa valeur d''initilaisation\n',destination);
                    else
		        eval(strcat(destination,'=',source,';'));
                    end    
		end
	end
% fin de la boucle while	   
end
             

% securite sur tdeb et tfin            
param.gene.tdeb = min(data.gene.temps);
param.gene.kmin = 1;
param.gene.tfin = max(data.gene.temps);
param.gene.nbt  = length(data.gene.temps);
param.gene.kmax = param.gene.nbt;
param.gene.t    = param.gene.tfin;
if length(data.gene.temps) >1
	param.gene.dt  = mean(diff(data.gene.temps));
end

% creation de la liste complete des modules externes
champ = fieldnames(param.fonction);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('param.fonction.',champ{k});
end
fnoms= {};
% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while (~isempty(champ))
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
            
		fnoms{end+1}      = champc;
	end
% fin de la boucle while	   
end

% verification de la compatibilite du parametrage des modules externes
warning on
[voids_liste,substitution] =zlist_module;

for k=1:length(fnoms)
   nomfc = fnoms{k};
   ind = max(find(nomfc=='.') + 1);
   if isempty(ind)
      nomf = deblank(nomfc);
   else
      nomf = deblank(nomfc(ind:end));
   end
   module  = eval(nomfc);
   
   % application des regles de substitutions
   substi = eval(strrep(nomfc,'param.fonction.','substitution.'),'[];');
   if ~isempty(substi)
   	for lzk = 1:length(substi)
		szk = substi{lzk};
		if strmatch(module,szk{1},'exact')
			fprintf('MODULE SUBSTITUTION: replacing %s by %s\n',module,szk{2});
			module = szk{2};
			eval(sprintf('%s = ''%s'';',nomfc,szk{2}));
		end
	end
   end
  
   % 1er test = existance du module
   rep = which(module);
   if isempty(module)
      % on ne fait rien
   elseif isempty(rep)
          warning(sprintf('unknown function %s for external module %s !',module,nomf));
   else 
     % connexion du module      
     try
	 nb    = param.nombre.(nomf);
     catch
	 nb    = 1;
     end
     if (nb ==0)|~isfinite(nb)
		sortie.valeur=[];
		sortie.type=[];
		sortie.borne=[];
		sortie.defaut=[];
		sortie.info='';
		sortie.interface.jet='';
		sortie.interface.ts='';
		sortie.description = 'vide';                 % description (une ligne) de la fonction
		sortie.help = '';			     % nom du fichier d'aide s'il existe, sinon aide de la fonction
		sortie.gui  ='';			     % nom de l'interface graphique specifique si elle existe
		sortie.controle = '';			     % nom de la fonction de controle des valeurs si elle existe
		sortie.resume ='';			     % nom de la fonction qui ecrit le resume du parametrage
     else
		eval(['sortie = ',module,'(',num2str(nb),');cr =0;'],'cr =1;');
     end
     eval([strrep(nomfc,'.fonction.','.cons.') ' = sortie.valeur;']);

     % comparaison des structure                          
     newpar     = sortie.valeur;
     oldpar     = eval(strrep(strrep(nomfc,'.fonction.','.cons.'),'param.','paramold.'),'[]');
     difference = 0;
     if ~ isempty(newpar)
       liste = fieldnames(newpar);
       for l = 1:length(liste)
         nomc = liste{l};
         if isfield(oldpar,nomc)
             newf = getfield(newpar,nomc);
             oldf = getfield(oldpar,nomc);
               % bug corrige le 16/06/2004
              if all(size(newf) == size(oldf))  | ischar(newf)
                if strcmp(class(newf),class(oldf))
                   % champ compatible 
                   newf = oldf;
                   %verification des bornes
                   % recupere les infos
                   type        = getfield(sortie.type,nomc);
                   borne       = getfield(sortie.borne,nomc);
                   defaut      = getfield(sortie.defaut,nomc);

                   % cas des entiers
                   if ~isempty(findstr(lower(type),'integer'))|~isempty(findstr(lower(type),'entier'))
	                   newf =fix(newf);
                   end

                  % les bornes
		  if isempty(borne)
			% pas de test
                  elseif iscell(borne)
	             flag = zeros(size(newf,1));
	             for kl =1:length(borne)
                        if ischar(newf) | iscell(newf)
                            % bug corrige le 16/06/2004
                            ind =strmatch(borne{kl},newf,'exact');
			    flag(ind) = 1;
                         else
		              ind = find(borne{kl} == newf);
			      flag(ind) = 1;
		        end
	             end
	             if any(flag == 0)
                        % flag difference a 1
                        difference =1;
                        
                        % juste un warning ecrit pour
                        warning(sprintf('%s.%s out of range (in  %s of %s)',strrep(nomfc,'.fonction.','.cons.'), ...
                                 nomc,module,nomf));
                        % on corrige
                        if iscell(newf)
                              ind = find(flag == 0);
                              newf{ind} = deal(defaut);
                        elseif ischar(newf)
                           if size(newf,2) <= length(defaut)
                              newf(find(flag == 0),:) =' ';
                              ind = find(flag == 0)
                              for kl = 1:length(ind)
                                 newf(ind(kl),1:length(defaut)) = defaut;
                              end
                           else
                              memf = newf;
                              newf = char(' ' .* ones(size(memf,1),length(defaut)));
                              for kl =1:size(memf,1)
                                    if flag(kl) == 0
                                       newf(kl,:) = defaut;
                                    else
                                       newf(kl,1:size(memf,1)) = memf(kl,:);
                                    end 
                              end     
                           end
                        else
		           newf(find(flag == 0)) = defaut;
                        end
	             end
                  else
	             if isempty(newf)
		           newf(:) = defaut;
	             end
                     ind = find(~ isfinite(newf));
                     if ~ isempty(ind)
                        newf(ind) = defaut
                     end
                     newf = max(min(borne),min(max(borne),newf));
                  % fin du if iscell
                  end
                  newpar = setfield(newpar,nomc,newf);

                else
                  difference = 1;
                  warning(sprintf('%s.%s class mismatch (in  %s of %s)',strrep(nomfc,'.fonction.','.cons.'), ...
                                 nomc,module,nomf));
               % fin du if sur le test de la classe de la variable
                end
             else
                  difference = 1;
                  warning(sprintf('%s.%s size mismatch (in  %s of %s)',strrep(nomfc,'.fonction.','.cons.'), ...
                                 nomc,module,nomf));
             % fin du if sur la verification de la taille
             end 
         else
            difference = 1;
            warning(sprintf('%s.%s new field (in  %s of %s)',strrep(nomfc,'.fonction.','.cons.'), ...
                                 nomc,module,nomf));
         %  fin du if sur l'existance d'un champ dans la structure de consigne pour la fonction concernee   
         end
       %fin de la boucle sur les champs de param.cons.fonction_courante
       end 
         eval([strrep(nomfc,'.fonction.','.cons.') ' = newpar;']);
      % fin du test de presence de parametres isempty(liste)
      end
      % si diference warning et edition
      if difference == 1
               warning(sprintf('bad parameters for the function %s from external module %s !',module,nomf));
        % fin gestion difference  
      end   
   % fin du if sur l'existance de la fonction  
   end
% fin de la boucle sur les champ de param.fonction
end 


% cas specifique de la variable data.mode.rotc
if all(~ isfinite(data.mode.rotc))
   data.mode.rotc = zeros(size(data.mode.rotc));
   disp('Warning : wrong definition of "data.mode.rotc"');
end
% cas de la variable c2c
if all(~ isfinite(data.equi.c2c)) & any(isfinite(data.equi.jmoy))
    ve                       = ones(1,size(data.equi.jmoy,2));
    data.equi.c2c            = - (data.equi.rhomax  * ve).* cumtrapz(param.gene.x,data.equi.jmoy .* param.phys.mu0 .* data.equi.ri .* data.equi.vpr,2) ./   ...
                                (rpdederive(param.gene.x,data.equi.psi,2,2,2,1) ./ (data.equi.rhomax * ve));
    data.equi.c2c            =   zcentre(data.equi.c2c);
    disp('Warning : C2C recompute"');
end


% cas particulier des donnees de l'equilibre
loc = eval(param.fonction.equi);
if ~all(size(data.equi.R) == [param.gene.nbt,loc.nbrho,loc.nbtheta])
   fprintf('Warning ZINEB_VUP : equilibrium data size mismatch : reset R,Z,BR, BZ, BPHI, rhoRZ, pisRZ, df2RZ,dprRZ\n');
   param.gene.nbrhorz   = loc.nbrho;
   param.gene.nbthetarz = loc.nbtheta;
   data.equi.R         = single(NaN .* ones(param.gene.nbt,loc.nbrho,loc.nbtheta));       % R des surfaces magentiques (m)
   data.equi.Z         = single(NaN .* ones(param.gene.nbt,loc.nbrho,loc.nbtheta));       % Z des surfaces magentiques (m)
   data.equi.BR        = single(NaN .* ones(param.gene.nbt,loc.nbrho,loc.nbtheta));       % Br des surfaces magentiques (T)
   data.equi.BZ        = single(NaN .* ones(param.gene.nbt,loc.nbrho,loc.nbtheta));       % Bz des surfaces magentiques (T)
   data.equi.BPHI      = single(NaN .* ones(param.gene.nbt,loc.nbrho,loc.nbtheta));       % Bphi des surfaces magentiques (T)
   data.equi.rhoRZ     = single(NaN .* ones(param.gene.nbt,loc.nbrho));                      % rho correspondant aux surfaces magentiques (m)
   data.equi.psiRZ     = single(NaN .* ones(param.gene.nbt,loc.nbrho)); 
   data.equi.df2RZ     = single(NaN .* ones(param.gene.nbt,loc.nbrho)); 
   data.equi.dprRZ     = single(NaN .* ones(param.gene.nbt,loc.nbrho)); 
end
if ~all(size(data.equi.frmode) == [param.gene.nbt,loc.nbmode])
   fprintf('Warning ZINEB_VUP : equilibrium data size mismatch : reset frmode\n');
   param.gene.nbmoderz   = loc.nbmode;
   data.equi.frmode      = single(NaN .* ones(param.gene.nbt,loc.nbmode)); 
end

% cas paticulier de device
if ~isempty(param.fonction.machine)
	loc = eval(param.fonction.machine);
	if ~all(size(data.cons.asser.pfcur) == [param.gene.nbt,loc.valeur.nbpfcoils])
		data.cons.asser.pfcur = NaN .* ones(param.gene.nbt,sortie.valeur.nbpfcoils)
		param.gene.nbpfcoils  = sortie.nbpfcoils;
	end
end


% mise a jour du numero de version
param.gene.version_zineb = zinebversion;

% on conserve les memoires intactes
param.memoire = paramold.memoire;
param.asser   = paramold.asser;


