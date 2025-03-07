% traducteur automatique en anglais de la structure zinfo de cronos 
function ztranslate_info

% lecture de la structure info courante
info = zinfo('noconnexion');

%
% chargement du fichier anglais
%
root     = getappdata(0,'root');
file     = fullfile(root,'anglais','info_anglais.mat');
if exist(file,'file') == 2
   var          = load(file);
   info_anglais = var.info_anglais;
else
   info_anglais =[];
end

%
% initialisation des variables deja
%
root     = getappdata(0,'root');
filedeja = fullfile(root,'anglais','deja_anglais.mat');
if exist(filedeja,'file') == 2
     load(filedeja);
else
   deja_fr     = {};
   deja_auto   = {};
   deja_us     = {};
   deja_date   = {};
   deja_afaire = {};
end

%
% boucle sur les champs
% 

% liste des champs la structures
champ = fieldnames(info);
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('info.',champ{k});
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
	   fprintf('traitement de : %s\n',champc);

		% recupere la fin du nom
		[racine,nomc] = strtok(champc,'.');
		champus       = strcat('info_anglais',nomc);
		% texte francais
		texte_fr  = eval(champc,'[]');
      if isempty(texte_fr)
		     eval(strcat(champus,'=[];'));
		  	  save(file,'info_anglais');
      else
		    % test si le champ a deja ete traduit
		    texte_us  = eval(champus,'[]');
          if isempty(texte_us)
		         fait = 0;
		    else
		         fait = 1;
		    end
		
		    % selon le cas
		    if fait == 0
			    % 1- recherche si le texte existe deja
				 dejafait = 0;
				 inddeja  = strmatch(texte_fr,deja_fr,'exact');
				 if ~isempty(inddeja)
				     for k = 1:length(inddeja)
				         if strcmp(texte_fr,deja_fr{inddeja(k)})
							   dejafait = 1;
								eval(strcat(champus,'=deja_us{inddeja(k)};'));
								save(file,'info_anglais');
								break
							end
				     end
				 end
				 % 2 - sinon appel au traducteur
				 if dejafait == 0
				    [repbabel,crbabel] = ztraduit_babel(texte_fr);
					 if crbabel == 0
						  eval(strcat(champus,'=repbabel;'));
			           inddeja  = length(deja_fr) +1;
				        deja_fr{inddeja}   = texte_fr;
				        deja_us{inddeja}   = repbabel;
				        deja_date{inddeja} = clock;
				        deja_auto{inddeja} = 1;
						  fprintf('%s -> %s \n',texte_fr,repbabel);
						  save(file,'info_anglais','info_anglais');
                    save(filedeja,'deja_fr','deja_us','deja_date','deja_afaire','deja_auto');
					 else
					     fprintf('Probleme de traduction avec babelfish : %s \n',texte_fr)
					 end
				 end
			 end
		 end
    end
end
%
% sauve la structure deja
%
save(filedeja,'deja_fr','deja_us','deja_date','deja_afaire','deja_auto');


% 
% suite de la traduction des ressoureces restantes
%
if ~isempty(deja_afaire)
   while (~isempty(deja_afaire))
	   texte_fr  = deja_afaire{1};
		if isempty(texte_fr)
		   texte_fr ='';
	   end
		deja_afaire(1) =[];
		%length(deja_afaire)
		dejafait = 0;
		inddeja  = strmatch(texte_fr,deja_fr,'exact');
		if ~isempty(inddeja)
			for k = 1:length(inddeja)
				 if strcmp(texte_fr,deja_fr{inddeja(k)})
					dejafait = 1;
					break
				 end
			end
		end
		if dejafait == 0
			 [repbabel,crbabel] = ztraduit_babel(texte_fr);
			 inddeja  = length(deja_fr) +1;
			 if crbabel == 0
				 deja_fr{inddeja}   = texte_fr;
				 deja_us{inddeja}   = repbabel;
				 deja_date{inddeja} = clock;
				 deja_auto{inddeja} = 1;
				 fprintf('%s -> %s \n',texte_fr,repbabel);
			 else
				 fprintf('Probleme de traduction avec babelfish : %s \n',texte_fr)
				 deja_fr{inddeja}   = texte_fr;
				 deja_us{inddeja}   = '';
				 deja_date{inddeja} = [];
				 deja_auto{inddeja} = -1;
			 end
      	 save(filedeja,'deja_fr','deja_us','deja_date','deja_afaire','deja_auto');
		end
   end
end


%
% nouvel essai avec les problemes
%
indpb  = find(zcell2mat(deja_auto) == -1);
if ~isempty(indpb)
   for k = 1:length(indpb)
	      inddeja  = indpb(k);
			[repbabel,crbabel] = ztraduit_babel(deja_fr{inddeja});
			 if crbabel == 0
				 deja_us{inddeja}   = repbabel;
				 deja_date{inddeja} = clock;
				 deja_auto{inddeja} = 1;
				 fprintf('%s -> %s \n',deja_fr{inddeja},repbabel);
			 else
				 fprintf('Probleme de traduction avec babelfish : %s \n',texte_fr)
			 end
      	 save(filedeja,'deja_fr','deja_us','deja_date','deja_afaire','deja_auto');
	end
end

%
% fin de la traduction
%



