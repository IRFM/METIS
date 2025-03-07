function [param,data]=zinfosetanglais(param,data)


if ~strcmp(getappdata(0,'langue_cronos'),'anglais')
     return
end

%
% si deja en memoire
%
info_anglais = getappdata(0,'info_anglais');
if ~isempty(info_anglais)
     param = info_anglais.param;
     data  = info_anglais.data;
     return	  
end

%
% chargement du fichier anglais
%
root     = getappdata(0,'root');
file     = fullfile(root,'anglais','info_anglais.mat');
if exist(file,'file') == 2
   var          = load(file);
   info_anglais = var.info_anglais;
else
   return;
end



%
% boucle sur les champs
% 
info.param = param;
info.data  = data;

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
		% recupere la fin du nom
		[racine,nomc] = strtok(champc,'.');
		champus       = strcat('info_anglais',nomc);
		% test si le champ a deja ete traduit
		texte_us  = eval(champus,'[]');
      if isempty(texte_us)
				eval(strcat(champus,'='''';'));
	   end
    end
end


% memorisation
setappdata(0,'info_anglais',info_anglais);

param = info_anglais.param;
data  = info_anglais.data;
