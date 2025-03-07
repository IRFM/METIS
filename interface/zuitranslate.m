function sortie = zuitranslate(entree)

global Zdictionnaire
global Zmots_reserves
global deja_anglais

% test de l'entree
if nargin == 0
	sortie ='';
	return
elseif isempty(entree)
	sortie ='';
	return
end

% mot reserve de zineb -> pas de traduction
reserve = Zmots_reserves;
if ~isempty(reserve)
	if strmatch(lower(entree),reserve,'exact')
		sortie = entree;
		return
	end
end		

    	
% lit le dictionnaire
dico = Zdictionnaire;    
if isempty(dico)
	sortie = entree;
	return
end
	
% recherche du texte a traduire
indexsum = find(dico.index == sum(abs(entree)));
if isempty(indexsum)
	sortie = entree;
	%
	% ajout dans la liste des ressources a traduire
	%
	if isempty(deja_anglais)
   	root     = getappdata(0,'root');
   	file     = fullfile(root,'anglais','deja_anglais.mat');
   	if exist(file,'file') == 2
      	 deja_anglais = load(file);
   	else
     	 	deja_anglais.deja_fr     = {};
      	deja_anglais.deja_us     = {};
      	deja_anglais.deja_date   = {};
      	deja_anglais.deja_afaire = {};
      	deja_anglais.deja_auto   = {};
   	end
	end
	deja_anglais.deja_afaire{end+1} = entree;
   %file     = fullfile(root,'anglais','deja_anglais.mat');
   %save(file,'deja_fr','deja_us','deja_date','deja_afaire','deja_auto');

	return
else
	indtexte = strmatch(entree,dico.fr(indexsum),'exact');
	if isempty(indtexte)
		sortie = entree;
		%
		% ajout dans la liste des ressources a traduire
		%
		if isempty(deja_anglais)
   		root     = getappdata(0,'root');
   		file     = fullfile(root,'anglais','deja_anglais.mat');
   		if exist(file,'file') == 2
      		 deja_anglais = load(file);
   		else
     	 		deja_anglais.deja_fr     = {};
      		deja_anglais.deja_us     = {};
      		deja_anglais.deja_date   = {};
      		deja_anglais.deja_afaire = {};
      		deja_anglais.deja_auto   = {};
   		end
		end
		deja_anglais.deja_afaire{end+1} = entree;
		return
	else
      	sortie = dico.us{indexsum(min(indtexte))};
      	return
	end
end
