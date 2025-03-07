% recupere les parametre d'un module pour les mettre dans un nouveau
%  cons  = parametres definis precedemment
%  nomf = nom du module module
%  nb   = nombre de coupleurs ...
function [par,nbok] = zrepicparam(cons,nomf,nb)

% valeur de depart
nbok   = 0;
par    = cons;

% lecture des donnees
par    = feval(nomf,nb);
if isfield(par,'valeur');
	valeur = par.valeur;
	if isempty(valeur)
	   return
	end
else
	% fonction sans parametre
	return
end

% si pas de parametres existant 
if isempty(cons)
	 par = valeur;
end

% boucle sur les champs
noms = fieldnames(valeur);
for k = 1:length(noms)
	nomc = noms{k};
	% si le champ existe deja
	if isfield(cons,nomc)
		valcons = getfield(cons,nomc);
		valval  = getfield(valeur,nomc);
		% si le type est commun
		if (iscell(valval) & iscell(valcons)) | ...
			(ischar(valval) & ischar(valcons)) | ...
			(isnumeric(valval) & isnumeric(valcons)) 
			% si memes dimensions
			if all(size(valval) == size(valcons)) & isfield(par.borne,nomc)
				borne  = getfield(par.borne,nomc);
				% si intervalle compatible
				ok = inborne(valcons,borne);
				if any(ok)
					valval(find(ok)) = valcons(find(ok));
					eval(sprintf('valeur.%s = valval;',nomc));
					nbok = nbok + 1;
				end
			end
		end
	end
end

% sortie
par = valeur;

function ok = inborne(val,borne)

% code retour
ok = ones(1,length(val));

for k = 1:length(val)
	% selon le type
	if iscell(val)
		vc  = val{k};		
		if iscell(borne)
			if ischar(vc)
				if ~any(strmatch(vc,borne,'exact'))
					ok(k) = 0;
				end
			else
				borne = zcell2mat(borne);
				if ~any(vc == borne)
					ok(k) = 0;
				end
			end	
		else
			if vc < min(borne)
				ok(k) = 0;
			end
			if vc > max(borne)
				ok(k) = 0;
			end
		end
	elseif isnumeric(val)
		vc  = val(k);		
		if iscell(borne)
			borne = zcell2mat(borne);
			if ~any(vc == borne)
				ok(k) = 0;
			end
		else
			if vc < min(borne)
				ok(k) = 0;
			end
			if vc > max(borne)
				ok(k) = 0;
			end
		end
	else
		ok(k) = 0;
	end
end
