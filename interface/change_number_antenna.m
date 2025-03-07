% script de changement de nombre d'antenne, ...
try 
	nb = param.nombre;
	cons = data.cons;
catch
	error('Matlab workspace must content a valid CRONOS data set');
end

% choix de la source
choice = menu('Select a source',{'ICRH','ECRH','LHCD','NBI','Pellet'});
number = -1;
while (~isfinite(number) | (number < 1))
	rep = inputdlg('number of antenna or injector ?');
	number = str2num(rep{1});
end

switch choice
case 1'
	field = 'fci';          
	alt   = 'nbfci';
case 2'
	field = 'fce';
        alt   = 'nbfce';
case 3
	field = 'hyb';
        alt   = 'nbhyb';
case 4
	field = 'idn';
        alt   = 'nbidn';
case 5
	field = 'glacon';
        alt   = 'nbglacon';
otherwise
	error('???');
end

%  pour impose la sauvegarde
zuisavenonok

% changement des nombres declares
memnb                = param.nombre.(field);
param.nombre.(field) = number;
param.gene.(alt)     = number;

% changement du parametrage 
info   = feval(param.fonction.(field),number);
valeur = info.valeur;
old    = param.cons.(field);
noms   = fieldnames(valeur);
for k=1:length(noms)
	if isnumeric(valeur.(noms{k}))
		if length(valeur.(noms{k})) ~= length(old.(noms{k}))
			if number > memnb
				valeur.(noms{k})(1:memnb) = old.(noms{k})(1:memnb);
				valeur.(noms{k})((memnb+1):number) = old.(noms{k})(memnb);
			elseif number < memnb
				valeur.(noms{k})(1:number) = old.(noms{k})(1:number);
			else
				valeur.(noms{k}) = old.(noms{k});
			end
		else
		    valeur.(noms{k}) = old.(noms{k});
		end
	elseif iscell(valeur.(noms{k}))
		if length(valeur.(noms{k})) ~= length(old.(noms{k}))
			if number > memnb
				valeur.(noms{k})(1:memnb) = old.(noms{k})(1:memnb);
				for l = (memnb+1):number
					valeur.(noms{k}){l} = old.(noms{k}){memnb};
				end
			elseif number < memnb
				for l = 1:number
					valeur.(noms{k}){l} = old.(noms{k}){l};
				end
			else
				valeur.(noms{k}) = old.(noms{k});
			end
		else
		    valeur.(noms{k}) = old.(noms{k});
		end
	elseif ischar(valeur.(noms{k}))
		if ~isempty(str2num(old.(noms{k})))
			res = str2num(old.(noms{k}));
			if length(res) == 1
 				valeur.(noms{k}) = old.(noms{k});
			elseif length(res) == memnb
				if number > memnb
					res_e = res;
					res_e((memnb+1):number) = res(end);
					valeur.(noms{k})(1:memnb) = num2str(res_e);
				elseif number < memnb
					res_e = res(1:number);
					valeur.(noms{k})(1:memnb) = num2str(res_e);
				else
					valeur.(noms{k}) = old.(noms{k});
				end
			else
 				valeur.(noms{k}) = old.(noms{k});
			end
		else
		    valeur.(noms{k}) = old.(noms{k});
		end
	else
		    valeur.(noms{k}) = old.(noms{k});
	end
end
param.cons.(field) = valeur;

% changement dans les consignes
% on conserve les valeurs existantes
if number > memnb
	data.cons.(field)(:,(memnb+1):number) = data.cons.(field)(:,end) * ones(1,number - memnb);
elseif number < memnb
	data.cons.(field)(:,number:end) = [];
end






