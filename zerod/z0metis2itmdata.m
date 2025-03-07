% fonction creant la liste des info pour les donnees 0D de metis dans l'ITM
function z0metis2itmdata

% acces a la liste des donnees metis
list_metis = zero1t;

% exemple 
pf= fullfile(fileparts(which('metis')),'certification','metis');
lf = dir(fullfile(pf,'*.mat'));
data = load(fullfile(pf,lf(1).name));
data = data.post;

% acces aux donnees de conversion deja cree
try
	tree = load('metis2itm0d.xml'); 
	list_itm = xml2struct(tree);
catch
	list_itm =[];
end

noms = fieldnames(list_metis);
for k=1:length(noms)
	nomc = noms{k};
	if ~isfield(list_itm,nomc)
		list_itm.(nomc).itm_name = nomc;	
	end
	info = list_metis.(nomc);
	[com,unit]  = strtok(info,'(');
	[unit,void] = strtok(unit,')');
	if ~isfield(list_itm.(nomc),'comment')
		list_itm.(nomc).('comment') = com;
	end
	if ~isfield(list_itm.(nomc),'unit') 
		if ~isempty(unit)
			list_itm.(nomc).('unit') = unit(2:end);
		else
			list_itm.(nomc).('unit') = '';
		end
	end
	if ~isfield(list_itm.(nomc),'type')
		if length(data.zerod.(nomc)) > 1
			list_itm.(nomc).('type') = 'time_dependant';
		else
			list_itm.(nomc).('type') = 'scalar';
		end
	end 
	
end

tree = struct2xml(list_itm);
save(tree,'metis2itm0d.xml')
