% tri des structure et affichage dans l'ordre
function tristruct(st)

names = fieldnames(st);
[names,indx] = sort(names);

for k=1:length(names)
	nom = names{k};
	txt = getfield(st,nom);
	fprintf('%s : %s\n',nom,txt);
end