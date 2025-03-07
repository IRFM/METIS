% copie l'indice de temps indnew de la structure snew dans l'indice de temps indin de la structure sin
function sout = zcopystruct(sin,snew,indin,indnew)

sout = sin;

noms = fieldnames(sin);
for k = 1:length(noms)
	nomc = noms{k};
	if isstruct(snew.(nomc))
		sout.(nomc) = zcopystruct(sin.(nomc),snew.(nomc),indin,indnew);	
	elseif length(size(sout.(nomc))) == 4
		sout.(nomc)(indin,:,:,:) = 	snew.(nomc)(indnew,:,:,:);	
	elseif length(size(sout.(nomc))) == 3
		sout.(nomc)(indin,:,:) = 	snew.(nomc)(indnew,:,:);	
	else
		sout.(nomc)(indin,:) = 	snew.(nomc)(indnew,:);
	end
end

 

