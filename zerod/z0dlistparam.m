% cree le fichier ascii des parametres de reference qui peu etre charge pour ecrase les parametres de la simulation
function z0dlistparam(post)

% sauvegarde 
fn =   sprintf('METIS_reference_parameter_%s@%d.txt',post.z0dinput.option.machine,fix(post.z0dinput.option.shot));

% ouverture du fichier 
[fid,mess] = fopen(fn,'w');
if fid < 0
	error(sprintf('Z0DLISTPARAM : %s (%s)',mess,fn));
end

blacklist = {'machine','shot'};

noms = fieldnames(post.z0dinput.option);
for k=1:length(noms)
	if isempty(strmatch(noms{k},blacklist,'exact'))
	    val = post.z0dinput.option.(noms{k});
	    if ischar(val)
		    fprintf(fid,'%s = %s\n',noms{k},val);
	    else
		    fprintf(fid,'%s = %g\n',noms{k},val);
	    end
        end
end
fclose(fid);
fprintf('METIS reference parameters file created : %s\n',fn);
