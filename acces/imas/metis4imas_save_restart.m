function metis4imas_save_restart(filename,post,diary_text)

% variable global importante pour les test
% utilisation du mexfile si disponible	 
if isappdata(0,'MEXSOLVER_IN_METIS')
	mexsolver = getappdata(0,'MEXSOLVER_IN_METIS');
else	
	mexsolver =  [];
end
if isempty(mexsolver)
	repmex=which(strcat('mexpde1dsolver.',mexext));
	if ~isempty(repmex)
		mexsolver = 1;
	else
		mexsolver = 0;	
	end
	setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
end  

langue      =  lower(getappdata(0,'langue_cronos'));
if isempty(post)
	warning('No data to be saved');
	return
end
if ~isfield(post,'zerod')
	warning('No data to be saved');
	return
end
if ~isfield(post,'z0dinput')
	warning('No data to be saved');
	return
end
post.profil0d = post.profil;
data.zerod    = post.zerod;
data.z0dinput = post.z0dinput;
data.profil0d = post.profil0d;
% sauvegarde des donnees externes
noms = fieldnames(getappdata(0));
for k = 1:length(noms)
	if findstr(noms{k},'_EXP')
		data.appdata.(noms{k}) = getappdata(0,noms{k});
	end
end

% sauvegarde de la sortie standard du model simulink
try
	simout = evalin('base','simout');
catch
	simout = [];
end
zassignin('base','post.simout',simout);
data.simout = simout;

% donnees LUKE  dans METIS
if isfield(post,'lukeinmetis')
	data.lukeinmetis = post.lukeinmetis;
end   


% ajout des donnees de certication
root = getappdata(0,'root');
post = data;
clear data
save(filename,'post','diary_text');
