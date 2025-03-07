function metis_save(filename,post)

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

if nargin < 2
  try
	post = evalin('base','post');
  catch
	post = [];
  end
end
langue      =  lower(getappdata(0,'langue_cronos'));
if isempty(post)
	switch langue
	case 'francais'
		warning('Pas de donnee a sauver')
	otherwise
		warning('No data to be saved');
	end
	return
end
if ~isfield(post,'zerod')
	switch langue
	case 'francais'
		warning('Pas de donnee a sauver')
	otherwise
		warning('No data to be saved');
	end
	return
end
if ~isfield(post,'z0dinput')
	switch langue
	case 'francais'
		warning('Pas de donnee a sauver')
	otherwise
		warning('No data to be saved');
	end
	return
end
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

% sauvegarde  des donnees de la simulation en mode evolution
try
	z0dstruct = evalin('base','z0dstruct');
catch
	z0dstruct = [];
end
zassignin('base','post.z0dstruct_inter',z0dstruct);
data.z0dstruct_inter = z0dstruct;

% donnees LUKE  dans METIS
if isfield(post,'lukeinmetis')
	data.lukeinmetis = post.lukeinmetis;
end   


% ajout des donnees de certication
root = getappdata(0,'root');
if isempty(root)
	zineb_path;
	root = getappdata(0,'root');
end
% structure d'information 
data.info_test_metis.date        = clock;
[s,t] = unix('uname -a');
if s == 0
	data.info_test_metis.machine = t;
else
	error(sprintf('error executing ''uname -a'' (%s)',t));
end
[s,t] = unix('env');
if s == 0
	data.info_test_metis.env = t;
else
	error(sprintf('error reading environnement variables  (%s)',t));
end
[data.info_test_metis.version,data.info_test_metis.date_version]  = zinebversion;
data.info_test_metis.matlab_version  = version;
data.info_test_metis.toolbox_version = ver;
data.info_test_metis.path = matlabpath;
data.info_test_metis.user = getenv('USER');
data.info_test_metis.root = root;
data.info_test_metis.mexsolver = mexsolver;

% save cooling ratse to be able to know which one has been used for the
% simulation.
try
     load('Lz_zave.mat','tabmat')
     data.tabmat = tabmat;
catch
    disp('metis_save: unable to read cooling rate');
end

post = data;
clear data
save(filename,'post');
fprintf('METIS data file : %s saved\n',filename);
z0dinterfacetitle(filename);
