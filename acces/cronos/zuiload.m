% ZUILOAD  charge un fichier de donnees
%-------------------------------------
% fichier zuiload.m ->  zuiload
%
% fonction Matlab 5 :
% Cette fonction  charge un fichier de donnees. 
%
% syntaxe  :
%	[argument_1,param,post]=zuiload(fichier)
%
%  charge les donnees dans le workspace :
%     zuiload            -> utilise uigetfile
%     zuiload fichier
%  
%  retourne les donnees :
%     [data,param,post]=zuiload;          -> utilise uigetfile
%     [data,param,post]=zuiload(fichier);
%  
%  retourne les donnees dans une structure :
%     jeux1=zuiload;          -> utilise uigetfile
%     jeux1=zuiload(fichier);
%
% entrees : 
%	variables suivant l'appel, voir la syntaxe
%
% sorties : 
%	variables suivant l'appel, voir la syntaxe
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 07/06/2005.
% 
% 
% liste des modifications : 
%
%   * 13/06/2001 -> ajout de la compression/decompression
%   * 18/07/2001 -> compatibilite des sorties avec le load
%   * 14/09/2001 -> ajout de la permutation des donnees (jeux1)
%   * 14/09/2001 -> ajout de la coherence des nom de fichiers
%   * 14/09/2001 -> ajout du warning sur la version
%   * 24/09/2001 -> modification du chemin par defaut + correction bug matlab
%   * 28/09/2001 -> le mode avec argument de sortie ne pose pas de question
%   * 02/10/2001 -> correction bug de gestion du chemin par defaut selon la session
%   * 09/11/2001 -> ajout de la securite pour le cas ou la structure param.edit n'existe pas
%   * 14/12/2004 -> ajout de la gestion des modele simulink
%   * 17/01/2005 -> remplace rm par rm -f
%   * 07/06/2005 -> correction chemin initial pour fichiers
%
%
%--------------------------------------------------------------
%
function [argument_1,param,post]=zuiload(fichier)

% test des entrees
if nargin < 1
	fichier = '' ;
end

% sortie
argument_1 = [];
param      = [];
post       = [];
data       = [];

% recupere le nom du fichier de sauvegarde
try
    filename = evalin('base','param.edit.currentfile','[]');
catch
    filename ='';
end
try
  saveok   = evalin('base','param.edit.saveok','0');
catch
    saveok = 0;
end
	
% securite anti ecrasement
% Changement pour eviter d'avoir a repondre au 1er chargement
% verification de l'existance des donnees
if (existbase('param')==1) & (existbase('data')==1) & (nargout < 1)
	if saveok == 0
		rep =questdlg('Do you want to save the data ?', ...
		              'Be careful -> unsaved data !', ...
		              'Yes','No','Cancel','Cancel');
		switch rep
		case 'Yes'
			disp('save in progress')
			zuisavedata;
		case 'Cancel'
			return
		case 'No'
			disp('data removed from workspace')
		end
	end
	
	% permutation avec jeux1
	rep =questdlg('Do you want to store the present data as the reference dataset ?', ...
		              'Reference dataset', ...
		              'Yes','No','Cancel','No');
	switch rep
	case 'Yes'
		disp('in progress')
		evalin('base','jeux1.data = data;','jeux1.data=[];');
		evalin('base','jeux1.param = param;','jeux1.param=[];');
		evalin('base','jeux1.post = post;','jeux1.post=[];');
		evalin('base','t1 = jeux1.data.gene.temps;','t1=[];');
	case 'Cancel'
		return
	case 'No'
		% rien
	end
	
end
if isempty(fichier)
	% extrait le chemin
	if ~isempty(filename)
	        [pp,ff] = fileparts(filename);
		if ~isempty(pp)
			name = pp;
		else
           		 if strcmp(getenv('USER'),'cgc')
				 %name = fullfile(getenv('HOME'),'matlab5','zineb','data'); %mis en commentaire FI 19/02/03
                 %name = '/usr/drfc/zineb/data';
                  name = '/home/sccp/gsem/cgc/zineb/data';
                 [s,t] = unix(sprintf('ls -a %s',name));
                 if s ~= 0
                    name = '/home/sccp/gsem/cgc/cgc_data/zineb/data';
                 end
	    		else
				 name = fullfile(getenv('HOME'),'zineb','data');
	    		end
		end
	else
            if strcmp(getenv('USER'),'cgc')
		     %name = fullfile(getenv('HOME'),'matlab5','zineb','data'); %mis en commentaire FI 19/02/03
                     %name = '/usr/drfc/zineb/data';
                     %name = '/home/sccp/gsem/cgc/zineb/data';
		     %name = '/home/sccp/gsem/cgc/cgc_data/zineb/data';
                     name = fullfile(getenv('HOME'),'zineb','data');
	    else
		 name = fullfile(getenv('HOME'),'zineb','data');
	    end
	end
	% dialogue
	pwd_mem = pwd;
	if isdir(name)
		cd(name);
	end
	[file,path]=uigetfile({'*.mat';'*.gz'},'file name to load ?');
	if isdir(pwd_mem)
		cd(pwd_mem);
	end
	drawnow
	if ~ischar(file)
		warning('loading canceled')
		return
	end
	if strcmp(fileparts(path),'cgc_data')    % nouveau bug apparu et corrige le 02/05/2005 F.I.
	   path=['/home/sccp/gsem/cgc',path];
	end
	filename = strcat(path,file);
	filename = strrep(filename,'/tmp_mnt','');
%	filename = strrep(filename,'/usr/deneb/cgc','/usr/drfc/cgc');
%	filename = strrep(filename,'/usr/deneb/gcgc/cgc','/usr/drfc/cgc');
	filename = strrep(filename,'/usr/deneb/cgc','/home/sccp/gsem/cgc');
	filename = strrep(filename,'/usr/deneb/gcgc/cgc','/home/sccp/gsem/cgc');
else
	filename = fichier;
end

% selon le nombre d'arguments
if nargout == 0
	% chargement
	[lfile,rmfile,cr]=zgzip(filename,'uncompress');
	evalin('base',strcat('load(''',lfile,''')')) ;
	if ~isempty(rmfile)
		[voids,voidt]=unix(['rm -f ',rmfile,' >& /dev/null']);
	end
	%evalin('base',strcat('load(''',filename,''')'));
	evalin('base','data = zreduit(param,data,''uncompact'');');

	% modifie les champs
	zassignin('base','param.edit.currentfile',filename);
	zassignin('base','param.edit.saveok',1);
	zassignin('base','void',[]);
	clear argument_1

	% verification de la coherence des noms
	cr = zveriffilename(filename);

	
	% verification de la coherence des versions
	zverifversion; 

	% check if the simulation was successfull (for a result file)
	zverifexecution; 

        % gestion du model simulink 
        simuok = evalin('base',['isfield(param,''simu'')']);
        if simuok == 1
               evalin('base','[crsimu,param.simu] = zdownloadmdl(param.simu);');
        end
	
	% ajout du repertoire de travail  
	evalin('base','setappdata(0,''CRONOS_WORK_DIR'',fileparts(param.gene.file))');

	
elseif nargout  == 1
	% chargement
	[lfile,rmfile,cr]=zgzip(filename,'uncompress');
	load(lfile) ;
	if ~isempty(rmfile)
		[voids,voidt]=unix(['rm -f ',rmfile,' >& /dev/null']);
	end
	%evalin('base',strcat('load(''',filename,''')'));
	data = zreduit(param,data,'uncompact');
	
	% modifie les champs
	param.edit.currentfile=filename;
	param.edit.saveok=1;

        % gestion du model simulink 
        simuok = isfield(param,'simu');
        if simuok == 1
              [crsimu,param.simu] = zdownloadmdl(param.simu);
        end
	
	% ajout du repertoire de travail  
	setappdata(0,'CRONOS_WORK_DIR',fileparts(param.gene.file));

        % affecte les sorties
	argument_1.param = param;
	argument_1.data = data;
	argument_1.post = post;
	
else
	% chargement
	[lfile,rmfile,cr]=zgzip(filename,'uncompress');
	load(lfile) ;
	if ~isempty(rmfile)
		[voids,voidt]=unix(['rm -f ',rmfile,' >& /dev/null']);
	end
	%evalin('base',strcat('load(''',filename,''')'));
	data = zreduit(param,data,'uncompact');
	
        % gestion du model simulink 
        simuok = isfield(param,'simu');
        if simuok == 1
              [crsimu,param.simu] = zdownloadmdl(param.simu);
        end

	% modifie les champs
	param.edit.currentfile=filename;
	param.edit.saveok=1;
	argument_1 = data;
	% ajout du repertoire de travail  
	setappdata(0,'CRONOS_WORK_DIR',fileparts(param.gene.file));
end

zassignin('base','param.edit.saveok',1);
[hfig,h] = zuiformhandle('direct');
if ishandle(hfig)
	set(h.radio_savefile,'foregroundcolor',[0 0.5 0]);
end


% verification de la coherence des noms
%cr = zveriffilename(filename);


% verifcation de la coherence des versions
%zverifversion; 

% fonction prise dans zuisavedata
% par C. Passeron le 29/05/2001
function cr = existbase(nom)

cr = evalin('base',strcat('exist(''',nom,''')'));

