% ZUISAVEDATA sauve les donnees de zineb
%--------------------------------------------------------------------------------------
% fichier zuisavedata.m ->  zuisavedata
%
% fonction Matlab 5 :r
%	Cette fonction sauve les donnees de zineb . 
%
% syntaxe  :
%	zuisavedata;
%
% entrees
%
% sorties
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 14/04/2005.
% 
% liste des modifications : 
%   * 25/06/2001 -> ajout de la compression
%   * 29/08/2001 -> ajout selected(off) sur le botton sauve
%   * 16/10/2001 -> ajout de la gestion du commentaire
%   * 06/11/2001 -> ajout de la question si deja sauvegarder
%   * 11/03/2002 -> modification de la gestion de "annulation" de la boite de dialogue pour les commentaires
%   * 25/03/2002 -> gestion du type de fichier source ou resultat
%   * 09/07/2002 -> ajout de la gestion automatique de rapauve
%   * 14/12/2004 -> ajout de la sauvegarde des donnees simulink
%   * 14/04/2005 -> sauvegardes en option '-V6' pour pouvoir recharger les fichiers sous Matlab5
%   * 03/05/2005 -> compatibilite multi matlab
%--------------------------------------------------------------
%
function zuisavedata(force)

% verification de l'existance des donnees
if ~(existbase('param')==1) | ~(existbase('data')==1)
	warning('no data to be saved')
	return
end
if ~(existbase('post')==1)
     assignin('base','post',[]);
end

% test du commentaire
com = evalin('base','param.from.creation.com');
if isempty(com)
	%  fenetre pour ajouter le commentaire
	user = evalin('base','param.from.creation.user');
	user = strrep(user(:)',sprintf('\n'),'');
	% dialogue
   prompt={'enter your name (ou user):','enter your comments :'};
   def={user,com};
   dlgTitle='comments link to a Cronos file';
   lineNo=2;
   ok = 0;
   while (ok == 0)
   	answer=zinputdlg(prompt,dlgTitle,lineNo,def);
   	if ~isempty(answer)
   		user = answer{1};
   		com  = answer{2};
   		if ~isempty(user) & ~isempty(com)
   			if ~all(user <= sprintf(' ')) & ~all(com <= sprintf(' '))
   				ok = 1;
   				zassignin('base','param.from.creation.com',com);
   				zassignin('base','param.from.creation.user',user);
    			end
   		end
   	   def={user,com};
	   else
	      warning('Canceled save')
	      return
   	end
   end
end

% gestion du model simulink
try
         onoffmdl = evalin('base','param.simu.onoff');
         if onoffmdl== 1
            evalin('base','param.simu = zuploadmdl(param.simu,0);');
         end
catch
 	warning('This file may be updated : structure for simulink missing')
end

% recupere le nom du fichier de sauvegarde
filename = evalin('base','param.edit.currentfile','[]');
if ~isempty(filename)
    % supression du .gz
    [pf,nf,ext,ver]=fileparts(filename);
    filename = fullfile(pf,nf);
    % supression du .mat
    [pf,nf,ext,ver]=fileparts(filename);
    filename = fullfile(pf,nf);
end
origine  = evalin('base','param.gene.origine','zineb*.mat.gz');
fileres  = evalin('base','param.gene.file','zineb*_resultat.mat.gz');
filetype = evalin('base','param.gene.filetype','source');
rapsauve = evalin('base','param.gene.rapsauve','');
saveok   = evalin('base','param.edit.saveok','0');

% choix du fichier si vide + changement de nom
if isempty(filename)
        cdmem = cd;
	if strcmp(filetype,'source');
	      [pf,nf,ext,ver]=fileparts(origine);
	      if isdir(pf)
	      	cd(pf);
	      end
	      [file,path]=uiputfile({'*.mat';'*.gz'},'Cronos datafile name ?');
	else
	      [pf,nf,ext,ver]=fileparts(fileres);
	      if isdir(pf)
	      	cd(pf);
	      end
	      [file,path]=uiputfile({'*.mat';'*.gz'},'Cronos datafile name ?');
   	end
	drawnow
	cd(cdmem);
	if ~ischar(file)
		warning('canceled save')
		return
	end
	filename = strcat(path,file);
	filename = strrep(filename,'/tmp_mnt','');
   % supression du .gz
   [pf,nf,ext,ver]=fileparts(filename);
   filename = fullfile(pf,nf);
   % supression du .mat
   [pf,nf,ext,ver]=fileparts(filename);
   filename = fullfile(pf,nf);
	% modification de la variable edit.currentfile
	zassignin('base','param.edit.currentfile',filename);
	% modification de la variable origine du fichier
	if strcmp(filetype,'source')
	    zassignin('base','param.gene.origine',filename);
	    % modification du champ file
       filec = strcat(filename,'_resultat');
	    zassignin('base','param.gene.file',filec);
	    % modification du champ rapsauve
		 rapsauve = creerapsauve(rapsauve,filename,filetype);
		 zassignin('base','param.gene.rapsauve',rapsauve);
	else
	    % fichier resultat (origine ne change pas)
	    zassignin('base','param.gene.file',filename);
	    % modification du champ rapsauve
		 rapsauve = creerapsauve(rapsauve,filename,filetype);
		 zassignin('base','param.gene.rapsauve',rapsauve);
	 end
else
	% securite anti resauvegarde
	filec = strcat(filename,'.mat');

	% si deja sauver information
	if (saveok == 1) & (exist(filec) == 2)
	     if nargin  == 0
			warndlg('data already saved !','Save :');
			return
	     else
			rep = questdlg('do you want to save again the datas ?','already saved data !', ...
			'Yes','No','No');
			if strcmp(rep,'No')
				return
	                end
	     end
	end
	filec = strcat(filename,'.mat.gz');

	% si deja sauver information
	if (saveok == 1) & (exist(filec) == 2)
	     if nargin  == 0
			warndlg('data already saved !','Save :');
			return
	     else
			rep = questdlg('do you want to save again the datas ?','already saved data !', ...
			'Yes','No','No');
			if strcmp(rep,'No')
				return
	                end
	     end
	end
end

% sauvegarde
% compactage des donnees
evalin('base','data=zreduit(param,data,''compact'');');
% sauvegarde
vermat  = version;
if str2num(vermat(1)) < 7
   evalin('base',strcat('savedb(''',filename,''',''param'',''data'',''post'')'));
   % compression du fichier
   zgzip(filename,'compress');

elseif strmatch('saturne',getenv('HOSTNAME'))
	if isappdata(0,'pid')
		pid = getappdata(0,'pid');
        else
		pid = [];
	end
	if isempty(pid)
		[pid,pid_sess,user]=getidprocess;
		setappdata(0,'pid',pid);
	end
	file_scratch = fullfile('/scratch',sprintf('zineb%d',pid));
        evalin('base',strcat('savedb(''',file_scratch,''',''param'',''data'',''post'',''-V6'')'));
        % compression du fichier
        zgzip(filename,'compress',file_scratch);
else
   evalin('base',strcat('savedb(''',filename,''',''param'',''data'',''post'')'));
end

% decompactage des donnees
evalin('base','data=zreduit(param,data,''uncompact'');');
% save ok
zassignin('base','param.edit.saveok',1);
[hfig,h] = zuiformhandle('direct');
if ishandle(hfig)
	set(h.radio_savefile,'foregroundcolor',[0 0.5 0]);
end

% ajout du repertoire de travail  
evalin('base','setappdata(0,''CRONOS_WORK_DIR'',fileparts(param.gene.file))');

% test l'existance d'une variable dans le workspace
function cr = existbase(nom)

cr = evalin('base',strcat('exist(''',nom,''')'));

% genere la variable rapsauve
function rapsauve = creerapsauve(rapsauve,filename,filetype)

% si vide pas de sauvegarde rapide
if isempty(rapsauve)
    return
end

% test si installation standard
[dirpath,name] = fileparts(filename);
std            = fullfile(dirpath,'rapsauve');
if exist(std,'dir')
	racine = std;
else
	[racine,void] = fileparts(filename);
end

if strcmp(filetype,'source')
    name = strcat(name,'_resultat');
end

rapsauve = fullfile(racine,name);

