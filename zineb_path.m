% ZINEB_PATH genere le path matlab pour zineb
%---------------------------------------------
% fichier zineb_path.m ->  zineb_path
%
%
% fonction Matlab 5 :
%
% Cette fonction cree ou update le path matlab pour zineb.
% Elle sauve le path utilisateur initiale. 
% Elle gere le path generale de zineb (root/zineb)
% et le path utilisateur de zineb (home/zineb).
% 
% la variable 'root' peut etre donnee comme argument, 
% sinon elle est lue dans la varibale d'environnement
% 'ZINEBROOT'. La valeur par defaut est  :
% '/usr/drfc/cgc'.
%
% syntaxe  :
%  
%       cr =zineb_path({local,root});
%
%
% entrees :
%
%     local     =  repertoire ou racine des repertoires contenants les fonctions de l'utilisateur
%     root      =  repertoire racine de l'arborescence des repertoire de 
%                  zineb ( les fonctions de zineb doivent etre sous 'root/zineb')
%                  
% sorties :
% 
%     cr         =  compte rendu d'execution
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 03/03/2005.
% 
% 
% liste des modifications : 
%
% * 18/09/2001 -> mise en place de la protection pour rigel
% * 25/09/2001 -> ajout du path pour rema
% * 17/10/2001 -> suppression des path html et png
% * 22/02/2002 -> suppression des path specifique rema (mis en local)
% * 28/05/2002 -> ajout du path utilisateur $HOME/zineb
% * 04/06/2002 -> mise en place de la protection pour hercule
% * 13/06/2002 -> gestion du path pour tcronos
% * 17/10/2002 -> ajout du path dynamique pour l'interface en fonction de la langue
% * 07/11/2002 -> positionne la langue par defaut pour utilisation en ligne de commande
% * 19/11/2002 -> addaptation pour linux
% * 17/07/2003 -> ajout de la variable editeur
% * 03/03/2005 -> path cronos en dur
% * 03/03/2005 -> modification de la variable editeur sur linux
% * 15/03/2005 -> suppression des protections hercule/rigel
% * 01/02/2006 -> ajout path mds+ pour zeus
%--------------------------------------------------------------
%
function cr =zineb_path(local,root)


if ~isappdata(0,'langue_cronos')
	setappdata(0,'langue_cronos','anglais');
end



% par defaut
cr =0;

% recherche du chemin vers la racine
if nargin == 0
	root ='';
	local='';
end
if nargin == 1
	root ='';
end

if isempty(root)
	root = getenv('ZINEBROOT');
end
if isempty(root)
        if isdeployed
              root = pwd;
        else
	      root  = fileparts(which('zineb_path'));
        end
end

% try to add python engine 
if exist('pyversion') && ~ispc
    try
        rep = pyversion;
        if isempty(rep)
           %search for python
           [s,t] = unix('which python');
           if s == 0
               pyversion(strtrim(t));
               % test python
               tw = py.textwrap.TextWrapper;
           else
               disp('No python interpreter in path');
           end
        end
    catch
        disp('Python engine is not available');
        disp(lasterr);
    end
end

% % add path to java class for MDS+ for compiled version
% if isdeployed
%     javaaddpath(pwd);
%     try
%         javaaddpath(fullfile(pwd,'java'));
%     catch
%         disp('unable to add path to java class for MDS+');
%         fullpath(pwd,'java')
%         lasterr
%     end
%     try
%         javaaddpath(fullfile(pwd,'javalib'));
%     catch
%         disp('unable to add path to java lib class for MDS+');
%         fullpath(pwd,'java')
%         lasterr
%     end
%     try
%         java.lang.System.load(fullfile(pwd,'javalib','JavaMds.dll'));
%     catch
%         disp('unable to load JavaMds.dll');
%         fullfile(pwd,'javalib','JavaMds.dll')
%         lasterr
%     end
%     try
%         java.lang.System.load(fullfile(pwd,'javalib','JavaMdsLib.dll'));
%     catch
%         disp('unable to load JavaMdsLib.dll');
%         fullfile(pwd,'javalib','JavaMdsLib.dll')
%         lasterr
%     end
% end

% supression des blancs
if ~isdir(root) && ~ispc
    ind=findstr(root,' ');
    if ~isempty(ind)
	    root(ind)=[];
    end
end
%root = strrep(root,'/usr/deneb/gcgc/cgc','/usr/drfc/cgc');

% sauvegarde de root
setappdata(0,'root',root);

% test si il s'agit d'un update du path de zineb
if isappdata(0,'userpath')  && ~isdeployed
	path(getappdata(0,'userpath'));
else
	% sauvegarde du path utilisateur
	setappdata(0,'userpath',path);
end

editeur = '/usr/bin/kate';
envediteur = getenv('EDITOR');

if ~isempty(envediteur) & ~strcmp(computer,'GLNX86')& ~strcmp(computer,'GLNXA64')
   editeur = envediteur;
else
   liste_editeur = {'kate','emacs','nedit','dtpad','kedit','gedit','kwrite','xemacs','textedit','notepad','write'};
   for k=1:length(liste_editeur)
         nomc = liste_editeur{k};
         [s,t]=unix(sprintf('which %s',nomc));
         if (s == 0) & ~isempty(t)
               editeur = t;
               editeur(editeur <= ' ') = [];
               break
         end
   end

end
eddiiid     = getappdata(0,'root');
if strcmp(eddiiid(2:3),'u3')
  setappdata(0,'editeur','gedit');
else
  setappdata(0,'editeur',editeur);
end


if isdeployed
	setappdata(0,'MPIFLAG',0);
else
    filename=fullfile(root,'arch.inc');
    file=textread(filename,'%s','delimiter','\n','whitespace','');

    x=strmatch('MPI=-DMPI',file,'exact');
    if isscalar(x)
	setappdata(0,'MPIFLAG',1);
    else
	setappdata(0,'MPIFLAG',0);
    end
end

if isdeployed
    setappdata(0,'TEMPDIR',tempdir);
    setappdata(0,'TEMPDIR_EXCHANGE',tempdir);
else
    x=strmatch('CRONOSTEMPDIR',file);
    if ~isempty(x)
      [token,remain]=strtok(file{x(1)},'=');
      remain=char(remain);
      remain=deblank(remain(2:end));
      remain(remain <= ' ') = [];	
      if (length(remain)>0)
	    if isdir(remain)
		    setappdata(0,'TEMPDIR',remain);
	    else
		    setappdata(0,'TEMPDIR','DirNoSpecified');
	    end
      else
	  setappdata(0,'TEMPDIR','DirNoSpecified');
      end
    else
      setappdata(0,'TEMPDIR','DirNoSpecified');
    end

    x=strmatch('CRONOSTEMPDIR_EXCHANGE',file);
    if ~isempty(x)
      [token,remain]=strtok(file{x},'=');
      remain=char(remain);
      remain=deblank(remain(2:end));
      remain(remain <= ' ') = [];	
      if (length(remain)>0)
	    if isdir(remain)
		    setappdata(0,'TEMPDIR_EXCHANGE',remain);
	    else
		    setappdata(0,'TEMPDIR_EXCHANGE','DirNoSpecified');
	    end
      else
	  setappdata(0,'TEMPDIR_EXCHANGE','DirNoSpecified');
      end
    else
      setappdata(0,'TEMPDIR_EXCHANGE','DirNoSpecified');
    end
end

% matlabpool
setappdata(0,'matlabpool_onoff',0);
% if verLessThan('matlab', '7.5.0')
%     % rien
%     setappdata(0,'matlabpool_onoff',0);
% elseif exist('matlabpool')
%     try
%          if matlabpool('size') == 0
%                 matlabpool open local    
%                 setappdata(0,'matlabpool_onoff',1);
%          end
%     catch
%     	setappdata(0,'matlabpool_onoff',0);
%     end
% end

% switch off legend automatic update
try
    set(0,'defaultlegendAutoUpdate','off');
catch
    % nothing
end

if isdeployed
    % pour la version compilee
    set(0,'defaultFigureToolBar','figure');
    % pas d'autre action si c'est une application standalone
    return
end

% creation du path zineb de base
warning off
	addpath(root);
warning on
% path en dur pour question de rapidite
addpath(fullfile(root,'acces/cronos'));
addpath(fullfile(root,'acces/imas'));
addpath(fullfile(root,'acces/imas/noimas_installed'),'-end');
addpath(fullfile(root,'acces/iter'));
addpath(fullfile(root,'acces/jet'));
addpath(fullfile(root,'acces/jet/sal'));
addpath(fullfile(root,'acces/ts'));
addpath(fullfile(root,'acces/ST40'));
addpath(fullfile(root,'acces/west'));
addpath(fullfile(root,'acces/JT-60SA'));
addpath(fullfile(root,'acces/tcv'));
addpath(fullfile(root,'certification'));
addpath(fullfile(root,'coef'));
addpath(fullfile(root,'fbe'));
addpath(fullfile(root,'graphe'));
addpath(fullfile(root,'import'));
addpath(fullfile(root,'import/chaine'));
addpath(fullfile(root,'import/divers'));
addpath(fullfile(root,'import/sampling'));
addpath(fullfile(root,'import/splash'));
addpath(fullfile(root,'import/tools'));
addpath(fullfile(root,'interface'));
addpath(fullfile(root,'init'));
addpath(fullfile(root,'op'));
addpath(fullfile(root,'op/polyroot'));
addpath(fullfile(root,'op/m2html'),'-end');
addpath(fullfile(root,'op/StructBrows/StructBrowser'));  % probleme graphique sous matlab 7 
% conflit avec les fonction elementaire matlab (abs, atan2, asin , erf ...)
addpath(fullfile(root,'op/psitbx'),'-end');
addpath(fullfile(root,'op/xml_io_tools'));  % probleme graphique sous matlab 7 
addpath(fullfile(root,'op/xml'));  % probleme graphique sous matlab 7 
addpath(fullfile(root,'op/xml/XML4MAT'));% % probleme graphique sous matlab 7
addpath(fullfile(root,'op/xml/xmltree'));%probleme graphique sous matlab 7 
addpath(fullfile(root,'op/expokit/matlab'),'-end');
addpath(fullfile(root,'simulink'));
addpath(fullfile(root,'solver'));
addpath(fullfile(root,'source/fus'));
addpath(fullfile(root,'source/icrh'));
if isdir(fullfile(root,'source/lh'))
	addpath(fullfile(root,'source/lh'));
end
addpath(fullfile(root,'test'));
addpath(fullfile(root,'trait'));
addpath(fullfile(root,'zerod'));


% if isdeployed 
%   warning off
%   splash('Metis_splash_screen','png');
%   drawnow
% end

% test availiability of signal toolbox
if exist('medfilt1') == 2
	try
		void = medfilt1(rand(1,11),3);
	catch
		addpath(fullfile(root,'compatibility','signal'),'-BEGIN');	
	end
else
	addpath(fullfile(root,'compatibility','signal'),'-BEGIN');
end

%
% --------- LUKE local and remote distributions ---------
%
if isdir('/Applications/LUKE') && isempty(which('startup_LUKE'))
    addpath('/Applications/LUKE')
    addpath('/Applications/LUKE/LUKE_LHCD')
    addpath(genpath(['/Applications/LUKE/ALOHA/libaloha'])); 
    
    if ~isempty(which('startup_LUKE'))
        disp('=========================================================')
        disp('Connnection of METIS and LUKE');
        disp('=========================================================')
        cdmem_ = pwd;
        try
            startup_LUKE;
        catch
            disp('An error appears during the startup of LUKE setting: LUKE is not yet available');
        end
        try
            if isdir(fullfile(root,'op/temp_mfile'))
                addpath(fullfile(root,'op/temp_mfile'));
            end
            select_LUKE;
            if isdir(fullfile(root,'op/temp_mfile'))
                rmpath(fullfile(root,'op/temp_mfile'));
            end
        catch
            try
                if isdir(fullfile(root,'op/temp_mfile'))
                    rmpath(fullfile(root,'op/temp_mfile'));
                end
            end
            disp('An error appears during the slection of LUKE version: LUKE version must be selected manually');
        end
        cd(cdmem_);
        disp('=========================================================')
        disp('Now starting METIS Graphical User Interface');
        disp('=========================================================')
        
    end
end

% gestion du path utilisateur de zineb
% test si la session d'installation de zineb est differente de la session utilisateur
home=getenv('HOME');
ll = 1:min(length(home),length(root));
if strcmp(home(ll),root(ll))
	return
end
% gestion utilisateur
if strcmp(getenv('USER'),'cgc') | strcmp(getenv('USER'),'trait') | strcmp(getenv('USER'),'devtrait')
	return
end

% repertoire local par defaut
if isempty(local)
    % gestion des directories utilisateur
    if ~exist(strcat(home,'/zineb'),'dir')
            [cr,text]=mkdir(home,'zineb');
            if cr==0
               disp('Erreur lors de la creation du repertoire utilisateur zineb:')
                    disp(text)
               path(getappdata(0,'userpath'));
               rmappdata(0,'userpath');
               disp('Path utilisateur restaure');
                    return
            end
    end
    local = strcat(home,'/zineb');
end

if ~exist(strcat(local,'/dynamique'),'dir')
	[cr,text]=mkdir(local,'dynamique');
	if cr==0
	   disp('Erreur lors de la creation du repertoire utilisateur pour les fonctions generes par zineb:')
		disp(text)
	   path(getappdata(0,'userpath'));
	   rmappdata(0,'userpath');
	   disp('Path utilisateur restaure');
		return
	end
end

% donnee appdata
setappdata(0,'local',local);

% creation du path zineb de lutilisateur
%  [cr,liste] = unix(['ls -FRL1 ',home, ...
%                       '/zineb | grep "/" | grep ":" ',...
%                       '| grep -v "private" | grep -v "@" | grep -v "maple"']);
%  if cr~=0
%  	disp('Erreur lors de la lecture du path de zineb (utilisateur):')
%  	disp(liste)
%  	path(getappdata(0,'userpath'));
%  	rmappdata(0,'userpath');
%  	disp('Path utilisateur restaure');
%  	return
%  end
%  ind=findstr(liste,sprintf('\n'));
%  liste(ind)=[];
%path(liste,path);
path(local,path);
path(fullfile(local,'dynamique'),path);

setappdata(0,'path_zineb',1);




