% ZUIEDIT_IMPORT_MODE formulaire d'importation de donnees zineb
%--------------------------------------------------------------
% fichier zuiedit_import_mode.m ->  zuiedit_import_mode
%		           zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	formulaire d'importation de donnees zineb
%	Il permet d'importer dans la structure de donn�s zineb des donn�s
%	provennant d'un fichier ou de l'espace de travail.
%	Il a 2 modes de fonctionnement:
%		- il importe des consignes
%		- il importe des profils
%
% syntaxe  :
%	hout=zuiedit_import_mode(nom_mode,mode_fct);
%
% entree :
%	nom_mode = nom de la variable import�
%	mode_fct = 'consigne' ou 'profil'
%
% sorties :
%	hout = handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 20/07/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function hout=zuiedit_import_mode(nom_mode,mode_fct)

% argument

if nargin < 1
	error('il faut donner le nom de la variable et le mode de fonctionnement') ;
end
if nargin < 2
	error ('il faut donner le mode de fonctionnement ''profil'' ou ''consigne''') ;
end

if strcmp(mode_fct,'profil')
	mode     = 1 ;
	titre = sprintf('profile loading - %s',nom_mode) ;
elseif strcmp(mode_fct,'consigne')
	mode     = 0 ;
	titre = sprintf('Reference loading - %s',nom_mode) ;
else
	error('le mode de fonctionnement ne peut etre que ''profil'' ou ''consigne''')
end

% nom du tag
% ----------
[tagc,reste] = strtok(nom_mode,'.') ;
while ~isempty(reste)
	[tagc,reste] = strtok(reste,'.') ;
end
% on supprime les caracteres : , ( ) du tagc
tagc=strrep(tagc,'(','') ;
tagc=strrep(tagc,':','') ;
tagc=strrep(tagc,',','') ;
tagc=strrep(tagc,')','') ;

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle([tagc '_import_mode']) ;
if ishandle(hform)
        zuiformvisible(hform) ;
	zuiuploadform(hform) ;
	return
end

% formulaire 
form={};

% Titr
sepa ={'separation_comm','frame','',4,''};
form{1} = {sepa};

colj = {'jump','jump','',[],''};
col1 = {'nom_mode','text@full',titre,[],''};
form{length(form)+1} = {colj,col1};

% S�aration
sepa ={'separation_comm','frame','',4,''};
form{length(form)+1} = {sepa};

% source, nom du fichier, format
col1 = {'text_srce','text','source',5,'select data source',''} ;
col2 = {'text_nom' ,'text','filename',30,'file name',''} ;
col3 = {'text_type','text','filetype',5,'file type',''} ;
form{length(form)+1} = {col1,colj,col2,colj,col3,colj,colj} ;

type_srce = {' File      ', ...
             ' Workspace ', ...
             ' Database  '} ;
valeur_srce = {1,0,-1} ;
col1 = {'pop_srce' ,'popup',type_srce,1,'select data source',valeur_srce,''} ;
col2 = {'edit_nom' ,'edit',' ',30,'file name',''} ;
type_file = {'Matlab', ...
             'ASCII '} ;
valeur_file = {0,1} ;
col3 = {'pop_file' ,'popup',type_file,1,'file type',valeur_file,''} ;
col4 = {'choix'    ,'radio','Choose',0,'open file chooser',[],''} ;
form{length(form)+1} = {col1,colj,col2,colj,col3,colj,col4} ;

% Separation
sepa ={'separation_comm','frame','',4,''};
form{length(form)+1} = {sepa};

% Nom des variables
col1 = {'sources','text@full','variable names',1,'',''};
form{length(form)+1} = {col1} ;

% temps, espace, donnee
if mode==1
	texte = 'space' ;
	toolt = 'name of the space variable associated to the data' ;
else
	texte = 'channel' ;
	toolt = 'selected channel in the loaded data : interger >0' ;
end
col1 = {'text_nom_tps'   ,'text','time',20,'name of the time variable',''} ;
col2 = {'text_nom_espace','text',texte,20,toolt,''} ;
col3 = {'text_nom_donnee','text','data',20,'name of the data variable',''} ;
form{length(form)+1} = {col1,colj,col2,colj,col3} ;

col1 = {'edit_nom_tps'   ,'edit',' ',20,'name of the time variable',''} ;
col2 = {'edit_nom_espace','edit',' ',20,toolt,''} ;
col3 = {'edit_nom_donnee','edit',' ',20,'name of the data variable',''} ;
form{length(form)+1} = {col1,colj,col2,colj,col3} ;

% Separation
sepa ={'separation_comm','frame','',4,''};
form{length(form)+1} = {sepa};

% Reechantillonage
col1 = {'sources','text@full','resample',1,'',''};
form{length(form)+1} = {col1} ;

% espace, valeurs, defini
col1 = {'text_espace','text','coordonate',20,'space profile coordinate',''} ;
col2 = {'text_valeur','text','default value',10,'default value out of the interval',''} ;
col3 = {'text_defini','text','always positive',10,'always positive (yes/no)',''} ;
form{length(form)+1} = {col1,colj,col2,colj,col3} ;

type_espace = {'          r/a          ', ...
               '      minor radius      ', ...
               '    sqrt(Phi/pi/b0)    ', ...
               '  sqrt(Phi) normalised '} ;
valeur_espace = {0,1,2,3} ;
%col1 = {'pop_espace' ,'popup',type_espace,3,'space',valeur_espace,''} ;
col1 = {'pop_espace' ,'text','sqrt(Phi/pi/b0)',20,'squareroot of the radial toroidal flux',''} ;
col2 = {'edit_valeur','edit',' ',10,'space / channel',''} ;
type_defini = {'yes', ...
               'no'} ;
valeur_defini = {1,0} ;
col3 = {'pop_defini' ,'popup',type_defini,1,'always positive (yes/no)',valeur_defini,''} ;
form{length(form)+1} = {col1,colj,col2,colj,col3} ;

% Separation
sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

% Formulaire
hout=zuicreeform(' ',[tagc '_import_mode'],'zuiedit_import_mode_fct','',form) ;

[hfig,h] = zuiformhandle([tagc '_import_mode']) ;
setappdata(hfig,'nom_mode',nom_mode)
setappdata(hfig,'mode',mode)


