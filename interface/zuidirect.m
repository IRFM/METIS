% ZUIDIRECT Creation du formulaire principal de Zineb en mode direct
%--------------------------------------------------------------
% fichier zuidirect.m ->  
%		    zuicreeform : creation du formulaire
% 		    zuidirect_fonction : fonction de gestion des callbacks des uicontrols
%					                du formulaire
%
% fonction Matlab 5 :
%	Cette fonction definie la feuille principale (interface graphique) sous 
%	la forme d'un formulaire. Le formulaire est constitue de lignes.
%	Chaque ligne peut avoir un nombre de colonnes quelconques. Les lignes
%	qui se suivent et qui ont un nombre de colonnes egal, voient leurs colonnes
%	alignees verticalement.
%
% syntaxe  :
%  
% entrees
%
% sorties 
%
% fonction ecrite par C. Passeron, poste 61 19
% version 3.0, du 13/10/2005.
%
% liste des modifications : 
%--------------------------------------------------------------

function zuidirect

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('direct');
if ishandle(hform)
        zuiformvisible(hform);
	return
end
% 1eres lignes : Nom du fichier de travail
% --------------------------------------
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full',' Filename',52,''};
form{length(form)+1} = {col1};

filename = evalin('base','param.edit.currentfile','[''  '']') ;
col1 = {'text_loadfile','text@full',filename,52,'Name of the presently edited file',[],'param.edit.currentfile'};
form{length(form)+1} = {col1} ;

% 2emes lignes : Machine, Numero de choc, occurence, temps
% ------------------------------------------------------
col1 = {'libel_nom_machine' ,'text'  ,'Tokamak',13,''} ;
col2 = {'libel_numchoc'     ,'text'  ,'Shot number',13,''} ;
col3 = {'libel_occurence'   ,'text'  ,'Occurence',13,''} ;
col4 = {'libel_temps'       ,'text'  ,'Time (s)',13,''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

machine = evalin('base','param.from.machine','[''  '']') ;
numshot = evalin('base','param.from.shot.num','[]') ;
tps1 = evalin('base','param.gene.tdeb','[]') ;
tps2 = evalin('base','param.gene.tfin','[]') ;
tps = [num2str(tps1,2) ' -> ' num2str(tps2,2)] ;

numchoc   = fix(numshot) ;
occurence = round((numshot-numchoc)*10) ;

col1 = {'text_nom_machine','text',machine,[],'Tokamak name',[],'param.from.machine'};
col2 = {'text_numchoc'    ,'text',numchoc,[],'shot number',[],[]};
col3 = {'text_occurence'  ,'text',occurence,[],'number of version in database',[],[]};
col4 = {'text_temps'      ,'text',tps,[],'time',[],[]};
form{length(form)+1} = {col1,col2,col3,col4};
%celldisp(form)

% separation
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% 3eme ligne : Chargement, Rename, Edition du fichier
% ---------------------------------------------------
col1 = {'radio_loadfile','radio@left' ,'Load',0,'Load a CRONOS file'};
col2 = {'radio_renamefile','radio@center','Save as',0,'Save CRONOS workspace in another file'};
col3 = {'radio_savefile','radio@right','Save',0,'Save workspace with same filename'};
form{length(form)+1} = {col1,col2,col3};

% separation
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% 4eme ligne : Cr�tion, Edition 
% ------------------------------
col1 = {'radio_createfile','radio@left',' New file ',0,'Create a new CRONOS file'};
col3 = {'radio_editfile','radio@right','   Edit ',0,'Edit data and simulation parameters'};
form{length(form)+1} = {col1,col3};

% separation
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% 5eme ligne : Execution Interactif ou Batch
% ------------------------------------------
col1 = {'radio_runinter','radio@left' ,'Interactive run',0,'Run CRONOS in the current workspace'};
col2 = {'radio_runbatch','radio@right','Batch run',0,'Run CRONOS in batch'};
form{length(form)+1} = {col1,col2};

% 6eme ligne : Rebuilt - Visualisation
% -------------------------------------------------
col1 = {'radio_rebuilt','radio@left' ,'Recover results',0,'Recover results from interrupted simulation'};
col2 = {'radio_visu'   ,'radio@right','Visualisation',0,'Display data'};
form{length(form)+1} = {col1,col2};

% 7eme ligne : Visualisation - Affichage du resume
% -------------------------------------------------
col1 = {'radio_resume'   ,'radio@left' ,'Summary',0,'Simulation summary'};
col2 = {'radio_assistant','radio@right','Assistant',0,'Assistant: data manipulation'};
form{length(form)+1} = {col1,col2};

% 8eme ligne : Mise �jour du code dynamique - Mise ajour du fichier de travail
% ------------------------------------------------------------------------------
col1 = {'radio_majcode','radio@left' ,'Update code',0,'Update code'};
col2 = {'radio_majfile','radio@right','Update data',0,'Update data structure'};
form{length(form)+1} = {col1,col2};

% 9eme ligne : Post-Traitement
% ----------------------------
col1 = {'radio_path'     ,'radio@left' ,'Change path',0,'Change path'};
col2 = {'radio_posttrait','radio@right','Post-Processing',0,'Post-Processing'};
form{length(form)+1} = {col1,col2};

% separation
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% 10eme ligne : Affichage du logfile
% ---------------------------------
col1 = {'radio_tracelog','radio@left' ,'Log file',0,'Display Log file'};
col2 = {'radio_suivilog','radio@center','monitor execution',0,'Display CRONOS messages during batch run'};
col3 = {'radio_loadresult','radio@right','Load result',0,'Load the related result file'};
form{length(form)+1} = {col1,col2,col3};

% separation
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% 11eme ligne : Bouton Quitter
% ----------------------------
comm{1}={'btn_quit','radio','Quit',0,''};
comm{2}={'aide','radio@right','Help',0,''};

hout=zuicreeform(sprintf('CRONOS (id %d)',getidprocess),'direct','zuidirect_fct','',form,comm,0,0,'resize','on','visible','on');

