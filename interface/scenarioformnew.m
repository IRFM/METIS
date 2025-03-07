function [hout,hui] = scenarioformnew(tdeb,tfin,numchoc,tmse,tpol)

form={};
% resevation de la place pour le plot
col1 = {'list_plot','list','Ceci est un texte tres long pour reserver de la place & Ceci est un texte tres long pour reserver de la place',1,''};
form{length(form)+1} = {col1};
col1 = {'jump_void','jump','Ceci est un texte tres long pour reserver de la place',[],''};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

col1 = {'jump_void','jump','place',[],''};
col2 = {'text_X','text','initial time',5,'first interval time'};
col3 = {'text_Y','text','final time  ',5,'last interval time'};
col4 = {'text_liste','text','analysis time  ',5,'first time'};
col5 = {'text_liste','text','time index  ',5,' MSE or Pol index'};
form{length(form)+1} = {col2,col1,col3,col1,col4,col5};

col2 = {'valeur_xdeb','edit',tdeb,5,'first interval time ','','prepare.tdeb'};
col3 = {'valeur_xfin','edit',tfin,5,'last interval time','','prepare.tfin'};

if ~isempty(tmse) & ~isempty(tpol)
  sliste =sprintf('Standard|Efit|MSE|Pol')
  vliste = {1,2,3,4};
elseif ~isempty(tmse)
  sliste =sprintf('Standard|Efit|MSE')
  vliste = {1,2,3};
elseif ~isempty(tpol)
  sliste =sprintf('Standard|Efit|Pol')
  vliste = {1,2,4};
else
  sliste =sprintf('Standard|Efit')
  vliste = {1,2};
end
col4   = {'select_listediag','popup',sliste,2,'first current profile origin (EFIT, MSE, ...)',vliste,'prepare.tmode'};
col5   = {'select_listetemps','popup','',0,'','','prepare.tmodev'};
    
form{length(form)+1} = {col2,col1,col3,col1,col4,col5};

col2 = {'select_xdeb','radio','Select',0,'Select on the graph the first time'};
col3 = {'select_xfin','radio','Select',0,'Sélect on the graph the last time'};
form{length(form)+1} = {col2,col1,col3,col1,col1,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(sprintf('time database for the shot %g',numchoc), ...
                 'tempsJET','zuifaitjetnew','zuicontroltjetnew',form) ;

[hform,hui] = zuiformhandle('tempsJET');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);
setappdata(hout,'tmse',tmse);
setappdata(hout,'tpol',tpol);

zuifaitjetnew('init');
