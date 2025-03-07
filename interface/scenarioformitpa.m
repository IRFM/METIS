function [hout,hui] = scenarioformitpa(tdeb,tfin,dt,numchoc)

form={};
% reservation de la place pour le plot
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
col2 = {'text_X','text','temps debut ',5,'temps de début de l''intervalle'};
col3 = {'text_Y','text','temps fin   ',5,'temps de fin de l''intervalle'};
col4 = {'text_Z','text','pas de temps',5,'pas de temps de sauvegarde'};
form{length(form)+1} = {col2,col1,col3,col1,col4};

col2 = {'valeur_xdeb','edit',tdeb,5,'temps de début de l''intervalle','','prepare.tdeb'};
col3 = {'valeur_xfin','edit',tfin,5,'temps de fin de l''intervalle','','prepare.tfin'};
col4 = {'valeur_xdt','edit',dt,5,'pas de temps de sauvegarde','','prepare.dt'};

form{length(form)+1} = {col2,col1,col3,col1,col4};

col2 = {'select_xdeb','radio','selectionner',0,'Sélectionner sur le graphique le temps de début de l''intervalle'};
col3 = {'select_xfin','radio','selectionner',0,'Sélectionner sur le graphique le temps de fin de l''intervalle'};
form{length(form)+1} = {col2,col1,col3,col1,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(sprintf('Intervalle de temps de la simulation pour le choc %g',numchoc), ...
                 'tempsITPA','zuifaitscenarioitpa','zuictrlscenarioitpa',form) ;

[hform,hui] = zuiformhandle('tempsITPA');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);
setappdata(hout,'dt',dt);

zuifaitscenarioitpa('init');
