function tana = plottuneluke(temps)

tdeb          = min(temps);
tfin          = max(temps);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
numchoc       = 0;
[hout,hui]    = scenariotune(tdeb,tfin,numchoc);
hui.axes_plot = zuiplotin(hui.list_plot);
data          = evalin('base','data');
conslh        = sum(abs(data.cons.hyb),2);
plot(data.gene.temps,data.gene.ip,data.gene.temps,data.gene.nbar/1e20,data.gene.temps,conslh/1e6)
legend('Ip(MA)','nbar','Plh')
title('Scheme');
xlabel('time (s)');
ylabel('');
% attente retour
zwaitfor(hout)
tana = evalin('base','prepare.tdeb');
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end







function [hout,hui] = scenariotune(tdeb,tfin,numchoc)

form={};
% resevation de la place pour le plot
col1 = {'list_plot','list','*************************************************************************************************************************',1,''};
form{length(form)+1} = {col1};
col1 = {'jump_void','jump','*********************************************************',[],''};
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

col1 = {'jump_void','jump','place place',[],''};
col2 = {'text_X','text','calculation time',5,'first interval time'};
form{length(form)+1} = {col1,col2,col1};

col2 = {'valeur_xdeb','edit',tdeb,5,'first interval time','','prepare.tdeb'};
form{length(form)+1} = {col1,col2,col1};

col2 = {'select_xdeb','radio','select',0,'choose on the graphic the first interval time'};
form{length(form)+1} = {col1,col2,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(sprintf('LUKE time interval %g',numchoc), ...
                 'tempsTS','zuifaittempsluke','zuicontroletime',form) ;

[hform,hui] = zuiformhandle('tempsTS');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);

zuifaittempsluke('init');




