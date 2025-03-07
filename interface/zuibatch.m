% ZUIBATCH formulaire de soumission en batch de zineb
%--------------------------------------------------------------
% fichier zuibatch.m ->  zuicreeform
%
% fonction Matlab 5 :
% 
% Cette fonction cree la feuille (interface graphique) sous 
% la forme d'un formulaire. Le formulaire est constitue de lignes.
% Chaque ligne peut avoir un nombre de colonnes quelconques. Les lignes
% qui se suivent et qui ont un nombre de colonnes egal, voient leurs colonnes
% alignees verticalement.
%
% syntaxe  :
%	zuibatch  
%
% entrees
%
% fonction �rite par C. Passeron , poste 61-19
% version  2.1  du 04/07/2003 
%
% sorties
% 
% liste des modifications : 
% * 28/11/2001 -> modification de la presentation
% * 18/03/2002 -> ajout du mode reprise sans reinitilisation de l'equilibre et des sources
% * 11/04/2002 -> remplacement hercule2 en hercule
% * 17/09/2002 -> ajout du pc cronos
% * 26/11/2002 -> ajout des jac du jet
% * 26/11/2002 -> nouveau batch pour TS
%
%--------------------------------------------------------------

function zuibatch

% On recupere sur une machine ouverte au batch LSF, la liste des machines
%[cr,bhosts]  = unix('rsh deneb bhosts  | cut -d'' '' -f1 | egrep -v ''HOST_NAME''  | tr -s ''\012'' ''|'' ');
bhosts = ['Batch TS|local'];
thosts ={'alcyone',''};
hostname = getenv('HOSTNAME');
if isempty(hostname)
     [s,hostname] = unix('/bin/hostname');
     if s~=0
 	    [s,hostname] = unix('/usr/bin/hostname');    
     end
end

if strmatch('pc-cronos',hostname)
   bhosts = 'pc-cronos';
   thosts ={'pc-cronos'};
elseif findstr(hostname,'jac-')
   bhosts = 'auto|locale';
   thosts ={'jac',getenv('HOSTNAME')};
elseif ~isempty(findstr(getenv('HOSTNAME'),'zeus')) | ~isempty(findstr(getenv('HOSTNAME'),'eos'))
   bhosts = 'zeus|eos';
   thosts ={'zeus','eos'};
elseif findstr(hostname,'saturne')
	bhosts = 'auto';
        thosts ={'saturne'};
end
bqueues = {'All'};

% On recupere sur une machine ouverte au batch LSF, la liste des queues batch
%[cr,bqueues] = unix('rsh deneb bqueues | cut -d'' '' -f1 | egrep -v ''QUEUE_NAME'' | tr -s ''\012'' ''|'' ');

filename = evalin('base','param.edit.currentfile','[]');
if isempty(filename)
	disp(' no Cronos input data')
	warndlg(' no Cronos input data','Zineb run batch');
	return
end

% formulaire
% ----------
form={};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

% 1ere ligne : Nom du fichier de la simulation
% --------------------------------------------
col1 = {'libel','text@full',' Cronos input data name',[],''};
form{length(form)+1} = {col1};

col1 = {'text','text@full',filename,35,'Cronos input data name',[],[]};
form{length(form)+1} = {col1};

% separation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% 2eme ligne : Nom des machines ouvertes au batch LSF
% ----------------------------------------------------
colj=  {'jump1','jump','jump',[],''};
%col1 = {'text_nom_machine_batch' ,'text'  ,' CPU name for batch ',[],''};
%col2 = {'popup_nom_machine_batch','popup' ,bhosts,1,' CPU available for batch ',thosts,'void'};
%form{length(form)+1} = {col1,colj,col2};

% 3eme ligne : liste des queues ouvertes au batch LSF
% ----------------------------------------------------
col1 = {'text_nom_queue_batch' ,'text'  ,'  batch queue name ',[],''};
col2 = {'popup_nom_queue_batch','popup' ,'Auto',1,'the choice of the queue is automatic',{'long','short','test'},'void'};
form{length(form)+1} = {col1,colj,col2};

% 3eme ligne : reprise
% --------------------
filerep  = evalin('base','param.gene.file','[]');
%
colj=  {'jump1','jump','jump',[],''};
colm1 = {'text_reprise_run','text','Cronos run',[],''};
if ~isempty(filerep)
	col1 = {'radio_reprise_batch','popup','Complete          |Restart + Init| Restart without init|Post-processing' ,1,'renstart of a partial Cronos run with/without initialization of the equilibrium and source term',{0,1,2,3}};

else
	col1 = {'radio_reprise_batch','jump',' ' ,0,''};
end
form{length(form)+1} = {colm1,colj,col1};

% separation
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

hout=zuicreeform('Interface batch Zineb','batch','zuibatch_fonction','zuibatch_control',form);
