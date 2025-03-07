%  ZUIACCES  formulaire d'acces aux donnees d'un fichier cronos
%------------------------------------------------------------------------------
% fichier : zuiacces.m
% 
% fonction Matlab 5 : 
%	creation du formulaire d'acces aux fichiers cronos
%  
% syntaxe :  
%	zuiacces
%  
% entrees :  
%  
% sorties :  
%  
% fonction �rite par J-F Artaud , poste 46-78
% version  1.7  du  29/09/2001  
%  
% liste des modifications :  
%	* 27/09/2001 -> on remplace les wardlg par des msgbox avec un icon 'help
%	* 16/06/2003 -> ajout de la lecture des donnees de tfce
%  
%------------------------------------------------------------------------------ 
%  

function [data,param,post] = zuiacces

[data,param,post] = zuiload;

file  = param.edit.currentfile;

% intervalle de temps plus scenario
% lecture des donnees du scenario

times = data.gene.temps;
ip    = data.gene.ip;
nbar  = data.gene.nbar;
plh   = data.gene.paddhyb;
pfci  = data.gene.paddfci;
pfce  = data.gene.paddfce;

% creation du formulaire
tdeb = min(times);
tfin = max(times);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
[hout,hui] = scenarioform(tdeb,tfin,file);

hui.axes_plot=zuiplotin(hui.list_plot);
plot(times,ip.*10,'r',times,nbar,'b',times,plh,'c',times,pfci,'m',times,pfce,'g');
legend('Ip','Nbar','Plh','Picrh','Pecrh')
title('Scenario');
xlabel('time (s)');
ylabel('');

% attente retour
zwaitfor(hout)

code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end


% creation de la base temps homogene
pas= mean(diff(times));
tdeb = evalin('base','prepare.tdeb');
tfin = evalin('base','prepare.tfin');
temps = [tdeb;tfin];
zassignin('base','prepare.vtemps',temps);
delta = [pas;pas];
zassignin('base','prepare.delta',delta);
zassignin('base','prepare.times',times);
zassignin('base','prepare.ip',ip.*10);
zassignin('base','prepare.nbar',nbar);
zassignin('base','prepare.plh',plh);
zassignin('base','prepare.pfci',pfci);

% dialogue pour la base temps
nom  = ['Edition de la basetemps pour le fichier ' file];
aide = 'Edition du pas de temps en fonctiuon du temps';
liste_ref = '     Ip    |Nbar|Plh|Picrh|empty';
var_ref   = {{'prepare.times','prepare.ip',':'}, ...
	             {'prepare.times','prepare.nbar',':'}, ...
	             {'prepare.times','prepare.plh',':'}, ...
	             {'prepare.times','prepare.pfci',':'}, ...
	             {'[]','[]',''}};
hout=zuieditcons(nom,aide,temps,delta,'temps (s)','pas de temps (s)', ...
                 'prepare.vtemps','prepare.delta',1,'',liste_ref,var_ref ...
                 ,'','','HandleVisibility','on');
% attente retour

zwaitfor(hout)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end


% calcul de la base temps
temps =[];
vtemps = evalin('base','prepare.vtemps');
delta = evalin('base','prepare.delta');
for k =1:(length(vtemps)-1)
	temps = cat(2,temps,vtemps(k):delta(k):vtemps(k+1));
end
ind = find(diff(temps) <= (min(delta)/10));
if ~isempty(ind)
	temps(ind)=[];
end
zassignin('base','prepare.temps',temps(:));

data = zselt(data,iround(times,temps));


function [hout,hui] = scenarioform(tdeb,tfin,file)

form={};
% resevation de la place pour le plot
col1 = {'list_plot','list','*********************************************************************************************************************',1,''};
form{length(form)+1} = {col1};
col1 = {'jump_void','jump','***********************************************************',[],''};
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
col2 = {'text_X','text','initial time',5,'initial interval time'};
col3 = {'text_Y','text','final time  ',5,'final interval time'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'valeur_xdeb','edit',tdeb,5,'initial interval time','','prepare.tdeb'};
col3 = {'valeur_xfin','edit',tfin,5,'final interval time','','prepare.tfin'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'select_xdeb','radio','select',0,'choose the initial time directly on the graphic'};
col3 = {'select_xfin','radio','select',0,'choose the final time directly on the graphic'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(['time interval of the CRONOS file ' file], ...
                 'tempsTS','zuifaittemps','zuicontroletemps',form, ...
                 {},0,0,'HandleVisibility','on') ;

[hout,hui] = zuiformhandle('tempsTS');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);

zuifaittemps('init');
