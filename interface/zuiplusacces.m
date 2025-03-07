function zuiplusacces

% ce progamme est en mode modale
%
% Modifications
% -------------
%  27/09/2001 -> on remplace les wardlg par des msgbox avec un icon 'help
%

%  1- on regarde s'il y a un jeux de donnee en memoire
% numchoc,chemin
try
   ok = evalin('base','isstruct(param)&isstruct(data)&~isempty(param.edit.currentfile)');
catch
   ok = 0;
end
if ok == 1
	rep = questdlg('Do you want to use the current file ?', ...
	               'Cronos simulation extension','Yes','No','Yes');
	if strcmp(rep,'No')
		ok = 0;
	end
end
if ~ok
	zuiload
end


% intervalle de temps plus scenario
% lecture des donnees du scenario
numchoc = evalin('base','param.from.shot.num');
machine = evalin('base','param.from.machine');
times   = evalin('base','data.gene.temps');
ip      = evalin('base','data.cons.ip');
if all(ip ==0) |all(~isfinite(ip))
	ip      = evalin('base','data.equi.ip');
end
pfci    = abs(evalin('base','data.cons.fci'));
pfce    = abs(evalin('base','data.cons.fce'));
phyb    = abs(evalin('base','data.cons.hyb'));
pidn    = abs(evalin('base','data.cons.idn'));
nbar    = abs(evalin('base','data.gene.nbar'));

% creation du formulaire
ButtonName=questdlg('minimum time allowed ?', ...
                         'Causality', ...
                         '0','Tmin','-3600','Tmin');
switch ButtonName,
   case '0'
      tdeb = 0;
  case 'Tmin'
      tdeb = min(times);
   case '-3600'
      tdeb = -3600;
end % switch

%tdeb = min(times);
tfin = 24*3600;
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
[hout,hui] = scenarioform(tdeb,tfin,numchoc,machine);
hui.axes_plot=zuiplotin(hui.list_plot);
plot(times,ip./1e5,'r',times,nbar./1e19,'b',times,sum(phyb,2)./1e6,'c',times,sum(pfci,2)./1e6,'m', ...
     times,sum(pidn,2)./1e6,'k',times,sum(pfce,2)./1e6,'g');
legend('Ip','Nbar','Phyb','Pfci','Pidn','Pfce')
title('Scheme');
xlabel('time (s)');
ylabel('');
zoom xon

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
zassignin('base','prepare.ip',ip./1e5);
zassignin('base','prepare.nbar',nbar./1e19);
zassignin('base','prepare.phyb',sum(phyb,2));
zassignin('base','prepare.pfci',sum(pfci,2));
zassignin('base','prepare.pidn',sum(pidn,2));
zassignin('base','prepare.pfce',sum(pfce,2));
zassignin('base','prepare.machine',machine);

% dialogue pour la base temps
nom  = sprintf('time database %s #%g',machine,numchoc);
aide = 'time database evolution';
liste_ref = '     Ip    |Nbar|Pnbi|Plh|Pfci|Pfce|vide';
var_ref   = {{'prepare.times','prepare.ip',':'}, ...
	             {'prepare.times','prepare.nbar',':'}, ...
	             {'prepare.times','prepare.pidn',':'}, ...
	             {'prepare.times','prepare.phyb',':'}, ...
	             {'prepare.times','prepare.pfci',':'}, ...
	             {'prepare.times','prepare.pfce',':'}, ...
	             {'[]','[]',''}};
hout=zuieditcons(nom,aide,temps,delta,'time (s)','time database (s)', ...
                 'prepare.vtemps','prepare.delta',1,'',liste_ref,var_ref);
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

% dialogue pour les parametres
% parametre par defaut
info         = zplusacces;
valeur       = info.valeur;
% recuperation des dernieres valeur
param         = evalin('base','param');
mode         = evalin('base','data.mode');
valeur.psimode     = round(mean(mode.psi));
valeur.pemode      = round(mean(mode.pe));
valeur.pionmode    = round(mean(mode.pion));
valeur.nelmode     = round(mean(mode.nel));
valeur.psilim      = round(mean(mode.cons.psi));
valeur.lambda      = param.gene.lambda;
valeur.modecoef    = param.gene.modecoef;
valeur.self        = param.gene.self;
valeur.fast        = param.gene.fast;
valeur.source_bord = param.gene.source_bord;
valeur.cn          = param.gene.cn;
valeur.plotonoff   = param.plot.onoff;
valeur.verbose     = param.gene.verbose;
valeur.rebuilt     = param.gene.rebuilt;
valeur.post        = param.gene.post;

zassignin('base','option',valeur);


h=zuicreefunform('zplusacces','option',1);
set(h,'name',sprintf('cronos simulation extension %s #%g',machine,numchoc));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% appel de la fonction d'acces
hdlg = msgbox('data structure in progress','Be patient ...','help');
drawnow;
evalin('base','[cr,data,param]=zplusacces(param,data,prepare.temps,option);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'zplusacces error');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('access data problem','zplusacces');
	return
end
% position flag d'edition
zuisavenonok;
evalin('base','param.edit.currentfile=param.gene.origine;');

% mise a jour menu principal de cronos
[hfig,h] = zuiformhandle('direct');
if ishandle(hfig)
	zuiuploadform(hfig) ;

	numshot=evalin('base','param.from.shot.num','[]');
	if ~isempty(numshot)
		numchoc = fix(numshot) ;
		zuidata(h.text_numchoc,num2str(numchoc));

		occurence = round((numshot-numchoc)*10);
		zuidata(h.text_occurence,num2str(occurence));

		tps1=evalin('base','param.gene.tdeb');
		tps2=evalin('base','param.gene.tfin');
		tps=[num2str(tps1,3) ' -> ' num2str(tps2,3)];
		zuidata(h.text_temps,tps);
	end
end

% consigne chauffages
hout=zuiedit_data_cons_chauf
zwaitfor(hout,'visible')

% condition aux limites
hout=zuiedit_data_cons_limites
zwaitfor(hout,'visible')

% injection de gaz
hout=zuiedit_data_cons_inj
zwaitfor(hout,'visible')


function [hout,hui] = scenarioform(tdeb,tfin,numchoc,machine)

form={};
% resevation de la place pour le plot
col1 = {'list_plot','list','**********************************************************************************************************',1,''};
form{length(form)+1} = {col1};
col1 = {'jump_void','jump','**********************************************************',[],''};
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
col2 = {'text_X','text','initial time',5,'initial time of the interval'};
col3 = {'text_Y','text','final time  ',5,'final time of the interval'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'valeur_xdeb','edit',tdeb,5,'initial time of the interval','','prepare.tdeb'};
col3 = {'valeur_xfin','edit',tfin,5,'final time of the interval','','prepare.tfin'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'select_xdeb','radio','select',0,'choose the initial time directly on the graph'};
col3 = {'select_xfin','radio','select',0,'choose the final time directly on the graph'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(sprintf('time interval, shot %s #%g',machine,numchoc), ...
                 'tempsTS','zuifaittemps','zuicontroletemps',form) ;

[hform,hui] = zuiformhandle('tempsTS');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);

zuifaittemps('init');
