%  ZUITSACCES  formulaire d'acces aux donnees d'un choc TS
%------------------------------------------------------------------------------
% fichier : zuitsacces.m
% 
% fonction Matlab 5 : 
%	creation du formulaire d'acces aux donnees d'un choc TS
%  
% syntaxe :  
%	zuitsacces
%  
% entrees :  
%  
% sorties :  
%  
% fonction write par J-F Artaud , poste 46-78
% version  3.0  du  14/10/2005
%  
% liste des modifications :  
%	* 27/09/2001 -> on remplace les wardlg par des msgbox avec un icon 'help
%  
%------------------------------------------------------------------------------ 
%  

function zuitsacces

% ce progamme est en mode modale
%
% dialogue pour le numero du choc et du chemin 
% numchoc,chemin
try
     chemin = evalin('base','prepare.chemin');
catch
    chemin ='';
end
if isempty(chemin)   
	if strcmp(getenv('USER'),'cgc')
		chemin = strcat(getenv('HOME'),'/cgc_data/zineb/data');
	else
		chemin = strcat(getenv('HOME'),'/zineb/data');
	end
	zassignin('base','prepare.chemin',chemin);
end
try
    numchoc = evalin('base','prepare.numchoc');
catch
    numchoc =[];
end
if isempty(numchoc)
	zassignin('base','prepare.numchoc',tsdernier_choc + 0.1);
end
hform=numform(chemin);
zuifaitnum('init');
zwaitfor(hform)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% intervalle de temps plus scenario
% lecture des donnees du scenario
numchoc = evalin('base','prepare.numchoc')
[ip,times,void1,void2] = tsbase(numchoc,'sprofip');
[nl,times,void1,void2] = tsbase(numchoc,'sprofnl');
[plh,times,void1,void2] = tsbase(numchoc,'sprofplh');
[pfci,times,void1,void2] = tsbase(numchoc,'gprofpfci');
[pecrh,polar,phi_pol,phi_tor,tecrh] = zecrh(fix(numchoc),times);
pecrh(tecrh <0) = [];
tecrh(tecrh <0) = [];
[tnbi,pnbi_cronos] = read_nbi_power(fix(numchoc));
indnbi = (tnbi >=  min(times)) &  (tnbi <=  max(times));
tnbi = tnbi(indnbi);
pnbi_cronos = pnbi_cronos(indnbi);
% creation du formulaire
tdeb = min(times);
tfin = max(times);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
[hout,hui] = scenarioform(tdeb,tfin,numchoc);
hui.axes_plot=zuiplotin(hui.list_plot);
plot(times,ip.*10,'r',times,nl,'b',times,sum(plh,2),'c',times,sum(pfci,2),'m',tecrh,sum(pecrh,2),'g',tnbi,pnbi_cronos/1e6,'k');
legend('Ip','Nl','Plh','Picrh','Pecrh','Pnbi')
title('Scheme');
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
zassignin('base','prepare.nl',nl);
zassignin('base','prepare.plh',sum(plh,2));
zassignin('base','prepare.pfci',sum(pfci,2));

% dialogue pour la base temps
nom  = sprintf('time database for the TS shot #%g',numchoc);
aide = 'time step edition';
liste_ref = '     Ip    |Nl|Plh|Picrh|empty';
var_ref   = {{'prepare.times','prepare.ip',':'}, ...
	             {'prepare.times','prepare.nl',':'}, ...
	             {'prepare.times','prepare.plh',':'}, ...
	             {'prepare.times','prepare.pfci',':'}, ...
	             {'[]','[]',''}};
hout=zuieditcons(nom,aide,temps,delta,'time (s)','time step (s)', ...
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
info         = ztsacces;
zassignin('base','option',info.valeur);                  
h=zuicreefunform('ztsacces','option',1);
set(h,'name',sprintf('shot TS %g',evalin('base','prepare.numchoc')));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% appel de la fonction d'acces
hdlg = msgbox('Creation of the Cronos structure','be patient ...','help');
drawnow;
evalin('base','[cr,data,param]=ztsacces(prepare.numchoc,prepare.chemin,prepare.temps,option);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'Erreur TSacces');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('data access problem','TSacces');
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
		tps=[num2str(tps1,2) ' -> ' num2str(tps2,2)];
		zuidata(h.text_temps,tps);
	end
end


% formulaire d'edition du numchoc et du chemin
function hout=numform(chemin)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesTS');
if ishandle(hfig)
	zuicloseone(hfig) ;
end

% creation du formulaire
form = {}; 
% ligne pour le numero du choc
col1 = {'text_numchoc','text','shot number :',16,'shot number ',''};
col2 = {'numchoc','edit','',26,'shot number','','prepare.numchoc'};
form{length(form)+1} = {col1,col2};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'text_chemin','text@full','complete path access :',[],'/usr/drfc/usr/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full',chemin,42,'complete path access ','','prepare.chemin'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'messages d''erreur',''};
form{length(form)+1} = {col1};

hout=zuicreeform('Tore Supra data access','accesTS','zuifaitnum','zuicontrolnum',form);

function [hout,hui] = scenarioform(tdeb,tfin,numchoc)

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
col2 = {'text_X','text','first time',5,'first interval time'};
col3 = {'text_Y','text','last time  ',5,'last inetrval time'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'valeur_xdeb','edit',tdeb,5,'first interval time','','prepare.tdeb'};
col3 = {'valeur_xfin','edit',tfin,5,'last interval time','','prepare.tfin'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'select_xdeb','radio','select',0,'choose on the graphic the first interval time'};
col3 = {'select_xfin','radio','select',0,'choose on the graphic the last interval time'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(sprintf('TS time interval %g',numchoc), ...
                 'tempsTS','zuifaittemps','zuicontroletemps',form) ;

[hform,hui] = zuiformhandle('tempsTS');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);

zuifaittemps('init');
