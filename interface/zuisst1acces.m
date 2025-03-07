function zuisst1acces

% ce progamme est en mode modale
%
% Modifications
% -------------
%  27/09/2001 -> on remplace les wardlg par des msgbox avec un icon 'help
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
	zassignin('base','prepare.numchoc',1);
end
try
    scenario = evalin('base','prepare.scenario');
catch
    scenario =[];
end
if isempty(scenario)
	zassignin('base','prepare.scenario',1);
end
hform=numsst1form(chemin);
zuifaitnumsst1('init');
zwaitfor(hform)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% intervalle de temps plus scenario
% lecture des donnees du scenario
numchoc   = evalin('base','prepare.numchoc');
nscenario = evalin('base','prepare.scenario');
sc        = zsst1scena(nscenario);
zassignin('base','scenario',sc);

% creation du formulaire
tdeb  = max(min(sc.t),0.1);
tfin  = max(sc.t);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
[hout,hui] = scenarioform(tdeb,tfin,numchoc);
hui.axes_plot=zuiplotin(hui.list_plot);
plot(sc.t,sc.ip,'r',sc.t,sc.nbar,'b',sc.t,sc.pidn,'c',sc.t,sc.pfci,'m', ...
     sc.t,sc.pfce,'k',sc.t,sc.plh,'g--');
legend('Ip','Nbar','Pnbi','Picrh','Pecrh','Plh')
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
pas= 1;
tdeb = evalin('base','prepare.tdeb');
tfin = evalin('base','prepare.tfin');
temps = [tdeb;tfin];
zassignin('base','prepare.vtemps',temps);
delta = [pas;pas];
zassignin('base','prepare.delta',delta);
zassignin('base','prepare.times',sc.t');
zassignin('base','prepare.ip',sc.ip');
zassignin('base','prepare.nl',sc.nbar');
zassignin('base','prepare.pidn',sc.pidn');
zassignin('base','prepare.pfci',sc.pfci');
zassignin('base','prepare.pfce',sc.pfce');
zassignin('base','prepare.plh',sc.plh');
zassignin('base','prepare.ntnd',sc.ftri');

% dialogue pour la base temps
nom  = sprintf('SST1 time database, shot #%g',numchoc);
aide = 'time step as a function of time';
liste_ref = '     Ip    |Nbar|Pnbi|Picrh|Pecrh|PLH|nT/nD|empty';
var_ref   = {{'prepare.times','prepare.ip',':'}, ...
	             {'prepare.times','prepare.nl',':'}, ...
	             {'prepare.times','prepare.pidn',':'}, ...
	             {'prepare.times','prepare.pfci',':'}, ...
	             {'prepare.times','prepare.pfce',':'}, ...
	             {'prepare.times','prepare.plh',':'}, ...
	             {'prepare.times','prepare.ntnd',':'}, ...
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


% dialogue pour la separtrice
zassignin('base','separatrice',[]);
separatrice = [];
while (isempty(separatrice))
	info         = zsst1sepa;
	zassignin('base','sepa_option',info.valeur);                  
	h=zuicreefunform('zsst1sepa','sepa_option',1);
	set(h,'name',sprintf('SST1 separatrix parameters, #%g',evalin('base','prepare.numchoc')));
	zwaitfor(h)
	code=getappdata(0,'coderetour');
	if ~strcmp(code.action,'validation')
		return
	end
	evalin('base','separatrice = zsst1sepa(prepare.temps,sepa_option);');
	separatrice = evalin('base','separatrice');
end
% dialogue pour les parametres
info         = zsst1acces;
zassignin('base','option',info.valeur);                  
h=zuicreefunform('zsst1acces','option',1);
set(h,'name',sprintf('SST1 shot preparation,  #%g',evalin('base','prepare.numchoc')));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% appel de la fonction d'acces
hdlg = msgbox('data structure in progress','be patient ...','help');
drawnow;
evalin('base','[cr,data,param]=zsst1acces(prepare.numchoc,prepare.chemin,prepare.temps,option,scenario,separatrice);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'Erreur SST1acces');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('data access problem','SST1acces');
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
		tps=[num2str(tps1) ' �' num2str(tps2)];
		zuidata(h.text_temps,tps);
	end
end


% formulaire d'edition du numchoc et du chemin
function hout=numsst1form(chemin)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesSST1');
if ishandle(hfig)
	zuicloseone(hfig) ;
end

% creation du formulaire
form = {}; 
% ligne pour le numero du choc
col1 = {'text_numchoc','text','shot number :',16,'shot number ',''};
col2 = {'numchoc','edit','',26,'shot number','','prepare.numchoc'};
form{length(form)+1} = {col1,col2};

col1 = {'text_scenario','text','scheme number :',16,'scheme number ',''};
col2 = {'scenario','popup','1',1,'scheme number',{1},'prepare.scenario'};
form{length(form)+1} = {col1,col2};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'text_chemin','text@full','path access to the input files :',[],'/usr/drfc/usr/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full',chemin,42,'path access to the input files ','','prepare.chemin'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'messages d''erreur',''};
form{length(form)+1} = {col1};

hout=zuicreeform('SST1 shot preparation','accesSST1','zuifaitnumsst1','zuictrlnumsst1',form);

function [hout,hui] = scenarioform(tdeb,tfin,numchoc)

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

col1 = {'jump_void','jump','place place',[],''};
col2 = {'text_X','text','first time',5,'first time of the simulation'};
col3 = {'text_Y','text','last time  ',5,'last time of the simulation'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'valeur_xdeb','edit',tdeb,5,'first time of the simulation','','prepare.tdeb'};
col3 = {'valeur_xfin','edit',tfin,5,'last time of the simulation','','prepare.tfin'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'select_xdeb','radio','selectionner',0,'choose on the graphic the first interval time'};
col3 = {'select_xfin','radio','selectionner',0,'choose on the graphic the last interval time'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(sprintf('SST1 time interval #%g',numchoc), ...
                 'tempsTS','zuifaittemps','zuicontroletemps',form) ;

[hform,hui] = zuiformhandle('tempsTS');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);

zuifaittemps('init');
