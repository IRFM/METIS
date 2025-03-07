function zuieastacces

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
     racinegfile = evalin('base','prepare.racine');
catch
     racinegfile ='';
end
if isempty(racinegfile)   
	if strcmp(getenv('USER'),'cgc')
		racinegfile = strcat(getenv('HOME'),'/cgc_data/zineb/data/EAST');
	else
		racinegfile = strcat(getenv('HOME'),'/zineb/data/EAST');
	end
	zassignin('base','prepare.racine',racinegfile);
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
hform=numeastform(chemin,racinegfile);
zuifaitnumeast('init');
zwaitfor(hform)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% intervalle de temps plus scenario
% lecture des donnees du scenario
numchoc       = evalin('base','prepare.numchoc');
nscenario     = evalin('base','prepare.scenario');
racinegfile   = evalin('base','prepare.racine') ;
%
% bug a corrige avec J.F.
%
nscenario  = nscenario(1);
if nscenario > 1
  sc        = zeastscena(nscenario);
else
  racineeq = fullfile(chemin,'EAST',int2str(numchoc),'easteq.mat');
  racinetemp = fullfile(chemin,'EAST',int2str(numchoc),'easttemp.mat');
  if ~exist(racineeq,'file')
    chemin  = '/local/basiuk';
    racineeq = fullfile(chemin,'EAST',int2str(numchoc),'easteq.mat');
    racinetemp = fullfile(chemin,'EAST',int2str(numchoc),'easttemp.mat');
    racinefit = fullfile(chemin,'EAST',int2str(numchoc),'eastfit.mat');
  end    
  load(racinetemp)
  load(racinefit);
  sc.t = easttemp.temps;
  sc.ip = easttemp.ip;
  sc.ne0 = eastfit.nex(:,1)/1e19;
  if length(sc.ne0) == 1
      sc.ne0=sc.ne0*ones(size(sc.t));
  end
  sc.pidn = 0*sc.t;
  sc.pfci = 0*sc.t;
  sc.phyb = 0*sc.t;
  sc.pfce = 0*sc.t;
  sc.zeff = 1*sc.t;
end
zassignin('base','scenario',sc);

% creation du formulaire
tdeb  = min(min(sc.t),0.1);
tfin  = max(sc.t);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
[hout,hui] = scenarioform(tdeb,tfin,numchoc);
hui.axes_plot=zuiplotin(hui.list_plot);
plot(sc.t,sc.ip,'r',sc.t,sc.ne0,'b',sc.t,sc.pidn,'c',sc.t,sc.pfci,'m', ...
     sc.t,sc.phyb,'k');
legend('Ip','ne0','Pnbi','Picrh','Plh')
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
zassignin('base','prepare.nl',sc.ne0');
zassignin('base','prepare.pidn',sc.pidn');
zassignin('base','prepare.pfci',sc.pfci');
zassignin('base','prepare.pfce',sc.phyb');


% dialogue pour la base temps
nom  = sprintf('EAST database edition #%g',numchoc);
aide = 'EAST dataabse construction';
liste_ref = '     Ip    |ne0|Pnbi|icrh|Plh|empty';
var_ref   = {{'prepare.times','prepare.ip',':'}, ...
	             {'prepare.times','prepare.nl',':'}, ...
	             {'prepare.times','prepare.pidn',':'}, ...
	             {'prepare.times','prepare.pfci',':'}, ...
	             {'prepare.times','prepare.phyb',':'}, ...
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

if isempty(scenario)
   scenario(1)=1;
end

% dialogue pour la separatrice
zassignin('base','separatrice',[]);
racineeq = fullfile(chemin,'EAST',int2str(numchoc),'easteq.mat');
racinetemp = fullfile(chemin,'EAST',int2str(numchoc),'easttemp.mat');
if ~exist(racineeq,'file')
  chemin  = '/local/basiuk';
  racineeq = fullfile(chemin,'EAST',int2str(numchoc),'easteq.mat');
  racinetemp = fullfile(chemin,'EAST',int2str(numchoc),'easttemp.mat');
end    
if exist(racineeq,'file')
  load(racineeq)  
  load(racinetemp)
  [teast,iteast]  = sort(easteq.temps);
  separatrice.r0  = (max(easteq.rext(iteast,:)')+min(easteq.rext(iteast,:)'))'/2;
  separatrice.z0  = mean(easteq.zext(iteast,:),2);
  separatrice.a   = (max(easteq.rext(iteast,:)')-min(easteq.rext(iteast,:)'))'/2; 
  separatrice.trh = zeros(size(teast));
  separatrice.trl = zeros(size(teast));
  separatrice.b0  = easttemp.b0;
  separatrice.e1  = (max(easteq.zext(iteast,:)')-min(easteq.zext(iteast,:)'))' ./ ...
                    (max(easteq.rext(iteast,:)')-min(easteq.rext(iteast,:)'))'; 
  separatrice.R   = double(easteq.rext(iteast,:));
  separatrice.Z   = double(easteq.zext(iteast,:));	
  separatrice.temps = teast;
  zassignin('base','separatrice',separatrice);
else
  separatrice = [];
end
while isempty(separatrice)
	info         = zeastsepanew2;
	zassignin('base','sepa_option',info.valeur);
	h=zuicreefunform('zeastsepanew2','sepa_option',1);
	set(h,'name',sprintf('separatrix parameter for  EAST shot #%g',evalin('base','prepare.numchoc')));
	zwaitfor(h)
	code=getappdata(0,'coderetour');
	if ~strcmp(code.action,'validation')
		return
    end
	evalin('base','separatrix = zeastsepanew2(prepare.temps,sepa_option);');
    separatrice = evalin('base','separatrix');
end

% dialogue pour les parametres
info         = zeastacces;
zassignin('base','option',info.valeur);                  
h=zuicreefunform('zeastacces','option',1);
set(h,'name',sprintf('inputfile for EAST shot (Predicitve mode) #%g',evalin('base','prepare.numchoc')));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% appel de la fonction d'acces
hdlg = msgbox('Cronos input file in progress','be patient ...','help');
drawnow;
evalin('base','[cr,data,param]=zeastacces(prepare.numchoc,prepare.chemin,prepare.temps,option,scenario,separatrice);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'Error EASTacces');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('data access problem','EASTacces');
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
function hout=numeastform(chemin,racinegfile)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesEAST');
if ishandle(hfig)
	zuicloseone(hfig) ;
end

% creation du formulaire
form = {};
% ligne pour le numero du choc
col1 = {'text_numchoc','text','shot number :',16,'shot number ',''};
col2 = {'numchoc','edit','',26,'shot number','','prepare.numchoc'};
form{length(form)+1} = {col1,col2};

col1 = {'text_scenario','text','scheme number :',16,'plasma scheme number ',''};
col2 = {'scenario','popup','existing shot| Ohmic Phase (10 s)| Ohmic Phase (60 s) | LHCD Phase | ICRH Phase | ICRH + LHCD | ICRH + LHCD + NB',1,'scheme number',{1,2,3,4,5,6},'prepare.scenario'};
form{length(form)+1} = {col1,col2};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'text_chemin','text@full','Complete path to access Cronos file :',[],'/cgc_data/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full',chemin,42,'Complete path to access Cronos file : ','','prepare.chemin'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'text_racine','text@full','Complete path to access gfile from EAST :',[],'/home/local/basiuk/EAST/',''};
form{length(form)+1} = {col1};
col1 = {'racine','edit@full',racinegfile,42,'Complete path to access gfile from EAST','','prepare.racine'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'error message',''};
form{length(form)+1} = {col1};

hout=zuicreeform('EAST data access','accesEAST','zuifaitnumeast','zuictrlnumeast',form);

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
col2 = {'text_X','text','temps debut',5,'first interval time'};
col3 = {'text_Y','text','temps fin  ',5,'last interval time'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'valeur_xdeb','edit',tdeb,5,'first interval time','','prepare.tdeb'};
col3 = {'valeur_xfin','edit',tfin,5,'last interval time','','prepare.tfin'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

col2 = {'select_xdeb','radio','choose',0,'Choose on the graph the first time'};
col3 = {'select_xfin','radio','choose',0,'Choose on the graph the last time'};
form{length(form)+1} = {col1,col2,col1,col3,col1};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text@full',' ',10,'',[],''} ;
form{length(form)+1} = {col3};

hout=zuicreeform(sprintf('time database fot the EAST shot #%g',numchoc), ...
                 'tempsTS','zuifaittemps','zuicontroletemps',form) ;

[hform,hui] = zuiformhandle('tempsTS');

setappdata(hout,'tdeb',tdeb);
setappdata(hout,'tfin',tfin);

zuifaittemps('init');
