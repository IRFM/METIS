function zuitsacces

% ce progamme est en mode modale
%
% Modifications
% -------------
%  27/09/2001 -> on remplace les wardlg par des msgbox avec un icon 'help
%  10/12/2002 -> interface en anglais
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

[pid,pid_sess,user]=getidprocess;
numchoc = pid;
machine = 'DEMO';
tdeb    = 0;
tfin    = 1;
pas     = 5e-2;
zassignin('base','prepare.numchoc',numchoc);
zassignin('base','prepare.machine',machine);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
zassignin('base','prepare.pas',pas);

hform=videform(chemin,machine,numchoc,tdeb,tfin,pas);
setappdata(hform,'pid',pid);
zuifaitvide('init');

zwaitfor(hform)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

numchoc = evalin('base','prepare.numchoc');
machine = evalin('base','prepare.machine');
tdeb    = evalin('base','prepare.tdeb');
tfin    = evalin('base','prepare.tfin');
pas     = evalin('base','prepare.pas');

% creation de la base temps homogene
temps = [tdeb;tfin];
zassignin('base','prepare.vtemps',temps);
delta = [pas;pas];
zassignin('base','prepare.delta',delta);

% dialogue pour la base temps
nom  = sprintf('time slice for  %s #%g',machine,numchoc);
aide = 'Edition of the saving time';
liste_ref = 'empty';
var_ref   = {{'[]','[]',''}};
hout=zuieditcons(nom,aide,temps,delta,'time (s)','dt (s)', ...
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

% dialogue pour la geometrie
info         = zvidegeo;
option       = info.valeur; 
zassignin('base','geovide',option);
h=zuicreefunform('zvidegeo','geovide',1);
set(h,'name',sprintf('shot geometry %s #%g',machine,numchoc));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% dialogure pour le scenatio
% edition de la consigne ip
nom     = 'Ip';
aide    = 'Plasma current';
x       = temps';
y       = 1e6 .* ones(size(x));
zassignin('base','scenario.ip',y);
texte_x = 'time (s)';
texte_y = 'Ip (A)'; 
var_x   = 'void';
var_y   = 'scenario.ip';
canal   = 1;
h=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'');
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% edition de la consigne Nbar/Ngr
nom     = 'Fn';
aide    = 'density value (compared to the Greenwald limit)';
x       = temps';
y       = 0.8 .* ones(size(x));
zassignin('base','scenario.nbngr',y);
texte_x = 'time (s)';
texte_y = 'Nbar/Ngr (su)'; 
var_x   = 'void';
var_y   = 'scenario.nbngr';
canal   = 1;
h=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'');
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% edition de la consigne zeffm
nom     = 'Zeff';
aide    = '<Zeff> value';
x       = temps';
y       = ones(size(x));
zassignin('base','scenario.zeffm',y);
texte_x = 'time (s)';
texte_y = '<Zeff>'; 
var_x   = 'void';
var_y   = 'scenario.zeffm';
canal   = 1;
h=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'');
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% edition de la consigne ploss
nom     = 'Ploss';
aide    = 'Ploss value';
x       = temps';
y       = 1e6 .* ones(size(x));
zassignin('base','scenario.ploss',y);
texte_x = 'time (s)';
texte_y = 'Ploss (W)'; 
var_x   = 'void';
var_y   = 'scenario.ploss';
canal   = 1;
h=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,'');
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% dialogue pour les parametres
info         = zvideacces;
option       = info.valeur; 
zassignin('base','option',option);
h=zuicreefunform('zvideacces','option',1);
set(h,'name',sprintf('Preparation of the shot %s #%g',machine,numchoc));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end



% appel de la fonction d'acces
hdlg = msgbox('Input data of CRONOS in progress','be patient ...','help');
drawnow;

evalin('base','[cr,data,param]=zvideacces(prepare.machine,prepare.numchoc,prepare.chemin,prepare.temps,option,scenario,geovide);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'Errorr JETacces');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('problem during reading input JET data','JETacces');
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

% appel des autres fenetres
% gaz
hout = zuiedit_param_comp;

% geomertie
hout = zuiedit_data_geom;

% consigne chauffages
hout=zuiedit_data_cons_chauf;

% condition aux limites
hout=zuiedit_data_cons_limites;

% injection de gaz
hout=zuiedit_data_cons_inj;

% profils initiaux
hout=zuiedit_data_prof;

% equilibre initial
zinitequi ;

% formulaire d'edition du numchoc  ....
function hout=videform(chemin,machine,numchoc,tdeb,tfin,pas);



% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesvide');
if ishandle(hfig)
	zuicloseone(hfig) ;
end

% creation du formulaire
form = {}; 
% ligne pour le numero du choc
col1 = {'text_numchoc','text','tokamak name',[],'tokamak name',''};
col2 = {'text_machine','text','shot number ',[],'shot number ',''};
form{length(form)+1} = {col1,col2};
col1 = {'machine','edit',machine,[],'tokamak name','','prepare.machine'};
col2 = {'numchoc','edit',num2str(numchoc),[],'shot number ','','prepare.numchoc'};
form{length(form)+1} = {col1,col2};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'text_tdeb','text','starting time (s)',16,'starting time (s)',''};
col2 = {'text_pas','text','dt (s)',16,'initial dt (s)',''};
col3 = {'text_tfin','text','end time (s)',16,'end time (s)',''};
form{length(form)+1} = {col1,col2,col3};
col1 = {'tdeb','edit',num2str(tdeb),[],'starting time (s)','','prepare.tdeb'};
col2 = {'pas','edit',num2str(tfin),[],'initial dt (s)','','prepare.pas'};
col3 = {'tfin','edit',num2str(pas),[],'end time (s)','','prepare.tfin'};
form{length(form)+1} = {col1,col2,col3};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'text_chemin','text@full','Path to CRONOS input file :',[],'/usr/drfc/usr/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full','chemin',42,'Path to CRONOS input file :','','prepare.chemin'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'error text',''};
form{length(form)+1} = {col1};

hout=zuicreeform('Preparation of an empty input file','accesvide','zuifaitvide','zuictrlvide',form);

