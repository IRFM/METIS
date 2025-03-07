%  ZUITCVACCES  formulaire d'acces aux donnees d'un choc TCV
%------------------------------------------------------------------------------
% fichier : zuitcvacces.m
% 
% fonction Matlab 5 : 
%	creation du formulaire d'acces aux donnees d'un choc TCV
%  
% syntaxe :  
%	zuitcvacces
%  
% entrees :  
%  
% sorties :  
%  
% fonction ecrite par J-F Artaud , poste 46-78
% version  3.0  du  13/12/2004  
%  
% liste des modifications :  
%  
%------------------------------------------------------------------------------ 
%  
function zuitcvacces

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
	zassignin('base','prepare.numchoc',8001);
end

try
    racine = evalin('base','prepare.racine');
catch
	if strcmp(getenv('USER'),'cgc')
		racine ='/usr/drfc/cgc/cgc_data/tcv/data';
	else
		racine = strcat(getenv('HOME'),'/zineb/data/tcv');
	end

end
zassignin('base','prepare.racine',racine);


hform=numform(chemin,racine);
zuifaitnumtcv('init');
zwaitfor(hform)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% intervalle de temps plus scenario
% lecture des donnees du scenario
numchoc  = evalin('base','prepare.numchoc') ;
racine   = evalin('base','prepare.racine') ;
if ~isempty(findstr(racine,int2str(numchoc)))
     racinec = racine;
else
     racinec = sprintf('%s/%d',racine,numchoc);
end
zassignin('base','prepare.racinec',racinec);

if exist(sprintf('%s/tcvtemp.mat',racinec),'file')
  load(sprintf('%s/tcvtemp',racinec));
else
  disp(' no temporal file, run ztcv for this shot')
  return
end
if exist(sprintf('%s/tcvprof.mat',racinec),'file')
  load(sprintf('%s/tcvprof',racinec));
else
  disp(' no profile file, run ztcv for this shot')
  return
end
if 1>2
load(sprintf('%s/tcvpsi',racinec));
load(sprintf('%s/tcveq',racinec));
load(sprintf('%s/tcvfit',racinec));
else
   tcvfit = []; 
end
% creation du formulaire
tdeb     = min(tcvprof.tTe);
tfin     = max(tcvprof.tTe);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
nbmse = 0;
nbpol = 0;
zassignin('base','prepare.nbmse',nbmse);
zassignin('base','prepare.mse',1);
zassignin('base','prepare.tmse',0);
zassignin('base','prepare.nbpol',nbpol);
zassignin('base','prepare.pol',1);
zassignin('base','prepare.tpol',0);
videt=[];
[hout,hui]    = scenarioformnew(tdeb,tfin,numchoc,videt,videt);
hui.axes_plot = zuiplotin(hui.list_plot);
if ~isempty(tcvfit)
   plot(tcvtemp.tip,abs(tcvtemp.ip)/1e5,'r',tcvfit.tnex,tcvfit.nex(:,1)/1e19,'b',...
     tcvtemp.tPec,sum(tcvtemp.Pec,2)/1e6);
else
   plot(tcvtemp.tip,abs(tcvtemp.ip)/1e5,'r',tcvprof.tTe,tcvprof.ne(:,1)/1e19,'b',...
     tcvtemp.tPec,sum(tcvtemp.Pec,2)/1e6);
end
legend('Ip*10','Ne(0)','Pecrh')

title(['Scenario'])
xlabel('time (s)');
ylabel('');
% attente retour
zwaitfor(hout)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end


% creation de la base temps homogene
pas   = mean(diff(tcvprof.tTe));
tdeb  = evalin('base','prepare.tdeb');
tfin  = evalin('base','prepare.tfin');
temps = [tdeb;tfin];
zassignin('base','prepare.vtemps',temps);
delta = [pas;pas];
zassignin('base','prepare.delta',delta);

zassignin('base','prepare.tip',tcvtemp.tip);
zassignin('base','prepare.ip',tcvtemp.ip/1e5);
if ~isempty(tcvfit)
  zassignin('base','prepare.tne',tcvfit.tnex);
  zassignin('base','prepare.ne0',tcvfit.nex(:,1)/1e19);
else
  zassignin('base','prepare.tne',tcvprof.tTe);
  zassignin('base','prepare.ne0',tcvprof.ne(:,1)/1e19);
    
end
zassignin('base','prepare.tpecrh',tcvtemp.tPec);
zassignin('base','prepare.pecrh',tcvtemp.Pec/1e6);

% dialogue pour la base tcv
nom  = sprintf('TCV time database, shot #%g',numchoc);
aide = 'time step edition';
liste_ref = '     Ip    |Ne0|Pecrh|empty';
var_ref   = {{'prepare.tefit','prepare.ip',':'}, ...
	     {'prepare.tnem','prepare.ne0',':'}, ...
	     {'prepare.tplh','prepare.pecrh',':'}, ...
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
temps  = [];                 
vtemps = evalin('base','prepare.vtemps');
delta  = evalin('base','prepare.delta');
for k  = 1:(length(vtemps)-1)
  temps = cat(2,temps,vtemps(k):delta(k):vtemps(k+1));
end
ind = find(diff(temps) <= (min(delta)/10));
if ~isempty(ind)
	temps(ind)=[];
end
zassignin('base','prepare.temps',temps(:));

% dialogue pour les parametres
info          = ztcvacces;
option        = info.valeur; 

% securite Te
temode ={};
if exist(fullfile(racine,sprintf('%d/tex%d.mat',numchoc,numchoc)))
   temode{end+1} = 1;
end 
if exist(fullfile(racine,sprintf('%d/teshthx%d.mat',numchoc,numchoc)))
   temode{end+1} = 2;
end 
if exist(fullfile(racine,sprintf('%d/teshx%d.mat',numchoc,numchoc)));
   temode{end+1} = 3;
end 
zassignin('base','prepare.fite',temode);

% composition 
zassignin('base','prepare.zcompoval',compovalid(numchoc,racine));



zassignin('base','option',option);
h=zuicreefunform('ztcvacces','option',1);
set(h,'name',sprintf('input file for TCV shot %g',evalin('base','prepare.numchoc')));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% appel de la fonction d'acces
hdlg = msgbox('Cronos input file in progress','be patient ...','help');
drawnow;
evalin('base','[cr,data,param]=ztcvacces(prepare.numchoc,prepare.chemin,prepare.temps,option,prepare.racinec);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'Error tcvacces');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('data access trouble','tcvacces');
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
		tps=[num2str(tps1) ' ->' num2str(tps2)];
		zuidata(h.text_temps,tps);
	end
end


% formulaire d'edition du numchoc et du chemin
function hout=numform(chemin,racine)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accestcv');
if ishandle(hfig)
	zuicloseone(hfig) ;
end

% creation du formulaire
form = {}; 
% ligne pour le numero du choc
col1 = {'text_numchoc','text','shot number :',16,'shot number ',''};
col2 = {'numchoc','edit','',26,'shot number ','','prepare.numchoc'};
form{length(form)+1} = {col1,col2};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'text_chemin','text@full','Complete path to access Cronos file :',[],'/usr/drfc/usr/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full',chemin,42,'Complete path to access Cronos file ','','prepare.chemin'};
form{length(form)+1} = {col1};

col1 = {'text_racine','text@full','Complete path to access TCV data file generated by ztcv :',[],'/usr/drfc/cgc/cgc_data/data/tcv/data/',''};
form{length(form)+1} = {col1};
col1 = {'racine','edit@full',racine,42,'Complete path to access TCV data file','','prepare.racine'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'error message',''};
form{length(form)+1} = {col1};

hout=zuicreeform('TCV data access','accesTCV','zuifaitnumtcv','zuictrlnumtcv',form);

