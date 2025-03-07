%  ZUIFTUACCES  formulaire d'acces aux donnees d'un choc ftu
%------------------------------------------------------------------------------
% fichier : zuiftuacces.m
% 
% fonction Matlab 5 : 
%	creation du formulaire d'acces aux donnees d'un choc ftu
%  
% syntaxe :  
%	zuiftuacces
%  
% entrees :  
%  
% sorties :  
%  
% fonction écrite par J-F Artaud , poste 46-78
% version  2.1  du  4/06/2003  
%  
% liste des modifications :  
%  
%------------------------------------------------------------------------------ 
%  
function zuiftuacces

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
	zassignin('base','prepare.numchoc',19103);
end

try
    racine = evalin('base','prepare.racine');
catch
	if strcmp(getenv('USER'),'cgc')
		racine ='/usr/drfc/cgc/cgc_data/ftu/data';
	else
		racine = strcat(getenv('HOME'),'/zineb/data/ftu');
	end

end
zassignin('base','prepare.racine',racine);


hform=numform(chemin,racine);
zuifaitnumftu('init');
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

load(sprintf('%s/ftutemp',racinec));
load(sprintf('%s/ftupsi',racinec));
load(sprintf('%s/ftuprof',racinec));
load(sprintf('%s/ftueq',racinec));
load(sprintf('%s/ftufit',racinec));

% creation du formulaire
tdeb     = min(ftuprof.tne);
tfin     = max(ftuprof.tne);
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
plot(ftutemp.tip,abs(ftutemp.ip)/1e5,'r',ftufit.tnex,ftufit.nex(:,1)/1e19,'b',...
     ftutemp.tplh,ftutemp.plh/1e6);
legend('Ip*10','Ne(0)','PLH')

title(['Scenario'])
xlabel('temps (s)');
ylabel('');
% attente retour
zwaitfor(hout)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end


% creation de la base temps homogene
pas   = mean(diff(ftuprof.tne));
tdeb  = evalin('base','prepare.tdeb');
tfin  = evalin('base','prepare.tfin');
temps = [tdeb;tfin];
zassignin('base','prepare.vtemps',temps);
delta = [pas;pas];
zassignin('base','prepare.delta',delta);

zassignin('base','prepare.tip',ftutemp.tip);
zassignin('base','prepare.ip',ftutemp.ip/1e5);
zassignin('base','prepare.tne',ftufit.tnex);
zassignin('base','prepare.ne0',ftufit.nex(:,1)/1e19);
zassignin('base','prepare.tplh',ftutemp.tplh);
zassignin('base','prepare.plh',ftutemp.plh/1e6);

% dialogue pour la base ftu
nom  = sprintf('Edition de la basetemps pour le choc FTU #%g',numchoc);
aide = 'Edition du pas de temps en fonction du temps';
liste_ref = '     Ip    |Ne0|Plh|vide';
var_ref   = {{'prepare.tefit','prepare.ip',':'}, ...
	     {'prepare.tnem','prepare.ne0',':'}, ...
	     {'prepare.tplh','prepare.plh',':'}, ...
	     {'[]','[]',''}};
hout=zuieditcons(nom,aide,temps,delta,'temps (s)','pas de temps (s)', ...
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
info          = zftuacces;
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
h=zuicreefunform('zftuacces','option',1);
set(h,'name',sprintf('Preparation du choc FTU %g',evalin('base','prepare.numchoc')));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% appel de la fonction d'acces
hdlg = msgbox('Creation de la structure de donnees Cronos en cours','Patience ...','help');
drawnow;
evalin('base','[cr,data,param]=zftuacces(prepare.numchoc,prepare.chemin,prepare.temps,option,prepare.racinec);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'Erreur ftuacces');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('Probleme lors de la lecture des donnees','ftuacces');
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
function hout=numform(chemin,racine)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesftu');
if ishandle(hfig)
	zuicloseone(hfig) ;
end

% creation du formulaire
form = {}; 
% ligne pour le numero du choc
col1 = {'text_numchoc','text','numero du choc :',16,'numero du choc ',''};
col2 = {'numchoc','edit','',26,'numero du choc ','','prepare.numchoc'};
form{length(form)+1} = {col1,col2};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'text_chemin','text@full','Chemin complet d''acces aux fichiers Cronos :',[],'/usr/drfc/usr/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full',chemin,42,'Chemin complet d''acces aux fichiers ','','prepare.chemin'};
form{length(form)+1} = {col1};

col1 = {'text_racine','text@full','Chemin complet d''acces aux fichiers de donnees ftu :',[],'/usr/drfc/cgc/cgc_data/data/ftu/data/',''};
form{length(form)+1} = {col1};
col1 = {'racine','edit@full',racine,42,'Chemin complet d''acces aux fichiers de donnees ftu','','prepare.racine'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'messages d''erreur',''};
form{length(form)+1} = {col1};

hout=zuicreeform('Acces aux donnees FTU','accesftu','zuifaitnumftu','zuictrlnumftu',form);

