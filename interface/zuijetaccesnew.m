%  ZUIJETACCES  formulaire d'acces aux donnees d'un choc JET
%------------------------------------------------------------------------------
% fichier : zuijetacces.m
% 
% fonction Matlab 5 : 
%	creation du formulaire d'acces aux donnees d'un choc JET
%  
% syntaxe :  
%	zuijetacces
%  
% entrees :  
%  
% sorties :  
%  
% fonction écrite par J-F Artaud , poste 46-78
% version  1.9  du  3/07/2002  
%  
% liste des modifications :  
%	* 27/09/2001 -> on remplace les wardlg par des msgbox avec un icon 'help
%	* 26/06/2002 -> test de l'existence des fichiers Te (utilisation de fullfile)
%	* 3/07/2002 ->  Par defaut, lecture des donnees JET sous le directoire home
%  
%------------------------------------------------------------------------------ 
%  
function zuijetacces

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
	zassignin('base','prepare.numchoc',51782);
end

try
    racine = evalin('base','prepare.racine');
catch
	if strcmp(getenv('USER'),'cgc')
		racine ='/usr/drfc/cgc/cgc_data/jet/data/';
	else
		racine = strcat(getenv('HOME'),'/zineb/data/JET/');
	end
%	racine ='/usr/drfc/cgc/matlab5/tcron/JET/data/';
end
zassignin('base','prepare.racine',racine);


hform=numform(chemin,racine);
zuifaitnumjet('init');
zwaitfor(hform)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% intervalle de temps plus scenario
% lecture des donnees du scenario
numchoc  = evalin('base','prepare.numchoc') ;
racine   = evalin('base','prepare.racine') ;
datajet  = load(sprintf('%s/%d/temp%d',racine,numchoc,numchoc));
dataefit = load(sprintf('%s/%d/efit%d',racine,numchoc,numchoc));
dataMSE.tqEFTMx  = [];
if exist(sprintf('%s/%d/qMSEx%d.mat',racine,numchoc,numchoc))
  dataMSE=load(sprintf('%s/%d/qMSEx%d.mat',racine,numchoc,numchoc));
end
dataPOL.tqPOLx  = [];
if exist(sprintf('%s/%d/qPOLx%d.mat',racine,numchoc,numchoc))
  dataPOL=load(sprintf('%s/%d/qPOLx%d.mat',racine,numchoc,numchoc));
end

% creation du formulaire
tdeb = min(dataefit.tefit);
tfin = max(dataefit.tefit);
zassignin('base','prepare.tdeb',tdeb);
zassignin('base','prepare.tfin',tfin);
if ~isempty(dataMSE.tqEFTMx)
    nbmse  = length(dataMSE.tqEFTMx);
else
    nbmse  = 0;
end
if ~isempty(dataPOL.tqPOLx)
    nbpol  = length(dataPOL.tqPOLx);
else
    nbpol  = 0;
end

zassignin('base','prepare.nbmse',nbmse);
zassignin('base','prepare.mse',1);
zassignin('base','prepare.tmse',0);
zassignin('base','prepare.nbpol',nbpol);
zassignin('base','prepare.pol',1);
zassignin('base','prepare.tpol',0);
    
[hout,hui]    = scenarioformnew(tdeb,tfin,numchoc,dataMSE.tqEFTMx,dataPOL.tqPOLx);
hui.axes_plot = zuiplotin(hui.list_plot);
plot(dataefit.tefit,abs(dataefit.Ip)/1e5,'r',datajet.tnem,datajet.nem/1e19,'b',...
     datajet.tplh,datajet.plh/1e6,'c',datajet.tpicrh,datajet.picrh/1e6,'m',...
     datajet.tpnbi,datajet.pnbi/1e6,'k');
legend('Ip*10','Nmoy','Plh','Pfci','Pnbi')
if ~isempty(dataMSE.tqEFTMx)  & ~isempty(dataPOL.tqPOLx)
  hold on
  v = axis;
  for kEFTM=1:length(dataMSE.tqEFTMx)
    plot([dataMSE.tqEFTMx(kEFTM) dataMSE.tqEFTMx(kEFTM)],[0 v(4)/2],'b:')
  end
  for kPOL=1:length(dataPOL.tqPOLx)
    plot([dataPOL.tqPOLx(kPOL) dataPOL.tqPOLx(kPOL)],[v(4)/2 v(4)],'r--')
  end
  title(['Scenario, + temps des reconstructions MSE et Pol'])
  hold off
elseif ~isempty(dataMSE.tqEFTMx) 
  hold on
  v = axis;
  for kEFTM=1:length(dataMSE.tqEFTMx)
    plot([dataMSE.tqEFTMx(kEFTM) dataMSE.tqEFTMx(kEFTM)],[0 v(4)/2],'b:')
  end
  title(['Scenario, + temps des reconstructions MSE '])
  hold off

elseif ~isempty(dataPOL.tqPOLx)
  hold on
  v = axis;
  for kPOL=1:length(dataPOL.tqPOLx)
    plot([dataPOL.tqPOLx(kPOL) dataPOL.tqPOLx(kPOL)],[0 v(4)/2],'r--')
  end
  title(['Scenario, + temps des reconstructions Pol'])
  hold off
else
  title(['Scenario'])
end
xlabel('temps (s)');
ylabel('');
% attente retour
zwaitfor(hout)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end


% creation de la base temps homogene
pas= mean(diff(dataefit.tefit));
tdeb = evalin('base','prepare.tdeb');
tfin = evalin('base','prepare.tfin');
temps = [tdeb;tfin];
zassignin('base','prepare.vtemps',temps);
delta = [pas;pas];
zassignin('base','prepare.delta',delta);

zassignin('base','prepare.tefit',dataefit.tefit);
zassignin('base','prepare.ip',dataefit.Ip/1e5);
zassignin('base','prepare.tnem',datajet.tnem);
zassignin('base','prepare.nem',datajet.nem/1e19);
zassignin('base','prepare.tplh',datajet.tplh);
zassignin('base','prepare.plh',datajet.plh/1e6);
zassignin('base','prepare.tpicrh',datajet.tpicrh);
zassignin('base','prepare.picrh',datajet.picrh/1e6);
zassignin('base','prepare.tpnbi',datajet.tpnbi);
zassignin('base','prepare.pnbi',datajet.pnbi/1e6);


% dialogue pour la base jet
nom  = sprintf('Edition de la basetemps pour le choc JET #%g',numchoc);
aide = 'Edition du pas de temps en fonction du temps';
liste_ref = '     Ip    |Nl|Plh|Pfci|Pidn|vide';
var_ref   = {{'prepare.tefit','prepare.ip',':'}, ...
	     {'prepare.tnem','prepare.nem',':'}, ...
	     {'prepare.tplh','prepare.plh',':'}, ...
	     {'prepare.tpicrh','prepare.picrh',':'}, ...
	     {'prepare.tpnbi','prepare.pnbi',':'}, ...
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
info          = zjetacces;
option        = info.valeur; 
valmse        = evalin('base','prepare.mse');
valpol        = evalin('base','prepare.pol');
if valmse == 1
  option.efit     = 1;
  option.tempsdeb = evalin('base','prepare.tmse');
elseif valpol == 1
  option.efit = 2;
  option.tempsdeb = evalin('base','prepare.tpol');
else
  option.efit     = valmse;
end

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
h=zuicreefunform('zjetacces','option',1);
set(h,'name',sprintf('Preparation du choc JET %g',evalin('base','prepare.numchoc')));
zwaitfor(h)
code=getappdata(0,'coderetour');
if ~strcmp(code.action,'validation')
	return
end

% appel de la fonction d'acces
hdlg = msgbox('Creation de la structure de donnees Cronos en cours','Patience ...','help');
drawnow;
evalin('base','[cr,data,param]=zjetacces(prepare.numchoc,prepare.chemin,prepare.temps,option,prepare.racine);');
if ishandle(hdlg)
	zuicloseone(hdlg);
end
% test creation donnees
ok = evalin('base','isstruct(param)&isstruct(data)');
if ok == 0
	warndlg(lasterr,'Erreur JETacces');
	return
end
cr = evalin('base','cr','NaN');
if cr ~= 0
	warndlg('Probleme lors de la lecture des donnees','JETacces');
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
		tps=[num2str(tps1) ' à ' num2str(tps2)];
		zuidata(h.text_temps,tps);
	end
end


% formulaire d'edition du numchoc et du chemin
function hout=numform(chemin,racine)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesJET');
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

col1 = {'text_racine','text@full','Chemin complet d''acces aux fichiers de donnees Jet :',[],'/usr/drfc/cgc/matlab5/tcron/JET/data/',''};
form{length(form)+1} = {col1};
col1 = {'racine','edit@full',racine,42,'Chemin complet d''acces aux fichiers de donnees Jet','','prepare.racine'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'messages d''erreur',''};
form{length(form)+1} = {col1};

hout=zuicreeform('Acces aux donnees JET','accesJET','zuifaitnumjet','zuictrlnumjet',form);

