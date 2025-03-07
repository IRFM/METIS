%  ZUIJETACCES  formulaire d'acces aux donnees d'un choc JET
%------------------------------------------------------------------------------
% fichier : zuijetacces.m
% 
% fonction Matlab 5 :
%	creation du formulaire d'acces aux donnees d'un choc JET
% version sans waitfor
%  
% syntaxe :  
%	zuijetacces({phase})
%  
% entrees :  
%
%   phase  = phase de deroulement de la porcedure 1 a 5
%  
% sorties :  
%  
% fonction ecrite par J-F Artaud , poste 46-78
% version  2.0  du  27/11/2002  
%  
% liste des modifications :  
%  
%------------------------------------------------------------------------------ 
%  
function zuijetacces2(phase);

% ce progamme n'est pas en mode modale
%
if nargin < 1
   phase =1;
end
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
		racine ='/usr/drfc/cgc/cgc_data/jet/data';
	else
		racine = strcat(getenv('HOME'),'/zineb/data/JET');
	end
%	racine ='/usr/drfc/cgc/matlab5/tcron/JET/data/';
end
zassignin('base','prepare.racine',racine);

if phase == 1
   hform=numform2(chemin,racine);
   zuifaitnumjet('init');
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

datajet  = load(sprintf('%s/temp%d',racinec,numchoc));
dataefit = load(sprintf('%s/efit%d',racinec,numchoc));
dataMSE.tqEFTMx  = [];
if exist(sprintf('%s/qMSEx%d.mat',racinec,numchoc))
  dataMSE=load(sprintf('%s/qMSEx%d.mat',racinec,numchoc));
end
dataPOL.tqPOLx  = [];
if exist(sprintf('%s/qPOLx%d.mat',racinec,numchoc))
  dataPOL=load(sprintf('%s/qPOLx%d.mat',racinec,numchoc));
end

% creation du formulaire
if phase == 2
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
    
   [hout,hui]    = scenarioformnew2(tdeb,tfin,numchoc,dataMSE.tqEFTMx,dataPOL.tqPOLx);
   hui.axes_plot = zuiplotin(hui.list_plot);
   
   if length(datajet.tplh)<1 
       datajet.tplh='n';
   end  
   if length(datajet.tpicrh)<1 
       datajet.tpicrh='n';
   end  
   if length(datajet.tpnbi)<1 
       datajet.tpnbi='n';
   end  
   
   if ~ischar(datajet.tplh(1)) & ~ischar(datajet.tpicrh(1))  & ~ischar(datajet.tpnbi(1))
      plot(dataefit.tefit,abs(dataefit.Ip)/1e5,'r',datajet.tnem,datajet.nem/1e19,'b',...
      datajet.tplh,datajet.plh/1e6,'c',datajet.tpicrh,datajet.picrh/1e6,'m',...
      datajet.tpnbi,datajet.pnbi/1e6,'k');
      legend('Ip*10','Nmoy','Plh','Picrh','Pnbi')
   elseif ~ischar(datajet.tplh(1)) & ~ischar(datajet.tpicrh(1))
      plot(dataefit.tefit,abs(dataefit.Ip)/1e5,'r',datajet.tnem,datajet.nem/1e19,'b',...
      datajet.tplh,datajet.plh/1e6,'c',datajet.tpicrh,datajet.picrh/1e6,'m');
      legend('Ip*10','Nmoy','Plh','Picrh')
   elseif ~ischar(datajet.tplh(1))
      plot(dataefit.tefit,abs(dataefit.Ip)/1e5,'r',datajet.tnem,datajet.nem/1e19,'b',...
      datajet.tplh,datajet.plh/1e6,'c');
      legend('Ip*10','Nmoy','Plh')
   else
      plot(dataefit.tefit,abs(dataefit.Ip)/1e5,'r',datajet.tnem,datajet.nem/1e19,'b');
      legend('Ip*10','Nmoy')
   end
   if ~isempty(dataMSE.tqEFTMx)  & ~isempty(dataPOL.tqPOLx)
      hold on
      v = axis;
      for kEFTM=1:length(dataMSE.tqEFTMx)
         plot([dataMSE.tqEFTMx(kEFTM) dataMSE.tqEFTMx(kEFTM)],[0 v(4)/2],'b:')
      end
      for kPOL=1:length(dataPOL.tqPOLx)
         plot([dataPOL.tqPOLx(kPOL) dataPOL.tqPOLx(kPOL)],[v(4)/2 v(4)],'r--')
      end
  title(['Scheme, + MSE and Polarimetry times'])
      hold off
   elseif ~isempty(dataMSE.tqEFTMx) 
      hold on
      v = axis;
      for kEFTM=1:length(dataMSE.tqEFTMx)
         plot([dataMSE.tqEFTMx(kEFTM) dataMSE.tqEFTMx(kEFTM)],[0 v(4)/2],'b:')
      end
   title(['Scheme, + MSE times'])
      hold off
   
   elseif ~isempty(dataPOL.tqPOLx)
      hold on
      v = axis;
%      for kPOL=1:length(dataPOL.tqPOLx)
%         plot([dataPOL.tqPOLx(kPOL) dataPOL.tqPOLx(kPOL)],[0 v(4)/2],'r--')
%      end
      sxmX=[dataPOL.tqPOLx';dataPOL.tqPOLx'];
      sxmY=[ones(1,length(dataPOL.tqPOLx))*v(4)/2;ones(1,length(dataPOL.tqPOLx))*v(4)];
      line(sxmX,sxmY,'Color','r','LineStyle','--');
      title(['Scheme, + Polarimetry times'])
      hold off
   else
      title(['Scheme'])
   end
   xlabel('time (s)');
   ylabel('');
	return
end


% creation de la base temps homogene
% dialogue pour la base jet
if phase == 3
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


   nom  = sprintf('JET time database, shot #%g',numchoc);
   aide = 'time step edition';
   liste_ref = '     Ip    |Nl|Plh|Picrh|Pnbi|empty';
   var_ref   = {{'prepare.tefit','prepare.ip',':'}, ...
   	     {'prepare.tnem','prepare.nem',':'}, ...
   	     {'prepare.tplh','prepare.plh',':'}, ...
   	     {'prepare.tpicrh','prepare.picrh',':'}, ...
   	     {'prepare.tpnbi','prepare.pnbi',':'}, ...
   	     {'[]','[]',''}};
   hout=zuieditcons(nom,aide,temps,delta,'time (s)','time step (s)', ...
                    'prepare.vtemps','prepare.delta',1,'zuijetacces2(4)', ...
                    liste_ref,var_ref);
	return
end

if phase  == 4   % debut phase 4
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
   h=zuicreefunform('zjetacces','option',1,0,'zuijetacces2(5);');
   set(h,'name',sprintf('input file creation for JET shot %g',evalin('base','prepare.numchoc')));
	return
end

% ici commence la phase 5
if phase  == 5
   % on affiche la fenetre specifique Ex-file, si celle-ci est demandee
   okex = evalin('base','option.exfile');
   if okex
     info = zjetexfile;
     optionex = info.valeur;
     zassignin('base','optionex',optionex);
     h=zuicreefunform('zjetexfile','optionex',1,0,'zuijetacces2(6);');
     set(h,'name',sprintf('JET ExFile data choice  #%g',evalin('base','prepare.numchoc')));
	  return
   else
     info = zjetexfile;    % si l'ex file n'est pas demandee
     optionex = info.valeur; % met a zero (valeur par defaut) toutes les valeurs d'optionex
     zassignin('base','optionex',optionex);
   end
end

% Ici commence la phase 6
% appel de la fonction d'acces aux donnees crees par zjet
%option = evalin('base','option');  % pour garder les choix effectues par l'utilisateur (remis a 0 par l'appel a zjetacces ligne 216)
hdlg = msgbox('Cronos input file in progress','Be Patient ...','help');
drawnow;
evalin('base','[cr,data,param]=zjetacces(prepare.numchoc,prepare.chemin,prepare.temps,option,optionex,prepare.racinec);');
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
   warndlg('Data access trouble','JETacces');
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
function hout=numform2(chemin,racine)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesJET');
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

col1 = {'text_chemin','text@full','Complete path to access Cronos file (input/output) :',[],'/usr/drfc/usr/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full',chemin,42,'Complete path to access Cronos file ','','prepare.chemin'};
form{length(form)+1} = {col1};

col1 = {'text_racine','text@full','Complete path to access JET data file generated by zjet :',[],'/usr/drfc/cgc/matlab5/tcron/JET/data/',''};
form{length(form)+1} = {col1};
col1 = {'racine','edit@full',racine,42,'Complete path to access JET data file Jet','','prepare.racine'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'messages d''erreur',''};
form{length(form)+1} = {col1};

hout=zuicreeform('Data JET access','accesJET','zuifaitnumjet2','zuictrlnumjet',form);
