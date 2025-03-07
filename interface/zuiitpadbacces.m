%  ZUIITPAACCES  formulaire d'acces aux donnees des bases ITPA Profils
%------------------------------------------------------------------------------
% fichier : zuiitpadbacces.m
% 
% fonction Matlab 5 : 
%	formulaire d'acces aux donnees des bases ITPA Profils
% version sans waitfor
%  
% syntaxe :  
%	zuiitpadbacces({phase})
%  
% entrees :  
%
%   phase  = phase de deroulement de la procedure 1 a 5
%  
% sorties :  
%  
% fonction ecrite par F. Imbeaux
% version  2.1  du  24/07/2003  
%  
% liste des modifications :  
%  
%------------------------------------------------------------------------------ 
%  
function zuiitpadbacces(phase);

% ce progamme n'est pas en mode modal
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
	zassignin('base','prepare.numchoc',53521);
end

%%%%%%%%%%%%%%%%%%%%%%% Phase 1 : selection du numero de choc, tokanak, base %%%%%%%%%%%%%%%%%%
if phase == 1        
   hform=numform2(chemin);  % definition du formulaire en fin de fichier !
   zuifaitnumitpa('init');
   return
end


%%%%%%%%%%%%%%%%%%%%%%% Phase > 1 %%%%%%%%%%%%%%%%%%

% lecture des donnees du scenario
numchoc  = evalin('base','prepare.numchoc') ;
base = evalin('base','prepare.base') ;
tokname = evalin('base','prepare.tokname') ;

%%%%%%%%%%%%%%%%%%%%%%% Phase 2 : selection de la base temps %%%%%%%%%%%%%%%%%%
% creation du formulaire
if phase == 2
   disp('Connecting to ITPA DB server ...')
   hdlg = msgbox('Connecting to ITPA DB server ...','Remote access','help');
   drawnow;
   %mdsconnect('tokamak-profiledb.ukaea.org.uk');  % connexion au serveur ITPA
   mdsconnect('194.81.223.118');        % connexion au serveur ITPA
   disp('Opening the requested shot number ...')
   eval(['[a,b]=mdsopen(''',base,tokname,''',',int2str(numchoc),');']) % ouverture de l'arbre correspondant au choc
   if ishandle(hdlg)
      zuicloseone(hdlg);
   end
      
   [dummy,time] = cgcgetitpa_1d('amin'); % lecture du temps pour avoir tmin et tmax des donnees de la base
   tdeb = min(time);
   tfin = max(time);
   dt = 0.05; % valeur par defaut du pas de temps (s)
   zassignin('base','prepare.tdeb',tdeb);
   zassignin('base','prepare.tfin',tfin);
   zassignin('base','prepare.dt',dt);    

   [hout,hui]    = scenarioformitpa(tdeb,tfin,dt,numchoc);
   hui.axes_plot = zuiplotin(hui.list_plot);

   % lecture des donnees correspondant au scenario (pour affichage seulement)
   ip = cgcgetitpa_1d('ip',time);
   pnbi = cgcgetitpa_1d('pnbi',time);
   if ischar(pnbi)
      pnbi = zeros(size(time));
   end   
   picrh = cgcgetitpa_1d('picrh',time);
   if ischar(picrh)
      picrh = zeros(size(time));
   end   
   pecrh = cgcgetitpa_1d('pech',time);
   if ischar(pecrh)
      pecrh = zeros(size(time));
   end   
   plh = cgcgetitpa_1d('plh',time);
   if ischar(plh)
      plh = zeros(size(time));
   end   
   nbar = cgcgetitpa_1d('nel',time);
   
   plot(time,ip/1e5,time,nbar/1e19,time,pnbi/1e6,time,picrh/1e6,time,plh/1e6,time,pecrh/1e6)
   legend('Ip*10 (MA)','Nbar/1e19 (m^-^3)','Pnbi (MW)','Picrh (MW)','Plh (MW)','Pecrh (MW)')

   title('Scenario')
   xlabel('time (s)');
   return
end

%%%%%%%%%%%%%%%%%%%%%%% Phase 3 : affichage du formulaire de preparation de la simulation %%%%%%%%%%%%%%%%%%
if phase == 3
   % formulaire pour le choix des parametres de la simulation
   info          = zitpadbacces;
   option        = info.valeur; 
   zassignin('base','option',option);
   h=zuicreefunform('zitpadbacces','option',1,0,'zuiitpadbacces(4);');
   set(h,'name',sprintf('Preparation du choc %g',evalin('base','prepare.numchoc')));
   return
end

%%%%%%%%%%%%%%%%%%%%%%% Phase 4 : lecture du formulaire de preparation et creation de la structure de donnees Cronos %%%%%%%%%%%%%%%%%%
if phase  == 4
   % creation de la base temps uniforme
   tdeb = evalin('base','prepare.tdeb');
   tfin = evalin('base','prepare.tfin');
   dt   = evalin('base','prepare.dt');
   
   temps = [tdeb:dt:tfin];  % base temps sur laquelle toutes les donnees vont etre reinterpolees
   zassignin('base','prepare.temps',temps);  % on recopie dans le base workspace 

   hdlg = msgbox('Downloading data from ITPA DB ...','Remote access','help');
   drawnow;
   % chargement des donnes, dans le base workspace -> remplissage des structures data et param
   evalin('base','[cr,data,param]=zitpadbacces(prepare.tokname,prepare.numchoc,prepare.chemin,prepare.temps,option);');
   if ishandle(hdlg)
      zuicloseone(hdlg);
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
		tps=[num2str(tps1) ' -> ' num2str(tps2)];
		zuidata(h.text_temps,tps);
	end
   end
end

% formulaire d'edition du numchoc et du chemin
function hout=numform2(chemin)

% on ferme le formulaire s'il existe
[hfig,h] = zuiformhandle('accesITPADB');
if ishandle(hfig)
	zuicloseone(hfig) ;
end

% creation du formulaire
form = {}; 
% ligne pour le numero du choc
col1 = {'text_tokname','text','tokname :',16,'Tokamak name ',''};
col2 = {'tokname','popup','JET|TFTR|DIII-D|RTP|ASDEX-U|C-MOD|FTU|ITER|JT60-U|MAST|T-10|TORE SUPRA|TEXTOR',1,'Tokamak name ',{'jet','tftr','d3d','rtp','aug','cmod','ftu','iter','jt60u','mast','t10','ts','txtr'},'prepare.tokname'};
form{length(form)+1} = {col1,col2};

col1 = {'text_numchoc','text','numero du choc :',16,'numero du choc ',''};
col2 = {'numchoc','edit','',26,'numero du choc ','','prepare.numchoc'};
form{length(form)+1} = {col1,col2};

col1 = {'text_base','text','base :',16,'base ',''};
col2 = {'base','popup','Profile DB|ITB DB',1,'nom de la base ',{'','itb_'},'prepare.base'};
form{length(form)+1} = {col1,col2};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};
 
col1 = {'text_chemin','text@full','Chemin complet d''acces aux fichiers Cronos :',[],'/usr/drfc/usr/zineb/data',''};
form{length(form)+1} = {col1};
col1 = {'chemin','edit@full',chemin,42,'Chemin complet d''acces aux fichiers ','','prepare.chemin'};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

col1 = {'etat','text@full',' ',10,'messages d''erreur',''};
form{length(form)+1} = {col1};

hout=zuicreeform('Acces aux donnees ITPA Profile DBs','accesITPADB','zuifaitnumitpa','zuictrlnumitpa',form);
