function [cr,data,param]=videacces(vparam,vdata,temps,option)

% declaration des parametres
if nargin <=1 
	valeur.psimode     = 2;       % mode Psi : 1 -> interpretatif, 2 -> predictif
	valeur.pemode      = 1;       % mode Pe : 1 -> interpretatif, 2 -> predictif
	valeur.pionmode    = 1;       % mode Pion : 1 -> interpretatif, 2 -> predictif
	valeur.nelmode     = 1;       % mode Ne : 1 -> interpretatif, 2 -> predictif
	valeur.psilim      = 0;       % condition au limite sur Psi : 0 -> Ip, 1 -> Vloop, 2 -> Psi(bord)
   valeur.lambda      = 0;       % facteur lambda dans les flux de chaleur provenant du flux de particules {0,3/.2,5/2}
  	valeur.modecoef    = 1;       % mode de calcul des coefficients des equations de transport  0 -> tous, 1 -> convectif +diagonaux
   %valeur.self        = 0;       % mode de fonctionnement completement  auto consistant si = 1, 0 -> coefficient +neoclassique auto consistant 
  	valeur.fast        = 0;       % 0 -> mode standart correct, 1-> pas de sources neoclassique, ni de coefficient neo self consistante, 2 -> mode optimiser pour la diffusion du courant
   valeur.source_bord = 0;       % controle du clacul de recyclage  : 1 = recyclage recalculer a chaque sous pas de temps, 0 = comme les sources
   valeur.cn          = 0.5;     % mode solveur : 0 -> implicite, 0.5 -> Cranck-Nickolson
   valeur.plotonoff   = 0;       % mode interractif avec plot si 1 
   valeur.verbose     = 1;       % 1-> commentaires , 0 -> pas de commentaires
	valeur.rebuilt     = 1;       % 1 -> reconstruction automatique du fichier resultat en fin d'execution 
	valeur.post        = 0;       % 1 -> postprocessing en fin d'execution
   
	type.psimode     = 'integer';       
	type.pemode      = 'integer';       
	type.pionmode    = 'integer';       
	type.nelmode     = 'integer';       
	type.psilim      = 'integer';       
	type.lambda      = 'integer';     
	type.modecoef    = 'integer';      
	%type.self        = 'integer';      
	type.fast        = 'integer';       
	type.source_bord = 'integer'; 
	type.cn          = 'float';
	type.plotonoff   = 'integer';       
	type.verbose     = 'integer'; 
	type.rebuilt     = 'integer';       
	type.post        = 'integer'; 
   
	borne.psimode     = {1,2};       
	borne.pemode      = {1,2};       
	borne.pionmode    = {1,2};      
	borne.nelmode     = {1,2};       
	borne.psilim      = {0,1,2};      
	borne.lambda      = {0,3/2,5/2};     
	borne.modecoef    = {0,1};      
	%borne.self        = {-1,0,1};       
	borne.fast        = {0,1,2};       
	borne.source_bord = {0,1};       
	borne.cn          = {0,1/2};
	borne.plotonoff   = {0,1};      
	borne.verbose     = {0,1};       
 	borne.rebuilt     = {0,1};      
	borne.post        = {0,1};       
   
	defaut.psimode     = 2;       
	defaut.pemode      = 1;       
	defaut.pionmode    = 1;       
	defaut.nelmode     = 1;       
	defaut.psilim      = 0;       
	defaut.lambda      = 0;     
	defaut.modecoef    = 1;      
	%defaut.self        = 0;      
	defaut.fast        = 0;       
	defaut.source_bord = 0;       
	defaut.cn          = 0.5;
 	defaut.plotonoff   = 1;       
	defaut.verbose     = 1;       
 	defaut.rebuilt     = 1;       
	defaut.post        = 0;       
  
	info.psimode     = 'mode Psi : 1 -> interpretatif, 2 -> predictif';
	info.pemode      = 'mode Pe : 1 -> interpretatif, 2 -> predictif';
	info.pionmode    = 'mode Pion : 1 -> interpretatif, 2 -> predictif';
	info.nelmode     = 'mode Ne : 1 -> interpretatif, 2 -> predictif';
	info.psilim      = 'condition au limite sur Psi : 0 -> Ip, 1 -> Vloop, 2 -> Psi(bord)';
	info.lambda      = 'facteur lambda dans les flux de chaleur provenant du flux de particules {0,3/.2,5/2}';
	info.modecoef    = 'mode de calcul des coefficients des equations de transport  0 -> tous, 1 -> convectif +diagonaux';
	%info.self        = 'mode de fonctionnement completement  auto consistant si = 1, 0 -> coefficient +neoclassique auto consistant';
	info.fast        = '0 -> mode standart correct, 1-> pas de sources neoclassique, ni de coefficient neo self consistante, 2 -> mode optimiser pour la diffusion du courant';
	info.source_bord = 'controle du clacul de recyclage  : 1 = recyclage recalculer a chaque sous pas de temps, 0 = comme les sources';
	info.cn          = 'mode solveur : 0 -> implicite, 0.5 -> Cranck-Nickolson';
	info.plotonoff   = 'mode interractif avec plot si = 1 ';
	info.verbose     = '1-> commentaires , 0 -> pas de commentaires';
	info.rebuilt     = '1 -> reconstruction automatique du fichier resultat en fin d''execution'; 
	info.post        = '1 -> postprocessing en fin d''execution';

	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'Module de parametrage pour le prolongement d''une simulation Cronos';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr = sortie;
	return
end

% compatibilite



% cr par defaut
cr =0;

% gestion des entrees
if nargin <4 
	cr = 1;
	disp('Nombre d''arguments incorrect !');
end

 
% valeurs par default
tdebut =min(temps);
tfin = max(temps);

% creation du nom du fichier
chemin  = fileparts(vparam.gene.file);
machine = vparam.from.machine;
numchoc = vparam.from.shot.num;
fichier = fullfile(chemin,strcat('zineb',int2str(fix(numchoc)),lower(machine),int2str(round(rem(numchoc,1)*10)), ...
                 'de',int2str(round(tdebut*1000)),'a',int2str(round(tfin*1000))));

% initialisation de la structure
[cr,data,param]=zinit('',temps,vparam.gene.nbrho,fichier,tdebut,tfin, ...
                         vparam.nombre.fci,vparam.nombre.fce,vparam.nombre.idn, ...
                         vparam.nombre.hyb,vparam.nombre.glacon,vparam.gene.nbg);
                         
if cr ~=0
	return
end

% connexion des modules externes
param.fonction = vparam.fonction;
[cr,data,param] = zconnexion(data,param);
if cr ~=0
	return
end

% gestion de la separatrice
data.geo.R = single(NaN .* ones(size(data.geo.R,1),size(vdata.geo.R,2)));
data.geo.Z = single(NaN .* ones(size(data.geo.Z,1),size(vdata.geo.Z,2)));


% parametre generaux
param.gene  = vparam.gene;

% parametres generaux
param.gene.modecoef     = option.modecoef;   % modele diffussif/convectif - pas de terme non diagonaux
%param.gene.self         = option.self;
param.gene.lambda       = option.lambda;     
param.gene.fast         = option.fast;       
param.gene.source_bord  = option.source_bord;       
param.gene.cn           = option.cn;       
param.gene.verbose      = option.verbose;       
param.gene.rebuilt      = option.rebuilt;     
param.gene.post         = option.post;     

if any([option.nelmode,option.pemode,option.pionmode] > 1)
	param.gene.nbeq_mode    = 4;        
else
	param.gene.nbeq_mode    = 1;       
	
end
% recopie des autres parametres
%param.gene.psiequi     = vparam.gene.psiequi;
%param.gene.adiabatic   = vparam.gene.adiabatic;
%param.gene.delta_adia  = vparam.gene.delta_adia;
%param.gene.nmax        = vparam.gene.nmax;
%param.gene.nequi_ini   = vparam.gene.nequi_ini;
%param.gene.nequi       = vparam.gene.nequi;
%param.gene.dpsi_ini    = vparam.gene.dpsi_ini;
%param.gene.djmoy       = vparam.gene.djmoy;
%param.gene.amorti      = vparam.gene.amorti;
%param.gene.mjmoy       = vparam.gene.mjmoy;
%param.gene.evx_inter   = vparam.gene.evx_inter;
%param.gene.critere     = vparam.gene.critere;
%param.gene.nbsauve     = vparam.gene.nbsauve;
%param.gene.rebuilt     = vparam.gene.rebuilt;
%param.gene.post        = vparam.gene.post;
 
%gene
param.compo      = vparam.compo;
param.fonction   = vparam.fonction;
param.cons       = vparam.cons;
param.plot       = vparam.plot;
param.plot.onoff = option.plotonoff;       
param.split      = vparam.split;
param.from       = vparam.from;
param.asser      = vparam.asser;
param.intervalle = vparam.intervalle; 
param.intervalle.temps_old = param.intervalle.temps; 
param.intervalle.temps(1)  = min(data.gene.temps);
param.intervalle.temps(2)  = max(data.gene.temps);
param.intervalle.calcul(1) = max(min(data.gene.temps),param.intervalle.calcul(1));
param.intervalle.calcul(2) = max(data.gene.temps);

% boucle sur les intervalles
noms = fieldnames(param.intervalle);
for kl = 1:length(noms)
      inter = getfield(param.intervalle,noms{kl});
      inter(1) = max(min(inter),min(data.gene.temps));
      inter(2) = min(max(inter),max(data.gene.temps));
      param.intervalle = setfield(param.intervalle,noms{kl},inter);
end


param.profile    = vparam.profile;


% 2 - remplissage des parametres

% les informations sur les donnees
param.from.creation.date =clock;
[s,whoami] = unix('whoami');
param.from.creation.user = whoami;
param.from.source.desc ={};
param.from.createur = 'zplusacces';
param.from.option  = option;

% indice de calcul
param.gene.kmin = 1;
param.gene.kmax = length(temps);
param.gene.k    = 1;
param.gene.nbt  = length(temps);
param.gene.tdeb = temps(1);
param.gene.t    = temps(1);
param.gene.tfin = temps(end);

% boucle sur les donnees
for k = 1 :length(temps)
	fprintf('.');
	tc = temps(k);
	l  = max(find(vdata.gene.temps <= tc));
        if isempty(l) 
           l = 1;
        elseif l < 1
           l = 1;
        end
	datal = zget1t(vdata,l);
	data  = zput1t(data,k,datal);
	data.gene.temps(k) = tc;
end
fprintf('\n');

% les modes
v1                             = ones(param.gene.nbt,1);
v0                             = zeros(param.gene.nbt,1);
ind                            = param.gene.kmin:1:param.gene.kmax;
von                            = v0;
voff                           = v0;            % mis a zeros
von(ind)                       = v1(ind); 
vlit                           = von;         	% donnees en entree
vcalc                          = 2 .* von;    	% calculees
vcopie                         = 3 .* von;    	% recopie du temps precedent
%
if option.psimode == 1
	data.mode.psi        = vlit    
else
	data.mode.psi        = vcalc;    
end
if option.nelmode == 1
	data.mode.nel        = vlit;                 
else
	data.mode.nel        = vcalc;                 
end
if option.pemode == 1
	data.mode.pe         = vlit;                
else
	data.mode.pe         = vcalc;                
end
if option.pionmode == 1
	data.mode.pion       = vlit;                 
else
	data.mode.pion       = vcalc;                 
end

if option.psilim == 0
	data.mode.cons.psi        = v0;           
elseif option.psilim == 1
	data.mode.cons.psi        = v1;           
else 
	data.mode.cons.psi        = 2 .* v1;           
end

% mode coef
if option.modecoef == 1
	data.mode.ee        =  vcalc;                 
	data.mode.ei        =  voff;                 
	data.mode.en        =  voff;                 
	data.mode.ej        =  voff;                 
	data.mode.ve        =  voff;                 
	data.mode.ep        =  voff;                 
	
	data.mode.ie        =  voff;                 
	data.mode.ii        =  vcalc;                 
	data.mode.in        =  voff;                 
	data.mode.ij        =  voff;                 
	data.mode.vi        =  voff;                 
	data.mode.ip        =  voff;                 
	
	data.mode.ne        =  voff;                 
	data.mode.ni        =  voff;                 
	data.mode.nn        =  vcalc;                 
	data.mode.nj        =  voff;                 
	data.mode.vn        =  vcalc;                 
	
	data.mode.fefe      =  voff;                 
	data.mode.fev       =  voff;                 
	
	data.mode.fifi      =  voff;                 
	data.mode.fiv       =  voff;                 
	
	data.mode.rotc       =  voff;                 
	data.mode.rotv      =  voff;                 
else
	data.mode.ee        =  vcalc;                 
	data.mode.ei        =  vcalc;                 
	data.mode.en        =  vcalc;                 
	data.mode.ej        =  vcalc;                 
	data.mode.ve        =  vcalc;                 
	data.mode.ep        =  vcalc;                 
	
	data.mode.ie        =  vcalc;                 
	data.mode.ii        =  vcalc;                 
	data.mode.in        =  vcalc;                 
	data.mode.ij        =  vcalc;                 
	data.mode.vi        =  vcalc;                 
	data.mode.ip        =  vcalc;                 
	
	data.mode.ne        =  vcalc;                 
	data.mode.ni        =  vcalc;                 
	data.mode.nn        =  vcalc;                 
	data.mode.nj        =  vcalc;                 
	data.mode.vn        =  vcalc;                 
	
	data.mode.fefe      =  vcalc;                 
	data.mode.fev       =  vcalc;                 
	
	data.mode.fifi      =  vcalc;                 
	data.mode.fiv       =  vcalc;                 
	
	data.mode.rotc       =  vcalc;                 
	data.mode.rotv      =  vcalc;                 
end	

% autres
if option.plotonoff == 1
	data.mode.plot          = von;  
else
	data.mode.plot          = voff;  
end

% modification du nom du fichier de sortie
param.gene.file = strcat(param.gene.origine,'_resultat');
param.plot.pause = 0;
[pp,fp,ep,vp]=fileparts(param.gene.file);
param.gene.rapsauve =fullfile(pp,'rapsauve',fp);

% changement du type de fichier 
param.gene.filetype ='source';

% sauvegarde du fichier
% sauvegarde du fichier
if nargin <3
  try
    % compactage des donnees
    data=zreduit(param,data,'compact');
    % sauvegarde
    post=[];
	if  verLessThan('matlab','7.0')
				savedb(param.gene.origine,'param','data','post','-V6');
	else
				savedb(param.gene.origine,'param','data','post');
	end
    %save(param.gene.origine,'param','data','post');
    % compression du fichier
    zgzip(param.gene.origine,'compress');
    %save(param.gene.origine,'param','data');
    data=zreduit(param,data,'uncompact');
  catch
    disp('Probleme lors de la sauvegarde du fichier :')
    disp(lasterr)
    cr =-11005;
    return
  end
end


