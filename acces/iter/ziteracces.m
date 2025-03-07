function [cr,data,param]=ziteracces(numchoc,chemin,temps,option,scenario,separatrice)
%
%
%
%
% version 2.1, 19 aout 2003
% dernieres modifications
% 19 aout 2003 -> consbord.fluxqe, qi et ge mis �1 par d�aut (utilisation de la fonction de bord)
% declaration des parametres
if nargin <=1 
	valeur.nbrho       = 101;     % nombre de points radiaux [101]
	
	valeur.compo       = 4;      % composition du plasma 1 -> H/H, 2 -> D/H, 3 -> D/D  et 4 -> D/T
	valeur.ip          = 17;     % Ip max (MA) [17]
	valeur.li          = 1;      % li de reference pendant la phase de chauffage [1]
	valeur.zeff        = 1.6;    % zeff moyen de reference pendant la phase de chauffage [1.6]
	valeur.ntnd        = 1;      % rapport maximun en T et D dans la decharge [1]
	valeur.pidn        = 33;     % puissance IDN par injecteur [<33]
	valeur.pfci        = 20;     % puissance FCI par antenne en MW [<20]
	valeur.pfce        = 20;     % puissance FCE par antenne en MW [<20]
	valeur.nbar        = 0.8;    % Nbar en fraction de la limite de Greenwald [0.8]
	valeur.forme       = 100;    % facteur de forme du profil de ne [100]
	valeur.phase       = 1;      % phase d'equipement de la machine (fixe la pussance de chauffage disponible) [1]
	
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
	valeur.verbose     = 1;       % 1 -> commentaires , 0 -> pas de commentaires
	valeur.rebuilt     = 1;       % 1 -> reconstruction automatique du fichier resultat en fin d'execution 
	valeur.post        = 0;       % 1 -> postprocessing en fin d'execution
	
	type.nbrho       = 'integer';    
	
	type.compo       = 'integer';
	type.ip          = 'float';
	type.li          = 'float';
	type.zeff        = 'float';
	type.ntnd        = 'float';
	type.pidn        = 'float';
	type.pfci        = 'float';
	type.pfce        = 'float';
	type.nbar        = 'float';
	type.forme       = 'float';
	type.phase       = 'integer';
	
	type.psimode     = 'integer';       
	type.pemode      = 'integer';       
	type.pionmode    = 'integer';       
	type.nelmode     = 'integer';       
	type.psilim      = 'integer';       
	type.lambda      = 'integer';     
	type.modecoef    = 'integer';      
	type.self        = 'integer';      
	type.fast        = 'integer';       
	type.source_bord = 'integer'; 
	type.cn          = 'float';
	type.plotonoff   = 'integer';       
	type.verbose     = 'integer'; 
	type.rebuilt     = 'integer';       
	type.post        = 'integer'; 
   
	borne.nbrho       = {21,51,101,201,501,1001};    
	
   borne.compo       = {1,2,3,4};
	borne.ip          = [0,17];
	borne.li          = [0.5,1.5];
	borne.zeff        = [1.1,2.5];
	borne.ntnd        = [0,1];
	borne.pidn        = [0,33];
	borne.pfci        = [0,20];
	borne.pfce        = [0,20];
	borne.nbar        = [0.1,1.2];
	borne.forme       = [3,1000];
	borne.phase       = {1};
	
	borne.psimode     = {1,2};       
	borne.pemode      = {1,2};       
	borne.pionmode    = {1,2};      
	borne.nelmode     = {1,2};       
	borne.psilim      = {0,1,2};      
	borne.lambda      = {0,3/2,5/2};     
	borne.modecoef    = {0,1};      
	borne.self        = {0,1};       
	borne.fast        = {0,1,2};       
	borne.source_bord = {0,1};       
	borne.cn          = {0,1/2};
	borne.plotonoff   = {0,1};      
	borne.verbose     = {0,1};       
 	borne.rebuilt     = {0,1};      
	borne.post        = {0,1};       
  
	defaut.nbrho       = 101;    
	
	defaut.compo       = 4;     
	defaut.ip          = 17;    
	defaut.li          = 1;     
	defaut.zeff        = 1.6;   
	defaut.ntnd        = 1;   
	defaut.pidn        = 33;    
	defaut.pfci        = 20;    
	defaut.pfce        = 20;    
	defaut.nbar        = 0.8;   
	defaut.forme       = 100;   
	defaut.phase       = 1;     
		
	defaut.psimode     = 2;       
	defaut.pemode      = 1;       
	defaut.pionmode    = 1;       
	defaut.nelmode     = 1;       
	defaut.psilim      = 0;       
	defaut.lambda      = 0;     
	defaut.modecoef    = 1;      
	defaut.self        = 0;      
	defaut.fast        = 0;       
	defaut.source_bord = 0;       
	defaut.cn          = 0.5;
 	defaut.plotonoff   = 0;       
	defaut.verbose     = 1;       
 	defaut.rebuilt     = 1;       
	defaut.post        = 0;       
  
	info.nbrho       = 'nombre de points radiaux [101]';
	
	info.compo       = 'composition du plasma 1 -> H/H, 2 -> D/H, 3 -> D/D  et 4 -> D/T';
	info.ip          = 'Ip max (MA) [17]';
	info.li          = 'li de reference pendant la phase de chauffage [1]';
	info.zeff        = 'zeff moyen de reference pendant la phase de chauffage [1.6]';
	info.ntnd        = 'rapport maximun en T et D dans la decharge [1]';
	info.pidn        = 'puissance IDN par injecteur [<33]';
	info.pfci        = 'puissance FCI par antenne en MW [<20]';
	info.pfce        = 'puissance FCE par antenne en MW [<20]';
	info.nbar        = 'Nbar en fraction de la limite de Greenwald [0.8]';
	info.forme       = 'facteur de forme du profil de ne [100]';
	info.phase       = 'phase d''equipement de la machine (fixe la pussance de chauffage disponible) [1]';
	
	info.psimode     = 'mode Psi : 1 -> interpretatif, 2 -> predictif';
	info.pemode      = 'mode Pe : 1 -> interpretatif, 2 -> predictif';
	info.pionmode    = 'mode Pion : 1 -> interpretatif, 2 -> predictif';
	info.nelmode     = 'mode Ne : 1 -> interpretatif, 2 -> predictif';
	info.psilim      = 'condition au limite sur Psi : 0 -> Ip, 1 -> Vloop, 2 -> Psi(bord)';
	info.lambda      = 'facteur lambda dans les flux de chaleur provenant du flux de particules {0,3/.2,5/2}';
	info.modecoef    = 'mode de calcul des coefficients des equations de transport  0 -> tous, 1 -> convectif +diagonaux';
	info.self        = 'mode de fonctionnement completement  auto consistant si = 1, 0 -> coefficient +neoclassique auto consistant ';
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
	
	sortie.description = 'Module d''acces au donnees de Tore Supra';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr = sortie;
	return
end
setappdata(0,'tokamak_name','ITER')

% compatibilite
nbrho        = option.nbrho;
tdebut       = [];
tfin         = [];
mode_compo   = option.compo;
ip0          = option.ip.*1e6;
li0          = option.li;
zeff0        = option.zeff;
ntnd0        = option.ntnd;
pidn0        = option.pidn.*1e6;
pfci0        = option.pfci.*1e6;
pfce0        = option.pfce.*1e6;
nbar0        = option.nbar; % en fraction de la limite de Greenwald
nfact        = option.forme; % facteur de forme du profil de ne
phase        = option.phase;
volume       = separatrice.vv;

% cr par defaut
cr =0;

% gestion des entrees
if nargin <6
	cr = 1;
	disp('Nombre d''arguments incorrect !');
end


 
 
% valeurs par default
if isempty(temps)
	temps=bile.times;
elseif size(temps,2)>1
	temps=temps';
end
if isempty(tdebut)
	tdebut =min(temps);
end
if isempty(tfin)
	tfin = max(temps);
end
if isempty(nbrho)
	nbrho=length(bile.rhofit);
end



% valeurs pour Iter
switch phase
	
case 1
	nbfci    = 1;
	nbfce    = 1;
	nbhyb    = 1;
	nbidn    = 16;
	nbglacon = 2;
	nbg      = 5;
	
otherwise
	error('cas non implante')
end

% creation du nom du fichier
if isempty(chemin)
	chemin=strcat(getenv('HOME'),'/zineb/data');
end
fichier = strcat(chemin,'/zineb',int2str(fix(numchoc)),'iter',int2str(round(rem(numchoc,1)*10)), ...
                 'de',int2str(round(tdebut*1000)),'a',int2str(round(tfin*1000)));

% initialisation de la structure
[cr,data,param]=zinit('',temps,nbrho,fichier,tdebut,tfin,nbfci,nbfce,nbidn,nbhyb,nbglacon,nbg);
if cr ~=0
	return
end
vt = ones(param.gene.nbt,1);
ve = ones(1,param.gene.nbrho);

% parametres generaux
param.gene.modecoef     = option.modecoef;   % modele diffussif/convectif - pas de terme non diagonaux
%param.gene.self         = option.self;
param.gene.lambda       = option.lambda;     
param.gene.fast         = option.fast;       
param.gene.source_bord  = option.source_bord;       
param.gene.cn           = option.cn;       
param.plot.onoff        = option.plotonoff;       
param.gene.verbose      = option.verbose;     
param.gene.rebuilt      = option.rebuilt;     
param.gene.post         = option.post;     
if any([option.nelmode,option.pemode,option.pionmode] > 1)
	param.gene.nbeq_mode    = 4;        
else
	param.gene.nbeq_mode    = 1;       
	
end

% connexion des modules externes
[cr,data,param] = zconnexion(data,param);
if cr ~=0
	return
end

% 1- preparation des donnees de l'equilibre

% 2 - remplissage des parametres
% les informations sur les donnees
param.from.machine = 'ITER';
param.from.shot.num = numchoc;
param.from.shot.date = [2013 08 15 22 23 23] ;
param.from.shot.info = '';
param.from.creation.date =clock;
[s,whoami] = unix('whoami');
param.from.creation.user = whoami;
param.from.source.desc ={};
% parametre d'echantillnage (cf. zsample.m)
%    1 - pour un signal simple :
signal.ondelette        = 0;  %1
signal.defaut.temps     = NaN;
signal.defaut.espace    = 0;
signal.defaut.inf       = [];
signal.plus             = 0;

%    2 - pour un groupe de signaux :
groupe.ondelette         = 0;
groupe.energie           = 1;  %0.01;
groupe.defaut.temps     = NaN;
groupe.defaut.espace    = 0;
groupe.defaut.inf       = [];
groupe.plus             = 0;
param.from.sample.signal   = signal;      % parametre d'echantillonage  pour un signal simple =s(temps,1)
param.from.sample.groupe   = groupe;      % parametre d'echantillonage  pour un groupe =g(temps,espace)
param.from.createur = 'ziteracces';
param.from.option  = option;

% la composition du plasma
% on commence par mettre des 0
param.compo.z = zeros(size(param.compo.z));
param.compo.a = zeros(size(param.compo.a));
switch mode_compo 
case 1
	% plasma H/H
	param.compo.z = [1,2,2,4,18];
	param.compo.a = [1,3,4,9,40];
	c1            = 0.01;
	c2            = 0.03;
	rimp          = 0.1;
case 2
	% plasma D/H
	param.compo.z = [1,1,2,4,18];
	param.compo.a = [2,1,4,9,40];
	c1            = 0.1;
	c2            = 0.04;
	rimp          = 0.1;
case 3
	% plasma D/D
	param.compo.z = [1,2,2,4,18];
	param.compo.a = [2,3,4,9,40];
	c1            = 0.01;
	c2            = 0.03;
	rimp          = 0.1;
case 4
	% plasma D/T
	param.compo.z = [1,1,2,4,18];
	param.compo.a = [2,3,4,9,40];
	c1            = 1;
	c2            = 0.04;
	rimp          = 0.1;
otherwise
	error('cas non implante ...')
end

if strcmp(param.fonction.impur,'zinebcompo')
	param.cons.impur.cmin1    =  c1;
	param.cons.impur.cmin2    =  c2;
	param.cons.impur.rimp     =  rimp;
	param.cons.impur.zeff     =  1;
end

%
% le signe est positif si le courant (ou le champ) est dans le sens trigonometrique
% lorsque le tokamak est regarde depuis le haut
%
param.gene.signe.ip  = -1;              % signe du courant plasma
param.gene.signe.b0  = -1;              % signe du champ toroidal

% 3 - remplissage des signaux simple (pas de rho)

% le temps
data.gene.temps = temps;

% la geometrie
data.geo.r0     = separatrice.r0;
data.geo.z0     = separatrice.z0;
data.geo.a      = separatrice.a;
data.geo.e1     = separatrice.e1;
data.geo.b0     = separatrice.b0;
data.geo.trh1   = separatrice.trh;
data.geo.trb1   = separatrice.trl;
data.geo.ind1   = zeros(size(temps));
data.geo.mode   = 2.*ones(size(temps));
data.geo.R      = single(separatrice.R);
data.geo.Z      = single(separatrice.Z);

% consignes
data.cons.ip      = ip0    .* medfilt1(interp1(scenario.t,scenario.ip,temps,'linear'),5);
data.cons.zeffm   = zeff0  .* medfilt1(interp1(scenario.t,scenario.zeff,temps,'linear'),5);
data.cons.nhnd    = ntnd0  .* medfilt1(interp1(scenario.t,scenario.ftri,temps,'linear'),5); 
% securite zeff
ind = find(data.cons.zeffm < min(param.compo.z));
if ~isempty(ind)
	data.cons.zeffm(ind) = min(param.compo.z) .* ones(1,length(ind));
end
ind = find(data.cons.zeffm > max(param.compo.z));
if ~isempty(ind)
	data.cons.zeffm(ind) = max(param.compo.z) .* ones(1,length(ind));
end
% consigne nhnd, utilise pour le T

% injection de neutres 
data.cons.idn              = (pidn0  .* medfilt1(interp1(scenario.t,scenario.pidn,temps,'linear'),5)) * ones(1,size(data.cons.idn,2))/size(data.cons.idn,2);
param.fonction.idn         = 'zsinbad2temps';
param.cons.idn.energie     = 900e3*ones(16,1);
param.cons.idn.align       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
param.cons.idn.fraction1   = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
param.cons.idn.fraction2   = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
param.cons.idn.fraction3   = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
param.cons.idn.masse       = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
param.cons.idn.charge      = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
param.cons.idn.type        = 3;
param.cons.idn.debut       = 'Avant';
param.cons.idn.machine     = 'ITER';



% injection de fci
data.cons.fci     = (pfci0  .* medfilt1(interp1(scenario.t,scenario.pfci,temps,'linear'),5)) * ones(1,size(data.cons.fci,2));
data.cons.fce     = (pfce0  .* medfilt1(interp1(scenario.t,scenario.pfce,temps,'linear'),5)) * ones(1,size(data.cons.fce,2));

% autre donnees de dimensionnement
ploss             = sum(data.cons.idn,2) ./ 1e6 + sum(data.cons.fci,2) ./ 1e6 +  ...
                    sum(data.cons.fce,2) ./ 1e6 + data.cons.ip ./ 1e7;
% nbar en 10^20 m^-3
nbar              = nbar0  .* (data.cons.ip ./ 1e6) ./ pi ./ data.geo.a .^ 2   .* interp1(scenario.t,scenario.nbar,temps,'linear');
% li evolution 
x                 =  (ploss- min(ploss)) ./ max(ploss);
li                =  li0  + (1 / 0.7 - 1) .* (1 - x);

% calcul de meff
[ae,nion,nmin1,nmin2,nimp1,nimp2]=itercompo(nbar .* 1e20 ,data.cons.zeffm,c1,c2,rimp, ...
                                  param.compo.z(1),param.compo.z(2),param.compo.z(3), ...
                                  param.compo.z(4),param.compo.z(5));
meff = (nion .* param.compo.a(1) + nmin1 .* param.compo.a(2) + nmin2 .* param.compo.a(3)  + ...
        nimp1 .* param.compo.a(4) + nimp2 .* param.compo.a(5)) ./  ...
       (nion + nmin1 + nmin2 + nimp1 + nimp2);                          

% loi d'echelles pour le calcul des  profil initiaux 
R    = data.geo.r0;
a    = data.geo.a;
ep   = a ./ R;
K    = data.geo.e1;
b    = a .* K ;
ip   = data.cons.ip ./ 1e6;
Bt   = data.geo.b0;
zeff = data.cons.zeffm;
ne   = 10 .* nbar;

% rlw
wrlw  = 2.6e-2 .* ne .^ (3/4) .* zeff .^ (1/4)  .*  Bt .^ (1/2) .* ip  .^ (1/2) .* (R .* a .* b) .^ (11/12) + ...
        1.2e-2 .* ip .* (R .* a .* b) .^ (1/2) .* ploss .* zeff .^ (-1/2); %MJ
% ITERH-98P(y)        
wthh  = (36e-3  .* ip .^ 0.97 .* Bt .^ 0.08 .* ne .^ 0.41 .* ploss .^ -0.63 .* ...
         R .^ 1.93 .* K .^ 0.67 .* ep .^ 0.23 .* meff .^ 0.2) .* ploss;    % MJ     
% ITERH-96P(th)        
wthl  = (23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* ploss .^ -0.73 .* ...
         R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2) .* ploss;   % Mj       
% seuil mode H loi LH99(1)        
plosslh   = 2.84 .* nbar .^ 0.58 .* Bt .^ 0.82 .* R .* a .^ 0.81 ./ meff; % MW
         
% pression au piedestale 
ppied    = 2.8 .* ip .* Bt .^0.9 ./ a .* R .^ -0.125 .* 1e3; % Pa

% piquage de j :
piqj = (exp(li)-1.65)./0.89;
ind  = find(piqj <1);
piqj(ind) = 1 * ones(1,length(ind));
ind  = find(piqj >10);
piqj(ind) = 10 * ones(1,length(ind));


% 4 - remplissage des donnees dependants de rho et du temps
% le profil de densite (initial)
% la densite de bord  doit etre au plus de 1/3 de la densite cenrale
%  nea =  a nbar ^ 2 
%pp  = [-0.0635    7.2139 -153.2643];
pp  = [ -0.2867   27.7006 -622.8086];
nea = min(nbar .* 1e20 ./ 3,exp(polyval(pp,log(nbar .* 1e20))));
xx  = 1 - param.gene.x;
ne  = tanh(nfact .* xx);
ne  = ne ./ ne(1);
nbr = nbar .* 1e20 - nea;
nen = trapz(param.gene.x,ne,2);
ne  = (nbr ./ nen) * ne + nea * ve;
ni  = (ae * ve) .* ne;

% le profil de courant (initial)
jmoy  = (1 - (vt * param.gene.x) .^ 2 ) .^ (piqj * ve);
% normalisation approchee (sans importance pour la suite, renormlise par l'equilibre)
inorm = (2 .* pi .* data.geo.a .^ 2 .* data.geo.e1) .* trapz(param.gene.x,jmoy .* (ones(size(jmoy,1),1)*param.gene.x),2);
j0    = data.cons.ip ./ inorm;
jmoy  = (j0 * ve) .* jmoy;

% le profil de psi (initial)
peri   = 4 .* data.geo.a .* sqrt(data.geo.e1)  .* ellie(sqrt((data.geo.e1 - 1) ./ data.geo.e1));
grhor  = (0.95 +( 1.22 .* peri ./ (2 .* pi .* data.geo.a)  -1) * (param.gene.x .^ 2)) ./ ( data.geo.r0 * ve);
rhomax = peri ./2 ./ pi;
bpol0  = param.phys.mu0 .* j0 .* data.geo.a .* data.geo.e1 ./ (1+ piqj)  ./ 2;
bpol   = (bpol0 * ve) .* ((1 - (1 - (vt * param.gene.x) .^ 2) .^ ((piqj * ve) + 1)) ./ (vt * param.gene.x +eps));
bpol(:,1) = 0 .* vt;
psi    = (rhomax * ve) .* cumtrapz(param.gene.x,bpol ./ grhor,2);
psi1   = (17/400) .* data.gene.temps + 0.45  .* param.phys.mu0 .* (data.geo.r0 + data.geo.a) .* data.cons.ip;
psi1   = max(psi1)/2 - psi1;
psi    = psi(:,end) * ve - psi + psi1 * ve ; % verifier le sens ...


% la temperature du bord est inconnue on choisi 35 eV en mode L et 150 en mode H
% pour le mode L
xx  = [0 0.05 0.2 0.7 1];
tx  = [1 0.995 0.9 0.2 0];
tprof = spline(xx,tx,param.gene.x);
% te mode L
te1   = 35;
we1   = (3/2) .* param.phys.e .* volume .* trapz(param.gene.x,te1 .* ne,2);
we0   = wrlw .* 1e6 - we1;
te0   = we0 ./ ((3/2) .* param.phys.e .* volume .* trapz(param.gene.x,(vt * tprof) .* ne,2));
tel   = (te0 * ve) .* (vt * tprof) + te1;
% ti mode L
ti1   = 35;
wi1   = (3/2) .* param.phys.e .* volume .* trapz(param.gene.x,ti1 .* ni,2);
wi0   = (wthl - wrlw) .* 1e6 - wi1;
ti0   = wi0 ./ ((3/2) .* param.phys.e .* volume .* trapz(param.gene.x,(vt * tprof) .* ni,2));
til   = (ti0 * ve) .* (vt * tprof) + ti1;
% pour le mode H (Ti = Te)
% bord
t1    = 150;  % temperature de bord 150 eV
% piedestale
rpied = 0.15; % largeur du piedestale
ipied = length(param.gene.x) - round(length(param.gene.x) .* rpied ./ rhomax) ;
xpied = param.gene.x(ipied);
wpied = 0 .* vt;
tprof_p = 0 .* (vt *ve);
for k = 1:length(data.gene.temps);
   tpied = ppied(k) ./ (ni(k,ipied) + ne(k,ipied)) ./ param.phys.e;
	xx            = [0 xpied(k) 1];
	tx            = [tpied(k) tpied(k)  t1];
	tprof_p(k,:)  = interp1(xx,tx,param.gene.x,'linear');
	wpied(k)      = (3/2) .* param.phys.e .* volume(k) .* trapz(param.gene.x,tprof_p(k,:) .* (ni(k,:)+ne(k,:)),2);
end
% core
xx  = [0 0.05 0.2  1];
tx  = [1 0.995 0.9  0];
tprof = spline(xx,tx,param.gene.x);
w0  = wthh .* 1e6 - wpied;
t0  = w0 ./ ((3/2) .* param.phys.e .* volume .* trapz(param.gene.x,(vt * tprof) .* (ni + ne),2));
teh = (t0 * ve) .* (vt * tprof) + tprof_p;
tih = teh;

% basullement L-> H et H -> L sans hysteresis
flh  = medfilt1( ploss > plosslh,5);
flhv = flh * ve;
te   = tel .* (1- flhv) + teh .* flhv;
ti   = til .* (1- flhv) + tih .* flhv;
wdia = wthl .* (1 - flh) + wthh .* flh;


% regalge du split dtmax
tau               = (wdia./ploss);
tau               = mean(tau(isfinite(tau)));
param.split.dtmax = max(tau ./ 10,param.split.dtmax);

% les profils
data.prof.psi          = psi;
data.prof.dpsidt       = zdxdt(data.prof.psi,data.gene.temps);
data.cons.flux         = data.prof.psi(:,end);
data.cons.vloop        = - 2 .* pi .* data.prof.dpsidt(:,end);
data.prof.jmoy         = jmoy;
data.prof.ne           = ne;
ind = find(data.prof.ne<3e17);
if ~isempty(ind)
        data.prof.ne(ind)   = 3e17 .* ones(1,length(ind));
end
ind = find(ni<3e17);
if ~isempty(ind)
        ni(ind)   = 3e17 .* ones(1,length(ind));
end
ind = find(ni < (data.prof.ne./10));
if ~isempty(ind)
        ni(ind)   = data.prof.ne(ind) ./10;
end
ind = find(ni > data.prof.ne);
if ~isempty(ind)
        ni(ind)   = data.prof.ne(ind);
end
% ni initiale peut changer dans zineb
data.prof.ni           = ni;
data.prof.te           = te;
ind = find(data.prof.te<13.6);
if ~isempty(ind)
        data.prof.te(ind)   = 13.6 .* ones(1,length(ind));
end
data.prof.pe           = data.prof.te .* data.prof.ne .* param.phys.e;
data.prof.ti           = ti;
ind = find(data.prof.ti<13.6);
if ~isempty(ind)
        data.prof.ti(ind)   = 13.6 .* ones(1,length(ind));
end
ind = find(data.prof.ti < (data.prof.te./5));
if ~isempty(ind)
        data.prof.ti(ind)   = data.prof.te(ind) ./5;
end
% attention ni et ti sont calculer dans zineb et peuvent etre different en sortie
data.prof.pion         =  data.prof.ti .* ni .*param.phys.e;
data.prof.ptot         =  data.prof.pion + data.prof.pe;

% donnees experimentales pour le plot
data.exp.ne0 = data.prof.ne(:,1);
data.exp.nea = data.prof.ne(:,end);
data.exp.ni0 = data.prof.ni(:,1);
data.exp.nia = data.prof.ni(:,end);
data.exp.te0 = data.prof.te(:,1);
data.exp.tea = data.prof.te(:,end);
data.exp.ti0 = data.prof.ti(:,1);
data.exp.tia = data.prof.ti(:,end);
data.exp.ip       = data.cons.ip;
data.exp.vloop    = data.cons.vloop;
data.exp.betadia  = wdia.*1e6./ ((3/8).*param.phys.mu0.*data.geo.r0.*data.exp.ip .^2);
data.exp.li       = li;

% les consignes suites
data.cons.ne1 = data.prof.ne(:,end);
data.cons.te1 = data.prof.te(:,end);
data.cons.ti1 = data.prof.ti(:,end);


% fci reglage de la phase par defaut
data.cons.fci = data.cons.fci .* exp(i.*pi);
data.cons.hyb=zeros(size(data.cons.hyb)).*1e6;


% donnees pour le bord
% a faire plus tard
 
% les modes
v1   = ones(param.gene.nbt,1);
v0   = zeros(param.gene.nbt,1);
ind  = param.gene.kmin:1:param.gene.kmax;
von  = v0;
voff = v0;            % mis a zeros
von(ind) = v1(ind); 
vlit   = von;         % donnees en entree
vcalc  = 2 .* von;    % calculees
vcopie = 3 .* von;    % recopie du temps precedent
%
data.mode.impur      = vcalc;  

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

data.mode.equi       = vcalc;                 
data.mode.neo        = vcalc;                 
data.mode.fluce      = voff;                 
data.mode.flucion    = voff;                 
data.mode.rot        = voff; 

if option.psilim == 0
	data.mode.cons.psi        = v0;           
elseif option.psilim == 1
	data.mode.cons.psi        = v1;           
else 
	data.mode.cons.psi        = 2 .* v1;           
end

data.mode.cons.ne         = v0;            
data.mode.cons.pe         = v0;            
data.mode.cons.pion       = v0;           
data.mode.cons.fluce      = v0;            
data.mode.cons.flucion    = v0;            
data.mode.cons.rot        = v0;            
data.mode.consbord.ne     = v0;            
data.mode.consbord.te     = v0;            
data.mode.consbord.ti     = v0;            
data.mode.consbord.fluxge = v1;            
data.mode.consbord.fluxqe = v1;            
data.mode.consbord.fluxqi = v1;            
data.mode.mhd.dds    = voff;                 
data.mode.mhd.elm    = voff;                 
data.mode.mhd.limite = voff;                 

% (sources)
data.mode.fci        = vcalc;                 
data.mode.fce        = vcalc;                 
data.mode.hyb        = vcalc;                 
data.mode.idn        = voff;                 
data.mode.n0         = vcalc;                 
data.mode.bord       = vcalc;                 
data.mode.glacon     = voff;                 
data.mode.fus        = vcalc;                 
data.mode.ohm        = vcalc;                 
data.mode.qneo       = voff;                 
data.mode.qei        = vcalc;                 
data.mode.prad       = vcalc;                 
data.mode.brem       = vcalc;                 
data.mode.cyclo      = vcalc;                 

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
	
	data.mode.rotc       =  vcalc;                 
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

%mode de variables neo
data.mode.eta       = vcalc;    
data.mode.jboot     = vcalc;    

% autres
if option.plotonoff == 1
	data.mode.plot          = von;  
else
	data.mode.plot          = voff;  
end

data.mode.zeff          = vcalc;    
data.mode.ae            = vcalc;    
data.mode.asser         = voff;    


% modification du nom du fichier de sortie
param.gene.file = strcat(param.gene.origine,'_resultat');
param.plot.pause = 0;
[pp,fp,ep]=fileparts(param.gene.file);
param.gene.rapsauve =fullfile(pp,'rapsauve',fp);

% sauvegarde du fichier
if nargin <3
  try
    % compactage des donnees
    data=zreduit(param,data,'compact');
    % sauvegarde
    post=[];
    savedb(param.gene.origine,'param','data','post');
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



% calcul de la compositions du plasma 
% en entree : 
%    ne      = densite electronique
%    zeff    = profil de zeff
%    c1      = rapport de la densite du 1er minoriatire sur la densite de l'espece principale
%    c2      = rapport de la densite du 2ieme minoriatire sur la densite de l'espece principale
%    rimp    = rapport de la densite de la 2ieme impurete sur la densite de la premiere impurete
%    zion    = numero atomic (ou charge moyenne) de l'espece principale
%    zmin1   = numero atomic (ou charge moyenne) du 1er minoritaire
%    zmin2   = numero atomic (ou charge moyenne) du 2ieme minoritaire
%    zimp1   = numero atomic (ou charge moyenne) de la 1ere impurete
%    zimp2   = numero atomic (ou charge moyenne) de la 2ieme impurete
%    
% en sortie : ae,nion,nmin1,nmin2,nimp1,nimp2 
% 
% equations :
% ae*ne =nion+nmin1+nmin2+nimp1+nimp2
% ni =ae*ne
% nmin1 = c1 *nion
% nmin2 = c2*nion
% nimp2 = rimp*nimp1
% ne = nion*zion +nmin1 *zmin1 +nmin2 *zmin2 +nimp1*zimp1 +nimp2*zimp2
% ne *zeff = nion*zion^2 +nmin1 *zmin1^2 +nmin2 *zmin2^2 +nimp1*zimp1^2 +nimp2*zimp2^2
% 
function [ae,nion,nmin1,nmin2,nimp1,nimp2]=itercompo(ne,zeff,c1,c2,rimp,zion,zmin1,zmin2,zimp1,zimp2)

% variables de calcul
de =  - zimp1 .^ 2 .* zion - zimp1 .^ 2 .* c1 .* zmin1 - zimp1 .^ 2 .* c2 .* zmin2 - ...
        rimp .* zimp2 .^ 2 .* zion - rimp .* zimp2 .^ 2 .* c1 .* zmin1 -  ...
        rimp .* zimp2 .^ 2 .* c2 .* zmin2 + zion .^ 2 .* zimp1 +  ...
        zion .^ 2 .* rimp .* zimp2 + c1 .* zmin1 .^ 2 .* zimp1 + ...
        c1 .* zmin1 .^ 2 .* rimp .* zimp2 + c2 .* zmin2 .^ 2 .* zimp1 + ...
        c2 .* zmin2 .^ 2 .* rimp .* zimp2;
       
% especes principales            
nion = -ne .* (zimp1 .^ 2 + rimp .* zimp2 .^ 2 - zeff .* zimp1 - zeff .* rimp .* zimp2) ./ de;
nmin1 = c1.*nion;
nmin2 = c2.*nion;

% impuretees:
nimp1 = (zion .^ 2 + c1 .* zmin1 .^ 2 + c2 .* zmin2 .^ 2 - zion .* zeff - ...
         c1 .* zmin1 .* zeff - c2 .* zmin2 .* zeff ) .* ne ./ de;
nimp2 = rimp.*nimp1;

% rapport somme(ni)/ne :
ae = (- zimp1 .^ 2 - rimp .* zimp2 .^ 2 + zeff .* zimp1 + zeff .* rimp .* zimp2 - ...
        zimp1 .^ 2 .* c1 - c1 .* rimp .* zimp2 .^ 2 + c1 .* zeff .* zimp1 + ...
        c1 .* zeff .* rimp .* zimp2 - c2 .* zimp1 .^ 2 - c2 .* rimp .* zimp2 .^ 2 + ...
        c2 .* zeff .* zimp1 + c2 .* zeff .* rimp .* zimp2 + zion .^ 2 + ...
        c1 .* zmin1 .^ 2 + c2 .* zmin2 .^ 2 - zion .* zeff - c1 .* zmin1 .* zeff - ...
        c2 .* zmin2 .* zeff + zion .^ 2 .* rimp + c1 .* zmin1 .^ 2 .* rimp + ...
        c2 .* zmin2 .^ 2 .* rimp - zion .* zeff .* rimp - c1 .* zmin1 .* zeff .* rimp - ...
        c2 .* zmin2 .* zeff .* rimp) ./ de;

