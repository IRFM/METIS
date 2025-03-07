function [cr,data,param]=ztcvacces(numchoc,chemin,temps,option,racine)
%
% Couplage machine TCV
% preparation d'une simulation CRONOS pour TCV, en utilisant l'ecriture des donnees de ztcv
%
% Auteur: V. Basiuk
% version 3.0, 30 novembre 2004
%
% dernieres modifications
%
langue = getappdata(0,'langue_cronos');

if nargin <=1 
	valeur.nbrho       = 101;     % nombre de points radiaux [101]
	
	valeur.compo       = 2;       % composition du plasma 1 -> H/H, 2 -> D/H, 3 -> D/D  et 4 -> He/D
	
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
  
	defaut=valeur ;   
	
   if strcmp(langue,'francais')
	  info.nbrho       = 'nombre de points radiaux [101]';
	  info.compo       = 'composition du plasma 1 -> H/H, 2 -> D/H, 3 -> D/D  et 4 -> D/T';
	  info.psimode     = 'mode Psi : 1 -> interpretative, 2 -> predictive';
	  info.pemode      = 'mode Pe  : 1 -> interpretative, 2 -> predictive';
	  info.pionmode    = 'mode Pion: 1 -> interpretative, 2 -> predictive';
	  info.nelmode     = 'mode Ne  : 1 -> interpretative, 2 -> predictive';
	  info.psilim      = 'frontiere (for Psi) : 0 -> Ip, 1 -> Vloop, 2 -> Psi(bord)';
	  info.lambda      = 'terme lambda des equations de transport {0,3/.2,5/2}';
	  info.modecoef    = 'coefficients equation de transport : 0 -> tout, 1 -> convective + diag.';
	  info.self        = 'mode de fonctionnement completement  auto consistant si = 1, 0 -> coefficient +neoclassique auto consistant ';
	  info.fast        = '0 -> mode standart correct, 1-> pas de sources neoclassique, ni de coefficient neo self consistante, 2 -> mode optimiser pour la diffusion du courant';
	  info.source_bord = 'controle du clacul de recyclage  : 1 = recyclage recalculer a chaque sous pas de temps, 0 = comme les sources';
	  info.cn          = 'mode du solveur : 0 -> implicite, 0.5 -> Cranck-Nickolson';
	  info.plotonoff   = '1 -> graphique en mode interactif ';
	  info.verbose     = '1-> commentaires , 0 -> pas de commentaires';
	  info.rebuilt     = '1 -> reconstruction en fin de run du fichier CRONOS'; 
	  info.post        = '1 -> postprocessing en fin de run';
   
   else
	  info.nbrho       = 'number of radial points [101]';
	  info.compo       = 'plasma composition 1 -> H/H, 2 -> D/H, 3 -> D/D  et 4 -> D/T';
	  info.psimode     = 'Psi  mode: 1 -> interpretative, 2 -> predictive';
	  info.pemode      = 'Pe   mode : 1 -> interpretative, 2 -> predictive';
	  info.pionmode    = 'Pion mode : 1 -> interpretative, 2 -> predictive';
	  info.nelmode     = 'Ne   mode  : 1 -> interpretative, 2 -> predictive';
	  info.psilim      = 'limit boundary (for Psi) : 0 -> Ip, 1 -> Vloop, 2 -> Psi(bord)';
	  info.lambda      = 'lambda factor including in transport equation  {0,3/.2,5/2}';
	  info.modecoef    = 'transport equation coeffecient calculation  0 -> whole, 1 -> convective + diag.';
	  info.self        = 'operating mode :1 -> purely self consistant, 0 -> coefficient + neoclassic self consistant ';
	  info.fast        = '0 -> standart mode, 1 -> no neoclassical source  2 -> optimized mode for current diffusion';
	  info.source_bord = 'recycling mode  : 1 -> at each internal time step, 0 -> like sources';
	  info.cn          = 'mode solveur : 0 -> implicite, 0.5 -> Cranck-Nickolson';
	  info.plotonoff   = 'interractive mode with plot if = 1 ';
	  info.verbose     = '1-> comments , 0 -> no comments';
	  info.rebuilt     = '1 -> rebuilt the output data file at the end'; 
	  info.post        = '1 -> postprocessing calculation at the end of CRNONOS';
   end
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

% compatibilite
nbrho        = option.nbrho;
tdebut       = [];
tfin         = [];
mode_compo   = option.compo;

% cr par defaut
cr =0;

% gestion des entrees
if nargin <5
	cr = 1;
	disp('wrong number of argument !');
        return
end
if strcmp(langue,'anglais')
  if nargin <5
	  cr = 1;
	  disp('wrong number of argument !');
          return
  end
end


 
 
% valeurs par default
if isempty(temps)
	temps=linspace(0,2,301);
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
	nbrho=101;
end



% valeurs pour TCV

nbfci    = 0;
nbfce    = 9;
nbhyb    = 0;
nbidn    = 0;
nbglacon = 0;
nbg      = 5;
	

directory = racine;
filetemp  = [directory,'/tcvtemp.mat'];
fileprof  = [directory,'/tcvprof.mat'];
filepsi   = [directory,'/tcvpsi.mat'];
fileeq    = [directory,'/tcveq.mat'];
filefit   = [directory,'/tcvfit.mat'];
load(filetemp)
load(fileprof)
load(fileeq)

if isempty(tdebut)
	tdebut =min(tcvprof.tne);
end
if isempty(tfin)
	tfin = max(tcvprof.tne);
end
if isempty(nbrho)
	nbrho=length(tcvprof.rho);
end


% creation du nom du fichier
if isempty(chemin)
	chemin=strcat(getenv('HOME'),'/zineb/data');
end
fichier = strcat(chemin,'/zineb',int2str(fix(numchoc)),'tcv_',int2str(round(rem(numchoc,1)*10)), ...
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
param.from.machine = 'TCV';
param.from.shot.num = numchoc;
param.from.shot.date = [2002 11 12 22 23 23] ;
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
param.from.createur = 'ztcvacces';
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
	% plasma He/D
	param.compo.z = [2,1,1,4,18];
	param.compo.a = [4,2,1,9,40];
	c1            = 0.1;
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

data.geo.r0     = interp1(tcveq.t,tcveq.R0,temps);
data.geo.z0     = 0*temps;
data.geo.a      = interp1(tcveq.t,tcveq.a,temps);
data.geo.e1     = interp1(tcveq.t,tcveq.e,temps);
data.geo.b0     = interp1(tcveq.t,tcveq.B0,temps);
data.geo.trh1   = 0*temps;
data.geo.trb1   = 0*temps;
data.geo.ind1   = zeros(size(temps));
data.geo.mode   = 2.*ones(size(temps));
data.geo.R      = single(interp1(tcvtemp.trext,tcvtemp.rext,temps));
data.geo.Z      = single(interp1(tcvtemp.trext,tcvtemp.zext,temps));

% consignes
data.cons.ip      = interp1(tcvtemp.tip,tcvtemp.ip,temps,'linear');
data.cons.vloop   = interp1(tcvtemp.tvl,tcvtemp.vl,temps,'linear');
if isempty(tcvtemp.zeffm)
  data.cons.zeffm   = ones(size(temps))+1;
else
  data.cons.zeffm   = interp1(tcvtemp.tzeffm,tcvtemp.zeffm(:,1),temps,'linear');
end
data.cons.nhnd    = 0  .* temps; 
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

% les profils
groupe.plus =0;
%
% elimination des nan au centre et au bord
% methode douteuse
%
npsia=interp1(tcveq.tpsia,tcveq.psia,tcveq.t);
psi1 = interp1(tcveq.t,tcveq.psiloc+npsia*ones(1,size(tcveq.psiloc,2)),temps);
vol1=interp1(tcveq.t,tcveq.volloc,temps);
rhotcv=tcveq.x;
vol2 =interp1(tcveq.t,tcveq.vol,temps);
indnan = find(temps > max(tcveq.t));
if ~isempty(indnan)
    vol1(indnan,:) = 0;
    vol2(indnan,:) = 0;
end    
for kpsi=1:length(temps)
  psi(kpsi,:)=interp1(vol1(kpsi,:)/vol1(kpsi,1),psi1(kpsi,:),vol2(kpsi,:)/vol2(kpsi,end),'linear');
end
xpsi = ones(size(temps))*rhotcv;
for k=1:length(temps)
  indpsi               = find(~isnan(xpsi(k,:)));
  nxpsi                = xpsi(k,indpsi);
  npsi                 = psi(k,indpsi);
  [nxpsi,ixpsi]        = sort(nxpsi);
  nxpsi(1)            = 0;
  nxpsi(end)         = 1;
  data.prof.psi(k,:)  = interp1(nxpsi,npsi(ixpsi),param.gene.x);
end

Te                     = interp1(tcvprof.tTe,tcvprof.Te,temps);
Te                     = interp1(tcvprof.rhotcv',Te',param.gene.x')';
Te(Te<0)               = 13.6;
ne                     = interp1(tcvprof.tTe,tcvprof.ne,temps);
ne                     = interp1(tcvprof.rhotcv',ne',param.gene.x')';
data.prof.te           = Te;
data.prof.ne           = ne;

data.prof.pe           = data.prof.te .* data.prof.ne .* param.phys.e;
rti                    = ones(size(temps));
if ~isfield(tcvtemp,'Ti0')
    tcvtemp.Ti0=[];
end
if ~isempty(tcvtemp.Ti0)
  indrti               = find( temps > min(tcvtemp.tTi0) & temps < max(tcvtemp.tTi0));
  nrti                 = interp1(tcvtemp.tTi0,tcvtemp.Ti0,temps(indrti))./data.prof.te(indrti,1)*1000;
  rti(indrti)         = nrti;
end
data.prof.ti           = data.prof.te .* (rti*ones(size(param.gene.x)));



% injection de fce
if ~isempty(tcvtemp.Pec)
  data.cons.fce        = interp1(tcvtemp.tPec,tcvtemp.Pec,temps,'linear');
  indval = find(~isnan(tcvtemp.ecfreq));
  param.cons.fce.freq_ghz(indval) = tcvtemp.ecfreq(indval)/1e9;

  param.cons.fce.angle_pol(indval) = tcvtemp.ectheta(indval);
  param.cons.fce.angle_tor(indval) = tcvtemp.ecphi(indval);
  param.cons.fce.rant = 1.227*ones(1,nbfce);
  param.cons.fce.zant = [0 0.515 0.515 0 0.515 0.515 0 0 0];

  %
  % case using mean value (on time) of poloidal and toroidal angle
  %
  fceok                = find(~isnan(tcvtemp.ecphi));
  data.cons.fce(isnan(data.cons.fce))=0;
  angle_mul            = zeros(1,9);
  angle_mul(fceok)     = (360 + rem(tcvtemp.ecphi(fceok),360)) .* 1e-10  + fix((360 + rem(tcvtemp.ectheta(fceok),360)) .* 1e4) .* 1e-7;
  data.cons.fce        = data.cons.fce .* exp(sqrt(-1) .* (ones(size(temps))*angle_mul));
 else
  data.cons.fce        = zeros(1,size(data.cons.fce,2));  
end

data.cons.hyb        = zeros(length(temps),nbhyb);

% autre donnees de dimensionnement
pfce                   = sum(abs(data.cons.fce'))';
pfce(isnan(pfce))     = 0;
phyb                   = sum(abs(data.cons.hyb'))';
phyb(isnan(phyb))     = 0;
ploss                  = pfce ./ 1e6 + phyb ./ 1e6 +  ...
						  data.cons.ip ./ 1e6;
% nbar en 10^20 m^-3

nbar                    = data.cons.ip*10/1e6/pi./data.geo.a;
%
% li lu en entree
%
indok = find(~isnan(tcvtemp.li));
li                      = interp1(tcvtemp.tli(indok),tcvtemp.li(indok),temps);
% calcul de meff
zeff  = data.cons.zeffm * ones(1,size(data.prof.ne,2));
[ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(data.prof.ne,zeff,c1,c2,rimp, ...
                                  param.compo.z(1),param.compo.z(2),param.compo.z(3), ...
                                  param.compo.z(4),param.compo.z(5));
meff = (nion .* param.compo.a(1) + nmin1 .* param.compo.a(2) + nmin2 .* param.compo.a(3)  + ...
        nimp1 .* param.compo.a(4) + nimp2 .* param.compo.a(5)) ./  ...
       (nion + nmin1 + nmin2 + nimp1 + nimp2);                          

% loi d'echelles pour le calcul des  profil initiaux 
R    = data.geo.r0;
a    = data.geo.a;
volume = 2*pi*pi*a.*a.*R;
ep   = a ./ R;
K    = data.geo.e1;
b    = a .* K ;
ip   = data.cons.ip ./ 1e6;
Bt   = data.geo.b0;
zeff = data.cons.zeffm;
ne   = 10 .* nbar;
nsat = 0.06e20*ip.*R.*sqrt(param.compo.a(1))./K./a.^(2.5)+1e17;
nneo = min(nsat,nbar*1e20)/1e20;
qcyl = 2*pi*a.*a.*Bt./(ip*1e6)./R./param.phys.mu0;
% rlw
wrlw  = 2.6e-2 .* ne .^ (3/4) .* zeff .^ (1/4)  .*  Bt .^ (1/2) .* ip  .^ (1/2) .* (R .* a .* b) .^ (11/12) + ...
        1.2e-2 .* ip .* (R .* a .* b) .^ (1/2) .* ploss .* zeff .^ (-1/2); %MJ
weneo = ploss.*0.07.*nneo.*a.*R.^2.*qcyl;


% piquage de j :
piqj = (exp(li)-1.65)./0.89;
ind  = find(piqj <1);
piqj(ind) = 1 * ones(1,length(ind));
ind  = find(piqj >10);
piqj(ind) = 10 * ones(1,length(ind));


% le profil de courant (initial)
jmoy  = (1 - (vt * param.gene.x) .^ 2 ) .^ (piqj * ve);
% normalisation approchee (sans importance pour la suite, renormlise par l'equilibre)
inorm = (2 .* pi .* data.geo.a .^ 2 .* data.geo.e1) .* trapz(param.gene.x,jmoy .* (ones(size(jmoy,1),1)*param.gene.x),2);
j0    = data.cons.ip ./ inorm;
jmoy  = (j0 * ve) .* jmoy;

if sum(tcveq.wdia) == 0
  wdia = 4/3*weneo;
else
  tcveq.wdia(isnan(tcveq.wdia))=0;
  wdia = interp1(tcveq.twdia,tcveq.wdia,temps,'linear')/1e6;
end

% regalge du split dtmax
tau               = (wdia./ploss);
tau               = mean(tau(isfinite(tau)));
param.split.dtmax = max(tau ./ 10,param.split.dtmax);

% les profils
tcvtemp.vl(isnan(tcvtemp.vl))=0;
data.prof.dpsidt       = zdxdt(data.prof.psi,data.gene.temps);
data.cons.flux         = interp1(tcveq.tpsia,tcveq.psia,temps,'linear');
data.cons.vloop        = interp1(tcvtemp.tvl,tcvtemp.vl,temps,'linear');
data.prof.jmoy         = jmoy;
ind = find(data.prof.ne<3e17);
if ~isempty(ind)
        data.prof.ne(ind)   = 3e17 .* ones(1,length(ind));
end
% ni initiale peut changer dans zineb
data.prof.ni           = nion;

% attention ni et ti sont calculer dans zineb et peuvent etre different en sortie
data.prof.pion         =  data.prof.ti .* data.prof.ni .*param.phys.e;
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


% lh reglage de la phase par defaut

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
data.mode.consbord.fluxge = v0;            
data.mode.consbord.fluxqe = v0;            
data.mode.consbord.fluxqi = v0;            
data.mode.mhd.dds    = voff;                 
data.mode.mhd.elm    = voff;                 
data.mode.mhd.limite = voff;                 

% (sources)
data.mode.fci        = voff;                 
data.mode.fce        = vcalc;                 
data.mode.hyb        = voff;                 
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
[pp,fp,ep,vp]=fileparts(param.gene.file);
param.gene.rapsauve =fullfile(pp,'rapsauve',fp);

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
   % save(param.gene.origine,'param','data','post');
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
function [ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(ne,zeff,c1,c2,rimp,zion,zmin1,zmin2,zimp1,zimp2)

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



function o = zfiltsvd2m(e,n)

[u,s,v]  = svd(e,0);
s        = diag(s);
er       = sqrt(sum(s .^ 2));
indmax   = min(max(find((s./er) >= 0.05)),3);
s((indmax+1):end) = 0;

if nargin >1 
	for k =1:indmax
		[v1,v2,v3,m] = meanfilt1(u(:,k),n); 
		u(:,k)            = m';
	end
end
o         = u * diag(s) * v';

