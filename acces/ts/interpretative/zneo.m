% ZNEO2 calcule les grandeurs neoclassiques avec Nclass
%--------------------------------------------------------------
% fichier zneo2.m ->  zneo2
%
% hello fred
%
% fonction Matlab 5 :
%
% Cette fonction calcule les grandeurs neoclassiques avec Nclass.
% Elle calcule aussi qei et qneo.
% 
% syntaxe  :
%  
%      [neo,coef_diff,coef_conv] = zneo(parametre,prof,geo,impur,equi,compo,phys,lambda,nbrho,nbg,x,etatcharge,forces,ergrho,memoire);
%
% entree :
%
%      parametre       =    parametres propre a la fonction (param.cons.neo)
%      prof            =    structure des profils (datak.prof)
%      geo             =    structure de la geometrie (datak.geo)
%      impur           =    structure des impuretes (datak.impur)
%      equi            =    structure de l'equilibre mhd (datak.equi)
%      compo           =    parametres de la composition du plasma (param.compo)
%      phys            =    constantes physiques (param.phys)
%      lambda          =    facteur multiplicatif des termes convectif des equation de la chaleur (param.gene.lambda)
%      nbrho           =    nombre de points radiaux (param.gene.nbrho)
%      nbg             =    nombre d'especes (param.gene.nbg)
%      x               =    coordonnee normalisee (param.gene.x)
%      etatcharge      =    densite de chaque espece pour chaque etat de charge  [nbrho,nbetat], doit suivre compo et comporte tous les etats de charges, y compris le neutre
%      force1          =    1er moment des forces externes  pour l'ion principal (datak.source.totale.q ?)
%      force2          =    2ieme moment des forces externes  (?)
%      force3          =    3ime moment des forces externes  pour l'ion principal (datak.source.totale.wb ?)
%      ergrho          =    Er/gradient(rho) (V/unite(rho))
%      memoire         =    memoire.neo
% 
% sortie :
% 
%     neo              =  structure neo de datak
%     coef_diff        =  coeficient de diffusion de la matiere associe a chaque masse et chaque etat de charge (m^2/s) [nbrho,nbetat]
%     coef_conv        =  vitesse de convection de la matiere associe a chaque masse et chaque etat de charge (m^2/s) [nbrho,nbetat]
%
% parametres :
% 
%  parametre.bootmodel    = model pour le bootstrap (defaut = 0) :
%                                  0  -> Houlberg
%                                  1  -> Hirshmann
%                                  2  -> Kessel
%  parametre.nbeq         = nombre de moments pour le calcul de l'equilibre neoclassique (defaut = 3) {2,3}:
%  parametre.banane       = largeur de banane finie au centre  -> 1 (defaut = 1 ) {0,1}
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 03/02/2005.
% 
% 
% liste des modifications : 
%
%  * 16/04/2002 -> securite sur epar pour le premier temps
%  * 16/04/2002 -> correction anti warning
%  * 22/10/2002 -> nclass fortran 90 + correction bug vpinch
%  * 11/12/2002 -> interface en anglais
%  * 11/02/2003 -> blindage erneo (si nions = 0 -> nions = 1e-18)
%  * 14/03/2003 -> parametre noncreux applique �ne, ni et coef.ee coef.ii
%  * 01/07/2003 -> qei calcule sur les profil natif (meme si creux)
%  * 09/09/2003 -> cr�tion des nouvelles variables (sans les remplir)
%  * 15/09/2003 -> matrice des etats de charge en entree
%  * 22/09/2003 -> passage de forces et champ electrique radial
%  * 05/03/2004 -> ajout des options de rotation simplifiee qui existaient auparavant dans zneorotform (refondu dans zneo)
%  * 05/03/2004 -> + correction de la position d une parenthese dans le calcul de erneo !!
%  * 06/03/2004 -> correction du calcul de vtheta pour l espece principale
%  * 23/03/2004 -> mise en place de l'appel de nclass_1 pour la prise en compte de la rotation
%  * 10/06/2004 -> blindage variable sm1 (division par bpol, nulle au centre)
%  * 21/06/2004 -> calcul du parametre G
%  * 13/09/2004 -> effet des runaways sur la resistivite, calcul simple 
%  * 08/12/2004 -> protection densite nul 
%  * 10/01/2005 -> suppression du mode old
%  * 10/01/2005 -> la rotation provient toujours de cronos et n'est plus claculee en interne.
%  * 21/01/2005 -> remplacement de spline par tsplinet.
%  * 24/01/2005 -> supression du parametre coefonly (pas utilier)
%  * 24/01/2005 -> ajout de memoire
%  * 31/01/2005 -> correction du signe de coef_conv
%  * 03/02/2005 -> jboot plat au centre
%  * 20/11/2006 -> Replacement de nclass_1 par neocall: version par fichiers
%  * 05/12/2008 -> reprise de interp1 a la place de tsplinet en mode pchip
%  (VB)
%--------------------------------------------------------------------------
%
function [neo,coef_diff,coef_conv,neodir]=zneo(parametre,prof,geo,impur,equi,compo,phys,lambda,nbrho,nbg,x,etatcharge,force1,force2,force3,ergrho,memoire)

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	langue           = getappdata(0,'langue_cronos');


	% nombre d'especes
	try
		nbg = evalin('base','param.gene.nbg');
	catch
		nbg = 5;
	end
	
	valeur.bootmodel = 0;                % model pour le bootstrap (defaut = 0)
	valeur.nbeq      = 3;                % nombre de moments pour le calcul de l'equilibre neoclassique (defaut = 3)
	valeur.banane    = 1;                % largeur de banane finie au centre  -> 1 (defaut = 0 )
	valeur.noncreux  = 1;                % profil non creux  au centre  -> 1 (defaut = 0 )
	valeur.sqlim     = 0.5;              % valeur maximale du facteur de deformation d'orbite
	valeur.injnbi    = 'D';              % gaz injecte par NBI
	%valeur.coefonly  = 0;                % si 1-> calcul uniquement les coefficients de transport pour le code d'impuretes
	valeur.runaway   = 0;             % si 1-> prise en compte des runaways pour le calcul de la resistivte effective
	valeur.er_effect = 'No';                %  prise en compte des effet du champ eletrique radial
	valeur.version      = 'f90';             % version de nclass
	valeur.save      = 0;                % si 1-> sauve le contexte apres appel du fortran
	valeur.source_vtor ='internal';      % source de la donnee de vitesse toroidal
	valeur.impurity_vtor_index = min(nbg,4); % index de l'espece correspondant a la mesure de vitesse
	
	type.bootmodel   = 'integer';       % type entier
	type.nbeq        = 'integer';       % type entier
	type.banane      = 'integer';       % type entier
	type.noncreux    = 'integer';       % type entier
	type.sqlim       = 'float';          % type reel
	type.injnbi      = 'string';          % type 
	%type.coefonly    = 'integer';          % type 
	type.runaway        = 'integer';          % type
	type.er_effect        = 'string';          % type
	type.version        = 'string';          % type
	type.save        = 'integer';          % type
	type.source_vtor = 'string';
	type.impurity_vtor_index = 'integer';
	
	borne.bootmodel  = {0,1,2};          % valeurs possible
	borne.nbeq       = {2,3};            % valeurs possible 
	borne.banane     = {0,1};            % valeurs possible 
	borne.noncreux   = {0,1};            % valeurs possible 
	borne.sqlim      = [0,1];             % valeurs possible 
	borne.injnbi     = {'H','D','T','He3','He4','D'};            % valeurs possible
	%borne.coefonly   = {0,1};            % valeurs possible
	borne.runaway       = {0,1};            % valeurs possible 
	borne.er_effect     = {'Yes','No'};            % valeurs possible
	borne.version     = {'f77','f90'};            % valeurs possible
	borne.save       = {0,1};            % valeurs possible
	borne.source_vtor ={'internal','measurement', 'grad_ti'};
	for k=1:nbg
		borne.impurity_vtor_index{k} = k;
	end

	defaut.bootmodel = 0;                 % valeurs par defaut
	defaut.nbeq      = 3;                 % valeurs par defaut 
	defaut.banane    = 1;                 % valeurs par defaut 
	defaut.noncreux    = 1;                 % valeurs par defaut 
	defaut.sqlim    = 0.5;                
	defaut.injnbi       = 'D';
	%defaut.coefonly       = 0;
	defaut.runaway       = 0;                
	defaut.er_effect       = 'No';
	defaut.version       = 'f90';
	defaut.save       = 0;
	defaut.source_vtor ='internal';
	defaut.impurity_vtor_index = min(nbg,4); % index de l'espece correspondant a la mesure de vitesse


	  info.bootmodel = 'bootstrap model (default = 0) 0  -> Houlberg; 1  -> Hirshmann; 2  -> Kessel';
	  info.nbeq      = 'moment number used in the resolution (default = 3) {2,3}';
	  info.banane    = 'banana finite width effect  -> 1 (default = 1 ) {0,1}';
	  info.noncreux  = 'negative gradient for the temperature and density profile  -> 1 (default = 0 )';
	  info.sqlim    = 'maximale diffence from unit of the squeezing factor';
	  info.injnbi    = 'name of gaz injected by NBI,for which forces are computed in NBI module(source.q et qource.wb) ';
	  %info.coefonly    = 'for coefonly=1, zneo compute only diffusivity and pinch fot the impurities transport module';
	  info.runaway      = 'for runaway=1, runaway effect on effective resistivity (simple model)';
	  info.er_effect = 'if = Yes, take into account the radial electric field to compute neoclassical equilibrium';
	  info.version      = 'nclass version, f77 (old version), f90 (new version with some problem)';
	  info.save      = 'for save=1, function stack variables are save after mexfile call';
	  info.source_vtor ='source of toroidal rotation data used in radial electric field computation : 1 = CRONOS toroidal moment equation, 2 = external toroidal velocity measurement & 3 = toroidal moment proportionnal to gradient Ti (used param.gene.grad_ti parameter)';
	  info.impurity_vtor_index = 'index in cronos composition of ion species used for toroidal velocity measurement'; % index de l'espece correspondant a la mesure de vitesse

	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';     % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	sortie.description = 'Calcul des grandeurs neoclassiques avec Nclass';   % description (une ligne) de la fonction
        if strcmp(langue,'anglais')	
	     sortie.description = 'Neoclassical calcul using NCLASS';   % description (une ligne) de la fonction
        end
   	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	
	neo=sortie;
	return
end

% reservation des sorties
mt = NaN.*ones(1,nbrho);
mt1 = NaN.*ones(1,nbrho,nbg);
neo.eta           = mt;    % resistivite (ohm * m)
neo.jboot         = mt;    % courant de bootstrap (A/m^2)
neo.coef.ee       = mt;    % diffusivite thermique electronique neo  (m^2/s)
neo.coef.ve       = mt;    % vitesse de convection thermique electronique neo
neo.coef.ii       = mt;    % diffusivite thermique ionique neo
neo.coef.vi       = mt;    % vitesse de convection thermique ionique neo
neo.coef.nn       = mt;    % diffusivite electronique neo
neo.coef.vn       = mt;    % vitesse de convection electronique neo
neo.coef.rotv     = mt;    % 
neo.coef.rot      = mt;    % 
neo.vtor          = mt1;    % profil de vitesse de rotation toroidale (m/s)
neo.vtheta        = mt1;    % profil de vitesse de rotation poloidale (m/s)
neo.mach          = mt1;    % profil de nombre de mach 
neo.er            = mt;    % profil de champ electrique radial 
neo.flux.ne       = mt;    % profil de flux eletronique neoclassique
neo.flux.nion     = mt;    % profil de flux ionique (somme sur les especes) neoclasique
neo.flux.qe       = mt;    % profil de flux de chaleur electronique neoclassique
neo.flux.qion     = mt;    % profil de flux de chaleur ionique (somme sur les especes) neoclassique
neo.qei           = mt;    % echange ion electron (equipartition)
neo.qneo          = mt;    % echange ion electron (effet neoclassique)
neo.w_ion         = mt1;
neo.w_e           = mt;
neo.utheta_i      = mt1;
neo.utheta_e      = mt;
neo.gammae        = mt;
neo.g             = mt;    % parametre G, Hellander, p 258. , le calcul NClass est correct si G <<1 , en particulier pour le calcul du transport des impuretes
neo.fail          = NaN;






% variable optionnelle
if ~isfield(parametre,'runaway')
    parametre.runaway = 0;
end
if ~isfield(parametre,'er_effect')
    parametre.er_effect = 'No';
    warning('ZNEO : missing parameter er_effect set to "No"');
end
if nargin < 12
    etatcharge =  NaN .* ones(nbrho,sum(compo.z) + nbg);
elseif ~all(isfinite(etatcharge(:)))
    % etatcharge existe mais n'a pas ete initialise par le module impur
    etatcharge =  NaN .* ones(nbrho,sum(compo.z) + nbg);    
elseif size(etatcharge,2) ~= (sum(compo.z) + nbg)
    % etatcharge existe mais n'a pas ete initialise par le module impur
    etatcharge =  NaN .* ones(nbrho,sum(compo.z) + nbg);    
    warning('ZNEO : etatcharge dimension mismatch');
end
% mode impurete ou non
etatcharge   = double(etatcharge);
modeimpur = all(isfinite(etatcharge(:)));


% creation des indices de demultiplexage de la matrice etatcharge
indcompo   = []; % vecteur reliant l'indice dans compo a l'indice dans etatcharge
aetat      = []; % nombre de masse pour chaque etat dans etatcharge
zetat      = []; % nombre de charge pour chaque etat dans etatcharge (y compris 0)
indinion   = []; % vecteur depassage des sortie de Nclass a la matrice etatcharge
indoution  = []; % vecteur depassage de la matrice etatcharge aux entrees de NClass
indinv.in  = {}; % passage de impur-> neo 
indinv.out = {}; % passage de neo -> impur
count      = 1;
for k =1:nbg
   nbz = compo.z(k) +1;
   indcompo      = cat(2,indcompo,k .* ones(1,nbz));
   aetat         = cat(2,aetat,compo.a(k) .* ones(1,nbz));
   zetat         = cat(2,zetat,0:compo.z(k));
   indinv.in{k}  = count+(1:compo.z(k));
   indinv.out{k} = find((indcompo == k) & (zetat >0));
   count         = count + compo.z(k);
end
indoution    = find(zetat >0);
indinion     = 1 + (1:length(indoution));

% sortie pour le transport d'impuretes
coef_diff        = zeros(size(etatcharge));
coef_conv        = zeros(size(etatcharge));


if nargin < 15
    force1 = [];
    force2 = [];
    force3 = [];
end
if nargin < 16
    ergrho = [];
end
 
 
% normalisation des donnees
if isempty(force1)
    force1 = zeros(1,nbrho);
else
    ind = find(~isfinite(force1));
    if ~ isempty(ind)
       force1(ind) = 0;
    end
end        
if isempty(force2)
    force2 = zeros(1,nbrho);
else
    ind = find(~isfinite(force2));
    if ~ isempty(ind)
       force2(ind) = 0;
    end
end        
if isempty(force3)
    force3 = zeros(1,nbrho);
else
    ind = find(~isfinite(force3));
    if ~ isempty(ind)
       force3(ind) = 0;
    end
end        
if isempty(ergrho)
    ergrho   = zeros(1,nbrho);
else
    ind = find(~isfinite(ergrho));
    if ~ isempty(ind)
       ergrho(ind) = 0;
    end
end
grphi    = - ergrho;
psidrho  = prof.psid1 ./ equi.rhomax;
warning off
inter    = zcentre(grphi ./ psidrho);
warning on
gr2phi   = psidrho .* zder5pts(x,inter) ./ equi.rhomax;

% calcul du parametre de  'squeezing'  : limitation au cas du transport local
if isfield(parametre,'sqlim')
   sm1max  = parametre.sqlim;
else
   sm1max  = 0.5;
   disp('undefined sm1max parameter : you must update your data set');
end
sm1 = zeros(size(prof.bpol));
sm1(2:end) =  max(compo.a) .*  phys.ua .* gr2phi(2:end) ./ phys.e ./ prof.bpol(2:end) .^ 2;
if any(sm1 < -1)
   zverbose('=> too negative gradient of the radial electrical field : non local problem is not address by NClass\n');
end

gr2phimax  = abs(sm1max .* phys.e .* prof.bpol .^ 2 ./  max(compo.a) ./ phys.ua);
if any(abs(gr2phi) > gr2phimax)
    zverbose(' => too high gradient of the radial electrical field : value corrected \n');
    indm = find(gr2phi < -gr2phimax);
    gr2phi(indm) =  -gr2phimax(indm);
    indp = find(gr2phi > gr2phimax);
    gr2phi(indp) =  gr2phimax(indp);

    % dans ce cas on recalcul un champ electrique constistant
    warning off
    grphinew = cumtrapz(x,zcentre(gr2phi .* equi.rhomax ./ psidrho),2) .* psidrho;
    warning on
    % la valeur au centre est nulle
    grphi = grphinew - grphinew(1);
    neo.er = - grphi;
    %disp('breakpoint in zneo : correction grphi');
    %keyboard
else
    % dans cette version le champ electrique est une donnee d'entree 
    % elle peut etre modifier afin de rester dans les limites de NClass
    neo.er = ergrho;
end
 
% securite
etatmin   = 2e13 ./ max(1,((ones(size(etatcharge,1),1)*zetat) .^ 2) .* (ones(size(etatcharge,1),1)*aetat));
etatmin((ones(size(etatcharge,1),1)*zetat) == 0) = 0;
indmin = find(isfinite(etatcharge) & (etatcharge < etatmin));
etatcharge(indmin) = etatmin(indmin);


% precalcul
% calcul des gradients utiles 
nions   = squeeze(impur.impur)';  % size(nions) = [nbg,nbrho] !
%
% protection densite nul au bord, VB, 8 decembre 2004
%
% valeur minimale
nions(nions<=0) = 1e13;
% verification electroneutralite 
% il peut y avoir une difference entre ne et sum(nions*zions) car le calcul du transport des impuretes est 
% fait a des endroits differents dans zsolver1t et de maniere asyncrhone. La valeur de ne doit etre concervee.
%ind = find(sum(nions,2) == 30*length(x));
rapneutre  = sum(nions .* (compo.z' * ones(size(x))),1)./prof.ne;
% correction si necessaire
nions      = nions ./ (ones(size(nions,1),1)*rapneutre);

% gradient negatif partout
nionsd1 = pdederive(x,nions,0,2,2,1) ;
nid1    = pdederive(x,prof.ni,0,2,2,1) ;
ted1    = pdederive(x,prof.te,0,2,2,1);
tid1    = pdederive(x,prof.ti,0,2,2,1);
ted1v   = pdederive(x,prof.te,0,2,2,1);
tid1v   = pdederive(x,prof.ti,0,2,2,1);
dd1     = pdederive(x,equi.d,0,2,2,1);
ned1    = pdederive(x,prof.ne,0,2,2,1);
%
% protection q negatif au centre
%
if sum(equi.q<0) > 0
   disp('warning, equi.q < 0')  
   equi.q(equi.q<0) = 0.1;
end

if parametre.noncreux == 1
   ted1       = ted1    .* (ted1<0);
   tid1       = tid1    .* (tid1<0);
   newte      = cumtrapz(x,ted1)-trapz(x,ted1)+max(30,prof.te(end));
   nionsd1    = nionsd1 .* (nionsd1<0);
   for knion = 1:5
     newnions(knion,:)   = cumtrapz(x,nionsd1(knion,:))-trapz(x,-abs(nionsd1(knion,:)))+max(30,nions(knion,end));
   end
   %
   % preservation eletroneutralite
   newne      = sum(newnions' .* (ones(size(newnions',1),1)*compo.z),2)';
   newni      = sum(newnions',2)';
   nionsd1 = pdederive(x,newnions,0,2,2,1) ;
   ned1    = pdederive(x,newne,0,2,2,1) ;
   nid1    = pdederive(x,newni,0,2,2,1) ;
   ned1       = ned1    .* (ned1<0);
   nid1       = nid1    .* (nid1<0);
   nionsd1    = nionsd1 .* (nionsd1<0);
   newti      = cumtrapz(x,tid1)-trapz(x,tid1)+max(30,prof.ti(end));
else
   newte      = prof.te;
   newti      = prof.ti;
   newne      = prof.ne;
   newni      = prof.ni;
   newnions   = nions;
end
% rexterieur
% remplacement de spline par tsplinet 21/01/2005 puis remise de interp1
% avec pchip et extrap (5/12/2008 VB
rext    = interp1(double(equi.rhoRZ),max(squeeze(double(equi.R)),[],2)',equi.rhomax .* x,'pchip','extrap')';
%rext     = tsplinet(double(equi.rhoRZ),max(squeeze(double(equi.R))'),equi.rhomax .* x)';


if modeimpur == 0
   % information sur les I/O de la fonction bootHoulbergMat5:
   %     Adaptation Tore Supra
   %     Version du 14 sept 1998
   %     V. Basiuk, W.A. Houlberg  
   %    syntaxe:
   %    tab = bootHoulbergMat41(p1,p2,p3,p4,p5,p6,p7,p8,p9);
   %entrees du Mexfile
   %p1   charge electron et ions
   %p2   masse electron et ion
   %p3   champ electrique radial (V/m) et sa derivee 
   %p4   temperature electron et des ions en keV
   %p5   densit electron et des ions en m-3
   %p6   derive par rapport au rayon normalise de te et ti
   %p7   derive par rapport au flux normalise de ne et ni
   %p8   rayon externe de la surface magnetique (m); 
   %     petit rayon normalis ; 
   %     qpsi; 
   %     fraction de particule piegee; 
   %     Bphi moyenne sur les surfaces de flux;
   %     dVdx
   %     a (m) si derivees rayon normalis/1 si non normalis;
   %     grand rayon de la dernire surface magntique (m); 
   %     moyenne de B.^2 sur une surface de flux;
   %     moyenne de 1/B.^2 sur une surface de flux;
   %     Qpsi(1)
   %     C3
   %     <gradrho2/B2>
   %     elongation au centre
   %     <E.B>   (VT/m)
   %p9(1)   = 0 formule de Houlberg;  = 1 Hirshmann; = 2 kessel
   %p9(2)   = 2 ou 3  -> nb de moments pour le calcul de l'equilibre
   %p9(3)   = 0 -> pas d'effet de largeur fini 
   %p9(4)   = 1 -> 
   %p9(5)   = 1 -> 
   %--------------------------------------------------------------
   %    sorties
   %       tab(1)  =  <Jbs.B> (A*T/m**2)
   %       tab(2)  = <Jbs.B> driven by unit p'/p of electrons(A*T*rho/m**3)
   %       tab(3)  = <Jbs.B> driven by unit p'/p of espece ionique 1 (A*T*rho/m**3)
   %       tab(4)  = <Jbs.B> driven by unit p'/p of espece ionique 2 (A*T*rho/m**3)
   %       tab(5)  = <Jbs.B> driven by unit T'/T of electrons(A*T*rho/m**3)
   %       tab(6)  = <Jbs.B> driven by unit T'/T of espece ionique 1 (A*T*rho/m**3)
   %       tab(7)  = <Jbs.B> driven by unit T'/T of espece ionique 2 (A*T*rho/m**3)
   %       tab(8)  = <Jex.B> current response to fexiz (A*T/m**2)
   %       tab(9)  = Parallel electrical resistivity (Ohm*m) 
   %       tab(10) = flux de chaleur electronique (J*rho/m**3/s)
   %       tab(11) = flux de chaleur ionique (J*rho/m**3/s)
   %       tab(12) = flux d'electrons (rho/m**3/s)
   %       tab(13) = flux d'ions (rho/m**3/s)
   %       tab(14) = diffusion coefficient (diag comp) of s (rho**2/s) 
   %       tab(15) = somme des differentes vitesses electroniques(rho/s)
   %       tab(16) = total heat electronic cond coefficient (rho**2/s)
   %       tab(17) = heat electronic convection velocity (rho/s)
   %       tab(18) = heat ionic cond coefficient (rho**2/s)
   %       tab(19) = heat ionic convection velocity (rho/s)
   %-------------------------------------------------------------------- 

   % charge
   p1    = [-1;compo.z(:)];
   % nombre de masse
   p2    = [1;compo.a(:)];
   % champ electrique radial et sa derivees ou indices du gaz injecte
   p3=[0;0];
   % parametres de la fonction       
   p9  = [parametre.bootmodel;parametre.nbeq;parametre.banane;1;1;length(x)];

   % correction epar <> 0
   epar    = prof.epar;
   ind     = find(~isfinite(epar));
   if ~isempty(ind)
    epar(ind) = 0;
   end
   eparnz  = 0.01 ./ 2 ./ pi ./ geo.r0;
   epar    = eparnz .* (epar == 0) + epar .* (epar ~= 0);
   
   % temperatures 
   p4  = [newte' newti' * ones(1,nbg)]  .* 1e-3;
	
   % densites
   p5  = [newne'  newnions']; 
	
   % derivees des temperatures
   p6  = [ted1(1,:)' tid1(1,:)' * ones(1,nbg)]  .* 1e-3;
	
   % derivees des densites
   p7  = [ned1(1,:)' nionsd1']; 
        
   % structure complexe
   p8  = [rext';                           ...
       x;                               ...
       equi.q;                          ...
       equi.ftrap;                      ...
       prof.bphi ;                      ...
       equi.vpr  .* equi.rhomax;        ...
       equi.rhomax*ones(size(x(1,:)));  ...
       geo.r0*ones(size(x(1,:)));       ...
       equi.b2 ;                        ...
       equi.b2i ;                       ...
       equi.q(1)*ones(size(x(1,:)));    ...
       equi.r2i .* equi.vpr .* equi.rhomax; ...
       equi.grho2b2;                    ...
       equi.e(1)*ones(size(x(1,:)));    ...
       epar .* geo.b0];
   %
   % dimension interne profil -> 200
   %
   ndim          = 200;
   p4(ndim,6)    = 0;
   p5(ndim,6)    = 0;
   p6(ndim,6)    = 0;
   p7(ndim,6)    = 0;
   p8(15,ndim)   = 0;
   rpts                  = length(x);

   % nouvel appel pour nclass_1
   %
   % decopue de p8 en 15 composants
   %       
   p10  =p8(1,:)';
   p11  =p8(2,:)';
   p12  =p8(3,:)';
   p13  =p8(4,:)';
   p14  =p8(5,:)';
   p15  =p8(6,:)';
   p16  =p8(7,:)';
   p17  =p8(8,:)';
   p18  =p8(9,:)';
   p19  =p8(10,:)';
   p20  =p8(11,:)';
   p21  =p8(12,:)';
   p22  =p8(13,:)';
   p23  =p8(14,:)';
   p24  =p8(15,:)';
else
   % cas des impuretes disponible par etat de charge
   newne = sum(etatcharge .* (ones(size(etatcharge,1),1) * zetat),2)';
   ned1    = pdederive(x,newne,0,2,2,1) ;
   newni   = sum(etatcharge,2)';
   nid1    = pdederive(x,newni,0,2,2,1) ;
   for k =1:nbg
      newnions(k,:) = sum(etatcharge(:,indinv.out{k}),2)';
   end
   nionsd1 = pdederive(x,newnions,0,2,2,1) ;

   % remplissage des matrices de Nclass
   p1(1) = -1;
   p2(1) = 1;
   p1(indinion) = zetat(indoution);
   p2(indinion) = aetat(indoution);
   % indices du gaz injecte par NBI (rempli apres)
   p3=[0;0];
   % parametres de la fonction       
   p9  = [parametre.bootmodel;parametre.nbeq;parametre.banane;1;1;length(x)];

   % correction epar <> 0
   epar    = prof.epar;
   ind     = find(~isfinite(epar));
   if ~isempty(ind)
    epar(ind) = 0;
   end
   eparnz  = 0.01 ./ 2 ./ pi ./ geo.r0;
   epar    = eparnz .* (epar == 0) + epar .* (epar ~= 0);
   
   % temperatures 
   p4  = cat(2,newte',newti' * ones(1,length(p2) -1))  .* 1e-3;
	
   % densites
   p5  = cat(2,newne',etatcharge(:,indoution)); 
	
   % derivees des temperatures
   p6  = cat(2,ted1(1,:)',tid1(1,:)' * ones(1,length(p2)-1))  .* 1e-3;
	
   % derivees des densites
   etatcharged1= pdederive(x,etatcharge,0,2,1,1);
   p7  = cat(2,ned1(1,:)',etatcharged1(:,indoution)); 
        
   % structure complexe
   p8  = [rext';                           ...
       x;                               ...
       equi.q;                          ...
       equi.ftrap;                      ...
       prof.bphi ;                      ...
       equi.vpr  .* equi.rhomax;        ...
       equi.rhomax*ones(size(x(1,:)));  ...
       geo.r0*ones(size(x(1,:)));       ...
       equi.b2 ;                        ...
       equi.b2i ;                       ...
       equi.q(1)*ones(size(x(1,:)));    ...
       equi.r2i .* equi.vpr .* equi.rhomax; ...
       equi.grho2b2;                    ...
       equi.e(1)*ones(size(x(1,:)));    ...
       epar .* geo.b0];
   %
   % dimension interne profil -> 200
   %
   ndim          = 200;
   nprof         = size(p4,2);
   p4(ndim,nprof)    = 0;
   p5(ndim,nprof)    = 0;
   p6(ndim,nprof)    = 0;
   p7(ndim,nprof)    = 0;
   p8(15,ndim)   = 0;
   rpts                  = length(x);

   % nouvel appel pour nclass_1
   %
   % decopue de p8 en 15 composants
   %       
   p10  =p8(1,:)';
   p11  =p8(2,:)';
   p12  =p8(3,:)';
   p13  =p8(4,:)';
   p14  =p8(5,:)';
   p15  =p8(6,:)';
   p16  =p8(7,:)';
   p17  =p8(8,:)';
   p18  =p8(9,:)';
   p19  =p8(10,:)';
   p20  =p8(11,:)';
   p21  =p8(12,:)';
   p22  =p8(13,:)';
   p23  =p8(14,:)';
   p24  =p8(15,:)';

end
        
%
%  p3 : p3(1) = indice de charge de l'espece injectee dans p1
%       p3(2) = p2(indice de charge de l'espece injectee dans p1)
switch parametre.injnbi
case 'H'
     indgaz = min(find(p1 == 1 & p2 == 1));
case 'D'
     indgaz = min(find(p1 == 1 & p2 == 2));
case 'T'
     indgaz = min(find(p1 == 1 & p2 == 3));
case 'He3'
     indgaz = min(find(p1 == 2 & p2 == 3));
case 'He4'
     indgaz = min(find(p1 == 2 & p2 == 4));
otherwise
     indgaz = 2;
end
if isempty(indgaz)
  indgaz = 2;
end
p3(1) = indgaz;
p3(2) =  p2(indgaz);


switch parametre.er_effect
case 'Yes'
	%  nouvelles entrees  liee au forces externes
	%  p25 = pgrphipr  : radial electric field Phi' (V/rho)
	%  p26 = pgr2phipr : radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
	%
	p25 =cat(1,grphi',zeros(length(p24) - length(grphi),1));
	p26 =cat(1,gr2phi',zeros(length(p24) - length(gr2phi),1));
	
	%  p27 = fexizpr1  : first moment of external parallel force on species i,z (T*j/m**3)
	%  p28 = fexizpr2  : second moment of external parallel force on species i,z (T*j/m**3)
	%  p29 = fexizpr3  : third moment of external parallel force on species i,z (T*j/m**3)
	%  le moment d'ordre 1 est la quantite de mouvement (v), 2 -> la chaleur (v^3), 3- > en v^5
	%  ref : W. A. Houlberg, Phys. PLasmas 4 (9), september 1997, p 3230- (eq 4 et suivantes)
	%  l'usage normale : fexizpr1 et fexizpr2 doivent etre fournis pour IDN et fexizpr3 = 0;
	%  fexizpr1 = datak.source.totale.wb
	%  fexizpr2 = datak.source.totale.q
	% 
	p27 =cat(1,force1',zeros(length(p24) - length(force1),1));
	p28 =cat(1,force2',zeros(length(p24) - length(force2),1));
	p29 =zeros(size(p24));
	%
	% pour emploie de bootHoulberg6
	%
	%pold = [p25 p26 p27 p28 p29]';
	%tab         = bootHoulbergMat6(p1,p2,p3,p4,p5,p6,p7,p8,p9,pold);
otherwise
	p25 = zeros(size(p24));
	p26 = zeros(size(p24));
	p27 = zeros(size(p24));
	p28 = zeros(size(p24));
	p29 = zeros(size(p24));
end
% test :
%  p1(indinion)-zetat(indoution)=0;
%  p2(indinion)-aetat(indoution)=0;
%
%
% nouvel appel de neocall (normalement compatible)
% save contexte_houlberg
%    <<<<<<<<<<<=============================>>>>>>>>>>>>>
% attention : correction provisoire en attente de modification de NClass1 : p25 ->  phys.e ./ 1e3 .* p25
%    <<<<<<<<<<<=============================>>>>>>>>>>>>>
%disp('with Pfirsch-Schluter viscosity') 
%p9(5) = 1;
%p9(3) = 0;
if strcmp(parametre.version,'f90')
  [tab,memoire]    = neocall(p1,p2,p3,p4,p5,p6,p7,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,memoire);
  neodir = memoire.data.neodir;
end
if strcmp(parametre.version,'f77')
  pold        = [p25 p26 p27 p28 p29]';
%  clear mex
  tab         = bootHoulbergMat6(p1,p2,p3,p4,p5,p6,p7,p8,p9,pold);
  neodir = [];
end
%--------------------
% Variables de sortie
%____________________
%tabtr(1,ipr) = <Jbs.B> (A*T/m**2)
%------------------------------------------------------------------
% <Jbs.B> driven by unit p'/p of species s (A*T*rho/m**3)
% 2 -> electrons, 3-> espece ionique 1, 4 -> espece ionique 2
%------------------------------------------------------------------
%tabtr(2,ipr) = <Jbs.B> driven by unit p'/p of electrons (A*T*rho/m**3)
%tabtr(3,ipr) = <Jbs.B> driven by unit p'/p of espece ionique 1 (A*T*rho/m**3)
%tabtr(4,ipr) = <Jbs.B> driven by unit p'/p of espece ionique 2 (A*T*rho/m**3)
%------------------------------------------------------------------
% <Jbs.B> driven by unit T'/T of s (A*T*rho/m**3)
% 5 -> electrons, 6 -> espece ionique 1, 7 -> espece ionique 2
%------------------------------------------------------------------
%tabtr(5,ipr) = <Jbs.B> driven by unit T'/T of electrons (A*T*rho/m**3)
%tabtr(6,ipr) = <Jbs.B> driven by unit T'/T of espece ionique 1 (A*T*rho/m**3)
%tabtr(7,ipr) = <Jbs.B> driven by unit T'/T of espece ionique 2 (A*T*rho/m**3)
%-----------------------------------------------
%tabtr(8,ipr) = <Jex.B> current response to fexiz (A*T/m**2)
%tabtr(9,ipr) =  Parallel electrical resistivity (Ohm*m)
%------------------------------------------------------------
%       1  : p' and T' driven banana-plateau flux
%       2  : Pfirsch-Schluter flux
%       3  : classical flux
%       4  : <E.B> driven flux
%       5  : external force driven flux
%------------------------------------------------------------
%tabtr(10,ipr)  = conduction heat flux of electron (J*rho/m**3/s) (1+2+3+4+5)
%tabtr(12,ipr)  = particle flux of electron(rho/m**3/s)
%tabtr(11,ipr)  = conduction heat flux of all ions (rho/m**3/s) (1+2+3+4+5)   
%tabtr(13,ipr)  = particle flux ofall ions(rho/m**3/s)
%
%  vns      -> convection velocity (off diag comps-p', T') of s (rho/s)
%  vebs     -> <E.B> particle convection velocity of s (rho/s)
%  vnex     -> external force particle convection velocity [rho/s]
%  DP_SS    -> diffusion coefficient of s2 on p'/p of s1 [rho**2/s]
%  DT_SS    -> diffusion coefficient of s2 on T'/T of s1 [rho**2/s]
%  qfl(5,1) -> external force driven flux
%
%tabtr(14,ipr)   = diffusion coefficient (diag comp) of s (rho**2/s) 
%tabtr(15,ipr)   = somme des differentes vitesses electroniques (vns+vebs+vnex)
%tabtr(16,ipr)   = diffusion coefficient (DT_SS+DP_SS)
%tabtr(17,ipr)   = heat electronic convection velocity (rho/s)
%tabtr(18,ipr)   = sum [heat ionic cond coefficient of s2 on p'/p of s1 (rho**2/s)
%                  + heat ionic cond coefficient of s2 on T'/T of s1 (rho**2/s)]
%tabtr(19,ipr)   = sum heat ionic convection velocity (rho/s)
%
%
%
%tabtr(20,ipr)   = densite  
%tabtr(21,ipr)   = gradient Te
%tabtr(22,ipr)   = gradient Ti
%
% donnees electoniques
%
%tabtr(23,ipr)   = Ti
%tabtr(24,ipr)   = densite ionique
%tabtr(25,ipr)   = gradientTi
%-----------------------------------------------
%tabtr(26,ipr)   = <J_OH.B> Ohmic current [A*T/m**2]
%--------------------------------------|
% 1, p', T', Phi'                      |
% 2, <E.B>                             |
% 3, fex_iz                            |
%--------------------------------------|
%tabtr(27:(27+nspec),ipr)   = parallel flow of s from force  [T*m/s] (1+2+3) {UPAR_S}
%tabtr(ind:(ind+nspec),ipr) = poloidal flow of s from force  [m/s/T] [T*m/s] (1+2+3) {UTHETA_S}
%tabtr(ind:(ind+nspec),ipr) = diffusion coefficients (diag comp) [rho**2/s]
%abtr(ind:(ind+nspec),ipr) = Vitesse totale de convection par espece et etat de charge [rho/s]

% separation
nspec      = prod(size(p1));
indmem     = 26;
upar_out   = tab((indmem+1):(indmem+nspec),1:rpts);   % le premier vecteur de la matrice correspond aux electrons
indmem     = indmem+nspec;
utheta_out = tab((indmem+1):(indmem+nspec),1:rpts);   % le premier vecteur de la matrice correspond aux electrons
%
% protection point centrale
% mise a zero
% VB, 14 septembre 2007
%
utheta_out(:,1) = 0;
indmem     = indmem+nspec;
diff_out   = tab((indmem+1):(indmem+nspec),1:rpts);   % le premier vecteur de la matrice correspond aux electrons  + il manque la normalisation a ce niveau la
indmem     = indmem+nspec;
pinch_out  = tab((indmem+1):(indmem+nspec),1:rpts);   % le premier vecteur de la matrice correspond aux electrons + il manque la normalisation a ce niveau la

if parametre.save == 1
   save contexte_houlberg
end


% cas specifique du transport d'impuretes
if modeimpur == 1
   % decodage des coefficients de transport
   norme_diff  = (equi.rhomax .^ 2 ./ equi.grho2)' * ones(1,nprof);
   norme_conv  = (equi.rhomax ./ equi.grho2)' * ones(1,nprof);
   coef_diff   = NaN .* ones(size(etatcharge));
   coef_conv   = NaN .* ones(size(etatcharge));
   coef_diff(:,indoution) = diff_out(indinion,:)' .* norme_diff(:,indinion);
   coef_conv(:,indoution) = - pinch_out(indinion,:)' .* norme_conv(:,indinion);   % attention a la convention de signe
   coef_diff(1,indoution) = 0;
   coef_conv(1,indoution) = 0;
end


% les grandeur neoclassique (neo)
neo.eta               = tab(9,1:rpts);              % valide
if neo.eta(1) < 0
%
%
%
 zverbose('probleme nclass resistivite centrale negative')
 zverbose(' correction manuelle (neo.eta(1) = neo.eta(2)) + sauvegarde du contexte')
 save contexte_neo_probleme
 neo.eta(1)= neo.eta(2);

end
if sum(neo.eta) == 0
  save houlberg_fail
  neo.fail = 1;
else
      neo.fail = 0;
end
%
% limitation de la resisitivite, effet des runaways
% calcul simple
%
if parametre.runaway == 1
%
% ref : R. Martin-Solis et al, Contr. Fusion and Plasma Physics, 22C (1998), 794-797
%
  lnei          =  14.9 - 0.5.*log(prof.ne(1) ./ 1e20) + log(newte(1) ./ 1e3);
  consb         = 4 .* pi .* phys.epsi0.^2 .* phys.me .* phys.c.^2;
  Erun          = (phys.e.^3 * lnei ./ consb .* prof.ne(1)).*(mean(prof.zeff)+2);
  Epar          = max(-Erun,min(Erun,prof.epar));
  Eparneo       = prof.epar;
  Eparneo(Eparneo == 0) = eps;
  etastar       = max(0.1,abs(Epar./Eparneo)).*neo.eta;
  neo.eta       = etastar;
end
% le bootstrap complet
neo.jboot             = tab(1,1:rpts) ./ geo.b0;    % valide     

   
% remarque :
% les flux de Nclass sont definit comme :
%  <G.grad(rho)> = rhomax*g = rhomax * ( -D * dNe/dx +Ne * V)
% nous avons definit :
%   <G.grad(rho)> =<grad(rho)^2> * g = <grad(rho)^2> * (-D/rhomax * dNe/dx - Ne * V)
% le lien entre G et g :
%   G = - <grho2> g ?
% Le * rhomax est necessaire pour avoir la bonne unite (utilisation de x a la place de rho)
neo.flux.ne       =   tab(12,1:rpts) .* equi.rhomax ./ equi.grho2;
neo.flux.nion     =   tab(13,1:rpts) .* equi.rhomax ./ equi.grho2;   
neo.flux.qe       =   tab(10,1:rpts) .* equi.rhomax ./ equi.grho2;
neo.flux.qion     =   tab(11,1:rpts) .* equi.rhomax ./ equi.grho2;
neo.coef.nn       =   tab(14,1:rpts) .* equi.rhomax .^ 2 ./ equi.grho2;                  % valide 
neo.coef.vn       = - tab(15,1:rpts) .* equi.rhomax ./ equi.grho2;
if parametre.noncreux == 1
  neo.coef.ee       =   tab(16,1:rpts) .* equi.rhomax .^ 2 ./ equi.grho2 .* newne;    % valide 
  neo.coef.ii       =   tab(18,1:rpts) .* equi.rhomax .^ 2 ./ equi.grho2 .* newni;    % valide 
else
  neo.coef.ee       =   tab(16,1:rpts) .* equi.rhomax .^ 2 ./ equi.grho2 .* prof.ne;    % valide 
  neo.coef.ii       =   tab(18,1:rpts) .* equi.rhomax .^ 2 ./ equi.grho2 .* prof.ni;    % valide 
end

neo.coef.ve       = - tab(17,1:rpts) .* equi.rhomax ./ equi.grho2;
neo.coef.vi       = - tab(19,1:rpts) .* equi.rhomax ./ equi.grho2;
%
% correction au centre
%
if neo.coef.ee(1) == 0
  neo.coef.ee       = zcentre(neo.coef.ee);
end
if neo.coef.ve(1) == 0
  neo.coef.ve       = zcentre(neo.coef.ve);
end
if neo.coef.ii(1) == 0
  neo.coef.ii       = zcentre(neo.coef.ii);
end
if neo.coef.vi(1) == 0
  neo.coef.vi       = zcentre(neo.coef.vi);
end
if neo.coef.vn(1) == 0
  neo.coef.vn       = zcentre(neo.coef.vn);
end
if neo.coef.nn(1) == 0
  neo.coef.nn       = zcentre(neo.coef.nn);
end
%
% creation de coef.rotv (par defaut on prend comme le transport ionique)
%
neo.coef.rotv = neo.coef.vi;
neo.coef.rot  = neo.coef.ii ./ prof.ni;

%% clacul de de qei
%lnei          =  15.2 - 0.5.*log(prof.ne ./ 1e20) + log(prof.te ./ 1e3);
lnei          =  14.9 - 0.5.*log(prof.ne ./ 1e20) + log(newte ./ 1e3);
ind = find(~isfinite(lnei) | (lnei <10));
if ~isempty(ind)
	lnei(ind) = 10 .* ones(1,length(ind));
end

taue          = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.me))  .* ...
                ((phys.e .* prof.te) .^ (3/2) ./ prof.ne ./ lnei);
ve            = ones(1,size(prof.te,2));                
neo.qei       = 3 .* phys.me ./ phys.mp ./ taue .* sum(nions .* (compo.z'*ve) .^ 2 ./ (compo.a'*ve)) .*  ...
                (phys.e .* prof.te - phys.e .* prof.ti) ;
                
ind = find(~isfinite(neo.qei));
if ~isempty(ind)                
	neo.qei(ind) = zeros(1,length(ind));
end

% calcul de qneo
if lambda == (5/2)
	% ref: Tokamaks, Wesson 2ieme edition p 160  
	piond1            = (prof.ni.*tid1v + nid1.*prof.ti) .* phys.e;
	neo.qneo          =  prof.flux.gi ./ prof.ni .* piond1 ./ equi.rhomax .* equi.grho2 - ...
	                     0.172 .* phys.e .* prof.flux.gi .* tid1v ./ equi.rhomax;
	ind = find(~isfinite(neo.qneo));
	if ~isempty(ind)                
		neo.qneo(ind) = zeros(1,length(ind));
	end
else
        neo.qneo          = zeros(1,nbrho); 		
end

% prolongement par continuite des grandeurs au centre
neo.eta            = zcentre(neo.eta);
% neo.jboot          = zcentre(neo.jboot);
% le courant de bootstrap est nul au centre
neo.jboot(1)       = 0;
% les flux sont nuls au centre
neo.flux.ne(1)     = 0;
neo.flux.qe(1)     = 0;
neo.flux.nion(1)   = 0;
neo.flux.qion(1)   = 0;
neo.coef.ee        = zneoprotect(neo.coef.ee,1,'ee');
neo.coef.ve        = zneoprotect(neo.coef.ve,0,'ve');
neo.coef.ii        = zneoprotect(neo.coef.ii,1,'ii');
neo.coef.vi        = zneoprotect(neo.coef.vi,0,'vi');
neo.coef.nn        = zneoprotect(neo.coef.nn,1,'nn');
neo.coef.vn        = zneoprotect(neo.coef.vn,0,'vn');
neo.coef.rot       = zneoprotect(neo.coef.rot,1,'rot');
neo.coef.rotv      = zneoprotect(neo.coef.rotv,0,'rotv');


if  modeimpur == 1
      neo.utheta_e      = utheta_out(1,:);
      neo.utheta_i      = NaN .* ones(1,rpts,nbg);
      % calcul de la moyenne ponderee des densites et des charge pour chaque especes
      for k = 1:nbg
            neo.utheta_i(1,:,k) = (sum(utheta_out(indinv.in{k},:) .* (zetat(indinv.out{k})' * ones(1,rpts)) .* ...
                                  etatcharge(:,indinv.out{k})',1) ./ sum((zetat(indinv.out{k})' * ones(1,rpts)) .* ...
                                  etatcharge(:,indinv.out{k})',1))';
      end
else
      neo.utheta_i      = shiftdim(utheta_out(2:end,:)',-1);
      neo.utheta_e      = utheta_out(1,:);
end


% calcul de neo.g
% ref : Collisional Transport in Magnetized Plasmas, P. Helander and D.J. Sigmar, Cambridge Monographics on plasma physics, p 258
% le parametre G mesure la grandeur des force de friction // par rapport au gradient de pression //.
% Le parametre G permet de valider le resultat de NClass. Ce parametre doit etre plus petit que 1; en particulier pour ce qui concerne les 
% flux et le coefficient de transport. Le calcul du transport des impuretes n'est valide que si abs(G) << 1. Dans le cas contraire, comme avec une
% forte rotation, il n'y a plus d'invariance sur une surface de flux des densites.
psidrho  = prof.psid1 ./ equi.rhomax;
psidrho(1) = eps;
lnii     = 17.3 - 0.5.*log(prof.ne ./ 1e20) + 3/2 .* log(prof.ti ./ 1e3); % pour H, D et T uniquement
% pour l'espece principale
% Tokamaks, Wesson p 663
tauii     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* compo.a(1))  .* ...
	                   (phys.e .* prof.ti) .^ (3/2) ./ sum(nions,1) ./ lnii ./ compo.z(1) .^ 4 ;
meff      =  sum(nions .* (compo.a' * ones(1,size(nions,2))),1)./ sum(nions,1) .* phys.ua;
omi       =  phys.e .* sqrt(prof.zeff ./ prof.ae) .* sqrt(equi.b2) ./ meff;
bpol      =  prof.bpol;
bpol(1)   =  eps;
tid1      = pdederive(x,prof.ti,0,2,2,1); 
neo.g     = - prof.zeff ./ prof.ae .* equi.F .* sqrt(equi.b2) ./ omi ./ tauii ./ bpol .*  ...
              (prof.piond1 ./ equi.rhomax ./ psidrho  ./ max(eps,prof.pion) - 3 ./ 2  .* tid1 ./ equi.rhomax ./ psidrho ./ max(13.6,prof.ti));

neo.g = zcentre(neo.g);


% memorisation des forces
neo.force1 = force1;
neo.force2 = force2;
neo.force3 = force3;



%disp('keyboard at the end of zneo')
%keyboard

% protection contre les NaN Inf et les Chi < 0
function coef = zneoprotect(coef,plus,name)

if nargin <2
   plus = 0;
elseif isempty(plus)
   plus =0;
end   

coef = zcentre(coef);

ind  = find(~isfinite(coef));
if ~isempty(ind)
   if nargin >2
      zverbose(sprintf('Warning ZNEO : NaN ou Inf dans les coefficients (%s)\n',name));
   else
      zverbose('Warning ZNEO : NaN ou Inf dans les coefficients\n');
   end
   coef(ind) = zeros(1,length(ind));
end

if plus > 0
    % la valeur au centre n'est pas definie
    coef(1) = coef(2) .* (coef(1) <=0) + coef(1) .* (coef(1) >0);

    if all(coef <=0)
      if nargin >2
         zverbose(sprintf('Warning ZNEO : le coefficient de diffusion %s est nul ou negatif\n',name));
      else
         zverbose('Warning ZNEO : 0 dans les coefficients de diffusion Chi ou D \n');
      end
    end
    indnn     = find(coef >0);
    if ~isempty(indnn)
	     vmin = mean(coef(indnn))/length(coef);
    else
	     vmin = mean(abs(coef))/length(coef); 
    end
    ind  = find(coef < 0);
    if ~isempty(ind)
      if nargin >2
         zverbose(sprintf('Warning ZNEO : 0 dans les coefficients de diffusion Chi ou D (%s)\n',name));
      else
         zverbose('Warning ZNEO : 0 dans les coefficients de diffusion Chi ou D \n');
      end
         coef(ind) = vmin .* ones(1,length(ind));
         
    end
end



% fonction qui rend les profils monotones
function s =znoncreux(x,e)

dedx   = pdederive(x,e,0,2,2,1);
if any(dedx >= 0)
	ind            = min(find(dedx < 0));
	if ~isempty(ind)
	  dedx(1:ind)    = dedx(ind) .* ones(1,ind);
	end
end
s     = cumtrapz(x,dedx);
%
% blindage profil negatif
%
s     = abs(s - s(end) + e(end) .* (e(end)>0));


% calcul de la derive sur 5 points
function s = zder5pts(x,y)

c  = size(y,2);

x = x(:);
y = y(:);
delta    = mean(diff(x));
xx       = cat(1,-2.*delta,-delta,x,delta,2.*delta);
yy       = cat(1,y(1),y(1),y,y(end),y(end));
s        = (2 .* yy(5:end) + yy(4:(end-1)) - yy(2:(end-3)) - 2.*  yy(1:(end-4))) ./ 10 ./ delta;

if c > 1
   s = s';
end
