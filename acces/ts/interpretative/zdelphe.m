% ZDELPHE profil de depot hybride calcule par DELPHINE (RT+FP)
%--------------------------------------------------------------
% fichier zdelphe.m -> fichiers source dans repertoire /delphine/
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule le depot pour l'onde hybride.
% Elle doit etre utilisee lorsqu'on est en asservissement sur le flux au bord
% 
% syntaxe  :
%  
%      [sortie,memoire] = zdelphe(parametre,proto,cons,cons_spec,geo,equi,injection, ...
%                             prof,neo,impur,phy,composition,gene,memoire);
%
% entree :
%
%      parametre       =    parametre propre a la fonction (param.cons.hyb)
%      proto           =    prototype de la structure pour les sources, valeurs a zeros (proto = zsourceproto;)
%      cons            =    consigne de puissance par coupleur (data.cons.hyb)
%      cons_spec       =    spectre en puissance par coupleur (data.cons.hybspec)
%      geo             =    geometrie du plasma (data.geo)
%      equi            =    donnees de l'equilibre plasma (data.equi)
%      injection       =    consigne d'injection de gaz (data.cons.c)
%      prof            =    profils des donnees calculees par le code (data.prof) 
%      neo             =    donnees neoclassiques (data.neo)
%      impur           =    sous strcuture des impurtes (data.impur)
%      phy             =    constantes physiques (param.phys)
%      composition     =    composition du plasma: charge et masse des atomes (param.compo)
%      gene            =    parametres generaux (param.gene)
%      memoire         =    structure des dernieres valeurs calculees
% 
% sortie :
% 
%     sortie           =  structure de type source remplie par le module (sortie === proto)
%     memoire          =  datak.memoire.hyb (valeur de reference pour le dernier calcul complet, 
%                         pas utiliser dans cette fonction, reserve pour d'autres modules)
% 
%
% fonction ecrite par J-F Artaud et F. Imbeaux, poste 46-78 ou 63-26
% version 2.2, du 28/01/2004.
% 
% 
% liste des modifications : 
%
%   * 06/07/2001 -> ajout du champ electrique local
%   * 10/07/2001 -> ajout des reflexions speculaires en cas de "blocage du rayon"
%   * 10/07/2001 -> passage des profils ne et Te en numerique, au lieu de l'ancien parametrage
%   * 05/09/2001 -> ajout d'une securite si boucle infinie sur "blocage du rayon" (boucleRT)
%   * 26/09/2001 -> bug corrige dans la correction du courant LH si depassement de 99%Ip-Iboot
%   * 02/10/2001 -> ajout de la structure memoire en sortie
%   * 11/10/2001 -> reduction du nombre de parametres a regler par l'utilisateur (mis par defaut)
%   * 11/10/2001 -> prise en compte du profil de Zeff (au lieu de zeff moyen) -> modifs dans initRT, launchFP, launchFPohm
%   * 16/10/2001 -> ajout de la fonction d'interpolation
%   * 17/10/2001 -> ajout de la sauvegarde a la fin
%   * 19/10/2001 -> bug corrige : affectation du coefficient de diffusion (appel a par.D non effectue avant)
%   * 19/10/2001 -> reajustement du courant pour compenser la puissance perdue (discretisation du spectre ou rayons non completement absorbes)
%   * 24/10/2001 -> ajout des codes d'erreur (sortie.err)
%   * 26/10/2001 -> blindage contre les champs Epar trop eleves (dans launchFP)
%   * 06/11/2001 -> blindage contre l'absence de structure memoire en cas d'erreur
%   * 08/11/2001 -> ajout d'une astuce plus efficace pour eviter la boucle infinie sur le blocage du rayon (boucleRT)
%   * 20/11/2001 -> ligne 242 : definition de la surface magnetique ultime par ap*1.4, pour etre sur d'avoir une scrape-off bien fermee
%   * 20/11/2001 -> ligne 264 : mise en commentaire (ca servait a quoi, avant ???)
%   * 21/11/2001 -> blindage si l'hybride fait trop de courant : etendu a 95% de (Ip-Iboot), pour eviter les sauts de Vloop et les depassements (99% trop juste ?)
%   * 18/02/2002 -> blindage contre les champs Epar trop eleves mis egalement dans launchFPohm (incorrect avant !)
%   * 20/02/2002 -> blindage contre les indices paralleles < 1 (apres prepa) et enleve le blindage sauvage mis par VB dans dampfi lorsque ipp=[]
%   * 20/02/2002 -> remplacement de rac = 5 par une formule dependant de drmin et donnant des valeurs plus faibles (afin que la reflexion se fasse vraiment dans la scrape off) -> modif boucleRT
%   * 08/03/2002 -> le code d'erreur -2 est transforme en erreur = 4, afin que le calcul puisse se poursuivre (err>0)
%   * 14/03/2002 -> ajout du parametre effmult, facteur correctif de l'efficacite de generation de courant
%   * 20/03/2002 -> ajout des parametres epar et depip
%   * 20/03/2002 -> ajout de la possibilite de recalculer le champ E// et jLH de facon coherente (par.epar = 2)
%   * 03/04/2002 -> blindage contre Epar=NaN au debut d'une simulation
%   * 05/06/2002 -> correction de la syntaxe du isfield ligne 400 (cas par.epar=2) 
%   * 10/07/2002 -> prise en compte de l'orientation du champ toroidal (gene.signe.b0)
%   * 22/07/2002 -> securites sur le nombre de rayons pour le second spectre, et sur les consignes de phase
%   * 16/09/2002 -> verification qu'il n'y a pas de variation aberrante de Ilh d'un calcul a l'autre (sortie.err = -5)
%   * 18/09/2002 -> correction du module de reflexion speculaire : mise en place de spec2.m dans boucleRT
%   * 07/10/2002 -> ajout du parametre Dlim, qui permet de borner la valeur du Dql et donc de limiter les pb de convergence du FP
%   * 11/12/2002 -> interface en anglais
%   * 15/01/2003 -> le vecteur abscisse est directement pris sur gene.x*equi.rhomax, afin d'eviter les erreurs d'arrondi dans le calcul de Vcouche et Scouche (absopmhu)
%   * 17/01/2003 -> ajout de la gestion des configurations ITER 5 GHz et 3.7 GHz
%   * 05/02/2003 -> version pour nouveau Fokker Planck DKE, en version normalisation et addressage locaux
%   * 05/02/2003 -> version capable de calculer l'absorption par les particules alpha
%   * 11/02/2003 -> nouvelle formule pour calculer le rapport d'aspect (absopmhu)
%   * 10/06/2003 -> sauvegarde du contexte sur demand (par.save)
%   * 28/07/2003 -> ajout de la configuration FTU
%   * 28/01/2004 -> ajout de la configuration SST1
%   * 21/07/2004 -> lecture du traitement GHYBSPEC pour avoir les spectres complets pour TS
%   * 29/11/2005 -> prise en compte rigoureuse des orientations de Bt et Ip donnes par Cronos (marchait par chance avant)
% NB : les consignes d'entree et le calcul FP supposent par convention n// positif pour la direction co-courant
% mais le calcul des rayons doit prendre en compte le signe exact de k// et des champs magnetiques d'equilibre
% par consequent la relation de signe entre k// et n// est : -gene.signe.ip.*gene.signe.b0 (si Bt et Ip sont dans le meme sens, k//<0 pour faire du co-courant)
% NB : dans Delphine, le sens positif de Ufi (direction toroidale) est comme dans Cronos, i.e. positif dans le sens trigonometrique vu depuis le haut du tokamak
%   * 05/12/2007 -> Modification de wpi0 dans initRT
%   * 05/12/2007 -> Changement de la routine de preinterpolation de la grille d'equilibre : cartes 2D rhobar, Bpr, Bpz (R,Z), en utilisant griddata et une nouvelle prolongation de Bp a l'exterieur de la DSMF
%   * 06/12/2007 -> ajout configuration ASDEX (etude)
%   * 31/10/2008 -> nouveau blindage en sortie du trace de rayon, pour eviter les rayons qui penetrent la zone de non-accessibilite (donne Mfinal < 0)
%-------------------------------------------------------------
%
function [sortie,memoire] = zdelphe(par,proto,cons,cons_spec,geo,equi,injection, ...
                             prof,neo,impur,phy,compo,gene,memoire)
                             
% mode initialisation 
% fonction auto declarante                             
langue                  = getappdata(0,'langue_cronos');
if nargin <=1
	if nargin ==0
		nbhyb=1;
	else
		nbhyb=par;   
	end
        valeur.conf              = 'MIXTE';                    % configuration antennes / machine	
        valeur.nraysh            = 40;                      % nb de rayons 1ere antenne
        valeur.nraysv            = 40;                       % nb de rayons 2eme antenne
        valeur.seuilspec         = 0.02;                    % seuil pour la prise en compte d'un pic du spectre continu (si conf = 'TSspec')
        valeur.multispec         = 1;                    % facteur de sur-echantillonnage du spectre SWAN (si conf = 'TSspec')
        valeur.effmult           = 1;                       % facteur phenomenologique qui multiplie l'efficacite de generation de courant
        valeur.epar              = 1;                       % 0 : ne tient pas compte de E//; 1 : tient compte de E//
        valeur.depip             = 0;                       % 0 : reduction du courant LH s'il depasse la consigne de Ip; 1 : pas d'effet
		  valeur.D                 = 1;                       % diffusion phenomenologique du depot LH (m2/s)
        valeur.itermax           = 1000;                    % nb de points stockes pour chaque rayon
        valeur.zinterpolation    = 'zdelphe_interpolation'; % nom de la fonction rapide d'interpolation pour k -> k+1, mettre a vide si pas d'interpolation
        valeur.Dlim              = 10;                      % valeur max du Dql (pour eviter les problemes de convergence du FP)
        valeur.parabs            = 1;                       % ne pas toucher
	valeur.alpha             = 0;                       % si 1, calcule l'absorption sur les particules alpha

        valeur.save             = 'No';	                  % sauvegarde des fichiers contextes


        type.conf              = 'string';                    % configuration antennes / machine	
        type.nraysh            = 'integer';                      % nb de rayons 1ere antenne
        type.nraysv            = 'integer';                       % nb de rayons 2eme antenne

        type.seuilspec         = 'float';                    % seuil pour la prise en compte d'un pic du spectre continu (si conf = 'TSspec')
        type.multispec         = 'float';                    % facteur de sur-echantillonnage du spectre SWAN (si conf = 'TSspec')
        type.effmult           = 'float';                           % facteur phenomenologique qui multiplie l'efficacite de generation de courant
        type.epar              = 'integer';                       % 0 : ne tient pas compte de E//; 1 : tient compte de E//
        type.depip             = 'integer';                       % 0 : reduction du courant LH s'il depasse la consigne de Ip; 1 : pas d'effet
        type.D                 = 'float';                       % diffusion phenomenologique du depot LH (m2/s)
        type.itermax           = 'integer';                    % nb de points stockes pour chaque rayon
        type.zinterpolation    = 'string';                    % nom de la fonction rapide d'interpolation pour k -> k+1, mettre a vide si pas d'interpolation
	     type.Dlim              = 'float';                      % valeur max du Dql (pour eviter les problemes de convergence du FP)
	     type.parabs            = 'float';
		  type.alpha             = 'integer';                    % si 1, calcule l'absorption sur les particules alpha
	type.save               = 'string';
	
        borne.conf              = {'TSspec','TS','JET','CIEL','MIXTE','FTU','SST1','ITER5GHz','ITER3.7GHz','DEMO3.7GHz','DEMO5.0GHz','AUG3.7GHz'};  % configuration antennes / machine	
        borne.nraysh            = [1,400];                      % nb de rayons 1ere antenne
        borne.nraysv            = [0,400];                       % nb de rayons 2eme antenne
        borne.seuilspec         = [0,1];                    % seuil pour la prise en compte d'un pic du spectre continu (si conf = 'TSspec')
        borne.multispec         = [0.1,10];                    % facteur de sur-echantillonnage du spectre SWAN (si conf = 'TSspec')
        borne.effmult           = [0,100];                           % facteur phenomenologique qui multiplie l'efficacite de generation de courant
        borne.epar              = {0,1,2};                       % 0 : ne tient pas compte de E//; 1 : tient compte de E//
        borne.depip             = {0,1};                       % 0 : reduction du courant LH s'il depasse la consigne de Ip; 1 : pas d'effet
        borne.D                 = [0,100];                       % diffusion phenomenologique du depot LH (m2/s)
        borne.itermax           = [1,10000];                    % nb de points stockes pour chaque rayon
        borne.zinterpolation    =  {'zdelphe_interpolation',''}; % nom de la fonction rapide d'interpolation pour k -> k+1, mettre a vide si pas d'interpolation
  	     borne.Dlim              = [0,100];                   % valeur max du Dql (pour eviter les problemes de convergence du FP) 
		  borne.parabs            = [0,100];
		  borne.alpha             = {0,1};                    % si 1, calcule l'absorption sur les particules alpha

	  borne.save              = {'Yes','No'};


        defaut.conf              = 'MIXTE';                 % configuration antennes / machine	
        defaut.nraysh            = 40;                      % nb de rayons 1ere antenne
        defaut.nraysv            = 40;                       % nb de rayons 2eme antenne
        defaut.seuilspec         = 0.02;                    % seuil pour la prise en compte d'un pic du spectre continu (si conf = 'TSspec')
        defaut.multispec         = 1;                    % facteur de sur-echantillonnage du spectre SWAN (si conf = 'TSspec')
        defaut.effmult           = 1;                        % facteur phenomenologique qui multiplie l'efficacite de generation de courant
        defaut.epar              = 1;                       % 0 : ne tient pas compte de E//; 1 : tient compte de E//
        defaut.depip             = 0;                       % 0 : reduction du courant LH s'il depasse la consigne de Ip; 1 : pas d'effet
        defaut.D                 = 1;
        defaut.itermax           = 1000;                    % nb de points stockes pour chaque rayon
        defaut.zinterpolation    = 'zdelphe_interpolation'; % nom de la fonction rapide d'interpolation pour k -> k+1, mettre a vide si pas d'interpolation
	     defaut.Dlim              = 10;                      % valeur max du Dql (pour eviter les problemes de convergence du FP)
        defaut.parabs            = 3;
		  defaut.alpha             = 0;                    % si 1, calcule l'absorption sur les particules alpha

	defaut.save          = 'No';


        if strcmp(langue,'francais')
          info.conf              = 'configuration antennes / machine';	
          info.nraysh            = 'nb de rayons 1ere antenne';
          info.nraysv            = 'nb de rayons 2eme antenne';
          info.seuilspec         = 'seuil pour la prise en compte d''un pic du spectre continu (si conf = ''TSspec'')';
          info.multispec         = 'facteur de sur-echantillonnage du spectre SWAN (si conf = ''TSspec'') (>1 pour rajouter des rayons)';
          info.effmult           = 'facteur correctif de l''efficacite de generation de courant';
          info.epar              = '0 : ne tient pas compte de E//; 1 : tient compte du E// donne; 2 : calcul coherent de jLH et E//';
          info.depip             = '0 : reduction du courant LH s''il depasse la consigne de Ip; 1 : pas d''effet';
          info.D                 = 'diffusion phenomenologique du depot LH (m2/s)';
          info.itermax           = 'nb de points stockes pour chaque rayon';
          info.zinterpolation    = 'nom de la fonction rapide d''interpolation pour k -> k+1, mettre a vide si pas d''interpolation';
	       info.Dlim              = 'valeur max du Dql (pour eviter les problemes de convergence du FP)';
		    info.parabs            = 'ne pas toucher';
		    info.alpha             = 'si 1, calcule l''absorption sur les particules alpha';
	  info.save                = 'sauvegarde du contexte pour test';
		  end
        if strcmp(langue,'anglais')
          info.conf              = 'configuration antennas / tokamak';	
          info.nraysh            = 'number of rays for the first antenna';
          info.nraysv            = 'number of rays for the second antenna';
          info.seuilspec         = 'threshold for taking into account a peak of the continuous spectrum (if conf = ''TSspec'')';
          info.multispec         = 'resampling factor of SWAN power spectrum (if conf = ''TSspec'') (>1 for adding new rays)';
          info.effmult           = 'multiplier of the current drive efficiency';
          info.epar              = '0 : does not take into account E//; 1 : takes into account E// ; 2 : calulates jLH and E// self-consistently';
          info.depip             = '0 : prevents the LH driven current to be greater than Ip; 1 : no change';
          info.D                 = 'LH phenomenological current diffusion coefficient(m2/s)';
          info.itermax           = 'number of points stored for each ray';
          info.zinterpolation    = 'name of interpolation function for k -> k+1; set to empty for no interpolation';
	       info.Dlim              = 'Dql max (to avoid FP convergence problem)';
		    info.parabs            = 'do not change';
		    info.alpha             = 'if 1, calculates absorption by alpha particles';
	  info.save                = 'save context for test';
		  end
	interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'Calcul du depot de puissance Hybride avec DELPHINE (RT/FP)';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	
	
	return
end

% test sur la puissance
if (all(abs(cons)< 2e5) | all(abs(angle(cons))<=1))
    sortie = proto;
    return
end

global rhoexp neexp teexp methode Nbpoints

% abs(cons(1,2)) contient les puissances des antennes 1 et 2
% angle(cons(1,2)) contient les n//0 des antennes 1 et 2

% debut du calcul
% petits vecteurs utils
va = ones(gene.nbhyb,1);
ve = ones(1,gene.nbrho);
sortie = proto;
depot=[];

% parametres internes
drmin = 5e-3; % distance entre deux points stockes par le RT (m)
parabs = par.parabs; % parametre du calcul de l'absorption
Dlim = par.Dlim;  % valeur max du Dql autorisee
dR = 2e-2; % maillage de l'equilibre en R (m)
dZ = 2e-2; % mallage de l'equilibre en Z (m)

% on definit le flag d'erreur a 0 par defaut (pas d'erreur)
sortie.err = 0;

% definition de la grille (R,Z) pour le RT
conf=par.conf;
if strcmp(conf,'JET')
   rmin=1.6;
   rmax=4.1;
   zmin=-2;
   zmax=2.3;
   flh = 3.7e9; % frequence de l'onde (Hz)
elseif strcmp(conf,'ITER3.7GHz')
   rmin = 3.5;
   rmax = 9;
   zmin = -4;
   zmax = 4.5;
   flh = 3.7e9; % frequence de l'onde (Hz)
elseif strcmp(conf,'DEMO3.7GHz')
   rmin = geo.r0-geo.a*1.3;
   rmax = geo.r0+geo.a*1.3;
   zmin = -geo.a*2.5;
   zmax = geo.a*2.5;
   flh = 3.7e9; % frequence de l'onde (Hz)
elseif strcmp(conf,'DEMO5.0GHz')
   rmin = geo.r0-geo.a*1.3;
   rmax = geo.r0+geo.a*1.3;
   zmin = -geo.a*2.5;
   zmax = geo.a*2.5;
   flh = 5.0e9; % frequence de l'onde (Hz)
elseif strcmp(conf,'ITER5GHz')
   rmin = 3.5;
   rmax = 9;
   zmin = -4;
   zmax = 4.5;
   flh = 5e9; % frequence de l'onde (Hz)
elseif strcmp(conf,'FTU')
   rmin = 0.6;
   rmax = 1.3;
   zmin = -0.4;
   zmax = 0.4;
   flh = 8e9; % frequence de l'onde (Hz)
elseif strcmp(conf,'SST1')
   rmin = 0.7;
   rmax = 1.3;
   zmin = -0.4;
   zmax = 0.4;
   flh = 3.7e9; % frequence de l'onde (Hz)
elseif strcmp(conf,'AUG3.7GHz')
   rmin = 0.9;  % 1 avant pour les simus du 13679
   rmax = 2.4;
   zmin = -1.05;
   zmax = 1.05;
   flh = 3.7e9; % frequence de l'onde (Hz)
else
   % pour TS :
   rmin = 1.35;
   rmax = 3.3;
   zmin = -0.95;
   zmax = 0.95;
   flh = 3.7e9; % frequence de l'onde (Hz)
end

% Construction de la grille 2D de DELPHINE
R   =  rmin:dR:rmax;
Z   =  zmin:dZ:zmax;

RR  = R'*ones(1,size(Z,2));
ZZ  = ones(size(RR,1),1)*Z;

% Construction des matrices 2D (rho,theta) a partir de l'equilibre CRONOS
RRin = squeeze(double(equi.R));
RRin = RRin(2:end,1:(end-1));
ZZin = squeeze(double(equi.Z));
ZZin = ZZin(2:end,1:(end-1));
BRin = squeeze(double(equi.BR));
BRin = BRin(2:end,1:(end-1));
BZin = squeeze(double(equi.BZ));
BZin = BZin(2:end,1:(end-1));
% Prolongation au dela de la separatrice de l'equilibre, a + 10 % pour decrire la SOL, a + 20 % pour avoir des interpolation a la limite exterieure de la SOL (normalement, aucun rayon ne sort de la SOL, a + 10 %)
RRplus=(RRin(end,:)-geo.r0)*1.1+geo.r0;
ZZplus=(ZZin(end,:)-geo.z0)*1.1+geo.z0;
RRplus2=(RRin(end,:)-geo.r0)*1.2+geo.r0;
ZZplus2=(ZZin(end,:)-geo.z0)*1.2+geo.z0;
RRes = [RRin' RRplus' RRplus2']';
ZZes = [ZZin' ZZplus' ZZplus2']';
% Prolongation du champ poloidal : Theoreme Ampere applique en geometrie circulaire --> 1/rho_geo  & 1/R de la geometrie torique
% Prolongation depuis la separatrice (BRin(end,:),RRin(end,:)
BRplus = BRin(end,:)/1.1 .* RRin(end,:) ./ RRplus;  % Approximation circulaire pour faire ce developpement, pas tres grave ...
BZplus = BZin(end,:)/1.1 .* RRin(end,:) ./ RRplus;  % Approximation circulaire pour faire ce developpement, pas tres grave ...
BRplus2 = BRin(end,:)/1.2 .* RRin(end,:) ./ RRplus2;  % Approximation circulaire pour faire ce developpement, pas tres grave ...
BZplus2 = BZin(end,:)/1.2 .* RRin(end,:) ./ RRplus2;  % Approximation circulaire pour faire ce developpement, pas tres grave ...
BRes = [BRin' BRplus' BRplus2']';
BZes = [BZin' BZplus' BZplus2']';
% coordonnee de flux toroidal pour le calcul de rhobar(R,Z)
rho = double(equi.rhoRZ(2:end)') * ones(1,size(RRin,2));
ap  = equi.rhomax; % coordonnee de flux toroidal de la derniere surface magnetique fermee (m)
av = ap * 1.1; % limite exterieure de la SOL placee a 10% au dela de la derniere surface magnetique fermee; aucun rayon au dela de cette zone
rhoes=[rho' ones(size(rho,2),1)*av ones(size(rho,2),1)*ap*1.4]'; %ap*1.4 cette derniere surface, afin d'etre sur de bien clore la scrape-off (ap*1.1)

% Interpolation de cet equilibre prolonge (+ SOL), en grille (rho,theta) sur la grille rectangulaire (R,Z) de DELPHINE  NEW 05/12/2007
rhobar = griddata(RRes,ZZes,rhoes,RR,ZZ,'cubic'); 
Bpr = griddata(RRes,ZZes,BRes,RR,ZZ,'cubic'); 
Bpz = griddata(RRes,ZZes,BZes,RR,ZZ,'cubic'); 

% calcul des derivees partielles de rhobar
[drhobardz drhobardr]=gradient(rhobar,Z,R);

% prise en compte du signe de Ip donne par Cronos
Bpr = - Bpr .* gene.signe.ip;        % car par convention, Helena donne Bpr et Bpz tels que Ip < 0
Bpz = - Bpz .* gene.signe.ip;

% position du centre magnetique (m)
Rm=double(equi.R(1,1,1));
Zm=double(equi.Z(1,1,1));

% te et ne sur la coordonnee de flux toroidal
methode = 1; % prend directement les profils de CRONOS, sans faire un nouveau fit
rhoexp = gene.x*equi.rhomax;
neexp = prof.ne;
teexp = prof.te/1e3;
Nbpoints = length(rhoexp);


% installation du ripple
if strcmp(conf,'JET')|strcmp(conf,'ITER3.7GHz')|strcmp(conf,'ITER5GHz')|strcmp(conf,'FTU')|strcmp(conf,'SST1')|strcmp(conf,'AUG3.7GHz')
   ripple = 0; % ripple desactive sur les autres machines
   N = 1; % nb de bobines toroidales (bidon)
else
   ripple = 1; % ripple pris en compte sur TS
   N = 18; % nb de bobines toroidales (sur TS)
end
[delta,ddeltadr,ddeltadz,a1,a2,a3,b1,b2]=mripple(ripple,R,Z);


% vecteur param (donnees de l'equilibre)
Bo      = geo.b0.*gene.signe.b0; % champ toroidal en Ro (T)
Ro      = geo.r0; % position en R du centre g�m�rique du plasma (m)
ti0     = prof.ti(1)./1e3; % Ti centrale (keV)
Zion    = compo.z(1); % charge de l'ion majoritaire
Aion    = compo.a(1); % masse de l'ion majoritaire (uma)


% parametres de la simulation
itermax = par.itermax;      

% parametres antennes / spectres
nraysh      =  par.nraysh;     
npar0h      =  angle(cons(1));
if (npar0h < 0) & strcmp(conf,'AUG3.7GHz')
   disp('Special case for AUG study : high positive n// for first spectrum : 2*pi added to npar0h')
   npar0h = npar0h + 2*3.14159;
end
if length(cons)>1
   npar0v      =  angle(cons(2));
   if (npar0v < 0) & strcmp(conf,'AUG3.7GHz')
      disp('Special case for AUG study : high positive n// for second spectrum : 2*pi added to npar0v')
      npar0v = npar0v + 2*3.14159;
   end
   if strcmp(conf,'ITER5GHz')
      disp('Special case for ITER : npar0v forced to -6')
      npar0v = -6;
   end
   if abs(npar0h-npar0v)<0.05 & ~strcmp(conf,'TSspec')
      if strcmp(langue,'anglais')
        disp('the two N// are close : the two launchers are coupled !!')
      else
        disp('LES DEUX N// SONT PROCHES : LES DEUX ANTENNES SONT GROUPEES !!')
      end
	   cons(1)= (abs(cons(1))+abs(cons(2)))*exp(i*npar0h);
	   cons(2)= 0;
	   nraysv = 0;
   end
   if abs(cons(2))<10000
      nraysv = 0; % s'il n'y a pas de puissance dans le deuxieme spectre, alors on ne calcule pas les rayons
      if strcmp(langue,'anglais')
      disp('NRAYSV = 0, as there is no injected power for the second launcher') 
      else
      disp('NRAYSV est mis a 0, car aucune puissance n''a ete entree pour le deuxieme spectre') 
      end
   else
      nraysv = par.nraysv;
	   if nraysv == 0
	      nraysv = 40;
      if strcmp(langue,'anglais')
		   disp('NRAYS is put to 40 as injected power is detected for the second launcher')
      else
		   disp('NRAYS force a 40 car presence de puissance pour le deuxieme spectre')
      end
	   end
	end
else
   nraysv=0;
end
nrays       =  nraysh + nraysv;
listeray    =  1:nrays;

if strcmp(conf,'TS')
%   poscoupleur = 3.1; % grand rayon de depart approximatif des rayons
   poscoupleur = Ro + ap;
   nparwdthh   =  0.4;
   nparwdthv   =  0.4;
   direch = 0.7; % directivite antenne
   direcv = 0.7;
elseif strcmp(conf,'MIXTE')
%   poscoupleur = 3.1; % grand rayon de depart approximatif des rayons
   poscoupleur = Ro + ap;
	poscoupleur = 3.05;
   nparwdthh   =  0.4;
   nparwdthv   =  0.3;
   direch = 0.7; % directivite antenne
   direcv = 0.6;
elseif strcmp(conf,'CIEL')
%   poscoupleur = 3.1; % grand rayon de depart approximatif des rayons (m)
   poscoupleur = Ro + ap;
   nparwdthh   =  0.3;
   nparwdthv   =  0.3;
   direch = 0.6; % directivite antenne
   direcv = 0.6;
elseif strcmp(conf,'TSspec')     % spectre continu pour TS
%   poscoupleur = 3.1; % grand rayon de depart approximatif des rayons (m)
   poscoupleur = Ro + ap;
   % nouveau 28/07/04 : on sur-echantillonne le spectre pour mettre plus de rayons
   npar1new = linspace(cons_spec.npar(1,1,1),cons_spec.npar(1,end,1),par.multispec*length(cons_spec.npar(1,:,1)));
   pow1new = interp1(cons_spec.npar(1,:,1),cons_spec.pow(1,:,1),npar1new);
   npar2new = linspace(cons_spec.npar(1,1,2),cons_spec.npar(1,end,2),par.multispec*length(cons_spec.npar(1,:,2)));
   pow2new = interp1(cons_spec.npar(1,:,2),cons_spec.pow(1,:,2),npar2new);

   % on ne garde que les rayons au dessus du seuil en puissance par.seuilspec
   ipar1 = find(pow1new./max(pow1new) > par.seuilspec);
   listenpar1 = npar1new(ipar1);
   listepow1  = pow1new(ipar1);
   ipar2 = find(pow2new./max(pow2new) > par.seuilspec);
   listenpar2 = npar2new(ipar2);
   listepow2  = pow2new(ipar2);
   % renormalisation a la puissance de chaque coupleur
   listepow1 = listepow1 .* abs(cons(:,1)) ./ sum(listepow1);
   listepow2 = listepow2 .* abs(cons(:,2)) ./ sum(listepow2);
   direch = 1; % directivite antenne
   direcv = 1;
   nraysh = length(ipar1);
   nraysv = length(ipar2);
   if abs(cons(2))<10000
      nraysv = 0; % s'il n'y a pas de puissance dans le deuxieme spectre, alors on ne calcule pas les rayons
      if strcmp(langue,'anglais')
         disp('NRAYSV = 0, as there is no injected power for the second launcher') 
      else
         disp('NRAYSV est mis a 0, car aucune puissance n''a ete entree pour le deuxieme spectre') 
      end
   end
   nrays = nraysh + nraysv;
   listeray =  1:nrays;
elseif strcmp(conf,'JET')
   poscoupleur = 0; % decalage arbitraire de l'antenne par rapport a sa position de reference (grand rayon, m)
   nparwdthh   =  0.46;
   nparwdthv   =  0.;
   direch = 0.7; % directivite antenne
   direcv = 0.7;
elseif strcmp(conf,'FTU')
   poscoupleur = 0.; % decalage arbitraire de l'antenne par rapport a sa position de reference (grand rayon, m)
   nparwdthh   =  0.5;
   nparwdthv   =  0.5;
   direch = 0.85; % directivite antenne
   direcv = 0.85;
elseif strcmp(conf,'SST1')
   poscoupleur = 0.; % decalage arbitraire de l'antenne par rapport a sa position de reference (grand rayon, m)
   nparwdthh   =  0.74;
   nparwdthv   =  0.74;
   direch = 0.75; % directivite antenne
   direcv = 0.75;
elseif strcmp(conf,'ITER3.7GHz')|strcmp(conf,'ITER5GHz')
   poscoupleur = 8.2; % grand rayon de depart approximatif des rayons (m)
   nparwdthh   =  0.4; % a preciser quand on en saura plus !
   nparwdthv   =  0.4;
   direch = 1; % directivite antenne   % 0.7 si un seul lobe
   direcv = 1; 
elseif strcmp(conf,'DEMO3.7GHz')
   poscoupleur = geo.r0+geo.a+0.2; % grand rayon de depart approximatif des rayons (m)
   nparwdthh   =  0.4; % a preciser quand on en saura plus !
   nparwdthv   =  0.4;
   direch = 0.7; % directivite antenne
   direcv = 0.7; % a preciser quand on en saura plus !
elseif strcmp(conf,'DEMO5.0GHz')
   poscoupleur = geo.r0+geo.a+0.2; % grand rayon de depart approximatif des rayons (m)
   nparwdthh   =  0.4; % a preciser quand on en saura plus !
   nparwdthv   =  0.4;
   direch = 0.7; % directivite antenne
   direcv = 0.7; % a preciser quand on en saura plus !
elseif strcmp(conf,'AUG3.7GHz')
   poscoupleur = 0.05; % avancee du coupleur par rapport a la position standard (coupleur avance poscoupleur > 0) (m)
   nparwdthh   =  0.3;    % largeur de spectre as C3 ! 
   nparwdthv   =  0.3;
   direch = 0.7; % directivite antenne as C2 !
   direcv = 0.7;
end   

phi0        =  pi/N;
nph         =  8;
if ~isempty(gene.rapsauve) & strcmp(par.save,'Yes')
   [chemin,void] = fileparts(gene.rapsauve);
   save(fullfile(chemin,sprintf('last_delphine@%s',int2str(fix(gene.t*1000)))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAY-TRACING
% prepare les parametres d'entree pour le RT
initRT
% effectue les calcul de RT
boucleRT
% calcule les parametres des rayons utiles pour l'absorption
[nparfinal,wpefinal,Mfinal,Jfinal,rhofinal,Malpha,nperpfinal]=prepa(listeray,long,rfinal,fifinal,zfinal,krfinal,nfinal,kzfinal,param,R,Z,Bpr,Bpz,delta,ddeltadr,ddeltadz,rhobar,drhobardr,drhobardz,-gene.signe.ip.*gene.signe.b0);

% quand npar < 1 (ca peut arriver apres une mauvaise reflexion ...) ou Mfinal < 0 (signifie que le rayon est entre dans la zone de non-accessibilite),
% c'est qu'on a eu un gros probleme dans le calcul de la propagation --> on enleve la fin du rayon
for iray=listeray
  aa = find((abs(nparfinal(1:long(iray),iray))<=1) | (Mfinal(1:long(iray),iray) < 0));
  if ~isempty(aa)
	long(iray) = aa(1)-1;
	disp(num2str(iray));  
	nparfinal(aa(1):itermax,iray)=zeros(itermax-aa(1)+1,1);
	Mfinal(aa(1):itermax,iray)=zeros(itermax-aa(1)+1,1);
	Malpha(aa(1):itermax,iray)=zeros(itermax-aa(1)+1,1);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOKKER-PLANCK
% preparation des spectres pour le FP
Ph=abs(cons(1))*direch/1e6;
if any(size(cons)>1)
   Pv=abs(cons(2))*direcv/1e6;
   Ptotal=Ph+Pv;
else
   Pv = 0;
   Ptotal = Ph;
end      

% option absorption par les particules alpha
alpha = par.alpha;

% Zeff profile : initialise dans absopmhu

% initialisation des variables pour la boucle absorption + FP
absopmhu

Epar=prof.epar; % pour le FP (launchFP)
aa=find(isnan(Epar));  % blindage si Epar = NaN, au debut d'une simulation
if ~isempty(aa)
	Epar(aa)=zeros(size(aa));
end	 

if par.epar == 2    % initialise jlhlast dans le cas jlh et E// calcules de facon coherente
	if ~isnan(memoire.t) & isfield(memoire.data,'sortie')
	   jlhlast = memoire.data.sortie.j;
   else	
	   jlhlast = zeros(1,gene.nbrho);  % Attention, cela pose -t-il un probleme quand on demarre la simulation au milieu d'un creneau LH ?
   end
end

% boucle entre le module d'absorption et le FP
boucoup

if par.epar == 1
   % on redonne un coup de FP sans hybride pour deduire le courant ohmique
   coupohm
   % on retranche le courant ohmique du courant total
   courantLH=couranttotal-courantohm;
	renorm;
elseif par.epar == 0
   courantLH=couranttotal;
	renorm;
end
% remarque : le calcul du courant ohmique et la renormalisation sont traites directement dans couplage pour par.epar = 2 

% REMPLISSAGE DE LA STRUCTURE DE SORTIE
sortie.j = jlh;
sortie.el = depotfinal;

if sortie.err>=0      % si le calcul est OK
   if par.depip & isfield(memoire.data,'sortie')
      % si la structure memoire existe, et si le resultat n'est pas borne par Ip,
		% on verifie que l'ecart de courant par rapport a la simulation d'avant n'est pas aberrant
      if any(size(cons)>1)
         courant_interp = memoire.data.sortie.j.*(abs(cons(1))+abs(cons(2)))./memoire.data.ptot./(prof.ne+1e13).*memoire.data.ne;
      else
         courant_interp = memoire.data.sortie.j.*(abs(cons(1)))./memoire.data.ptot./(prof.ne+1e13).*memoire.data.ne;     
      end
      Ilh_interp = trapz(equi.rhomax*gene.x,courant_interp.*equi.spr);
      Ilhfinal = trapz(equi.rhomax*gene.x,sortie.j.*equi.spr);
		if abs((Ilhfinal-Ilh_interp)/Ilh_interp) > 1
      if strcmp(langue,'anglais')
		   disp('LH current too different from the previous DELPHINE run, the previous value is kept')
      else
		   disp('Courant trop different de l''iteration precedente : report de l''ancien resultat')
      end

         sortie.j = courant_interp;
         if any(size(cons)>1)
           sortie.el = memoire.data.sortie.el.*(abs(cons(1))+abs(cons(2)))./memoire.data.ptot;
         else
           sortie.el = memoire.data.sortie.el.*(abs(cons(1)))./memoire.data.ptot;        
         end
         sortie.err = -5;
		end			
   end
   if ~par.depip
	   % correction du courant si ihyb > ip_equi
      iboot=trapz(equi.rhomax*gene.x,neo.jboot.*equi.spr);
      % il faudrait peut-etre egalement tenir compte des autres sources de courant non-inductives ...
      Ilhfinal = trapz(equi.rhomax*gene.x,sortie.j.*equi.spr);
      if Ilhfinal > (0.95*equi.ip-iboot)
      if strcmp(langue,'anglais')
         disp('LH current correction as   Ilh greater than  (95% de Ip-Iboot)')
      else
         disp('Correction du courant LH car Ilh superieur a (95% de Ip-Iboot)')
      end

         sortie.j=sortie.j/Ilhfinal*(0.95.*equi.ip-iboot);
         sortie.err = (1-(0.95.*equi.ip-iboot)/Ilhfinal)*100; % le courant LH a ete reduit de sortie.err % 
		end
	end	
   % remplissage de la structure memoire
   mem.sortie   = sortie;
   if any(size(cons) > 1) 
     mem.ptot     = abs(cons(1))+abs(cons(2));
   else
     mem.ptot     = abs(cons(1));
   end
   mem.ne       = prof.ne; 
   memoire.t    = gene.t;
   memoire.data = mem;
else
   % s'il y a eu une erreur dans le calcul, on reprend le temps precedent avec les memes formules que dans zdelphe_interpolation
   if isfield(memoire.data,'sortie')
     if any(size(cons) > 1) 
       sortie.j = memoire.data.sortie.j.*(abs(cons(1))+abs(cons(2)))./memoire.data.ptot./(prof.ne+1e13).*memoire.data.ne;
       sortie.el = memoire.data.sortie.el.*(abs(cons(1))+abs(cons(2)))./memoire.data.ptot;
     else
       sortie.j = memoire.data.sortie.j.*(abs(cons(1)))./memoire.data.ptot./(prof.ne+1e13).*memoire.data.ne;
       sortie.el = memoire.data.sortie.el.*(abs(cons(1)))./memoire.data.ptot;    
     end
   else
      % si la structure memoire est vide, alors on sort des zeros
      sortie.j=zeros(size(gene.x));
      sortie.el=zeros(size(gene.x));
   end
end




if ~isempty(gene.rapsauve) & strcmp(par.save,'Yes')
   [chemin,void] = fileparts(gene.rapsauve);
   save(fullfile(chemin,sprintf('last_delphine@%s',int2str(fix(gene.t*1000)))));
end

