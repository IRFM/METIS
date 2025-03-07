% ZREMA interface cronos pour REMA (code de calcul des sources ECRH)
%------------------------------------------------------------------
% fichier zrema.m ->  zrema
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule les sources du chauffage ECRH.
% Elle utilise le mexfile rema.mexaxp
% 
% syntaxe  :
%  
%      [sortie,memoire] = zfcesimple(parametre,proto,cons,geo,equi,concentration, ...
%                          prof,neo,impur,phy,composition,gene)
%
% entree :
%
%      parametre       =    parametre propre a la fonction (param.cons.fce)
%      proto           =    prototype de la structure pour les sources, valeurs a zeros (proto = zsourceproto;)
%      cons            =    consigne de puissance par coupleur (data.cons.fce)
%      geo             =    geometrie du plasma (data.geo)
%      equi            =    donnees de l'equilibre plasma (data.equi)
%      concentration   =    profils de concentration des differentes espece d'ions/atomes (data.compo)
%      prof            =    profils des donnees calculees par le code (data.prof) 
%      neo             =    donnees neoclassiques (data.neo)
%      impur           =    sous strcuture des impurtes (data.impur)
%      phy             =    constantes physiques (param.phys)
%      composition     =    composition du plasma: charge et masse des atomes ( param.compo)
%      gene            =    parametres generaux (param.gne)
%      memoire         =    structure des dernieres valeurs calculees
%
% 
% sortie :
% 
%     sortie           =  structure de type source remplie par le module (sortie === proto)
%     memoire          =  datak.memoire.fce (valeur de reference pour le dernier calcul complet, 
%                         pas utiliser dans cette fonction, reserve pour d'autres modules)
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 26/09/2005
%
%
% liste des modifications :
%
%   * 02/10/2001 -> ajout de la structure memoire en sortie
%   * 09/11/2001 -> petit bug corrige ligne 302 (pa->ppa)
%   * 02/02/2002 -> ajout de l'equilibre complet
%   * 27/09/2002 -> on remplace b0 = geo.b0 par b0 = ge0.b0*signe.b0
%   * 30/09/2002 -> prise en compte complete des signes
%   * 01/09/2002 -> correction bug normalisation courant
%   * 24/01/2003 -> synergie simple en presence d'hybride
%   * 10/06/2003 -> sauvegarde du contexte sur demande pour test
%   * 13/06/2003 -> ajout de la possibilite de passer les angles des miroirs en consignes
%   * 26/09/2005 -> mise au propre de la gestion du nombre d'antennes
%   * 29/12/2005 -> passage possible avec 9 antennes
%--------------------------------------------------------------
%
function [sortie,memoire] = zremafile(parametre,proto,cons,geo,equi,concentration, ...
                             prof,neo,impur,phy,composition,gene,memoire)

% mode initialisation
% fonction auto declarante
if nargin <=1
	if nargin ==0
		nbfce=1;
	else
		nbfce=parametre;
	end
	langue                  = getappdata(0,'langue_cronos');

	valeur.nbcouronne   = 8;                        % nombre de couronnes de rayon (1 a 8) [8]
	valeur.nbray        = 2;                        % nombre de rayon par couronne * 7 (1 a 5) [2]
	valeur.modpolar     = 2 .* ones(1,nbfce);       % mode polarisation : 1 -> X, 2 -> O [2]
	valeur.freq_ghz     = 118 .* ones(1,nbfce);     % frequence en  Ghz [TS = 118]
	valeur.harm_min     = 1 .* ones(1,nbfce);       % harmonique min prise en compte [1]
	valeur.harm_max     = 2 .* ones(1,nbfce);       % harmonique max prise en compte [2]
	valeur.rant         = 3.53 .* ones(1,nbfce);    % position en R du dernier miroir (m) [TS = 3.53]
	if nbfce == 3
	   valeur.zant = [0,0.2,-0.2];                  % position en Z du dernier miroir (m) [valeurs par defaut antenne TS]
	else
           valeur.zant = zeros(1,nbfce);
	end
	valeur.angle_pol    = 0 .* ones(1,nbfce);       % angle d'injection dans le plan poloidal par rapport a l'horizontal, en degres (Ts = -20 a 20, >0 vers le haut) [0]
	valeur.angle_tor    = 0 .* ones(1,nbfce);      % angle d'injection dans le plan toroidal par rapport a la normal a Btor en degres (Ts = -30 a 30, >0 dans le sens de ip) [10]
        valeur.synergie     = 1 .* ones(1,nbfce);        % synergie sur le courant
        valeur.angle_var    = 0 .* ones(1,nbfce);        % angle des miroirs passes en consigne si =1
        valeur.equi         = 'cronos';	                  % sauvegarde des fichiers contextes
        valeur.save         = 'No';	                  % sauvegarde des fichiers contextes
        valeur.logout       = 'No';	                  % sauvegarde des fichiers contextes
        valeur.fastmode     = 'No';	                  % mode de cacul rapide

        valeur.curfor       = 1;	                  % sauvegarde des fichiers contextes

	type.nbcouronne     = 'integer';     % type de la donnnee
	type.nbray          = 'integer';     % type de la donnnee
	type.modpolar       = 'integer';     % type de la donnnee
	type.freq_ghz       = 'float';       % type de la donnnee
	type.harm_min       = 'integer';     % type de la donnnee
	type.harm_max       = 'integer';     % type de la donnnee
	type.rant           = 'float';       % type de la donnnee
	type.synergie       = 'float';       % type de la donnnee
	type.zant           = 'float';       % type de la donnnee
	type.angle_pol      = 'float';       % type de la donnnee
	type.angle_tor      = 'float';       % type de la donnnee
	type.angle_var      = 'integer';     % type de la donnnee
	type.equi           = 'string';     % type de la donnnee
	type.save               = 'string';
	type.logout               = 'string';

	type.fastmode               = 'string';

        type.curfor       = 'integer';	                  % sauvegarde des fichiers contextes
	
	borne.nbcouronne    = [1,32];         % bornes
	borne.nbray         = [1,8];
	borne.modpolar      = {1,2};
	borne.freq_ghz      = [1,1000];
	borne.harm_min      = [1,3];
	borne.harm_max      = [1,7];
	borne.rant          = [0,30];
	borne.synergie      = [-10,10];
	borne.zant          = [-20,20];
	borne.angle_pol     = [-90,90];
	borne.angle_tor     = [-90,90];
	borne.angle_var      = {0,1};
	borne.equi              = {'cronos','simple1','simple2'};
	borne.save              = {'Yes','No'};
	borne.logout              = {'Yes','No'};

	borne.fastmode              = {'Yes','No'};

        borne.curfor = [0 1];

	
	defaut.nbcouronne   = 8;           % valeur par defaut
	defaut.nbray        = 2;
	defaut.modpolar     = 2;
	defaut.freq_ghz     = 118;
	defaut.harm_min     = 1;
	defaut.harm_max     = 2;
	defaut.rant         = 3.53; 
	defaut.synergie     = 1; 
	defaut.zant         = 0;
	defaut.angle_pol    = 0;
	defaut.angle_tor    = 10;
	defaut.angle_var    = 0;
 	defaut.equi          = 'cronos';
 	defaut.save          = 'No';
 	defaut.logout          = 'No';

 	defaut.fastmode          = 'No';
 	

 	defaut.curfor = 1;


	  info.nbcouronne   = 'number of corrona (1 -> 8) [8]';
	  info.nbray        = 'number of rays per corrona * 7 (1 a 5) [2]';
	  info.modpolar     = 'polarisation mode : 1 -> X, 2 -> O [2]';
	  info.freq_ghz     = 'frequency [Ghz, TS = 118]';
	  info.harm_min     = 'harmonique min  [1]';
	  info.harm_max     = 'harmonique max  [2]';
	  info.rant         = 'last mirror major radius (m) [TS = 3.53]';
	  info.synergie     = 'synergy factor (/MW of LH) on EC current during LH';
	  info.zant         = 'last mirror altitude  (m) [0]' ;            
	  info.angle_pol    = 'poloidal injection angle [degree]  (Ts = -20 a 20, >0 vers le haut) [0]';
	  info.angle_tor    = 'toroidal injection angle [degree]  (Ts = -30 a 30, >0 dans le sens de ip) [10]';
	  info.angle_var    = 'if set at  1, antena mirrors angles are given as reference in input';
	  info.equi         = 'equilibrium : cronos, simple 1 (use edge and central safety factor), simple 2 (use current peaking factor';
	  info.save         = 'save context for test';
	  info.logout       = 'logout file (fort.9)';
          info.curfor       = 'current drive calculation (1 -> mode complet))';
	  info.fastmode     = 'tune rema for fast computation mode';
	
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'ECRH source term, using REMA';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	
	
	return
end


% debut du calcul
nba = gene.nbfce;
if nba > 9
	error('antenna number > 9, only 9 antennas !');
end

if ~isfield(parametre,'synergie')

  parametre.synergie = ones(1,nba);

end

if ~isfield(parametre,'angle_var')

  parametre.angle_var = 0 .* ones(1,nba);

end

if ~isfield(parametre,'fastmode')

  parametre.fastmode = 'No';

end


% seuil de puissance 
if all(abs(cons) < 1e4)
   sortie    = proto;
	return
end

% petits vecteurs utils


va = ones(gene.nbfce,1);
ve = ones(1,gene.nbrho);
betap      = zintsurf(prof.ptot,gene.x,equi.spr,equi.rhomax) ./  ...
             zintsurf(ones(size(gene.x)),gene.x,equi.spr,equi.rhomax) ./ ...
                      (prof.bpol(end) .^2 ./ 2 ./ phy.mu0);
[pja,pjq]=piquage(gene.x,equi.jmoy');
% preparation des donnes pour rema
indant     = find(abs(cons) > 1e4);
ipar       = sum(abs(cons) > 1e4);
switch parametre.fastmode
case 'Yes'
	ifa0       = 1;
	infa0      = 1;
	forcecor   = 1;
        parametre.curfor =0;
	
otherwise
	ifa0       = parametre.nbcouronne;
	infa0      = parametre.nbray;
	forcecor   = 0;
end
iout       = 0; % le mode iout = 1 pour debug, mais plante matlab5
a0         = equi.a(end);
r0         = equi.raxe(end);
beli       = betap+equi.li/2;
rcha       = geo.r0;
elon       = equi.e(1);
iq         = gene.signe.ip .* gene.signe.b0;
b0         = geo.b0;    %.* gene.signe.b0;
qp0        = equi.q(1);
qpa        = equi.q(end);
d0         = equi.d(1);
imodv      = comp(parametre.modpolar(indant));
frv        = comp(parametre.freq_ghz(indant));
Ipla       = equi.ip;
% decodage de la consigne
pfce       = abs(cons);
angle_mul  = angle(cons);
angle_toro = fix(abs(angle_mul) .* 1e7) .* 1e-4 - 360;
angle_polo = (abs(angle_mul) .* 1e7 - fix(abs(angle_mul) .* 1e7)) .* 1e3 - 360;
potv       = comp(abs(pfce(indant)));

% commutateur 
angle_toro = angle_toro .* (parametre.angle_var == 1) + parametre.angle_tor .* (parametre.angle_var == 0);
angle_polo = angle_polo .* (parametre.angle_var == 1) + parametre.angle_pol .* (parametre.angle_var == 0);

% l'angle toroidal positif correspond a du co courant (absorbtion downshift sans okawa)
angle_toro = -gene.signe.b0 .* angle_toro;

iar1v      = comp(parametre.harm_min(indant));
iar2v      = comp(parametre.harm_max(indant));
rav        = comp(parametre.rant(indant));
zav        = comp(parametre.zant(indant));
arzv       = comp(angle_polo(indant));
dphk0v     = comp(angle_toro(indant));
freq0      = 1; 
dfr2       = frv(1)+2;
rifle1     = 0.75;

%
% couplage et interpolation avec l'equilibre
%
rho        = double(equi.rhoRZ);
x          = rho ./ rho(end);

R          = squeeze(double(equi.R));
Z          = squeeze(double(equi.Z));
%
% sortie d'HELENA de BR BZ corresponde a un courant dans le sens des aiguilles d'une montre vue du dessus du tokamak
%
BR         = -gene.signe.ip .* squeeze(double(equi.BR));
BZ         = -gene.signe.ip .* squeeze(double(equi.BZ));
BPHI       =  gene.signe.b0 .* squeeze(double(equi.BPHI));
%
% changement de grille pour les profils
%
te         = interp1(gene.x,prof.te,x,'linear') ;
ne         = interp1(gene.x,prof.ne,x,'linear') ;
ze         = interp1(gene.x,prof.zeff,x,'linear') ;
spr        = interp1(gene.x,equi.spr,x,'linear') ;
vpr        = interp1(gene.x,equi.vpr,x,'linear') ;
vpr(end)  = 2.* vpr(end-1) - vpr(end-2);
spr(end)  = 2.* spr(end-1) - spr(end-2);
%save newprofiles te ne vpr -V4
% securite anti nan
te(~isfinite(te)) = 0;
ne(~isfinite(ne)) = 0;
ze(~isfinite(ze)) = 0;
ze(~isfinite(ze)) = 0;
ze = abs(ze);


%C		INPUT DATA DESCRIPTION
%c	
%c	IPAR:	Number of wave beams launched into the plasma (1 - 3)
%c		Each wave beam is described by many rays, distributed
%c		on IFA0 concentric circles, 7*INFA0 equally spaced rays
%c		per circle, plus the central ray of the beam
%c	IFA0:	Number of concentric circles used to describe the beam
%c		Suggested values: 1 - 32
%c	INFA0:	7*INFA0 is the number of rays per circle
%c		Suggested values: 1 - 2
%c	IOUTW:	IOUTW=1 -> the matlab file with all the graphic outputs
%c		is written; IOUTW=0 -> the file is not written
%c
%c	Plasma parameters
%c
%c	A0W:	plasma minor radius in cm
%c	R0W:	plasma major radius in cm
%c	RCHAW:	major radius (in cm) entering the definition of the central
%c		magnetic field: B0=B(RCHAW)
%c	ELONW:	Elongation of the elliptical cross-section
%c	IQW:	Sign of the current with respect to the magnetic field
%c		(+1 or -1)
%c	FB:	for FB=1 the equilibrium is determined by the central
%c		magnetic field B0, central and edge safety factor QP0W
%c		and QPAW, maximum Shafranov shift D0.
%c		for FB=0 the values of IBTORw,IPLw,BELIw,ALPHAP1w are used
%c	B0W:	Central magnetic field in Tesla: B0W=B(RCHAW)
%c	QP0W:	Central value of the safety factor (the profile is parabolic)
%c	QPAW:	Edge value of the safety factor
%c	D0W:	Shafranov shift (in cm) of the magnetic axis (the profile
%c		is linear in the normalised magnetic flux)
%c	FI:	for FI=2 the equilibrium files are used, i.e., the normalized
%c		radius rhov (array dimension ic) and the quantities 
%c		RV,ZV,BPHIV,BRV,BZV (array dimension ic, it).
%c               for FI=1 the equilibrium is determined by the current 
%c		flowing in the toroidal magnetic field coils IBTORW, the
%c		plasma current IPLW, beta poloidal + li/2 BELIW, and the
%c		current peaking factor ALPHAP1W
%c		for FI=0 the values of B0w,QP0w,QPAw,D0w are used
%c	IBTORW:	current flowing in the toroidal magnetic field coils in Amp.
%c	IPLW:	Plasma current in MA
%c	BELIW:	beta poloidal + li/2
%c	ALPHAP1W:	current peaking factor
%c
%c	ECE parameters
%c
%c	IFRE0:	number of frequencies used in the calculation of the ECE
%c		spectrum <= IFF (array dimension)
%c	DFR2:	extension of the frequency range (in GHz) for the calculation
%c		of the spectrum (the same for the three cases)
%c	RIFLE1:	wall reflection coefficient (between 0 and 1)
%c
%c	IC:	Number of points for the density, temperature and Zeff profiles
%c       IC2M1:  2*IC-1
%c	RHOV:	Array of normalised radius for the profiles
%c	NEV:	Array of density profile
%c	TEV:	Array of temperature profile
%c	ZEV:	Array of Zeff profile
%c	
%c	Wave parameters
%c	
%c	IMODV:	Wave polarisation mode:	IMODV=1 -> extraordinary mode
%c					IMODV=2 -> ordinary mode
%c	FRV:	Wave frequency in GHz
%c	POTV:	Wave power in kW
%c	IAR1V:	Lower boundary of the range of harmonics included in
%c		the calculation.  Suggested value: 1.
%c	IAR2V:	Upper boundary of the range of harmonics included in
%c		the calculation.  Suggested value: 2.
%c	RAV:	Major radius of the antenna location (in cm)
%c	ZAV:	Vertical antenna location with respect to the equatorial
%c		plane (in cm)
%c	ARZV:	Projection in a poloidal plane of the wave launching
%c		angle, with respect to the horizontal direction.
%c		> 0: wave launched upwards; < 0: wave launched downwards
%c	DPHK0V:	Projection in a toroidal plane of the wave launching
%c		angle, with respect to the normal to the magnetic field.
%c		> 0: in the direction of the plasma current.
%c
%c
%c		OUTPUT DATA DESCRIPTION
%c
%c	PPA:	Power deposition profile (in kW/cm^3) for the 1st wave beam
%c	JA:	Driven current profile (in kA/cm^2) for the 1st wave beam
%c		The current is computed using a linear formula including
%c		trapped electrons effects (R. Cohen, Phys. Fluids 30, 2442 (1987))
%c
%c	PPB:	Power deposition profile (in kW/cm^3) for the 2nd wave beam
%c	JB:	Driven current profile (in kA/cm^2) for the 2nd wave beam
%c
%c	PPC:	Power deposition profile (in kW/cm^3) for the 3rd wave beam
%c	JC:	Driven current profile (in kA/cm^2) for the 3rd wave beam
%c
%c	ECE outputs
%c	
%c	FREV1:	Frequency array (in GHz) of the 1st ECE spectrum
%c	TRADV1:	Radiation temperature (in keV) of the 1st ECE spectrum
%c	TAFV1:	Optical depth associated to the 1st ECE spectrum
%c	FREV2:	Frequency array (in GHz) of the 2nd ECE spectrum
%c	TRADV2:	Radiation temperature (in keV) of the 2nd ECE spectrum
%c	TAFV2:	Optical depth associated to the 2nd ECE spectrum
%c	FREV3:	Frequency array (in GHz) of the 3rd ECE spectrum
%c	TRADV3:	Radiation temperature (in keV) of the 3rd ECE spectrum
%c	TAFV3:	Optical depth associated to the 3rd ECE spectrum
%c
%c
if ~isfield(parametre,'equi')
    parametre.equi = 'cronos';
end
if ~isfield(parametre,'logout')
    parametre.logout = 'No';
end
if strcmp(parametre.equi,'cronos')
  FI(1) = 2;
elseif strcmp(parametre.equi,'simple1')
  FI(1) = 0;
  b0    = gene.signe.b0 * b0;
else
  FI(1) = 1;
end
%
% formule pour le calcul du courant
%
FI(2) = parametre.curfor;
remacallgen

if strcmp(parametre.save,'Yes')

  save contexte_rema

end
nbr = length(ind);
% appel de rema
%[xo,ppa,ja,ppb,jb,ppc,jc]   = remamex(ipar,ifa0,infa0,iout, ...
%				      a0,r0,rcha,elon,iq, ...
%				      b0,qp0,qpa,d0, ...
%				      x,te,ne,ze, ...
%				      imodv,frv,potv,iar1v,iar2v, ...
%				      rav,zav,arzv,dphk0v);

%changement de grille pour les sortie equilibre -> cronos
wa          = equi.rhomax .* trapz(xx,vpr .* psauve(1,:));
ppa         = interp1(xx,psauve(1,:),gene.x,'linear') ;
wb          = equi.rhomax .* trapz(xx,vpr .* psauve(2,:));
ppb         = interp1(xx,psauve(2,:),gene.x,'linear') ;
wc          = equi.rhomax .* trapz(xx,vpr .* psauve(3,:));
ppc         = interp1(xx,psauve(3,:),gene.x,'linear') ;
ia          = equi.rhomax.* trapz(xx,spr .* jsauve(1,:));
ja          = interp1(xx,jsauve(1,:),gene.x,'linear') ;
ib          = equi.rhomax .* trapz(xx,spr .* jsauve(2,:));
jb          = interp1(xx,jsauve(2,:),gene.x,'linear') ;
ic          = equi.rhomax .* trapz(xx,spr .* jsauve(3,:));
jc          = interp1(xx,jsauve(3,:),gene.x,'linear') ;
if nbr > 3 & nbr < 7
 wd          = equi.rhomax .* trapz(xx,vpr .* psauve(4,:));
 ppd         = interp1(xx,psauve(4,:),gene.x,'linear') ;
 we          = equi.rhomax .* trapz(xx,vpr .* psauve(5,:));
 ppe         = interp1(xx,psauve(5,:),gene.x,'linear') ;
 wf          = equi.rhomax .* trapz(xx,vpr .* psauve(6,:));
 ppf         = interp1(xx,psauve(6,:),gene.x,'linear') ;
 id          = equi.rhomax.* trapz(xx,spr .* jsauve(4,:));
 jd          = interp1(xx,jsauve(4,:),gene.x,'linear') ;
 ie          = equi.rhomax .* trapz(xx,spr .* jsauve(5,:));
 je          = interp1(xx,jsauve(5,:),gene.x,'linear') ;
 iff         = equi.rhomax .* trapz(xx,spr .* jsauve(6,:));
 jf          = interp1(xx,jsauve(6,:),gene.x,'linear') ;
end
if nbr > 7
 wg          = equi.rhomax .* trapz(xx,vpr .* psauve(7,:));
 ppg         = interp1(xx,psauve(7,:),gene.x,'linear') ;
 wh          = equi.rhomax .* trapz(xx,vpr .* psauve(8,:));
 pph         = interp1(xx,psauve(8,:),gene.x,'linear') ;
 wi          = equi.rhomax .* trapz(xx,vpr .* psauve(9,:));
 ppi         = interp1(xx,psauve(9,:),gene.x,'linear') ;
 ig          = equi.rhomax.* trapz(xx,spr .* jsauve(7,:));
 jg          = interp1(xx,jsauve(7,:),gene.x,'linear') ;
 ih          = equi.rhomax .* trapz(xx,spr .* jsauve(8,:));
 jh          = interp1(xx,jsauve(8,:),gene.x,'linear') ;
 ifi         = equi.rhomax .* trapz(xx,spr .* jsauve(9,:));
 ji          = interp1(xx,jsauve(9,:),gene.x,'linear') ;
    
end


% gestion des sorties
sortie    = proto;
switch ipar
case 1
	sortie.el = ppa(:)';
	sortie.j = ja(:)';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1))';
	iece     = ia;
	wece     = wa;
case 2
	sortie.el = (ppa(:) + ppb(:))';
	sortie.j = (ja(:) + jb(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1))';
	iece     = ia + ib;
	wece     = wa + wb;
case 3
	sortie.el = (ppa(:) + ppb(:) + ppc(:))';
	sortie.j = (ja(:) + jb(:) + jc(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1) + ...
                     jc(:) * (parametre.synergie(indant(3))-1))';
	iece     = ia + ib + ic;
	wece     = wa + wb + wc;
case 4
	sortie.el = (ppa(:) + ppb(:) + ppc(:) + ppd(:))';
	sortie.j = (ja(:) + jb(:) + jc(:) + jd(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1) + ...
                     jc(:) * (parametre.synergie(indant(3))-1) + ...
                     jd(:) * (parametre.synergie(indant(4))-1))';
	iece     = ia + ib + ic + id;
	wece     = wa + wb + wc + wd;
case 5
	sortie.el = (ppa(:) + ppb(:) + ppc(:) + ppd(:) + ppe(:))';
	sortie.j = (ja(:) + jb(:) + jc(:) + jd(:) + je(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1) + ...
                     jc(:) * (parametre.synergie(indant(3))-1) + ...
                     jd(:) * (parametre.synergie(indant(4))-1) + ...
                     je(:) * (parametre.synergie(indant(5))-1))';
	iece     = ia + ib + ic + id + ie;
	wece     = wa + wb + wc + wd + we;
case 6
	sortie.el = (ppa(:) + ppb(:) + ppc(:) + ppd(:) + ppe(:) + ppf(:))';
	sortie.j = (ja(:) + jb(:) + jc(:) + jd(:) + je(:) + jf(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1) + ...
                     jc(:) * (parametre.synergie(indant(3))-1) + ...
                     jd(:) * (parametre.synergie(indant(4))-1) + ...
                     je(:) * (parametre.synergie(indant(5))-1) + ...
                     jf(:) * (parametre.synergie(indant(6))-1))';
	iece     = ia + ib + ic + id + ie + iff;
	wece     = wa + wb + wc + wd + we + wf;
case 7
	sortie.el = (ppa(:) + ppb(:) + ppc(:) + ppd(:) + ppe(:) + ppf(:) + ppg(:))';
	sortie.j = (ja(:) + jb(:) + jc(:) + jd(:) + je(:) + jf(:) + jg(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1) + ...
                     jc(:) * (parametre.synergie(indant(3))-1) + ...
                     jd(:) * (parametre.synergie(indant(4))-1) + ...
                     je(:) * (parametre.synergie(indant(5))-1) + ...
                     jf(:) * (parametre.synergie(indant(6))-1) + ...
                     jg(:) * (parametre.synergie(indant(7))-1))';
	iece     = ia + ib + ic + id + ie + iff + ig;
	wece     = wa + wb + wc + wd + we + wf + wg;
case 8
	sortie.el = (ppa(:) + ppb(:) + ppc(:) + ppd(:) + ppe(:) + ppf(:) + ppg(:) + pph(:))';
	sortie.j = (ja(:) + jb(:) + jc(:) + jd(:) + je(:) + jf(:) + jg(:) + jh(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1) + ...
                     jc(:) * (parametre.synergie(indant(3))-1) + ...
                     jd(:) * (parametre.synergie(indant(4))-1) + ...
                     je(:) * (parametre.synergie(indant(5))-1) + ...
                     jf(:) * (parametre.synergie(indant(6))-1) + ...
                     jg(:) * (parametre.synergie(indant(7))-1) + ...
                     jh(:) * (parametre.synergie(indant(8))-1))';
	iece     = ia + ib + ic + id + ie + iff + ig + ih;
	wece     = wa + wb + wc + wd + we + wf + wg + wh;
case 9
	sortie.el = (ppa(:) + ppb(:) + ppc(:) + ppd(:) + ppe(:) + ppf(:) + ppg(:) + pph(:) + ppj(:))';
	sortie.j = (ja(:) + jb(:) + jc(:) + jd(:) + je(:) + jf(:) + jg(:) + jh(:) + jj(:))';
   sortie.synergie = (ja(:) * (parametre.synergie(indant(1))-1) + ...
                     jb(:) * (parametre.synergie(indant(2))-1) + ...
                     jc(:) * (parametre.synergie(indant(3))-1) + ...
                     jd(:) * (parametre.synergie(indant(4))-1) + ...
                     je(:) * (parametre.synergie(indant(5))-1) + ...
                     jf(:) * (parametre.synergie(indant(6))-1) + ...
                     jg(:) * (parametre.synergie(indant(7))-1) + ...
                     jh(:) * (parametre.synergie(indant(8))-1) + ...
                     jj(:) * (parametre.synergie(indant(9))-1))';
	iece     = ia + ib + ic + id + ie + iff + ig + ih + ij;
	wece     = wa + wb + wc + wd + we + wf + wg + wh + wj;

otherwise
   sortie.synergie = zeros(size(gene.x));


end

rap = trapz(gene.x,sortie.el .* equi.vpr).*equi.rhomax ./ sum(abs(cons));
iece_out  = equi.rhomax .* trapz(gene.x,equi.spr .* sortie.j);
fprintf('Pfce/cons = %g\n',rap);
fprintf('Prema/cons = %g\n',wece ./ sum(abs(cons)));
fprintf('Iout/Irema = %g\n',iece_out ./ iece);

if ((rap > 1.2) | (rap < 0.8)) | forcecor
        fprintf('renormalisation to the reference power (MW) : %g\n',cons/1e6);
	sortie.el       = sortie.el ./ rap;
	sortie.j        = sortie.j ./ rap;
        sortie.synergie = sortie.synergie ./ rap;
end


% compelete les parametres
function s = comp(e)

s = zeros(1,3);
s(1:length(e)) = e(:)';


% integralle surfacique (d'une grandeur independante de theta)
%  s = integrale de surface
%  e = valeur a integree
%  x = coordonnees normalisee
%  sp = datak.equi.sp
%  rhomax = datak.equi.rhomax   
function s=zintsurf(e,x,sp,rhomax)   

    s = rhomax .* trapz(x,sp .* e,2);
  
  
% integralle volumique 
%  s = integrale de volume
%  e = valeur a integree
%  x = coordonnees normalisee
%  vpr = datak.equi.vpr
%  rhomax = datak.equi.rhomax   
function s=zintvol(e,x,vpr,rhomax)   

  s = rhomax.*trapz(x,vpr .* e,2);
