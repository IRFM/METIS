% ZMURGAZ calcule simple le bord et l'etat du mur 
%------------------------------------------------
% fichier zmurgaz.m ->  zmurgaz, zintvol
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule (de maniere simple) les donnees du bord du plasma et l'etat du mur.
% elle gere uniquement le recyclage, le pompage (efficacite 100%) et l'injection direct de gaxz 
%  
% syntaxe  :
%  
%     [murkp1,memoire]=zrecycle(par,murk,pomp,gaz,prof,impur,source,compo,memoire,vpr,spr,rhomax,grho2,x,phys,nbg,dt)
%
% entrees :
%
%     par              = struture des consignes de calcul de la fonction (param.cons.bord)
%     murk             = structure "bord" au temps precedent (datak.bord)
%     pomp             = consigne du pompage (datakp1.cons.pomp)
%     gaz              = consigne d'injection de gaz au bord (datakp1.cons.c)
%     prof             = structure des profils (datakp1.prof)
%     impur            = structure des impuretes et du rayonnement (datakp1.impur)
%     source           = structure des sources (datakp1.source)
%     compo            = composition du plasma (param.compo)
%     memoire          = datak.memoire.bord (valeur de reference pour le dernier calcul complet, 
%                        pas utiliser dans cette fonction, reserve pour d'autres modules)
%     vpr              = surface du plasma (datakp1.equi.vpr)
%     spr              = circonference du plasma (datakp1.equi.spr)
%     rhomax           = rayon moyen du plasma (datakp1.equi.rhomax)
%     grho2            = <|gradient(rho)|^2> (datakp1.equi.grho2)
%     x                = param.gene.x
%     phys             = constante physiques (param.phys)
%     nbg              = nombre d'especes considerees (param.gene.nbg)
%     dt               = pas de temps local
%
% sorties :
% 
%     murkp1     = structure d bord au temps k+1 (datakp1.bord) 
%     memoire    =  datak.memoire.bord (valeur de reference pour le dernier calcul complet, 
%                   pas utiliser dans cette fonction, reserve pour d'autres modules)
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0 08/12/2004.
% 
% 
% liste des modifications : 
%
%   * 26/06/2001 -> temperature de bord
%   * 27/06/2001 -> (ne,ni,Te,Ti) self consistant
%   * 02/10/2001 -> ajout de la structure memoire en sortie
%   * 25/01/2002 -> correction du calcul de Te et Ti bord
%   * 04/02/2002 -> changement de la boucle de convergence pour Te et Ti bord
%   * 13/02/2002 -> nouveau calcul des donnees de bord
%   * 22/02/2002 -> ajout de parametres (fixer avant dans le code)
%   * 20/03/2002 -> correction du lissage si dt > dtlissage
%   * 05/04/2003 -> securite sur les consignes (anti NaN et Inf)
%   * 07/04/2004 -> securite sur le flux sortant
%   * 08/04/2004 -> ajout du recyclage (parametre)
%   * 08/12/2004 -> protection densite nul au bord
%--------------------------------------------------------------
%
function [murkp1,memoire]=zrecycle(par,murk,pomp,gaz,prof,impur,source,compo,memoire,vpr,spr,rhomax,grho2,x,phys,nbg,dt)

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 

	valeur.gain_te   = 40;  % facteur multiplicatif de la loi donnant te de bord (eV)
	valeur.gain_ti   = 40;  % facteur multiplicatif de la loi donnant ti de bord (eV)
	valeur.te_max    = 200; % valeur maximale de te de bord (eV)
	valeur.ti_max    = 200; % valeur maximale de ti de bord (eV)
	valeur.lissage   = 0.1;   % constante de temps d'evolution des temperature de bord et de la densite de bord (s)
	valeur.recycle   = 1;   % facteur de recyclage reel

	type.gain_te     = 'float';     % type reel
	type.gain_ti     = 'float';     % type reel
	type.te_max      = 'float';     % type reel
	type.ti_max      = 'float';     % type reel
	type.lissage     = 'float';     % type reel
	type.recycle      = 'float';     % type reel

	borne.gain_te    = [1,1000];  % valeurs possible
	borne.gain_ti    = [1,1000];  % valeurs possible
	borne.te_max     = [13.6,1000];  % valeurs possible
	borne.ti_max     = [13.6,1000];  % valeurs possible
	borne.lissage   = [0,1];       % valeurs possible 
	borne.recycle   = [0,1.5];       % valeurs possible 

	defaut.gain_te   = 40;   % valeurs par defaut
	defaut.gain_ti   = 40;   % valeurs par defaut
	defaut.te_max    = 200; 
	defaut.ti_max    = 200; 
	defaut.lissage   = 0.1;  % valeurs par defaut
	defaut.recycle   = 1;  % valeurs par defaut

	% informtions
	info.gain_te    = 'facteur multiplicatif de la loi donnant te de bord (eV)';
	info.gain_ti    = 'facteur multiplicatif de la loi donnant ti de bord (eV)';
	info.te_max     = 'valeur maximale de te de bord (eV)';
	info.ti_max     = 'valeur maximale de ti de bord (eV)';
	info.lissage    = 'constante de temps d''evolution des temperature de bord et de la densite de bord (s)';
	info.recycle    = 'facteur de recyclage effectif tel que la source au bord = R * flux sortant';

	interface.ts        = 'a faire ts';   % nom de la fonction d'interfacage avec les donnees TS
	interface.jet       = 'a faire jet';  % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'Calcul tres simple du bord et de l''etat du mur';   % description (une ligne) de la fonction
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	murkp1 = sortie;
	return

end

% le proto en sortie
murkp1 = murk;

% non calculer
murkp1.pref       = 400;
murkp1.temp       = 400;
murkp1.reflex = 1;
murkp1.nb     = 1;
murkp1.pchauffe = 0;

% flux en provenance du plasma
fluxout = prof.flux.sortant.ge(end) ;
% le flux doit etre positif
% par default on calcul pour 1 particule/s
fluxout = max(1,fluxout .* (fluxout >0));

% securite
gaz(~isfinite(gaz)) =0;
% securite
if ~isfinite(pomp)
   pomp = 0;
end
if fluxout == 1
   fluxout = abs(prof.flux.sortant.ge(end)) ;
   zverbose('flux sortant de matiere <= 0 !\n');
   if fluxout == 1
      fluxout = pomp + murk.fluxplasma;
   end
   if fluxout <= 1
      fluxout =  prof.ne(end) ./ 10;
   end
end
if isfield(par,'recycle')
   fluxout = fluxout .* par.recycle;
end

% flux pompe en equivalent electrons
murkp1.fluxpompe  = min(fluxout,pomp);
% flux reentrant en equivalent electrons
murkp1.fluxplasma = max(1,fluxout-pomp);
% rapport du flux reentrant au flux sortant
rapflux   = murkp1.fluxplasma ./ (fluxout+eps);
if fluxout == 0
	rapflux = 1;
end

% flux sortant par especes 
% hypothese proportionnalite au flux electronique equivalent
%for k =1:nbg
%	murkp1.fluxgazout(k) = fluxout .* squeeze(impur.impur(1,end,k))' ./  ...
%	                       sum(squeeze(impur.impur(1,end,1:nbg)).*compo.z');
%end 
%
% protection densite nul au bord, VB 8 decembre 2004
%
for k =1:nbg
	murkp1.fluxgazout(k) = fluxout .* squeeze(impur.impur(1,end,k))' ./  ...
	                       sum(squeeze(impur.impur(1,end,1:nbg)+1e16).*compo.z');
end 

% renormalisation pour conserver les charges
%
% protection quand murkp1.fluxgazout est nul VB, 8 decembre 2004
% ind = find(~isfinite(murkp1.fluxgazout));

murkp1.fluxgazout(~isfinite(murkp1.fluxgazout))=0;
indz = find(murkp1.fluxgazout~=0);
if ~isempty(indz)
  murkp1.fluxgazout(indz)  = murkp1.fluxgazout(indz) .* (murkp1.fluxplasma ./ sum( murkp1.fluxgazout(indz).*compo.z(indz)));
else
  murkp1.fluxgazout         = zeros(1,nbg);
end

% flux de gaz (fueling)                     
murkp1.fluxgaz    = gaz;
ind  = find(~isfinite(murkp1.fluxgaz));
if ~isempty(ind)
     murkp1.fluxgaz(ind) = 1; % 1 particule par s
end
% flux entrant dans le plasma (en equivalent electrons)
% neutres chauds
murkp1.fluxmur_c    = murkp1.fluxplasma;
%neutres froids
murkp1.fluxmur_f    = sum(murkp1.fluxgaz .* compo.z);
% le flux total entrant dans le plasma
murkp1.fluxgebord   = murkp1.fluxmur_c + murkp1.fluxmur_f;

% flux entrant par espece 
% neutres chauds
murkp1.fluxgazin_c    = murkp1.fluxgazout;
% neutres froids
murkp1.fluxgazin_f    = murkp1.fluxgaz;

% flux sortant comme conditions au limites
% energie des neutres froids (eV)
e0   = murkp1.temp .* phys.k ;

qebord            =  max(prof.flux.sortant.qe(end),0);
qebordmax         =  max(rhomax .* trapz(x,vpr .* abs(source.totale.el),2),0);
qebordmin         =  max(rhomax .* trapz(x,vpr .* prof.pe,2)./vpr(end),0);
if isfinite(qebordmax)
    qebord            =  min(max(qebord,qebordmin),qebordmax);
else
    qebord            = qebordmin;
end
qibord            =  max(prof.flux.sortant.qi(end),0);
qibordmax         =  max(rhomax .* trapz(x,vpr .* abs(source.totale.ion),2),0);
qibordmin         =  max(rhomax .* trapz(x,vpr .* prof.pion,2)./vpr(end),0);
if isfinite(qibordmax)
    qibord            =  min(max(qibord,qibordmin),qibordmax);
else
    qibord            = qibordmin;
end

murkp1.fluxqebord = max(0,murkp1.fluxmur_c ./ max(1,fluxout) .* qebord - 25  .* phys.e .* murkp1.fluxmur_f);
murkp1.fluxqibord = murkp1.fluxmur_c ./ max(1,fluxout) .* qibord + e0 .* murkp1.fluxmur_f;
if ~isfinite(murkp1.fluxqebord )
	murkp1.fluxqebord = 0;
end
if ~isfinite(murkp1.fluxqibord )
	murkp1.fluxqibord = 0;
end

% initialisation des variables
murkp1.nebord   = min(prof.ne(prof.ne > 0));
nebord          = murkp1.nebord;
murkp1.nibord   = squeeze(min(impur.impur,[],2))';
nibord          = murkp1.nibord;
murkp1.tebord   = min(prof.te(prof.te > 0));
murkp1.tibord   = min(prof.ti(prof.ti >0));
if ~isfinite(murkp1.nebord)
	murkp1.nebord = 3e17;
end
ind =   find(~isfinite(murkp1.nibord));
if ~isempty(ind)
	 murkp1.nibord(ind)   = 0;
end
if murkp1.nibord(1) == 0
	murkp1.nibord(1) = murkp1.nebord ./ compo.z(1);
end

if ~isfinite(murkp1.tebord)
	murkp1.tebord = 30;
end
if ~isfinite(murkp1.tibord)
	murkp1.tibord = 30;
end

% calcul de la densite de bord autonome pour calculer les Te et Ti de bord(modele du Wesson)
nbar          = trapz(x,prof.ne);
nea           = znebord(nbar);
nia           = max(1e17,prof.ae(end) .* nea);
% calcul des densite avec le flux
nebord_new = nea;
%
nibord_new      = (murkp1.fluxgazin_f(:) + murkp1.fluxgazin_c(:))';
ind                = find(~isfinite(nibord_new));
nibord_new(ind) = zeros(1,length(ind));
neeq             = sum(nibord_new(:).*compo.z(:));
if neeq == 0
	neeq = nea;
end
if all(nibord_new == 0)
	nibord_new(1) = nia;
end
% normalisation
nibord_new = nibord_new ./ neeq .* nea;

% lissage temporel
if dt > par.lissage
  	murkp1.nebord = nebord_new;
  	murkp1.nibord = nibord_new;
else
  	murkp1.nebord = (murkp1.nebord .* (par.lissage - dt) + nebord_new * dt) ./ par.lissage;
  	murkp1.nibord = (murkp1.nibord .* (par.lissage - dt) + nibord_new * dt) ./ par.lissage;
end

% calcul des temperatures temperature avec une loi simple
%
% nibord = nibord + 1e16, VB, 8 decembre 2004
%
tebord_new   = min(max(par.gain_te .* (qebordmax ./ 1e6) .^(2/3) ./ (nebord ./ 1e18),30),par.te_max); 
tibord_new   = min(max(par.gain_ti .* (qibordmax ./ 1e6) .^(2/3) ./ (sum(nibord+1e16) ./ 1e18),30),par.ti_max); 
if ~isfinite(tebord_new)
	tebord_new = 30;
end
if ~isfinite(tibord_new)
	tibord_new = 30;
end
% lissage temporel
if dt > par.lissage
  	murkp1.tebord = tebord_new;
  	murkp1.tibord = tibord_new;
else
  	murkp1.tebord = (murkp1.tebord .* (par.lissage - dt) + tebord_new * dt) ./ par.lissage;
  	murkp1.tibord = (murkp1.tibord .* (par.lissage - dt) + tibord_new * dt) ./ par.lissage;
end

% integralle volumique 
%  s = integrale de volume
%  e = valeur a integree
%  x = coordonnees normalisee
%  vpr = datak.equi.vpr
%  rhomax = datak.equi.rhomax   
function s=zintvol(e,x,vpr,rhomax)   

  s = rhomax.*trapz(x,vpr .* e,2);
