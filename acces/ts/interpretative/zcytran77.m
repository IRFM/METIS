% ZCYTRAN77 profil de depot fci calcule en utilisant Pion et / ou absor
%--------------------------------------------------------------------
% fichier zcytran77.m ->  mexfile cytran_77
%
% fonction Matlab 7 :
%
% Cette fonction calcule les pertes (profil) de rayonnement cyclotronique en utilisant Cytran
% version fortran 77 
% 
% syntaxe  :
%  
%      pcyclo = zcytran(parametre,geo,equi,prof,neo,impur,phys,composition,gene);
%
% entree :
%
%      parametre       =    parametre propre a la fonction (param.cons.sync)
%      geo             =    geometrie du plasma (data.geo)
%      equi            =    donnees de l'equilibre plasma (data.equi)
%      injection       =    consigne d'injection de gaz (data.cons.c)
%      prof            =    profils des donnees calculees par le code (data.prof) 
%      neo             =    donnees neoclassiques (data.neo)
%      impur           =    sous strcuture des impurtes (data.impur)
%      phys             =    constantes physiques (param.phys)
%      composition     =    composition du plasma: charge et masse des atomes ( param.compo)
%      gene            =    parametres generaux (param.gne)
%      memoire         =    structure des dernieres valeurs calculees
% 
% sortie :
% 
%     pcyclo           =  puissance cyclotronique perdue (W/m^3)
%                         ! attention au signe !
%      memoire         =    structure des dernieres valeurs calculees
% 
%
% fonction ecrite par V. Basiuk , poste 61_26
% version vcvs.
% 
%--------------------------------------------------------------
%
function [pcyclo,memoire] = zcytran(parametre,geo,equi,prof,neo,impur,phys,composition,gene,memoire)

% fonction auto declarante                             
if nargin <=1 

	  valeur.ree        = 0.76;                     % refl of X mode from incident X mode [0-1]
	  valeur.reo        = 0.04;                      % refl of X mode from incident O mode [0-1]
  	  valeur.roo        = 0.76;                   % refl of O mode from incident O mode [0-1]
	  valeur.roe        = 0.04;                   % refl of O mode from incident X mode [0-1]
	  valeur.albajar    = 1;                   % renormalisation sur le fit de albajar
	  valeur.reflex        = 0.8;         % coefficient de reflexion du rayonnement cyclo par les parois

	  type.ree        = 'float';                               % type reel
	  type.reo        = 'float';                               % type reel
	  type.roo        = 'float';                               % type reel
	  type.roe        = 'float';                               % type reel
	  type.albajar    = 'integer';                               %
	  type.reflex     = 'float';     % type reel

	  borne.ree       = [0,1];                              % valeurs possible
	  borne.reo       = [0,1];                              % valeurs possible
	  borne.roo       = [0,1];                               % valeurs possibles
	  borne.roe       = [0,1];                           % valeurs possibles
	  borne.albajar   = {0,1};                           % valeurs possibles
	  borne.reflex    = [0,1];       % valeurs possible
	
	  defaut.ree      = 0.76;                                   % valeurs par defaut
	  defaut.reo      = 0.04;                                     % valeurs par defaut
	  defaut.roo      = 0.76;                                  % valeur par defaut
	  defaut.roe      = 0.04;                                  % valeur par defaut
	  defaut.albajar  = 1;                                  % valeur par defaut
	  defaut.reflex   = 0.8;         % valeurs par defaut
      
     langue = getappdata(0,'langue_cronos');
 
     if strcmp(langue,'francais')
	    info.ree        = 'refl of X mode from incident X mode [0-1]';
	    info.reo        = 'refl of X mode from incident O mode [0-1]';
	    info.roo        = 'refl of O mode from incident O mode [0-1]';
	    info.roe        = 'refl of O mode from incident X mode [0-1]';
	    info.albajar    = 'if =1, normalization of total power on Albajar scaling law';
            info.reflex     = 'edge reflection coefficient of the cyclotronic radiation';
     else
	    info.ree        = 'refl of X mode from incident X mode [0-1]';
	    info.reo        = 'refl of X mode from incident O mode [0-1]';
	    info.roo        = 'refl of O mode from incident O mode [0-1]';
	    info.roe        = 'refl of O mode from incident X mode [0-1]';
	    info.albajar    = 'if =1, normalization of total power on Albajar scaling law';
            info.reflex     = 'edge reflection coefficient of the cyclotronic radiation';
     end

	  interface.ts           = '';                            % nom de la fonction d'interfacage avec les donnees TS
	  interface.jet          = '';                            % nom de la fonction d'interfacage avec les donnees Jet

	  sortie.valeur          = valeur;
	  sortie.type            = type;
	  sortie.borne           = borne;
	  sortie.defaut          = defaut;
	  sortie.info            = info;
	  sortie.interface       = interface;
	  sortie.description     = 'Compute the cycltronic radiation with CYTRAN code';   % description (une ligne) de la fonction
	  sortie.help            = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	  sortie.gui             = '';                            % nom de l'interface graphique specifique si elle existe
	  sortie.controle        = '';                            % nom de la fonction de controle des valeurs si elle existe

     pcyclo = sortie;	
     return
end

% selon le nombre de points
ne   = prof.ne;
te   = prof.te;
bmod = sqrt(equi.b2);
vpr  = equi.vpr;
x    = gene.x;
% mise au pas
xn   = linspace(0,1,101);
ne   = interp1(x,ne,xn,'spline'); 
te   = interp1(x,te,xn,'spline'); 
bmod = interp1(x,bmod,xn,'spline'); 
vpr  = interp1(x,vpr,xn,'spline'); 



% preparation des donnees de cytran
ne   = (ne(1:end-1) + ne(2:end)) ./ 2;
te   = (te(1:end-1) + te(2:end)) ./ 2e3;
bmod = (bmod(1:end-1) + bmod(2:end)) ./ 2;
xc   = (xn(1:end-1) + xn(2:end)) ./ 2;
vol  =  equi.rhomax .* diff(cumtrapz(xn,vpr,2));
surf = vpr;

% consigne
cons(1) = parametre.ree;
cons(2) = parametre.reo;
cons(3) = parametre.roo;
cons(4) = parametre.roe;

% appel de cytran
%save contexte_cytran_SIZE
%disp('avant cytran_1')
[pcyclo,ptot] = cytran_77(te,ne,surf,vol,bmod,cons);
%disp('apres cytran_1')
pcyclo(isnan(pcyclo))=0;
if isnan(ptot)
  ptot=0.01;
end

indnok = find(xc>0.95);

if ~isempty(indnok)
pcyclo(indnok) = [];
	
end
pcyclo(100)=0;



% changement d'unite
pcyclo = pcyclo .* 1e6;
% puissance totale
%ptot   = ptot(1) .* 1e6;
ptot   = ptot(1);

% changement de grilles
pcyclo = pchip(xc,pcyclo,gene.x);

% conservation de la puissance + signe cronos
%pcyclo =  pcyclo .* ptot ./ equi.rhomax ./ trapz(gene.x,equi.vpr .* pcyclo,2);
%ptot ./ equi.rhomax ./ trapz(gene.x,equi.vpr .* pcyclo,2)
% renormalisation sur le fit de Albajar
if parametre.albajar == 1
   [pvoid,palbajar] = zcyclo(gene.x,prof.te,prof.ne,equi.rmoy,equi.a,equi.e,equi.F,equi.b2, ...
                              equi.vpr,equi.spr,equi.rhomax,parametre.reflex,phys);

   %disp(-palbajar ./ ptot)
   pcyclo = - pcyclo ./ ptot .* palbajar;
end



[pvoid,palbajar] = zcyclo(gene.x,prof.te,prof.ne,equi.rmoy,equi.a,equi.e,equi.F,equi.b2, ...
                              equi.vpr,equi.spr,equi.rhomax,parametre.reflex,phys);
%figure(20);plot(gene.x,pcyclo,'b')
%drawnow
