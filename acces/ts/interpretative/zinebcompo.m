% ZINEBCOMPO calcul de la puissance rayonnee et du rayonnement de freinage
%-------------------------------------------------------------------------
% fichier zinebcompo.m ->  zinebcompo, zcompo
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule le rapport alphae, les profils d'impuretes correspondant,
% et les densite de puissance rayonnees. Elle peut aussi calculer un profil de zeff
% (de la forme Ne^a, a= -0.4 par defaut) a partir du zeff lineique moyen. 
%  
% syntaxe  :
%  
%     [impur,memoire]=zinebcompo(cons,geo,equi,bord,zeffm,prof,neo,phys,compo,gene,memoire,source,impurmem,fonctionneo,consneo,coef,neomem);
%
% entrees :
%
%     cons    =  struture des consignes de calcul de la fonction
%     geo     =  datak.geo 
%     equi    =  datak.equi
%     bord    =  datak.mur.fluxgazin (pas utiliser dans cette fonction, reserve pour d'autres modules)
%     zeffm   =  datak.cons.zeffm (zeff lineique moyen utiliser pour calculer le profil de zeff)
%     nhnd    =  datak.cons.nhnd  (rapport nh/nD si disponible)
%     prof    =  datak.prof
%     neo     =  datak.neo  (pas utiliser dans cette fonction, reserve pour d'autres modules)
%     phys    =  param.phys
%     compo   =  param.compo (numero atomique et nombre de masses des gaz)
%     gene    =  param.gene
%     memoire =  datak.memoire.impur (valeur de reference pour le dernier calcul complet, 
%                pas utiliser dans cette fonction, reserve pour d'autres modules)
%     source  =  datak.source.n (source d'impurte bord +fusion+idn,pas utiliser 
%                dans cette fonction, reserve pour d'autres modules)
%     impurmem = datak.impur au temps precedent (pas utiliser dans cette fonction)
%     fonctionneo = nom de la fonction de calcul neoclassique
%     consneo     = parametres associes.
%     coef        = coefficients de transport .
%     neomem      = structure memoire de neo
%
%
% sorties :
% 
%     impur   = structure impur :
%                   impur.prad   = densite  de puissance rayonnee (raies+recombinaison)
%                   impur.brem   = densite de  de puissance de rayonnement de frainage (valeur indicative)
%                   impur.cyclo  = densite  de puissance de rayonnement cyclotron
%                   impur.zeff   = profil de zeff calculer dans la fonction
%                   impur.impur  = profils de densite d'impuretes par espece
%                   impur.ae     = profil de facteur reliant ne a ni(eff)
%                   impur.conv   = nombre de boucles necessaires a la convergence du code de transport des impuretes avec le transport neoclassique
%                   impur.fail   = non nul si le  code de transport des impuretes n'a pas converge ou que le resultat n'est pas consitent avec le transport neoclassique
%                   impur.neofail= non nul si le code neoclassique appele en meme temps que le code de transport des impuretes n'a pas converge
%
%     memoire = datak.memoire.impur (valeur de reference pour le dernier calcul complet, 
%               pas utiliser dans cette fonction, reserve pour d'autres modules)
% 
% parametres :
% 
%     cons.zeff     =  0 => zeff donnee (prof.zeff)
%                      1 => zeff calcule par la fonction en utilisant la consigne  zeffm
%                      2 => zeff calcule par la fonction en utilisant la loi d'echelle TS  
%                           (dans ce cas zeffm doit etre donne pour nbar = 1e20 m-3)     
%     cons.exposant =  exposant tel que zeff = cte * Ne ^exposant
%     cons.cmin1    =  rapport de la densite du 1er minoriatire sur la densite de l'espece principale
%     cons.cmin2    =  rapport de la densite du 2ieme minoriatire sur la densite de l'espece principale
%     cons.rimp     =  rapport de la densite de la 2ieme impurete sur la densite de la premiere impurete
%
% remarque : la puissance rayonnee correcte est donnees par prad + brem .
%      
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 31/01/2005.
% 
% 
% liste des modifications : 
%
%  * 07/09/2001 -> ajout de l'utilisation de la consigne nhnd pour le control du rapport isotopique D/T
%  * 02/10/2001 -> ajout de la structure memoire en sortie
%  * 16/04/2002 -> changement des securite sur zeffm et zeff
%  * 17/07/2002 -> correction cas ou zeffm n'est pas defini
%  * 28/08/2002 -> ajout de la normalistion du rayonnement
%  * 14/11/2002 -> ajout de l'entree impurmem pour le code mariobea
%  * 11/12/2002 -> interface anglais
%  * 06/03/2003 -> correction d'une division par zero de "pc" au centre
%  * 10/06/2003 -> securite electroneutralite exacte (1e-16)
%  * 10/09/2003 -> utilisation du 0D de Albajar poour le cacul de cyclo
%  * 11/09/2003 -> ajout des profils d'impuretes neoclassique
%  * 11/09/2003 -> ajout de la gestion de la fraction ntnd
%  * 23/02/2004 -> lignes autour de 276 enlevees (rajoutaient des NaN dans les concentration de minoritaire en absence de nhnd ?)
%  * 10/01/2005 -> suppression du calcul de la puissance cyclotronique
%  * 11/01/2005 -> ajout du nom de la fonction neo et des parametres associees
%  * 11/01/2005 -> ajout des coefficients de transports en entree
%  * 14/01/2005 -> ajout de impur.fail, impur.conv et impur.neofail
%  * 24/01/2005 -> ajout de memoire.neo en entree
%  * 31/01/2005 -> ajout de impur.pradsol
%  * 31/01/2005 -> ajout de impur.alpha
%  * 05/12/2007 -> Correction des commentaires version anglaise
%--------------------------------------------------------------
%
function [impur,memoire] = zinebcompo(cons,geo,equi,bord,zeffm,nhnd,prof,neo,phys,compo,gene,memoire,source,impurmem,fonctionneo,consneo,coef,neomem)

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	langue                  = getappdata(0,'langue_cronos');

	valeur.zeff          = 0;           % zeff donnee (prof.zeff) (0-> zeff donne | 1-> zeff calcule (Ne^exposant)| 2 -> zeffm dependant de nbar et zeff calcule comme en 1)
	valeur.exposant      = -0.4;        % exposant tel que zeff = cte * Ne ^exposant
	valeur.cmin1         = 0.1;         % rapport de la densite du 1er minoriatire sur la densite de l'espece principale (He3/D)
	valeur.cmin2         = 0;           % rapport de la densite du 2ieme minoriatire sur la densite de l'espece principale
	valeur.rimp          = 0.1;         % rapport de la densite de la 2ieme impurete sur la densite de la premiere impurete
	valeur.neoimp        = 0;           % profil impurete : 0- > zeff + parametre,1 - > neoclassique + zeff au centre
	%valeur.reflex        = 0.8;         % coefficient de reflexion du rayonnement cyclo par les parois
	valeur.nhndonoff     = 1;           % 0 -> n'utilise pas la consigne nhnd et nhnt, 1 -> utilise la consigne nhnd (pour D/H et D/T)
	valeur.norme         = 1;           % 0 -> pas de normalisation, 1 -> normalisation de la puissance rayonnee 
	valeur.crad          = 1;           % position du maximum de puissance rayonnee due a la reconbinaison (en x)
	valeur.lrad          = 0.1;         % largeur de la source de puissance rayonnee due a la reconbinaison (en fraction de rhomax)
	valeur.frad          = 1;           % facteur correctif a appliquer a la  puissance rayonnee due a la reconbinaison 

	type.zeff            = 'logical';   % type logique
	type.exposant        = 'float';     % type reel
	type.cmin1           = 'float';     % type reel
	type.cmin2           = 'float';     % type reel
	type.rimp            = 'float';     % type reel
	type.neoimp          = 'integer';   % type entier
	%type.reflex          = 'float';     % type reel
	type.nhndonoff       = 'logical';   % type logique
	type.norme           = 'logical';   % type logique
	type.crad            = 'float';     % type reel
	type.lrad            = 'float';     % type reel
	type.frad            = 'float';     % type reel

	borne.zeff           = {0,1,2};       % valeurs possible
	borne.exposant       = [0,-2];      % valeurs possible 
	borne.cmin1          = [0,1];       % valeurs possible 
	borne.cmin2          = [0,1];       % valeurs possible 
	borne.rimp           = [0,1];       % valeurs possible 
	borne.neoimp         = {0,1,2};      % valeurs possible
	%borne.reflex         = [0,1];       % valeurs possible 
	borne.nhndonoff      = {0,1};       % valeurs possible
	borne.norme          = {0,1};       % valeurs possible
	borne.crad           = [0,1];       % valeurs possible 
	borne.lrad           = [0,1];       % valeurs possible 
	borne.frad           = [0.1,10];     % valeurs possible 

	defaut.zeff          = 0;           % valeurs par defaut
	defaut.exposant      = -0.4;        % valeurs par defaut 
	defaut.cmin1         = 0.1;         % valeurs par defaut 
	defaut.cmin2         = 0;           % valeurs par defaut 
	defaut.rimp          = 0.1;         % valeurs par defaut
	defaut.neoimp        = 0;           % valeurs par defaut
	%defaut.reflex        = 0.8;         % valeurs par defaut
	defaut.nhndonoff     = 1;           % valeurs par defaut
	defaut.norme         = 1;           % valeurs par defaut
	defaut.crad          = 1;           % valeurs par defaut
	defaut.lrad          = 0.1;         % valeurs par defaut
	defaut.frad          = 1;           % valeurs par defaut
   if strcmp(langue,'francais')
	  info.zeff           = '0-> zeff donne | 1-> zeff calcule (Ne^exposant) normalise a zeffm| 2 -> zeffm dependant de nbar et zeff calcule comme en 1';       % informations
	  info.exposant       = 'exposant de la loi zeff = a * Ne ^exposant';
	  info.cmin1          = 'rapport de la densite du 1er minoritaire sur la densite de l''espece principale';
	  info.cmin2          = 'rapport de la densite du 2ieme minoritaire sur la densite de l''espece principale';
	  info.rimp           = 'rapport de la densite de la 2ieme impurete sur la densite de la premiere impurete';
	  info.neoimp         = 'profils impuretes : 0- > zeff + parametre (comme cronos < v2.2),1 - > neoclassique + zeff central'; 
	  %info.reflex         = 'coefficient de reflexion du rayonnement cyclo par les parois';
	  info.nhndonoff      = '0 -> n''utilise pas la consigne nhnd et ntnd, 1 -> utilise la consigne nhnd et ntnd';
	  info.norme          = '0 -> pas de normalisation, 1 -> normalisation de la puissance rayonnee  {1}';
	  info.crad           = 'position du maximum de puissance rayonnee due a la recombinaison (en x)';
	  info.lrad           = 'largeur de la source de puissance rayonnee due a la recombinaison (en fraction de rhomax)';
	  info.frad           = 'facteur correctif a appliquer a la  puissance rayonnee due a la reconbinaison'; 
  else
	  info.zeff           = '0-> zeff as profile input (data.prof.zeff) | 1-> zeff profile deduced from average zeff (data.cons.zeffm) with parametric radial shape (see exposant) | 2 -> average zeff deduced from scaling law + parametric radial shape (see exposant) ';       % informations
	  info.exposant       = 'parametric radial shape of zeff profile -> zeff = a * Ne ^ exposant';
	  info.cmin1          = 'density ratio of the first minority species over main species';
	  info.cmin2          = 'density ratio of the second minority species over main species';
	  info.rimp           = 'density ratio of the second impurity over the first one';
	  info.neoimp         = 'impurities profiles : 0- > zeff + parameters (like cronos < V2.2), 1 - > neoclassic + zeff(x=0) '; 
	  %info.reflex         = 'edge reflection coefficient of the cyclotronic radiation';
	  info.nhndonoff      = '0 -> use cmin1 and cmin2 to calculate nH, nD, nT; 1 -> use input nh/nd and nt/nd ratios (data.cons.nhnd) to calculate the ions species densities';
	  info.norme          = '1 -> normalisation on the radiative power {1}';
	  info.crad           = 'radial location of the radiative power due to CXR (in normalized radius)';
	  info.lrad           = 'width of the CXR radiative power (in normalized radius)';
	  info.frad           = 'CXR correctif factor to the radiative power'; 
   end
	interface.ts        = 'a faire ts';   % nom de la fonction d'interfacage avec les donnees TS
	interface.jet       = 'a faire jet';  % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'Calcul simple de alphae, du rayonnement et de la densite d''impurtes (+zeff)';   % description (une ligne) de la fonction
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	impur = sortie;
	return

end

% provisoire 
if ~isfield(cons,'norme')
   cons.norme = 0;
end
if ~isfield(cons,'neoimp')
   cons.neoimp = 0;
end
if ~isfield(cons,'crad') | ~isfield(cons,'lrad')| ~isfield(cons,'frad')
	cons.crad  = 1;
	cons.lrad  = 0.1;
	cons.frad  = 1;
end

% calcul de zmax
zmax   = compo.z(end - 1) .* (1 - cons.rimp) + compo.z(end) .* cons.rimp;
zmax   = min(max(compo.z) .* 0.9,zmax);  
zmax   = max(min(compo.z) .* 2.1 , zmax);  

% calcul du zeff	                         
if cons.zeff ==0
	% le profil de zeff est une donnee
	impur.zeff = prof.zeff;
else	
	% si cons.zeff == 2 => zeffm est calcule avec la loi d'echelle TS
	% utiliser pour l'injection de glacon
	if (cons.zeff == 2) | ~isfinite(zeffm)
		nbar = trapz(gene.x,prof.ne)./1e19;
		if compo.z(1) == 1
		     zeffm = 4.35 .* nbar .^ -0.4; 	
		elseif compo.z(1) ==2
		     zeffm = 10 .* nbar .^ -0.7;      
		end
		zeffm = min(max(zeffm,compo.z(1) .* 1.1),zmax - 0.2);    
	end
	
	
	% le profil de zeff est calcule 
	% zeff est une fonction de ne
	% test de compatibilite
	ind  = find(zeffm < compo.z);
	if isempty(ind)
		disp('Zeff references not comptible with the set of species !')
	end
	if cons.exposant == 0
		impur.zeff = zeffm .* ones(size(prof.ne));
	else	
		zeff = prof.ne .^ cons.exposant;
		zeff(end)  = 2 .* zeff(end-1) - zeff(end-2);
		g = trapz(gene.x,zeff,2);
		alpha = (zeffm - compo.z(1)) ./ (g - 1e21.^cons.exposant);
		gamma = (g .* compo.z(1)  - zeffm .* 1e21.^cons.exposant) ./(g - 1e21.^cons.exposant); 
		impur.zeff = alpha .* zeff + gamma;
	end
	
end

% securite zeff
ind = find(impur.zeff < min(compo.z));
if ~isempty(ind)
	impur.zeff(ind) = min(compo.z) .* ones(1,length(ind));
end
ind = find(impur.zeff > zmax);
if ~isempty(ind)
	impur.zeff(ind) = zmax .* ones(1,length(ind));
end

% utilisation de la consigne nhnd pour H et T
% real(nhnd) == nH/nD
% imag(nhnd) = nT/nD
if isfinite(nhnd)
  if  cons.nhndonoff == 1
     if isfinite(real(nhnd))
       ind = find((compo.z == 1) & (compo.a == 1));
       if ~isempty(ind)
         if ind == 2
           cons.cmin1 = real(nhnd);
         elseif ind == 3
           cons.cmin2 = real(nhnd);
         end
       end
     end
     if isfinite(imag(nhnd))
       ind = find((compo.z(2) == 1) & (compo.a(2) == 3));
       if ~isempty(ind)
         if ind == 2
           cons.cmin1 = imag(nhnd);
         elseif ind == 3
           cons.cmin2 = imag(nhnd);
         end
       end
     end
  end
%else     %%%%%%%%%%%%%%% ENLEVE le 23/02/2004 (F.I.), ca fait bugger zcompo de mettre des NaN dans les minoritaires !!!
% AUTANT LAISSER LES PARAMETRES RENTRES PAR L'UTILISATEUR !!!
%  disp('presence de NaN dans nhnd, zinebcompo')
%  cons.cmin1 = NaN;
%  cons.cmin2 = NaN;
end

% choix du mode des impuretees
% calcul de la composition du plasma
[ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(prof.ne,impur.zeff,cons.cmin1,cons.cmin2,cons.rimp, ...
                                         compo.z(1),compo.z(2),compo.z(3),compo.z(4),compo.z(5));

% securite (le plasma ne peut pas contenir que des impuretes)
while any(nion <= 0)
      zmax = zmax/1.1;
      impur.zeff(impur.zeff > (zmax)) = zmax;  
      [ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(prof.ne,impur.zeff,cons.cmin1,cons.cmin2,cons.rimp, ...
                                             compo.z(1),compo.z(2),compo.z(3),compo.z(4),compo.z(5));
end
% mise des resultats dans la structure impur
impur.ae    = ae .* (ae > 0.01) + 0.01 .* (ae <= 0.01);
nion        = nion  .* (nion > 0);
nmin1       = nmin1 .* (nmin1 > 0);
nmin2       = nmin2 .* (nmin2 > 0);
nimp2       = nimp2 .* (nimp2 > 0);
nimp1       = nimp1 .* (nimp1 > 0);
impur.impur = cat(3,nion,nmin1,nmin2,nimp1,nimp2);   

% choix du mode des impuretees
if cons.neoimp == 1
      % calcul des formes des profils d'impuretes neoclassiques
      impneo  = ((prof.ne ./prof.ne(1))'*ones(1,length(compo.a))) .^ (ones(size(prof.ne'))* compo.z) .* ...
                ((prof.ti ./prof.ti(1))'*ones(1,length(compo.a))) .^ (ones(size(prof.ne'))* (compo.z ./ 2 - 1));
     
       % normalisation 
      impneo   = impneo ./ (ones(size(prof.ne')) * max(impneo,[],1)); 
      
      % recherche des constantes
      indne      = 1;
      nmainneo   = impneo(indne,1) .* compo.z(1) + cons.cmin1 .* impneo(indne,2) .* compo.z(2) + cons.cmin2 .* impneo(indne,3) .* compo.z(3);
      nimpneo    = impneo(indne,4) .* compo.z(4) + cons.rimp  .* impneo(indne,5) .* compo.z(5);
      nmainneo2  = impneo(indne,1) .* compo.z(1) .^ 2 + cons.cmin1 .* impneo(indne,2) .* compo.z(2) .^ 2 + cons.cmin2 .* impneo(indne,3) .* compo.z(3) .^ 2;
      nimpneo2   = impneo(indne,4) .* compo.z(4) .^ 2 + cons.rimp  .* impneo(indne,5) .* compo.z(5) .^ 2;
     
      matneo     = [nmainneo2,nimpneo2;nmainneo,nimpneo];
      cneo       = matneo \ [prof.ne(indne).* impur.zeff(indne);prof.ne(indne)];

      nions    = cneo(1) .* impneo(:,1);
      nmin1s   = cneo(1) .* cons.cmin1 .* impneo(:,2);
      nmin2s   = cneo(1) .* cons.cmin2 .* impneo(:,3);
      nimp1s   = cneo(2) .* impneo(:,4);
      nimp2s   = cneo(2) .* cons.rimp  .* impneo(:,5);
      % mixte des solutions
      xi  = gene.x;
      xn  = 1 - gene.x;
      nion    = xn .* nions'  + xi .* nion;
      nmin1   = xn .* nmin1s' + xi .* nmin1;
      nmin2   = xn .* nmin2s' + xi .* nmin2;
      nimp1   = xn .* nimp1s' + xi .* nimp1;
      nimp2   = xn .* nimp2s' + xi .* nimp2;
      
      % recalcul
      ae      = (nion +  nmin1 + nmin2 + nimp1 + nimp2) ./ prof.ne;
      zeff    = (nion .* compo.z(1) .^ 2 + nmin1 .* compo.z(2) .^ 2 +  ...
                 nmin2 .* compo.z(3) .^ 2 + nimp1 .* compo.z(4) .^ 2 +  ...
                 nimp2 .* compo.z(5) .^ 2) ./ prof.ne;

      % mise des resultats dans la structure impur
      impur.ae    = ae .* (ae > 0.01) + 0.01 .* (ae <= 0.01);
      nion        = nion  .* (nion > 0);
      nmin1       = nmin1 .* (nmin1 > 0);
      nmin2       = nmin2 .* (nmin2 > 0);
      nimp2       = nimp2 .* (nimp2 > 0);
      nimp1       = nimp1 .* (nimp1 > 0);
      impur.impur = cat(3,nion,nmin1,nmin2,nimp1,nimp2);   
      impur.zeff  = zeff;
      
end

% calcul des puissances rayonnees
impur.prad = 0.* ae;
impur.brem = 0.* ae;

% boucle sur les especes
for k=1:gene.nbg
	[prad,brem,lz] = zradiation(compo.z(k),compo.a(k),prof.ne,impur.impur(1,:,k),prof.te);
	if ~all(prad == 0);
		impur.prad = impur.prad + prad -brem;
	end
	impur.brem = impur.brem + brem;
end

% normalisation du rayonnement
if (cons.norme == 1) && isfinite(cons.lrad) && isfinite(cons.crad) && (cons.lrad >0)
		nb20  = trapz(gene.x,prof.ne,2)./1e20;
                prtot = cons.frad .* 1e6 .* (trapz(gene.x,impur.zeff,2) - compo.z(1)) ./ 7 .* equi.vpr(end) .* nb20 .^ 2;
		prcor = prtot - equi.rhomax .* trapz(gene.x,equi.vpr .* (impur.brem + impur.prad),2);
		if prcor > 0
			recom = exp(-(gene.x-cons.crad) .^ 2 ./ cons.lrad .^2);
 			prcal = equi.rhomax .* trapz(gene.x,equi.vpr .* recom,2);
			recom = prcor ./ prcal .* recom;
			impur.prad = impur.prad + recom;
		end
end

% correction relativiste de Pbrem
% P.E. Stott, Plasma. Physics and Controlled Fusion 47 (2005) 1305-1338 ref [15]
tec2   = 511e3;% me c^2
xrel   = (1 + 2 .* prof.te ./ tec2) .* (1 + (2 ./ impur.zeff) .* (1 - (1 + prof.te ./ tec2) .^ (-1)));
xrel(imag(xrel) ~= 0) = 1;
impur.brem = xrel .* impur.brem;


% nouvelle variables pour le code de transport d'impuretes
impur.conv       = 0 ;    % nombre de boucles necessaires a la convergence du code de transport des impuretes avec le transport neoclassique
impur.fail       = 0;    % non nul si le  code de transport des impuretes n'a pas converge
impur.neofail    = 0;    % non nul si le code neoclassique appele en mem temps que le code de transport des impuretes n'a pas converge
impur.pradsol    = 2 .* trapz(gene.x,impur.prad .* equi.vpr,2);   % par convention pour le ploss
impur.alpha      = ones(size(gene.x));

% calcul de la puissance dissipee par le rayonnement cyclotronique
%[pcyclo,ptot] = zcyclo(gene.x,prof.te,prof.ne,equi.rmoy,equi.a,equi.e,equi.F,equi.b2, ...
%                              equi.vpr,equi.spr,equi.rhomax,cons.reflex,phys);
%impur.cyclo  = pcyclo;
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

