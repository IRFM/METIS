% ZBGBS modele de coefficients de transport Bohm/Gyro Bohm 
%--------------------------------------------------------------
% fichier zjetto.m ->  zjetto
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule les coefficients de transport du modele Bohm/Gyro Bohm
% Elle inclue un effet de bord pour le mode H et un effet de scisaillemnt inverse
% pour l'amelioration du confinement.
% Elle permet de rajouter au transport de la matiere, le pinch de Ware  et une vitesse
% de convection anormale en 1/Ne^2.
% 
% syntaxe  :
%  
%      coef = zjetto(parametre,coef,geo,equi,impur,prof,neo,phys,compo,nbrho,nbg,x,source);
%
% entree :
%
%      parametre       =    parametre propre a la fonction (param.cons.fci)
%      coef            =    structure pour les coefficients ;
%      geo             =    geometrie du plasma (data.geo)
%      equi            =    donnees de l'equilibre plasma (data.equi)
%      impur           =    sous structure impuretes
%      prof            =    profils des donnees calculees par le code (data.prof) 
%      neo             =    donnees neoclassiques (data.neo)
%      phys            =    constantes physiques (param.phys)
%      compo           =    composition du plasma: charge et masse des atomes ( param.compo)
%      nbrho           =    nombre de points radiaux
%      nbg             =    nombre de gaz
%      x               =    coordonnees normalisee (gene.x)
%      source          =    structure source pour les modeles cale sur les lois d'echelles
%
% 
% sortie :
% 
%     coef             =  structure coef modifiee
% 
%
% fonction ecrite par J-F Artaud et F. Imbeaux , poste 46-78
% version 2.2 du 02/06/2004
% 
% liste des modifications : 
%
%  * 20/09/2001 -> remplacement de prof.shear par equi.shear
%  * 08/11/2001 -> ajout de l'entree source
%  * 02/04/2003 -> reecriture transport densite selon la reference Garzotti et al, Particle transport ... to appear in Nuclear Fusion 2003
%  * 02/06/2004 -> ajout des coefficient pour la rotation
%--------------------------------------------------------------
%
function coef = zbgbs(parametre,coef,geo,equi,impur,prof,neo,phys,compo,nbrho,nbg,x,source)

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	langue                  = getappdata(0,'langue_cronos');
	
	valeur.inverse      = 1;                        % effet scisaillement inverse (oui -> 1, non -> 0) [1]
	valeur.modeh        = 1;                        % effet mode H (oui -> 1, non -> 0) [0]
	valeur.seuil        = 0.05;                     % seuil pour le scisaillement inverse [0.05]
	valeur.facteur      = 20;                       % facteur de l'effet du scisaillement inverse [1]
	valeur.aeb          = 8e-5;                     % coefficient du terme bohm de Xie [8e-5]
	valeur.aegb         = 7e-2;                     % coefficient du terme gyro bohm de Xie [7e-2]
	valeur.aib          = 2 .* 8e-5;                % coefficient du terme bohm de Xii [1.6e-4]
	valeur.aigb         = 0.25 .* 7e-2;             % coefficient du terme gyro bohm de Xii [1.75e-2]
	valeur.typeconvection = 'Grad q/q';             % type de convection de la matiere 
	valeur.convection   = 0.5;                      % vitesse de convection de la matiere  [0.5]
	valeur.ximax        = 20;                       % valeur maximale des Xie et Xii
	
	type.inverse        = 'integer';                % type entier
	type.modeh          = 'integer';                % type entier
	type.seuil          = 'float';                  % type reel
	type.facteur        = 'float';                  % type reel
	type.aeb            = 'float';                  % type reel
	type.aegb           = 'float';                  % type reel
	type.aib            = 'float';                  % type reel
	type.aigb           = 'float';                  % type reel
	type.typeconvection = 'string';                 % type string
	type.convection     = 'float';                  % type reel
	type.ximax          = 'float';                  % type reel
	
	borne.inverse       = {0,1};                    % valeurs possibles
	borne.modeh         = {0,1};                    % valeurs possibles
	borne.seuil         = [0,inf];                  % valeurs possibles
	borne.facteur       = [0,inf];                  % valeurs possibles
	borne.aeb           = [0,inf];                  % valeurs possibles
	borne.aegb          = [0,inf];                  % valeurs possibles
	borne.aib           = [0,inf];                  % valeurs possibles
	borne.aigb          = [0,inf];                  % valeurs possibles
	borne.typeconvection= {'Grad q/q','Grad Te/Te'};      % valeurs possibles
	borne.convection    = [0,inf];                  % valeurs possibles
	borne.ximax         = [10,inf];                 % valeurs possibles
	
	defaut.inverse      = 1;                        % valeurs par defaut
	defaut.modeh        = 1;                        % valeurs par defaut
	defaut.seuil        = 0.05;                     % valeurs par defaut
	defaut.facteur      = 20;                       % valeurs par defaut
	defaut.aeb          = 8e-5;                     % valeurs par defaut
	defaut.aegb         = 7e-2;                     % valeurs par defaut
	defaut.aib          = 2 .* 8e-5;                % valeurs par defaut
	defaut.aigb         = 0.25 .* 7e-2;             % valeurs par defaut
	defaut.typeconvection = 'Grad q/q';                        % valeurs par defaut
	defaut.convection   = 0.5;                        % valeurs par defaut
	defaut.ximax        = 20;                        % valeurs par defaut

   if strcmp(langue,'francais')	
	  info.inverse      = 'reduction du transport si cisaillement magnetique negatif (oui -> 1, non -> 0) [1]';   %informations
	  info.modeh        = 'facteur mode H (facteur Grad Te/Te|bord dans terme Bohm) (oui -> 1, non -> 0) [1]';
	  info.seuil        = 'seuil de l''effet du cisaillement negatif [0.05]';
	  info.facteur      = 'facteur de l''effet du cisaillement negatif [20]';
	  info.aeb          = 'coefficient du terme Bohm de Xie [8e-5]';
	  info.aegb         = 'coefficient du terme gyro-Bohm de Xie [7e-2]';
	  info.aib          = 'coefficient du terme Bohm de Xii [1.6e-4]';
	  info.aigb         = 'coefficient du terme gyro-Bohm de Xii [1.75e-2]';
	  info.typeconvection = 'expression de la vitesse de convection de la matiere [Grad q/q]';
	  info.convection   = 'coefficient du terme de vitesse de convection de la matiere [0.5]';
	  info.ximax        = 'valeur maximale des Xie et Xii';;
	end
   if strcmp(langue,'anglais')	
	  info.inverse      = 'transport reduction for negative magnetic shear (yes -> 1, no -> 0) [1]';   %informations
	  info.modeh        = 'H mode factor (Grad Te/Te|edge factor in Bohm term)(yes -> 1, no -> 0) [1]';
	  info.seuil        = 'shear threshold [0.05]';
	  info.facteur      = 'shear factor [20]';
	  info.aeb          = 'Bohm coefficient for Xie [8e-5]';
	  info.aegb         = 'Gyro-Bohm coefficient for  Xie [7e-2]';
	  info.aib          = 'Bohm coefficient for Xii [1.6e-4]';
	  info.aigb         = 'Gyro-Bohm coefficient for Xii [1.75e-2]';
	  info.typeconvection = 'expression convective speed for density [Grad q/q]';
	  info.convection   = 'convective speed coefficient for the density';
	  info.ximax        = 'max values of Xie and Xii';;
	end
   
	interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	coef.valeur     = valeur;
	coef.type       = type;
	coef.borne      = borne;
	coef.defaut     = defaut;
	coef.info       = info;
	coef.interface  = interface;
	
	coef.description = 'Modele de coefficients de transport Bohm/Gyro Bohm ';   % description (une ligne) de la fonction
	
	coef.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	coef.gui      ='';                             % nom de l'interface graphique specifique si elle existe
	coef.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	return
end

% les division par 0 sont traitees
warning off

% calcul de rho_star (de l'ion majoritaire)
rho_star = 1.02e-4 .* sqrt(compo.a(1) .* prof.te) ./ (compo.z(1) .* equi.rhomax .* prof.bphi);

% gradient mode H
if parametre.modeh  == 1
    tebord  =prof.te(end);
    ind = min(find(x>=0.8));
    gh  = abs(prof.te(ind) - tebord)./tebord;
else
    gh  = abs(800 -100)./100;
end
if ~isfinite(gh)| (gh <0)
  gh = 1;
end

% effet d scisaillement inverse :
if parametre.inverse  == 1
	frs = 1 ./ (1 + exp(parametre.facteur .* (parametre.seuil - prof.shear)));
else
	frs=ones(1,nbrho);
end

% Xie bohm
terme1 = abs(prof.ped1) ./ prof.pe;
xib  = frs .* prof.te ./ prof.bphi .* terme1 .* prof.q .^ 2 .* gh;

% Xie gyro bohm
ted1    = abs(pdederive(x,prof.te,0,2,2,1) ./ equi.rhomax );
xigb    = rho_star  ./ prof.bphi .* ted1 .* (ted1 > 0);

% les coefficients de transport
xie  = parametre.aegb .* xigb + parametre.aeb .* xib;
xie(xie > parametre.ximax) = parametre.ximax;
xii  = parametre.aigb .* xigb + parametre.aib .* xib;
xii(xii > parametre.ximax) = parametre.ximax;

% pour eviter les gradients de chi trop forts au niveau de la barriere
if parametre.inverse  == 1
	aa = find(xie(1:50)>5);
	if ~isempty(aa)
		xie(aa) = 5;
	end
	bb = find(xii(1:50)>6);
	if ~isempty(bb)
		xii(bb) = 6;
	end
end

d    = linspace(1,0.3,nbrho) .* xie .* xii ./ (xie + xii);
ind = find(~isfinite(d));
if ~isempty(ind)
	d(ind) = zeros(1,length(ind));
end


kie = xie  .* prof.ne;
kii = xii .* prof.ni;

% Convection densite
if strcmp(parametre.typeconvection,'Grad q/q')
   qd1   =  pdederive(x,prof.q,0,2,2,1) ./ equi.rhomax;
   vn    =  - parametre.convection .* d .*qd1 ./prof.q;     % turbulence equipartition
elseif strcmp(parametre.typeconvection,'Grad Te/Te')
   vn    =  parametre.convection .* d ./prof.lte;     % thermodiffusion
end

ind = find(~isfinite(vn));
if ~isempty(ind)
	vn(ind)    = zeros(1,length(ind));
end



% remplissage de la structure coef
coef.ee = kie;
coef.ii = kii;
coef.nn = d;
coef.vn = vn;
coef.rot = xii;
coef.rotv = 0 .* vn;

% fin de la conction
warning on

%return

if 1 > 2
h = findobj(0,'type','figure','tag','zbgbs');
if isempty(h)
    h = figure('tag','zbgbs');
    hc=uicontrol(h,'style','radio','tag','onpause','value',0);
else
    figure(h);
    hc = findobj(h,'style','radio','tag','onpause');
end
set(h,'color',[1 1 1])
subplot(3,2,1)
plot(x,coef.ee ./ prof.ne,'r',x,coef.ii ./ prof.ni,'b',x,neo.coef.ee./prof.ne,'m',x,neo.coef.ii./prof.ni,'c');
ylabel('Chi')
subplot(3,2,2)
plot(x,prof.q,'r',x,prof.shear,'b');
ylabel(' q & s');
subplot(3,2,3)
plot(x,prof.ne,'r',x,prof.ni,'b');
ylabel('density')
set(gca,'ylim',[0,inf]);
subplot(3,2,4)
plot(x,prof.te,'r',x,prof.ti,'b');
set(gca,'ylim',[0,max(30,max(max(prof.ti),max(prof.te)))]);
ylabel('Te & Ti');
subplot(3,2,5)
plot(x,prof.jmoy ./ 1e6,x,source.jboot ./ 1e6);
ylabel('Jmoy & Jboot')
end
