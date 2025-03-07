% ZBGBS modele de coefficients de transport Bohm/Gyro Bohm
% avec les coefficients optimises pour TS [Erba et al, Nucl. Fus. 38 (1998) 1013]
%--------------------------------------------------------------
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule les coefficients de transport du modele Bohm/Gyro Bohm
% Elle inclue un effet de bord pour le mode H et un effet de cisaillemnt inverse
% pour l'amelioration du confinement.
% Elle permet de rajouter au transport de la matiere, le pinch de Ware  et une vitesse
% de convection anormale en 1/Ne^2.
% 
% syntaxe  :
%  
%      coef = zbgbs_TS(parametre,coef,geo,equi,impur,prof,neo,phys,compo,nbrho,nbg,x,source);
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
% fonction ecrite par J-F Artaud et F. Imbeaux,
% version 2.0, du 21/11/2002.
% 
% liste des modifications : 
%
%  * 20/09/2001 -> remplacement de prof.shear par equi.shear
%  * 08/11/2001 -> ajout de l'entree source
%  * 22/11/2002 -> optimisation des coefficients par defaut pour TS (valeurs donnees dans Erba et al, Nucl. Fus. 38 (1998) 1013)
%  * 22/11/2002 -> suppression du "mode H" (revenait a multiplier par 7 aeb et aib pour TS)
%  * 22/11/2002 -> suppression du terme de pinch pour la chaleur electronique (systematiquement pris en compte avant, mais negligeable ~2e-4 m/s)
%  * 25/11/2002 -> remplacement de la thermodiffusion par une diffusion des particules en grad q/q 
%--------------------------------------------------------------
%
function coef = zbgbs_TS(parametre,coef,geo,equi,impur,prof,neo,phys,compo,nbrho,nbg,x,source)

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	
	valeur.inverse_TS      = 1;                        % effet cisaillement magnetique (oui -> 1, non -> 0) [1]
	valeur.seuil_TS        = 0.05;                     % seuil pour le cisaillement magnetique [0.05]
	valeur.facteur_TS      = 20;                       % intensite de l'effet du scisaillement inverse [1]
	valeur.aeb_TS          = 2.5e-4;                   % coefficient du terme bohm de Xie [2.5e-4]
	valeur.aegb_TS         = 3.5e-2;                   % coefficient du terme gyro bohm de Xie [3.5e-2]
	valeur.aib_TS          = 2 .* valeur.aeb_TS;          % coefficient du terme bohm de Xii [5e-4]
	valeur.aigb_TS         = 3.5e-2;                   % coefficient du terme gyro bohm de Xii [3.5e-2]
	valeur.d_TS            = 0.1;                      % coefficient multiplicateur de D [0.1]
	valeur.convection_TS   = 0.6;                      % vitesse de convection de la matiere  [0.6]
	valeur.ximax_TS        = 20;                       % valeur maximale des Xie et Xii
	
	type.inverse_TS        = 'integer';                % type entier
	type.seuil_TS          = 'float';                  % type reel
	type.facteur_TS        = 'float';                  % type reel
	type.aeb_TS            = 'float';                  % type reel
	type.aegb_TS           = 'float';                  % type reel
	type.aib_TS            = 'float';                  % type reel
	type.aigb_TS           = 'float';                  % type reel
	type.d_TS              = 'float';                  % type reel
	type.convection_TS     = 'float';                  % type reel
	type.ximax_TS          = 'float';                  % type reel
	
	borne.inverse_TS       = {0,1};                    % valeurs possibles
	borne.seuil_TS         = [0,inf];                  % valeurs possibles
	borne.facteur_TS       = [0,inf];                  % valeurs possibles
	borne.aeb_TS           = [0,inf];                  % valeurs possibles
	borne.aegb_TS          = [0,inf];                  % valeurs possibles
	borne.aib_TS           = [0,inf];                  % valeurs possibles
	borne.aigb_TS          = [0,inf];                  % valeurs possibles
	borne.d_TS             = [0,inf];                  % valeurs possibles
	borne.convection_TS    = [0,inf];                  % valeurs possibles
	borne.ximax_TS         = [10,inf];                 % valeurs possibles
	
	defaut.inverse_TS      = 1;                        % valeurs par defaut
	defaut.seuil_TS        = 0.05;                     % valeurs par defaut
	defaut.facteur_TS      = 20;                       % valeurs par defaut
	defaut.aeb_TS          = 2.5e-4;                     % valeurs par defaut
	defaut.aegb_TS         = 3.5e-2;                     % valeurs par defaut
	defaut.aib_TS          = 2 .* defaut.aeb_TS;                % valeurs par defaut
	defaut.aigb_TS         = 3.5e-2;             % valeurs par defaut
	defaut.d_TS            = 0.1;                        % valeurs par defaut
	defaut.convection_TS   = 0.6;                        % valeurs par defaut
	defaut.ximax_TS        = 20;                        % valeurs par defaut
	
	info.inverse_TS      = 'reduction du transport si cisaillement magnetique < seuil (oui -> 1, non -> 0) [1]';   %informations
	info.seuil_TS        = 'seuil de la reduction du transport liee au cisaillement magnetique [0.05]';
	info.facteur_TS      = 'intensite de la reduction du transport liee au cisaillement magnetique [20]';
	info.aeb_TS          = 'coefficient du terme bohm de Xie [2.5e-4]';
	info.aegb_TS         = 'coefficient du terme gyro bohm de Xie [3.5e-2]';
	info.aib_TS          = 'coefficient du terme bohm de Xii [5e-4]';
	info.aigb_TS         = 'coefficient du terme gyro bohm de Xii [3.5e-2]';
	info.d_TS            = 'coefficient multiplicateur de D [0.1]';
	info.convection_TS   = 'coefficient du terme de vitesse de convection de la matiere';
	info.ximax_TS        = 'valeur maximale des Xie et Xii';;
	
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

% effet d scisaillement inverse :
if parametre.inverse_TS  == 1
	frs = 1 ./ (1 + exp(parametre.facteur_TS .* (parametre.seuil_TS - prof.shear)));
   %aa = find(prof.shear>=1);
   %frs(aa)=ones(size(aa));
else
	frs=ones(1,nbrho);
end

% Xie bohm
terme1 = abs(prof.ped1) ./ prof.pe;
xib  = frs .* prof.te ./ prof.bphi .* terme1 .* prof.q .^ 2;

% Xie gyro bohm
ted1    = abs(pdederive(x,prof.te,0,2,2,1) ./ equi.rhomax );
xigb    = rho_star  ./ prof.bphi .* ted1 .* (ted1 > 0);

% les coefficients de transport
xie  = parametre.aegb_TS .* xigb + parametre.aeb_TS .* xib;
xie(xie > parametre.ximax_TS) = parametre.ximax_TS;
xii  = parametre.aigb_TS .* xigb + parametre.aib_TS .* xib;
xii(xii > parametre.ximax_TS) = parametre.ximax_TS;
d    = parametre.d_TS .* xie .* xii ./ (xie + xii);
ind = find(~isfinite(d));
if ~isempty(ind)
	d(ind) = zeros(1,length(ind));
end

%xie = zregular(x,xie);
%xii = zregular(x,xii);
%d   = zregular(x,d);


kie = xie  .* prof.ne;
kii = xii .* prof.ni;

% le pinch suplementaire :
% vn    =  parametre.convection_TS .* d ./prof.lte;
vn    =  parametre.convection_TS .* d .*(prof.shear./(x.*equi.rhomax+1e-3));
vn    =  zregular(x,vn);

ind = find(~isfinite(vn));
if ~isempty(ind)
	vn(ind)    = zeros(1,length(ind));
end

% remplissage de la structure coef
coef.ee = kie;
coef.ii = kii;
coef.nn = d;
coef.vn = vn;

%disp([num2str(max(vn)),'   ',num2str(min(vn))])

% fin de la conction
warning on

