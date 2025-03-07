% ZNEUTRES calcule la distribution de neutres aux bord et les sources associees
%-------------------------------------------------------------------------------
% fichier zneutres.m ->  zneutres
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule (de maniere simple et fausse) la distribution de neutres au bord 
% du plasma. Elle calcule aussi les puits et sources associees. 
% 
% Nous attendons impatiament une fonction complete ... 
%  
% syntaxe  :
%  
%     [source,nn0,n0th,memoire] = zneutres(par,source,geo,equi,gazbord,prof,neo, ...
%                                  impur,bord,phys,compo,gene,memoire);
%
% entrees :
%
%     par        = struture des consignes de calcul de la fonction
%     source     = portotype a 0 de la structure source remplie par la fonction
%     geo        = geometrie de la derniere surface magnetique (datak.geo)
%     equi       = structure equilibre (datak.equi)
%     gazbord    = consigne d'injection de gaz au bord (datak.cons.c)
%     prof       = structure des profils (datak.prof)
%     neo        = structure des grandeurs neoclassiques (datak.neo)
%     impur      = structure des impuretes et du rayonnement (datak.impur)
%     bord       = structure decrivant le bord et le mur (datak.bord)
%     phys       = constante physiques (param.phys)
%     compo      = composition du plasma (param.compo)
%     gene       = parametre generiques (param.gene)
%     memoire    = datak.memoire.n0 (valeur de reference pour le dernier calcul complet, 
%                  pas utiliser dans cette fonction, reserve pour d'autres modules)
%     premiertemps = si 1, force le flux sortant a une valeur realiste
%     n0interp   = dernier calcul de n0 par le module pour interpolation
%     nn0interp  = dernier calcul de n.n0 par le module pour interpolation
%     n0thinterp = dernier calcul de n.n0 par le module pour interpolation
%
% sorties :
% 
%     source     = structure source  remplie, pour les neutres du bord
%     nn0        = profils de sources d'impuretes par espece du aux neutres du bord
%     n0th       = profil  de neutres thermiques au bord (1er echange de charges) 
%     memoire    = datak.memoire.n0 (valeur de reference pour le dernier calcul complet, 
%                  pas utiliser dans cette fonction, reserve pour d'autres modules)
% 
% parametres :
% 
%    cette fonction n'a pas de parametres
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 22/06/2004.
% 
% 
% liste des modifications : 
%
% * 02/10/2001 -> ajout de la structure memoire en sortie
% * 29/11/2002 -> ajout du mode premier temps
% * 22/09/2003 -> ajout du terme de friction pour la rotation
% * 22/06/2004 -> modification de l'energie sur le canal ionique
%
%--------------------------------------------------------------
%
function [source,nn0,n0th,memoire] = zneutres(par,source,geo,equi,gazbord,prof,neo, ...
                                      impur,bord,phys,compo,gene,memoire,premiertemps,n0interp,nn0interp,n0thinterp)

% fonction auto declarante                             
if nargin <=1 
	
	% cette fonction utilise aucun parametre
	valeur          = []; 
	type            = []; 
	borne           = [];      
	defaut          = [];           
	info            = [];
	interface.ts        = 'a faire ts';   
	interface.jet       = 'a faire jet';  
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'Calcul de l''equilibre';   % description (une ligne) de la fonction
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	sortie.resume ='';                           % nom de la fonction qui ecrit le resume du parametrage
	
	source = sortie;
	return

elseif nargin < 14 
   premiertemps = 0;
elseif isempty(premiertemps)
   premiertemps = 0;
end

if premiertemps == 1
  % calcul du flux min
  tau = equi.rmoy(end) .^ 2 /30;
  fluxmin = equi.rhomax .* trapz(gene.x,prof.ne .* equi.vpr) ./ tau;
  if ~isfinite(bord.fluxmur_c)
      bord.fluxmur_c = fluxmin;
  else
      bord.fluxmur_c = max(fluxmin,bord.fluxmur_c);
  end
end

% vpr modifier
vpr =equi.vpr;
vpr(1)=eps;
% volume plasma
vol = equi.rhomax .* trapz(gene.x,vpr,2);

% calcul des sections efficaces 
% 1- ionisation (par les electrons) 
te=[0.03 2 10 100 5e3 1e5];
sigmav=[1e-25 1e-17 7e-14 5e-14 1e-14 1e-15]; % m^3*s^-1
svi=exp(spline(log(te),log(sigmav),log(prof.te)));
% echange de charge
svcx = 7e-14;
% energie des neutres froids (eV)
if ~isfinite(bord.temp)
	bord.temp = 300;
elseif bord.temp <= 0
	bord.temp = 300;
end
e0   = bord.temp .* phys.k ./ phys.e; 
% les vitesse de penetration des neutres
v0i  = sqrt(2 .* phys.e .* max(13.6,prof.ti(end)) ./ phys.ua ./ compo.a(1));
v0th = sqrt(2 .* phys.e .* e0 ./ phys.ua ./ compo.a(1));

% calcul de la variable d'atenuation pour la ionisation (cf. Wesson p432 9.5)
sic = equi.rhomax .* cumtrapz(gene.x,prof.ne .* svi ./ v0i,2);
sic = 1  - (sic(end) -sic);
lic = 1 - gene.x(max(find(sic <=0)));
if isempty(lic)
	lic =0.2;
end

sif = equi.rhomax .* cumtrapz(gene.x,prof.ne .* svi ./ v0th,2);
sif = 1  - (sif(end) -sif);
lif  = 1 - gene.x(max(find(sif <=0)));
if isempty(lif)
	lif =0.2;
end

% calcul de l'attenuation pour l'echange de charge (cf. Wesson p 432)
scxc = equi.rhomax .* cumtrapz(gene.x,prof.ne .* svcx ./ v0i,2);
scxc = 1  - (scxc(end) -scxc);
lcxc  = 1 - gene.x(max(find(scxc <=0)));
if isempty(lcxc)
	lcxc =0.2;
end

scxf = equi.rhomax .* cumtrapz(gene.x,prof.ne .* svcx ./ v0th,2);
scxf = 1  - (scxf(end) -scxf);
lcxf  = 1 - gene.x(max(find(scxf <=0)));
if isempty(lcxf)
	lcxf =0.2;
end

% source associee + normalisation a 1
sc    = exp(- sqrt(1-gene.x.^2)./lcxc) .* svi;
intc  = equi.rhomax .*  trapz(gene.x,sc .* vpr);
sc    = sc ./ intc;

sf    = exp(- sqrt(1-gene.x.^2)./(lif + 5 .*lcxf)) .* svi;
%sf(end)  = sf(end-1);
%intf  = equi.rhomax .*  trapz(gene.x(1:(end-1)),sf(1:(end-1)) .* vpr(1:(end-1)));
intf  = equi.rhomax .*  trapz(gene.x,sf .* vpr);
sf    = sf  ./ intf ;

% source de matiere en equivalent electrons
% par default caclule pour 1 electron par seconde
source.ne      = sc .* max(bord.fluxmur_c,0.5) + sf .* max(bord.fluxmur_f,0.5);   % atomes*m^-3*s^-1
source.ne(1)   = 0;
% puit sur les electrons (25 eV par ionisation)
source.el      = - source.ne .* 25 .* phys.e;  % W*m^-3

% source sur les ions (W/m^3)
source.ion     =  sc .* max(bord.fluxmur_c,0.5) .* sqrt(e0 .* prof.ti(end)) .* phys.e +  ...
                  sf .* max(bord.fluxmur_f,0.5) .* e0 .* phys.e ; 
source.ion(1)  =  0;

% profil de neutre thermique
n0th           = (sf .*  max(bord.fluxmur_f,0.5) .* (lif + 5 .*lcxf) ./ v0th +  ...
                  sc .*  max(bord.fluxmur_c,0.5) .* lcxc ./ v0i) .* equi.rhomax;
n0th(1)        = 0;

% profil de neutres par especes n'est pas implante dans ce module
% proportionnel au flux entrant
nn0= NaN .* ones(1,gene.nbrho,gene.nbg);
for k=1:gene.nbg
   n0       =  bord.fluxgazin_c(k) .* sc + bord.fluxgazin_f(k) .* sf;
   n0(1)    =  0; 
   nn0(1,:,k) = n0;
end

% calcul du terme de friction sur les neutres pour la rotation
% calcul de la section efficace
svcx       = sigvcxH(prof.te,e0)';
% terme de friction selon phi
source.w   = - svcx .* n0th .* prof.rot;
% projection sur B pour le module neoclassique
% le module NClass ne prend pas en compte ce genre de force
% pour le moment.
%source.wb   = ?
%source.qb   = ?

