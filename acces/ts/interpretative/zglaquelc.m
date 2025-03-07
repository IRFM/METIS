% ZGLAQUELC interface de GLAQUELC (module de calcul du depot des glacons avec derive)
%--------------------------------------------------------------------------------------
% fichier zglaquelc.m ->  zglaquelc
%
%
% fonction Matlab 5 :
%
% Cette fonction interface le module de calcul de depot de glacons avec derive GLAQUELC.
% 
% syntaxe  :
%  
%     [dne_out,dni_out] = zglaquelc(par,geo,equi,cons,prof,neo,impur,bord,source,phys,compo,gene)
%
% entree :
%
%      par             =    parametre propre a la fonction (param.cons.fci)
%      cons            =    consigne de puissance par coupleur (data.cons.fci)
%      geo             =    geometrie du plasma (data.geo)
%      equi            =    donnees de l'equilibre plasma (data.equi)
%      cons            =    consigne de puissance par coupleur (data.cons.fci)
%      prof            =    profils des donnees calculees par le code (data.prof) 
%      neo             =    donnees neoclassiques (data.neo)
%      impur           =    sous strcuture des impurtes (data.impur)
%      bord            =    sous structure decrivant le bord (data.bord) 
%      source          =    sous structure des sources (data.source)
%      phys            =    constantes physiques (param.phys)
%      compo           =    composition du plasma: charge et masse des atomes ( param.compo)
%      gene            =    parametres generaux (param.gene)
% 
% sortie :
% 
%     dne_out          =    depot de matiere  en nombre d'electrons (m^-3) [1,nbrho]
%     dni_out          =    depot de matiee par espece (m^-3) [1,nbrho,nbg]
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.21, du 29/11/2004.
% 
% 
% liste des modifications : 
% * 29/11/2004 --> correction dans inhomog.m a la fin (pb d'interpolation, a valider par concepteur
%--------------------------------------------------------------
%
function [dne_out,dni_out,dpe_out,dpion_out,duration] = zglaquelc(par,geo,equi,cons,prof,neo,impur,bord,source,phys,compo,gene)

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	 if nargin ==0
	 	nb=1;
	 else
	 	nb=par;   
	 end
    valeur.ablat       = 1;                   % modele d'ablation 1 = NGS-TS
    valeur.deriv       = 0;                   % calcul de la derive si = 1
    valeur.graph       = 0;                   % si >0 trace le graphe d'injection
    valeur.rapide      = 0;                   % si =1 geometrie simplifiee, accelere le calcul
    valeur.mix         = 0.5;                 % melange neutre/ion [0..1]
    valeur.nu          = 0;                   % correction de la dependance du modele de derive en Te (0)
    valeur.rinj        = 10 .* ones(1,nb);    % rayon du point initial de la trajectoire d'injection (m)
    valeur.zinj        = zeros(1,nb);         % azimut du point initial de la trajectoire d'injection (m)
    valeur.angleinj    = 180.*ones(1,nb);     % angle de la trajectoire d'injection en degres (sens trigo, 0 = HFS ; 180 = LFS) 
    valeur.phiglac     = [];                  % diametre du glacon (m) { preciser soit le diametre soit le nombre d'atome}
    valeur.np          = 5e20;                % nombre d'atome dans le glacon (su) { preciser soit le diametre soit le nombre d'atome}
    valeur.vitesse     = 500 .* ones(1,nb);   % vitesse d'injection (m/s)
    valeur.amass       = 2 .* ones(1,nb);     % nombre de masse de l'isotope de l'hydrogene utilise pour le glacon
    valeur.prob        = ones(1,nb);          % probabilite d'injection [0..1]
       
    type.ablat       = 'integer';             % type des donnees
    type.deriv       = 'integer';
    type.graph       = 'integer';
    type.rapide      = 'integer';
    type.mix         = 'float';              
    type.nu          = 'float';              
    type.rinj        = 'float';
    type.zinj        = 'float'; 
    type.angleinj    = 'float';
    type.phiglac     = 'float';
    type.np          = 'float';
    type.vitesse     = 'float';
    type.amass       = 'integer';
    type.prob        = 'float';
       
    borne.ablat       = {1};                 % bornes
    borne.deriv       = {0,1};         
    borne.graph       = {0,1};         
    borne.rapide      = {0,1};         
    borne.mix         = [0,1];
    borne.nu          = [-10,10];            
    borne.rinj        = [-inf,inf];
    borne.zinj        = [-inf,inf];
    borne.angleinj    = [0,360];
    borne.phiglac     = [0,0.1]; 
    borne.np          = [1e17,1e24]; 
    borne.vitesse     = [1,100e3];
    borne.amass       = {1,2,3}; 
    borne.prob        = [0,1];
       
    defaut.ablat       = 1;                   %valeur par defaut
    defaut.deriv       = 0;
    defaut.graph       = 0;
    defaut.rapide      = 0;
    defaut.mix         = 0.5;
    defaut.nu          = 0;
    defaut.rinj        = 10;
    defaut.zinj        = 0;
    defaut.angleinj    = 180;
    defaut.phiglac     = [];
    defaut.np          = 5e20;
    defaut.vitesse     = 500;
    defaut.amass       = 2;
    defaut.prob        = 1;
       
    info.ablat       = 'modele d''ablation 1 = NGS-TS';
    info.deriv       = 'calcul de la derive si = 1';
    info.graph       = 'si >0, trace le graphe d''injection';
    info.rapide      = 'si =1 geometrie simplifiee, accelere le calcul';
    info.mix         = 'melange neutre/ion [0..1]';
    info.nu          = 'correction de la dependance du modele de derive en Te (0)';
    info.rinj        = 'rayon du point initial de la trajectoire d''injection (m)';
    info.zinj        = 'azimut du point initial de la trajectoire d''injection (m)';
    info.angleinj    = 'angle de la trajectoire d''injection en degres (sens trigo, 0 = HFS ; 180 = LFS)'; 
    info.phiglac     = 'diametre du glacon (m) { preciser soit le diametre soit le nombre d''atome}';
    info.np          = 'nombre d''atome dans le glacon (su) { preciser soit le diametre soit le nombre d''atome}';
    info.vitesse     = 'vitesse d''injection (m/s)';
    info.amass       = 'nombre de masse de l''isotope de l''hydrogene utilise pour le glacon';
    info.prob        = 'probabilite d''injection [0..1]';
	
    valeur.frequency   = zeros(1,nb);
    type.frequency      = 'float';
    borne.frequency      = [0,100];
    defaut.frequency   = 0;
    info.frequency   = 'pellet injection frequency in Hz; if = 0, injection is controled by cronos referencences (param.gene.evx_inter must be set to zero if pellet injection frequency is used, otherwise only injection controled by cronos works)';

    interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
    interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
    sortie.valeur=valeur;
    sortie.type=type;
    sortie.borne=borne;
    sortie.defaut=defaut;
    sortie.info=info;
    sortie.interface=interface;
	
    sortie.description = 'Calcul du depot de matiere de l''injection de glacon utilisant GLAQUELC (Parks/Fois)';   % description (une ligne) de la fonction
	
    sortie.help = '';             % nom du fichier d'aide s'il existe, sinon aide de la fonction
    sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
    sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
    dne_out =sortie;
	
    return
end

% rho local
rhoin   = equi.a ./ equi.a(end);
dne_out = 0 .* prof.ne;
dpe_out = 0 .* prof.ne;
dni_out = 0 .* impur.impur;
dpion_out = 0 .* prof.ne;
duration = 1e-6;

% boucle sur les injecteurs
for k = 1:gene.nbglacon
   % test sur la consigne et la probabilite d'injection
   if (cons(k) > 0) & (par.prob(k) >= rand(1))
       % variables d'entrees
       graph    = par.graph;
       fich     = 0;
       repli    = 1;
       codes    = [graph,fich,par.ablat,repli,2 .* (par.deriv == 1),par.mix,par.nu];
       id       = 'XXL';
       Nech     = [gene.nbrho,1001];
       Xo       = par.rinj(k);
       Zo       = par.zinj(k);
       teta     = par.angleinj(k) ./ 180 .* pi;
       if isempty(par.np)
       	rp       = par.phiglac(k) ./ 2;
       	Np       = [];
       else
       	rp       = [];
       	Np       = par.np(k);
       end
       Vp       = par.vitesse(k);
       Ap       = [par.amass(k),compo.a(1)];
       rhoneo   = rhoin;
       neo      = prof.ne;
       rhoTeo   = rhoin;
       Teo      = prof.te;
       Ro       = geo.r0;
       a        = geo.a;
       xx       = linspace(0,1,length(rhoin));
       if par.rapide == 1
       	Do       = equi.d(1);
       	ko       = equi.e(end);
       	delo     = 0.5 .* (equi.trh(end) + equi.trl(end));
       else
       	Do       = interp1(rhoin,equi.d,xx)';
       	ko       = interp1(rhoin,equi.e,xx)';
       	delo     = interp1(rhoin,0.5 .* (equi.trh + equi.trl),xx)';
       end
       Bt       = geo.b0;
       Ip       = equi.ip;

       % appel de la fonction de calcul du depot des glacons
      [codes,chcodes,pivot,rhop,rhom,rhocm,rhoL,effinj,errdep,rho,DN,dne,ne,dTe,Te,ABL,EQM]=...
       GLAQUELC(codes,id,Nech,Xo,Zo,teta,rp,Np,Vp,Ap,rhoneo,neo,rhoTeo,Teo,Ro,a,Do,ko,delo,Bt,Ip);

       %save last_glacon
       
       % reechantillonage de dne sur x
       dne_r   = interp1(rho',dne(:,end)',rhoin,'spline');
       if any(dne_r < 0)
              dne_r   = interp1(rho',dne(:,end)',rhoin,'linear');
       end
       if any(dne_r < 0)
              dne_r   = interp1(rho',dne(:,end)',rhoin,'nearest');
       end
       if any(dne_r < 0)
              dne_r   = dne_r .* (dne_r >0);
       end
       % coefficient de renormalisation pour tenir compte de la vrai geometrie
       coefnorme = effinj .* ABL(1,7) ./ trapz(gene.x.*equi.rhomax,equi.vpr.*dne_r);

       
       dne_out = dne_out + dne_r .* coefnorme;          
       % calcul de dni
       dni_out = shiftdim((dne_r'.*coefnorme) * ((compo.a == par.amass(k))./compo.z),-1) + dni_out;
   end
end
