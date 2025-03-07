% ZHYBSIMPLE profil de depot hybride simple
%--------------------------------------------------------------
% fichier zhybsimple.m ->  zhybsimple, zintvol, cur_eff, polys,
%                          cur_plh, cur_eff_lh
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule de maniere simple un depot pour l'Hybride.
% 
% syntaxe  :
%  
%      [sortie,memoire] = zhybsimple(parametre,proto,cons,geo,equi,injection, ...
%                             prof,neo,impur,phy,composition,gene,memoire);
%
% entree :
%
%      parametre       =    parametre propre a la fonction (param.cons.hyb)
%      proto           =    prototype de la structure pour les sources, valeurs a zeros (proto = zsourceproto;)
%      cons            =    consigne de puissance par coupleur (data.cons.hyb)
%      geo             =    geometrie du plasma (data.geo)
%      equi            =    donnees de l'equilibre plasma (data.equi)
%      injection       =    consigne d'injection de gaz (data.cons.c)
%      prof            =    profils des donnees calculees par le code (data.prof) 
%      neo             =    donnees neoclassiques (data.neo)
%      impur           =    sous strcuture des impurtes (data.impur)
%      phy             =    constantes physiques (param.phys)
%      composition     =    composition du plasma: charge et masse des atomes ( param.compo)
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
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 08/02/2005.
% 
% 
% liste des modifications : 
% 
%  * 20/03/2001 -> modulation de l'efficacite de generation de courant avec l'activite MHD
%  * 02/10/2001 -> ajout de la structure memoire en sortie
%  * 24/01/2003 -> synergie FCE et LH
%  * 11/03/2003 -> suppression de la proportionnalite a 1/ne du profil de courant LH lorsqu'on utilise les Xdurs (j proportionnel a P/ne de facon systematique avant)
%  * 29/07/2003 -> efficacite par defaut a 5.8e18 
%  * 28/08/2003 -> possibilite de prendre en compte les pertes ripple pour TS (remodifie 02/09/03)  
%  * 30/10/2003 -> suppression du parametre gammastar, non utilise 
%  * 10/03/2004 -> ajout efficacite variable au cours du temps (data.cons.asser.etalh) 
%  * 27/04/2004 -> suppression des parametre npar, frequence, fraction, non utilises
%  * 13/09/2004 -> ajout de l'effet de la conductivite chaude (correction d'ordre 1)
%  * 08/02/2005 -> ajout des scaling d'efficacite de generation de courant 
%--------------------------------------------------------------
%
function [sortie,memoire] = zhybsimple(parametre,proto,cons,geo,equi,injection, ...
                             prof,neo,impur,phy,composition,gene,memoire)
                             
% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	if nargin ==0
		nbhyb=1;
	else
		nbhyb=parametre;   
	end
	
        valeur.machine           = 'TS';                     % tokamak	
	valeur.centre            = zeros(1,nbhyb);           % centre de la gaussienne 
	valeur.largeur           = 0.2.*ones(1,nbhyb);       % largeur de la gaussienne 
	valeur.efficacite        = 0.58e19*ones(1,nbhyb);    % efficacite de generation de courant ( A/W/m^2)
        valeur.scaling           = 'Input';
	valeur.xdur              = 1;                        % si 1 utilise prof.xdur a la place de la gaussienne pour la forme du depot
	valeur.effet_mhd         = 0;                        % prise en compte de la mhd sur le profil de j
	valeur.ripple            = 0;                        % prise en compte du ripple sur TS
	valeur.etavar            = 'No';                        % prise en compte d'une efficacite variable
	valeur.etachaud            = 'No';                        % prise en compte de la correction pour la conductivite chaude
	
        type.machine             = 'string';                 % configuration machine	
	type.centre              = 'float';                  % type reel
	type.largeur             = 'float';                  % type reel
	type.fraction            = 'float';                  % type reel
	type.efficacite          = 'float';                  % type reel
        type.scaling             = 'string';                 % configuration machine	
	type.xdur                = 'integer';                % type entier
	type.effet_mhd           = 'integer';                % type entier
	type.ripple              = 'integer';                % type entier
        type.etavar              = 'string';                 % configuration machine	
        type.etachaud              = 'string';                 % configuration machine	
	
        borne.machine            = {'TS','Autres'};         % configuration machine	
	borne.centre             = [0,1];                    % valeurs possibles 
	borne.largeur            = [0,0.5];                  % valeurs possibles 
	borne.efficacite          = [-1e20,1e20];            % valeurs possibles 
        borne.scaling            = {'Input','Goniche','SimulTS'};
	borne.xdur               = {0,1};                    % valeurs possibles
	borne.effet_mhd          = {0,1};                    % valeurs possibles
	borne.ripple             = {0,1};                    % valeurs possibles
        borne.etavar             = {'Yes','No'};         % configuration machine	
        borne.etachaud             = {'Yes','No'};         % configuration machine	
	
        defaut.machine           = 'TS';                     % configuration machine	
	defaut.centre            = 0;                        % valeurs par defaut 
	defaut.largeur           = 0.2;                      % valeurs par defaut 
	defaut.efficacite        = 0.58e19;                  % valeurs par defaut
        defaut.scaling           = 'Input';
	defaut.xdur              = 1;                        % valeurs par defaut
	defaut.effet_mhd         = 0;                        % valeurs par defaut
	defaut.ripple            = 0;                        % valeurs par defaut
        defaut.etavar            = 'No';                     % configuration machine	
        defaut.etachaud            = 'No';                     % configuration machine	
	
	
        info.machine             = 'configuration machine (pour le ripple)';	
	info.centre              = 'centre de la gaussienne' ;
	info.largeur             = 'largeur de la gaussienne'; 
	info.efficacite          = 'efficacite de generation de courant';
        info.scaling             = 'Input -> donnee par le scalaire efficacite ou par data.cons.asser.etalh selon la valeur de etavar, Goniche -> loi d''echelle de M. Goniche, SimulTS -> loi d''echelle SimulTS';
	info.xdur                = 'si 1 utilise prof.xdur a la place de la gaussienne pour la forme du depot';
	info.effet_mhd           = 'prise en compte de la mhd sur le profil de j';
	info.ripple              = 'prise en compte des pertes ripple pour TS';
	info.etavar              = 'efficacite LH variable au cours du temps (data.cons.asser.etalh)';
	info.etachaud              = 'prise en compte de la correction pour la conductivite chaude';
	
	interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'Calcul simple du depot de puisance Hybride';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	
	
	return
end

if ~isfield(parametre,'effet_mhd')
  parametre.effet_mhd=0;
end
if ~isfield(parametre,'scaling')
  parametre.scaling='Input';
end

% debut du calcul
% petits vecteurs utils
va = ones(gene.nbhyb,1);
ve = ones(1,gene.nbrho);

% Parametres du Ripple (modif 02/09/03)  (loi d'???helle pour le diagnostic DRIPPLE)
alpha = [0.6662639631967 -0.503759394462734 0.235509558581612;0.942854084818663 -0.212233746620332 -0.606046698438498;1.17449344302278 0.406956986489841 -1.55647263246365;1.26490295439731 0.807360991632855 -1.77323814523948;1.34391626203478 0.829563873473516 -1.67640737283367;1.27263426767301 0.445323609124571 -1.06744934260616;0.886779152061443 -0.4217528815393 -0.069881000409826;0.146970707366509 -0.758932144063924 0.917539880313246];
dalpha = [0.00155681373475998 0.00518733657582794 0.00765538869327886;0.00104073495105042 0.00346775105896577 0.00511765177752616;0.000831525873545265 0.00277066194964786 0.00408889877342242;0.000773591595767396 0.00257762370017632 0.00380401599961882;0.000934569549592637 0.00311400567647508 0.00459559997660044;0.00123139407027466 0.00410303131156422 0.00605518825539183;0.00152775646084226 0.00509051711925342 0.0075125040814332;0.00186667541971748 0.00621980232040562 0.0091790851937274];
c0 = [2.85482634661752;2.24535175106449;8.06969332838054;14.8139712278853;10.4570952588182;3.32534745291915;0.533723654679241;0.0812120305892144];

% selon les parametres le lieu du maximum du depot est donne ou calcule
centre = parametre.centre;
nplow  = NaN .* ones(1,gene.nbhyb);
npup   = NaN .* ones(1,gene.nbhyb);
floss  = NaN .* ones(1,gene.nbhyb);

% Dans le passe, on pouvait faire appel a une fonction de calcul du centre de la gaussienne, de type First pass / Wave diffusion (?)
% option supprimee
%if parametre.calcul_centre == 1
%	for k = 1:gene.nbhyb
%	    [centre(1,k) , nplow(k) , npup(k) ,floss(k)] =cur_plh(gene.x',equi.rhomax,geo.r0,geo.b0,composition.z(1),composition.a(1), ...
%	                                                          prof.te',prof.ne',squeeze(impur.impur(1,:,1)),prof.q', ...
%	                                                          parametre.npar(k),parametre.freq_ghz(k),gamstar);
%	end
%end

% puissance
p = exp( - ((va*gene.x - centre' * ve)./ (parametre.largeur' * ve)).^2);
if (parametre.xdur == 1)
	if all(isfinite(prof.xdur))
		if any(prof.xdur >0)
			p =(prof.xdur' * ones(1,size(p,1)))';
			disp('utilisation du profil xdur');
		end
	end	
end 

% normalisation a la puissance de consigne
pint    = trapz(equi.rhomax .* gene.x , p .* (va * equi.vpr),2);
nbar    = trapz(gene.x,prof.ne);
plhcons = sum((abs(cons)'));


if parametre.ripple == 1 & strcmp(parametre.machine,'TS')
  %
  % premire version en attendant une loi d'???helle sur le pertes par calorim???rie
  % loi d'???helle pour l'???ergie
  escal         = [0.6 0.4];
  %
  % loi d'???helle pour l'??????ation en temp???ature
  %
  cT            = 3.45;
  aT            = [1.74 -1.92 1.06];
  %  
  % estimation nl voie 2 de TS
  %
  nl            = nbar/1e19 * 2 * geo.a * 0.72;
  rip.profil    = c0 .* (plhcons/1e6)    .^  alpha(:,1) .* ...
                    (equi.ip/1e6)    .^  alpha(:,2).* ...
		    (nl)             .^  alpha(:,3);
  rip.Erip      = 130*(equi.ip/1e6) .^ 0.6 .* (nl).^0.4;
  rip.ploss     = (ones(8,1)*rip.Erip) .* rip.profil * 18;
  %
  % facteur 1.4 tenant compte du fait que DRIPPLE ne voit que 70 % des pertes 
  %
  rip.ptot      = sum(rip.ploss)*1.4;
  rip.fr        = rip.ptot / plhcons;
  p             = p .* ((abs(cons)' ./ pint) * ve) * (1-rip.fr);
else
  p             = p .* ((abs(cons)' ./ pint) * ve);
end
% repartition de la puissance et ecriture des sorties
sortie      = proto;
sortie.el   = sum(p,1); % somme sur les antennes

% generation de courant 
%nemoy = zintvol(prof.ne,gene.x,equi.vpr,equi.rhomax) ./ ...
%        zintvol(ones(size(equi.vpr)),gene.x,equi.vpr,equi.rhomax);

zeffm  = trapz(gene.x,prof.zeff);
% efficacite LH (=parametre.efficacite si non nul, scaling sinon)
%--------------
switch parametre.scaling
case 'SimulTS'
    
   % scaling simults
   npar  = abs(angle(cons.'));
   phase = max(0,(npar - 1.8) .* 230);
   %cst    = 0.75e19.*sqrt(2.1+3) ./ 4 .* 7;
   %etalh = cst ./ sqrt(nbar ./ 1e19+3).* bt ./ (zeff+5) .* (1-phase./230) .* sign(option.etalh); % scaling etaLH
   cst    = 0.75*sqrt(2.1+3)/4*7;
   eta    = ((cst./sqrt(nbar./1e19+3)*geo.b0./(zeffm+5)) .* va).*(1 - phase ./ 230) .* 1e19 ;  % scaling etaLH

case 'Goniche'

   %phase =  angle(cons') .* 360;
   %npar  =  phase ./ 230 + 1.8;
   Dn    = 2.01 - 0.63 .*  abs(angle(cons.')); 
   %tem   = trapz(gene.x,prof.te .* equi.vpr,2) ./ trapz(gene.x,equi.vpr,2) ./1e3; 
   eta   = min(2.4e20 ./ (5+zeffm),max(1e17,1.18e19 .* Dn .^ 0.55 .* (zeffm * va) .^ (-0.24) .*  ((equi.ip./1e6) * va) .^ 0.43 .* (geo.b0 * va) ./ 3.6));
   
otherwise
   eta    =  parametre.efficacite';   
   switch parametre.etavar
   case 'Yes'

      eta = parametre.machine_cons.etalh;

   end
end
% courant genere (en fonction de l'efficacite)
icd  = eta .* abs(cons') ./ geo.r0 ./ nbar;
iboot=trapz(equi.rhomax*gene.x,neo.jboot.*equi.spr);

if sum(icd(:)) > (0.95*equi.ip-iboot)
   disp('LH current correction as Ilh is greater than 95% of (Ip-Iboot)')
   sortie.err = (1-(0.95.*equi.ip-iboot)/sum(icd(:)))*100; % le courant LH a ete reduit de sortie.err
   icd = (0.95.*equi.ip-iboot) .* icd ./sum(icd(:));  % maintains the distribution of current between the antennas
   disp(['ILH current has been reduced of ', num2str(sortie.err), ' %'])
end






warning off
if (parametre.xdur == 1)
   fact =  p;             % j proportionnel aux Xdurs
else
   fact =  p ./ (va * prof.ne);    % j proportionnel a P/n
end   
ind = find(~isfinite(fact));
if ~isempty(ind)
	fact(ind) = zeros(1,length(ind));
end
inte =  fact .* (va * equi.spr);
jcd  = ((icd ./ trapz(equi.rhomax .* gene.x ,inte,2))* ve) .* fact;  % normalisation jcd
ind = find(~isfinite(jcd)); 
if ~isempty(jcd)
	jcd(ind)  = zeros(1,length(ind));
end
%
% effet conductivite chaude
%
if strcmp(parametre.etachaud,'Yes')

  cons = pi^2 / 2;
  sigmahot = cons .* p ./ (va*prof.ne).^2 .* (5+va*prof.zeff).^2 ./ (3+va*prof.zeff) .* eta.^2;
  if isfield(proto,'sigmasyn')
	 sortie.sigmasyn = sum(sigmahot,1);
  else
  	jhot = (va*prof.epar) .* sigmahot;
  	jcd = jcd + jhot;   
  end
elseif isfield(proto,'sigmasyn')
 	sortie.sigmasyn(:) = 0; 
end

% modulation de jcd avec la mhd
if (parametre.effet_mhd == 1) & all(all(isfinite(prof.mhd_cd)))
	%disp('prise en compte de l''activite mhd');
	jcd = jcd .* (va * prof.mhd_cd);
end

warning on
sortie.j = sum(jcd,1);
%
% synergie FCE, LH, renormalisation sur la puissance LH injectee
%
if isfield(proto,'synergie')
%  sortie.j = sortie.j + proto.synergie .* sum(abs(cons))/1e6;
%
% synergie ind???endante de Phyb
%
  sortie.j = sortie.j + proto.synergie ;
end


