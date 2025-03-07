% ZLIMMHD calcul primaire des limites MHD
%--------------------------------------------------------------
% fichier zlimmhd.m ->  zlimmhd
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule le seuil de declenchement des DDS (critere du shear critique)
% ref :  F. porcelli et all, Plasma Phys. Control Fusion, vol 38, p 2163-2186, 1996
%    
%
% syntaxe  :
%  
%      [cr,evx] = zlimmhd(cons,datak,mem,gene,compo,phys)
%
% entree :
%
%      cons            =    parametre propre a la fonction (param.cons.mhd.limite)
%      datak           =    structure de donnees data au temps d'interet.
%      mem             =    structure de donnees au temps precedent
%      gene            =    parametres generaux (param.gene)
%      compo           =    composition du plasma: charge et masse des atomes ( param.compo)
%      phys            =    constantes physiques (param.phys)
% 
% sortie :
% 
%     cr               =  compte rendu d'execution (0 = ok)
%     evx              =  structure des evenemnts 
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.9, du 30/04/2002.
% 
% 
% liste des modifications : 
%
%   * 30/05/2002 -> changement de definition de s1
%   * 30/05/2002 -> changement de definition de rp et rn
%   * 30/05/2002 -> suppression de la dependance en rn dans scrit
%
%--------------------------------------------------------------
%
function [cr,evx] = zlimtemps(cons,datak,mem,gene,compo,phys)

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 
	valeur.cstab   = 10;                % coefficient de stabilite des dds
	valeur.beta    = 0;                  % choix du beta : 0 -> beta ionique , 1 -> beta ionique + supra thermique
	
	type.cstab     = 'float';                % type reel
	type.beta      = 'integer';     
	
	borne.cstab    = [0.1,10];               % valeurs possible 
	borne.beta     = {0,1};               % valeurs possible 
	
	defaut.cstab   = 0.025;                  % valeurs par defaut
	defaut.beta    = 0;                  % valeurs par defaut
	
	info.cstab     = 'coefficient de stabilite des dds';
	info.beta     = 'choix du beta : 0 -> beta ionique , 1 -> beta total (y compris supra thermique)';
	
	interface.ts = '';      % nom de la fonction d'interfacage avec les donnees TS
	interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;
	
	sortie.description = 'limite de stabilite des DDS (modele du shear critique F.Porcelli)';   % description (une ligne) de la fonction
	
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr = sortie;
	return
end



% initialisation des sortie
evx.dds = 0;
evx.elm = 0;
cr      = 0;

%
% calcul de l'indice de la surface q = 1
%
ind1 = max(find(datak.prof.q <= 1));
if isempty(ind1)
	% pas de dds
	return
elseif gene.x(ind1) < 0.1
   % top petite pas prise en compte
	return
end

% shear en q =1 
r1     = datak.equi.a(ind1) .* sqrt(datak.equi.e(ind1));
inds1  = min(find(datak.equi.a >= r1 ));
s1 = datak.prof.shear(inds1);

% calcul du seuil
mu       = sum(compo.a .* squeeze(datak.impur.impur(1,1,:))') ./ sum(squeeze(datak.impur.impur(1,1,:)));

tau      = datak.prof.te(ind1) ./ datak.prof.ti(1);
alpha    = 1.5 .* 3 .^ -(7/6) .* (tau ./ (1+tau)) .^ (7/12);
S        = 5.7e9 ./ sqrt(mu ./ 2.5) .* ( 8.1 .* datak.equi.ri(1)) .* (sqrt(datak.equi.b2(1)) ./ 5.7) ./ ...
           sqrt(datak.prof.ne(1) ./ 1e20) .* (datak.prof.te(1) ./ 20e3) .^ (3/2) .* (r1 ./ 1.6) .^ 2; 
rhoi     = 0.4e-2 .* sqrt(mu ./ 2.5) .* sqrt(datak.prof.te(1) ./ 20e3) ./ (sqrt(datak.equi.b2(1)) ./ 5.7);
rhohat   = rhoi ./ r1;

% ca semble mieux marcher
if cons.beta == 0
   betai1   = 2 .* phys.mu0 .* datak.prof.pion(ind1) ./ datak.equi.b2(ind1);
else	
   betai1   = 2 .* phys.mu0 .* (datak.prof.ptot(ind1) - datak.prof.pe(ind1)) ./ datak.equi.b2(ind1);
end

lmax     = 2.* datak.geo.r0;
lmin     = datak.equi.rhomax ./ length(gene.x);

%rn       = datak.prof.lne(ind1);
rn       = datak.prof.ne(1)  .* r1 ./max(1,abs(datak.prof.ne(1) - datak.prof.ne(inds1)));
rn       = min(lmax,max(lmin,rn));
gpt      = rpdederive(gene.x,datak.prof.ptot,0,2,2,1) ./ datak.equi.rhomax .* datak.equi.grho;
%rp       = lcalc(datak.prof.ptot,gpt,lmin,lmax);
%rp       = rp(ind1);
rp       = datak.prof.ptot(1)  .* r1 ./max(eps,abs(datak.prof.ptot(1) - datak.prof.ptot(inds1)));
rp       = min(lmax,max(lmin,rp));
%scrit    = cons.cstab .* alpha .* sqrt(S .^ (1/3) .* rhohat) .* (betai1 .* datak.equi.r2(ind1) ./r1 .^ 2) .^ (7/12) .* ...
%           (r1 ./ rn) .* (r1 ./ rp) .^ (1/6)
scrit    = cons.cstab .* alpha .* sqrt(S .^ (1/3) .* rhohat) .*  ...
           (betai1 .* datak.equi.r2(ind1) ./r1 .^ 2) .^ (7/12) .* ...
           (r1 ./ rp) .^ (1/6);

fprintf('x1 =  %g, s1 = %g, scrit = %g \n',gene.x(ind1),s1,scrit);
% trigger des DDS	
if	s1 > scrit	
	evx.dds  =1 ;
end
% fin de la fonction


% calcul de la longueur de gradient
function out = lcalc(v,gv,lmin,lmax)

ind = find( abs(gv) < eps);
if ~isempty(ind)
	 gv(ind) = eps;
end
out = abs(v ./ gv);
if ~isempty(ind)
	 out(ind) = lmax;
end
out = min(lmax,max(lmin,out));


