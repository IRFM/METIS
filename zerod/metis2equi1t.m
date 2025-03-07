% fonction pour calcul d'equilibre sans convergence
function  [profil,deuxd,moment,scalaire,facteur,Ip_lcfs] = metis2equi1t(profil,geo,phys,nbth,facteur,exposhape,regularisation,factor_two_pi)


if nargin < 8
  factor_two_pi = 1;
end
% donnees de metis 
%       profil.x;
%	profil.kx     
%	profil.dx     
%	profil.rmx  
%       profil.raxe

% vecteur utiles
if isfield(profil,'x')
    x  = profil.x;
else
    x  = profil.xli;
    profil.x = x;
end
    
ve = ones(size(x));
% dpsidx s'annule au centre
psid1    = min(-eps,pdederive(x,profil.psi,0,2,2,1));
% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
psid2    = pdederive(x,profil.psi,1,2,2,2);
if any(psid1 == -eps)
	psi = cumtrapz(x,psid1,2);
	psid2    = pdederive(x,psi,1,2,2,2);
end

profil.deltax = profil.Raxe - profil.Raxe(end);
Raxe          = profil.Raxe;
epsi          = max(eps,(geo.a * x)) ./ Raxe;
profil.epsi   = epsi;

% calcul de la geometrie 2d
deuxd = ggeo2d(profil.x,profil.Raxe,geo.a,profil.kx,profil.dx,geo.vp,geo.sp,geo.Rsepa,geo.Zsepa,nbth,profil.psi,exposhape,factor_two_pi);

% calcul des valeurs moyennes
fdia = profil.fdia;

[profil.volume,profil.surface,profil.phi,dvdx,dsdx,dphidx] = volphi(x,deuxd,fdia,factor_two_pi);
% calcul directe de q
profil.q = - pdederive(profil.psi,profil.phi,2,2,2,1) ./ 2 ./pi;
profil.q(1:6) = pchip(x(6:end),profil.q(6:end),x(1:6));
% rho cronos et vpr
profil.rho = sqrt(abs(profil.phi ./ pi ./ geo.b0));
profil.vpr = pdederive(profil.rho,profil.volume,0,2,2,1);
profil.spr = pdederive(profil.rho,profil.surface,0,2,2,1);
% calcul de gradient rho
deuxd = grho_comp(x,deuxd,profil.rho,factor_two_pi);

% moyenne geometrique pour le calcul de F et q
[profil.r2i,profil.ri,profil.rmoy,profil.r2,profil.rmx,profil.b2,profil.b2i] = griri2(x,deuxd,fdia,factor_two_pi);
% C3 (compatibilite ancienne version
profil.C3  = profil.vpr .* profil.r2i;
% BPHI
deuxd.BPHI = (profil.fdia'*ones(1,size(deuxd.R,2)))./deuxd.R;
%
deuxd.drhodpsi       = pdederive(profil.psi,profil.rho,0,2,2,1)' * ones(1,size(deuxd.R,2));
%deuxd.dRdPSI         = pdederive(profil.psi',deuxd.R,0,2,1,1);
%deuxd.dZdPSI         = pdederive(profil.psi',deuxd.Z,0,2,1,1);
warning off
deuxd.dRdPSI         = deuxd.dRdx ./ deuxd.dPSIdx;
deuxd.dZdPSI         = deuxd.dZdx ./ deuxd.dPSIdx;
warning on
deuxd.dRdPSI(1,:) = 0;
deuxd.dZdPSI(1,:) = 0;

% valeur moyenne avec grad rho
pass        = profil.q ./ max(eps,profil.rho) ./ geo.b0;
pass(1)     = pchip(profil.rho(2:end),pass(2:end),profil.rho(1));
profil.pass = pass;    
[profil.vpr_ctr,profil.grho2r2,profil.grho,profil.grho2, ...
 profil.grhor,profil.grho2b2,profil.bpolm] = ggg(profil.x,deuxd,fdia,regularisation);
% calcul de C2
profil.C2   = profil.grho2r2 .* profil.vpr;
profil.vpr_ctr = profil.vpr_ctr ./ pdederive(x,profil.rho,2,2,2,1);
% calcul de li
profil.bpol  = sqrt(profil.grho2r2) .* abs(pdederive(profil.rho,profil.psi,0,2,2,1));    
  
		     
% calcul de la fraction de pieges
if isfield(profil,'bpol') && isfield(profil,'fdia')
	b2m       = profil.bpol .^ 2  + profil.fdia .^ 2 .* profil.r2i;
	bm        = profil.fdia .* profil.ri .* (1 + 0.5 .* profil.bpol ./ sqrt(profil.grho2r2) ./ ...
	            profil.fdia .* sqrt(profil.grho2));
	% calcul de bmax
	btor         = (profil.fdia ./ (Raxe - geo.a * x));
%	grho         = abs((profil.rmx(:,end) * ve) ./ max(eps,abs(pdederive(x,profil.Raxe - geo.a * x,0,2,2,1))));
%	grho(:,1)    = 2 .* grho(:,2) - grho(:,3);
%	bpol         = abs(pdederive(x,profil.psi,0,2,2,1))./ (profil.Raxe - geo.a * x) .* ...
%	                grho ./ (profil.rmx(:,end) * ve);
	grho         = deuxd.gradrho(:,1)';
	bpol         = abs(pdederive(profil.rho,profil.psi,0,2,2,1))./ (profil.Raxe - geo.a * x) .* grho;
	bmax         = sqrt(btor .^ 2 + bpol .^ 2);
	% variable
	h  = min(1,bm ./ bmax);
	h2 = min(1,b2m ./ bmax .^ 2);
	 
	% expression de ftrap  Lin-Liu and Miller Phys. of Plasmas 2 (5) 1995
	ftu          = 1 - h2 ./ h .^ 2 .* ( 1 - sqrt(1 - h) .* (1 + 0.5 .* h));
	ftl          = 1 - h2 .* (1./ h2 - sqrt(1-sqrt(h2)) ./ h2  - sqrt(1 - sqrt(h2)) ./ 2 ./h);
	ftc          = 1 - (1-profil.epsi) .^ 2 ./ sqrt(1 - profil.epsi .^ 2) ./ (1 + 1.46 .* sqrt(profil.epsi));  % Wesson , Tokamak
	fte          = 0.75 .* ftu + 0.25 .* ftl;
	profil.ftrap = ftc .* (1 - x) + fte .* x;
	
else
	profil.ftrap     = 1 - (1-profil.epsi) .^ 2 ./ sqrt(1 - profil.epsi .^ 2) ./ (1 + 1.46 .* sqrt(profil.epsi));  % Wesson , Tokamak
end

% moments
R = squeeze(deuxd.R);
Z = squeeze(deuxd.Z);
vu = ones(1,size(Z,2));
Rmin = min(R,[],2);
Rmax = max(R,[],2);
[Zmin,lmin] = min(Z,[],2);
[Zmax,lmax] = max(Z,[],2);
dmin      = exp(-(Z - Zmin * ones(1,size(Z,2))) .^ 2 ./ (geo.a ./ 16) .^ 2);  
dmax      = exp(-(Z - Zmax * ones(1,size(Z,2))) .^ 2 ./ (geo.a ./ 16) .^ 2);
Rzmin     = sum(R .* dmin,2) ./ sum(dmin,2);
Rzmax    = sum(R .* dmax,2) ./ sum(dmax,2);
maskRmax = (R == (Rmax * vu));
moment.zaxe = sum(Z .* maskRmax,2) ./ max(1,sum(maskRmax,2));
moment.raxe  = (Rmax + Rmin) ./ 2;
moment.a     = (Rmax - Rmin) ./ 2;
moment.e     = (Zmax - Zmin) ./ max(eps,Rmax - Rmin);
moment.e(1)  = pchip(profil.x(2:end),moment.e(2:end),0);
moment.trh   = (moment.raxe - Rzmax) ./ max(eps,moment.a);
moment.trh(1) = 0;
indp = find(moment.trh >= 0);
if length(indp) > 7
	moment.trh   = pchip(x(indp),moment.trh(indp),x);
else
	moment.trh   = zeros(size(x));
end
moment.trl   = (moment.raxe - Rzmin) ./ max(eps,moment.a);
moment.trl(1) = 0;
indp = find(moment.trl >= 0);
if length(indp) > 7
	moment.trl   = pchip(x(indp),moment.trl(indp),x);
else
	moment.trl  = zeros(size(x));
end
% autres scalaires
scalaire.ipout   = trapz(profil.rho,profil.spr .* profil.jmoy,2);
scalaire.liout     = 2 .* trapz(profil.rho,profil.bpol .^ 2  .* profil.vpr,2) ./  ...
                     ((phys.mu0 .*scalaire.ipout ) .^ 2  .* geo.R);
% definition de CRONOS (ce n'est pas celle de METIS qui correspond au betap(1))
scalaire.betap  = trapz(profil.rho,profil.ptot .* profil.spr,2) ./ trapz(profil.rho,profil.spr,2) ./ (profil.bpol(end) .^2 ./ 2 ./ phys.mu0);   
scalaire.betat  = 2 .* phys.mu0 .* trapz(profil.rho,profil.ptot .* profil.vpr,2) ./ ...
                                   trapz(profil.rho,(profil.fdia .* profil.ri) .^2 .* profil.vpr,2);   
scalaire.betan  = 3./2 .* trapz(profil.rho,profil.ptot .* profil.vpr,2) .* (1.6.*pi./3) .* moment.a(end) ./...
                          trapz(profil.rho,profil.vpr,2) ./ geo.b0 ./ scalaire.ipout;			  
			  
% calcul de jmoy pour la mesure de convergence
jmoy       = -pdederive(profil.rho,profil.vpr .* profil.grho2r2 .* ...
             pdederive(profil.rho,profil.psi,0,2,2,1),0,2,2,1) ./  ...
	     (phys.mu0 .* max(eps,profil.vpr) .* profil.ri);
jmoy(1)    = pchip(profil.rho(2:end),jmoy(2:end),0);
profil.jmoy = jmoy;


% calcul de dthetadR et dZ
deuxd.dthdR = - deuxd.dZdx ./ max(eps,deuxd.RZjac);
deuxd.dthdR(1,:) = 0;
deuxd.dthdZ =  deuxd.dRdx ./ max(eps,deuxd.RZjac);
deuxd.dthdZ(1,:) = 0;

% valeur en circulaire
%  rc = (deuxd.R - ones(size(deuxd.R,1),1) * deuxd.R(1,:));
%  zc = (deuxd.Z - ones(size(deuxd.Z,1),1) * deuxd.Z(1,:));
%  dthdr = - zc ./ rc .^ 2 ./ (1 + zc .^ 2 ./ rc .^ 2); 
%  dthdz =  1 ./ rc  ./ (1 + zc .^ 2 ./ rc .^ 2); 
% computing circulation of B on LCFS for diagnostic purpose
for k=1:size(deuxd.R,1)
  Ip_lcfs(k) = sign(factor_two_pi) .* sum(((deuxd.BR(k,1:end-1) + deuxd.BR(k,2:end)) .* diff(deuxd.R(k,:)) ./ 2) + ...
          ((deuxd.BZ(k,1:end-1) + deuxd.BZ(k,2:end)) .* diff(deuxd.Z(k,:)) ./ 2)) ./ (4.*pi.*1e-7);
end

% coefficient de geometrie
function deuxd = ggeo2d(x,Raxe,a,kx,d,Vp,Sp,Rsepa,Zsepa,nbth,psi,expo,factor_two_pi)

% securite sur les valeurs d'entree
Raxe = max(0.01,min(1e3,Raxe));
a    = max(0.01,min(0.9 .* Raxe(:,end),a));
kx   = max(0.5,min(10,kx));
d    = max(0,min(0.999,d));

% donnees utiles
ve = ones(size(x));
amin  = a  * ve;
t  =  asin(d);
R0    = Raxe(end) * ve;

% angle theta
th   = linspace(0,2.*pi,nbth);
vh   = ones(size(th));
% separatrice donnees par les moments
Rmom   = R0(end) * vh + (a * vh) .* cos(th + t(end) * sin(th));
Zmom   = (a .* kx(end)) * sin(th);


% donnees 2 D en sortie
deuxd.th        = ones(length(x),1) * th;
deuxd.R        = NaN .* ones(length(x),length(th));
deuxd.Z        = NaN .* ones(length(x),length(th));
deuxd.rhog     = NaN .* ones(length(x),length(th));
deuxd.dRdx     = NaN .* ones(length(x),length(th));
deuxd.dZdx     = NaN .* ones(length(x),length(th));
deuxd.dRdth    = NaN .* ones(length(x),length(th));
deuxd.dZdth    = NaN .* ones(length(x),length(th));
deuxd.cs       = NaN .* ones(length(x),length(th));
deuxd.ss       = NaN .* ones(length(x),length(th));
deuxd.PSI      = psi' * ones(1,length(th));
deuxd.dPSIdx   = NaN .* ones(length(x),length(th));
%  deuxd.dPSIdth  = NaN .* ones(length(x),length(th));
deuxd.BR       = NaN .* ones(length(x),length(th));
deuxd.BZ       = NaN .* ones(length(x),length(th));
deuxd.RZjac    = NaN .* ones(length(x),length(th));
deuxd.dPSIdR   = NaN .* ones(length(x),length(th)); 
deuxd.dPSIdZ   = NaN .* ones(length(x),length(th)); 

% boucle pour la correction de forme
for l = 1:length(x)

%  	if l > 1
%  		Rm = R;
%  		Zm = Z;
%  		R  = Rp;
%  		Z  = Zp;
%  		rho2m = rho2;
%  		rho2  = rho2p;
%  	else
		R = Raxe(:,l) * vh + (amin(:,l) * vh) .* (x(:,l) * vh) .*  ...
			cos(th + (t(:,l) * vh) .* sin(th));
		Z = (kx(:,l) * vh) .*(amin(:,l) * vh) .* (x(:,l) * vh) .* sin(th);
		% morphing
		if ~isempty(Rsepa) && ~isempty(Zsepa) && (expo > 0)
			[R,Z]   = zemorph(R,Z,Rsepa,Zsepa,Rmom,Zmom,x(1,l),expo);
		end
		Raxel       = 0.5 .* (min(R,[],2) + max(R,[],2));
		rho2    = (R - Raxel * vh) .^2 +  Z .^ 2;
%  	end
%  	if l < size(x,2)
%  		Rp = Raxe(:,l+1) * vh + (amin(:,l+1) * vh) .* (x(:,l+1) * vh) .*  ...
%  			cos(th + (t(:,l+1) * vh) .* sin(th));
%  		Zp = (kx(:,l+1) * vh) .*(amin(:,l+1) * vh) .* (x(:,l+1) * vh) .* sin(th);
%  		% morphing
%  		if expo > 0
%  			[Rp,Zp]   = zemorph(Rp,Zp,Rsepa,Zsepa,Rmom,Zmom,x(1,l+1),expo);
%  		end
%  		Raxelp       = 0.5 .* (min(Rp,[],2) + max(Rp,[],2));
%  		rho2p = (Rp - Raxelp * vh) .^2 + Zp .^ 2;
%  	end


	% donnees 2D (grille droite)
	Raxel       = Raxe(1);
	cx   = (R - Raxel) + sqrt(-1) .* Z;
	thx  = unwrap(angle(cx),[],2);
	rhox = abs(cx);
	thx(thx<0) = thx(thx<0) + 2 .* pi;
	[thx,indx] = sort(thx,2);
	rhox       = rhox(indx);
	rhox = cat(2,rhox,rhox,rhox);
	thx = cat(2,thx -2.*pi,thx,thx+2.*pi);
	indnok = find(any(diff(thx,1,2)<=0,1));
	thx(:,indnok) =[];
	rhox(:,indnok)  = [];
	% rhoo   = pchip(thx,rhox,th);
	rhoo   = spline(thx,rhox,th);
	deuxd.cs(l,:)    = cos(th);
	deuxd.R(l,:)     = Raxel + rhoo .* deuxd.cs(l,:);
	deuxd.ss(l,:)    = sin(th);
	deuxd.Z(l,:)     = rhoo .* deuxd.ss(l,:);
	deuxd.rhog(l,:)  = rhoo;
	
	
end

% la coordonnee inetrne est a*x
%deuxd.dRdx           =  pdederive(x',deuxd.R,2,2,1,1);
%deuxd.dZdx           =  pdederive(x',deuxd.Z,2,2,1,1);
deuxd.dRdx           =  0.5 .* (pdederive(x',deuxd.R,1,2,1,1) + pdederive(x',deuxd.R,1,1,1,1));
deuxd.dZdx           =  0.5 .* (pdederive(x',deuxd.Z,1,2,1,1) + pdederive(x',deuxd.Z,1,1,1,1));
%drhogdx              =  pdederive(x',deuxd.rhog,2,2,1,1);
deuxd.dPSIdx         =  0.5 .* (pdederive(x',deuxd.PSI,0,2,1,1) +  pdederive(x',deuxd.PSI,0,1,1,1));
deuxd.dRdth          =  pdederive(th,deuxd.R,2,2,2,1);
deuxd.dRdth(:,1)     = 0.5 .* (deuxd.dRdth(:,1) + deuxd.dRdth(:,end));
deuxd.dRdth(:,end)   = deuxd.dRdth(:,1);
deuxd.dZdth          =  pdederive(th,deuxd.Z,2,2,2,1);
deuxd.dZdth(:,1)     = 0.5 .* (deuxd.dZdth(:,1) + deuxd.dZdth(:,end));
deuxd.dZdth(:,end)   = deuxd.dZdth(:,1);
%deuxd.RZjac          = abs(deuxd.dRdx .* deuxd.dZdth - deuxd.dZdx .* deuxd.dRdth);
deuxd.RZjac          = deuxd.dRdx .* deuxd.dZdth - deuxd.dZdx .* deuxd.dRdth;
if any(any(deuxd.RZjac(2:end,:) <= eps))
  indbad = find(any(deuxd.RZjac <= eps,2));
  dRdx_alt           =  pdederive(x',deuxd.R,1,1,1,1);
  dZdx_alt           =  pdederive(x',deuxd.Z,1,1,1,1);
  RZjac_alt          =  dRdx_alt .* deuxd.dZdth - dZdx_alt .* deuxd.dRdth;
  deuxd.dRdx(indbad,:)   =  dRdx_alt(indbad,:);
  deuxd.dZdx(indbad,:)   =  dZdx_alt(indbad,:);
  deuxd.RZjac(indbad,:)  =  RZjac_alt(indbad,:);
end
deuxd.RZjac          = abs(deuxd.RZjac);

% PSI ne depend pas de theta
deuxd.dPSIdR   =  deuxd.dPSIdx .* deuxd.dZdth  ./ max(eps,deuxd.RZjac);
deuxd.dPSIdZ   = -deuxd.dPSIdx .* deuxd.dRdth ./ max(eps,deuxd.RZjac);
% champ poloidal
deuxd.BR      =   deuxd.dPSIdZ ./ deuxd.R ./ abs(factor_two_pi);
deuxd.BZ      = - deuxd.dPSIdR ./ deuxd.R ./ abs(factor_two_pi);
deuxd.BR(1,:) = 0;
deuxd.BZ(1,:) = 0;

% coefficient de geometrie
function deuxd = grho_comp(x,deuxd,rho,factor_two_pi)

% passage a 2D
th        = deuxd.th(end,:);
deuxd.rho = rho' * ones(size(th));
deuxd.dRdrho = 0.5 .* (pdederive(rho',deuxd.R,1,2,1,1) + pdederive(rho',deuxd.R,1,1,1,1));
deuxd.dZdrho = 0.5 .* (pdederive(rho',deuxd.Z,1,2,1,1) + pdederive(rho',deuxd.Z,1,1,1,1));

% la coordonnee inetrne est a*x
%deuxd.drhodx              =  pdederive(x',deuxd.rho,2,2,1,1);
deuxd.drhodx              =  0.5 .* (pdederive(x',deuxd.rho,1,2,1,1) + pdederive(x',deuxd.rho,1,1,1,1));
deuxd.drhodth             =  pdederive(th,deuxd.rho,2,2,2,1);
deuxd.drhodth(:,1)        = 0.5 .* (deuxd.drhodth(:,1) +deuxd.drhodth(:,end));
deuxd.drhodth(:,end)      = deuxd.drhodth(:,1);

% calcul de gradrho
deuxd.drhodR   = (deuxd.drhodx .* deuxd.dZdth -  deuxd.drhodth .* deuxd.dZdx)  ./ max(eps,deuxd.RZjac);
deuxd.drhodZ   = (-deuxd.drhodx .* deuxd.dRdth + deuxd.drhodth .* deuxd.dRdx) ./ max(eps,deuxd.RZjac);
deuxd.gradrho = sqrt(deuxd.drhodR .^ 2 + deuxd.drhodZ .^ 2);



% fonction de morphing des surface de flux
% Rin, Zin : coordonnees de la surface non deformee
% Rsepa,Zsepa : coordonnees de la vrai separatrice
% Rmom,Zmom ; coordonnees de la separtrice donnees par les moments (ou de la dsmf)
% xin : coordonnees radiale de la surface de flux
% expo : exposant de la fonction de morphing
function [Rx,Zx] = zemorph(Rin,Zin,Rsepa,Zsepa,Rmom,Zmom,xin,expo)

% vecteur auxiliaire
vu = ones(1,size(Rin,2));
vs = ones(1,size(Rsepa,2));

% la separatrice est centree verticalement 
% calul de l'offset vertical
%maskrmax  = (Rsepa == (max(Rsepa,[],2) * ones(1,size(Rsepa,2))));
%z0         = sum(Zsepa .* maskrmax,2) ./ sum(maskrmax,2);
%Zsepa      = Zsepa - z0 * vs;
Raxe       = 0.5 .* (min(Rin,[],2) + max(Rin,[],2));

% deformation
cx   = (Rin - Raxe * vu) + sqrt(-1) .* Zin;
thx  = unwrap(angle(cx),[],2);
rhox = abs(cx);
thx(thx<0) = thx(thx<0) + 2 .* pi;
[thx,indx] = sort(thx,2);
for k = 1:size(indx,1)
	rhox(k,:)       = rhox(k,indx(k,:));
end

% separtrice a ce temps version analytique
cl   = (Rmom - Raxe * vu) + sqrt(-1) .* Zmom;
thl  = unwrap(angle(cl),[],2);
rhol = abs(cl);
thl(thl<0) = thl(thl<0) + 2 .* pi;
[thl,indl] = sort(thl,2);
for k = 1:size(indl,1)
	rhol(k,:)       = rhol(k,indl(k,:));
end
rhol = cat(2,rhol,rhol,rhol);
thl = cat(2,thl -2.*pi,thl,thl+2.*pi);
indnok = find(any(diff(thl,1,2)<=0,1));
thl(:,indnok) =[];
rhol(:,indnok)  = [];

% separtrice a ce temps complete
cc   = (Rsepa - Raxe * vs) + sqrt(-1) .* Zsepa;
thc  = unwrap(angle(cc),[],2);
rhoc = abs(cc);
thc(thc<0) = thc(thc<0) + 2 .* pi;
[thc,indc] = sort(thc,2);
for k = 1:size(indc,1)
	rhoc(k,:)       = rhoc(k,indc(k,:));
end
rhoc = cat(2,rhoc,rhoc,rhoc);
thc = cat(2,thc -2.*pi,thc,thc+2.*pi);
indnok = find(any(diff(thc,1,2)<=0,1));
thc(:,indnok) =[];
rhoc(:,indnok)  = [];


% regle de trois, mise a l'echelle
%fx      = 1 - erf(expo .* (1 - xin));
fx      = xin .^ expo;
morf    = (1 - fx) + fx .* pchip(thc,rhoc,thx) ./ pchip(thl,rhol,thx);
morf(:,1) = (morf(:,end) + morf(:,1)) ./ 2;
morf(:,end) = morf(:,1);
rhox = rhox .* morf;
Rx   = rhox .* cos(thx) + Raxe * vu;
Zx   = rhox .* sin(thx);


% coefficient de geometrie
function [vpr,grho2r2,grho,grho2,grhor,grho2b2,bpolm] = ggg(x,deuxd,fdia,regularisation,factor_two_pi)


% B2
BPHI = (fdia' * ones(1,size(deuxd.R,2))) ./ deuxd.R;
B2   = deuxd.BR .^ 2  + deuxd.BZ .^ 2 + BPHI .^ 2; 
%gg,
warning off
factinte = 2 .* pi .* deuxd.R .* deuxd.RZjac;
warning on
factinte(1,:) = 0;

% donnees finale
th       = deuxd.th(end,:);
vpr      = trapz(th,factinte,2)';
C2       = trapz(th,deuxd.gradrho .^2 ./ deuxd.R .^2 .* factinte,2)';
Cgrho    = trapz(th,deuxd.gradrho .* factinte,2)'; 
CgrhoR   = trapz(th,deuxd.gradrho ./ deuxd.R .* factinte,2)'; 
Cgrho2   = trapz(th,deuxd.gradrho .^2 .* factinte,2)'; 
Cgrho2b2 = trapz(th,deuxd.gradrho .^2 ./ B2 .* factinte,2)';
Cbpol2   = trapz(th,(deuxd.BR .^ 2  + deuxd.BZ .^ 2).* factinte,2)';

rmx(1) = 0;
vpr(1) = 0;
C2(1) = 0;
Cgrho(1) = 0;
Cgrho2(1) = 0;
CgrhoR(1) = 0;
Cgrho2b2(1) = 0;
Cbpol2(1)   = 0;

% securite
vpr   = max(0,vpr);

if regularisation == 0
  % donnees derivees
	grho2r2      = C2 ./ max(eps,vpr);
	grho2r2(1)   = pchip(x(2:end),grho2r2(2:end),0);
	grho         = Cgrho ./ max(eps,vpr);
	grho(1)      = pchip(x(2:end),grho(2:end),0);
	grhor        = CgrhoR ./ max(eps,vpr);
	grhor(1)     = pchip(x(2:end),grhor(2:end),0);
	grho2        = Cgrho2 ./ max(eps,vpr);
	grho2(1)     = pchip(x(2:end),grho2(2:end),0);
	grho2b2      = Cgrho2b2 ./ max(eps,vpr);
	grho2b2(1)   = pchip(x(2:end),grho2b2(2:end),0);
	bpolm        = sqrt(Cbpol2 ./ max(eps,vpr));
	bpolm(1)     = 0;
else
	%  % donnees derivees
	ll = max(2,find((deuxd.rho(:,1)./deuxd.rho(end,1)) > 0.05,1));
	grho2r2      = C2 ./ max(eps,vpr);
	grho2r2(1:ll)   = pchip(x(ll:end),grho2r2(ll:end),x(1:ll));
	grho         = Cgrho ./ max(eps,vpr);
	grho(1:ll)      = pchip(x(ll:end),grho(ll:end),x(1:ll));
	grhor        = CgrhoR ./ max(eps,vpr);
	grhor(1:ll)     = pchip(x(ll:end),grhor(ll:end),x(1:ll));
	grho2        = Cgrho2 ./ max(eps,vpr);
	grho2(1:ll)     = pchip(x(ll:end),grho2(ll:end),x(1:ll));
	grho2b2      = Cgrho2b2 ./ max(eps,vpr);
	grho2b2(1:ll)   = pchip(x(ll:end),grho2b2(ll:end),x(1:ll));
	bpolm        = sqrt(Cbpol2 ./ max(eps,vpr));
	bpolm(1:ll)     = 0;
end

% coefficient de geometrie
function [volume,surface,phi,dvdx,dsdx,dphidx] = volphi(x,deuxd,fdia,factor_two_pi)

%
dvdx      = zeros(size(x));
dsdx      = zeros(size(x));
dphidx    = zeros(size(x));
volume    = zeros(size(x));
surface   = zeros(size(x));
phi       = zeros(size(x));
%
BPHI = (fdia' * ones(1,size(deuxd.R,2))) ./ deuxd.R;
%  dv   = 2 .* pi .* deuxd.R .* deuxd.rhog;
%  ds   = deuxd.rhog;
%  dphi = BPHI .* deuxd.rhog;
dx   = mean(diff(x));

for k = 2:length(x)
	% dvdx(k)   = trapz(deuxd.th(end,:),(dv(k - 1,:)   + dv(k,:))   .* (deuxd.rhog(k,:) -deuxd.rhog(k-1,:)) ./ 2,2) ./ dx;
	% dsdx(k)   = trapz(deuxd.th(end,:),(ds(k - 1,:)   + ds(k,:))   .* (deuxd.rhog(k,:) -deuxd.rhog(k-1,:)) ./ 2,2) ./ dx;
	% dphidx(k) = trapz(deuxd.th(end,:),(dphi(k - 1,:) + dphi(k,:)) .* (deuxd.rhog(k,:) -deuxd.rhog(k-1,:)) ./ 2,2) ./ dx;
	dvdx(k)     = trapz(deuxd.th(end,:), pi .* (deuxd.R(k,:) + deuxd.R(k-1,:)) .* (deuxd.rhog(k,:) .^ 2 - deuxd.rhog(k-1,:) .^ 2) ./ 2,2) ./ dx;
	dsdx(k)     = trapz(deuxd.th(end,:), deuxd.rhog(k,:) .^ 2 - deuxd.rhog(k-1,:) .^ 2,2) ./ dx ./ 2;
	dphidx(k)   = trapz(deuxd.th(end,:),(BPHI(k - 1,:) + BPHI(k,:)) .* (deuxd.rhog(k,:).^ 2 -deuxd.rhog(k-1,:) .^ 2) ./ 4,2) ./ dx;
	volume(k)   = trapz(deuxd.th(end,:),pi .* deuxd.R(k,:) .* deuxd.rhog(k,:) .^ 2,2);
	surface(k)  = trapz(deuxd.th(end,:),deuxd.rhog(k,:) .^ 2 / 2,2);
end

%volume2   = dx .* cumsum(dvdx,2);
%surface2  = dx .* cumsum(dsdx,2);
phi        = dx .* cumsum(dphidx,2);

% ri et ri2 pour calculer F et q
function [r2i,ri,rmoy,r2,rmx,b2,b2i] = griri2(x,deuxd,fdia,factor_two_pi)

% B2
BPHI = (fdia' * ones(1,size(deuxd.R,2))) ./ deuxd.R;
B2   = deuxd.BR .^ 2  + deuxd.BZ .^ 2 + BPHI .^ 2; 
%gg,
warning off
factinte = 2 .* pi .* deuxd.R .* deuxd.RZjac;
warning on
factinte(1,:) = 0;

% donnees finale
th       = deuxd.th(end,:);
vpr      = trapz(th,factinte,2)';
C3       = trapz(th,factinte ./ deuxd.R .^ 2,2)';
C4       = trapz(th,factinte ./ deuxd.R,2)';
C5       = trapz(th,factinte .* deuxd.R,2)';
C6       = trapz(th,factinte .* deuxd.R .^ 2,2)';
Cb2      = trapz(th,B2.* factinte,2)';
Cb2i     = trapz(th,factinte./ B2,2)';
rmx      = sqrt(trapz(th,deuxd.rhog .^ 2,2)  ./ 2 ./ pi)';
		
rmx(1) = 0;
vpr(1) = 0;
C3(1) = 0;
C4(1) = 0;
C5(1) = 0;
C6(1) = 0;
Cb2(1)      = 0;
Cb2i(1)     = 0;

% securite
vpr = max(0,vpr);

% donnees derivees
r2i          = C3 ./ max(eps,vpr);
r2i(1)       = pchip(x(2:end),r2i(2:end),0);
ri           = C4 ./ max(eps,vpr);
ri(1)        = pchip(x(2:end),ri(2:end),0);
rmoy         = C5 ./ max(eps,vpr);
rmoy(1)      = pchip(x(2:end),rmoy(2:end),0);
r2           = C6 ./ max(eps,vpr);
r2(1)        = pchip(x(2:end),r2(2:end),0);
b2           = Cb2 ./ max(eps,vpr);
b2(1)        = pchip(x(2:end),b2(2:end),0);
b2i          = Cb2i ./ max(eps,vpr);
b2i(1)       = pchip(x(2:end),b2i(2:end),0);
		

% ri et ri2 pour calculer F et q
function jmoy = gjmoy(x,deuxd)

warning off
factinte = 2 .* pi .* deuxd.R .* deuxd.RZjac;
warning on
factinte(:,1) = 0;

% donnees finale
th       = deuxd.th(end,:);
vpr      = trapz(th,factinte,2)';
Cj       = trapz(th,deuxd.jphi .* factinte,2)';

vpr(1) = 0;
Cj(1) = 0;

% securite
vpr = max(0,vpr);

% donnees derivees
jmoy         = Cj ./ max(eps,vpr);
jmoy(1)      = pchip(x(2:end),jmoy(2:end),0);

