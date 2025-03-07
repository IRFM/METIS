% fonction pour calcul d'equilibre sans convergence
function  [profil,deuxd,moment,Ip_lcfs] = metis2equi1t_imas(profil,geo,phys,nbth,exposhape,factor_two_pi)

x = profil.x;
ve = ones(size(x));
% dpsidx s'annule au centre
psid1    = profil.psid1;
% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
psid2    = profil.psid2;
profil.deltax = profil.Raxe - profil.Raxe(end);
Raxe          = profil.Raxe;
epsi          = profil.epsi;

% calcul de la geometrie 2d
deuxd = ggeo2d(profil.x,profil.Raxe,geo.a,profil.kx,profil.dx,geo.vp,geo.sp,geo.Rsepa,geo.Zsepa,nbth,profil.psi,exposhape,factor_two_pi);
% calcul de gradient rho
deuxd = grho_comp(x,deuxd,profil.rho,factor_two_pi);
% moyenne geometrique pour le calcul de F et q
[profil.r2i_ctr,profil.ri_ctr,profil.rmoy,profil.r2,profil.rmx_ctr,profil.b2,profil.b2i] = griri2(x,deuxd,profil.fdia,factor_two_pi);
% C3 (compatibilite ancienne version
profil.C3_ctr  = profil.vpr .* profil.r2i_ctr;
% BPHI
deuxd.BPHI = (profil.fdia'*ones(1,size(deuxd.R,2)))./deuxd.R;
%
deuxd.drhodpsi       = pdederive(profil.psi,profil.rho,0,2,2,1)' * ones(1,size(deuxd.R,2));
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
[profil.vpr_tor_ctr,profil.grho2r2_ctr,profil.grho_ctr,profil.grho2_ctr, ...
 profil.grhor,profil.grho2b2,profil.bpolm] = ggg(profil.x,deuxd,profil.fdia,factor_two_pi);
% calcul de C2
profil.C2_ctr   = profil.grho2r2_ctr .* profil.vpr;
profil.vpr_tor_ctr = profil.vpr_tor_ctr ./ pdederive(x,profil.rho,2,2,2,1);
% calcul de li
profil.bpol_ctr  = abs(sqrt(profil.grho2r2_ctr) .* abs(pdederive(profil.rho,profil.psi,0,2,2,1)) ./ factor_two_pi);    
  
		     
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
% calcul de dthetadR et dZ
deuxd.dthdR = - deuxd.dZdx ./ max(eps,deuxd.RZjac);
deuxd.dthdR(1,:) = 0;
deuxd.dthdZ =  deuxd.dRdx ./ max(eps,deuxd.RZjac);
deuxd.dthdZ(1,:) = 0;

% computing circulation of B on LCFS for diagnostic purpose
for k=1:size(deuxd.R,1)
  Ip_lcfs(k) = sign(factor_two_pi) .* sum(((deuxd.BR(k,1:end-1) + deuxd.BR(k,2:end)) .* diff(deuxd.R(k,:)) ./ 2) + ...
          ((deuxd.BZ(k,1:end-1) + deuxd.BZ(k,2:end)) .* diff(deuxd.Z(k,:)) ./ 2)) ./ (4.*pi.*1e-7);
end
%  figure(118);clf
%  plot(deuxd.R(end,:),deuxd.Z(end,:),'r');
%  hold on
%  quiver(deuxd.R(end,:),deuxd.Z(end,:),deuxd.BR(end,:),deuxd.BZ(end,:),'b')
%  quiver((deuxd.R(end,1:end-1) + deuxd.R(end,2:end)) ./2, ...
%         (deuxd.Z(end,1:end-1) + deuxd.Z(end,2:end)) ./2, ...
%         (deuxd.BR(end,1:end-1) + deuxd.BR(end,2:end)) ./2, ...
%         (deuxd.BZ(end,1:end-1) + deuxd.BZ(end,2:end)) ./2,'k');
%  drawnow

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
if nbth < 0
    th   = linspace(2.*pi,0,abs(nbth));
else
    th   = linspace(0,2.*pi,nbth);
end
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
deuxd.BR       = NaN .* ones(length(x),length(th));
deuxd.BZ       = NaN .* ones(length(x),length(th));
deuxd.RZjac    = NaN .* ones(length(x),length(th));
deuxd.dPSIdR   = NaN .* ones(length(x),length(th)); 
deuxd.dPSIdZ   = NaN .* ones(length(x),length(th)); 

% boucle pour la correction de forme
for l = 1:length(x)

	R = Raxe(:,l) * vh + (amin(:,l) * vh) .* (x(:,l) * vh) .*  ...
		cos(th + (t(:,l) * vh) .* sin(th));
	Z = (kx(:,l) * vh) .*(amin(:,l) * vh) .* (x(:,l) * vh) .* sin(th);
	% morphing
	if ~isempty(Rsepa) && ~isempty(Zsepa) && (expo > 0)
		[R,Z]   = zemorph(R,Z,Rsepa,Zsepa,Rmom,Zmom,x(1,l),expo);
	end
	Raxel       = 0.5 .* (min(R,[],2) + max(R,[],2));
	rho2    = (R - Raxel * vh) .^2 +  Z .^ 2;

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
	rhoo   = spline(thx,rhox,th);
	deuxd.cs(l,:)    = cos(th);
	deuxd.R(l,:)     = Raxel + rhoo .* deuxd.cs(l,:);
	deuxd.ss(l,:)    = sin(th);
	deuxd.Z(l,:)     = rhoo .* deuxd.ss(l,:);
	deuxd.rhog(l,:)  = rhoo;
		
end

% la coordonnee inetrne est a*x
deuxd.dRdx           =  0.5 .* (pdederive(x',deuxd.R,1,2,1,1) + pdederive(x',deuxd.R,1,1,1,1));
deuxd.dZdx           =  0.5 .* (pdederive(x',deuxd.Z,1,2,1,1) + pdederive(x',deuxd.Z,1,1,1,1));
deuxd.dPSIdx         =  0.5 .* (pdederive(x',deuxd.PSI,0,2,1,1) +  pdederive(x',deuxd.PSI,0,1,1,1));
deuxd.dRdth          =  pdederive(th,deuxd.R,2,2,2,1);
deuxd.dRdth(:,1)     = 0.5 .* (deuxd.dRdth(:,1) + deuxd.dRdth(:,end));
deuxd.dRdth(:,end)   = deuxd.dRdth(:,1);
deuxd.dZdth          =  pdederive(th,deuxd.Z,2,2,2,1);
deuxd.dZdth(:,1)     = 0.5 .* (deuxd.dZdth(:,1) + deuxd.dZdth(:,end));
deuxd.dZdth(:,end)   = deuxd.dZdth(:,1);
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
deuxd.dPSIdZ   = -deuxd.dPSIdx .* deuxd.dRdth  ./ max(eps,deuxd.RZjac);
% champ poloidal
deuxd.BR      =    deuxd.dPSIdZ ./ deuxd.R ./ abs(factor_two_pi);
deuxd.BZ      =  - deuxd.dPSIdR ./ deuxd.R ./ abs(factor_two_pi);
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
deuxd.drhodx              =  0.5 .* (pdederive(x',deuxd.rho,1,2,1,1) + pdederive(x',deuxd.rho,1,1,1,1));
deuxd.drhodth             =  pdederive(th,deuxd.rho,2,2,2,1);
deuxd.drhodth(:,1)        = 0.5 .* (deuxd.drhodth(:,1) +deuxd.drhodth(:,end));
deuxd.drhodth(:,end)      = deuxd.drhodth(:,1);

% calcul de gradrho
deuxd.drhodR   = (deuxd.drhodx .* deuxd.dZdth -  deuxd.drhodth .* deuxd.dZdx)  ./ max(eps,deuxd.RZjac);
deuxd.drhodZ   = (-deuxd.drhodx .* deuxd.dRdth + deuxd.drhodth .* deuxd.dRdx)  ./ max(eps,deuxd.RZjac);
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
function [vpr,grho2r2,grho,grho2,grhor,grho2b2,bpolm] = ggg(x,deuxd,fdia,factor_two_pi)


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
% added abs to protect against orientation of LCFS points
vpr      = abs(trapz(th,factinte,2))';
C2       = abs(trapz(th,deuxd.gradrho .^2 ./ deuxd.R .^2 .* factinte,2))';
Cgrho    = abs(trapz(th,deuxd.gradrho .* factinte,2))'; 
CgrhoR   = abs(trapz(th,deuxd.gradrho ./ deuxd.R .* factinte,2))'; 
Cgrho2   = abs(trapz(th,deuxd.gradrho .^2 .* factinte,2))'; 
Cgrho2b2 = abs(trapz(th,deuxd.gradrho .^2 ./ B2 .* factinte,2))';
Cbpol2   = abs(trapz(th,(deuxd.BR .^ 2  + deuxd.BZ .^ 2).* factinte,2))';
rmx      = sqrt(abs(trapz(th,deuxd.rhog .^ 2,2))  ./ 2 ./ pi)';

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
% added abs to protect against orientation of LCFS points
vpr      = abs(trapz(th,factinte,2))';
C3       = abs(trapz(th,factinte ./ deuxd.R .^ 2,2))';
C4       = abs(trapz(th,factinte ./ deuxd.R,2))';
C5       = abs(trapz(th,factinte .* deuxd.R,2))';
C6       = abs(trapz(th,factinte .* deuxd.R .^ 2,2))';
Cb2      = abs(trapz(th,B2.* factinte,2))';
Cb2i     = abs(trapz(th,factinte./ B2,2))';
rmx      = sqrt(abs(trapz(th,deuxd.rhog .^ 2,2))  ./ 2 ./ pi)';
		
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

% security 
b2  = max(0,b2);
b2i = max(0,b2i);

