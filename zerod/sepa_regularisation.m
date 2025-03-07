% regularisation de la separatrice par fit d'une solution analytique dans le vide de GS
function [Rout,Zout,coef] = sepa_regularisation(Rin,Zin,dpsi,nbout,weight_x,plotonoff)

if nargin < 3
	dpsi = 1;
end
if nargin < 4
	nbout =201;
end
if nargin < 5
	weight_x = 5;
end
if nargin < 6
	plotonoff = 0;
end

% mise en forme
smem = size(Rin);
Rin = Rin(:);
Zin = Zin(:); 
% ajout de l'axe

R = cat(1,Rin,mean(Rin));
Z = cat(1,Zin,mean(Zin));


% fit des coefficients
mpol = NaN * ones(length(R),16);
for k = 1:16
  coef = zeros(16);
  coef(k) = 1;
  [psi_out,BR_out,BZ_out] = GS_vide_updown_16(coef,R,Z);
  mpol(:,k) = psi_out;
end
ypol = cat(1,-dpsi * ones(size(Rin)),0);
% augmentation du poid de points pres du point X
indx= cat(1,find(Zin == max(Zin),1),find(Zin == min(Zin),1));
mpol(indx,:) = weight_x .* mpol(indx,:);
ypol(indx)   = weight_x .* ypol(indx);
coef =  mpol \ ypol;


% mise en forme pour sortie uniforme
zerod_void.d0 = (max(Rin) - min(Rin)) / 20;
[Rsepa,Zsepa] = z0reshape_sepa(fillgeo(Rin,Zin),zerod_void,Rin',Zin',5,0.95,nbout);


% test erreur initiale
[psi_sepa,BR_sepa,BZ_sepa] = GS_vide_updown_16(coef,Rsepa,Zsepa); 
errpsi  =  sum((psi_sepa + dpsi) .^ 2);

% passage en complexe (plus compact)
p_sepa_rf   = Rsepa + sqrt(-1) .* Zsepa;
% initilisation du Newton
dpsidr_sepa  =   real(p_sepa_rf) .* BZ_sepa;
dpsidz_sepa  = - real(p_sepa_rf) .* BR_sepa; 
grad_psi     = dpsidr_sepa + sqrt(-1) .* dpsidz_sepa;


errpsi_mem = 0;
nbloop     = 100;  % nombre maximal de loop
nstep      = 1 ./ sqrt(2);  % pas du Newton

% boucle de convergence
while (errpsi > (eps .* dpsi)) && (nbloop > 0)
        % correction de la position des points
	p_sepa_rf = p_sepa_rf + nstep .* grad_psi ./ abs(grad_psi) .^ 2 .* (psi_sepa + dpsi);
        % nb
	nbloop = nbloop - 1;

        % mise a jour 
        [psi_sepa,BR_sepa,BZ_sepa] = GS_vide_updown_16(coef,real(p_sepa_rf),imag(p_sepa_rf)); 
	dpsidr_sepa  =   real(p_sepa_rf) .* BZ_sepa;
	dpsidz_sepa  = - real(p_sepa_rf) .* BR_sepa;
	grad_psi     = dpsidr_sepa + sqrt(-1) .* dpsidz_sepa;
	errpsi  =  sum((psi_sepa + dpsi) .^ 2);

end

Rout	= real(p_sepa_rf);
Zout	= imag(p_sepa_rf);


if plotonoff == 1
	% moment de la separatrice
	a = (max(Rin) - min(Rin)) / 2;
	b = (max(Zin) - min(Zin)) / 2;
	% calcul sur une grille pour retrouver la nouvelle separatrice
	nbeqdsk = 101;
	% (R,Z) GRID
	r = linspace(min(Rin)- a ./ 4 ,max(Rin) + a ./ 4,nbeqdsk);
	z = linspace(min(Zin)- b ./ 4 ,max(Zin) + b  ./ 4 ,ceil(nbeqdsk .* max(1,b/a)));
	[r2d,z2d] = meshgrid(r,z);
	[psi2d,BR2d,BZ2d] = GS_vide_updown_16(coef,r2d,z2d); 
	
	h = findobj(0,'type','figure','tag','sepa_regularisation');
	if isempty(h)
	h=figure('tag','sepa_regularisation');
	else
	figure(h);
	end   
	clf
	set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
		'defaultlinelinewidth',1,'color',[1 1 1])
	contour(r2d,z2d,psi2d,101)
	hold on
	plot(Rin,Zin,'r.')
	contour(r2d,z2d,psi2d,[-dpsi , -dpsi],'g')
	plot(Rout,Zout,'k')
	drawnow
end
if nbloop <= 0
	figure(203)
	clf
	plot(real(p_sepa_rf),imag(p_sepa_rf),'b',Rsepa,Zsepa,'.r')
	hold on
	quiver(real(p_sepa_rf),imag(p_sepa_rf),dpsidr_sepa ,dpsidz_sepa);
	drawnow
        errpsi
        nb = 100 - nbloop
end

% fonction appele par zfitvide
function [psi_out,BR_out,BZ_out] = GS_vide_updown_16(coef,Rin,Zin)


% separation
s1   = coef(1);
s2   = coef(2);
A    = coef(3);
B    = coef(4);
cn   = coef(5:16);

% dimension
nrz = size(Rin);
R  = Rin(:);
Z  = Zin(:);

% init (la partie en a + b * r^2) est dans le developpment
% solution trivial de Psi
rr   = sqrt(R .^ 2 + Z .^ 2);
psiv = s1 .* rr;
BRv  = - s1 .* Z ./ rr ./ R;
BZv  = s1  ./ rr;

% nouvelle contribution (groupe de symetrie de l'equation)
% ref : Y.E. Litvinenko  PoP 17, 074502 (2010). 
psiv  = psiv + s2 .* (Z .* sqrt(R .^ 2 + Z .^ 2) + R .^ 2 .* asinh(Z./R));
BRv   = BRv  - s2 .* 2 .* sqrt(Z .^ 2 + R .^ 2) ./ R;
BZv   = BZv  + s2 .* 2.* asinh(Z ./ R);


% calcul (partie courant)
% solution soloviev
% ref : P. J. Mc Carthy POP 1999 p 3554-...
psiv = psiv +  A .* R .^ 4 ./ 8 + B .* Z .^2 ./ 2;
BRv  = BRv  - B .* Z ./ R; 
BZv  = BZv +  A .* R .^ 2 ./ 2;


% calcul (partie contribution du vide)
% ref :
% Analytical tokamak equilibrium for shaped plasmas
% S. B. Zheng,a) A. J. Wootton, and Emilia R. Solano
% Phys. Plasmas 3 (3), March 1996

% A.J. Cerfon and J.P. Freidberg, POP 17, 032502 (2010)
% choix R0 = 1
% 1
psi_1 = 1;
BR_1  = 0;
BZ_1  = 0;

% 2 
psi_2 = R .^ 2; 
BR_2  = 0;
BZ_2  = 2;

% 3
psi_3 = Z .^ 2  - R .^ 2 .* log(R);
BR_3  = - 2 .* Z ./ R;
BZ_3  = - 2 .* log(R) - 1;

%4 
psi_4 = R .^ 4 - 4 .* R .^ 2 .* Z .^ 2;
BR_4  = 8 .* R .* Z;
BZ_4  = 4 .* R .^ 2 - 8 .* Z .^ 2;

%5
psi_5 = 2 .* Z .^ 4  - 9 .* R .^ 2 .* Z .^ 2  + 3 .* R .^ 4  .* log(R)  - 12 .* R .^ 2 .* Z .^ 2 .* log(R);
BR_5  = - (8 .* Z .^ 3 + (-24 .* R .^ 2 .* log(R) - 18 .* R .^ 2) .* Z) ./ R;
BZ_5  =  -(24 .* log(R) + 30) .* Z .^ 2 + 12 .* R .^ 2 .* log(R) + 3 .* R .^ 2;

%6
psi_6 = R .^ 6  - 12 .* Z .^ 2 .* R .^ 4 + 8 .* Z .^ 4 .* R .^ 2;
BR_6  = -(32 .* R .* Z .^ 3 - 24 .* R .^ 3 .* Z);
BZ_6  = 16 .* Z .^ 4 - 48 .* R .^ 2 .* Z .^ 2 + 6 .* R .^ 4;

%
psi_7 = 8 .* Z .^ 6  - 140 .* Z .^ 4 .* R .^ 2 + 75 .* Z .^ 2 .* R .^ 4  - 15  .* R .^ 6 .* log(R) + 180 .* R .^ 4 .* Z .* 2 .* log(R) - 120 .* R .^ 2 .* Z .^ 4 .* log(R); 
BR_7  = -(48 .* Z .^ 5 + (-480 .* R .^ 2 .* log(R) - 560 .* R .^ 2) .* Z .^ 3 + 150 .* R .^ 4 .* Z + 360 .* R .^ 4 .* log(R)) ./ R;
BZ_7  = (-240 .* log(R) - 400) .* Z .^ 4 + 300 .* R .^ 2 .* Z .^ 2 + (1440 .* R .^ 2 .* log(R) + 360 .* R .^ 2) .* Z - 90 .* R .^ 4 .* log(R) - 15 .* R .^ 4;


% 8
psi_8 = Z;
BR_8  = - 1./ R;
BZ_8  = 0;

% 9
psi_9 = Z .* R .^ 2;
BR_9  = - R;
BZ_9  = 2.* Z;

% 10
psi_10 = Z .^ 3 - 3 .* Z .* R .^ 2 .* log(R);
BR_10  = -(3 .* Z .^ 2 - 3 .* R .^ 2 .* log(R)) ./ R;
BZ_10  = (-6 .* log(R) - 3) .* Z;

% 11
psi_11 = 3 .* Z .* R .^ 4  - 4 .* Z .^ 3 .* R .^ 2;
BR_11  = - 3 .* R .^ 3 + 12 .* R .* Z .^ 2;
BZ_11 = 12 .* R .^ 2 .* Z - 8 .* Z .^ 3;

%12
psi_12 = 8 .* Z .^ 5 - 45 .* Z .* R .^ 4 - 80 .* Z .^ 3 .* R .^ 2 .* log(R) + 60 .* Z .* R .^ 4 .* log(R); 
BR_12  = - (40 .* Z .^ 4 - 240 .* R .^ 2 .* log(R) .* Z .^ 2 + 60 .* R .^ 4 .* log(R) - 45 .* R .^ 4) ./ R;
BZ_12  = (-160 .* log(R) - 80) .* Z .^ 3 + (240 .* R .^ 2 .* log(R) - 120 .* R .^ 2) .* Z;

% sommation
psiv = psiv + cn(1) .* psi_1 + cn(2)  .* psi_2 +  cn(3)  .* psi_3  + cn(4)  .* psi_4 + ...
              cn(5) .* psi_5 + cn(6)  .* psi_6 +  cn(7)  .* psi_7  + cn(8)  .* psi_8 + ...
              cn(9) .* psi_9 + cn(10) .* psi_10 + cn(11) .* psi_11 + cn(12) .* psi_12;

BRv = BRv + cn(1) .* BR_1 + cn(2)  .* BR_2 +  cn(3)  .* BR_3  + cn(4)  .* BR_4 + ...
            cn(5) .* BR_5 + cn(6)  .* BR_6 +  cn(7)  .* BR_7  + cn(8)  .* BR_8 + ...
            cn(9) .* BR_9 + cn(10) .* BR_10 + cn(11) .* BR_11 + cn(12) .* BR_12;

BZv = BZv + cn(1) .* BZ_1 + cn(2)  .* BZ_2 +  cn(3)  .* BZ_3  + cn(4)  .* BZ_4 + ...
            cn(5) .* BZ_5 + cn(6)  .* BZ_6 +  cn(7)  .* BZ_7  + cn(8)  .* BZ_8 + ...
            cn(9) .* BZ_9 + cn(10) .* BZ_10 + cn(11) .* BZ_11 + cn(12) .* BZ_12;


% probleme de convention de signe 
BRv = - BRv;
BZv = - BZv;

% mise en forme
psi_out   = reshape(psiv,nrz);
BR_out    = reshape(BRv,nrz);
BZ_out    = reshape(BZv,nrz);


function geo = fillgeo(Rsepa,Zsepa)

if size(Rsepa,2) == 1
	Rsepa = Rsepa';
	Zsepa = Zsepa';
end

geo.z0 = (max(Zsepa) + min(Zsepa)) ./ 2;
Zsepa  = Zsepa - geo.z0;

% calcul des moments pour assurer la coherence dans tous les cas avec les moments
% centre pour angle d'integration
rc = mean(Rsepa,2);
zc = mean(Zsepa,2);
vc = ones(1,size(Rsepa,2));
%uc = unwrap(angle((Rsepa-rc*vc) + sqrt(-1) .* (Zsepa  -zc*vc)));
uc = unwrap(angle((Rsepa-rc*vc) + sqrt(-1) .* Zsepa ));
uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
uc(:,1)   = uc(:,end) + 2 .* pi;
xu    = linspace(0,1,length(vc));
dRdx  = pdederive(xu,Rsepa,2,2,2,1);
dZdx  = pdederive(xu,Zsepa,2,2,2,1);
% calcul de R0 et Z0
maskrmax  = (Rsepa == (max(Rsepa,[],2) * vc));
% recalcul des parametres sur le vecteur final
rmin  = min(Rsepa,[],2);
rmax  = max(Rsepa,[],2);
geo.a = 0.5 .* (rmax - rmin);
geo.R = 0.5 .* (rmax + rmin);
zmin  = min(Zsepa,[],2);
zmax  = max(Zsepa,[],2);
% geo.z0   = (zmax + zmin + sum(Zsepa .* maskrmax,2) ./ sum(maskrmax,2)) ./ 3;
geo.K    = (abs(trapz(xu,Zsepa .*  dRdx,2) ./ pi ./ geo.a) + (zmax - zmin)) ./ 3 ./ geo.a;

rzmax = geo.R;
rzmin = geo.R;
for k = 1:size(Zsepa,1)
	rzmax(k) = Rsepa(k,min(find(Zsepa(k,:) == zmax(k))));
	rzmin(k) = Rsepa(k,min(find(Zsepa(k,:) == zmin(k))));
end
uu   =  angle(rzmax - geo.R + sqrt(-1) .* zmax);
ul   =  angle(rzmin - geo.R + sqrt(-1) .* zmin);
tu   =  abs((acos((rzmax - geo.R) ./ geo.a) - acos(cos(uu))) ./ sin(uu));
tl   =  abs((acos((rzmin - geo.R) ./ geo.a) - acos(cos(ul))) ./ sin(ul));
tm   =  (tl + tu) ./ 2;
d    =   abs(rzmax + rzmin -  2 .* geo.R) ./ 2 ./ geo.a;
geo.d =  0.6 .* d + 0.4  .* sin(tm);

	      