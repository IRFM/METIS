% modele de coefficients de transport CDBM
function [kie,kie_itb,kii,kii_itb,chie,chii,d,vn] = z0cdbm(meff,profil)

% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)


% fonction d'amelioration de metis
frs  = max(0.1,min(10,profil.xieshape_itb ./ max(eps,profil.xieshape)));

% calcul du  parametre alpha
b2m    = profil.bpol .^ 2  + profil.fdia .^ 2 .* profil.r2i;
beta   = profil.ptot ./ (b2m ./ 2 ./ phys.mu0);
gbeta  = pdederive(profil.xli,beta,0,2,2,1) ./ max(eps,pdederive(profil.xli,profil.rmx,0,2,2,1)) .* profil.grho;
alpha  = - profil.qjli .^ 2 ./ profil.ri .* gbeta;
alpha  = max(-2,min(2,alpha));
% magnetic shear
s = pdederive(profil.xli,profil.qjli,0,2,2,1) ./ max(eps,pdederive(profil.xli,profil.rmx,0,2,2,1)) .* profil.rmx ./ profil.qjli;
% for this model
s = max(-0.5,min(1.5,s));
%
c = phys.c;
%
q  = profil.qjli;
ve = ones(size(profil.xli));
R  = profil.Raxe(:,end) * ve;
a  = (profil.Raxe(:,end) .* profil.epsi(:,end));
ax = a * profil.xli;
a  = a * ve;
b0 = (profil.fdia(:,end) ./ profil.Raxe(:,end)) * ve;
%
pth          = phys.e .* (profil.nep .* profil.tep + profil.nip .* profil.tip);
betathermal  = 2 .* phys.mu0 .* pth ./ b0 .^ 2;
betat        = 2 .* phys.mu0 .* profil.ptot ./ b0 .^ 2;
alphathermal = R .* (-profil.qjli .^2 ) .* pdederive(ax,betathermal,0,2,2,1);
%
we   = ((profil.nep .* phys.e .^ 2) ./ (phys.me .* phys.epsi0)) .^ 0.5;
bphi = profil.fdia ./ profil.Raxe;
vA2  = (bphi .^ 2) ./ (phys.mu0 .* phys.mp .* profil.nep);
vA   = vA2 .^ 0.5;
%
che = (12 .* ((abs(alphathermal)) .^ 1.5) .* vA * c .^ 2) ./ ((we .^ 2) .* q.* R);
che = 0.1 + 30 * tanh(che/30);
s2  = s - alpha;
s2(alpha < 0) = -s2(alpha < 0);
% condition s2 << 1/2
s2  = min(0.45,s2);
f   = ones(size(q));
neg = find(s2<0);
pos = find(s2>0);
f(neg) = ( 2 .* (1 - 2 .* s2(neg)) .* (1 - 2 .* s2(neg) + 3 .* s2(neg) .^ 2)) .^ (-0.5);
f(pos) = (( 2 .^ 0.5 .* (1 - 2 .* s2(pos) + 3 .* s2(pos) .^ 2 + 2 .* s2(pos) .^ 3)) .^ (-1)) .* (1 + 9 .* 2 .^ 0.5 .* s2(pos) .^ 2.5);
% neglected in inital paper
% kim = - profil.epsi.* (1 - 1 ./ profil.qjli .^ 2);
% fim = (kim ./ profil.epsi) .^ (3/2) ./ s .^ 2;
% fim(kim <=0) = 0;
% f   = max(f,fim);
% %
elong = profil.kx;
kappa =(( 2 .* elong .^ 0.5) ./ (elong .^ 2 + 1)) .^ 1.5;
chie  = che .* f .* kappa;
chii  = che .* f .* kappa;
d     = 0.2 .* chie;
cq    = 0.25 + (2/3) .* s .* (s >= 0);
vq    = d .* 2 .* cq .* profil.ri;

% thermodiffusion (hypothese thermodiffusion quasi nulle)
% les gradients 
gte = pdederive(profil.xli,profil.tep,0,2,2,1) ./ max(eps,pdederive(profil.xli,profil.rmx,0,2,2,1)) .* profil.grho;
gti = pdederive(profil.xli,profil.tip,0,2,2,1) ./ max(eps,pdederive(profil.xli,profil.rmx,0,2,2,1)) .* profil.grho;
% borne pour les longueurs de gradient
% maximum
lmax = 2./ profil.ri;
% minimum
lmin = (profil.rmx(:,end) ./ length(profil.xli)) * ve;
% longueur de gradient
lte = lcalc(profil.tep,gte,lmin,lmax);
lti = lcalc(profil.tip,gti,lmin,lmax);
vt   = 0.1 .* (chie ./ lte - chie ./ lti);
% sommation
vn   = (vq + vt);
ind = find(~isfinite(vn));
if ~isempty(ind)
    vn(ind)    = zeros(1,length(ind));
end

% donnees pour la forme du profil
xie      = chie;
xii      = chii;
%xii_itb  =  frs .* xii;
%xie_itb  =  frs .* xie;
% ITB is already included in the model
xii_itb  =  1 .* xii;
xie_itb  =  1 .* xie;
% valeur centrale
xie(:,1)     = xie(:,2);
xii(:,1)     = xii(:,2);
xie_itb(:,1) = xie_itb(:,2);
xii_itb(:,1) = xii_itb(:,2);

% sorties en K
ne_shape   = profil.nep ./ max(1,max(profil.nep,[],2) * ve);
ni_shape   = profil.nip ./ max(1,max(profil.nip,[],2) * ve);
kie        = profil.xieshape .* xie     .* ne_shape;
kie_itb    = profil.xieshape .* xie_itb .* ne_shape; 
kii        = profil.xieshape .* xii     .* ni_shape;
kii_itb    = profil.xieshape .* xii_itb .* ni_shape;

return
figure(123);
subplot(3,4,1)
plot(profil.xli,profil.xieshape)
title('xieshape');
subplot(3,4,2)
plot(profil.xli,xie)
title('xie');
subplot(3,4,3)
plot(profil.xli,ne_shape)
title('ne_shape');
subplot(3,4,4)
plot(profil.xli,kie)
title('kie');
subplot(3,4,5)
plot(profil.xli,che)
title('che');
subplot(3,4,6)
plot(profil.xli,f)
title('f');
subplot(3,4,7)
plot(profil.xli,kappa)
title('kappa');
subplot(3,4,8)
plot(profil.xli,s2)
title('s - alpha');
subplot(3,4,9)
plot(profil.xli,alpha)
title('alpha');
subplot(3,4,10)
plot(profil.xli,s)
title('s');
subplot(3,4,11)
plot(profil.xli,gbeta)
title('gbeta');
% subplot(3,4,11)
% plot(profil.xli,fim)
% title('fim');
drawnow


% calcul de la longueur de gradient
function out = lcalc(v,gv,lmin,lmax)

ind = find( abs(gv) < eps);
if ~isempty(ind)
	 gv(ind) =eps;
end
out = abs(v ./ gv);
if ~isempty(ind)
	 out(ind) = lmax(ind);
end
out = min(lmax,max(lmin,out));
