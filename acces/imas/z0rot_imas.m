% dans metis la direction toroidal est dans le sens du courant plasma, 
% de telle sorte que Btheta est toujours positif
% le moment injecte par l'IDN est posistif si l'injection est co courant
% le champs toroidal est positif s'il est dans le sens du courant
% calcul de la vitesse de rotation moyenne
% cette fonction calcul la rotation toroidal pour chaque especes 
% ainsi que le rotation poloidal
function [rtor,vtor,vpol,omega,mtor] = z0rot_imas(zs,profli,option,frhe3,geo,cons)

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

% calcul du profil de rotation toroidal
% palsma de fond
% switch option.gaz
% case 1
%    zj = 1;
%    aj = 1;
% case 2
%    zj = 1;
%    aj = 2;
% case 3
%    zj = 1;
%    aj = mean(2 .* (1 - cons.iso) +  3 .* cons.iso);
% case 4
%    zj = 2;
%    aj = 4;
% end

% isotopic composition for option.gaz == 5
if option.gaz == 5
    nHe3onD = real(cons.iso);
    nTonD   = imag(cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(cons.iso));
    nTonD   = real(cons.iso);
end
cons.iso = real(cons.iso);

% case of Hydrogen NBI in DT plasma
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
   gas_nbi = -1; 
else
   gas_nbi = option.gaz; 
end

% Sn
if ~isfield(option,'Sn_fraction')
    option.Sn_fraction = 0;
end

% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - option.zimp);
mask = (dd == min(dd));
aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(aimp)
    aimp = 7/3 .* option.zimp;
end
dd   = abs(Z_el - option.zmax);
mask = (dd == min(dd));
amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(amax)
    amax = 7/3 .* option.zmax;
end
% % impurete principale
zimp = option.zimp;
% aimp = ceil(zimp .* (7/3));
% 
% % 2ieme impurete
zmax = option.zmax;
% amax = ceil(zmax .* (7/3));


% pour chaque espece d'ions
x     = profli.xli;
ve    = ones(size(x));
vt    = ones(size(profli.n1p,1),1);
nDm   = interp1_imas(zs.temps,zs.nDm,profli.temps,'pchip','extrap');
nTm   = interp1_imas(zs.temps,zs.nTm,profli.temps,'pchip','extrap');
nDp   = max(1e13,profli.n1p .* ((nDm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
nTp   = max(1e13,profli.n1p .* ((nTm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
switch option.gaz
    case 11
        nbp   = nTp;
        nTp   = zeros(size(nTp));
        nHp   = max(1e13,profli.n1p - nDp);
    otherwise
        nbp   = zeros(size(nTp));
        nHp   = max(1e13,profli.n1p - nTp - nDp);
end
switch option.gaz
    case 5
        nhep3  = profli.nhep;
        nhep   = option.frhe0 .* profli.nep;
    otherwise
        nhep  = max(1e13,profli.nhep .* (1 - frhe3));
        nhep3 = max(1e13,profli.nhep .* frhe3);
end
nz1p  = max(1e13,profli.nzp);
nz2p  = max(1e11,profli.nzp .* option.rimp);  
nwp   = max(1, profli.nwp);
% masse
if option.Sn_fraction > 0
    Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p +  ...
              (1 - option.Sn_fraction) .*183.84 .* nwp) + option.Sn_fraction .*118.71 .* nwp + 3.02 .* nhep3 + 11 .* nbp;
else
    Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p + 183.84 .* nwp + 3.02 .* nhep3 + 11 .* nbp);
end
% omega homothetic a Ti (cas Jet avec NBI)
% calcul du moment (definition de cronos sum(n_i*m_i*<R*V_phi>)
rtor = Mtor .* profli.omega ./ profli.r2i;
% calcul de la rottaion poloidal 
warning off


% formulaire ORNL
lnii    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(1); 
lnhe    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(2); 
lnhe3   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(2); 
lnz1    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zimp); 
lnz2    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zmax); 
lnw     = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(z0wavez(profli.tep)); 
lnb    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(5);

% pour l'espece principale
% Tokamaks, Wesson p 663
taui_h     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 1)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nHp ./ lnii ./ 1 .^ 4 ;
taui_d     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 2)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nDp ./ lnii ./ 1 .^ 4 ;
taui_t     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 3)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nTp ./ lnii ./ 1 .^ 4 ;
taui_he     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 4)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nhep ./ lnhe ./ 2 .^ 4 ;
taui_he3     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 3)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nhep ./ lnhe3 ./ 2 .^ 4 ;
taui_z1     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* aimp)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nz1p ./ lnz1 ./ zimp .^ 4 ;
taui_z2     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* amax)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nz2p ./ lnz2 ./ zmax .^ 4 ;
taui_w      =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 183.84)  .* ...
	                    (phys.e .* profli.tip) .^ (3/2) ./ nwp ./ lnw ./ max(1,z0wavez(profli.tep)) .^ 4 ; 
                    
taui_b      =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 11)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ nbp ./ lnb ./ 5 .^ 4 ;
                   
if option.Sn_fraction > 0
    lnsn    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(z0snavez(profli.tep));
    taui_sn     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 118.71)  .* ...
        (phys.e .* profli.tip) .^ (3/2) ./ (option.Sn_fraction .* nwp) ./ lnsn ./ max(1,z0snavez(profli.tep)) .^ 4 ;
    if option.Sn_fraction < 1
        taui_w      =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 183.84)  .* ...
            (phys.e .* profli.tip) .^ (3/2) ./ ((1 - option.Sn_fraction) .* nwp) ./ lnw ./ max(1,z0wavez(profli.tep)) .^ 4 ;
    end
end




% Plasma rotation in Tokamaks, V. Rozhansky and M. Tendler,Review of plasma physics, tome 19, p 163
vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 1 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_h;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkh     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 2 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi./ taui_d;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkd     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 3 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_t;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkt     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 4 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_he;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkhe     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 3 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_he;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkhe3     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ aimp ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z1;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkz1     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ amax ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z1;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkz2     = max(-2.1,min(1.7,fk));

if option.Sn_fraction  == 1
    vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 183.84 ./ phys.mp);
    nuis   = zeros(size(profli.Raxe));
else
    vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 183.84 ./ phys.mp);
    nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_w;
    nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
end
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
    (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkw     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 11 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_b;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
    (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkb    = max(-2.1,min(1.7,fk));
fkb(~isfinite(fkb)) = 0;

if option.Sn_fraction > 0
    vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 118.71 ./ phys.mp);
    nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_sn;
    nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
    fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
        (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
    fksn     = max(-2.1,min(1.7,fk));
end


warning on


% changement de repere
dpsidrho = abs(pdederive(x,profli.psi,0,2,2,1) ./ (profli.rmx(:,end) * ve));
dpsidrho(:,1) = NaN;

% vitesse poloidal 
b2 = (profli.fdia .* profli.ri) .^ 2; 
% standard  + orbit squeezing
gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .*nHp,0,2,2,1) ./ nHp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_h  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_h(:,1) = 2 .* utheta_h(:,2) - utheta_h(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nDp,0,2,2,1) ./ nDp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_d  = - fkd .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_d(:,1) = 2 .* utheta_d(:,2) - utheta_d(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nTp,0,2,2,1) ./ nTp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_t  = - fkt .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_t(:,1) = 2 .* utheta_t(:,2) - utheta_t(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nhep,0,2,2,1) ./ nhep;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_he  = - fkhe .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho; 
utheta_he(:,1) = 2 .* utheta_he(:,2) - utheta_he(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nhep3,0,2,2,1) ./ nhep3;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_he3  = - fkhe3 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho; 
utheta_he3(:,1) = 2 .* utheta_he3(:,2) - utheta_he3(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nz1p,0,2,2,1) ./ nz1p;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_z1  = - fkz1 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zimp .* gtheta ./ dpsidrho;
utheta_z1(:,1) = 2 .* utheta_z1(:,2) - utheta_z1(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nz2p,0,2,2,1) ./ nz2p;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_z2  = - fkz2 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zmax .* gtheta ./ dpsidrho;
utheta_z2(:,1) = 2 .* utheta_z2(:,2) - utheta_z2(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nbp,0,2,2,1) ./ max(1,nbp);
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_b  = - fkb .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho;
utheta_b(:,1) = 2 .* utheta_b(:,2) - utheta_b(:,3);


if option.Sn_fraction > 0
    gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
    gtheta2 = pdederive(x,phys.e .* profli.tip .* nwp,0,2,2,1) ./ nwp;
    stheta  = min(sign(gtheta1),sign(gtheta2));
    gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
    utheta_w  = - fkw .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ max(1,z0wavez(profli.tep)) .* gtheta ./ dpsidrho;
    utheta_w(:,1) = 2 .* utheta_w(:,2) - utheta_w(:,3);
    utheta_sn  = - fksn .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ max(1,z0snavez(profli.tep)) .* gtheta ./ dpsidrho;
    utheta_sn(:,1) = 2 .* utheta_sn(:,2) - utheta_sn(:,3);
else
    gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
    gtheta2 = pdederive(x,phys.e .* profli.tip .* nwp,0,2,2,1) ./ nwp;
    stheta  = min(sign(gtheta1),sign(gtheta2));
    gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
    utheta_w  = - fkw .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ max(1,z0wavez(profli.tep)) .* gtheta ./ dpsidrho;
    utheta_w(:,1) = 2 .* utheta_w(:,2) - utheta_w(:,3);
end


% les champ en Rmax
a = interp1_imas(cons.temps,geo.a,profli.temps,'pchip','extrap');
rmax         = profli.Raxe + a * x;
btor         = option.signe .* (profli.fdia ./rmax);
grho         = abs((profli.rmx(:,end) * ve) ./ max(eps,pdederive(x,rmax,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,profli.psi,0,2,2,1)./ rmax .* grho ./ (profli.rmx(:,end) * ve);
btot         = sqrt(btor .^ 2 + bpol .^ 2);

switch option.mino
    case 'He3'
        ag = 3;
        zg = 2;
    case 'T'
        ag = 3;
        zg = 1;
    case 'He4'
        ag = 4;
        zg = 2;
    case 'D'
        ag = 2;
        zg = 1;
    case 'B'
        ag = 11;
        zg = 5;
    otherwise
        ag = 1;
        zg = 1;
end

switch option.gaz
    case 1
        zj = 1;
        aj = 1;
    case 2
        zj = 1;
        aj = 2;
    case 3
        zj = 1;
        aj = mean((2  + 3 .*  cons.iso) ./ (1 + cons.iso));
    case 4
        zj = 2;
        aj = 4;
    case 5
        zj = mean((1  + 4 .* cons.iso)   ./  (1 + 2 .* cons.iso));
        aj = mean((2  + 3 .* cons.iso)   ./  (1 + cons.iso));
    case 11
        zj = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
        aj = mean((1  + 11 .* cons.iso)  ./  (1 + cons.iso));
        
end


% calcul de la rotation poloidal pour l'impurete principale
% Y. B. Kim et all Phys. Fluids. B 3  (8) 1991  p 2050-
switch option.gaz
case 4
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.nhep) .* 4);
	nii   = max(1e13,profli.nhep);
otherwise
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.n1p) .* 1); 
	nii   = max(1e13,profli.n1p);
end
% ontraite les impureptes a l'etat de trace ...
alpha  = min(zimp,alpha);
beta   = (27/4) .^ 2 .* (aj ./ aimp) .^ 2 ./ (15/2 + sqrt(2 .* alpha) .* sqrt(aimp ./ aj));
g      = profli.ftrap ./ max(0.01,1 - profli.ftrap);
mui_00 = g .* (          alpha + sqrt(2)             -           log(1 + sqrt(2)));
mui_10 = g .* ((3/2)  .* alpha + 4  ./ sqrt(2)       - (5/2)  .* log(1 + sqrt(2)));
mui_01 = mui_10;
mui_11 = g .* ((13/4) .* alpha + 39 ./ (4 * sqrt(2)) - (25/4) .* log(1 + sqrt(2)));
D      = mui_00 .* (mui_11 + sqrt(2) + alpha - alpha .* beta) - mui_01 .^ 2;
D(D==0) = 1e38;
K1     = mui_01 ./ D .* (sqrt(2) + alpha - alpha .* beta);
K2     = (mui_00 .* mui_11 - mui_01 .* mui_10) ./ D;
vth    = sqrt(2 .* profli.tip .* phys.e ./ (phys.mp .* aj));
b2     = sqrt(profli.bpol .^ 2 + (profli.fdia .* profli.ri) .^ 2);
rhoi   = 4.57e-3 .* sqrt(aj .* profli.tip ./ 1e3 ./ b2);
ltim1  = pdederive(x,profli.tip,0,2,2,1) ./ profli.tip ./ (profli.rmx(:,end) * ve); 
lpiim1 = pdederive(x,profli.tip .* nii,0,2,2,1)  ./ (profli.tip .* nii)  ./ (profli.rmx(:,end) * ve); 
lpiIm1 = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.tip .* nz1p) ./ (profli.rmx(:,end) * ve); 
%figure(21);plot(x,1./ltim1,'b',x,1./lpiim1,'r',x,1./lpiIm1,'g');drawnow
%
vtehta_z1 = option.signe .*  0.5 .* vth .* rhoi .* ((K1 + (3/2) .* K2) .* ltim1 - lpiim1 + (zj ./ zimp) .* 1 .* lpiIm1) .* (profli.fdia .* profli.ri) ./ sqrt(b2);
utheta_z1 = vtehta_z1 ./ max(eps,profli.bpol);
utheta_z1(:,1) = 0;

% calcul de la rotation poloidal pour l'impurete principale
% Y. B. Kim et all Phys. Fluids. B 3  (8) 1991  p 2050-
switch option.gaz
case 4
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.nhep) .* 4);
	nii   = max(1e13,profli.nhep);
otherwise
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.n1p) .* 1); 
	nii   = max(1e13,profli.n1p);
end
% ontraite les imurepte a l'etat de trace ...
alpha  = min(zmax,alpha);
beta   = (27/4) .^ 2 .* (aj ./ aimp) .^ 2 ./ (15/2 + sqrt(2 .* alpha) .* sqrt(amax ./ aj));
g      = profli.ftrap ./ max(0.01,1 - profli.ftrap);
mui_00 = g .* (          alpha + sqrt(2)             -           log(1 + sqrt(2)));
mui_10 = g .* ((3/2)  .* alpha + 4  ./ sqrt(2)       - (5/2)  .* log(1 + sqrt(2)));
mui_01 = mui_10;
mui_11 = g .* ((13/4) .* alpha + 39 ./ (4 * sqrt(2)) - (25/4) .* log(1 + sqrt(2)));
D      = mui_00 .* (mui_11 + sqrt(2) + alpha - alpha .* beta) - mui_01 .^ 2;
D(D==0) = 1e38;
K1     = mui_01 ./ D .* (sqrt(2) + alpha - alpha .* beta);
K2     = (mui_00 .* mui_11 - mui_01 .* mui_10) ./ D;
vth    = sqrt(2 .* profli.tip .* phys.e ./ (phys.mp .* aj));
b2     = sqrt(profli.bpol .^ 2 + (profli.fdia .* profli.ri) .^ 2);
rhoi   = 4.57e-3 .* sqrt(aj .* profli.tip ./ 1e3 ./ b2);
ltim1  = pdederive(x,profli.tip,0,2,2,1) ./ profli.tip ./ (profli.rmx(:,end) * ve); 
lpiim1 = pdederive(x,profli.tip .* nii,0,2,2,1)  ./ (profli.tip .* nii)  ./ (profli.rmx(:,end) * ve); 
lpiIm1 = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.tip .* nz1p) ./ (profli.rmx(:,end) * ve); 
%figure(21);plot(x,1./ltim1,'b',x,1./lpiim1,'r',x,1./lpiIm1,'g');drawnow
%
vtehta_z2 = option.signe .*  0.5 .* vth .* rhoi .* ((K1 + (3/2) .* K2) .* ltim1 - lpiim1 + (zj ./ zimp) .* 1 .* lpiIm1) .* (profli.fdia .* profli.ri) ./ sqrt(b2);
utheta_z2 = vtehta_z2 ./ max(eps,profli.bpol);
utheta_z2(:,1) = 0;



% changement de repere
dpsidrho =  pdederive(x,profli.psi,0,2,2,1) ./ (profli.rmx(:,end) * ve);
dpsidrho(:,1) = 0;
% calcul du champ electrique radial (Er gradient(rho))
if option.Sn_fraction > 0
    Ptor    = phys.mp .* pdederive(x,profli.tip .* nHp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nDp,0,2,2,1) + ...
        3 .* phys.mp .* pdederive(x,profli.tip .* nTp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nhep,0,2,2,1) + ...
        (3/2) .* phys.mp .* pdederive(x,profli.tip .* nhep3,0,2,2,1) + ...
        aimp ./ zimp .* phys.mp .* pdederive(x,profli.tip .* nz1p,0,2,2,1) + ...
        amax ./ zmax .* phys.mp .* pdederive(x,profli.tip .* nz2p,0,2,2,1) + ...
        11/5 .* phys.mp .* pdederive(x,profli.tip .* nbp,0,2,2,1) + ...
        (1 - option.Sn_fraction) .* 183.84  ./ max(1,z0wavez(profli.tep)) .* phys.mp .* pdederive(x,profli.tip .* nwp,0,2,2,1) + ...
        option.Sn_fraction .* 118.71  ./ max(1,z0snavez(profli.tep)) .* phys.mp .* pdederive(x,profli.tip .* nwp,0,2,2,1);
    
else
    Ptor    = phys.mp .* pdederive(x,profli.tip .* nHp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nDp,0,2,2,1) + ...
        3 .* phys.mp .* pdederive(x,profli.tip .* nTp,0,2,2,1) + ...
        2 .* phys.mp .* pdederive(x,profli.tip .* nhep,0,2,2,1) + ...
        (3/2) .* phys.mp .* pdederive(x,profli.tip .* nhep3,0,2,2,1) + ...
        aimp ./ zimp .* phys.mp .* pdederive(x,profli.tip .* nz1p,0,2,2,1) + ...
        amax ./ zmax .* phys.mp .* pdederive(x,profli.tip .* nz2p,0,2,2,1) + ...
        11/5 .* phys.mp .* pdederive(x,profli.tip .* nbp,0,2,2,1) + ...
        183.84  ./ max(1,z0wavez(profli.tep)) .* phys.mp .* pdederive(x,profli.tip .* nwp,0,2,2,1);
end


Ptor    = Ptor ./ (profli.rmx(:,end) * ve);



% sorties
vtor         = ones(length(profli.temps),length(profli.xli),8);
vpol         = ones(length(profli.temps),length(profli.xli),8);

switch option.mode_vtheta
case 'same v_tor'

      vtor_H   = profli.vtor;
      vtor_D   = profli.vtor;
      vtor_T   = profli.vtor;
      vtor_he3 = profli.vtor;
      vtor_he  = profli.vtor;
      vtor_z1  = profli.vtor;
      vtor_z2  = profli.vtor;
      vtor_w   = profli.vtor;
      vtor_b   = profli.vtor;
      if option.Sn_fraction > 0
         vtor_sn   = profli.vtor;
      end

      % calculde la roration poloidale pour imp
      inter          = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zimp .* nz1p);
      warning off
      omega_z1       = (profli.er - inter) ./ dpsidrho; 
      utheta_z1      = (vtor_z1    - omega_z1 .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_z1(:,1) = 0;
      warning on
      vtheta_z1   = utheta_z1 .* bpol;

      % calculde la roration poloidale pour max
      inter          = pdederive(x,profli.tip .* nz2p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zmax .* nz2p);
      warning off
      omega_z2       = (profli.er - inter) ./ dpsidrho; 
      utheta_z2      = (vtor_z2    - omega_z2 .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_z2(:,1) = 0;
      warning on
      vtheta_z2   = utheta_z2 .* bpol;

      % calculde la roration poloidale pour W & Sn
      inter          = pdederive(x,profli.tip .* nwp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (max(1,z0wavez(profli.tep)) .* nwp);
      warning off
      omega_w        = (profli.er - inter) ./ dpsidrho;
      utheta_w       = (vtor_w    - omega_w .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_w(:,1)  = 0;
      warning on
      vtheta_w       = utheta_w .* bpol;
      if option.Sn_fraction > 0
          inter          = pdederive(x,profli.tip .* nwp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (max(1,z0snavez(profli.tep)) .* nwp);
          warning off
          omega_sn        = (profli.er - inter) ./ dpsidrho;
          utheta_sn       = (vtor_sn    - omega_sn .* rmax) ./ profli.fdia .* rmax .* option.signe;
          utheta_sn(:,1)  = 0;
          warning on
          vtheta_sn       = utheta_sn .* bpol;
      end
      
      % calculde la roration poloidale pour he
      inter          = pdederive(x,profli.tip .* nhep,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nhep);
      warning off
      omega_he       = (profli.er - inter) ./ dpsidrho; 
      utheta_he      = (vtor_he    - omega_he .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_he(:,1) = 0;
      warning on
      vtheta_he      = utheta_he .* bpol;

      % calculde la roration poloidale pour he3
      inter           = pdederive(x,profli.tip .* nhep3,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nhep3);
      warning off
      omega_he3       = (profli.er - inter) ./ dpsidrho; 
      utheta_he3      = (vtor_he3    - omega_he3 .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_he3(:,1) = 0;
      warning on
      vtheta_he3      = utheta_he3 .* bpol;
      
       % calculde la roration poloidale pour B11
      inter           = pdederive(x,profli.tip .* nbp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nbp);
      warning off
      omega_b       = (profli.er - inter) ./ dpsidrho; 
      utheta_b      = (vtor_b    - omega_b .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_b(:,1) = 0;
      warning on
      vtheta_b      = utheta_b .* bpol;
      
      % calculde la roration poloidale pour H
      inter           = pdederive(x,profli.tip .* nHp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nHp);
      warning off
      omega_H         = (profli.er - inter) ./ dpsidrho; 
      utheta_h        = (vtor_H    - omega_H .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_h(:,1)   = 0;
      warning on
      vtheta_h        = utheta_h .* bpol;
      
      % calculde la roration poloidale pour D
      inter           = pdederive(x,profli.tip .* nDp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nDp);
      warning off
      omega_D         = (profli.er - inter) ./ dpsidrho; 
      utheta_d        = (vtor_D    - omega_D .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_d(:,1)   = 0;
      warning on
      vtheta_d        = utheta_d .* bpol;
      
      % calculde la roration poloidale pour D
      inter           = pdederive(x,profli.tip .* nTp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nTp);
      warning off
      omega_T         = (profli.er - inter) ./ dpsidrho; 
      utheta_t        = (vtor_T    - omega_T .* rmax) ./ profli.fdia .* rmax .* option.signe;
      utheta_t(:,1)   = 0;
      warning on
      vtheta_t        = utheta_t .* bpol;
      
otherwise

      % calul de la vitessse toroidal en Rmax pour l'impurete principale
      vtheta_z1   = utheta_z1 .* bpol;
      vtheta_z2   = utheta_z2 .* bpol;
      vtheta_he   = utheta_he  .* bpol;
      vtheta_he3  = utheta_he3 .* bpol;
      vtheta_h    = utheta_h .* bpol;
      vtheta_d    = utheta_d .* bpol;
      vtheta_t    = utheta_t .* bpol;
      vtheta_b    = utheta_b .* bpol;
      vtheta_w    = utheta_w .* bpol;
      if option.Sn_fraction > 0
          vtheta_sn    = utheta_sn .* bpol;      
      end
      
      % calculde la roration toroidale
      warning off
      inter          = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zimp .* nz1p);
      omega_z1       = (profli.er - inter) ./ dpsidrho; 
      vtor_z1        = omega_z1 .* rmax + utheta_z1 .* option.signe .* profli.fdia ./ rmax;  
      vtor_z1(:,1)   = 2 .* vtor_z1(:,2) - vtor_z1(:,3);

      inter          = pdederive(x,profli.tip .* nz2p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zmax .* nz2p);
      omega_z2       = (profli.er - inter) ./ dpsidrho; 
      vtor_z2        = omega_z2 .* rmax + utheta_z2 .* option.signe .* profli.fdia ./ rmax;  
      vtor_z2(:,1)   = 2 .* vtor_z2(:,2) - vtor_z2(:,3);

      inter         = pdederive(x,profli.tip .* nwp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (max(1,z0wavez(profli.tep)) .* nwp);
      omega_w       = (profli.er - inter) ./ dpsidrho; 
      vtor_w        = omega_w .* rmax + utheta_w .* option.signe .* profli.fdia ./ rmax;  
      vtor_w(:,1)   = 2 .* vtor_w(:,2) - vtor_w(:,3);
      if option.Sn_fraction > 0
          inter         = pdederive(x,profli.tip .* nwp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (max(1,z0snavez(profli.tep)) .* nwp);
          omega_sn       = (profli.er - inter) ./ dpsidrho;
          vtor_sn        = omega_sn .* rmax + utheta_sn .* option.signe .* profli.fdia ./ rmax;
          vtor_sn(:,1)   = 2 .* vtor_sn(:,2) - vtor_sn(:,3);
      end

      inter          = pdederive(x,profli.tip .* nhep,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nhep);
      omega_he       = (profli.er - inter) ./ dpsidrho; 
      vtor_he        = omega_he .* rmax + utheta_he .* option.signe .* profli.fdia ./ rmax;  
      vtor_he(:,1)   = 2 .* vtor_he(:,2) - vtor_he(:,3);

      inter          = pdederive(x,profli.tip .* nhep3,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nhep3);
      omega_he3      = (profli.er - inter) ./ dpsidrho; 
      vtor_he3       = omega_he3 .* rmax + utheta_he3 .* option.signe .* profli.fdia ./ rmax;  
      vtor_he3(:,1)  = 2 .* vtor_he3(:,2) - vtor_he3(:,3);

      inter          = pdederive(x,profli.tip .* nhep3,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nbp);
      omega_b      = (profli.er - inter) ./ dpsidrho; 
      vtor_b       = omega_b .* rmax + utheta_b .* option.signe .* profli.fdia ./ rmax;  
      vtor_b(:,1)  = 2 .* vtor_b(:,2) - vtor_b(:,3);

      inter         = pdederive(x,profli.tip .* nHp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nHp);
      omega_H       = (profli.er - inter) ./ dpsidrho; 
      vtor_H        = omega_H .* rmax + utheta_h .* option.signe .* profli.fdia ./ rmax;  
      vtor_H(:,1)   = 2 .* vtor_H(:,2) - vtor_H(:,3);

      inter         = pdederive(x,profli.tip .* nDp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nDp);
      omega_D       = (profli.er - inter) ./ dpsidrho; 
      vtor_D        = omega_D .* rmax + utheta_d .* option.signe .* profli.fdia ./ rmax;  
      vtor_D(:,1) = 2 .* vtor_D(:,2) - vtor_D(:,3);

      inter         = pdederive(x,profli.tip .* nTp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nTp);
      omega_T       = (profli.er - inter) ./ dpsidrho; 
      vtor_T        = omega_T .* rmax + utheta_t .* option.signe .* profli.fdia ./ rmax;  
      vtor_T(:,1)   = 2 .* vtor_T(:,2) - vtor_T(:,3);
						      
      warning on
end		 
%
% mise en forme
%
vtor(:,:,1) = vtor_H;
vtor(:,:,2) = vtor_D;
vtor(:,:,3) = vtor_T;
vtor(:,:,4) = vtor_he3;
vtor(:,:,5) = vtor_he;
vtor(:,:,6) = vtor_z1;
vtor(:,:,7) = vtor_z2;
vtor(:,:,8) = vtor_w;
if option.Sn_fraction > 0
    vtor(:,:,9)  = vtor_sn;
else
    vtor(:,:,9) = zeros(size(vtor_w));   
end
if all(isfinite(vtor_b(:)))
    vtor(:,:,10) = vtor_b;
else
    vtor(:,:,10) = 0;    
end

vpol(:,:,1) = vtheta_h;
vpol(:,:,2) = vtheta_d;
vpol(:,:,3) = vtheta_t;
vpol(:,:,4) = vtheta_he3;
vpol(:,:,5) = vtheta_he;
vpol(:,:,6) = vtheta_z1;
vpol(:,:,7) = vtheta_z2;
vpol(:,:,8) = vtheta_w;
if option.Sn_fraction > 0
    vpol(:,:,9) = vtheta_sn;
else
    vpol(:,:,9) = zeros(size(vtheta_w));
end
if all(isfinite(vtheta_b(:)))
    vpol(:,:,10) = vtheta_b;
else
    vpol(:,:,10) = 0;
end

omega(:,:,1) = omega_H;
omega(:,:,2) = omega_D;
omega(:,:,3) = omega_T;
omega(:,:,4) = omega_he3;
omega(:,:,5) = omega_he;
omega(:,:,6) = omega_z1;
omega(:,:,7) = omega_z2;
omega(:,:,8) = omega_w;
if option.Sn_fraction > 0
    omega(:,:,9) = omega_sn;  
else
    omega(:,:,9) = zeros(size(omega_w));      
end
if all(isfinite(omega_b(:)))
    omega(:,:,10) = omega_b;
else
    omega(:,:,10) = 0;
end

iso   = interp1_imas(cons.temps,cons.iso,profli.temps,'pchip','extrap');
ftnbi = (real(cons.pnbi) .* real(cons.ftnbi) + imag(cons.pnbi) .* imag(cons.ftnbi)) ./ ...
        max(1,real(cons.pnbi) + imag(cons.ftnbi));
ftnbi = interp1_imas(cons.temps,ftnbi,profli.temps,'pchip','extrap');
mtor_H = 0 .* omega_H;
mtor_D = 0 .* omega_D;
mtor_T = 0 .* omega_T;
mtor_he = 0 .* omega_he;
mtor_he3 = 0 .* omega_he3;
mtor_b= 0 .* omega_b;
%
switch option.gaz
    case 1
        mtor_H = mtor_H + profli.rot_n0;
    case 2
        mtor_D = mtor_D + profli.rot_n0;
    case 3
        mtor_D = mtor_D + profli.rot_n0 ./ (1 + iso * ones(1,size(mtor_D,2)));
        mtor_T = mtor_T + profli.rot_n0 .* (iso * ones(1,size(mtor_T,2))) ./ (1 + iso * ones(1,size(mtor_D,2)));
    case 5
        mtor_D = mtor_D + profli.rot_n0 ./ (1 + iso * ones(1,size(mtor_D,2)));
        mtor_he3 = mtor_he3 + profli.rot_n0 .* (iso * ones(1,size(mtor_he3,2))) ./ (1 + iso * ones(1,size(mtor_D,2)));
    case 11
        mtor_H = mtor_H + profli.rot_n0 ./ (1 + iso * ones(1,size(mtor_D,2)));
        mtor_b = mtor_b + profli.rot_n0 .* (iso * ones(1,size(mtor_b,2))) ./ (1 + iso * ones(1,size(mtor_D,2)));
    case 4
        mtor_he = mtor_he + profli.rot_n0;
end

switch gas_nbi
    case -1
        %mtor_D = mtor_D;
        mtor_H = mtor_H + profli.rot_nbi;
    case 3
        mtor_D = mtor_D + profli.rot_nbi .* (1 - ftnbi * ones(1,size(mtor_D,2)));
        mtor_T = mtor_T + profli.rot_nbi .* (ftnbi * ones(1,size(mtor_T,2)));
    case 5
        mtor_D   = mtor_D + profli.rot_nbi .* (1 - ftnbi * ones(1,size(mtor_D,2)));
        mtor_he3 = mtor_he3 + profli.rot_nbi .* (ftnbi * ones(1,size(mtor_he3,2))); 
    case 11
        mtor_H   = mtor_H + profli.rot_nbi .* (1 - ftnbi * ones(1,size(mtor_H,2)));
        mtor_b = mtor_b + profli.rot_nbi .* (ftnbi * ones(1,size(mtor_b,2))); 
    otherwise
        mtor_D = mtor_D + profli.rot_nbi .* (1 - ftnbi * ones(1,size(mtor_D,2)));
        mtor_H = mtor_H + profli.rot_nbi .* (ftnbi * ones(1,size(mtor_H,2)));
end

mtor(:,:,1) = mtor_H;
mtor(:,:,2) = mtor_D;
mtor(:,:,3) = mtor_T;
mtor(:,:,4) = mtor_he3;
mtor(:,:,5) = mtor_he;
mtor(:,:,6) = 0;
mtor(:,:,7) = 0;
mtor(:,:,8) = 0;
mtor(:,:,9) = 0;
if all(isfinite(mtor_b(:)))
    mtor(:,:,10) = mtor_b;
else
    mtor(:,:,10) = 0;
end

% effect of toroidal field orientation
if isfield(option,'orientation')

  rtor   = option.orientation .* rtor;
  vtor   = option.orientation .* vtor;
  vpol   = - option.orientation .* vpol;
  omega  = option.orientation .* omega;
  mtor   = option.orientation .* mtor;

end