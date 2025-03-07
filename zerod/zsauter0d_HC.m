% fonction bootstrap Sauter corrigee (j*b)
% avec modification pour le piedestal 
% references :
%  O. sauter et al, PoP vol 6 (1999) p 2834-2839	
%  O. Sauter et al, PoP vol 9 (2002) p 5140-5141 (errat)
%  S. Koh et al, PoP vol 19 (2012)   p 072505-1 to 072505-12 
%  R. Hager and C;S; Chang PoP vol 23 (2016) p 042503- 
%  C.S Chang, Physics of Plasmas vol 4 (1197) # 6 p 2241-
%
%  jboot must be devided by b0 to have the current density

function jboot = zsauter0d_HC(x,te,ti,ne,ni,q,Zeff,Zion,Raxe,ft,epsi,dpsidx,F,modeh,fgg,model_zi,rmx,r2i,grho2r2,force_sauter)

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


% le facteur est un fit avec NClass (l'integrale n'est pas precise au bord)
if nargin < 15
    fgg = 1;
end
if nargin < 16
    model_zi = 1;
end
if nargin < 20
  force_sauter = 0;
end
ve   = ones(size(x));
vt   = ones(size(te,1),1);

% various effective charges
Z_bar   = ne ./ ni;
Z_alpha = Zion;
Z       = (Zion .^ 2 .* Z_bar .* Zeff) .^ (1/4);
% rules for Z
% for electron Z = Zeff
% for ions collisionality Z = Z
% elsewhere, including Coulomb logarithm : Z = Z_bar;
% Z_alpha = Zion
if force_sauter == 1
    Z_bar = Zeff;
    Z     = Zion;
end
% collisionality
lni            = 30.0-log(Z_bar.^3 .* sqrt(ni) ./ (ti.^(1.5)));
lne            = 31.3-log(sqrt(ne)./te);
fact           = Raxe ./ (epsi.^(1.5));
nustari        = fact .* (phys.e .^ 2 ./ (12 .* pi .^ (3 ./ 2) .* phys.epsi0 .^ 2)) .* ni .* lni .* q ./ (ti .^ 2) .* (Z .^ 4);
nustare        = fact .* (sqrt(2) .* phys.e .^ 2 ./ (12 .* pi .^ (3 ./ 2) .* phys.epsi0 .^ 2)) .* ne .* lne .* q ./ (te .^ 2) .* Zeff;

% Sauter part
alpha_0 = - 1.17 .* ( 1 - ft) ./ (1 - 0.22 .* ft  - 0.19 .* ft .^ 2);
alpha_i = ((alpha_0 + 0.25 .* (1 - ft .^ 2) .* sqrt(nustari)) ./ (1 + 0.5 .* sqrt(nustari)) + 0.315 .* nustari .^ 2 .* ft .^ 6) ./  ...
          (1 + 0.15 .* nustari .^ 2 .* ft .^ 6);

% effective trapped particles fractions
% Z = Z_bar
% f_t_eff_31 is modified from A9 to formula 25
f_t_eff_31    = ft ./ (1 + (1 - 0.1 .* ft).* sqrt(nustare) + 0.5 .* (1 - 0.99 .* ft) .* (nustare ./ Z_bar));
f_t_eff_32_ee = ft ./ (1 + 0.26 .* (1 - ft) .* sqrt(nustare) + 0.18 .* (1 - 0.37 .* ft) .* (nustare ./ sqrt(Z_bar)));
f_t_eff_32_ei = ft ./ (1 + (1 + 0.6 .* ft) .* sqrt(nustare) + 0.85 .* (1 - 0.37 .* ft) .* nustare .* (1 + Z_bar));
f_t_eff_34    = ft ./ (1 + (1 - 0.1 .* ft).* sqrt(nustare) + 0.5 .* (1 - 0.5 .* ft) .* (nustare ./ Z_bar));
%
L_31 = F_31(f_t_eff_31,Z_bar);
L_32 = F_32(f_t_eff_32_ee,f_t_eff_32_ei,Z_bar);
L_34 = F_31(f_t_eff_34,Z_bar);

% new terms
% geometical effect:
% best match with full equilibrium computation is the formula with dpsidrho
dpsidrho = - rmx  ./ q .* (F(:,end) * ve) ./ (Raxe(:,end) * ve);
c_B = 3 ./ 2 .* dpsidrho .^ 2 ./ F .^ 2 .* r2i;
%c_B_bis =  3.*  epsi .^ 2 ./ 2 ./ q .^ 2 .* grho2r2;
%c_B_helander =      3 .* epsi .^ 2 ./ 2 ./ q .^ 2 ./ Raxe .^ 2;
% from appendix B
% ion coef
ci1 = 1.357;
ci2 = 2.192;
ci3 = 6.919;
% electron coef
ce1 = (10./12) .* (289.914 + 408 .* Zeff .^ 2) ./ (178 + 425.678 .* Zeff + 192 .* Zeff .^ 2); 
ce2 = (10./12) .* (468.105 + 912 .* Zeff .^ 2) ./ (178 + 425.678 .* Zeff + 192 .* Zeff .^ 2); 
ce3 = (10./12) .* (1477.85 + 2592 .* Zeff .^ 2) ./ (178 + 425.678 .* Zeff + 192 .* Zeff .^ 2); 

% L_h
fh2     = c_B .* q .^ 2 .* (Raxe(:,end) * ve) .^ 2 .* Z_bar .^ 2 ./ 4 ./ (sqrt(2) + Z_bar) ./ epsi .^ 3 ./ nustare .^ 2; 
L_31_h2 = fh2 .* ((4 .* sqrt(2) + 13 .* Z_bar) .* ce1 + 6 .* Z_bar .* ce2);
L_32_h2 = fh2 .* ((4 .* sqrt(2) + 13 .* Z_bar) .* ce2 + 6 .* Z_bar .* ce3);
% this is not in the original paper !
L_34_h2 = L_31_h2; % from Helander 12.51

% alpha_i_h
alpha_i_h_0 = ci2 ./ ci1;


% L_1
L_31_1    = 2 .* ft .* Z_bar .* (2.4 + Z_bar) ./ (1 + Z_bar) ./ (1 - 0.99 .* ft) ./ nustare;
L_32_1    = ft .* ((0.05 + 0.62 .* Z_bar) ./ (0.18 - 0.0666 .* ft) .* sqrt(Z_bar) - (0.56 + 1.93 .* Z_bar) ./ (0.85 - 0.3145 .* ft) ./ (1 + Z_bar)) ./  ...
                  (1 + 0.44 .* Z_bar) ./ Z_bar ./ nustare;
L_34_1    = ft .* Z_bar .* (2.4 + Z_bar) ./ (1 + Z_bar) ./ (0.5 - 0.25 .* ft) ./ nustare;
alpha_i_0 = 2.1;

% critical collisionality
a_c       = 32 .* ce1 .^ 2 + 8 .* sqrt(2) .* (13 .* ce1 .^ 2 + 9 .* ce1 .* ce2 + 2 .* ce2 .^ 2) .* Z_bar + ...
            ((13 .* ce1  + 9 .* ce2) .^ 2 + 7 .* ce2 .^ 2 + 6 .* ce3  .* (6 .* ce1 + 4 .* ce2)) .* Z_bar .^ 2;
b_c       = 6 .* ce2 .* Z_bar .* (sqrt(2) + Z_bar) + ce1 .* (8 + 17 .* sqrt(2) .* Z_bar + 13 .* Z_bar .^ 2);
nustare_c = 5 .* sqrt(c_B) .* q .* (Raxe(:,end) * ve) .* Z_bar ./ epsi .^ (3 ./ 2) .* sqrt(a_c ./ b_c);
% c5 coefficient:
c5        = (L_31_h2 .* nustare .^ 2 ./ nustare_c .^ 2 ) ./ (L_31_1 .* nustare ./ nustare_c) ./ nustare_c;

% choosen coefficients


% fitted coefficients
c2  = 0.2947;
a1  = 0.0488;
a2  = 0.0488;
b1  = 0;
b2  = 0.0086;
epsi_c_1 = 0.15;
epsi_c_2 = 0.5099;
w1  = 0.1;
w2  = 0.15;
ld_1 = 19.1702;
ld_2 = 1.9056;
ld_3 = 0.0106;
ld_4 = 2.994;
ld_5 = 0.9958;

% derived coefficients
% possible typo !
c1   = 0.5 .* ((a1 - b1) .* tanh(2 .* (epsi - epsi_c_1)  ./ w1) + (a1 + b1)) - ...
       0.5 .* ((a2 - b2) .* tanh(2 .* (epsi - epsi_c_2)  ./ w2) + (a2 - b2)); 

% gamma functions
gamma_31 = (1 + c5 .* nustare) ./ (1 + c5 .* L_31_1 ./ L_31_h2 .* nustare);
gamma_32 = (1 + c5 .* nustare) ./ (1 + c5 .* L_32_1 ./ L_32_h2 .* nustare);
gamma_34 = (1 + c5 .* nustare) ./ (1 + c5 .* L_34_1 ./ L_34_h2 .* nustare);

% gamma alpha 
gamma_alpha =(1 + c5 .* nustari) ./ (1 + c5 .* alpha_i_0 ./ alpha_i_h_0 .* nustari);


% beta correction
beta_col = (1 + c1 .* c2 .* nustare .^ 2) ./ (1 + c1 .* nustare .^ 2);

% gradients 
pel      = phys.e .* te .* ne;
pion     = phys.e .* ti .* ni;
ptot     = pel + pion;
rpel     = pel ./ max(ptot,1);
dpdpsi   = pdederive(x,ptot,0,2,2,1);
dqdpsi   = pdederive(x,q,0,2,2,1);
dtedpsi  = pdederive(x,te,0,2,2,1);
dnedpsi  = pdederive(x,ne,0,2,2,1);
dtidpsi  = pdederive(x,ti,0,2,2,1);
% correction piedestal (derivee 2 points + air du rectangle / aire triangle)
indh = find(modeh);
if ~isempty(indh) && (fgg ~= 0)
	dx            = x(end) - x(end-1);
	dpdpsi(indh,end-1)  = fgg .* (ptot(indh,end) - ptot(indh,end-1)) ./ dx;
	dqdpsi(indh,end-1)  = fgg .* (q(indh,end) - q(indh,end-1)) ./ dx;
	dtedpsi(indh,end-1) = fgg .* (te(indh,end) - te(indh,end-1)) ./ dx;
	dnedpsi(indh,end-1) = fgg .* (ne(indh,end) - ne(indh,end-1)) ./ dx;
	dtidpsi(indh,end-1) = fgg .* (ti(indh,end) - ti(indh,end-1)) ./ dx;
end
%
dpdpsi   = dpdpsi  ./ dpsidx;
dqdpsi   = dqdpsi  ./ dpsidx;
dtedpsi  = dtedpsi ./ dpsidx;
dnedpsi  = dnedpsi ./ dpsidx;
dtidpsi  = dtidpsi ./ dpsidx;


% orbite correction
% estimation de delta_psi
rext   = Raxe .* (1 + epsi);
ral    = 4.55e-3.* sqrt(ti ./ 1e3) ./ (F ./ rext); % en m
%qp     = qmin * ve + (qa -qmin) * (x .^ 2);
% patato
dp1     = Raxe .* (2 .* q .* ral ./ Raxe) .^ (2/3);
% banana
dp2     = sqrt(epsi) .* ral .* q;
dp      = dp2 .* (dp2 < (Raxe .* epsi)) + dp1 .* (dp2 >= (Raxe .* epsi));
%
delta_psi = abs(dp ./max(dp,Raxe .* epsi) .* dpsidx);
% delta_R
rint =  Raxe .* (1 - epsi);
delta_R = pdederive(x,Raxe,0,2,2,1) ./ pdederive(x,rint,0,2,2,1);
% correction du to orbit width in grad_ti
beta_grad_ti = - ld_1 .* nustari ./ (1 + ld_2 .* nustari .^ 2 + ld_3 .* nustari .^ 4) .* ...
                 (1 - epsi) .^ ld_4 .* abs(delta_R) .^ ld_5 .* delta_psi .* dtidpsi ./ ti;


% F_orbit
% sign of these terms ?
% if density decrease from center to edge Olne must be > 0
Olti = dtidpsi ./ ti;
Olne = dnedpsi ./ ne;
Olq  = dqdpsi ./ q;
F_orbit = 1 - 2 .* delta_psi .* (-(3./2) .* Olti + Olne + Olq);

if force_sauter
  % modified formula
  disp('Sauter mode')
  term_1 = L_31 .* dpdpsi  ./ pel;
  term_2 = L_32 .* dtedpsi ./ te;
  term_3 = L_34  .* alpha_i  .* dtidpsi .* (1 - rpel) ./ rpel ./ ti;
  jboot  = F .* pel .* (term_1 + term_2 + term_3);
else
  % modified formula
  term_1 = gamma_31 .* L_31 .* dpdpsi  ./ pel;
  term_2 = gamma_32 .* L_32 .* dtedpsi ./ te;
  term_3 = gamma_34 .* L_34 .* gamma_alpha .* alpha_i .* F_orbit .* dtidpsi ./ (te .* Z_bar);
  jboot  = (beta_col + beta_grad_ti) .* F .* pel .* (term_1 + term_2 + term_3);
end
% les gradients s'annulent au centre
jboot(:,1)     = 0;

%save contexte_boot_hc

% L functions
function out = F_31(x,Z)

out = (1 + 1.4 ./ ( Z + 1)) .* x - 1.9 ./ (Z + 1) .* x .^ 2 + 0.3 ./ (Z + 1) .* x .^ 3 + 0.2 ./ (Z + 1) .* x .^ 4;

%
function out = F_32(x_ee,y_ei,Z)

out = F_32_ee(x_ee,Z) + F_32_ei(y_ei,Z);

%
function out = F_32_ee(x_ee,Z)

out = (0.05 + 0.62 .* Z) ./ (Z .* (1 + 0.44 .* Z)) .* (x_ee - x_ee .^ 4) +  ...
      1 ./ (1 + 0.22 .* Z) .* (x_ee .^ 2 - x_ee .^ 4 - 1.2 .* (x_ee .^ 3 - x_ee .^ 4)) + 1.2 ./ (1 + 0.5 .* Z) .* x_ee .^ 4;

%
function out = F_32_ei(x_ei,Z)

out = -(0.56 + 1.93 .* Z) ./ (Z .* (1 + 0.44 .* Z)) .* (x_ei - x_ei .^ 4) +  ...
      4.95 ./ (1 + 2.48 .* Z) .* (x_ei .^ 2 - x_ei .^ 4 - 0.55 .* (x_ei .^ 3 - x_ei .^ 4)) - 1.2 ./ (1 + 0.5 .* Z) .* x_ei .^ 4;




