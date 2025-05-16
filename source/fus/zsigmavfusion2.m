% ZSIGMAVFUSION calcul les reactivites thermique des reactions de fusion
%-------------------------------------------------------------------------
% fichier zsigmavfusion.m ->  zsigmavfusion
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule les reactivites thermique des reactions de fusion  
% et donne les energies des produits de fusion. 
%  
% syntaxe  :
%  
%     [dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(ti);
%
% entrees :
%
%     ti = tabelau des temperatures ioniques (eV) 
%
% sorties :
% 
%     dd_p      = structure de donnees pour la reaction  D + D   -> T + p
%     dd_n      = structure de donnees pour la reaction  D + D   -> He3 + n
%     dt        = structure de donnees pour la reaction  D + T   -> He4 + n
%     dhe3      = structure de donnees pour la reaction  D + He3 -> He4 + p
%     tt        = structure de donnees pour la reaction  T + T   -> He4 + 2 n
%     the3_pn   = structure de donnees pour la reaction  T + He3 -> He4 + p + n
%     the3_d    = structure de donnees pour la reaction  T+ He3  -> He4 + D
% 
% les structures ont le format suivant :
%  
%     x.sv      = reactivites thermique <sigma*v>|maxwellienne (m^3 * s^-1)
%     x.delta   = precision relative d(sv)/sv
%     x.timin   = temperature minimale pour la validite des donnees (eV)
%     x.timax   = temperature maximale pour la validite des donnees (eV)
%     x.p       = energie des protons produits (eV, 0 si aucun)
%     x.n       = energie des neutrons produits (eV, somme de tous les termes, 0 si aucun)
%     x.d       = energie des noyaux de deuterium produits (eV, 0 si aucun)
%     x.t       = energie des noyaux de tritium produits  (eV, 0 si aucun)
%     x.he3     = energie des noyaux d'helium 3 produits (eV, 0 si aucun)
%     x.he4     = energie des noyaux d'helium 4 produits (eV, 0 si aucun)
% 
% references :
% 
% [1] la fusion thermonucleaire controlee par confinement magnetique, 
%     Collection CEA, Masson, 1987. 
%
% [2] Etude des produits de fusion charges dans un tokamak,
%     G. Martin, these de l'Universite de Paris-Sud, 1985.
%     
% [3] Improved formulas for fusion cross-sections and thermal reactivities,
%     Nuclear Fusion, Vol. 32, N� 4, 1992
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 17/02/2000.
% 
% 
% liste des modifications : 
% 
%     [dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion2(ti);
%
% This function was added by Dehong Chen for Te > 100keV with the formula in
% Ref[4] for DDp, DDn, DT and D3He fusions on 01/12/2024:
%
% [4] Aleksey A Mavrin 2018 New analytic representation of the thermonuclear
%     reaction rates. Plasma Phys. Control. Fusion 60 092001
%--------------------------------------------------------------
%
function [dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion2(ti)

% ti doit etre en keV
tikev    = max(1e-3,real(ti) ./ 1e3);

% pour D + D -> T + p  : D(d,p)T
dd_p.sv    = zformsv_NR(tikev, 1.46e-14, 18.81, -0.105, 0.0819, -2.34e-3, 0.0, 1/eps);
dd_p.delta = 2e-2;
dd_p.timin = 1e3;
dd_p.timax = 1000e3;
dd_p.p     = 3.02e6; % (ev)
dd_p.n     = 0;
dd_p.d     = 0;
dd_p.t     = 1.01e6; % (eV)
dd_p.he3   = 0;
dd_p.he4   = 0;
dd_p.sv    = zbornesv(dd_p.timin,dd_p.timax,ti,dd_p.sv);

% pour D + D -> He3 + n : D(d,n)He3 
dd_n.sv    = zformsv_NR(tikev, 1.39e-14, 18.81, -0.102, 0.111, -4.06e-3, 0.0, 1/eps);
dd_n.delta = 2e-2;
dd_n.timin = 1e3;
dd_n.timax = 1000e3;
dd_n.p     = 0;
dd_n.n     = 2.45e6; % (eV)
dd_n.d     = 0;
dd_n.t     = 0;
dd_n.he3   = 0.82e6; % (eV)
dd_n.he4   = 0;
dd_n.sv    = zbornesv(dd_n.timin,dd_n.timax,ti,dd_n.sv);

% pour D + T -> He4 + n : T(d,n)He4
dt.sv      =  zformsv_NR(tikev, 1.46e-11, 19.98, -2.278, 2.194, -7.98e-1, 9.34e-2, 20) + ...
              zformsv_R(tikev, 7.35e-13, 65, 1.09e-2, -5.42e-7);
dt.delta = 6e-2;
dt.timin = 1e3;
dt.timax = 1000e3;
dt.p     = 0;
dt.n     = 14.03e6; % (eV)
dt.d     = 0;
dt.t     = 0;
dt.he3   = 0;
dt.he4   = 3.56e6; % (eV)
dt.sv    = zbornesv(dt.timin,dt.timax,ti,dt.sv);

% pour D + He3 -> He4 + p : He3(d,p)He4
dhe3.sv    =  zformsv_NR(tikev, 2.52e-12, 31.72, -0.389, 0.117, 3.35e-2, 0.0, 41) + ...
              zformsv_R(tikev, 1.12e-12, 240, 6.07e-3, -5.51e-7);
dhe3.delta = 4e-2;
dhe3.timin = 1e3;
dhe3.timax = 1000e3;
dhe3.p     = 14.64e6; % (eV)
dhe3.n     = 0;
dhe3.d     = 0;
dhe3.t     = 0;
dhe3.he3   = 0;
dhe3.he4   = 3.71e6; % (eV) 
dhe3.sv    = zbornesv(dhe3.timin,dhe3.timax,ti,dhe3.sv);

% pour T + T -> He4 + 2n : T(t,2n)He4
% donnees pour le spline
t  = [0.2  ,   1.22,  3.3,   17.4,     60,       600,    1000]; % (keV)
sv = [3e-33,   1e-27, 4e-26, 1.79e-24, 1.2e-23, 9e-23,  8e-23]; % (m^3 * s^-1)
tt.sv    = 10 .^ interp1(log10(t),log10(sv),log10(tikev),'spline');
tt.delta = 0.25;
tt.timin = 0.2e3;
tt.timax = 1000e3;
tt.p     = 0; 
tt.n     = 11.332e6; % (eV -> 2 neutrons);
tt.t     = 0;
tt.d     = 0;
tt.he3   = 0;
tt.he4   = 0e6; % presque rien
tt.sv    = zbornesv(tt.timin,tt.timax,ti,tt.sv);

% new cross section from D; Coster
[sigma_,sv] = tt_cross_section(tikev .* 1e3);
%hf =gcf;figure(21);clf;loglog(tikev,sv,'.r',tikev,tt.sv);drawnow;figure(hf);
tt.sv    = sv;
tt.delta = 0.25;
tt.timin = 0.2e3;
tt.timax = 10000e3;
tt.p     = 0; 
tt.n     = 11.332e6; % (eV -> 2 neutrons);
tt.t     = 0;
tt.d     = 0;
tt.he3   = 0;
tt.he4   = 0e6; % presque rien
tt.sv    = zbornesv(tt.timin,tt.timax,ti,tt.sv);

% total T + He3 -> X + Y (+ Z)
t  = [1  ,   6.5,   10,    17.4,     60,       1000]; % (keV)
sv = [1e-33, 1e-27, 1e-26, 1.25e-25, 1.48e-23, 3e-22]; % (m^3 * s^-1)
the3.sv    = 10 .^ interp1(log10(t),log10(sv),log10(tikev),'spline');
the3.delta = 0.25;
the3.timin = 1e3;
the3.timax = 1000e3;
the3.p     = 0;
the3.n     = 0;
the3.t     = 0;
the3.d     = 0;
the3.he3   = 0;
the3.he4   = 0;
the3.sv    = zbornesv(the3.timin,the3.timax,ti,the3.sv);
% pour T + He3 -> He4 + p + n : He3(t,n+p)He4
the3_pn =the3;
the3_pn.sv   = 0.57 * the3.sv;
the3_pn.p    = (12.96e6-12.1e6) .* (4/5) ; % approximation grossiere
the3_pn.n    = 12.1e6; % (eV)
the3_pn.he4  = (12.96e6-12.1e6) .* (1/5) ; % approximation grossiere
% pour T + He3 -> He4 + D : He3(t,d)He4
the3_d      = the3;
the3_d.sv   = 0.43 * the3.sv;
the3_d.he4  = 4.8e6; % (eV)
the3_d.d    = 9.52e6; % (eV)

% formule de fit des section efficace
function sv_NR =zformsv_NR(Ti, A, EG, C1, C2, C3, C4, Tmix)

sv_NR = A ./ Ti .^ (2/3) .* exp(- EG ./ Ti .^ (1/3) - (Ti ./ Tmix) .^ 2) .* ...
        (1 + C1 .* Ti .^ (1/3) + C2 .* Ti .^ (2/3) + C3 .* Ti .^ (3/3) + C4 .*Ti .^ (4/3) );
    
% convertion cm^3 * s^-1 -> m^3 * s^-1
sv_NR   = sv_NR ./ 1e6;

function sv_R = zformsv_R(Ti, B, E0, C1, C2)

sv_R = B ./ Ti .^ (3/2) .* exp(- E0 ./ Ti) .* ...
       (1 + C1 .*Ti + C2 .* Ti.^ 2) ;

% convertion cm^3 * s^-1 -> m^3 * s^-1
sv_R   = sv_R ./ 1e6;


% fonction pour les bornes  sur ti (sv=0 hors intervalle)
function sv=zbornesv(timin,timax,ti,sv)

ind = find(ti < min(timin,timax) | ti > max(timin,timax));
sv(ind)=zeros(1,length(ind));

