% la vitesse de pompage est un consigne externe (Pa m^3 s^-1)
if isfield(post.z0dinput.cons,'pump')
	edge.pump = post.z0dinput.cons.pump;
else
	edge.pump = input('pumping speed (m^3 s^-1) ?') .* ones(size(post.z0dinput.cons.temps));
end
edge.temps = post.z0dinput.cons.temps;
% script d'estiamtion de la pression de neutre dans la chambre
phys.c           =   2.99792458e8;             % speed of light in vacuum (m/s)  (definition)
phys.h           =   6.62606876e-34;           % Planck constant (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivity of vacuum (F/m)  (definition)
phys.g           =   6.673e-11;                % gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % Boltzmann constant (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % fine structure constant (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % electron mass (at rest) (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % proton mass (at rest) (kg)
phys.ua          =   1.66053873e-27;           % Atomic mass unit (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % Avogadro number (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)
% valeur pour le bord
edge.te_sep = post.zerod.tebord;
edge.ti_sep = post.zerod.tibord;
edge.ne_sep = post.zerod.nebord;
edge.ni_sep = post.zerod.nibord;
edge.pion = max(0,post.zerod.pion + post.zerod.pei);
edge.pel  = max(0,post.zerod.pel - post.zerod.pei - post.zerod.pradsol);
edge.psol = edge.pion + edge.pel;
% flux en Pa m^3 s^-1 (equivalent electron)
edge.gnbi = post.zerod.pnbi ./ (post.z0dinput.option.einj .* phys.e) ./ (4.41e-4 .* phys.avo);
edge.gpellet  = post.zerod.n0a./ max(eps,1 - post.zerod.frac_pellet) .* post.zerod.frac_pellet ./ (4.41e-4 .* phys.avo);
edge.grecycle = post.zerod.n0a./ max(eps,1 - post.zerod.frac_pellet) ./ (4.41e-4 .* phys.avo);
edge.gtot     = edge.grecycle + edge.gpellet + edge.gnbi;
edge.gdt_sep  = pchip(post.profil0d.temps,post.profil0d.n1p(:,end),post.zerod.temps) ./ edge.ne_sep .* edge.gtot;
edge.gcore    = edge.gpellet + edge.gnbi;
edge.eta_c    = edge.gcore ./ edge.gtot; 
% cette definition n'est pas celle du papier !
edge.epsi_ei  = edge.pel ./ edge.pion; 
% cette definition n'est pas celle du papier !
%edge.epsi_ei  = post.zerod.pel ./ post.zerod.pion; 
edge.tau_he_star  = post.zerod.nhem .* post.zerod.vp ./  post.zerod.salpha;
edge.taue         = post.zerod.taue;
edge.ghe_sep      = post.zerod.nhem .* post.zerod.vp ./  (post.zerod.tauhe - post.zerod.taue) ./ (4.41e-4 .* phys.avo);
edge.nhe_sep      = pchip(post.profil0d.temps,post.profil0d.nhep(:,end),post.zerod.temps);
edge.nimp_sep     = pchip(post.profil0d.temps,post.profil0d.nzp(:,end),post.zerod.temps);
edge.qpk          = post.zerod.peakdiv ./1e6;
edge.qpk(1)       = edge.qpk(2);
edge.qpk          = edge.qpk .* post.zerod.xpoint;
% pression dans la sol
% vitesse moyenne des neutres
% 3 eV = dissotiation
% wall = 500K
edge.en0 = (post.zerod.telim +  3 + 600 .* phys.k ./ phys.e) ./ 3;
edge.vn0 = sqrt(2 .* edge.en0 .* phys.e ./ (post.zerod.meff .* phys.ua));
% densite de neutre (m^-3)
edge.n0 = edge.gtot .* (4.41e-4 .* phys.avo) ./ post.zerod.sext ./ edge.vn0;
% pression de neutre equivalente (Pa)
edge.p0_sol = edge.en0 .* phys.e .* edge.n0;
% zonne prive du limiteur
% la vitesse de pompage est suppose constante et capable de retirer 20 Pa m^3 s^-1 dans ITER
edge.p0_div  = edge.gtot ./ edge.pump;
edge.p0      = post.zerod.xpoint  .* edge.p0_div + (~post.zerod.xpoint) .* edge.p0_sol;


% calcul utilisant les scaling pour ITER 
% ref : Kukushkin et al NF 43 (2003) p 716-723;
edge.sc.ff    = 1 + 0.18 .* edge.eta_c;
edge.sc.fhe   = 30 .* post.zerod.nhem ./ post.zerod.nem; 
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = post.z0dinput.cons.picrh + post.z0dinput.cons.pecrh + post.z0dinput.cons.pnbi + post.zerod.pohm + post.zerod.plh;
Q    = rfan .* post.zerod.pfus ./ padd;
frad = (post.zerod.prad + post.zerod.pbrem + post.zerod.pcyclo) ./ post.zerod.pin;
edge.sc.fhe   = 0.21 .* (5 .* Q ./ (Q + 5)) ./ max(eps,1-frad);
edge.sc.gstar = edge.gtot ./ 124;
edge.sc.sstar = edge.pump ./ 20;
edge.sc.pstar = edge.psol ./ 100e6;
edge.sc.p0stat = edge.p0 ./ 6.2; 
edge.sc.mu    = edge.sc.gstar ./ edge.sc.sstar .* edge.sc.pstar .^ 0.87 .* edge.sc.ff .^ 2;

% valeur non sature
edge.ns.qpk      = 7.55 .* edge.sc.pstar .^ 2 .* edge.sc.p0stat .^ -0.85;
edge.ns.ne_sep   = 3.89e19 .* edge.sc.ff .^ 0.53 .* edge.sc.pstar .^ 0.25 .* edge.epsi_ei .^ 0.05 .* edge.sc.p0stat .^ 0.36;
edge.ns.gdt_sep  = 16.4 .* edge.sc.ff .^ -3 .* edge.sc.sstar .^ 0.3 .* edge.sc.pstar .^ -0.22 .* edge.sc.p0stat .^ 0.25;
edge.ns.nhe_sep  = 3.06e17 .* edge.sc.fhe .^ 1 .* edge.sc.ff .^ -1 .* edge.sc.sstar .^ -1 .* edge.sc.pstar .^ 2.44 .* edge.epsi_ei .^ -0.1  .* edge.sc.p0stat .^ -2;
edge.ns.ghe_sep  = 0.512 .* edge.sc.fhe .^ 1 .* edge.sc.ff .^ -1 .* edge.sc.sstar .^ -1 .* edge.sc.pstar .^ 2.44 .* edge.epsi_ei .^ 0  .* edge.sc.p0stat .^ - 2.21;
edge.ns.te_sep   = 162 .* edge.sc.ff .^ -0.06 .* edge.sc.sstar .^ -0.02 .* edge.sc.pstar .^ 0.47 .* edge.epsi_ei .^ 0.05  .* edge.sc.p0stat .^ - 0.17;
edge.ns.ti_sep   = 270 .* edge.sc.ff .^ -0.32 .* edge.sc.sstar .^ -0.04 .* edge.sc.pstar .^ 0.61 .* edge.epsi_ei .^ -0.116  .* edge.sc.p0stat .^ - 0.29;


% valeur sature
edge.sv.qpk      = 7.55 .* edge.sc.ff .^ -1.7 .* edge.sc.pstar .^ 1.26 ;
edge.sv.ne_sep   = 3.89e19 .* edge.sc.ff .^ 1.25 .* edge.sc.pstar .^ 0.55 .* edge.epsi_ei .^ 0.05;
% ce n'est pas calculable avec metis (flux de neutre provenant de la SOL et non pas du limiteur)
edge.sv.gdt_sep  = 16.4 .* edge.sc.ff .^ -2.5 .* edge.sc.sstar .^ 0.3;
edge.sv.nhe_sep  = 3.06e17 .* edge.sc.fhe .^ 1 .* edge.sc.ff .^ -5 .* edge.sc.sstar .^ -1 .* edge.sc.pstar .^ 0.7 .* edge.epsi_ei .^ -0.1;
edge.sv.ghe_sep  = 0.512 .* edge.sc.fhe .^ 1 .* edge.sc.ff .^ -5.42 .* edge.sc.sstar .^ -1 .* edge.sc.pstar .^ 0.52;
edge.sv.te_sep   = 162 .* edge.sc.ff .^ -0.4 .* edge.sc.sstar .^ -0.02 .* edge.sc.pstar .^ 0.32 .* edge.epsi_ei .^ 0.049;
edge.sv.ti_sep   = 270 .* edge.sc.ff .^ -0.9 .* edge.sc.sstar .^ -0.04 .* edge.sc.pstar .^ 0.36 .* edge.epsi_ei .^ -0.115;
edge.sv.gtot     = 124 .* edge.sc.ff .^ 2 .* edge.sc.sstar .* edge.sc.pstar .^ 0.87;

% sortie
m  = (edge.sc.mu >= 1.1);
cm = ~m;
edge.out.qpk     = cm .* edge.ns.qpk + m .* edge.sv.qpk;
edge.out.ne_sep  = cm .* edge.ns.ne_sep + m .* edge.sv.ne_sep;
edge.out.gdt_sep = cm .* edge.ns.gdt_sep + m .* edge.sv.gdt_sep;
edge.out.nhe_sep = cm .* edge.ns.nhe_sep + m .* edge.sv.nhe_sep;
edge.out.ghe_sep = cm .* edge.ns.ghe_sep + m .* edge.sv.ghe_sep;
edge.out.te_sep  = cm .* edge.ns.te_sep + m .* edge.sv.te_sep;
edge.out.ti_sep  = cm .* edge.ns.ti_sep + m .* edge.sv.ti_sep;
edge.out.gtot    = cm .* edge.gtot + m .* edge.sv.gtot;



% figure
h = findobj(0,'type','figure','tag','z0plotedge');
if isempty(h)
       h=figure('tag','z0plotedge');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(4,2,1)
plot(edge.temps,edge.qpk,'b',edge.temps,edge.out.qpk,'r');
ylabel('Q_{peak} (MW)')
set(gca,'Ylim',[0,50]);
title('ITER edge')

subplot(4,2,2)
plot(edge.temps,edge.ne_sep,'b',edge.temps,edge.out.ne_sep,'r');
ylabel('n_{e,sep} (m^-^3)')
title('blue -> METIS, red -> Kukushkin scaling')

subplot(4,2,3)
plot(edge.temps,edge.sc.mu,'b',edge.temps,1.1 .* ones(size(edge.temps)),'r-.');
ylabel('\mu')

subplot(4,2,4)
plot(edge.temps,edge.nhe_sep,'b',edge.temps,edge.out.nhe_sep,'r');
ylabel('n_{He,sep} (m^-^3)')

subplot(4,2,5)
plot(edge.temps,edge.ghe_sep,'b',edge.temps,edge.out.ghe_sep,'r');
ylabel('\Gamma _{He,sep} (Pa m^3 s^-^1)')

subplot(4,2,6)
plot(edge.temps,edge.te_sep,'b',edge.temps,edge.out.te_sep,'r');
ylabel('T_{e,sep} (eV)')

subplot(4,2,7)
plot(edge.temps,edge.ti_sep,'b',edge.temps,edge.out.ti_sep,'r');
ylabel('T_{i,sep} (eV)')
xlabel('time (s)')

subplot(4,2,8)
plot(edge.temps,edge.gtot,'b',edge.temps,edge.out.gtot,'-.r');
ylabel('\Gamma _{sep} (Pa m^3 s^-^1)')
xlabel('time (s)')


