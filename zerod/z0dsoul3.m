% fonction d'apppel de SOUL3 depuis METIS
% post -> metis data
% tc -> temps (s)
% plotonoff -> flag pour le graphes
function [souldata,restart,texte] = z0dsoul3(post,tc,x_x_fraction,angle_target,plotonoff,restart,dtsoul)

if nargin < 5
	plotonoff = 1;
end
if nargin < 6 
	restart = [];
end
if nargin < 7 
	dtsoul = 0;
end

% constante physique (phys)
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


zs   = post.zerod;
op0d = post.z0dinput.option;
cons = post.z0dinput.cons;
geo  = post.z0dinput.geo;
profli = post.profil0d;

indcp = min(find(profli.temps >= tc));
indc = min(find(zs.temps >= profli.temps(indcp)));
t    = zs.temps(indc);
x    = profli.xli;

gammai = 2.51;
gammae = 5;

gamma      = 2.5 .* zs.tibord(indc)./zs.tebord(indc) + 2 ./ (1 - op0d.de)  -  ...
             0.5 .* log((2 .* pi.* phys.me ./ phys.mp) .* (1 + zs.tibord(indc)./zs.tebord(indc)) .* (1 - op0d.de) .^ -2); % 25.46

gammai     = max(2.51,2.5 .* zs.tibord(indc)./zs.tebord(indc) -0.5);
gammae     = gamma - gammai;

disp([gammai,gammae,gammae+gammai,gamma])

% Preparation appel de SOUL

% calcul des flux + limitations
win = (3/2)  .* trapz(x,profli.ptot(indcp,:) .* profli.vpr(indcp,:));
ntot = trapz(x,profli.nep(indcp,:) .* profli.vpr(indcp,:));
% le temps de confinement doit etre entre 1e-3 s et 1e3s
qmax = win ./ 1e-3;
qmin = win ./ 1e3;
gmax = ntot ./ 1e-3;
gmin = ntot ./ 1e3;
q_edge_e = min(qmax,max(qmin,profli.qe(indcp,end)));
q_edge_i = min(qmax,max(qmin,profli.qi(indcp,end)));
Q_perp_e = q_edge_e./  profli.vpr(indcp,end) ./ zs.dsol(indc);       % perp. input power [W/m^3]
Q_perp_i = q_edge_i./  profli.vpr(indcp,end) ./ zs.dsol(indc);       % perp. input power [W/m^3]
% les donnees ne sont pas valides
if ~isfinite(Q_perp_e) || ~isfinite(Q_perp_i)
	    save contexte_NaN_Q_soul3
	    error('NaN in Q soul3');
end

% fuelling dans le divertor
% calcul plus detailler de Smin
ge_edge = profli.ge(indcp,end) .* profli.vpr(indcp,end) .* profli.grho2r2(indcp,end);
S_perp	= min(gmax,max(1 ,max(gmin,ge_edge))) ./  profli.vpr(indcp,end) ./  zs.dsol(indc);	% perp. input particles [/m^3/s]


% these E. Tsitrone
lclim = pi .* geo.R(indc) .* zs.qa(indc);
lcpol = pi .* geo.R(indc);
%lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* 5 .* zs.qa) .^ 2);   % le 5 pour tenir compte du point X
%lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* 7 .* zs.qa) .^ 2);   % le 7 pour tenir compte du point X (valeur etalonnee sur A.S Kukushkin et al NF 43 p 716-723)
lcx = sqrt(zs.peri(indc) .^ 2  + (pi .* geo.R(indc) .* op0d.lcx .* zs.qa(indc)) .^ 2);  
switch op0d.configuration
case 0
	lc = lcpol;
case 1
	lc = lclim;
case 2
	lc  = zs.xpoint(indc) .* lcx + (~zs.xpoint(indc)) .* lcpol;
case 3
	lc  = zs.xpoint(indc) .* lcx + (~zs.xpoint(indc)) .* lclim;
otherwise
	lc  = lcx;
end


% parametre d'entree de SOUL
param.n         = 127;    % number of point along the field line
param.xinit     = 0; 		% value of normalized coordinate at the stagnation point (0)
param.xend      = 1;   		% value of normalized coordinate at the plate (1) 
param.beta      = 1;            % grid stretching parameter (beta=0->Uniform grid)	
% il n'y a pas le /2, c'est déja la longueur du point de stagnation jusqu'a la plaque du divertor.			
param.Lpar      = lc;			% length of the magnetic field line between stagnation point and the plate (m)
param.L_xpt     = x_x_fraction .* lc;  % length of the magnetic field line between stagnation point and the x-point (m)
param.kb        = phys.e;			% Boltzmann constant (1.6022d-19)   
param.atm	= zs.meff(indc); %At Mass No[H=1,D=2,T=3]
param.mp        = phys.mp;       % proton mass [kg]
param.me        = phys.me;       % electron mass [kg]
if op0d.eioniz == 0
    param.E_i      = z0eioniz_div(zs.telim(indc),zs.nelim(indc))
else
    param.E_i      = op0d.eioniz; % electron energy loss by ioniz. & excitation (eV)
end
    
param.gama_e    = gammae;      % electron (sheath) heat transmission coefficient at mach number = 0
param.gama_i    = gammai;      % ion (sheath) heat transmission coefficient at mach number = 0
param.ang       = angle_target;       % field line angle with plate (degree])
if op0d.zimp == 6
   param.xi_i      = zs.nimpm(indc) ./ zs.nem(indc);
elseif op0d.zmax == 6
   param.xi_i      = op0d.rimp .* zs.nimpm(indc) ./ zs.nem(indc);
else 
   param.xi_i      = 0;
end
param.S_perp    = S_perp ;                            % volumetric ion source from core to SOL [/m^3/s]
param.Q_perp_e  = op0d.fpower .* Q_perp_e ;   % electron volumetric heat srce from core->SOL [W/m^3]
param.Q_perp_i  = op0d.fpower .* Q_perp_i ;   % ion volumetric heat srce from core->SOL [W/m^3]
param.limiter = 'y';
param.Rcyl      = op0d.Recycling;                  % fraction of incident ion flux on the plate that are recycling in the neutral flux.
param.verbosity = plotonoff; 		% selection of verbosity level (0= no trace)
if length(profli.temps) == length(zs.temps)
    param.RTOL      = 1e-3; 		% relative tolerance for timestepping
else
    param.RTOL      = 1e-2; 		% relative tolerance for timestepping
end
% appel de la fonction
[x_out,grid_out,ne_out,vflow_out,te_out,ti_out,neutral,mach,flux,qcde,qcve,qcdi,qcvi,texte,restart] = interface_soul3(param,dtsoul,restart);
save contexte_SOUL3_from_METIS

% graphes
if plotonoff == 1    
	figure;
	subplot(4,2,1)
	plot(x_out,te_out,'b',x_out,ti_out,'r');
        hold on
        plot(min(x_out),zs.tebord(indc),'oc',max(x_out),zs.telim(indc),'ok');
        plot(min(x_out),zs.tibord(indc),'om');
	xlabel('L//')
	ylabel('eV');
	legend('Te SOUl','Ti SOUL','Te METIS LCFS','T METIS divertor','Ti METIS LCFS');
	title(sprintf('Q_{perp,e} = %g  & Q_{perp,i} = %g ',Q_perp_e,Q_perp_i));
	subplot(4,2,2)
	plot(x_out,ne_out);
        hold on
        plot(min(x_out),zs.nebord(indc),'or',max(x_out),zs.nelim(indc),'or');
	xlabel('L//')
	ylabel('density (m^-3)');
	title(sprintf('S_{perp} = %g',S_perp));
        legend('N SOUL','N Metis');
	subplot(4,2,3)
	plot(x_out,neutral);
	xlabel('L//')
	ylabel('n_0 (m^-3)');
	%title(sprintf('Gaspuff_target = %g',gaspuff_div));
	subplot(4,2,4)
	plot(x_out,mach);
	xlabel('L//')
	ylabel('mach');
	%title(sprintf('n_{scale} = %g',param.n0));
	subplot(4,2,5)
	plot(x_out,flux);
	xlabel('L//')
	ylabel('density flux (/m^2/s)');
	%title(sprintf('t_{scale} = %g',param.t0));
	subplot(4,2,6)
	plot(x_out,qcde,'b',x_out,qcve,'r',x_out,qcdi,'c',x_out,qcvi,'m');
	xlabel('L//')
	legend('e-heat cnd flux','e-heat cnv flux','i-heat cnd flux','i-heat cnv flux');
	ylabel('Cnd+Cnv flux (W/m^2)');
	title(sprintf('power_{ratio2div} * Q_{perp,e} = %g & power_{ratio2div} * Q_{perp,i} = %g',op0d.fpower .* Q_perp_e,op0d.fpower .* Q_perp_i));
	subplot(4,2,2)
        edition2
	drawnow
end
	
tau_sol = NaN;
if any(~isfinite(te_out)) ||  any(~isfinite(ti_out)) || any(~isfinite(ne_out)) || any(~isfinite(mach)) || any(te_out < 0)|| any(ti_out < 0)
		disp('ZSOUL3: non convergence');
else
	tau_sol = trapz(x_out(2:end),1./vflow_out(2:end));
	
	% structure pour la memoire
	souldata.x_out    = x_out;  
	souldata.te_out   = te_out; 
	souldata.ti_out   = ti_out; 
	souldata.ne_out   = ne_out; 
	souldata.neutral  = neutral;
	souldata.mach     = mach;   
	souldata.flux     = flux;   
	souldata.qcde     = qcde;    
	souldata.qcve     = qcve;    
	souldata.qcdi     = qcdi;    
	souldata.qcvi     = qcvi;    
	souldata.tau_sol  = tau_sol;    
	
end	
if (te_out(end) >= te_out(1))
    disp('ZSOUL3: uncorrect normalisation of density or temperature'); 	
end	

